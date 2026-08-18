/* ****************************************************************** **
**    OpenSees - Open System for Earthquake Engineering Simulation    **
**          Pacific Earthquake Engineering Research Center            **
**                                                                    **
**                                                                    **
** (C) Copyright 1999, The Regents of the University of California    **
** All Rights Reserved.                                               **
**                                                                    **
** Commercial use of this program without express permission of the   **
** University of California, Berkeley, is strictly prohibited.  See   **
** file 'COPYRIGHT'  in main directory for information on usage and   **
** redistribution,  and for a DISCLAIMER OF ALL WARRANTIES.           **
**                                                                    **
** Developed by:                                                      **
**   Frank McKenna (fmckenna@ce.berkeley.edu)                         **
**   Gregory L. Fenves (fenves@ce.berkeley.edu)                       **
**   Filip C. Filippou (filippou@ce.berkeley.edu)                     **
**                                                                    **
** ****************************************************************** */

// Description: 2-node displacement-based Timoshenko beam with uniform
// reduced integration (single mid-span Gauss point). See the header for
// the formulation notes.

#include <ASDTimoshenkoBeam3d.h>
#include <Node.h>
#include <SectionForceDeformation.h>
#include <CrdTransf.h>
#include <Damping.h>
#include <Matrix.h>
#include <Vector.h>
#include <ID.h>
#include <Renderer.h>
#include <Domain.h>
#include <string.h>
#include <Information.h>
#include <Parameter.h>
#include <Channel.h>
#include <FEM_ObjectBroker.h>
#include <ElementResponse.h>
#include <ElementalLoad.h>
#include <elementAPI.h>
#include <math.h>
#include <stdlib.h>

Matrix ASDTimoshenkoBeam3d::K(12,12);
Vector ASDTimoshenkoBeam3d::P(12);

// anonymous namespace for utilities
namespace
{
    // options for the calculateAll method
    constexpr int OPT_NONE = 0;
    constexpr int OPT_UPDATE = (1 << 0);
    constexpr int OPT_LHS = (1 << 1);
    constexpr int OPT_RHS = (1 << 2);
    constexpr int OPT_LHS_IS_INITIAL = (1 << 3);

    // uniform reduced integration: a single Gauss point at mid-span.
    // xi is the normalized coordinate in [0, 1], the weight is normalized
    // so that the physical weight is GAUSS_WEIGHT * L.
    constexpr double GAUSS_XI = 0.5;
    constexpr double GAUSS_WEIGHT = 1.0;

    // max supported section order for the preallocated workspaces
    constexpr int MAX_SECTION_ORDER = 20;

    // number of basic DOFs:
    // (axial elongation, theta_z1, theta_z2, theta_y1, theta_y2, twist)
    constexpr int NBASIC = 6;

    /** \brief ASDTimoshenkoBeam3dGlobals
     *
     * This singleton class stores the workspaces for the beam calculations,
     * statically instantiated to avoid useless re-allocations
     *
     */
    class ASDTimoshenkoBeam3dGlobals
    {
    private:
        ASDTimoshenkoBeam3dGlobals() = default;

    public:

        double B_buffer[MAX_SECTION_ORDER * NBASIC];  // storage for the strain-displacement matrix
        double e_buffer[MAX_SECTION_ORDER];           // storage for the section deformations

        Matrix kb = Matrix(NBASIC, NBASIC); // basic stiffness

    public:
        static ASDTimoshenkoBeam3dGlobals& instance() {
            static ASDTimoshenkoBeam3dGlobals _instance;
            return _instance;
        }
    };

    // Fills the L-scaled strain-displacement matrix B (order x 6) at the
    // normalized coordinate xi, mapping the basic displacements
    // v = (u, theta_z1, theta_z2, theta_y1, theta_y2, phi) to the section
    // deformations e = (1/L) B v:
    //   axial strain    = u / L                              (constant)
    //   curvature (z)   = (theta_z2 - theta_z1) / L          (constant)
    //   shear (y)       = (1-xi) theta_z1 + xi theta_z2      (linear, sampled
    //   curvature (y)   = (theta_y2 - theta_y1) / L           only at mid-span
    //   shear (z)       = (1-xi) theta_y1 + xi theta_y2       by the reduced
    //   twist           = phi / L                              integration)
    void computeBMatrix(const ID& code, double xi, double L, Matrix& B)
    {
        B.Zero();
        for (int j = 0; j < code.Size(); j++) {
            switch (code(j)) {
            case SECTION_RESPONSE_P:
                B(j, 0) = 1.0;
                break;
            case SECTION_RESPONSE_MZ:
                B(j, 1) = -1.0;
                B(j, 2) = 1.0;
                break;
            case SECTION_RESPONSE_VY:
                B(j, 1) = L * (1.0 - xi);
                B(j, 2) = L * xi;
                break;
            case SECTION_RESPONSE_MY:
                B(j, 3) = -1.0;
                B(j, 4) = 1.0;
                break;
            case SECTION_RESPONSE_VZ:
                B(j, 3) = L * (1.0 - xi);
                B(j, 4) = L * xi;
                break;
            case SECTION_RESPONSE_T:
                B(j, 5) = 1.0;
                break;
            default:
                break;
            }
        }
    }

    // The reduced-integration Timoshenko formulation needs the full
    // 3D beam response of the section: without the VY/VZ terms the
    // symmetric bending modes would be zero-energy modes, and without
    // MY/T the out-of-plane and torsional DOFs would be singular.
    bool sectionIsCompatible(SectionForceDeformation& section, bool print_error)
    {
        const ID& code = section.getType();
        bool has_P = false, has_MZ = false, has_VY = false;
        bool has_MY = false, has_VZ = false, has_T = false;
        for (int j = 0; j < code.Size(); j++) {
            switch (code(j)) {
            case SECTION_RESPONSE_P:  has_P = true; break;
            case SECTION_RESPONSE_MZ: has_MZ = true; break;
            case SECTION_RESPONSE_VY: has_VY = true; break;
            case SECTION_RESPONSE_MY: has_MY = true; break;
            case SECTION_RESPONSE_VZ: has_VZ = true; break;
            case SECTION_RESPONSE_T:  has_T = true; break;
            default: break;
            }
        }
        bool ok = has_P && has_MZ && has_VY && has_MY && has_VZ && has_T
            && (code.Size() <= MAX_SECTION_ORDER);
        if (!ok && print_error) {
            opserr << "ASDTimoshenkoBeam3d - section " << section.getTag()
                << " is not compatible: it must provide the P, MZ, MY, VY, VZ and T response codes "
                << "(use, e.g., 'section Elastic' with shear, or a section aggregated with Vy/Vz/T materials). "
                << "A section without shear stiffness would make this element singular.\n";
        }
        return ok;
    }
}

void* OPS_ASDTimoshenkoBeam3d()
{
    int ndm = OPS_GetNDM();
    int ndf = OPS_GetNDF();
    if (ndm != 3 || ndf != 6) {
        opserr << "ASDTimoshenkoBeam3d - ndm must be 3 and ndf must be 6\n";
        return 0;
    }

    if (OPS_GetNumRemainingInputArgs() < 5) {
        opserr << "insufficient arguments: element ASDTimoshenkoBeam $tag $iNode $jNode $transfTag $secTag <-mass $m> <-cMass> <-damp $dampingTag>\n";
        return 0;
    }

    // mandatory: tag, iNode, jNode, transfTag, secTag
    int iData[5];
    int numData = 5;
    if (OPS_GetIntInput(&numData, &iData[0]) < 0) {
        opserr << "ASDTimoshenkoBeam3d - invalid integer inputs\n";
        return 0;
    }

    // options
    double mass = 0.0;
    int cmass = 0;
    int dampingTag = 0;
    Damping* theDamping = 0;
    numData = 1;
    while (OPS_GetNumRemainingInputArgs() > 0) {
        const char* type = OPS_GetString();
        if (strcmp(type, "-cMass") == 0) {
            cmass = 1;
        }
        else if (strcmp(type, "-mass") == 0) {
            if (OPS_GetNumRemainingInputArgs() > 0) {
                if (OPS_GetDoubleInput(&numData, &mass) < 0) {
                    opserr << "ASDTimoshenkoBeam3d - invalid mass\n";
                    return 0;
                }
            }
        }
        else if (strcmp(type, "-damp") == 0) {
            if (OPS_GetNumRemainingInputArgs() > 0) {
                if (OPS_GetIntInput(&numData, &dampingTag) < 0) return 0;
                theDamping = OPS_getDamping(dampingTag);
                if (theDamping == 0) {
                    opserr << "ASDTimoshenkoBeam3d - damping not found\n";
                    return 0;
                }
            }
        }
    }

    CrdTransf* theTransf = OPS_getCrdTransf(iData[3]);
    if (theTransf == 0) {
        opserr << "ASDTimoshenkoBeam3d - coord transformation " << iData[3] << " not found\n";
        return 0;
    }

    SectionForceDeformation* theSection = OPS_getSectionForceDeformation(iData[4]);
    if (theSection == 0) {
        opserr << "ASDTimoshenkoBeam3d - section " << iData[4] << " not found\n";
        return 0;
    }
    if (!sectionIsCompatible(*theSection, true))
        return 0;

    return new ASDTimoshenkoBeam3d(iData[0], iData[1], iData[2],
        *theSection, *theTransf, mass, cmass, theDamping);
}

ASDTimoshenkoBeam3d::ASDTimoshenkoBeam3d(int tag, int nd1, int nd2,
                                         SectionForceDeformation &section,
                                         CrdTransf &coordTransf, double r, int cm,
                                         Damping *damping)
    :Element(tag, ELE_TAG_ASDTimoshenkoBeam3d),
     theSection(0), crdTransf(0), theDamping(0),
     connectedExternalNodes(2),
     Q(12), q(6), rho(r), cMass(cm), parameterID(0)
{
    theSection = section.getCopy();
    if (theSection == 0) {
        opserr << "ASDTimoshenkoBeam3d::ASDTimoshenkoBeam3d - failed to get a copy of section model\n";
        exit(-1);
    }

    crdTransf = coordTransf.getCopy3d();
    if (crdTransf == 0) {
        opserr << "ASDTimoshenkoBeam3d::ASDTimoshenkoBeam3d - failed to copy coordinate transformation\n";
        exit(-1);
    }

    if (damping) {
        theDamping = damping->getCopy();
        if (theDamping == 0) {
            opserr << "ASDTimoshenkoBeam3d::ASDTimoshenkoBeam3d - failed to copy damping\n";
            exit(-1);
        }
    }

    connectedExternalNodes(0) = nd1;
    connectedExternalNodes(1) = nd2;

    theNodes[0] = 0;
    theNodes[1] = 0;

    for (int i = 0; i < 5; i++) {
        q0[i] = 0.0;
        p0[i] = 0.0;
    }
}

ASDTimoshenkoBeam3d::ASDTimoshenkoBeam3d()
    :Element(0, ELE_TAG_ASDTimoshenkoBeam3d),
     theSection(0), crdTransf(0), theDamping(0),
     connectedExternalNodes(2),
     Q(12), q(6), rho(0.0), cMass(0), parameterID(0)
{
    theNodes[0] = 0;
    theNodes[1] = 0;

    for (int i = 0; i < 5; i++) {
        q0[i] = 0.0;
        p0[i] = 0.0;
    }
}

ASDTimoshenkoBeam3d::~ASDTimoshenkoBeam3d()
{
    if (theSection)
        delete theSection;

    if (crdTransf)
        delete crdTransf;

    if (theDamping)
        delete theDamping;
}

int
ASDTimoshenkoBeam3d::getNumExternalNodes() const
{
    return 2;
}

const ID&
ASDTimoshenkoBeam3d::getExternalNodes()
{
    return connectedExternalNodes;
}

Node **
ASDTimoshenkoBeam3d::getNodePtrs()
{
    return theNodes;
}

int
ASDTimoshenkoBeam3d::getNumDOF()
{
    return 12;
}

void
ASDTimoshenkoBeam3d::setDomain(Domain *theDomain)
{
    // Check Domain is not null - invoked when object removed from a domain
    if (theDomain == 0) {
        theNodes[0] = 0;
        theNodes[1] = 0;
        return;
    }

    int Nd1 = connectedExternalNodes(0);
    int Nd2 = connectedExternalNodes(1);

    theNodes[0] = theDomain->getNode(Nd1);
    theNodes[1] = theDomain->getNode(Nd2);

    if (theNodes[0] == 0 || theNodes[1] == 0) {
        opserr << "WARNING ASDTimoshenkoBeam3d (tag: " << this->getTag() << "), node not found in domain\n";
        return;
    }

    int dofNd1 = theNodes[0]->getNumberDOF();
    int dofNd2 = theNodes[1]->getNumberDOF();

    if (dofNd1 != 6 || dofNd2 != 6) {
        opserr << "WARNING ASDTimoshenkoBeam3d (tag: " << this->getTag() << "), needs 6 DOFs at each node\n";
        return;
    }

    // check the section provides the shear response (the broker path
    // bypasses the check done in the OPS_ function)
    if (theSection != 0)
        sectionIsCompatible(*theSection, true);

    if (crdTransf->initialize(theNodes[0], theNodes[1])) {
        opserr << "WARNING ASDTimoshenkoBeam3d (tag: " << this->getTag() << "), failed to initialize coordinate transformation\n";
        return;
    }

    // initialize the damping
    if (theDamping && theDamping->setDomain(theDomain, 6)) {
        opserr << "ASDTimoshenkoBeam3d::setDomain - error initializing damping\n";
        exit(0);
    }

    double L = crdTransf->getInitialLength();
    if (L == 0.0) {
        opserr << "WARNING ASDTimoshenkoBeam3d (tag: " << this->getTag() << "), zero length\n";
        return;
    }

    this->DomainComponent::setDomain(theDomain);

    this->update();
}

int
ASDTimoshenkoBeam3d::setDamping(Domain *theDomain, Damping *damping)
{
    if (theDomain && damping) {
        if (theDamping)
            delete theDamping;

        theDamping = damping->getCopy();
        if (!theDamping) {
            opserr << "ASDTimoshenkoBeam3d::setDamping - failed to get copy of damping\n";
            return -1;
        }
        if (theDamping->setDomain(theDomain, 6)) {
            opserr << "ASDTimoshenkoBeam3d::setDamping - error initializing damping\n";
            return -2;
        }
    }

    return 0;
}

int
ASDTimoshenkoBeam3d::commitState()
{
    int retVal = 0;

    // call element commitState to do any base class stuff
    if ((retVal = this->Element::commitState()) != 0) {
        opserr << "ASDTimoshenkoBeam3d::commitState () - failed in base class";
    }

    retVal += theSection->commitState();
    retVal += crdTransf->commitState();
    if (theDamping)
        retVal += theDamping->commitState();

    return retVal;
}

int
ASDTimoshenkoBeam3d::revertToLastCommit()
{
    int retVal = 0;

    retVal += theSection->revertToLastCommit();
    retVal += crdTransf->revertToLastCommit();
    if (theDamping)
        retVal += theDamping->revertToLastCommit();

    return retVal;
}

int
ASDTimoshenkoBeam3d::revertToStart()
{
    int retVal = 0;

    retVal += theSection->revertToStart();
    retVal += crdTransf->revertToStart();
    if (theDamping)
        retVal += theDamping->revertToStart();

    return retVal;
}

int
ASDTimoshenkoBeam3d::calculateAll(Matrix &kb, Vector &qb, int options)
{
    int result = 0;

    auto& globals = ASDTimoshenkoBeam3dGlobals::instance();

    double L = crdTransf->getInitialLength();
    double oneOverL = 1.0 / L;

    // the L-scaled strain-displacement matrix at the single Gauss point
    int order = theSection->getOrder();
    const ID &code = theSection->getType();
    Matrix B(globals.B_buffer, order, NBASIC);
    computeBMatrix(code, GAUSS_XI, L, B);

    // impose the trial section deformations: e = (1/L) B v
    if (options & OPT_UPDATE) {
        crdTransf->update();
        const Vector &v = crdTransf->getBasicTrialDisp();
        Vector e(globals.e_buffer, order);
        e.addMatrixVector(0.0, B, v, oneOverL);
        result += theSection->setTrialSectionDeformation(e);
    }

    // basic stiffness: kb = (w/L) B^T ks B
    if (options & OPT_LHS) {
        const Matrix &ks = (options & OPT_LHS_IS_INITIAL) ?
            theSection->getInitialTangent() : theSection->getSectionTangent();
        kb.addMatrixTripleProduct(0.0, B, ks, GAUSS_WEIGHT * oneOverL);
    }

    // basic forces: q = w B^T s + q0
    if (options & OPT_RHS) {
        const Vector &s = theSection->getStressResultant();
        qb.addMatrixTransposeVector(0.0, B, s, GAUSS_WEIGHT);
        qb(0) += q0[0];
        qb(1) += q0[1];
        qb(2) += q0[2];
        qb(3) += q0[3];
        qb(4) += q0[4];
    }

    return result;
}

int
ASDTimoshenkoBeam3d::update(void)
{
    auto& kb = ASDTimoshenkoBeam3dGlobals::instance().kb;
    int err = calculateAll(kb, q, OPT_UPDATE);
    if (err != 0) {
        opserr << "ASDTimoshenkoBeam3d::update() - failed setTrialSectionDeformation()\n";
        return err;
    }
    return 0;
}

const Matrix&
ASDTimoshenkoBeam3d::getTangentStiff()
{
    auto& kb = ASDTimoshenkoBeam3dGlobals::instance().kb;

    // the basic forces are needed by the transformation for the
    // geometric stiffness terms (PDelta, Corotational)
    calculateAll(kb, q, OPT_LHS | OPT_RHS);

    // Transform to global stiffness
    K = crdTransf->getGlobalStiffMatrix(kb, q);

    return K;
}

const Matrix&
ASDTimoshenkoBeam3d::getInitialStiff()
{
    auto& kb = ASDTimoshenkoBeam3dGlobals::instance().kb;

    calculateAll(kb, q, OPT_LHS | OPT_LHS_IS_INITIAL);
    if (theDamping)
        kb *= theDamping->getStiffnessMultiplier();

    // Transform to global stiffness
    K = crdTransf->getInitialGlobalStiffMatrix(kb);

    return K;
}

const Matrix&
ASDTimoshenkoBeam3d::getMass()
{
    K.Zero();

    if (rho == 0.0)
        return K;

    double L = crdTransf->getInitialLength();
    if (cMass == 0) {
        // lumped mass matrix
        double m = 0.5 * rho * L;
        K(0,0) = K(1,1) = K(2,2) = K(6,6) = K(7,7) = K(8,8) = m;
    }
    else {
        // consistent mass matrix, based on the same linear interpolation
        // used for the displacements (no rotary/torsional inertia)
        static Matrix ml(12,12);
        ml.Zero();
        double m = rho * L / 6.0;
        for (int i = 0; i < 3; i++) {
            ml(i,i) = ml(i+6,i+6) = 2.0 * m;
            ml(i,i+6) = ml(i+6,i) = m;
        }

        // transform local mass matrix to global system
        K = crdTransf->getGlobalMatrixFromLocal(ml);
    }

    return K;
}

void
ASDTimoshenkoBeam3d::zeroLoad(void)
{
    Q.Zero();

    for (int i = 0; i < 5; i++) {
        q0[i] = 0.0;
        p0[i] = 0.0;
    }

    return;
}

int
ASDTimoshenkoBeam3d::addLoad(ElementalLoad *theLoad, double loadFactor)
{
    int type;
    const Vector &data = theLoad->getData(type, loadFactor);
    double L = crdTransf->getInitialLength();

    // Note: p0 (reactions) comes from statics and is the same as in any
    // other beam; q0 (fixed-end basic forces) is consistent with the
    // linear interpolation of this element: a transverse load does no
    // work on the basic rotations, so the moment entries are zero.
    if (type == LOAD_TAG_Beam3dUniformLoad) {
        double wy = data(0) * loadFactor;  // Transverse y
        double wz = data(1) * loadFactor;  // Transverse z
        double wx = data(2) * loadFactor;  // Axial (+ve from node I to J)

        double Vy = 0.5 * wy * L;
        double Vz = 0.5 * wz * L;
        double Pax = wx * L;

        // Reactions in basic system
        p0[0] -= Pax;
        p0[1] -= Vy;
        p0[2] -= Vy;
        p0[3] -= Vz;
        p0[4] -= Vz;

        // Fixed end forces in basic system
        q0[0] -= 0.5 * Pax;
    }
    else if (type == LOAD_TAG_Beam3dPointLoad) {
        double Py = data(0) * loadFactor;
        double Pz = data(1) * loadFactor;
        double N  = data(2) * loadFactor;
        double aOverL = data(3);

        if (aOverL < 0.0 || aOverL > 1.0)
            return 0;

        // Reactions in basic system
        p0[0] -= N;
        p0[1] -= Py * (1.0 - aOverL);
        p0[2] -= Py * aOverL;
        p0[3] -= Pz * (1.0 - aOverL);
        p0[4] -= Pz * aOverL;

        // Fixed end forces in basic system
        q0[0] -= N * aOverL;
    }
    else {
        opserr << "ASDTimoshenkoBeam3d::addLoad() - load type unknown for element with tag: "
            << this->getTag() << "\n";
        return -1;
    }

    return 0;
}

int
ASDTimoshenkoBeam3d::addInertiaLoadToUnbalance(const Vector &accel)
{
    // Check for a quick return
    if (rho == 0.0)
        return 0;

    // Get R * accel from the nodes
    const Vector &Raccel1 = theNodes[0]->getRV(accel);
    const Vector &Raccel2 = theNodes[1]->getRV(accel);

    if (6 != Raccel1.Size() || 6 != Raccel2.Size()) {
        opserr << "ASDTimoshenkoBeam3d::addInertiaLoadToUnbalance matrix and vector sizes are incompatible\n";
        return -1;
    }

    // want to add ( - fact * M R * accel ) to unbalance
    if (cMass == 0) {
        // take advantage of lumped mass matrix
        double L = crdTransf->getInitialLength();
        double m = 0.5 * rho * L;

        Q(0) -= m * Raccel1(0);
        Q(1) -= m * Raccel1(1);
        Q(2) -= m * Raccel1(2);
        Q(6) -= m * Raccel2(0);
        Q(7) -= m * Raccel2(1);
        Q(8) -= m * Raccel2(2);
    }
    else {
        // use matrix vector multip. for consistent mass matrix
        static Vector Raccel(12);
        for (int i = 0; i < 6; i++) {
            Raccel(i)   = Raccel1(i);
            Raccel(i+6) = Raccel2(i);
        }
        Q.addMatrixVector(1.0, this->getMass(), Raccel, -1.0);
    }

    return 0;
}

const Vector&
ASDTimoshenkoBeam3d::getResistingForce()
{
    auto& kb = ASDTimoshenkoBeam3dGlobals::instance().kb;

    calculateAll(kb, q, OPT_RHS);

    if (theDamping)
        theDamping->update(q);

    // Vector for reactions in basic system
    Vector p0Vec(p0, 5);

    P = crdTransf->getGlobalResistingForce(q, p0Vec);

    // Subtract other external nodal loads ... P_res = P_int - P_ext
    if (rho != 0)
        P.addVector(1.0, Q, -1.0);

    return P;
}

const Vector&
ASDTimoshenkoBeam3d::getDampingForce(void)
{
    crdTransf->update();

    return crdTransf->getGlobalResistingForce(theDamping->getDampingForce(), Vector(5));
}

const Vector&
ASDTimoshenkoBeam3d::getResistingForceIncInertia()
{
    P = this->getResistingForce();

    if (theDamping)
        P += this->getDampingForce();

    if (rho != 0.0) {
        const Vector &accel1 = theNodes[0]->getTrialAccel();
        const Vector &accel2 = theNodes[1]->getTrialAccel();

        if (cMass == 0) {
            // take advantage of lumped mass matrix
            double L = crdTransf->getInitialLength();
            double m = 0.5 * rho * L;

            P(0) += m * accel1(0);
            P(1) += m * accel1(1);
            P(2) += m * accel1(2);
            P(6) += m * accel2(0);
            P(7) += m * accel2(1);
            P(8) += m * accel2(2);
        }
        else {
            // use matrix vector multip. for consistent mass matrix
            static Vector accel(12);
            for (int i = 0; i < 6; i++) {
                accel(i)   = accel1(i);
                accel(i+6) = accel2(i);
            }
            P.addMatrixVector(1.0, this->getMass(), accel, 1.0);
        }

        // add the damping forces if rayleigh damping
        if (alphaM != 0.0 || betaK != 0.0 || betaK0 != 0.0 || betaKc != 0.0)
            P.addVector(1.0, this->getRayleighDampingForces(), 1.0);
    }
    else {
        // add the damping forces if rayleigh damping
        if (betaK != 0.0 || betaK0 != 0.0 || betaKc != 0.0)
            P.addVector(1.0, this->getRayleighDampingForces(), 1.0);
    }

    return P;
}

int
ASDTimoshenkoBeam3d::sendSelf(int commitTag, Channel &theChannel)
{
    int dbTag = this->getDbTag();

    static Vector data(15);
    data(0) = this->getTag();
    data(1) = connectedExternalNodes(0);
    data(2) = connectedExternalNodes(1);
    data(3) = crdTransf->getClassTag();
    int crdTransfDbTag = crdTransf->getDbTag();
    if (crdTransfDbTag == 0) {
        crdTransfDbTag = theChannel.getDbTag();
        if (crdTransfDbTag != 0)
            crdTransf->setDbTag(crdTransfDbTag);
    }
    data(4) = crdTransfDbTag;
    data(5) = theSection->getClassTag();
    int sectDbTag = theSection->getDbTag();
    if (sectDbTag == 0) {
        sectDbTag = theChannel.getDbTag();
        if (sectDbTag != 0)
            theSection->setDbTag(sectDbTag);
    }
    data(6) = sectDbTag;
    data(7) = rho;
    data(8) = cMass;
    data(9) = alphaM;
    data(10) = betaK;
    data(11) = betaK0;
    data(12) = betaKc;

    data(13) = 0;
    data(14) = 0;
    if (theDamping) {
        data(13) = theDamping->getClassTag();
        int dampingDbTag = theDamping->getDbTag();
        if (dampingDbTag == 0) {
            dampingDbTag = theChannel.getDbTag();
            if (dampingDbTag != 0)
                theDamping->setDbTag(dampingDbTag);
        }
        data(14) = dampingDbTag;
    }

    if (theChannel.sendVector(dbTag, commitTag, data) < 0) {
        opserr << "ASDTimoshenkoBeam3d::sendSelf() - failed to send data Vector\n";
        return -1;
    }

    // send the coordinate transformation
    if (crdTransf->sendSelf(commitTag, theChannel) < 0) {
        opserr << "ASDTimoshenkoBeam3d::sendSelf() - failed to send crdTransf\n";
        return -1;
    }

    // send the section
    if (theSection->sendSelf(commitTag, theChannel) < 0) {
        opserr << "ASDTimoshenkoBeam3d::sendSelf() - failed to send section\n";
        return -1;
    }

    // ask the Damping to send itself
    if (theDamping && theDamping->sendSelf(commitTag, theChannel) < 0) {
        opserr << "ASDTimoshenkoBeam3d::sendSelf() - could not send Damping\n";
        return -1;
    }

    return 0;
}

int
ASDTimoshenkoBeam3d::recvSelf(int commitTag, Channel &theChannel,
                              FEM_ObjectBroker &theBroker)
{
    int dbTag = this->getDbTag();

    static Vector data(15);

    if (theChannel.recvVector(dbTag, commitTag, data) < 0) {
        opserr << "ASDTimoshenkoBeam3d::recvSelf() - failed to recv data Vector\n";
        return -1;
    }

    this->setTag((int)data(0));
    connectedExternalNodes(0) = (int)data(1);
    connectedExternalNodes(1) = (int)data(2);
    int crdTransfClassTag = (int)data(3);
    int crdTransfDbTag = (int)data(4);
    int sectClassTag = (int)data(5);
    int sectDbTag = (int)data(6);

    rho = data(7);
    cMass = (int)data(8);

    alphaM = data(9);
    betaK = data(10);
    betaK0 = data(11);
    betaKc = data(12);

    // create a new crdTransf object if one needed
    if (crdTransf == 0 || crdTransf->getClassTag() != crdTransfClassTag) {
        if (crdTransf != 0)
            delete crdTransf;

        crdTransf = theBroker.getNewCrdTransf(crdTransfClassTag);

        if (crdTransf == 0) {
            opserr << "ASDTimoshenkoBeam3d::recvSelf() - failed to obtain a CrdTransf object with classTag "
                << crdTransfClassTag << "\n";
            return -2;
        }
    }
    crdTransf->setDbTag(crdTransfDbTag);

    // invoke recvSelf on the crdTransf object
    if (crdTransf->recvSelf(commitTag, theChannel, theBroker) < 0) {
        opserr << "ASDTimoshenkoBeam3d::recvSelf() - failed to recv crdTransf\n";
        return -3;
    }

    // create a new section object if one needed
    if (theSection == 0 || theSection->getClassTag() != sectClassTag) {
        if (theSection != 0)
            delete theSection;

        theSection = theBroker.getNewSection(sectClassTag);

        if (theSection == 0) {
            opserr << "ASDTimoshenkoBeam3d::recvSelf() - Broker could not create Section of class type "
                << sectClassTag << "\n";
            return -4;
        }
    }
    theSection->setDbTag(sectDbTag);

    // invoke recvSelf on the section object
    if (theSection->recvSelf(commitTag, theChannel, theBroker) < 0) {
        opserr << "ASDTimoshenkoBeam3d::recvSelf() - section failed to recv itself\n";
        return -5;
    }

    // check if the Damping is null; if so, get a new one
    int dmpTag = (int)data(13);
    if (dmpTag) {
        if (theDamping == 0 || theDamping->getClassTag() != dmpTag) {
            if (theDamping != 0)
                delete theDamping;
            theDamping = theBroker.getNewDamping(dmpTag);
            if (theDamping == 0) {
                opserr << "ASDTimoshenkoBeam3d::recvSelf() - could not get a Damping\n";
                return -6;
            }
        }

        // now, receive the Damping
        theDamping->setDbTag((int)data(14));
        if (theDamping->recvSelf(commitTag, theChannel, theBroker) < 0) {
            opserr << "ASDTimoshenkoBeam3d::recvSelf() - could not receive Damping\n";
            return -7;
        }
    }
    else {
        if (theDamping) {
            delete theDamping;
            theDamping = 0;
        }
    }

    return 0;
}

void
ASDTimoshenkoBeam3d::Print(OPS_Stream &s, int flag)
{
    if (flag == OPS_PRINT_CURRENTSTATE) {
        s << "\nASDTimoshenkoBeam3d, element id:  " << this->getTag() << "\n";
        s << "\tConnected external nodes:  " << connectedExternalNodes;
        s << "\tCoordTransf: " << crdTransf->getTag() << "\n";
        s << "\tmass density:  " << rho << ", cMass: " << cMass << "\n";

        double L = crdTransf->getInitialLength();
        double oneOverL = 1.0 / L;
        double N  = q(0);
        double Mz1 = q(1);
        double Mz2 = q(2);
        double Vy = (Mz1 + Mz2) * oneOverL;
        double My1 = q(3);
        double My2 = q(4);
        double Vz = (My1 + My2) * oneOverL;
        double T  = q(5);

        s << "\tEnd 1 Forces (P Mz Vy My Vz T): "
            << -N + p0[0] << " " << Mz1 << " " << Vy + p0[1] << " "
            << My1 << " " << -Vz + p0[3] << " " << -T << "\n";
        s << "\tEnd 2 Forces (P Mz Vy My Vz T): "
            << N << " " << Mz2 << " " << -Vy + p0[2] << " "
            << My2 << " " << Vz + p0[4] << " " << T << "\n";

        theSection->Print(s, flag);
    }

    if (flag == OPS_PRINT_PRINTMODEL_JSON) {
        s << "\t\t\t{";
        s << "\"name\": " << this->getTag() << ", ";
        s << "\"type\": \"ASDTimoshenkoBeam3d\", ";
        s << "\"nodes\": [" << connectedExternalNodes(0) << ", " << connectedExternalNodes(1) << "], ";
        s << "\"section\": \"" << theSection->getTag() << "\", ";
        s << "\"massperlength\": " << rho << ", ";
        s << "\"crdTransformation\": \"" << crdTransf->getTag() << "\"}";
    }
}

int
ASDTimoshenkoBeam3d::displaySelf(Renderer &theViewer, int displayMode, float fact,
                                 const char **displayModes, int numModes)
{
    static Vector v1(3);
    static Vector v2(3);

    theNodes[0]->getDisplayCrds(v1, fact, displayMode);
    theNodes[1]->getDisplayCrds(v2, fact, displayMode);

    return theViewer.drawLine(v1, v2, 1.0, 1.0, this->getTag());
}

Response*
ASDTimoshenkoBeam3d::setResponse(const char **argv, int argc,
                                 OPS_Stream &output)
{
    Response *theResponse = 0;

    output.tag("ElementOutput");
    output.attr("eleType", "ASDTimoshenkoBeam3d");
    output.attr("eleTag", this->getTag());
    output.attr("node1", connectedExternalNodes[0]);
    output.attr("node2", connectedExternalNodes[1]);

    // global force -
    if (strcmp(argv[0],"forces") == 0 || strcmp(argv[0],"force") == 0
        || strcmp(argv[0],"globalForce") == 0 || strcmp(argv[0],"globalForces") == 0) {

        output.tag("ResponseType","Px_1");
        output.tag("ResponseType","Py_1");
        output.tag("ResponseType","Pz_1");
        output.tag("ResponseType","Mx_1");
        output.tag("ResponseType","My_1");
        output.tag("ResponseType","Mz_1");
        output.tag("ResponseType","Px_2");
        output.tag("ResponseType","Py_2");
        output.tag("ResponseType","Pz_2");
        output.tag("ResponseType","Mx_2");
        output.tag("ResponseType","My_2");
        output.tag("ResponseType","Mz_2");

        theResponse = new ElementResponse(this, 1, P);
    }

    // local force -
    else if (strcmp(argv[0],"localForce") == 0 || strcmp(argv[0],"localForces") == 0) {

        output.tag("ResponseType","N_1");
        output.tag("ResponseType","Vy_1");
        output.tag("ResponseType","Vz_1");
        output.tag("ResponseType","T_1");
        output.tag("ResponseType","My_1");
        output.tag("ResponseType","Mz_1");
        output.tag("ResponseType","N_2");
        output.tag("ResponseType","Vy_2");
        output.tag("ResponseType","Vz_2");
        output.tag("ResponseType","T_2");
        output.tag("ResponseType","My_2");
        output.tag("ResponseType","Mz_2");

        theResponse = new ElementResponse(this, 2, P);
    }

    // basic force -
    else if (strcmp(argv[0],"basicForce") == 0 || strcmp(argv[0],"basicForces") == 0) {

        output.tag("ResponseType","N");
        output.tag("ResponseType","Mz_1");
        output.tag("ResponseType","Mz_2");
        output.tag("ResponseType","My_1");
        output.tag("ResponseType","My_2");
        output.tag("ResponseType","T");

        theResponse = new ElementResponse(this, 9, Vector(6));
    }

    // global damping force -
    else if (theDamping && (strcmp(argv[0],"globalDampingForce") == 0 || strcmp(argv[0],"globalDampingForces") == 0)) {

        output.tag("ResponseType","Px_1");
        output.tag("ResponseType","Py_1");
        output.tag("ResponseType","Pz_1");
        output.tag("ResponseType","Mx_1");
        output.tag("ResponseType","My_1");
        output.tag("ResponseType","Mz_1");
        output.tag("ResponseType","Px_2");
        output.tag("ResponseType","Py_2");
        output.tag("ResponseType","Pz_2");
        output.tag("ResponseType","Mx_2");
        output.tag("ResponseType","My_2");
        output.tag("ResponseType","Mz_2");

        theResponse = new ElementResponse(this, 21, P);
    }

    // local damping force -
    else if (theDamping && (strcmp(argv[0],"localDampingForce") == 0 || strcmp(argv[0],"localDampingForces") == 0)) {

        output.tag("ResponseType","N_1");
        output.tag("ResponseType","Vy_1");
        output.tag("ResponseType","Vz_1");
        output.tag("ResponseType","T_1");
        output.tag("ResponseType","My_1");
        output.tag("ResponseType","Mz_1");
        output.tag("ResponseType","N_2");
        output.tag("ResponseType","Vy_2");
        output.tag("ResponseType","Vz_2");
        output.tag("ResponseType","T_2");
        output.tag("ResponseType","My_2");
        output.tag("ResponseType","Mz_2");

        theResponse = new ElementResponse(this, 22, P);
    }

    // basic damping force -
    else if (theDamping && (strcmp(argv[0],"basicDampingForce") == 0 || strcmp(argv[0],"basicDampingForces") == 0)) {

        output.tag("ResponseType","N");
        output.tag("ResponseType","Mz_1");
        output.tag("ResponseType","Mz_2");
        output.tag("ResponseType","My_1");
        output.tag("ResponseType","My_2");
        output.tag("ResponseType","T");

        theResponse = new ElementResponse(this, 23, Vector(6));
    }

    // basic stiffness -
    else if (strcmp(argv[0],"basicStiffness") == 0) {

        output.tag("ResponseType","N");
        output.tag("ResponseType","Mz_1");
        output.tag("ResponseType","Mz_2");
        output.tag("ResponseType","My_1");
        output.tag("ResponseType","My_2");
        output.tag("ResponseType","T");

        theResponse = new ElementResponse(this, 19, Matrix(6,6));
    }

    // chord rotation -
    else if (strcmp(argv[0],"chordRotation") == 0 || strcmp(argv[0],"chordDeformation") == 0
             || strcmp(argv[0],"basicDeformation") == 0) {

        output.tag("ResponseType","eps");
        output.tag("ResponseType","thetaZ_1");
        output.tag("ResponseType","thetaZ_2");
        output.tag("ResponseType","thetaY_1");
        output.tag("ResponseType","thetaY_2");
        output.tag("ResponseType","thetaX");

        theResponse = new ElementResponse(this, 3, Vector(6));
    }

    // plastic rotation -
    else if (strcmp(argv[0],"plasticRotation") == 0 || strcmp(argv[0],"plasticDeformation") == 0) {

        output.tag("ResponseType","epsP");
        output.tag("ResponseType","thetaZP_1");
        output.tag("ResponseType","thetaZP_2");
        output.tag("ResponseType","thetaYP_1");
        output.tag("ResponseType","thetaYP_2");
        output.tag("ResponseType","thetaXP");

        theResponse = new ElementResponse(this, 4, Vector(6));
    }

    else if (strcmp(argv[0],"RayleighForces") == 0 ||
             strcmp(argv[0],"rayleighForces") == 0 ||
             strcmp(argv[0],"dampingForces") == 0) {

        theResponse = new ElementResponse(this, 12, P);
    }

    // section response -
    else if (strcmp(argv[0],"sectionX") == 0) {
        if (argc > 2) {
            output.tag("GaussPointOutput");
            output.attr("number", 1);
            output.attr("eta", GAUSS_XI * crdTransf->getInitialLength());
            theResponse = theSection->setResponse(&argv[2], argc-2, output);
            output.endTag();
        }
    }
    else if (strcmp(argv[0],"section") == 0) {
        if (argc > 1) {
            int sectionNum = atoi(argv[1]);
            if (sectionNum == 1 && argc > 2) {
                // by number: only 1 is valid (single Gauss point)
                output.tag("GaussPointOutput");
                output.attr("number", 1);
                output.attr("eta", GAUSS_XI * crdTransf->getInitialLength());
                theResponse = theSection->setResponse(&argv[2], argc-2, output);
                output.endTag();
            }
            else if (sectionNum == 0) {
                // argv[1] was not an int: forward everything to the section
                output.tag("GaussPointOutput");
                output.attr("number", 1);
                output.attr("eta", GAUSS_XI * crdTransf->getInitialLength());
                theResponse = theSection->setResponse(&argv[1], argc-1, output);
                output.endTag();
            }
        }
    }

    else if (strcmp(argv[0],"integrationPoints") == 0)
        theResponse = new ElementResponse(this, 7, Vector(1));

    else if (strcmp(argv[0],"integrationWeights") == 0)
        theResponse = new ElementResponse(this, 8, Vector(1));

    else if (strcmp(argv[0],"sectionTags") == 0)
        theResponse = new ElementResponse(this, 110, ID(1));

    else if (strcmp(argv[0],"energy") == 0)
        theResponse = new ElementResponse(this, 10, 0.0);

    if (theResponse == 0)
        theResponse = crdTransf->setResponse(argv, argc, output);

    output.endTag();

    if (theResponse == 0)
        return Element::setResponse(argv, argc, output);
    else
        return theResponse;
}

int
ASDTimoshenkoBeam3d::getResponse(int responseID, Information &eleInfo)
{
    double N, V, M1, M2, T;
    double L = crdTransf->getInitialLength();
    double oneOverL = 1.0 / L;

    if (responseID == 1)
        return eleInfo.setVector(this->getResistingForce());

    else if (responseID == 12) {
        P.Zero();
        P.addVector(1.0, this->getRayleighDampingForces(), 1.0);
        return eleInfo.setVector(P);
    }

    else if (responseID == 2) {
        // Axial
        N = q(0);
        P(6) =  N;
        P(0) = -N + p0[0];

        // Torsion
        T = q(5);
        P(9) =  T;
        P(3) = -T;

        // Moments about z and shears along y
        M1 = q(1);
        M2 = q(2);
        P(5)  = M1;
        P(11) = M2;
        V = (M1 + M2) * oneOverL;
        P(1) =  V + p0[1];
        P(7) = -V + p0[2];

        // Moments about y and shears along z
        M1 = q(3);
        M2 = q(4);
        P(4)  = M1;
        P(10) = M2;
        V = (M1 + M2) * oneOverL;
        P(2) = -V + p0[3];
        P(8) =  V + p0[4];

        return eleInfo.setVector(P);
    }

    else if (responseID == 9)
        return eleInfo.setVector(q);

    else if (responseID == 19) {
        static Matrix kb(6,6);
        calculateAll(kb, q, OPT_LHS);
        return eleInfo.setMatrix(kb);
    }

    else if (responseID == 21)
        return eleInfo.setVector(this->getDampingForce());

    else if (responseID == 22) {
        Vector Sd(6);
        Sd = theDamping->getDampingForce();

        // Axial
        N = Sd(0);
        P(6) =  N;
        P(0) = -N;

        // Torsion
        T = Sd(5);
        P(9) =  T;
        P(3) = -T;

        // Moments about z and shears along y
        M1 = Sd(1);
        M2 = Sd(2);
        P(5)  = M1;
        P(11) = M2;
        V = (M1 + M2) * oneOverL;
        P(1) =  V;
        P(7) = -V;

        // Moments about y and shears along z
        M1 = Sd(3);
        M2 = Sd(4);
        P(4)  = M1;
        P(10) = M2;
        V = (M1 + M2) * oneOverL;
        P(2) = -V;
        P(8) =  V;

        return eleInfo.setVector(P);
    }

    else if (responseID == 23)
        return eleInfo.setVector(theDamping->getDampingForce());

    // Chord rotation
    else if (responseID == 3)
        return eleInfo.setVector(crdTransf->getBasicTrialDisp());

    // Plastic rotation
    else if (responseID == 4) {
        static Vector vp(6);
        static Vector ve(6);
        static Matrix kb(6,6);
        calculateAll(kb, q, OPT_LHS | OPT_LHS_IS_INITIAL);
        kb.Solve(q, ve);
        vp = crdTransf->getBasicTrialDisp();
        vp -= ve;
        return eleInfo.setVector(vp);
    }

    else if (responseID == 7) {
        Vector locs(1);
        locs(0) = GAUSS_XI * L;
        return eleInfo.setVector(locs);
    }

    else if (responseID == 8) {
        Vector weights(1);
        weights(0) = GAUSS_WEIGHT * L;
        return eleInfo.setVector(weights);
    }

    else if (responseID == 110) {
        ID tags(1);
        tags(0) = theSection->getTag();
        return eleInfo.setID(tags);
    }

    else if (responseID == 10)
        return eleInfo.setDouble(theSection->getEnergy() * GAUSS_WEIGHT * L);

    else
        return Element::getResponse(responseID, eleInfo);
}

int
ASDTimoshenkoBeam3d::setParameter(const char **argv, int argc, Parameter &param)
{
    if (argc < 1)
        return -1;

    // If the parameter belongs to the element itself
    if (strcmp(argv[0],"rho") == 0) {
        param.setValue(rho);
        return param.addObject(1, this);
    }

    // damping
    if (strstr(argv[0],"damp") != 0) {
        if (argc < 2 || !theDamping)
            return -1;
        return theDamping->setParameter(&argv[1], argc-1, param);
    }

    // section response
    if (strstr(argv[0],"sectionX") != 0) {
        if (argc < 3)
            return -1;
        return theSection->setParameter(&argv[2], argc-2, param);
    }

    if (strstr(argv[0],"section") != 0) {
        if (argc < 3)
            return -1;
        // the section number must be 1 (single Gauss point)
        int sectionNum = atoi(argv[1]);
        if (sectionNum == 1)
            return theSection->setParameter(&argv[2], argc-2, param);
        return -1;
    }

    // default: delegate to the section
    return theSection->setParameter(argv, argc, param);
}

int
ASDTimoshenkoBeam3d::updateParameter(int parameterID, Information &info)
{
    if (parameterID == 1) {
        rho = info.theDouble;
        return 0;
    }
    return -1;
}

int
ASDTimoshenkoBeam3d::activateParameter(int passedParameterID)
{
    parameterID = passedParameterID;
    return 0;
}
