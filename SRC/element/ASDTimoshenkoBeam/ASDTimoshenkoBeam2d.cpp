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

#include <ASDTimoshenkoBeam2d.h>
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

Matrix ASDTimoshenkoBeam2d::K(6,6);
Vector ASDTimoshenkoBeam2d::P(6);

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

    // number of basic DOFs: (axial elongation, theta_1, theta_2)
    constexpr int NBASIC = 3;

    /** \brief ASDTimoshenkoBeam2dGlobals
     *
     * This singleton class stores the workspaces for the beam calculations,
     * statically instantiated to avoid useless re-allocations
     *
     */
    class ASDTimoshenkoBeam2dGlobals
    {
    private:
        ASDTimoshenkoBeam2dGlobals() = default;

    public:

        double B_buffer[MAX_SECTION_ORDER * NBASIC];  // storage for the strain-displacement matrix
        double e_buffer[MAX_SECTION_ORDER];           // storage for the section deformations

        Matrix kb = Matrix(NBASIC, NBASIC); // basic stiffness

    public:
        static ASDTimoshenkoBeam2dGlobals& instance() {
            static ASDTimoshenkoBeam2dGlobals _instance;
            return _instance;
        }
    };

    // Fills the L-scaled strain-displacement matrix B (order x 3) at the
    // normalized coordinate xi, mapping the basic displacements
    // v = (u, theta_1, theta_2) to the section deformations e = (1/L) B v:
    //   axial strain     = u / L                          (constant)
    //   curvature        = (theta_2 - theta_1) / L        (constant)
    //   shear strain     = (1-xi) theta_1 + xi theta_2    (linear, sampled
    //                      only at mid-span by the reduced integration)
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
            default:
                break;
            }
        }
    }

    // The reduced-integration Timoshenko formulation needs the shear
    // response of the section: without the VY term the symmetric bending
    // mode (theta_1 = theta_2) would be a zero-energy mode.
    bool sectionIsCompatible(SectionForceDeformation& section, bool print_error)
    {
        const ID& code = section.getType();
        bool has_P = false, has_MZ = false, has_VY = false;
        for (int j = 0; j < code.Size(); j++) {
            switch (code(j)) {
            case SECTION_RESPONSE_P:  has_P = true; break;
            case SECTION_RESPONSE_MZ: has_MZ = true; break;
            case SECTION_RESPONSE_VY: has_VY = true; break;
            default: break;
            }
        }
        bool ok = has_P && has_MZ && has_VY && (code.Size() <= MAX_SECTION_ORDER);
        if (!ok && print_error) {
            opserr << "ASDTimoshenkoBeam2d - section " << section.getTag()
                << " is not compatible: it must provide the P, MZ and VY response codes "
                << "(use, e.g., 'section Elastic' with shear, or a section aggregated with a Vy material). "
                << "A section without shear stiffness would make this element singular.\n";
        }
        return ok;
    }
}

void* OPS_ASDTimoshenkoBeam2d()
{
    int ndm = OPS_GetNDM();
    int ndf = OPS_GetNDF();
    if (ndm != 2 || ndf != 3) {
        opserr << "ASDTimoshenkoBeam2d - ndm must be 2 and ndf must be 3\n";
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
        opserr << "ASDTimoshenkoBeam2d - invalid integer inputs\n";
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
                    opserr << "ASDTimoshenkoBeam2d - invalid mass\n";
                    return 0;
                }
            }
        }
        else if (strcmp(type, "-damp") == 0) {
            if (OPS_GetNumRemainingInputArgs() > 0) {
                if (OPS_GetIntInput(&numData, &dampingTag) < 0) return 0;
                theDamping = OPS_getDamping(dampingTag);
                if (theDamping == 0) {
                    opserr << "ASDTimoshenkoBeam2d - damping not found\n";
                    return 0;
                }
            }
        }
    }

    CrdTransf* theTransf = OPS_getCrdTransf(iData[3]);
    if (theTransf == 0) {
        opserr << "ASDTimoshenkoBeam2d - coord transformation " << iData[3] << " not found\n";
        return 0;
    }

    SectionForceDeformation* theSection = OPS_getSectionForceDeformation(iData[4]);
    if (theSection == 0) {
        opserr << "ASDTimoshenkoBeam2d - section " << iData[4] << " not found\n";
        return 0;
    }
    if (!sectionIsCompatible(*theSection, true))
        return 0;

    return new ASDTimoshenkoBeam2d(iData[0], iData[1], iData[2],
        *theSection, *theTransf, mass, cmass, theDamping);
}

ASDTimoshenkoBeam2d::ASDTimoshenkoBeam2d(int tag, int nd1, int nd2,
                                         SectionForceDeformation &section,
                                         CrdTransf &coordTransf, double r, int cm,
                                         Damping *damping)
    :Element(tag, ELE_TAG_ASDTimoshenkoBeam2d),
     theSection(0), crdTransf(0), theDamping(0),
     connectedExternalNodes(2),
     Q(6), q(3), rho(r), cMass(cm), parameterID(0)
{
    theSection = section.getCopy();
    if (theSection == 0) {
        opserr << "ASDTimoshenkoBeam2d::ASDTimoshenkoBeam2d - failed to get a copy of section model\n";
        exit(-1);
    }

    crdTransf = coordTransf.getCopy2d();
    if (crdTransf == 0) {
        opserr << "ASDTimoshenkoBeam2d::ASDTimoshenkoBeam2d - failed to copy coordinate transformation\n";
        exit(-1);
    }

    if (damping) {
        theDamping = damping->getCopy();
        if (theDamping == 0) {
            opserr << "ASDTimoshenkoBeam2d::ASDTimoshenkoBeam2d - failed to copy damping\n";
            exit(-1);
        }
    }

    connectedExternalNodes(0) = nd1;
    connectedExternalNodes(1) = nd2;

    theNodes[0] = 0;
    theNodes[1] = 0;

    q0[0] = 0.0;
    q0[1] = 0.0;
    q0[2] = 0.0;

    p0[0] = 0.0;
    p0[1] = 0.0;
    p0[2] = 0.0;
}

ASDTimoshenkoBeam2d::ASDTimoshenkoBeam2d()
    :Element(0, ELE_TAG_ASDTimoshenkoBeam2d),
     theSection(0), crdTransf(0), theDamping(0),
     connectedExternalNodes(2),
     Q(6), q(3), rho(0.0), cMass(0), parameterID(0)
{
    theNodes[0] = 0;
    theNodes[1] = 0;

    q0[0] = 0.0;
    q0[1] = 0.0;
    q0[2] = 0.0;

    p0[0] = 0.0;
    p0[1] = 0.0;
    p0[2] = 0.0;
}

ASDTimoshenkoBeam2d::~ASDTimoshenkoBeam2d()
{
    if (theSection)
        delete theSection;

    if (crdTransf)
        delete crdTransf;

    if (theDamping)
        delete theDamping;
}

int
ASDTimoshenkoBeam2d::getNumExternalNodes() const
{
    return 2;
}

const ID&
ASDTimoshenkoBeam2d::getExternalNodes()
{
    return connectedExternalNodes;
}

Node **
ASDTimoshenkoBeam2d::getNodePtrs()
{
    return theNodes;
}

int
ASDTimoshenkoBeam2d::getNumDOF()
{
    return 6;
}

void
ASDTimoshenkoBeam2d::setDomain(Domain *theDomain)
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
        opserr << "WARNING ASDTimoshenkoBeam2d (tag: " << this->getTag() << "), node not found in domain\n";
        return;
    }

    int dofNd1 = theNodes[0]->getNumberDOF();
    int dofNd2 = theNodes[1]->getNumberDOF();

    if (dofNd1 != 3 || dofNd2 != 3) {
        opserr << "WARNING ASDTimoshenkoBeam2d (tag: " << this->getTag() << "), needs 3 DOFs at each node\n";
        return;
    }

    // check the section provides the shear response (the broker path
    // bypasses the check done in the OPS_ function)
    if (theSection != 0)
        sectionIsCompatible(*theSection, true);

    if (crdTransf->initialize(theNodes[0], theNodes[1])) {
        opserr << "WARNING ASDTimoshenkoBeam2d (tag: " << this->getTag() << "), failed to initialize coordinate transformation\n";
        return;
    }

    // initialize the damping
    if (theDamping && theDamping->setDomain(theDomain, 3)) {
        opserr << "ASDTimoshenkoBeam2d::setDomain - error initializing damping\n";
        exit(0);
    }

    double L = crdTransf->getInitialLength();
    if (L == 0.0) {
        opserr << "WARNING ASDTimoshenkoBeam2d (tag: " << this->getTag() << "), zero length\n";
        return;
    }

    this->DomainComponent::setDomain(theDomain);

    this->update();
}

int
ASDTimoshenkoBeam2d::setDamping(Domain *theDomain, Damping *damping)
{
    if (theDomain && damping) {
        if (theDamping)
            delete theDamping;

        theDamping = damping->getCopy();
        if (!theDamping) {
            opserr << "ASDTimoshenkoBeam2d::setDamping - failed to get copy of damping\n";
            return -1;
        }
        if (theDamping->setDomain(theDomain, 3)) {
            opserr << "ASDTimoshenkoBeam2d::setDamping - error initializing damping\n";
            return -2;
        }
    }

    return 0;
}

int
ASDTimoshenkoBeam2d::commitState()
{
    int retVal = 0;

    // call element commitState to do any base class stuff
    if ((retVal = this->Element::commitState()) != 0) {
        opserr << "ASDTimoshenkoBeam2d::commitState () - failed in base class";
    }

    retVal += theSection->commitState();
    retVal += crdTransf->commitState();
    if (theDamping)
        retVal += theDamping->commitState();

    return retVal;
}

int
ASDTimoshenkoBeam2d::revertToLastCommit()
{
    int retVal = 0;

    retVal += theSection->revertToLastCommit();
    retVal += crdTransf->revertToLastCommit();
    if (theDamping)
        retVal += theDamping->revertToLastCommit();

    return retVal;
}

int
ASDTimoshenkoBeam2d::revertToStart()
{
    int retVal = 0;

    retVal += theSection->revertToStart();
    retVal += crdTransf->revertToStart();
    if (theDamping)
        retVal += theDamping->revertToStart();

    return retVal;
}

int
ASDTimoshenkoBeam2d::calculateAll(Matrix &kb, Vector &qb, int options)
{
    int result = 0;

    auto& globals = ASDTimoshenkoBeam2dGlobals::instance();

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
    }

    return result;
}

int
ASDTimoshenkoBeam2d::update(void)
{
    auto& kb = ASDTimoshenkoBeam2dGlobals::instance().kb;
    int err = calculateAll(kb, q, OPT_UPDATE);
    if (err != 0) {
        opserr << "ASDTimoshenkoBeam2d::update() - failed setTrialSectionDeformation()\n";
        return err;
    }
    return 0;
}

const Matrix&
ASDTimoshenkoBeam2d::getTangentStiff()
{
    auto& kb = ASDTimoshenkoBeam2dGlobals::instance().kb;

    // the basic forces are needed by the transformation for the
    // geometric stiffness terms (PDelta, Corotational)
    calculateAll(kb, q, OPT_LHS | OPT_RHS);

    // Transform to global stiffness
    K = crdTransf->getGlobalStiffMatrix(kb, q);

    return K;
}

const Matrix&
ASDTimoshenkoBeam2d::getInitialStiff()
{
    auto& kb = ASDTimoshenkoBeam2dGlobals::instance().kb;

    calculateAll(kb, q, OPT_LHS | OPT_LHS_IS_INITIAL);
    if (theDamping)
        kb *= theDamping->getStiffnessMultiplier();

    // Transform to global stiffness
    K = crdTransf->getInitialGlobalStiffMatrix(kb);

    return K;
}

const Matrix&
ASDTimoshenkoBeam2d::getMass()
{
    K.Zero();

    if (rho == 0.0)
        return K;

    double L = crdTransf->getInitialLength();
    if (cMass == 0) {
        // lumped mass matrix
        double m = 0.5 * rho * L;
        K(0,0) = K(1,1) = K(3,3) = K(4,4) = m;
    }
    else {
        // consistent mass matrix, based on the same linear interpolation
        // used for the displacements (no rotary inertia)
        static Matrix ml(6,6);
        ml.Zero();
        double m = rho * L / 6.0;
        ml(0,0) = ml(1,1) = ml(3,3) = ml(4,4) = 2.0 * m;
        ml(0,3) = ml(3,0) = ml(1,4) = ml(4,1) = m;

        // transform local mass matrix to global system
        K = crdTransf->getGlobalMatrixFromLocal(ml);
    }

    return K;
}

void
ASDTimoshenkoBeam2d::zeroLoad(void)
{
    Q.Zero();

    q0[0] = 0.0;
    q0[1] = 0.0;
    q0[2] = 0.0;

    p0[0] = 0.0;
    p0[1] = 0.0;
    p0[2] = 0.0;

    return;
}

int
ASDTimoshenkoBeam2d::addLoad(ElementalLoad *theLoad, double loadFactor)
{
    int type;
    const Vector &data = theLoad->getData(type, loadFactor);
    double L = crdTransf->getInitialLength();

    // Note: p0 (reactions) comes from statics and is the same as in any
    // other beam; q0 (fixed-end basic forces) is consistent with the
    // linear interpolation of this element: a transverse load does no
    // work on the basic rotations, so the moment entries are zero.
    if (type == LOAD_TAG_Beam2dUniformLoad) {
        double wt = data(0) * loadFactor;  // Transverse (+ve upward)
        double wa = data(1) * loadFactor;  // Axial (+ve from node I to J)

        double V = 0.5 * wt * L;
        double Pax = wa * L;

        // Reactions in basic system
        p0[0] -= Pax;
        p0[1] -= V;
        p0[2] -= V;

        // Fixed end forces in basic system
        q0[0] -= 0.5 * Pax;
    }
    else if (type == LOAD_TAG_BeamUniformMoment) {
        double mz = data(2) * loadFactor;  // About z, per unit length

        // Reactions in basic system
        p0[1] += mz;
        p0[2] -= mz;

        // Fixed end forces in basic system: the distributed moment is
        // work-conjugate to the (linearly interpolated) rotation field
        q0[1] -= 0.5 * mz * L;
        q0[2] -= 0.5 * mz * L;
    }
    else if (type == LOAD_TAG_Beam2dPointLoad) {
        double Pt = data(0) * loadFactor;
        double N = data(1) * loadFactor;
        double aOverL = data(2);

        if (aOverL < 0.0 || aOverL > 1.0)
            return 0;

        // Reactions in basic system
        p0[0] -= N;
        p0[1] -= Pt * (1.0 - aOverL);
        p0[2] -= Pt * aOverL;

        // Fixed end forces in basic system
        q0[0] -= N * aOverL;
    }
    else {
        opserr << "ASDTimoshenkoBeam2d::addLoad() - load type unknown for element with tag: "
            << this->getTag() << "\n";
        return -1;
    }

    return 0;
}

int
ASDTimoshenkoBeam2d::addInertiaLoadToUnbalance(const Vector &accel)
{
    // Check for a quick return
    if (rho == 0.0)
        return 0;

    // Get R * accel from the nodes
    const Vector &Raccel1 = theNodes[0]->getRV(accel);
    const Vector &Raccel2 = theNodes[1]->getRV(accel);

    if (3 != Raccel1.Size() || 3 != Raccel2.Size()) {
        opserr << "ASDTimoshenkoBeam2d::addInertiaLoadToUnbalance matrix and vector sizes are incompatible\n";
        return -1;
    }

    // want to add ( - fact * M R * accel ) to unbalance
    if (cMass == 0) {
        // take advantage of lumped mass matrix
        double L = crdTransf->getInitialLength();
        double m = 0.5 * rho * L;

        Q(0) -= m * Raccel1(0);
        Q(1) -= m * Raccel1(1);
        Q(3) -= m * Raccel2(0);
        Q(4) -= m * Raccel2(1);
    }
    else {
        // use matrix vector multip. for consistent mass matrix
        static Vector Raccel(6);
        for (int i = 0; i < 3; i++) {
            Raccel(i)   = Raccel1(i);
            Raccel(i+3) = Raccel2(i);
        }
        Q.addMatrixVector(1.0, this->getMass(), Raccel, -1.0);
    }

    return 0;
}

const Vector&
ASDTimoshenkoBeam2d::getResistingForce()
{
    auto& kb = ASDTimoshenkoBeam2dGlobals::instance().kb;

    calculateAll(kb, q, OPT_RHS);

    if (theDamping)
        theDamping->update(q);

    // Vector for reactions in basic system
    Vector p0Vec(p0, 3);

    P = crdTransf->getGlobalResistingForce(q, p0Vec);

    // Subtract other external nodal loads ... P_res = P_int - P_ext
    if (rho != 0)
        P.addVector(1.0, Q, -1.0);

    return P;
}

const Vector&
ASDTimoshenkoBeam2d::getDampingForce(void)
{
    crdTransf->update();

    return crdTransf->getGlobalResistingForce(theDamping->getDampingForce(), Vector(3));
}

const Vector&
ASDTimoshenkoBeam2d::getResistingForceIncInertia()
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
            P(3) += m * accel2(0);
            P(4) += m * accel2(1);
        }
        else {
            // use matrix vector multip. for consistent mass matrix
            static Vector accel(6);
            for (int i = 0; i < 3; i++) {
                accel(i)   = accel1(i);
                accel(i+3) = accel2(i);
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
ASDTimoshenkoBeam2d::sendSelf(int commitTag, Channel &theChannel)
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
        opserr << "ASDTimoshenkoBeam2d::sendSelf() - failed to send data Vector\n";
        return -1;
    }

    // send the coordinate transformation
    if (crdTransf->sendSelf(commitTag, theChannel) < 0) {
        opserr << "ASDTimoshenkoBeam2d::sendSelf() - failed to send crdTransf\n";
        return -1;
    }

    // send the section
    if (theSection->sendSelf(commitTag, theChannel) < 0) {
        opserr << "ASDTimoshenkoBeam2d::sendSelf() - failed to send section\n";
        return -1;
    }

    // ask the Damping to send itself
    if (theDamping && theDamping->sendSelf(commitTag, theChannel) < 0) {
        opserr << "ASDTimoshenkoBeam2d::sendSelf() - could not send Damping\n";
        return -1;
    }

    return 0;
}

int
ASDTimoshenkoBeam2d::recvSelf(int commitTag, Channel &theChannel,
                              FEM_ObjectBroker &theBroker)
{
    int dbTag = this->getDbTag();

    static Vector data(15);

    if (theChannel.recvVector(dbTag, commitTag, data) < 0) {
        opserr << "ASDTimoshenkoBeam2d::recvSelf() - failed to recv data Vector\n";
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
            opserr << "ASDTimoshenkoBeam2d::recvSelf() - failed to obtain a CrdTransf object with classTag "
                << crdTransfClassTag << "\n";
            return -2;
        }
    }
    crdTransf->setDbTag(crdTransfDbTag);

    // invoke recvSelf on the crdTransf object
    if (crdTransf->recvSelf(commitTag, theChannel, theBroker) < 0) {
        opserr << "ASDTimoshenkoBeam2d::recvSelf() - failed to recv crdTransf\n";
        return -3;
    }

    // create a new section object if one needed
    if (theSection == 0 || theSection->getClassTag() != sectClassTag) {
        if (theSection != 0)
            delete theSection;

        theSection = theBroker.getNewSection(sectClassTag);

        if (theSection == 0) {
            opserr << "ASDTimoshenkoBeam2d::recvSelf() - Broker could not create Section of class type "
                << sectClassTag << "\n";
            return -4;
        }
    }
    theSection->setDbTag(sectDbTag);

    // invoke recvSelf on the section object
    if (theSection->recvSelf(commitTag, theChannel, theBroker) < 0) {
        opserr << "ASDTimoshenkoBeam2d::recvSelf() - section failed to recv itself\n";
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
                opserr << "ASDTimoshenkoBeam2d::recvSelf() - could not get a Damping\n";
                return -6;
            }
        }

        // now, receive the Damping
        theDamping->setDbTag((int)data(14));
        if (theDamping->recvSelf(commitTag, theChannel, theBroker) < 0) {
            opserr << "ASDTimoshenkoBeam2d::recvSelf() - could not receive Damping\n";
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
ASDTimoshenkoBeam2d::Print(OPS_Stream &s, int flag)
{
    if (flag == OPS_PRINT_CURRENTSTATE) {
        s << "\nASDTimoshenkoBeam2d, element id:  " << this->getTag() << "\n";
        s << "\tConnected external nodes:  " << connectedExternalNodes;
        s << "\tCoordTransf: " << crdTransf->getTag() << "\n";
        s << "\tmass density:  " << rho << ", cMass: " << cMass << "\n";

        double L = crdTransf->getInitialLength();
        double Pax = q(0);
        double M1 = q(1);
        double M2 = q(2);
        double V = (M1 + M2) / L;

        s << "\tEnd 1 Forces (P V M): " << -Pax + p0[0]
            << " " << V + p0[1] << " " << M1 << "\n";
        s << "\tEnd 2 Forces (P V M): " << Pax
            << " " << -V + p0[2] << " " << M2 << "\n";

        theSection->Print(s, flag);
    }

    if (flag == OPS_PRINT_PRINTMODEL_JSON) {
        s << "\t\t\t{";
        s << "\"name\": " << this->getTag() << ", ";
        s << "\"type\": \"ASDTimoshenkoBeam2d\", ";
        s << "\"nodes\": [" << connectedExternalNodes(0) << ", " << connectedExternalNodes(1) << "], ";
        s << "\"section\": \"" << theSection->getTag() << "\", ";
        s << "\"massperlength\": " << rho << ", ";
        s << "\"crdTransformation\": \"" << crdTransf->getTag() << "\"}";
    }
}

int
ASDTimoshenkoBeam2d::displaySelf(Renderer &theViewer, int displayMode, float fact,
                                 const char **displayModes, int numModes)
{
    static Vector v1(3);
    static Vector v2(3);

    theNodes[0]->getDisplayCrds(v1, fact, displayMode);
    theNodes[1]->getDisplayCrds(v2, fact, displayMode);

    return theViewer.drawLine(v1, v2, 1.0, 1.0, this->getTag());
}

Response*
ASDTimoshenkoBeam2d::setResponse(const char **argv, int argc,
                                 OPS_Stream &output)
{
    Response *theResponse = 0;

    output.tag("ElementOutput");
    output.attr("eleType", "ASDTimoshenkoBeam2d");
    output.attr("eleTag", this->getTag());
    output.attr("node1", connectedExternalNodes[0]);
    output.attr("node2", connectedExternalNodes[1]);

    // global force -
    if (strcmp(argv[0],"forces") == 0 || strcmp(argv[0],"force") == 0
        || strcmp(argv[0],"globalForce") == 0 || strcmp(argv[0],"globalForces") == 0) {

        output.tag("ResponseType","Px_1");
        output.tag("ResponseType","Py_1");
        output.tag("ResponseType","Mz_1");
        output.tag("ResponseType","Px_2");
        output.tag("ResponseType","Py_2");
        output.tag("ResponseType","Mz_2");

        theResponse = new ElementResponse(this, 1, P);
    }

    // local force -
    else if (strcmp(argv[0],"localForce") == 0 || strcmp(argv[0],"localForces") == 0) {

        output.tag("ResponseType","N1");
        output.tag("ResponseType","V1");
        output.tag("ResponseType","M1");
        output.tag("ResponseType","N2");
        output.tag("ResponseType","V2");
        output.tag("ResponseType","M2");

        theResponse = new ElementResponse(this, 2, P);
    }

    // basic force -
    else if (strcmp(argv[0],"basicForce") == 0 || strcmp(argv[0],"basicForces") == 0) {

        output.tag("ResponseType","N");
        output.tag("ResponseType","M1");
        output.tag("ResponseType","M2");

        theResponse = new ElementResponse(this, 9, Vector(3));
    }

    // global damping force -
    else if (theDamping && (strcmp(argv[0],"globalDampingForce") == 0 || strcmp(argv[0],"globalDampingForces") == 0)) {

        output.tag("ResponseType","Px_1");
        output.tag("ResponseType","Py_1");
        output.tag("ResponseType","Mz_1");
        output.tag("ResponseType","Px_2");
        output.tag("ResponseType","Py_2");
        output.tag("ResponseType","Mz_2");

        theResponse = new ElementResponse(this, 21, P);
    }

    // local damping force -
    else if (theDamping && (strcmp(argv[0],"localDampingForce") == 0 || strcmp(argv[0],"localDampingForces") == 0)) {

        output.tag("ResponseType","N1");
        output.tag("ResponseType","V1");
        output.tag("ResponseType","M1");
        output.tag("ResponseType","N2");
        output.tag("ResponseType","V2");
        output.tag("ResponseType","M2");

        theResponse = new ElementResponse(this, 22, P);
    }

    // basic damping force -
    else if (theDamping && (strcmp(argv[0],"basicDampingForce") == 0 || strcmp(argv[0],"basicDampingForces") == 0)) {

        output.tag("ResponseType","N");
        output.tag("ResponseType","M1");
        output.tag("ResponseType","M2");

        theResponse = new ElementResponse(this, 23, Vector(3));
    }

    // basic stiffness -
    else if (strcmp(argv[0],"basicStiffness") == 0) {

        output.tag("ResponseType","N");
        output.tag("ResponseType","M1");
        output.tag("ResponseType","M2");

        theResponse = new ElementResponse(this, 19, Matrix(3,3));
    }

    // chord rotation -
    else if (strcmp(argv[0],"chordRotation") == 0 || strcmp(argv[0],"chordDeformation") == 0
             || strcmp(argv[0],"basicDeformation") == 0) {

        output.tag("ResponseType","eps");
        output.tag("ResponseType","theta1");
        output.tag("ResponseType","theta2");

        theResponse = new ElementResponse(this, 3, Vector(3));
    }

    // plastic rotation -
    else if (strcmp(argv[0],"plasticRotation") == 0 || strcmp(argv[0],"plasticDeformation") == 0) {

        output.tag("ResponseType","epsP");
        output.tag("ResponseType","theta1P");
        output.tag("ResponseType","theta2P");

        theResponse = new ElementResponse(this, 4, Vector(3));
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
ASDTimoshenkoBeam2d::getResponse(int responseID, Information &eleInfo)
{
    double V;
    double L = crdTransf->getInitialLength();

    if (responseID == 1)
        return eleInfo.setVector(this->getResistingForce());

    else if (responseID == 12) {
        P.Zero();
        P.addVector(1.0, this->getRayleighDampingForces(), 1.0);
        return eleInfo.setVector(P);
    }

    else if (responseID == 2) {
        P(3) =  q(0);
        P(0) = -q(0) + p0[0];
        P(2) = q(1);
        P(5) = q(2);
        V = (q(1) + q(2)) / L;
        P(1) =  V + p0[1];
        P(4) = -V + p0[2];
        return eleInfo.setVector(P);
    }

    else if (responseID == 9)
        return eleInfo.setVector(q);

    else if (responseID == 19) {
        static Matrix kb(3,3);
        calculateAll(kb, q, OPT_LHS);
        return eleInfo.setMatrix(kb);
    }

    else if (responseID == 21)
        return eleInfo.setVector(this->getDampingForce());

    else if (responseID == 22) {
        Vector Sd(3);
        Sd = theDamping->getDampingForce();
        P(3) =  Sd(0);
        P(0) = -Sd(0);
        P(2) = Sd(1);
        P(5) = Sd(2);
        V = (Sd(1) + Sd(2)) / L;
        P(1) =  V;
        P(4) = -V;
        return eleInfo.setVector(P);
    }

    else if (responseID == 23)
        return eleInfo.setVector(theDamping->getDampingForce());

    // Chord rotation
    else if (responseID == 3)
        return eleInfo.setVector(crdTransf->getBasicTrialDisp());

    // Plastic rotation
    else if (responseID == 4) {
        static Vector vp(3);
        static Vector ve(3);
        static Matrix kb(3,3);
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
ASDTimoshenkoBeam2d::setParameter(const char **argv, int argc, Parameter &param)
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
ASDTimoshenkoBeam2d::updateParameter(int parameterID, Information &info)
{
    if (parameterID == 1) {
        rho = info.theDouble;
        return 0;
    }
    return -1;
}

int
ASDTimoshenkoBeam2d::activateParameter(int passedParameterID)
{
    parameterID = passedParameterID;
    return 0;
}
