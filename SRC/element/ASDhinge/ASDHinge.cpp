/* ****************************************************************** **
**    OpenSees - Open System for Earthquake Engineering Simulation    **
**          Pacific Earthquake Engineering Research Center            **
** ****************************************************************** */

// $Revision: 1.0 $
// $Date: 2026/08/23 $

// Original implementation: Massimo Petracca (ASDEA)

#include "ASDHinge.h"
#include "ASDHingeCorotationalTransformation.h"

#include <Domain.h>
#include <Node.h>
#include <Channel.h>
#include <FEM_ObjectBroker.h>
#include <Renderer.h>
#include <Information.h>
#include <Parameter.h>
#include <ElementResponse.h>
#include <ElementalLoad.h>
#include <UniaxialMaterial.h>
#include <Damping.h>
#include <DummyStream.h>
#include <elementAPI.h>
#include <classTags.h>

#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <vector>

// tolerance on the element length, relative to the node coordinates
#define ASDHINGE_LENTOL 1.0e-6

Matrix ASDHinge::s_K(12, 12);
Vector ASDHinge::s_P(12);

namespace {

const char* const DOF_NAMES[6] = { "Ux", "Uy", "Uz", "Rx", "Ry", "Rz" };
// chosen among the names already in the component_map of
// beam_force_deformation_traits in STKO/odb_rules/OpenSees.json, so that the
// post-processing side is a handful of lines
const char* const FORCE_NAMES[6] = { "N", "Vy", "Vz", "T", "My", "Mz" };

/// same construction as ZeroLength::setUp: z = x cross yp, y = z cross x
bool buildOrientation(const double* vx, const double* vyp, double R0[3][3])
{
    using namespace ASDHingeUtils;
    double x[3] = { vx[0], vx[1], vx[2] };
    double nx = norm3(x);
    if (nx < 1.0e-12)
        return false;
    for (int i = 0; i < 3; ++i) x[i] /= nx;

    double yp[3] = { vyp[0], vyp[1], vyp[2] };
    double z[3];
    cross3(x, yp, z);
    double nz = norm3(z);
    if (nz < 1.0e-12)
        return false;
    for (int i = 0; i < 3; ++i) z[i] /= nz;

    double y[3];
    cross3(z, x, y);

    for (int i = 0; i < 3; ++i) {
        R0[i][0] = x[i];
        R0[i][1] = y[i];
        R0[i][2] = z[i];
    }
    return true;
}

} // namespace


// --------------------------------------------------------------------------
// the command
// --------------------------------------------------------------------------

void* OPS_ASDHinge(void)
{
    static bool first_done = false;
    if (!first_done) {
        opserr << "Using ASDHinge - Developed by: ASDEA Software Technology\n";
        first_done = true;
    }

    const char* usage =
        "element ASDHinge $tag $iNode $jNode "
        "<-mat $m1 ... -dir $d1 ...> <-K $k1 ... $k6> "
        "<-orient $x1 $x2 $x3 $y1 $y2 $y3> "
        "<-corotational> <-doRayleigh> <-damp $tag>\n";

    if (OPS_GetNumRemainingInputArgs() < 3) {
        opserr << "ASDHinge: too few arguments.\nWant: " << usage;
        return 0;
    }

    int idata[3];
    int numdata = 3;
    if (OPS_GetIntInput(&numdata, idata) < 0) {
        opserr << "ASDHinge: failed to read tag, iNode, jNode.\n";
        return 0;
    }

    std::vector<int> matTags;
    std::vector<int> dirs;
    double kdiag[6] = { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 };
    bool haveK[6] = { false, false, false, false, false, false };
    double R0[3][3];
    ASDHingeUtils::identity3(R0);
    bool corotational = false;
    int doRayleigh = 0;
    int dampingTag = 0;

    while (OPS_GetNumRemainingInputArgs() > 0) {
        const char* type = OPS_GetString();

        if (strcmp(type, "-mat") == 0) {
            while (OPS_GetNumRemainingInputArgs() > 0) {
                int mtag;
                int one = 1;
                int nleft = OPS_GetNumRemainingInputArgs();
                if (OPS_GetIntInput(&one, &mtag) < 0) {
                    if (nleft > OPS_GetNumRemainingInputArgs())
                        OPS_ResetCurrentInputArg(-1);
                    break;
                }
                matTags.push_back(mtag);
            }
        }
        else if (strcmp(type, "-dir") == 0 || strcmp(type, "-dof") == 0) {
            while (OPS_GetNumRemainingInputArgs() > 0) {
                int d;
                int one = 1;
                int nleft = OPS_GetNumRemainingInputArgs();
                if (OPS_GetIntInput(&one, &d) < 0) {
                    if (nleft > OPS_GetNumRemainingInputArgs())
                        OPS_ResetCurrentInputArg(-1);
                    break;
                }
                dirs.push_back(d);
            }
        }
        else if (strcmp(type, "-K") == 0) {
            int six = 6;
            if (OPS_GetNumRemainingInputArgs() < 6 ||
                OPS_GetDoubleInput(&six, kdiag) < 0) {
                opserr << "ASDHinge: -K wants exactly 6 values, "
                          "one per local dof (0 = free).\n";
                return 0;
            }
            for (int i = 0; i < 6; ++i)
                haveK[i] = (kdiag[i] != 0.0);
        }
        else if (strcmp(type, "-orient") == 0) {
            double v[6];
            int six = 6;
            if (OPS_GetNumRemainingInputArgs() < 6 ||
                OPS_GetDoubleInput(&six, v) < 0) {
                opserr << "ASDHinge: -orient wants 6 values, "
                          "the local x and the local xy plane vectors.\n";
                return 0;
            }
            if (!buildOrientation(v, v + 3, R0)) {
                opserr << "ASDHinge: -orient vectors are null or parallel.\n";
                return 0;
            }
        }
        else if (strcmp(type, "-corotational") == 0) {
            corotational = true;
        }
        else if (strcmp(type, "-doRayleigh") == 0) {
            // zeroLength accepts an optional integer here, and HARD FAILS if
            // the next token is not one -- so "-doRayleigh -orient ..." is a
            // parse error there.  Take the integer only if it really is one.
            doRayleigh = 1;
            if (OPS_GetNumRemainingInputArgs() > 0) {
                int flag;
                int one = 1;
                int nleft = OPS_GetNumRemainingInputArgs();
                if (OPS_GetIntInput(&one, &flag) < 0) {
                    if (nleft > OPS_GetNumRemainingInputArgs())
                        OPS_ResetCurrentInputArg(-1);
                }
                else {
                    doRayleigh = flag;
                }
            }
        }
        else if (strcmp(type, "-damp") == 0) {
            int one = 1;
            if (OPS_GetNumRemainingInputArgs() < 1 ||
                OPS_GetIntInput(&one, &dampingTag) < 0) {
                opserr << "ASDHinge: -damp wants a damping tag.\n";
                return 0;
            }
        }
        else {
            opserr << "ASDHinge: unknown option '" << type << "'.\nWant: " << usage;
            return 0;
        }
    }

    if (matTags.size() != dirs.size()) {
        opserr << "ASDHinge: -mat has " << (int)matTags.size()
               << " entries and -dir has " << (int)dirs.size()
               << ": they must match one to one.\n";
        return 0;
    }

    // build the kernel: six slots, free unless told otherwise
    ASDHingeUncoupledKernel* kernel = new ASDHingeUncoupledKernel();
    for (int i = 0; i < 6; ++i) {
        if (haveK[i])
            kernel->setLinear(i, kdiag[i]);
        else
            kernel->setFree(i);
    }

    for (size_t n = 0; n < matTags.size(); ++n) {
        int d = dirs[n];
        if (d < 1 || d > 6) {
            opserr << "ASDHinge: -dir " << d << " is out of range. The local "
                      "dofs are 1..6 (Ux Uy Uz Rx Ry Rz).\n";
            delete kernel;
            return 0;
        }
        int slot = d - 1;
        if (haveK[slot]) {
            opserr << "ASDHinge: local dof " << d << " got both a material and "
                      "a non zero -K. Pick one.\n";
            delete kernel;
            return 0;
        }
        if (kernel->getSlotState(slot) == ASDHingeUncoupledKernel::Material) {
            opserr << "ASDHinge: local dof " << d
                   << " was given a material twice.\n";
            delete kernel;
            return 0;
        }
        UniaxialMaterial* mat = OPS_getUniaxialMaterial(matTags[n]);
        if (mat == 0) {
            opserr << "ASDHinge: no uniaxialMaterial with tag "
                   << matTags[n] << ".\n";
            delete kernel;
            return 0;
        }
        if (mat->getClassTag() != MAT_TAG_ASDHysteretic1DMaterial) {
            opserr << "ASDHinge: material " << matTags[n] << " is a "
                   << mat->getClassType() << ". ASDHinge only accepts "
                      "ASDHysteretic1D: the IMPL-EX error codes and the hinge "
                      "responses (damage, freeEnergy, limitStateRatio) are "
                      "that material's. Use a zeroLength for anything else.\n";
            delete kernel;
            return 0;
        }
        kernel->setMaterial(slot, mat);
    }
    kernel->initialize();

    Damping* theDamping = 0;
    if (dampingTag != 0) {
        theDamping = OPS_getDamping(dampingTag);
        if (theDamping == 0) {
            opserr << "ASDHinge: damping " << dampingTag << " not found.\n";
            delete kernel;
            return 0;
        }
    }

    return new ASDHinge(idata[0], idata[1], idata[2], R0, kernel,
                        corotational, doRayleigh, theDamping);
}


// --------------------------------------------------------------------------
// construction
// --------------------------------------------------------------------------

ASDHinge::ASDHinge(int tag, int Nd1, int Nd2, const double R0[3][3],
                   ASDHingeKernel* kernel, bool corotational,
                   int doRayleigh, Damping* damping)
    : Element(tag, ELE_TAG_ASDHinge)
    , m_connectedExternalNodes(2)
    , m_kernel(kernel)
    , m_transf(0)
    , m_doRayleigh(doRayleigh)
    , m_damping(0)
    , m_lsProbeDone(false)
    , m_initialize(1)
{
    m_connectedExternalNodes(0) = Nd1;
    m_connectedExternalNodes(1) = Nd2;
    m_nodes[0] = 0;
    m_nodes[1] = 0;
    for (int i = 0; i < 6; ++i)
        m_lsProbe[i] = 0;

    m_transf = corotational
        ? static_cast<ASDHingeTransformation*>(new ASDHingeCorotationalTransformation())
        : new ASDHingeTransformation();
    m_transf->setOrientation(R0);

    if (damping) {
        m_damping = damping->getCopy();
        if (m_damping == 0)
            opserr << "ASDHinge::ASDHinge - failed to copy the damping\n";
    }
}

ASDHinge::ASDHinge()
    : Element(0, ELE_TAG_ASDHinge)
    , m_connectedExternalNodes(2)
    , m_kernel(0)
    , m_transf(0)
    , m_doRayleigh(0)
    , m_damping(0)
    , m_lsProbeDone(false)
    , m_initialize(0)
{
    m_nodes[0] = 0;
    m_nodes[1] = 0;
    for (int i = 0; i < 6; ++i)
        m_lsProbe[i] = 0;
}

ASDHinge::~ASDHinge()
{
    this->releaseLimitStateProbes();
    if (m_kernel) delete m_kernel;
    if (m_transf) delete m_transf;
    if (m_damping) delete m_damping;
}

void ASDHinge::releaseLimitStateProbes()
{
    for (int i = 0; i < 6; ++i) {
        if (m_lsProbe[i]) {
            delete m_lsProbe[i];
            m_lsProbe[i] = 0;
        }
    }
    m_lsProbeDone = false;
}


// --------------------------------------------------------------------------
// connectivity
// --------------------------------------------------------------------------

int ASDHinge::getNumExternalNodes() const { return 2; }
const ID& ASDHinge::getExternalNodes() { return m_connectedExternalNodes; }
Node** ASDHinge::getNodePtrs() { return m_nodes; }
int ASDHinge::getNumDOF() { return 12; }

void ASDHinge::setDomain(Domain* theDomain)
{
    if (theDomain == 0) {
        m_nodes[0] = 0;
        m_nodes[1] = 0;
        return;
    }

    m_nodes[0] = theDomain->getNode(m_connectedExternalNodes(0));
    m_nodes[1] = theDomain->getNode(m_connectedExternalNodes(1));
    for (int k = 0; k < 2; ++k) {
        if (m_nodes[k] == 0) {
            opserr << "ASDHinge::setDomain - node "
                   << m_connectedExternalNodes(k) << " does not exist (element "
                   << this->getTag() << ")\n";
            exit(-1);
        }
    }

    // six dofs per node are required.  With coincident nodes and translations
    // only there is no frame to extract: every lever arm is zero, so the
    // corotational branch would be undefined.  Hinges live on beams.
    for (int k = 0; k < 2; ++k) {
        int ndf = m_nodes[k]->getNumberDOF();
        if (ndf != 6) {
            opserr << "ASDHinge::setDomain - element " << this->getTag()
                   << ": node " << m_connectedExternalNodes(k) << " has " << ndf
                   << " dofs, ASDHinge needs 6 (3D beam nodes). Use a "
                      "zeroLength for a 2D or a solid model.\n";
            exit(-1);
        }
        if (m_nodes[k]->getCrds().Size() != 3) {
            opserr << "ASDHinge::setDomain - element " << this->getTag()
                   << ": node " << m_connectedExternalNodes(k)
                   << " is not a 3D node.\n";
            exit(-1);
        }
    }

    // the two nodes must be coincident
    const Vector& X1 = m_nodes[0]->getCrds();
    const Vector& X2 = m_nodes[1]->getCrds();
    double L = 0.0, s1 = 0.0, s2 = 0.0;
    for (int i = 0; i < 3; ++i) {
        double d = X2(i) - X1(i);
        L += d * d;
        s1 += X1(i) * X1(i);
        s2 += X2(i) * X2(i);
    }
    L = sqrt(L);
    double vm = sqrt(s1 > s2 ? s1 : s2);
    if (L > ASDHINGE_LENTOL * (vm > 1.0 ? vm : 1.0)) {
        opserr << "ASDHinge::setDomain - element " << this->getTag()
               << ": the two nodes are " << L
               << " apart. ASDHinge is a zero length element.\n";
        exit(-1);
    }

    m_transf->setDomain(m_nodes, m_initialize == 0);

    if (m_damping) {
        if (m_damping->setDomain(theDomain, 6) < 0) {
            opserr << "ASDHinge::setDomain - failed to set the damping domain\n";
            exit(-1);
        }
    }

    this->DomainComponent::setDomain(theDomain);
}

int ASDHinge::setDamping(Domain* theDomain, Damping* damping)
{
    if (theDomain && damping) {
        if (m_damping) delete m_damping;
        m_damping = damping->getCopy();
        if (m_damping == 0) {
            opserr << "ASDHinge::setDamping - failed to copy the damping\n";
            return -1;
        }
        if (m_damping->setDomain(theDomain, 6) < 0) {
            opserr << "ASDHinge::setDamping - failed to set the damping domain\n";
            return -2;
        }
    }
    return 0;
}


// --------------------------------------------------------------------------
// state
// --------------------------------------------------------------------------

int ASDHinge::commitState()
{
    int retVal = this->Element::commitState();
    if (retVal < 0) {
        opserr << "ASDHinge::commitState - base class failed\n";
        return retVal;
    }
    int res = m_kernel->commitState();
    if (res < retVal) retVal = res;
    m_transf->commit();
    if (m_damping) {
        res = m_damping->commitState();
        if (res < retVal) retVal = res;
    }
    return retVal;
}

int ASDHinge::revertToLastCommit()
{
    int retVal = m_kernel->revertToLastCommit();
    m_transf->revertToLastCommit();
    if (m_damping) {
        int res = m_damping->revertToLastCommit();
        if (res < retVal) retVal = res;
    }
    return retVal;
}

int ASDHinge::revertToStart()
{
    int retVal = m_kernel->revertToStart();
    m_transf->revertToStart();
    if (m_damping) {
        int res = m_damping->revertToStart();
        if (res < retVal) retVal = res;
    }
    return retVal;
}

int ASDHinge::update()
{
    ASDHingeGlobals& g = ASDHingeGlobals::instance();

    m_transf->computeGlobalDisplacements(g.UG);
    // only here: the corotational frame is accumulated INCREMENTALLY, calling
    // update twice in a step would apply the increment twice
    m_transf->update(g.UG);

    m_transf->computeDeformations(g.UG, g.e);

    // deformation rate: the frame rotation contributes nothing at first order
    // that a uniaxial rate law can use, so the local rate is B * VG
    m_transf->computeGlobalVelocities(g.VG);
    m_transf->computeB(g.UG, g.B);
    g.edot.addMatrixVector(0.0, g.B, g.VG, 1.0);

    int res = m_kernel->setTrialDeformation(g.e, g.edot);

    if (m_damping)
        m_damping->update(m_kernel->getForce());

    return res;
}


// --------------------------------------------------------------------------
// matrices
// --------------------------------------------------------------------------

const Matrix& ASDHinge::getTangentStiff()
{
    ASDHingeGlobals& g = ASDHingeGlobals::instance();
    m_transf->computeGlobalDisplacements(g.UG);

    static Matrix kL(6, 6);
    kL = m_kernel->getTangent();
    if (m_damping)
        kL *= m_damping->getStiffnessMultiplier();

    s_K.Zero();
    s_P.Zero();
    m_transf->transformToGlobal(g.UG, m_kernel->getForce(), kL, s_K, s_P, true);
    return s_K;
}

const Matrix& ASDHinge::getInitialStiff()
{
    ASDHingeGlobals& g = ASDHingeGlobals::instance();
    // the reference configuration, corotational or not
    m_transf->computeLinearB(g.B);
    g.BtK.addMatrixTransposeProduct(0.0, g.B, m_kernel->getInitialTangent(), 1.0);
    s_K.Zero();
    s_K.addMatrixProduct(0.0, g.BtK, g.B, 1.0);
    return s_K;
}

const Matrix& ASDHinge::getDamp()
{
    if (m_doRayleigh == 1)
        return this->Element::getDamp();

    ASDHingeGlobals& g = ASDHingeGlobals::instance();
    m_transf->computeGlobalDisplacements(g.UG);
    m_transf->computeB(g.UG, g.B);
    g.BtK.addMatrixTransposeProduct(0.0, g.B, m_kernel->getDampTangent(), 1.0);
    s_K.Zero();
    s_K.addMatrixProduct(0.0, g.BtK, g.B, 1.0);
    return s_K;
}

const Matrix& ASDHinge::getMass()
{
    s_K.Zero();
    return s_K;
}


// --------------------------------------------------------------------------
// loads
// --------------------------------------------------------------------------

void ASDHinge::zeroLoad()
{
}

int ASDHinge::addLoad(ElementalLoad* theLoad, double loadFactor)
{
    (void)theLoad;
    (void)loadFactor;
    opserr << "ASDHinge::addLoad - element " << this->getTag()
           << " takes no elemental loads\n";
    return -1;
}

int ASDHinge::addInertiaLoadToUnbalance(const Vector& accel)
{
    (void)accel;
    return 0; // no mass
}


// --------------------------------------------------------------------------
// forces
// --------------------------------------------------------------------------

const Vector& ASDHinge::getResistingForce()
{
    ASDHingeGlobals& g = ASDHingeGlobals::instance();
    m_transf->computeGlobalDisplacements(g.UG);
    static Matrix dummy(6, 6);
    s_P.Zero();
    m_transf->transformToGlobal(g.UG, m_kernel->getForce(), dummy, s_K, s_P, false);
    if (m_damping)
        this->assembleDampingForce(s_P);
    return s_P;
}

void ASDHinge::assembleDampingForce(Vector& P)
{
    ASDHingeGlobals& g = ASDHingeGlobals::instance();
    m_transf->computeB(g.UG, g.B);
    P.addMatrixTransposeVector(1.0, g.B, m_damping->getDampingForce(), 1.0);
}

const Vector& ASDHinge::getResistingForceIncInertia()
{
    // no mass; only the Rayleigh damping forces, when asked for
    const Vector& P = this->getResistingForce();
    if (m_doRayleigh == 1 &&
        (alphaM != 0.0 || betaK != 0.0 || betaK0 != 0.0 || betaKc != 0.0))
        s_P.addVector(1.0, this->getRayleighDampingForces(), 1.0);
    return P;
}


// --------------------------------------------------------------------------
// staged construction
// --------------------------------------------------------------------------

void ASDHinge::onActivate()
{
    // The baseline moves to the current configuration, so a hinge born inside
    // a construction stage has zero deformation.  For the corotational
    // transformation this also re-bases the nodal triads to the identity: it
    // is not enough to reset the translational offset, or the hinge is born
    // with a finite relative rotation.
    m_transf->forceCaptureInitialDisp();
    this->update();
}

void ASDHinge::onDeactivate()
{
}


// --------------------------------------------------------------------------
// serialization
// --------------------------------------------------------------------------

int ASDHinge::sendSelf(int commitTag, Channel& theChannel)
{
    int dataTag = this->getDbTag();

    static ID idata(9);
    idata(0) = this->getTag();
    idata(1) = m_connectedExternalNodes(0);
    idata(2) = m_connectedExternalNodes(1);
    idata(3) = m_transf->isLinear() ? 0 : 1;
    idata(4) = m_kernel ? m_kernel->getTypeId() : -1;
    idata(5) = m_doRayleigh;
    idata(6) = is_this_element_active ? 1 : 0;
    if (m_damping) {
        idata(7) = m_damping->getClassTag();
        int dbTag = m_damping->getDbTag();
        if (dbTag == 0) {
            dbTag = theChannel.getDbTag();
            if (dbTag != 0)
                m_damping->setDbTag(dbTag);
        }
        idata(8) = dbTag;
    }
    else {
        idata(7) = 0;
        idata(8) = 0;
    }
    if (theChannel.sendID(dataTag, commitTag, idata) < 0) {
        opserr << "ASDHinge::sendSelf - failed to send the ID\n";
        return -1;
    }

    // ONE Vector: the orientation, then whatever the transformation owns.
    // Keeping it to a single message means its size cannot be confused with
    // another one (the trap ZeroLength::sendSelf documents).
    int nt = m_transf->internalDataSize();
    Vector vdata(9 + nt);
    double axis[3];
    for (int j = 0; j < 3; ++j) {
        m_transf->getLocalAxis(j, axis);
        for (int i = 0; i < 3; ++i)
            vdata(3 * j + i) = axis[i];
    }
    m_transf->saveInternalData(vdata, 9);
    if (theChannel.sendVector(dataTag, commitTag, vdata) < 0) {
        opserr << "ASDHinge::sendSelf - failed to send the Vector\n";
        return -1;
    }

    if (m_kernel && m_kernel->sendSelf(commitTag, theChannel) < 0) {
        opserr << "ASDHinge::sendSelf - the kernel failed\n";
        return -1;
    }

    if (m_damping && m_damping->sendSelf(commitTag, theChannel) < 0) {
        opserr << "ASDHinge::sendSelf - the damping failed\n";
        return -1;
    }

    return 0;
}

int ASDHinge::recvSelf(int commitTag, Channel& theChannel,
                       FEM_ObjectBroker& theBroker)
{
    int dataTag = this->getDbTag();

    static ID idata(9);
    if (theChannel.recvID(dataTag, commitTag, idata) < 0) {
        opserr << "ASDHinge::recvSelf - failed to recv the ID\n";
        return -1;
    }
    this->setTag(idata(0));
    m_connectedExternalNodes(0) = idata(1);
    m_connectedExternalNodes(1) = idata(2);
    bool corotational = (idata(3) == 1);
    int kernelType = idata(4);
    m_doRayleigh = idata(5);

    if (m_transf) delete m_transf;
    m_transf = corotational
        ? static_cast<ASDHingeTransformation*>(new ASDHingeCorotationalTransformation())
        : new ASDHingeTransformation();

    int nt = m_transf->internalDataSize();
    Vector vdata(9 + nt);
    if (theChannel.recvVector(dataTag, commitTag, vdata) < 0) {
        opserr << "ASDHinge::recvSelf - failed to recv the Vector\n";
        return -1;
    }
    double R0[3][3];
    for (int j = 0; j < 3; ++j)
        for (int i = 0; i < 3; ++i)
            R0[i][j] = vdata(3 * j + i);
    m_transf->setOrientation(R0);
    m_transf->restoreInternalData(vdata, 9);

    this->releaseLimitStateProbes();
    if (m_kernel) delete m_kernel;
    m_kernel = ASDHingeKernel::create(kernelType);
    if (m_kernel == 0)
        return -1;
    if (m_kernel->recvSelf(commitTag, theChannel, theBroker) < 0) {
        opserr << "ASDHinge::recvSelf - the kernel failed\n";
        return -1;
    }

    if (m_damping) {
        delete m_damping;
        m_damping = 0;
    }
    if (idata(7) != 0) {
        m_damping = theBroker.getNewDamping(idata(7));
        if (m_damping == 0) {
            opserr << "ASDHinge::recvSelf - no damping of class tag "
                   << idata(7) << "\n";
            return -1;
        }
        m_damping->setDbTag(idata(8));
        if (m_damping->recvSelf(commitTag, theChannel, theBroker) < 0) {
            opserr << "ASDHinge::recvSelf - the damping failed\n";
            return -1;
        }
    }

    // an element deactivated before the transfer must come back deactivated
    is_this_element_active = (idata(6) == 1);

    // the baseline came in with the transformation: setDomain must not
    // overwrite it
    m_initialize = 1;

    return 0;
}


// --------------------------------------------------------------------------
// print
// --------------------------------------------------------------------------

void ASDHinge::Print(OPS_Stream& s, int flag)
{
    (void)flag;
    s << "Element: " << this->getTag() << " type: ASDHinge  iNode: "
      << m_connectedExternalNodes(0) << "  jNode: "
      << m_connectedExternalNodes(1) << "\n";
    s << "  frame: " << (m_transf->isLinear() ? "linear" : "corotational")
      << "\n";
    double axis[3];
    for (int j = 0; j < 3; ++j) {
        m_transf->getLocalAxis(j, axis);
        s << "  local " << (char)('x' + j) << ": (" << axis[0] << ", "
          << axis[1] << ", " << axis[2] << ")\n";
    }
    s << "  slots:\n";
    if (m_kernel)
        m_kernel->Print(s);
}


// --------------------------------------------------------------------------
// responses
// --------------------------------------------------------------------------

Response* ASDHinge::setResponse(const char** argv, int argc, OPS_Stream& output)
{
    if (argc < 1)
        return 0;

    Response* theResponse = 0;

    output.tag("ElementOutput");
    output.attr("eleType", "ASDHinge");
    output.attr("eleTag", this->getTag());
    output.attr("node1", m_connectedExternalNodes(0));
    output.attr("node2", m_connectedExternalNodes(1));

    // -- the documented names: Hinge.Force, Hinge.Deformation, ... -----------
    // MPCO splits an -E key on the dot, so -E "Hinge.Force" arrives here as
    // {"Hinge", "Force"} and goes back to STKO under the same name.
    const char* sub = 0;
    if (argc > 1 && (strcmp(argv[0], "Hinge") == 0 || strcmp(argv[0], "hinge") == 0))
        sub = argv[1];

    if (sub && (strcmp(sub, "Force") == 0 || strcmp(sub, "force") == 0)) {
        // The Hinge.* responses live under a GaussPointOutput tag: the MPCO
        // recorder derives the result location from the tag stack, and STKO
        // lists a result in the gauss-point plots only if its components sit
        // under a gauss point.  This element is registered in the recorder as
        // Line_2N with 1 gauss point, so "number" is always 1.
        output.tag("GaussPointOutput");
        output.attr("number", 1);
        for (int i = 0; i < 6; ++i)
            output.tag("ResponseType", FORCE_NAMES[i]);
        theResponse = new ElementResponse(this, 2, Vector(6));
        output.endTag();
    }
    else if (sub && (strcmp(sub, "Deformation") == 0 ||
                     strcmp(sub, "deformation") == 0)) {
        output.tag("GaussPointOutput");
        output.attr("number", 1);
        for (int i = 0; i < 6; ++i)
            output.tag("ResponseType", DOF_NAMES[i]);
        theResponse = new ElementResponse(this, 3, Vector(6));
        output.endTag();
    }
    else if (sub && (strcmp(sub, "Stiffness") == 0 ||
                     strcmp(sub, "stiffness") == 0)) {
        theResponse = new ElementResponse(this, 13, Matrix(6, 6));
    }
    else if (sub && (strcmp(sub, "LimitStateRatio") == 0 ||
                     strcmp(sub, "limitStateRatio") == 0 ||
                     strcmp(sub, "LS") == 0)) {
        // Only exists if at least one slot provides it.  ASDHysteretic1D
        // refuses the response when it has no -limitStates, and it is right to:
        // a permanent zero in a recorder looks like a model that never
        // yielded.  Slots that do not provide it read 0, and that is stated in
        // the documentation and printed by Print.
        if (this->buildLimitStateProbes()) {
            output.tag("GaussPointOutput");
            output.attr("number", 1);
            for (int i = 0; i < 6; ++i)
                output.tag("ResponseType", DOF_NAMES[i]);
            theResponse = new ElementResponse(this, 7, Vector(6));
            output.endTag();
        }
        else {
            opserr << "ASDHinge::setResponse - element " << this->getTag()
                   << ": no slot defines -limitStates, "
                      "Hinge.LimitStateRatio is not available\n";
        }
    }
    else if (sub && (strcmp(sub, "SlotState") == 0 ||
                     strcmp(sub, "slotState") == 0 ||
                     strcmp(sub, "SlotStates") == 0)) {
        // What each slot IS - Free (0), Linear (1) or Material (2) - as a
        // gauss-point result, so a post-processor can tell an elastic spring
        // from a slot that has no spring at all.  Hinge.LimitStateRatio cannot:
        // it reads 0 in both cases.  The classic "slotStates" below returns the
        // same numbers as an ID and stays for the text recorders; this one is
        // tagged like the other Hinge.* results, so it reaches STKO named after
        // the six local dofs instead of C1..C6.
        output.tag("GaussPointOutput");
        output.attr("number", 1);
        for (int i = 0; i < 6; ++i)
            output.tag("ResponseType", DOF_NAMES[i]);
        theResponse = new ElementResponse(this, 8, Vector(6));
        output.endTag();
    }

    // -- the classic element names ------------------------------------------
    else if (strcmp(argv[0], "force") == 0 || strcmp(argv[0], "forces") == 0 ||
             strcmp(argv[0], "globalForce") == 0 ||
             strcmp(argv[0], "globalForces") == 0) {
        char buf[16];
        for (int n = 0; n < 2; ++n)
            for (int i = 0; i < 6; ++i) {
                sprintf(buf, "P%d_%d", n + 1, i + 1);
                output.tag("ResponseType", buf);
            }
        theResponse = new ElementResponse(this, 1, Vector(12));
    }
    else if (strcmp(argv[0], "localForce") == 0 ||
             strcmp(argv[0], "localForces") == 0 ||
             strcmp(argv[0], "basicForce") == 0 ||
             strcmp(argv[0], "basicForces") == 0) {
        for (int i = 0; i < 6; ++i)
            output.tag("ResponseType", FORCE_NAMES[i]);
        theResponse = new ElementResponse(this, 2, Vector(6));
    }
    else if (strcmp(argv[0], "deformation") == 0 ||
             strcmp(argv[0], "deformations") == 0 ||
             strcmp(argv[0], "defo") == 0 ||
             strcmp(argv[0], "basicDeformation") == 0) {
        for (int i = 0; i < 6; ++i)
            output.tag("ResponseType", DOF_NAMES[i]);
        theResponse = new ElementResponse(this, 3, Vector(6));
    }
    else if (strcmp(argv[0], "basicStiffness") == 0) {
        theResponse = new ElementResponse(this, 13, Matrix(6, 6));
    }
    else if (strcmp(argv[0], "dampingForces") == 0 ||
             strcmp(argv[0], "rayleighForces") == 0) {
        theResponse = new ElementResponse(this, 15, Vector(12));
    }
    else if (strcmp(argv[0], "xaxis") == 0) {
        theResponse = new ElementResponse(this, 20, Vector(3));
    }
    else if (strcmp(argv[0], "yaxis") == 0) {
        theResponse = new ElementResponse(this, 21, Vector(3));
    }
    else if (strcmp(argv[0], "zaxis") == 0) {
        theResponse = new ElementResponse(this, 22, Vector(3));
    }
    else if (strcmp(argv[0], "materials") == 0) {
        theResponse = new ElementResponse(this, 23, ID(6));
    }
    else if (strcmp(argv[0], "slotStates") == 0) {
        theResponse = new ElementResponse(this, 24, ID(6));
    }

    // -- forwarding to one slot ---------------------------------------------
    // The index is the LOCAL DOF, 1..6, always.  In a zeroLength it is the
    // position in a -mat list whose length varies from hinge to hinge, which
    // is exactly the ambiguity this element removes.
    else if (strcmp(argv[0], "material") == 0 && argc > 2) {
        int slot = atoi(argv[1]) - 1;
        if (slot >= 0 && slot < 6) {
            UniaxialMaterial* mat = m_kernel ? m_kernel->getMaterial(slot) : 0;
            if (mat)
                theResponse = mat->setResponse(&argv[2], argc - 2, output);
        }
    }

    output.endTag();
    return theResponse;
}

bool ASDHinge::buildLimitStateProbes()
{
    if (m_lsProbeDone) {
        for (int i = 0; i < 6; ++i)
            if (m_lsProbe[i]) return true;
        return false;
    }
    m_lsProbeDone = true;
    bool any = false;
    const char* arg = "limitStateRatio";
    for (int i = 0; i < 6; ++i) {
        UniaxialMaterial* mat = m_kernel ? m_kernel->getMaterial(i) : 0;
        if (mat == 0) continue;
        // a throwaway stream: the material tags its own output, and they must
        // not end up in the element descriptor
        DummyStream dummy;
        m_lsProbe[i] = mat->setResponse(&arg, 1, dummy);
        if (m_lsProbe[i]) any = true;
    }
    return any;
}

int ASDHinge::getResponse(int responseID, Information& eleInfo)
{
    ASDHingeGlobals& g = ASDHingeGlobals::instance();
    static Vector v6(6);
    static Vector v3(3);
    static Vector v12(12);
    static ID id6(6);

    switch (responseID) {

    case 1:
        return eleInfo.setVector(this->getResistingForce());

    case 2:
        return eleInfo.setVector(m_kernel->getForce());

    case 3:
        m_transf->computeGlobalDisplacements(g.UG);
        m_transf->computeDeformations(g.UG, g.e);
        return eleInfo.setVector(g.e);

    case 7: {
        v6.Zero();
        for (int i = 0; i < 6; ++i) {
            if (m_lsProbe[i] == 0) continue;
            if (m_lsProbe[i]->getResponse() < 0) continue;
            Information& mi = m_lsProbe[i]->getInformation();
            if (mi.theVector != 0 && mi.theVector->Size() > 0)
                v6(i) = (*mi.theVector)(0);
            else
                v6(i) = mi.theDouble;
        }
        return eleInfo.setVector(v6);
    }

    case 8: {
        for (int i = 0; i < 6; ++i)
            v6(i) = m_kernel ? static_cast<double>(m_kernel->getSlotState(i)) : 0.0;
        return eleInfo.setVector(v6);
    }

    case 13:
        return eleInfo.setMatrix(m_kernel->getTangent());

    case 15: {
        v12.Zero();
        if (m_doRayleigh == 1)
            v12 = this->getRayleighDampingForces();
        if (m_damping) {
            m_transf->computeGlobalDisplacements(g.UG);
            this->assembleDampingForce(v12);
        }
        return eleInfo.setVector(v12);
    }

    case 20:
    case 21:
    case 22: {
        double axis[3];
        m_transf->getLocalAxis(responseID - 20, axis);
        for (int i = 0; i < 3; ++i) v3(i) = axis[i];
        return eleInfo.setVector(v3);
    }

    case 23: {
        for (int i = 0; i < 6; ++i) {
            UniaxialMaterial* mat = m_kernel ? m_kernel->getMaterial(i) : 0;
            id6(i) = mat ? mat->getTag() : 0;
        }
        return eleInfo.setID(id6);
    }

    case 24: {
        for (int i = 0; i < 6; ++i)
            id6(i) = m_kernel ? m_kernel->getSlotState(i) : 0;
        return eleInfo.setID(id6);
    }

    default:
        return -1;
    }
}

int ASDHinge::setParameter(const char** argv, int argc, Parameter& param)
{
    if (argc < 1)
        return -1;

    // material $i ... , with i the LOCAL DOF 1..6
    if ((strcmp(argv[0], "material") == 0 || strcmp(argv[0], "-material") == 0)
        && argc > 2) {
        int slot = atoi(argv[1]) - 1;
        if (slot < 0 || slot > 5)
            return -1;
        UniaxialMaterial* mat = m_kernel ? m_kernel->getMaterial(slot) : 0;
        if (mat == 0)
            return -1;
        return mat->setParameter(&argv[2], argc - 2, param);
    }

    // otherwise offer it to every slot
    int result = -1;
    for (int i = 0; i < 6; ++i) {
        UniaxialMaterial* mat = m_kernel ? m_kernel->getMaterial(i) : 0;
        if (mat == 0) continue;
        int res = mat->setParameter(argv, argc, param);
        if (res != -1)
            result = res;
    }
    return result;
}
