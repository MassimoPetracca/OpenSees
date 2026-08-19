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

// Original implementation: Massimo Petracca (ASDEA)
//
// See ASDTet.h for the design notes.

#include <ASDTet.h>
#include <ASDSolidTet4CorotationalTransformation.h>
#include <NDMaterial.h>
#include <Domain.h>
#include <Node.h>
#include <ElementalLoad.h>
#include <Channel.h>
#include <FEM_ObjectBroker.h>
#include <MaterialResponse.h>
#include <ElementResponse.h>
#include <Information.h>
#include <Parameter.h>
#include <elementAPI.h>
#include <cmath>
#include <cstring>

namespace {

    constexpr int OPT_UPDATE = 1;
    constexpr int OPT_LHS = 2;
    constexpr int OPT_RHS = 4;
    constexpr int OPT_LHS_IS_INITIAL = 8;

}

void*
OPS_ASDTet(void)
{
    static bool first_done = false;
    if (!first_done) {
        opserr << "Using ASDTet - Developed by: Massimo Petracca, Guido Camata, ASDEA Software Technology\n";
        first_done = true;
    }

    int numArgs = OPS_GetNumRemainingInputArgs();
    if (numArgs < 6) {
        opserr << "Want: element ASDTet $tag $Node1 $Node2 $Node3 $Node4 $matTag "
            "<-corotational> <-b $bx $by $bz>\n";
        return 0;
    }

    int iData[6];
    int numData = 6;
    if (OPS_GetInt(&numData, iData) != 0) {
        opserr << "WARNING invalid integer input: element ASDTet\n";
        return 0;
    }

    NDMaterial* mat = OPS_getNDMaterial(iData[5]);
    if (mat == nullptr) {
        opserr << "ERROR: element ASDTet " << iData[0]
            << " NDMaterial " << iData[5] << " not found\n";
        return 0;
    }

    bool corotational = false;
    double body[3] = { 0.0, 0.0, 0.0 };
    while (OPS_GetNumRemainingInputArgs() > 0) {
        const char* type = OPS_GetString();
        if (strcmp(type, "-corotational") == 0) {
            corotational = true;
        }
        else if (strcmp(type, "-b") == 0) {
            if (OPS_GetNumRemainingInputArgs() < 3) {
                opserr << "Error: element ASDTet: -b needs 3 components\n";
                return 0;
            }
            int nd = 3;
            if (OPS_GetDoubleInput(&nd, body) < 0) {
                opserr << "Error: element ASDTet: invalid -b components\n";
                return 0;
            }
        }
        else {
            opserr << "Error: element ASDTet: unknown option '" << type << "'\n";
            return 0;
        }
    }

    return new ASDTet(iData[0], iData[1], iData[2], iData[3], iData[4],
        mat, corotational, body);
}

ASDTet::ASDTet()
    : Element(0, ELE_TAG_ASDTet)
    , m_node_ids(4)
    , m_g(4, 3)
    , m_B(6, 12)
    , m_P0(12)
    , m_U0(12)
{
}

ASDTet::ASDTet(int tag, int node1, int node2, int node3, int node4,
    NDMaterial* material, bool corotational, const double* body)
    : Element(tag, ELE_TAG_ASDTet)
    , m_node_ids(4)
    , m_use_corotational(corotational)
    , m_g(4, 3)
    , m_B(6, 12)
    , m_P0(12)
    , m_U0(12)
{
    m_node_ids(0) = node1;
    m_node_ids(1) = node2;
    m_node_ids(2) = node3;
    m_node_ids(3) = node4;
    m_material = material->getCopy("ThreeDimensional");
    if (m_material == nullptr) {
        opserr << "ASDTet ERROR: element " << tag << " - material " << material->getTag()
            << " cannot return a ThreeDimensional copy.\n";
        exit(-1);
    }
    for (int i = 0; i < 3; ++i)
        m_body[i] = body ? body[i] : 0.0;
    if (m_use_corotational)
        m_transformation = new ASDSolidTet4CorotationalTransformation();
}

ASDTet::~ASDTet()
{
    if (m_material)
        delete m_material;
    if (m_transformation)
        delete m_transformation;
}

void ASDTet::setDomain(Domain* theDomain)
{
    if (theDomain == nullptr) {
        for (int i = 0; i < 4; ++i)
            m_nodes[static_cast<std::size_t>(i)] = nullptr;
        DomainComponent::setDomain(theDomain);
        return;
    }

    for (int i = 0; i < 4; ++i) {
        Node* node = theDomain->getNode(m_node_ids(i));
        if (node == nullptr) {
            opserr << "ASDTet ERROR in setDomain: node " << m_node_ids(i)
                << " does not exist in the domain\n";
            exit(-1);
        }
        if (node->getCrds().Size() != 3 || node->getNumberDOF() != 3) {
            opserr << "ASDTet ERROR in setDomain: node " << m_node_ids(i)
                << " should have 3 coordinates and 3 DOFs\n";
            exit(-1);
        }
        m_nodes[static_cast<std::size_t>(i)] = node;
    }

    // constant gradients and volume from the reference geometry
    static Matrix J(3, 3);
    for (int j = 0; j < 3; ++j) {
        const Vector& X0 = m_nodes[0]->getCrds();
        const Vector& Xj = m_nodes[static_cast<std::size_t>(j + 1)]->getCrds();
        for (int i = 0; i < 3; ++i)
            J(i, j) = Xj(i) - X0(i);
    }
    double det =
        J(0, 0) * (J(1, 1) * J(2, 2) - J(1, 2) * J(2, 1)) -
        J(0, 1) * (J(1, 0) * J(2, 2) - J(1, 2) * J(2, 0)) +
        J(0, 2) * (J(1, 0) * J(2, 1) - J(1, 1) * J(2, 0));
    m_V = det / 6.0;
    if (m_V <= 0.0) {
        opserr << "ASDTet ERROR in setDomain: element " << getTag()
            << " has non-positive volume " << m_V
            << " (check the node ordering: 1-2-3 counterclockwise seen from 4)\n";
        exit(-1);
    }
    static Matrix invJ(3, 3);
    J.Invert(invJ);
    // rows of invJ are dN_a/dX for a = 1..3; node 0 closes the partition
    for (int j = 0; j < 3; ++j) {
        double s = 0.0;
        for (int a = 0; a < 3; ++a) {
            m_g(a + 1, j) = invJ(a, j);
            s += invJ(a, j);
        }
        m_g(0, j) = -s;
    }
    // B (6x12), Voigt [exx eyy ezz gxy gyz gzx], engineering shear
    m_B.Zero();
    for (int a = 0; a < 4; ++a) {
        int c = 3 * a;
        double gx = m_g(a, 0), gy = m_g(a, 1), gz = m_g(a, 2);
        m_B(0, c) = gx;
        m_B(1, c + 1) = gy;
        m_B(2, c + 2) = gz;
        m_B(3, c) = gy; m_B(3, c + 1) = gx;
        m_B(4, c + 1) = gz; m_B(4, c + 2) = gy;
        m_B(5, c) = gz; m_B(5, c + 2) = gx;
    }

    if (m_transformation)
        m_transformation->setReference(m_nodes, m_g);

    // initial displacements (activation-time), for staged construction
    if (!m_initialized) {
        for (int a = 0; a < 4; ++a) {
            const Vector& iu = m_nodes[static_cast<std::size_t>(a)]->getTrialDisp();
            for (int i = 0; i < 3; ++i)
                m_U0(3 * a + i) = iu(i);
        }
        m_initialized = true;
    }

    DomainComponent::setDomain(theDomain);
}

void ASDTet::Print(OPS_Stream& s, int flag)
{
    if (flag == OPS_PRINT_PRINTMODEL_JSON) {
        s << "\t\t\t{";
        s << "\"name\": " << this->getTag() << ", ";
        s << "\"type\": \"ASDTet\", ";
        s << "\"nodes\": [" << m_node_ids(0) << ", " << m_node_ids(1) << ", "
            << m_node_ids(2) << ", " << m_node_ids(3) << "], ";
        s << "\"material\": \"" << m_material->getTag() << "\"}";
        return;
    }
    s << "ASDTet tag: " << this->getTag() << endln;
    s << "   nodes: " << m_node_ids(0) << " " << m_node_ids(1) << " "
        << m_node_ids(2) << " " << m_node_ids(3) << endln;
    s << "   material: " << m_material->getTag() << ", volume: " << m_V << endln;
    s << "   kinematics: " << (m_use_corotational ? "corotational" : "linear") << endln;
}

int ASDTet::getNumExternalNodes() const
{
    return 4;
}

const ID& ASDTet::getExternalNodes()
{
    return m_node_ids;
}

Node** ASDTet::getNodePtrs()
{
    return m_nodes.data();
}

int ASDTet::getNumDOF()
{
    return 12;
}

void ASDTet::computeDisplacements(Vector& d) const
{
    for (int a = 0; a < 4; ++a) {
        const Vector& iu = m_nodes[static_cast<std::size_t>(a)]->getTrialDisp();
        for (int i = 0; i < 3; ++i)
            d(3 * a + i) = iu(i) - m_U0(3 * a + i);
    }
}

int ASDTet::calculateAll(Matrix& K, Vector& f, int options)
{
    int result = 0;

    // deformational displacements
    static Vector d(12);
    static Vector UL(12);
    computeDisplacements(d);
    if (m_transformation)
        m_transformation->calculateLocalDisplacements(d, UL);
    else
        UL = d;

    // strain update at the single integration point
    if (options & OPT_UPDATE) {
        static Vector eps(6);
        eps.addMatrixVector(0.0, m_B, UL, 1.0);
        result += m_material->setTrialStrain(eps);
    }

    // local stiffness: K_L = V * B^T D B
    bool lhs = (options & OPT_LHS);
    if (lhs) {
        const Matrix& D = (options & OPT_LHS_IS_INITIAL) ?
            m_material->getInitialTangent() : m_material->getTangent();
        K.addMatrixTripleProduct(0.0, m_B, D, m_V);
    }

    // local internal force: f_L = V * B^T sigma
    if (options & OPT_RHS) {
        const Vector& sig = m_material->getStress();
        f.addMatrixTransposeVector(0.0, m_B, sig, m_V);
    }

    // corotational transformation to global
    if (m_transformation && (options & (OPT_LHS | OPT_RHS))) {
        // the transformation needs the RHS even when only the LHS is
        // requested (the geometric terms are built from the local forces)
        if (!(options & OPT_RHS)) {
            const Vector& sig = m_material->getStress();
            f.addMatrixTransposeVector(0.0, m_B, sig, m_V);
        }
        m_transformation->transformToGlobal(d, K, f, lhs);
    }

    return result;
}

int ASDTet::update()
{
    static Matrix K(12, 12);
    static Vector f(12);
    return calculateAll(K, f, OPT_UPDATE);
}

int ASDTet::commitState()
{
    return m_material->commitState();
}

int ASDTet::revertToLastCommit()
{
    return m_material->revertToLastCommit();
}

int ASDTet::revertToStart()
{
    return m_material->revertToStart();
}

const Matrix& ASDTet::getTangentStiff()
{
    static Matrix K(12, 12);
    static Vector f(12);
    calculateAll(K, f, OPT_LHS | OPT_RHS);
    return K;
}

const Matrix& ASDTet::getInitialStiff()
{
    static Matrix K(12, 12);
    static Vector f(12);
    calculateAll(K, f, OPT_LHS | OPT_LHS_IS_INITIAL);
    return K;
}

const Matrix& ASDTet::getMass()
{
    // lumped mass: rho V / 4 per node
    static Matrix M(12, 12);
    M.Zero();
    double rho = m_material->getRho();
    if (rho > 0.0) {
        double m = rho * m_V * 0.25;
        for (int i = 0; i < 12; ++i)
            M(i, i) = m;
    }
    return M;
}

void ASDTet::zeroLoad()
{
    m_P0.Zero();
}

int ASDTet::addLoad(ElementalLoad* theLoad, double loadFactor)
{
    int type;
    const Vector& data = theLoad->getData(type, loadFactor);

    if (type == LOAD_TAG_BrickSelfWeight) {
        // uses the element -b components (body force per unit mass)
        double rho = m_material->getRho();
        double m = rho * m_V * 0.25 * loadFactor;
        for (int a = 0; a < 4; ++a)
            for (int i = 0; i < 3; ++i)
                m_P0(3 * a + i) += m * m_body[i];
        return 0;
    }
    else if (type == LOAD_TAG_SelfWeight) {
        // data = (xf, yf, zf): gravity direction factors
        double rho = m_material->getRho();
        double m = rho * m_V * 0.25 * loadFactor;
        for (int a = 0; a < 4; ++a)
            for (int i = 0; i < 3; ++i)
                m_P0(3 * a + i) += m * data(i);
        return 0;
    }

    opserr << "ASDTet::addLoad() - element " << getTag()
        << " - load type " << type << " unknown\n";
    return -1;
}

int ASDTet::addInertiaLoadToUnbalance(const Vector& accel)
{
    double rho = m_material->getRho();
    if (rho == 0.0)
        return 0;
    double m = rho * m_V * 0.25;
    for (int a = 0; a < 4; ++a) {
        const Vector& Raccel = m_nodes[static_cast<std::size_t>(a)]->getRV(accel);
        for (int i = 0; i < 3; ++i)
            m_P0(3 * a + i) -= m * Raccel(i);
    }
    return 0;
}

const Vector& ASDTet::getResistingForce()
{
    static Matrix K(12, 12);
    static Vector f(12);
    calculateAll(K, f, OPT_RHS);
    f.addVector(1.0, m_P0, -1.0);
    return f;
}

const Vector& ASDTet::getResistingForceIncInertia()
{
    static Vector f(12);
    f = getResistingForce();
    double rho = m_material->getRho();
    if (rho > 0.0) {
        double m = rho * m_V * 0.25;
        for (int a = 0; a < 4; ++a) {
            const Vector& accel = m_nodes[static_cast<std::size_t>(a)]->getTrialAccel();
            for (int i = 0; i < 3; ++i)
                f(3 * a + i) += m * accel(i);
        }
    }
    return f;
}

int ASDTet::sendSelf(int commitTag, Channel& theChannel)
{
    int res = 0;
    int dataTag = getDbTag();

    // INT data: tag, nodes, flags, material (classTag, dbTag)
    static ID idData(10);
    idData(0) = getTag();
    for (int i = 0; i < 4; ++i)
        idData(1 + i) = m_node_ids(i);
    idData(5) = m_use_corotational ? 1 : 0;
    idData(6) = m_initialized ? 1 : 0;
    idData(7) = m_material->getClassTag();
    int matDbTag = m_material->getDbTag();
    if (matDbTag == 0) {
        matDbTag = theChannel.getDbTag();
        if (matDbTag != 0)
            m_material->setDbTag(matDbTag);
    }
    idData(8) = matDbTag;
    idData(9) = 0; // reserved
    res = theChannel.sendID(dataTag, commitTag, idData);
    if (res < 0) {
        opserr << "WARNING ASDTet::sendSelf() - " << getTag() << " failed to send ID\n";
        return res;
    }

    // DOUBLE data: body force, initial displacements
    static Vector vectData(15);
    for (int i = 0; i < 3; ++i)
        vectData(i) = m_body[i];
    for (int i = 0; i < 12; ++i)
        vectData(3 + i) = m_U0(i);
    res = theChannel.sendVector(dataTag, commitTag, vectData);
    if (res < 0) {
        opserr << "WARNING ASDTet::sendSelf() - " << getTag() << " failed to send Vector\n";
        return res;
    }

    // the material
    res = m_material->sendSelf(commitTag, theChannel);
    if (res < 0) {
        opserr << "WARNING ASDTet::sendSelf() - " << getTag() << " failed to send its material\n";
        return res;
    }

    return res;
}

int ASDTet::recvSelf(int commitTag, Channel& theChannel, FEM_ObjectBroker& theBroker)
{
    int res = 0;
    int dataTag = getDbTag();

    static ID idData(10);
    res = theChannel.recvID(dataTag, commitTag, idData);
    if (res < 0) {
        opserr << "WARNING ASDTet::recvSelf() - failed to receive ID\n";
        return res;
    }
    setTag(idData(0));
    for (int i = 0; i < 4; ++i)
        m_node_ids(i) = idData(1 + i);
    m_use_corotational = idData(5) == 1;
    m_initialized = idData(6) == 1;
    int matClassTag = idData(7);
    int matDbTag = idData(8);

    if (m_use_corotational && m_transformation == nullptr)
        m_transformation = new ASDSolidTet4CorotationalTransformation();
    else if (!m_use_corotational && m_transformation) {
        delete m_transformation;
        m_transformation = nullptr;
    }

    static Vector vectData(15);
    res = theChannel.recvVector(dataTag, commitTag, vectData);
    if (res < 0) {
        opserr << "WARNING ASDTet::recvSelf() - failed to receive Vector\n";
        return res;
    }
    for (int i = 0; i < 3; ++i)
        m_body[i] = vectData(i);
    for (int i = 0; i < 12; ++i)
        m_U0(i) = vectData(3 + i);

    if (m_material == nullptr || m_material->getClassTag() != matClassTag) {
        if (m_material)
            delete m_material;
        m_material = theBroker.getNewNDMaterial(matClassTag);
        if (m_material == nullptr) {
            opserr << "WARNING ASDTet::recvSelf() - could not get an NDMaterial with classTag "
                << matClassTag << "\n";
            return -1;
        }
    }
    m_material->setDbTag(matDbTag);
    res = m_material->recvSelf(commitTag, theChannel, theBroker);
    if (res < 0) {
        opserr << "WARNING ASDTet::recvSelf() - the material failed to recvSelf\n";
        return res;
    }

    return res;
}

Response* ASDTet::setResponse(const char** argv, int argc, OPS_Stream& output)
{
    Response* theResponse = nullptr;

    output.tag("ElementOutput");
    output.attr("eleType", "ASDTet");
    output.attr("eleTag", getTag());
    for (int i = 0; i < 4; ++i) {
        char nodeTag[32];
        snprintf(nodeTag, 32, "node%d", i + 1);
        output.attr(nodeTag, m_node_ids(i));
    }

    if (strcmp(argv[0], "force") == 0 || strcmp(argv[0], "forces") == 0 ||
        strcmp(argv[0], "globalForce") == 0 || strcmp(argv[0], "globalForces") == 0) {
        theResponse = new ElementResponse(this, 1, Vector(12));
    }
    else if (strcmp(argv[0], "material") == 0 && argc > 1) {
        // accept both "material <resp> ..." and "material 1 <resp> ..."
        // (the single integration point)
        if (argc > 2 && strcmp(argv[1], "1") == 0)
            theResponse = m_material->setResponse(&argv[2], argc - 2, output);
        else
            theResponse = m_material->setResponse(&argv[1], argc - 1, output);
    }
    else if (strcmp(argv[0], "stress") == 0 || strcmp(argv[0], "stresses") == 0) {
        theResponse = m_material->setResponse(argv, argc, output);
    }
    else if (strcmp(argv[0], "strain") == 0 || strcmp(argv[0], "strains") == 0) {
        theResponse = m_material->setResponse(argv, argc, output);
    }

    output.endTag();
    return theResponse;
}

int ASDTet::getResponse(int responseID, Information& eleInfo)
{
    if (responseID == 1)
        return eleInfo.setVector(getResistingForce());
    return -1;
}
