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
                                                                        
// $Revision: 1.10 $
// $Date: 2021/04/28 22:51:21 $

// Original implementation: Massimo Petracca (ASDEA)
//
//

#include <ASDConstraintEquationElement.h>

#include <Domain.h>
#include <Node.h>
#include <ErrorHandler.h>
#include <Channel.h>
#include <FEM_ObjectBroker.h>
#include <elementAPI.h>
#include <Renderer.h>
#include <analysis/dof_grp/DOF_Group.h>

#include <stdio.h>
#include <stdlib.h>
#include <cmath>
#include <string>
#include <map>

void *
OPS_ASDConstraintEquationElement(void)
{
    static bool first_done = false;
    if (!first_done) {
        opserr << "Using ASDConstraintEquationElement - Developed by: Massimo Petracca, Guido Camata, ASDEA Software Technology\n";
        first_done = true;
    }

    const char* descr = "Want: element ASDConstraintEquationElement $tag $K $cNode $cDof <$rNode1 $rDof1 $rFact1 ... $rNodeN $rDofN $rFactN>\n";

    int numArgs = OPS_GetNumRemainingInputArgs();
    if (numArgs < 4) {
        opserr << "ASDConstraintEquationElement ERROR : Few arguments:\n" << descr;
        return 0;
    }
    
    // mandatory parameters
    int numData = 1;
    int tag;
    if (OPS_GetInt(&numData, &tag) != 0) {
        opserr << "ASDConstraintEquationElement ERROR invalid integer for $tag.\n";
        return 0;
    }
    double K;
    if (OPS_GetDouble(&numData, &K) != 0) {
        opserr << "ASDConstraintEquationElement ERROR invalid floating point number for $K.\n";
        return 0;
    }
    int cNode;
    if (OPS_GetInt(&numData, &cNode) != 0) {
        opserr << "ASDConstraintEquationElement ERROR invalid integer for $cNode.\n";
        return 0;
    }
    int cDof;
    if (OPS_GetInt(&numData, &cDof) != 0) {
        opserr << "ASDConstraintEquationElement ERROR invalid integer for $cDof.\n";
        return 0;
    }

    // retained nodes data
    numArgs = OPS_GetNumRemainingInputArgs();
    if (numArgs % 3 != 0) {
        opserr << "ASDConstraintEquationElement ERROR optional parameters for retained nodes "
            "<$rNode1 $rDof1 $rFact1 ... $rNodeN $rDofN $rFactN> should be multiple of 3 (found " << numArgs << ")\n" << descr;
        return 0;
    }
    int numTuples = numArgs / 3;
    // put the constrained node at the beginning
    ID rNodes(numTuples + 1);
    ID rDofs(numTuples + 1);
    Vector rFactors(numTuples + 1);
    rNodes(0) = cNode;
    rDofs(0) = cDof - 1;
    rFactors(0) = -1.0;
    // read retained nodes data
    int ivalue;
    double dvalue;
    for (int i = 0; i < numTuples; ++i) {
        if (OPS_GetInt(&numData, &ivalue) != 0) {
            opserr << "ASDConstraintEquationElement ERROR invalid integer for $rNode.\n";
            return 0;
        }
        rNodes(i + 1) = ivalue;
        if (OPS_GetInt(&numData, &ivalue) != 0) {
            opserr << "ASDConstraintEquationElement ERROR invalid integer for $rDof.\n";
            return 0;
        }
        rDofs(i + 1) = ivalue - 1;
        if (OPS_GetDouble(&numData, &dvalue) != 0) {
            opserr << "ASDConstraintEquationElement ERROR invalid floating point number for $rFactor.\n";
            return 0;
        }
        rFactors(i + 1) = dvalue;
    }

    // done
    return new ASDConstraintEquationElement(tag, rNodes, rDofs, rFactors, K);
}

ASDConstraintEquationElement::ASDConstraintEquationElement() 
    : Element(0, ELE_TAG_ASDConstraintEquationElement)
{
}

ASDConstraintEquationElement::ASDConstraintEquationElement(int tag, const ID& theNodes, const ID& theDofs, const Vector& theFactors, double K)
    : Element(tag, ELE_TAG_ASDConstraintEquationElement)
    , m_local_nodes(theNodes)
    , m_local_dofs(theDofs)
    , m_factors(theFactors)
    , m_K(K)
{
}

ASDConstraintEquationElement::~ASDConstraintEquationElement( )
{
}

const char* ASDConstraintEquationElement::getClassType(void) const
{
    return "ASDConstraintEquationElement";
}

void ASDConstraintEquationElement::setDomain(Domain* theDomain)
{
    // order nodes and keep unique values, because in the equations the retained nodes may apper multiple times.
    // key = node id
    // value = num dofs
    std::map<int, int> aux_map;
    for (int i = 0; i < m_local_nodes.Size(); ++i) {
        int node_id = m_local_nodes(i);
        // check node
        Node* node = theDomain->getNode(node_id);
        if (node == nullptr) {
            opserr << "ASDConstraintEquationElement Error in setDomain: node " << node_id << " does not exit in the domain\n";
            exit(-1);
        }
        // get number of dofs
        int ndofs = node->getNumberDOF();
        aux_map[node_id] = ndofs;
    }

    // count the total number of dofs
    m_num_dofs = 0;
    for (const auto& item : aux_map) {
        int ndofs = item.second;
        m_num_dofs += ndofs;
    }

    // now that they are ordered and unique, we can compute another map
    // with the starting local dof and the number of dofs for each node
    std::map<int, std::pair<int, int>> node_dof_map;
    {
        int pos = 0;
        for (const auto& item : aux_map) {
            int node_id = item.first;
            int ndof = item.second;
            node_dof_map[node_id] = std::make_pair(pos, ndof);
            pos += ndof;
        }
    }

    // compute mapping
    m_mapping.resize(m_local_dofs.Size());
    for (int i = 0; i < m_local_dofs.Size(); ++i) {
        int node_id = m_local_nodes(i);
        int node_dof = m_local_dofs(i);
        const auto& node_data = node_dof_map.at(node_id);
        int pos = node_data.first;
        int count = node_data.second;
        if (node_dof < 0 || node_dof >= count) {
            opserr << "ASDConstraintEquationElement Error in setDomain: DOF " << node_dof + 1 << " of node " << node_id << " is out of bounds.\n";
            exit(-1);
        }
        m_mapping(i) = node_dof + pos;
    }

    // store node pointers
    {
        std::size_t counter = 0;
        m_nodes.resize(node_dof_map.size());
        m_node_ids.resize(static_cast<int>(node_dof_map.size()));
        for (const auto& item : node_dof_map) {
            m_nodes[counter] = theDomain->getNode(item.first);
            m_node_ids(static_cast<int>(counter)) = item.first;
            ++counter;
        }
    }

    // compute initial displacement vector
    if (!m_U0_computed) {
        m_U0.resize(m_local_dofs.Size());
        m_U0 = getLocalDisplacements();
        m_U0_computed = true;
    }

    // call base class implementation
    DomainComponent::setDomain(theDomain);
}

void ASDConstraintEquationElement::Print(OPS_Stream& s, int flag)
{
    if (flag == -1) {
        int eleTag = this->getTag();
        s << "EL_ASDConstraintEquationElement\t" << eleTag << " :";
        for (int i = 0; i < m_node_ids.Size(); ++i)
            s << "\t" << m_node_ids(i);
        s << endln;
    }

    if (flag == OPS_PRINT_PRINTMODEL_JSON) {
        s << "\t\t\t{";
        s << "\"name\": " << this->getTag() << ", ";
        s << "\"type\": \"ASDConstraintEquationElement\", ";
        s << "\"nodes\": [";
        for (int i = 0; i < m_node_ids.Size(); ++i) {
            if (i > 0)
                s << ", ";
            s << m_node_ids(i);
        }
        s << "]}";
    }
}

int ASDConstraintEquationElement::getNumExternalNodes() const
{
    return m_node_ids.Size();
}

const ID& ASDConstraintEquationElement::getExternalNodes()
{
    return m_node_ids;
}

Node**
ASDConstraintEquationElement::getNodePtrs(void)
{
    return m_nodes.data();
}

int ASDConstraintEquationElement::getNumDOF()
{
    return m_num_dofs;
}

int ASDConstraintEquationElement::update()
{
    return 0;
}

int ASDConstraintEquationElement::revertToLastCommit()
{
    return 0;
}

const Matrix& ASDConstraintEquationElement::getTangentStiff()
{
    // output local matrix
    const Matrix& KL = getLocalMatrix();

    // output global matrix
    static Matrix K;
    K.resize(m_num_dofs, m_num_dofs);
    K.Zero();

    // copy in global dofset
    for (int i = 0; i < KL.noRows(); ++i) {
        int ig = m_mapping(i);
        for (int j = 0; j < KL.noCols(); ++j) {
            int jg = m_mapping(j);
            K(ig, jg) = KL(i, j);
        }
    }

    // done
    return K;
}

const Matrix& ASDConstraintEquationElement::getInitialStiff()
{
    return getTangentStiff();
}

const Matrix& ASDConstraintEquationElement::getMass()
{
    static Matrix M;
    M.resize(m_num_dofs, m_num_dofs);
    M.Zero();
    return M;
}

const Matrix& ASDConstraintEquationElement::getDamp()
{
    static Matrix C;
    C.resize(m_num_dofs, m_num_dofs);
    C.Zero();
    return C;
}

int
ASDConstraintEquationElement::addInertiaLoadToUnbalance(const Vector& accel)
{
    return 0;
}

const Vector& ASDConstraintEquationElement::getResistingForce()
{
    // output residual vector
    static Vector F;
    F.resize(m_num_dofs);

    // skip if computing reactions
    if (m_is_computing_reactions) {
        F.Zero();
        return F;
    }

    static Vector FL;
    FL.resize(m_local_dofs.Size());
    const Matrix& KL = getLocalMatrix();
    const Vector& UL = getLocalDisplacements();
    FL.addMatrixVector(0.0, KL, UL, 1.0);

    // copy in global dofset
    for (int i = 0; i < FL.Size(); ++i) {
        int ig = m_mapping(i);
        F(ig) = FL(i);
    }

    // done
    return F;
}

const Vector& ASDConstraintEquationElement::getResistingForceIncInertia()
{
    return getResistingForce();
}

int ASDConstraintEquationElement::addResistingForceToNodalReaction(int flag)
{
    m_is_computing_reactions = true;
    int result = Element::addResistingForceToNodalReaction(flag);
    m_is_computing_reactions = false;
    return result;
}

int ASDConstraintEquationElement::sendSelf(int commitTag, Channel& theChannel)
{
    int res = 0;
    int dataTag = getDbTag();

    // INT data 1
    // tag
    // number of local nodes (NL)
    // U0_computed
    int NL = m_local_nodes.Size();
    static ID idData1(3);
    idData1(0) = getTag();
    idData1(1) = NL;
    idData1(2) = static_cast<int>(m_U0_computed);
    res = theChannel.sendID(dataTag, commitTag, idData1);
    if (res < 0) {
        opserr << "WARNING ASDConstraintEquationElement::sendSelf() - " << this->getTag() << " failed to send ID 1\n";
        return res;
    }

    // INT data 2
    // NL local nodes
    // NL local dofs
    ID idData2(2 * NL);
    int pos = 0;
    for (int i = 0; i < NL; ++i) idData2(pos++) = m_local_nodes(i);
    for (int i = 0; i < NL; ++i) idData2(pos++) = m_local_dofs(i);
    res = theChannel.sendID(dataTag, commitTag, idData2);
    if (res < 0) {
        opserr << "WARNING ASDConstraintEquationElement::sendSelf() - " << this->getTag() << " failed to send ID 2\n";
        return res;
    }

    // DOUBLE data
    // K
    // NL factors
    // NL initial displacement
    Vector vectData(1 + 2 * NL);
    pos = 0;
    vectData(pos++) = m_K;
    for (int i = 0; i < NL; ++i) vectData(pos++) = m_factors(i);
    for (int i = 0; i < NL; ++i) vectData(pos++) = m_U0(i);
    res = theChannel.sendVector(dataTag, commitTag, vectData);
    if (res < 0) {
        opserr << "WARNING ASDConstraintEquationElement::sendSelf() - " << this->getTag() << " failed to send Vector\n";
        return res;
    }

    // done
    return res;
}

int ASDConstraintEquationElement::recvSelf(int commitTag, Channel& theChannel, FEM_ObjectBroker& theBroker)
{
    int res = 0;
    int dataTag = this->getDbTag();

    // INT data 1
    // tag
    // number of local nodes (NL)
    // U0_computed
    static ID idData1(3);
    res = theChannel.recvID(dataTag, commitTag, idData1);
    if (res < 0) {
        opserr << "WARNING ASDConstraintEquationElement::recvSelf() - " << this->getTag() << " failed to receive ID 1\n";
        return res;
    }
    setTag(idData1(0));
    int NL = idData1(1);
    m_U0_computed = static_cast<bool>(idData1(2));

    // INT data 2
    // NL local nodes
    // NL local dofs
    ID idData2(2 * NL);
    m_local_nodes.resize(NL);
    m_local_dofs.resize(NL);
    res = theChannel.recvID(dataTag, commitTag, idData2);
    if (res < 0) {
        opserr << "WARNING ASDConstraintEquationElement::recvSelf() - " << this->getTag() << " failed to receive ID 2\n";
        return res;
    }
    int pos = 0;
    for (int i = 0; i < NL; ++i) m_local_nodes(i) = idData2(pos++);
    for (int i = 0; i < NL; ++i) m_local_dofs(i) = idData2(pos++);

    // DOUBLE data
    // K
    // NL factors
    // NL initial displacement
    Vector vectData(1 + 2 * NL);
    m_factors.resize(NL);
    m_U0.resize(NL);
    res = theChannel.recvVector(dataTag, commitTag, vectData);
    if (res < 0) {
        opserr << "WARNING ASDConstraintEquationElement::recvSelf() - " << this->getTag() << " failed to receive Vector\n";
        return res;
    }
    pos = 0;
    m_K = vectData(pos++);
    for (int i = 0; i < NL; ++i) m_factors(i) = vectData(pos++);
    for (int i = 0; i < NL; ++i) m_U0(i) = vectData(pos++);

    // done
    return res;
}

const Vector& ASDConstraintEquationElement::getLocalDisplacements() const
{
    // global displacement vector
    static Vector U;
    U.resize(m_num_dofs);
    int counter = 0;
    for (Node* node : m_nodes) {
        const Vector& iu = node->getTrialDisp();
        for (int i = 0; i < iu.Size(); ++i) {
            U(counter++) = iu(i);
        }
    }

    // local displacement vector
    static Vector UL;
    UL.resize(m_local_dofs.Size());
    for (int i = 0; i < m_local_dofs.Size(); ++i)
        UL(i) = U(m_mapping(i));

    // remove initial displacements
    if (m_U0_computed) {
        UL.addVector(1.0, m_U0, -1.0);
    }
    return UL;
}

const Matrix& ASDConstraintEquationElement::getLocalMatrix() const
{
    // output local matrix
    static Matrix KL;
    KL.resize(m_local_dofs.Size(), m_local_dofs.Size());

    // compute
    static Matrix B;
    B.resize(1, m_local_dofs.Size());
    for (int i = 0; i < m_local_dofs.Size(); ++i)
        B(0, i) = m_factors(i);
    KL.addMatrixTransposeProduct(0.0, B, B, m_K);

    // done
    return KL;
}


