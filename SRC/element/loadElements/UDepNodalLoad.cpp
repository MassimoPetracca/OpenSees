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
// $Date: 2025/03/06 22:51:21 $

// Original implementation: Massimo Petracca (ASDEA)
// A nodal load whose components depend linearly on the displacement
// via user-defined coefficients
// 

#include <UDepNodalLoad.h>
#include <Domain.h>
#include <Node.h>
#include <ErrorHandler.h>
#include <Channel.h>
#include <FEM_ObjectBroker.h>
#include <elementAPI.h>
#include <Renderer.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <map>

void *
OPS_UDepNodalLoad(void)
{
    static bool first_done = false;
    if (!first_done) {
        opserr << "Using UDepNodalLoad - Developed by: Massimo Petracca, Guido Camata, ASDEA Software Technology\n";
        first_done = true;
    }

    int numArgs = OPS_GetNumRemainingInputArgs();
    if (numArgs < 1) {
        opserr << "Want: element UDepNodalLoad $tag <-F $FN $FD $UN $UD $UV>";
        return 0;
    }

    int tag;
    int numData = 1;
    if (OPS_GetInt(&numData, &tag) != 0) {
        opserr << "WARNING invalid integer tag: element UDepNodalLoad \n";
        return 0;
    }

    // add your logic to parse the tcl or python command
    std::vector<UDepNodalLoad::item_t> items;

    while (OPS_GetNumRemainingInputArgs() > 0) {
        const char* type = OPS_GetString();
        // do something to build the ith item_t ...
        if (strcmp(type, "-F") == 0) {
            if (OPS_GetNumRemainingInputArgs() < 5) {
                opserr << "UDepNodalLoad Error: 5 args required for each -F keyword!\n";
                return nullptr;
            }
            UDepNodalLoad::item_t item;
            if (OPS_GetIntInput(&numData, &item.fnode) < 0) {
                opserr << "UDepNodalLoad Error: cannot parse FN\n";
                return nullptr;
            }
            if (OPS_GetIntInput(&numData, &item.fdof) < 0) {
                opserr << "UDepNodalLoad Error: cannot parse FD\n";
                return nullptr;
            }
            if (OPS_GetIntInput(&numData, &item.unode) < 0) {
                opserr << "UDepNodalLoad Error: cannot parse UN\n";
                return nullptr;
            }
            if (OPS_GetIntInput(&numData, &item.udof) < 0) {
                opserr << "UDepNodalLoad Error: cannot parse UD\n";
                return nullptr;
            }
            if (OPS_GetDoubleInput(&numData, &item.uval) < 0) {
                opserr << "UDepNodalLoad Error: cannot parse UV\n";
                return nullptr;
            }
            item.fdof--;
            item.udof--;
            items.push_back(item);
        }
    }
    
    // create a new instance of this element
    return new UDepNodalLoad(tag, items);
}

UDepNodalLoad::UDepNodalLoad() 
    : Element(0, ELE_TAG_UDepNodalLoad)
{
    // nothing to do here
}

UDepNodalLoad::UDepNodalLoad(int tag, const std::vector<item_t>& items)
    : Element(tag, ELE_TAG_UDepNodalLoad)
    , m_items(items)
{
    // nothing to do here
}

UDepNodalLoad::~UDepNodalLoad( )
{
    // no memory to clean for this element
}

void  UDepNodalLoad::setDomain(Domain* theDomain)
{
    // if domain is null
    if (theDomain == nullptr) {
        m_nodes.clear();
        // call base class implementation
        DomainComponent::setDomain(theDomain);
        return;
    }

    // collect a unique set of nodes and map them to a continuous 0-based index
    std::map<int, std::size_t> node_map;
    std::size_t node_counter = 0;
    for (const auto& item : m_items) {
        if (node_map.count(item.fnode) == 0)
            node_map[item.fnode] = node_counter++;
        if (node_map.count(item.unode) == 0)
            node_map[item.unode] = node_counter++;
    }

    // get node pointers
    m_nodes.resize(node_counter);
    m_node_ids.resize(node_counter);
    for (const auto& item : node_map) {
        int node_id = item.first;
        std::size_t node_loc = item.second;
        Node* node = theDomain->getNode(node_id);
        if (node == 0) {
            opserr << "UDepNodalLoad Error: node " << node_id << " not in domain\n";
            exit(-1);
        }
        m_nodes[node_loc] = node;
        m_node_ids[node_loc] = node_id;
    }

    // only if not already initialized from recvSelf
    if (!m_initialized) {

        // compute number of DOFs
        m_num_dofs = 0;
        std::vector<int> node_dof_offset(node_counter);
        for (std::size_t i = 0; i < node_counter; ++i) {
            Node* node = m_nodes[i];
            node_dof_offset[i] = m_num_dofs;
            m_num_dofs += node->getNumberDOF();
        }

        // allocate the stiffness matrix
        m_K.resize(m_num_dofs, m_num_dofs);
        m_K.Zero();

        // fill the stiffness matrix
        for (const auto& item : m_items) {
            // check fnode and fdof: obtain the row index
            std::size_t fnode_loc = node_map[item.fnode];
            Node* fnode = m_nodes[fnode_loc];
            if (item.fdof < 0 || item.fdof >= fnode->getNumberDOF()) {
                opserr << "UDepNodalLoad Error: F DOF " << item.fdof + 1 << " out of range (1, " << fnode->getNumberDOF() << ")\n";
                exit(-1);
            }
            int row_offset = node_dof_offset[fnode_loc];
            int row = row_offset + item.fdof;
            // check unode and udof: obtain the column index
            std::size_t unode_loc = node_map[item.unode];
            Node* unode = m_nodes[unode_loc];
            if (item.udof < 0 || item.udof >= unode->getNumberDOF()) {
                opserr << "UDepNodalLoad Error: U DOF " << item.udof + 1 << " out of range (1, " << unode->getNumberDOF() << ")\n";
                exit(-1);
            }
            int col_offset = node_dof_offset[unode_loc];
            int col = col_offset + item.udof;
            // add the UVal term to the stiffness matrix (use minus sign, because it's an external force!)
            m_K(row, col) -= item.uval;
            opserr << m_K << "\n";
        }

        // initialized
        m_initialized = true;
    }

    // call base class implementation
    DomainComponent::setDomain(theDomain);
}

void UDepNodalLoad::Print(OPS_Stream& s, int flag)
{
    if(flag == -1) {
        int eleTag = this->getTag();
        s << "EL_UDepNodalLoad\t" << eleTag << " :";
        for (int i = 0; i < m_node_ids.Size(); ++i)
            s << "\t" << m_node_ids(i);
        s << endln;
    }

    if (flag == OPS_PRINT_PRINTMODEL_JSON) {
        s << "\t\t\t{";
        s << "\"name\": " << this->getTag() << ", ";
        s << "\"type\": \"UDepNodalLoad\", ";
        s << "\"nodes\": [";
        for (int i = 0; i < m_node_ids.Size(); ++i) {
            if (i > 0)
                s << ", ";
            s << m_node_ids(i);
        }
        s << "]}";
    }
}

int  UDepNodalLoad::getNumExternalNodes() const
{
    return m_node_ids.Size();
}

const ID& UDepNodalLoad::getExternalNodes()
{
    return m_node_ids;
}

Node**
UDepNodalLoad::getNodePtrs(void)
{
    return m_nodes.data();
}

int UDepNodalLoad::getNumDOF()
{
    return m_num_dofs;
}

int  UDepNodalLoad::commitState()
{
    // nothing to do here, the element is linear
    return 0;
}

int  UDepNodalLoad::revertToLastCommit()
{
    // nothing to do here, the element is linear
    return 0;
}

int  UDepNodalLoad::revertToStart()
{
    // nothing to do here, the element is linear
    return 0;
}

int UDepNodalLoad::update()
{
    // nothing to do here, K and F are computed in the associated methods
    return 0;
}

const Matrix& UDepNodalLoad::getTangentStiff()
{
    // K has been precomputed in setDomain
    return m_K;
}

const Matrix& UDepNodalLoad::getInitialStiff()
{
    // K has been precomputed in setDomain
    // same as tangent, the element is linear
    return m_K;
}

const Matrix& UDepNodalLoad::getMass()
{
    // no mass
    static Matrix M;
    M.resize(m_num_dofs, m_num_dofs);
    M.Zero();
    return M;
}

const Matrix& UDepNodalLoad::getDamp()
{
    // no damping
    static Matrix D;
    D.resize(m_num_dofs, m_num_dofs);
    D.Zero();
    return D;
}

void  UDepNodalLoad::zeroLoad()
{
    // nothing to do here
}

int
UDepNodalLoad::addLoad(ElementalLoad* theLoad, double loadFactor)
{
    // nothing to do here
    opserr << "UDepNodalLoad::addLoad - load type unknown for ele with tag: " << this->getTag() << endln;
    return -1;
}

int
UDepNodalLoad::addInertiaLoadToUnbalance(const Vector& accel)
{
    // nothing to do here
    return 0;
}

const Vector& UDepNodalLoad::getResistingForce()
{
    // calculate RHS contribution of this element
    static Vector RHS;
    RHS.resize(m_num_dofs);

    // todo: add contribution
    static Vector U;
    U.resize(m_num_dofs);
    int counter = 0;
    for (Node* node : m_nodes) {
        const Vector& iu = node->getTrialDisp();
        for (int i = 0; i < iu.Size(); ++i) {
            U(counter++) = iu(i);
        }
    }
    RHS.addMatrixVector(0.0, m_K, U, 1.0);

    // done
    return RHS;
}

const Vector& UDepNodalLoad::getResistingForceIncInertia()
{
    // no mass and no damp, same as getResistingForce
    return getResistingForce();
}

int  UDepNodalLoad::sendSelf(int commitTag, Channel& theChannel)
{
    int res = 0;

    // note: we don't check for dataTag == 0 for Element
    // objects as that is taken care of in a commit by the Domain
    // object - don't want to have to do the check if sending data
    int dataTag = this->getDbTag();

    // count the number of double values we need
    int num_items = m_items.size();
    int num_nodes = m_node_ids.Size();
    int num_double =
        num_items * 5 + // items
        m_num_dofs * m_num_dofs + // K
        num_nodes; // node_ids

    // a counter
    int counter;

    // INT data
    // 1 tag +
    // 1 initialization flag +
    // 1 num dofs + 
    // 1 num items +
    // 1 num nodes
    static ID idData(5);
    counter = 0;
    idData(counter++) = getTag();
    idData(counter++) = static_cast<int>(m_initialized);
    idData(counter++) = m_num_dofs;
    idData(counter++) = num_items;
    idData(counter++) = num_nodes;

    res = theChannel.sendID(dataTag, commitTag, idData);
    if (res < 0) {
        opserr << "WARNING UDepNodalLoad::sendSelf() - " << this->getTag() << " failed to send ID\n";
        return res;
    }

    // DOUBLE data
    Vector vectData(num_double);
    counter = 0;
    // items
    for (const auto& item : m_items) {
        vectData(counter++) = static_cast<double>(item.fnode);
        vectData(counter++) = static_cast<double>(item.fdof);
        vectData(counter++) = static_cast<double>(item.unode);
        vectData(counter++) = static_cast<double>(item.udof);
        vectData(counter++) = item.uval;
    }
    // K
    for (int i = 0; i < m_num_dofs; ++i)
        for (int j = 0; j < m_num_dofs; ++j)
            vectData(counter++) = m_K(i, j);
    // node IDs
    for (int i = 0; i < m_node_ids.Size(); ++i)
        vectData(counter++) = static_cast<double>(m_node_ids(i));

    res = theChannel.sendVector(dataTag, commitTag, vectData);
    if (res < 0) {
        opserr << "WARNING UDepNodalLoad::sendSelf() - " << this->getTag() << " failed to send Vector\n";
        return res;
    }

    // done
    return res;
}

int  UDepNodalLoad::recvSelf(int commitTag, Channel& theChannel, FEM_ObjectBroker& theBroker)
{
    int res = 0;

    int dataTag = this->getDbTag();

    // a counter
    int counter;

    // INT data
    // 1 tag +
    // 1 initialization flag +
    // 1 num dofs + 
    // 1 num items +
    // 1 num nodes
    static ID idData(5);
    res = theChannel.recvID(dataTag, commitTag, idData);
    if (res < 0) {
        opserr << "WARNING UDepNodalLoad::recvSelf() - " << this->getTag() << " failed to receive ID\n";
        return res;
    }
    counter = 0;
    setTag(idData(counter++));
    m_initialized = static_cast<bool>(idData(counter++));
    m_num_dofs = idData(counter++);
    int num_items = idData(counter++);
    int num_nodes = idData(counter++);
    int num_double =
        num_items * 5 + // items
        m_num_dofs * m_num_dofs + // K
        num_nodes; // node_ids
    
    // DOUBLE data
    Vector vectData(num_double);
    res = theChannel.recvVector(dataTag, commitTag, vectData);
    if (res < 0) {
        opserr << "WARNING UDepNodalLoad::recvSelf() - " << this->getTag() << " failed to receive Vector\n";
        return res;
    }
    counter = 0;
    // items
    m_items.resize(static_cast<std::size_t>(num_items));
    for (std::size_t i = 0; i < m_items.size(); ++i) {
        auto& item = m_items[i];
        item.fnode = static_cast<int>(vectData(counter++));
        item.fdof = static_cast<int>(vectData(counter++));
        item.unode = static_cast<int>(vectData(counter++));
        item.udof = static_cast<int>(vectData(counter++));
        item.uval = vectData(counter++);
    }
    // K
    for (int i = 0; i < m_num_dofs; ++i)
        for (int j = 0; j < m_num_dofs; ++j)
            m_K(i, j) = vectData(counter++);
    // node IDs
    m_node_ids.resize(num_nodes);
    for (int i = 0; i < m_node_ids.Size(); ++i)
        m_node_ids(i) = static_cast<int>(vectData(counter++));
    
    // done
    return res;
}

Response*
UDepNodalLoad::setResponse(const char **argv, int argc, OPS_Stream &output)
{
    // nothing to record now
    return nullptr;
}

int
UDepNodalLoad::getResponse(int responseID, Information &eleInfo)
{
    // nothing to record now
    return -1;
}

int UDepNodalLoad::setParameter(const char** argv, int argc, Parameter& param)
{
    // no param to update now
    return -1;
}

int UDepNodalLoad::updateParameter(int parameterID, Information& info)
{
    // no param to update now
    return 0;
}
