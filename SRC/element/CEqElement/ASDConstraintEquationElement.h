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
// $Date: 2022/02/24 22:51:21 $

// Original implementation: Massimo Petracca (ASDEA)
//
//
// Notes:
//
// 

#ifndef ASDConstraintEquationElement_h
#define ASDConstraintEquationElement_h

#include <Element.h>
#include <ID.h>
#include <Vector.h>
#include <Matrix.h>
#include <vector>

class ASDConstraintEquationElement : public Element
{

public:

    // life cycle
    ASDConstraintEquationElement();
    ASDConstraintEquationElement(int tag, const ID &theNodes, const ID &theDofs, const Vector &theFactors, double K, double rhs);
    virtual ~ASDConstraintEquationElement();

    // domain
    const char* getClassType(void) const;
    void setDomain(Domain* theDomain);

    // print
    void Print(OPS_Stream& s, int flag);

    // methods dealing with nodes and number of external dof
    int getNumExternalNodes() const;
    const ID& getExternalNodes();
    Node** getNodePtrs();
    int getNumDOF();

    // methods dealing with committed state and update
    int update();
    int revertToLastCommit();

    // methods to return the current linearized stiffness,
    // damping and mass matrices
    const Matrix& getTangentStiff();
    const Matrix& getInitialStiff();
    const Matrix& getMass();
    const Matrix& getDamp();

    // methods for applying loads
    int addInertiaLoadToUnbalance(const Vector& accel);

    // methods for obtaining resisting force (force includes elemental loads)
    const Vector& getResistingForce();
    const Vector& getResistingForceIncInertia();

    // computation of reactions
    int addResistingForceToNodalReaction(int flag);

    // public methods for element output
    int sendSelf(int commitTag, Channel& theChannel);
    int recvSelf(int commitTag, Channel& theChannel, FEM_ObjectBroker& theBroker);

    // parameters
    int setParameter(const char** argv, int argc, Parameter& param);
    int updateParameter(int parameterID, Information& info);

private:
    const Vector& getLocalDisplacements() const;
    const Matrix& getLocalMatrix() const;

private:
    // nodes involved in the constraint equation, first one is the constrained
    ID m_local_nodes;
    // local dofs involved in the constraint equation
    ID m_local_dofs;
    // factors
    Vector m_factors;
    // unique set of nodes
    ID m_node_ids;
    // the ndoes
    std::vector<Node*> m_nodes;
    // total number of dofs of the global matrix
    int m_num_dofs = 0;
    // a vector containing the local id mapping for assembling
    // into the element matrix and vectors
    ID m_mapping;
    // stiffness penalty value to impose the constraint
    double m_K = 1.0e18;
    // rhs for non-homogeneous constraints
    double m_rhs = 0.0;
    // initial displacements
    Vector m_U0;
    bool m_U0_computed = false;
    // flag for computation of reactions
    bool m_is_computing_reactions = false;

};

#endif // ASDConstraintEquationElement_h
