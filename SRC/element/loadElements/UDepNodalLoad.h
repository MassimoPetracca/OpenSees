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

#ifndef UDepNodalLoad_h
#define UDepNodalLoad_h

#include <Element.h>
#include <ID.h>
#include <Vector.h>
#include <Matrix.h>
#include <vector>

class UDepNodalLoad : public Element
{
public:
    class item_t {
    public:
        int fnode = 0;
        int fdof = 0;
        int unode = 0;
        int udof = 0;
        double uval = 0.0;
    };

public:

    // life cycle
    UDepNodalLoad();
    UDepNodalLoad(int tag, const std::vector<item_t>& items);
    virtual ~UDepNodalLoad();

    const char *getClassType(void) const {return "UDepNodalLoad";}
    
    // domain
    void setDomain(Domain* theDomain);

    // print
    void Print(OPS_Stream& s, int flag);

    // methods dealing with nodes and number of external dof
    int getNumExternalNodes() const;
    const ID& getExternalNodes();
    Node** getNodePtrs();
    int getNumDOF();

    // methods dealing with committed state and update
    int commitState();
    int revertToLastCommit();
    int revertToStart();
    int update();

    // methods to return the current linearized stiffness,
    // damping and mass matrices
    const Matrix& getTangentStiff();
    const Matrix& getInitialStiff();
    const Matrix& getMass();
    const Matrix& getDamp();

    // methods for applying loads
    void zeroLoad();
    int addLoad(ElementalLoad* theLoad, double loadFactor);
    int addInertiaLoadToUnbalance(const Vector& accel);

    // methods for obtaining resisting force (force includes elemental loads)
    const Vector& getResistingForce();
    const Vector& getResistingForceIncInertia();

    // public methods for element output
    int sendSelf(int commitTag, Channel& theChannel);
    int recvSelf(int commitTag, Channel& theChannel, FEM_ObjectBroker& theBroker);

    // public methods for response
    Response* setResponse(const char** argv, int argc, OPS_Stream& output);
    int getResponse(int responseID, Information& eleInfo);

    // public methods for parameters
    int setParameter(const char** argv, int argc, Parameter& param);
    int updateParameter(int parameterID, Information& info);

private:
    // the input items
    std::vector<item_t> m_items;
    // the coefficient matrix
    Matrix m_K;
    // initialization flag
    bool m_initialized = false;
    // nodal ids and node pointers
    ID m_node_ids;
    std::vector<Node*> m_nodes;
    // the total number of DOFs
    int m_num_dofs = 0;
};

#endif // UDepNodalLoad_h
