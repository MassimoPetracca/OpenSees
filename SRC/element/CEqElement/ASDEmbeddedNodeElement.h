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
// Notes:
//
// 

#ifndef ASDEmbeddedNodeElement_h
#define ASDEmbeddedNodeElement_h

#include <Element.h>
#include <ID.h>
#include <Vector.h>
#include <Matrix.h>
#include <vector>

class ASDEmbeddedNodeElement : public Element
{

public:

    // the geometric family of the embedding (host) element, deduced in
    // setDomain from the number of retained nodes and from ndm
    enum HostFamily {
        Fam_Unknown = 0,
        Fam_Tri = 1,     // 3 retained nodes (2D or 3D)
        Fam_Tet = 2,     // 4 retained nodes in 3D
        Fam_Quad = 3,    // 4 retained nodes in 2D
        Fam_Hexa = 4,    // 8 retained nodes in 3D
        Fam_Quad3D = 5   // 4 retained nodes in 3D (a shell or solid face)
    };

    // which dofs of the constrained node the element ties to the host field
    enum ConstraintMode {
        Mode_U = 0,   // translations only
        Mode_UR = 1,  // translations + rotations (skew part of the gradient; on a
                      // surface host the two bending rotations come from the slope
                      // of the transverse displacement, or, with -shearDeformable,
                      // from the interpolated nodal rotations of the host)
        Mode_UP = 2   // translations + pressure (u-p nodes)
    };

    // life cycle
    ASDEmbeddedNodeElement();
    ASDEmbeddedNodeElement(int tag, int cNode, const ID& rNodes, bool rot_flag, bool p_flag, double K, double KP, int shape_request = Fam_Unknown, bool shear_flag = false);
    virtual ~ASDEmbeddedNodeElement();

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

    // public methods for element output
    int sendSelf(int commitTag, Channel& theChannel);
    int recvSelf(int commitTag, Channel& theChannel, FEM_ObjectBroker& theBroker);

private:
    const Vector& getGlobalDisplacements() const;
    const Matrix& TRI_2D_U();
    const Matrix& TRI_2D_UP();
    const Matrix& TRI_2D_UR();
    const Matrix& TRI_3D_U();
    const Matrix& TRI_3D_UR();
    const Matrix& TRI_3D_UP();
    const Matrix& TET_3D_U();
    const Matrix& TET_3D_UR();
    const Matrix& TET_3D_UP();
    // isoparametric (non-simplex) hosts: the natural coordinate of the
    // constrained node is obtained with a Newton inversion of the map, and the
    // three constraint modes share one generic assembler
    const Matrix& QUAD_2D(int mode);
    const Matrix& HEX_3D(int mode);
    // a quadrilateral host in 3D is a surface, not a volume: the constraint is
    // written in the local frame of the (possibly warped) face, following the
    // same scheme as TRI_3D_UR
    const Matrix& QUAD_3D(int mode);
    // resolves m_family from the retained node count, from m_ndm and, when
    // those are not enough, from the user's -shape request
    int resolveFamily() const;

private:

    // the nodal ids, the first one is the constrained node,
    // the others are the retained nodes: 3 (triangle in 2D or shell triangle
    // in 3D), 4 (quadrilateral in 2D or 3D, or tetrahedron in 3D) or 8
    // (hexahedron)
    ID m_node_ids;
    // the nodes
    std::vector<Node*> m_nodes;
    // store the number of dimensions (2 or 3 are allowed)
    int m_ndm = 0;
    // geometric family of the host element (see HostFamily)
    int m_family = Fam_Unknown;
    // the host family requested with -shape, or Fam_Unknown to let the element
    // deduce it. It is only ever needed to tell a quadrilateral face from a
    // tetrahedron, the one case where 4 retained nodes in 3D are ambiguous
    int m_shape_request = Fam_Unknown;
    // total number of dofs
    int m_num_dofs = 0;
    // user input to constrain, if necessary, the rotation of the constrained node
    // if the constrained node has rotational DOFs
    bool m_rot_c_flag = false;
    // user input to constrain, if necessary, the pressure of the constrained node
    // if all nodes are U-P
    bool m_p_flag = false;
    // true if the constrained node has rotational DOFs and the user flag is true
    bool m_rot_c = false;
    // user input (-shearDeformable): on a surface host in 3D, tie the two
    // bending rotations of the constrained node to the interpolated nodal
    // rotations of the host instead of the slope of the transverse
    // displacement. The slope equals the rotation only for thin (Kirchhoff)
    // plates; a shear-deformable (Mindlin) shell carries theta = slope + gamma,
    // and its nodal rotations are the physical rotation of the fiber.
    // The drilling rotation keeps the in-plane skew gradient in either case:
    // it is not a director rotation, and hosts with a penalty drilling dof
    // would otherwise leak a weakly-constrained dof into the constraint.
    bool m_shear_flag = false;
    // true when m_shear_flag is accepted: -rot active, 3D surface host
    // (triangle, or quadrilateral face), every retained node with 6 dofs
    bool m_shear = false;
    // true if both constrained and retained nodes are U-P nodes
    bool m_up = false;
    // a vector containing the local id mapping for assembling
    // into the element matrix and vectors
    ID m_mapping;
    // stiffness penalty value to impose the constraint
    double m_K = 1.0e18;
    double m_KP = 1.0e18;
    // initial displacements
    Vector m_U0;
    bool m_U0_computed = false;

};

#endif // ASDEmbeddedNodeElement_h
