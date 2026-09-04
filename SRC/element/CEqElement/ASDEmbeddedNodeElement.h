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

class UniaxialMaterial;

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
    ASDEmbeddedNodeElement(int tag, int cNode, const ID& rNodes, bool rot_flag, bool p_flag, double K, double KP, int shape_request = Fam_Unknown, bool shear_flag = false, bool corot_flag = false,
        UniaxialMaterial* slip_mat = nullptr, int slip_node = 0, double KS = 0.0, const Vector* slip_x = nullptr,
        double slip_area = 1.0, const Vector* rot_axis = nullptr);
    virtual ~ASDEmbeddedNodeElement();

    // domain
    const char* getClassType(void) const;
    void setDomain(Domain* theDomain);
    // staged construction (elementActivate/elementDeactivate): a reactivated
    // element is reborn strain-free at the CURRENT configuration
    void onActivate();
    void onDeactivate();

    // print
    void Print(OPS_Stream& s, int flag);

    // methods dealing with nodes and number of external dof
    int getNumExternalNodes() const;
    const ID& getExternalNodes();
    Node** getNodePtrs();
    int getNumDOF();

    // methods dealing with committed state and update
    int update();
    int commitState();
    int revertToLastCommit();
    int revertToStart();

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

    // -slip: recorder access to the bond response (the element has no other
    // responses: without -slip it is a pure penalty constraint)
    Response* setResponse(const char** argv, int argc, OPS_Stream& output);
    int getResponse(int responseID, Information& eleInfo);

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
    // number of retained (host) nodes: the first node is the constrained one
    // and, with -slip, the last one is the real rebar node, not a host node
    int numRetained() const;
    // non-corotational embedding stiffness in the reduced local dofset
    // (dispatch to the family kernels), WITHOUT the -slip block: it is the
    // elastic part reused by both getTangentStiff and getResistingForce
    const Matrix& embedLocalStiffness();
    // -slip: current relative kinematics real->AUX, resolved on the local
    // triad (frozen m_slip_T0, or R*T0 with -corotational). Fills m_slip_g
    // and m_slip_axis, and returns the slip (first component).
    double slipComputeGap();
    // -slip, non-corotational path: add the zeroLength-equivalent block to
    // the FULL element matrix/vector (positions of the AUX and real node dofs
    // are known without the reduced mapping)
    void slipAddLinear(Matrix* K, Vector* F);
    // -corotational: lazy setup of the reference (activation-time) quantities,
    // and evaluation of the exact constraint g and its exact first variation B
    // in the reduced dofset. Gated by the numpy oracles in
    // OpenSees-Testing/new-asd-elements/ASDEmbeddedNodeElement/corot/.
    void corotSetup();
    void corotComputeBg(Matrix& B, Vector& g);
    // -rotAxis, corotational path: the axis resolved in the frame the
    // rotational rows of corotComputeBg are written in (global components on
    // a volume host, E0-local on a surface host). The LINEAR kernels resolve
    // it per call instead, each on its own face frame.
    void rotAxisRowFrame(double* kr) const;

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
    // -rotAxis (a modifier of -rot): tie ONLY the rotation about this axis
    // (unit, reference configuration - e.g. the axis of an embedded bar).
    // Implemented as a rank-1 weight on the 3 rotational rows of every kernel,
    // kU*(k (x) k), which equals the single projected row exactly:
    //     B_rot^T (k k^T) B_rot == (k^T B_rot)^T (k^T B_rot)
    // (same identity for the force), so the five UR kernels keep their rows
    // and only the weight changes. No bending clamp anywhere - the full -rot
    // reads the skew part of the HOST element's displacement gradient, which
    // is element-wise and invents bending on a bar between two cells - while
    // the rigid twist of the bar, otherwise a zero-energy mode (every bar
    // node lies on its own axis), stays held. Gated by verify_rot_axis.py in
    // OpenSees-Testing/new-asd-elements/ASDEmbeddedNodeElement/rotaxis/.
    bool m_rot_axis_flag = false;                  // user input
    bool m_rot_axis = false;                       // accepted: -rot active, 3D
    double m_rot_axis_v[3] = { 0.0, 0.0, 0.0 };
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
    // (Re)captures m_U0 at the current configuration. Clears m_U0_computed
    // first, because getGlobalDisplacements() subtracts m_U0 while that flag
    // is set.
    void captureInitialDisp();

    // -corotational: make the UR constraint exact under finite rotations of
    // the host by writing the linear kernel on deformational quantities in
    // the frame R = polar(F) of the host patch (the ASDhexa frame rule).
    // The frame is history-free; the only state is the slave TOTAL rotation,
    // tracked as a quaternion exactly the way the ASD shells track their
    // nodal quaternions (additive rotation dofs -> incremental composition).
    bool m_corot_flag = false;   // user input
    bool m_corot = false;        // accepted: -rot active, 3D host
    bool m_corot_init = false;   // reference data below is filled
    double m_qs[4] = { 1.0, 0.0, 0.0, 0.0 };      // slave quaternion (w,x,y,z)
    double m_rv[3] = { 0.0, 0.0, 0.0 };           // last additive rotation vector
    double m_qs_conv[4] = { 1.0, 0.0, 0.0, 0.0 };
    double m_rv_conv[3] = { 0.0, 0.0, 0.0 };
    // -corotational + -shearDeformable (surface hosts): the bending rows read
    // the deformational nodal rotations of the host, so each retained node
    // gets the same quaternion bookkeeping as the slave (4 + 3 doubles per
    // node, trial and converged). Sized in setDomain, preserved by recvSelf.
    std::vector<double> m_qa, m_rva, m_qa_conv, m_rva_conv;
    // reference (activation-time) data: shape functions at the material point,
    // cartesian gradients there, center gradients (frame rule), centroid,
    // local positions. Computed on X + U0, where F = I so R0 = I.
    Vector m_cN;    // 8
    Matrix m_cD;    // 8x3, dN/dX at xi_s
    Matrix m_cgc;   // 8x3, center gradients g_a
    Vector m_cc0;   // 3
    Matrix m_cY0;   // 8x3, X_a + U0_a - c0
    Vector m_cY0s;  // 3
    double m_ciK = 0.0; // penalty, m_K * cbrt(V) (solids) or m_K * sqrt(A) (surfaces)
    // surface hosts only: reference face frame (columns e1,e2,e3) and the 2D
    // cartesian gradients at the material point live in m_cE0 / m_cD(:,0:1);
    // for the Kabsch frame m_cgc rows hold the reference local positions Xh_a
    bool m_corot_surf = false;
    Matrix m_cE0;

    // -slip: absorb the rebar-slip zeroLength into the element. The real
    // rebar node is appended as the LAST external node; the AUX (constrained)
    // node stays embedded in the host. The uniaxial material acts on the
    // relative displacement along the local X (the bar axis), a stiff elastic
    // tie (m_KS, a raw [F/L] stiffness like the zeroLength penalty, NOT the
    // host-scaled -K pressure) acts on every other relative dof. With
    // -corotational the local X rotates with the host frame R.
    bool m_slip = false;            // -slip given
    bool m_slip_rot = false;        // both AUX and real node carry rotations
    UniaxialMaterial* m_slip_mat = nullptr; // owned copy of the tau-slip law
    double m_KS = 0.0;              // rigid-tie stiffness, used raw
    // bond area of the embedded node (the rebar node's tributary length times
    // the bar circumference): the slip material's stress AND tangent are
    // multiplied by it wherever they enter the mechanics, so the material can
    // be the tau-slip (bond STRESS vs slip) law itself, queried raw through
    // the 'bondStress' response. 1.0 by default = the material is a
    // force-slip law (the pre--slipArea scripts, which wrapped the tau law in
    // a 'Parallel ... -factors <area>', keep meaning what they meant).
    double m_slip_area = 1.0;
    Vector m_slip_x0;               // bar axis in the reference configuration
    // local triad (rows x0,y0,z0): y0/z0 are an arbitrary stable completion
    // built with the same rule STKO's frame_from_x uses, so the reported
    // transverse components match the legacy zeroLength assembly
    Matrix m_slip_T0;
    // -corotational rotational tie: the real node gets the same quaternion
    // bookkeeping as the slave (additive dofs -> incremental composition)
    double m_qr[4] = { 1.0, 0.0, 0.0, 0.0 };
    double m_rvr[3] = { 0.0, 0.0, 0.0 };
    double m_qr_conv[4] = { 1.0, 0.0, 0.0, 0.0 };
    double m_rvr_conv[3] = { 0.0, 0.0, 0.0 };
    // last computed local gap [slip, t1, t2, (rx, ry, rz)] and current axis,
    // kept for the recorder responses
    Vector m_slip_g;
    Vector m_slip_axis;

};

#endif // ASDEmbeddedNodeElement_h
