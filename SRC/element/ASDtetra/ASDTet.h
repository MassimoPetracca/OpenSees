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
// A 4-node tetrahedron companion to ASDHex: classic constant-strain kernel
// (one integration point, constant B, V = detJ/6) with linear or
// corotational (EICR) kinematics.
//
// No Petrov-Galerkin / EAS here, on purpose: the P1 tetrahedron is affine
// for any shape, so it keeps full linear completeness at any distortion
// and has no parasitic strains for PG to cure (unlike the isoparametric
// hexa). Its real limits are constant strain (poor in bending: refine, or
// use the 10-node tet / ASDHex) and volumetric locking near
// incompressibility (real cures need cross-element patches or an internal
// bubble; deliberately out of scope). Intended use: filler elements in
// hexa-dominant meshes.
//
// Design notes and the numpy red gate (which measured the Felippa
// consistent tangent to be the EXACT Hessian for this element):
// OpenSees-Testing/new-asd-elements/ASDTet/

#ifndef ASDTet_h
#define ASDTet_h

#include <Element.h>
#include <ID.h>
#include <Vector.h>
#include <Matrix.h>
#include <array>

class NDMaterial;
class ASDSolidTet4CorotationalTransformation;

class ASDTet : public Element
{

public:

    // life cycle
    ASDTet();
    ASDTet(int tag, int node1, int node2, int node3, int node4,
        NDMaterial* material, bool corotational, const double* body);
    virtual ~ASDTet();

    // domain
    const char* getClassType(void) const { return "ASDTet"; }
    void setDomain(Domain* theDomain);
    // staged construction: re-capture the initial displacement offset (m_U0)
    void onActivate();
    void onDeactivate();

    // print
    void Print(OPS_Stream& s, int flag);

    // nodes and dofs
    int getNumExternalNodes() const;
    const ID& getExternalNodes();
    Node** getNodePtrs();
    int getNumDOF();

    // state
    int update();
    int commitState();
    int revertToLastCommit();
    int revertToStart();

    // stiffness, mass, damping
    const Matrix& getTangentStiff();
    const Matrix& getInitialStiff();
    const Matrix& getMass();

    // loads
    void zeroLoad();
    int addLoad(ElementalLoad* theLoad, double loadFactor);
    int addInertiaLoadToUnbalance(const Vector& accel);

    // resisting force
    const Vector& getResistingForce();
    const Vector& getResistingForceIncInertia();

    // output
    int sendSelf(int commitTag, Channel& theChannel);
    int recvSelf(int commitTag, Channel& theChannel, FEM_ObjectBroker& theBroker);
    Response* setResponse(const char** argv, int argc, OPS_Stream& output);
    int getResponse(int responseID, Information& eleInfo);

private:

    // displacement since activation (U - U0), 12 components
    void computeDisplacements(Vector& d) const;
    // fills m_U0 from the trial displacements (see m_U0)
    void captureInitialDisp();
    // local kernel: strain update (OPT_UPDATE), local stiffness (tangent or
    // initial) and local internal force, on the deformational displacements
    int calculateAll(Matrix& K, Vector& f, int options);

    ID m_node_ids;
    std::array<Node*, 4> m_nodes = { nullptr, nullptr, nullptr, nullptr };
    NDMaterial* m_material = nullptr;      // owned ThreeDimensional copy
    bool m_use_corotational = false;
    ASDSolidTet4CorotationalTransformation* m_transformation = nullptr;

    // reference geometry: constant gradients (rows: dN_a/dX), volume,
    // strain-displacement matrix (6x12, Voigt exx eyy ezz gxy gyz gzx)
    Matrix m_g;
    Matrix m_B;
    double m_V = 0.0;

    // body force per unit mass (used by selfWeight-type element loads)
    double m_body[3] = { 0.0, 0.0, 0.0 };
    Vector m_P0;                            // accumulated element load vector

    // initial displacements (activation-time), for staged construction
    Vector m_U0;
    bool m_initialized = false;
};

#endif // ASDTet_h
