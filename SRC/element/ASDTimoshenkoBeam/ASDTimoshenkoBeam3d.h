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

// Original implementation:
//
// A 2-node displacement-based Timoshenko beam element for 3D analysis.
// Both the transverse displacements and the cross-section rotations are
// interpolated with linear shape functions, and all the section
// deformations (axial strain, curvatures, shear strains, torsional
// twist) are evaluated at a single mid-span Gauss point (uniform
// reduced integration) to avoid shear locking.
//
// The element works with any SectionForceDeformation providing the
// P, MZ, MY, VY, VZ and T response codes (e.g. Elastic with shear, a
// fiber section aggregated with shear/torsion responses, ...). A
// section missing one of those codes is rejected: without the shear
// terms this formulation would have zero-energy modes.
//
// The element uses the standard CrdTransf (basic system) machinery, so
// geometric nonlinearity is available today through the transformation:
// PDelta, or Corotational for large displacements (verified through full
// 2*pi rotations: rigid-body objectivity, roll-up, torsion, buckling).
// An element-embedded corotational formulation (Felippa EICR, as in
// ASDShellQ4/T3 and ASDHex) remains possible future work for uniformity,
// not the repair of a correctness gap.

#ifndef ASDTimoshenkoBeam3d_h
#define ASDTimoshenkoBeam3d_h

#include <Element.h>
#include <Matrix.h>
#include <Vector.h>
#include <ID.h>

class Node;
class SectionForceDeformation;
class CrdTransf;
class Response;
class Damping;

class ASDTimoshenkoBeam3d : public Element
{
  public:
    ASDTimoshenkoBeam3d(int tag, int nd1, int nd2,
                        SectionForceDeformation &section,
                        CrdTransf &coordTransf,
                        double rho = 0.0, int cMass = 0,
                        Damping *damping = 0);
    ASDTimoshenkoBeam3d();
    ~ASDTimoshenkoBeam3d();

    const char *getClassType(void) const {return "ASDTimoshenkoBeam3d";};

    int getNumExternalNodes(void) const;
    const ID &getExternalNodes(void);
    Node **getNodePtrs(void);

    int getNumDOF(void);
    void setDomain(Domain *theDomain);
    // staged construction: re-capture the initial displacement offset,
    // which lives in the coordinate transformation
    void onActivate(void);
    void onDeactivate(void);

    // public methods to set the state of the element
    int commitState(void);
    int revertToLastCommit(void);
    int revertToStart(void);

    // public methods to obtain stiffness, mass, damping and residual information
    int update(void);
    const Matrix &getTangentStiff(void);
    const Matrix &getInitialStiff(void);
    const Matrix &getMass(void);

    void zeroLoad();
    int addLoad(ElementalLoad *theLoad, double loadFactor);
    int addInertiaLoadToUnbalance(const Vector &accel);

    int setDamping(Domain *theDomain, Damping *damping);
    const Vector &getDampingForce(void);

    const Vector &getResistingForce(void);
    const Vector &getResistingForceIncInertia(void);

    // public methods for element output
    int sendSelf(int commitTag, Channel &theChannel);
    int recvSelf(int commitTag, Channel &theChannel, FEM_ObjectBroker &theBroker);
    int displaySelf(Renderer &theViewer, int displayMode, float fact,
                    const char **modes = 0, int numModes = 0);
    void Print(OPS_Stream &s, int flag = 0);

    Response *setResponse(const char **argv, int argc, OPS_Stream &s);
    int getResponse(int responseID, Information &eleInfo);

    int setParameter(const char **argv, int argc, Parameter &param);
    int updateParameter(int parameterID, Information &info);
    int activateParameter(int parameterID);

  protected:

  private:
    // All the element computations in the basic system live here,
    // controlled by the OPT_ flags (see the anonymous namespace in the
    // .cpp): update of the section state, basic stiffness (tangent or
    // initial) and basic forces.
    int calculateAll(Matrix &kb, Vector &q, int options);

    SectionForceDeformation *theSection; // the section at the single Gauss point (owned copy)
    CrdTransf *crdTransf;                // coordinate transformation (owned copy)
    Damping *theDamping;                 // damping object (owned copy, may be null)

    ID connectedExternalNodes;
    Node *theNodes[2];

    static Matrix K;  // element stiffness, damping and mass matrix
    static Vector P;  // element resisting force vector

    Vector Q;      // applied nodal loads
    Vector q;      // basic forces
    double q0[5];  // fixed end forces in basic system (no torsion)
    double p0[5];  // reactions in basic system (no torsion)

    double rho;    // mass density per unit length
    int cMass;     // consistent mass flag

    int parameterID;
};

#endif
