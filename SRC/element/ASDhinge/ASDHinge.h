/* ****************************************************************** **
**    OpenSees - Open System for Earthquake Engineering Simulation    **
**          Pacific Earthquake Engineering Research Center            **
** ****************************************************************** */

// $Revision: 1.0 $
// $Date: 2026/08/23 $

// Original implementation: Massimo Petracca (ASDEA)
//
// ASDHinge: a zero-length hinge for the plastic hinges of STKO PRO.
//
// It does what a zeroLength with six ASDHysteretic1D materials does, and three
// things it cannot:
//
//   1. an optional COROTATIONAL frame, so the hinge forces follow the
//      structure under large rotations instead of staying in the directions of
//      the reference configuration;
//   2. SIX FIXED SLOTS.  With zeroLength the -dir/-mat lists have variable
//      length, so recorder column 3 is Vz on one hinge and Mz on another;
//      here slot k is always local dof k, whatever the user assigned;
//   3. output under the Hinge keyword -- Hinge.Force and Hinge.Deformation,
//      instead of material.stress and material.strain, which for a uniaxial
//      material are generic words that say nothing about a hinge.
//
// And it leaves room for the fourth: the constitutive core is an
// ASDHingeKernel, whose tangent is 6x6.  Today it is diagonal; a coupled law
// (N-M domain, biaxial shear) plugs in without the element changing.
//
//   element ASDHinge $tag $iNode $jNode
//           <-mat $m1 $m2 ... -dir $d1 $d2 ...>
//           <-K $k1 $k2 $k3 $k4 $k5 $k6>
//           <-orient $x1 $x2 $x3 $y1 $y2 $y3>
//           <-corotational> <-doRayleigh> <-damp $tag>
//
// Six dofs per node are REQUIRED: with coincident nodes and translations only
// there is no frame to extract (every lever arm is zero, the polar
// decomposition of the fitted deformation gradient is degenerate), so
// -corotational would be meaningless.  Hinges live on beams, so this costs
// nothing.

#ifndef ASDHinge_h
#define ASDHinge_h

#include <Element.h>
#include <Matrix.h>
#include <Vector.h>
#include <ID.h>

#include "ASDHingeKernel.h"
#include "ASDHingeTransformation.h"

class Node;
class Channel;
class Response;
class Damping;

class ASDHinge : public Element
{
public:
    ASDHinge(int tag, int Nd1, int Nd2,
             const double R0[3][3],
             ASDHingeKernel* kernel,
             bool corotational,
             int doRayleigh = 0,
             Damping* damping = 0);
    ASDHinge();
    ~ASDHinge();

    const char* getClassType() const { return "ASDHinge"; }

    // connectivity
    int getNumExternalNodes() const;
    const ID& getExternalNodes();
    Node** getNodePtrs();
    int getNumDOF();
    void setDomain(Domain* theDomain);
    int setDamping(Domain* theDomain, Damping* damping);

    // state
    int commitState();
    int revertToLastCommit();
    int revertToStart();
    int update();

    // matrices
    const Matrix& getTangentStiff();
    const Matrix& getInitialStiff();
    const Matrix& getDamp();
    const Matrix& getMass();

    // loads
    void zeroLoad();
    int addLoad(ElementalLoad* theLoad, double loadFactor);
    int addInertiaLoadToUnbalance(const Vector& accel);

    // forces
    const Vector& getResistingForce();
    const Vector& getResistingForceIncInertia();

    // i/o
    int sendSelf(int commitTag, Channel& theChannel);
    int recvSelf(int commitTag, Channel& theChannel, FEM_ObjectBroker& theBroker);
    void Print(OPS_Stream& s, int flag = 0);

    Response* setResponse(const char** argv, int argc, OPS_Stream& s);
    int getResponse(int responseID, Information& eleInformation);

    int setParameter(const char** argv, int argc, Parameter& param);

    /// A zero-length element has no size, and the base class would return 0
    /// from the internodal distances -- which ASDHysteretic1D reads on its
    /// first setTrialStrain for the -autoRegularization divisor.
    double getCharacteristicLength() { return 1.0; }

    // staged construction
    void onActivate();
    void onDeactivate();

private:
    void assembleDampingForce(Vector& P);
    void releaseLimitStateProbes();
    /// true if at least one slot answers limitStateRatio
    bool buildLimitStateProbes();

private:
    ID m_connectedExternalNodes;
    Node* m_nodes[2];

    ASDHingeKernel* m_kernel;
    ASDHingeTransformation* m_transf;

    /// 0 = material damp tangent, 1 = Rayleigh
    int m_doRayleigh;
    Damping* m_damping;

    /// lazily built forwarders for the Hinge.LimitStateRatio response
    Response* m_lsProbe[6];
    bool m_lsProbeDone;

    /// 0 when the broker built us: recvSelf brings the baseline in and
    /// setDomain must not overwrite it.  Same role as ZeroLength::mInitialize.
    int m_initialize;

    // class wide working storage.  Vector and Matrix are dynamically sized in
    // OpenSees: owning them per element would still be fine, but the assembly
    // loop touches one element at a time and the family (ZeroLength,
    // ASDShellQ4) shares statics for exactly this reason.
    static Matrix s_K;
    static Vector s_P;
};

#endif // ASDHinge_h
