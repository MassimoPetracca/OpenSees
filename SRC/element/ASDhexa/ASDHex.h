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
// $Date: 2020/05/18 22:51:21 $

// Original implementation: 
//
// An 8-node solid hexahedral element for three-dimensional continuum analysis,
// based on an enhanced assumed strain / Petrov-Galerkin formulation (PG-EAS).
//
// It supports both linear and corotational kinematics.

#ifndef ASDHex_h
#define ASDHex_h

// ======================================================================
// ASDHEX_EAS_IMPERF -- EXPERIMENT, off by default. Define it to enable the
// -imperfection option, which multiplies each gauss point's stress AND tangent
// by 1 + a*r(tag, gp), with r a deterministic pseudo-random number in [-1, 1).
//
// Its purpose is to DISCRIMINATE, not to fix. With a softening material the
// element-internal equilibrium int(G_test^T sigma dV) = 0 admits a localized
// solution, and which gauss point localizes is decided by round-off. An
// imperfection removes that tie. So:
//   - if a = 1e-6 makes the analysis converge, the mechanism was the round-off
//     tie, and the enhanced-block regularization (see the A-op comment in the
//     condensation) is the proper cure;
//   - if a = 1e-6 only changes WHICH gauss point fails, the mechanism is the
//     definiteness of k_qq and this branch is dead weight.
//
// It is not shipped enabled, and it should not be. Scaling sigma and C by the
// same factor is exactly a material with E, ft, fc and Gf all scaled -- so an
// amplitude of 5% injects 5% of scatter into the ELASTIC stiffness -- the
// 'stresses' / 'strains' responses report the material's raw sigma rather than
// f*sigma, and the response depends on the element TAG, which any remesh
// renumbers. A strength-only perturbation would need a strengthScale parameter
// on the material plus Brick-style per-point setParameter routing here, which is
// not a change worth making to chase a stability problem.
// ======================================================================
// #define ASDHEX_EAS_IMPERF

#include <Element.h>
#include <ID.h>
#include <Vector.h>
#include <Matrix.h>
#include <nDmaterial.h>
#include <ASDSolidHexCorotationalTransformation.h>
 


class NDMaterial;
class Node;
class ElementalLoad;
class Renderer;
class Channel;
class FEM_ObjectBroker;
class OPS_Stream;
class Information;
class Parameter;

class ASDSolidHexCorotationalTransformation;
class Damping;
// per-element cache of the reference-geometry metric basis, defined in ASDHex.cpp
struct ASDSolidHexRefMetric;

class ASDSolidHex : public Element
{
public:
    // EAS Class
    class EASData {
    public:
		// enhanced parameters: trial, last committed, and the residual of the
		// enhanced equilibrium equation h = -int(G^T sigma dV)
		Vector alpha = Vector(12);
		Vector alpha_commit = Vector(12);
		Vector alpha_residual = Vector(12);

		// local displacements: U is the previous ITERATION's trial value (the
		// one updatePG_EAS differentiates to get dU), U_converged is the last
		// committed one. Same roles as ASDShellQ4's EASData::U/U_converged.
		Vector U = Vector(24);
		Vector U_converged = Vector(24);

		Matrix Kqq_inv = Matrix(12, 12);
		Matrix Kqu = Matrix(12, 24);
		Matrix Kuq = Matrix(24, 12);

		// ----------------------------------------------------------------
		// diagnostics of the enhanced solve, snapshotted PER ELEMENT.
		//
		// k_qq and its inverse live in the ASDSolidHexGlobals singleton, which
		// is process wide: the next element evaluated overwrites them. A
		// recorder reading them there would report whichever element ran last,
		// the same trap the display-values comment in calculateAll warns about.
		// These are copies taken while the data is still ours.
		//
		// TRANSIENT: recomputed at every evaluation, so they are deliberately
		// NOT counted in ASDSolidHex_EAS_DATA_SIZE and NOT serialized. Do not
		// add them there -- there is nothing to restore, recvSelf's first
		// evaluation refills them.
		// ----------------------------------------------------------------
		double rcond = 0.0;        // 1 / (||k_qq||_inf * ||k_qq^-1||_inf)
		double min_pivot = 0.0;    // min LDL^T pivot of sym(k_qq) / mean diagonal
		double stab_used = 0.0;    // the regularization factor s actually applied
    };
public:

    // life cycle
    ASDSolidHex();
    ASDSolidHex(
        int tag,
        int node1,
        int node2,
        int node3,
        int node4,
        int node5,
        int node6,
        int node7,
        int node8,
        NDMaterial* mat,
        bool corotational = false,
        Damping* damping = nullptr,
        const double* body = nullptr,
        double eas_stab = 0.0,
        double eas_penalty = 0.0,
        bool eas_auto = true
#ifdef ASDHEX_EAS_IMPERF
        , double imperfection = 0.0
#endif
    );
    virtual ~ASDSolidHex();

    // this element owns m_material[8], m_eas and m_transformation as raw
    // pointers, so the implicitly generated copy operations would double-free
    ASDSolidHex(const ASDSolidHex&) = delete;
    ASDSolidHex& operator = (const ASDSolidHex&) = delete;

    const char* getClassType(void) const { return "ASDSolidHex"; }

    // domain
    void setDomain(Domain* theDomain);
    // staged construction: re-capture the initial displacement offset
    void onActivate();
    void onDeactivate();

    // damping
    int setDamping(Domain* theDomain, Damping* damping);


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
    const Matrix& getTangentStiff();
    const Matrix& getInitialStiff();
    const Matrix& getMass();

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

    Response* setResponse(const char** argv, int argc, OPS_Stream& output);
    int getResponse(int responseID, Information& eleInfo);

    int setParameter(const char** argv, int argc, Parameter& param);

    // calculate the characteristic length for this element
    double getCharacteristicLength(void);

    // display 
    int displaySelf(Renderer&, int mode, float fact, const char** displayModes = 0, int numModes = 0);

private:

    // internal method to compute everything using switches...
    int calculateAll(Matrix& LHS, Vector& RHS, int options);
    void updatePG_EAS(const Vector& U);
    void initializePG_EAS();
    // fills m_U0 from the trial displacements (see m_U0)
    void captureInitialDisp();

private:

    static const int NumNodes = 8;
    static const int NDOF = 3 * NumNodes; // 24

    // nodal ids
    ID m_node_ids = ID(8);
    Node* nodePtrs[NumNodes] = { nullptr, nullptr, nullptr, nullptr,nullptr, nullptr, nullptr, nullptr }; //pointers to four nodes 

    NDMaterial* m_material[8];

    // one Damping per gauss point, as in ASDShellQ4. nullptr when the element was
    // created without -damp.
    Damping* m_damping[8];

    // body force per unit mass (e.g. gravity), set by the -b option and used by
    // addLoad for -brickSelfWeight / -selfWeight
    double m_body[3] = { 0.0, 0.0, 0.0 };

    // A-op: user floor on the regularization of the enhanced block, set by
    // -easStab. The gate in calculateAll raises s above this on its own whenever
    // k_qq is not positive definite, so 0 -- the default -- does NOT mean "no
    // regularization ever", it means "none unless the enhanced solve needs it".
    //
    // Only the OPERATOR is regularized: the enhanced residual stays exact, so the
    // converged answer is independent of this value. It is a knob on convergence,
    // not on the response. Initialized here so the FEM_ObjectBroker constructor
    // gets it right too.
    double m_eas_stab = 0.0;

    // A-pen: TRUE penalty on the enhanced modes, set by -easPenalty. Default 0.
    //
    // Unlike m_eas_stab this one enters the RESIDUAL, not just the operator, so it
    // CHANGES THE CONVERGED ANSWER. It adds an artificial elastic energy
    //     Pi_stab = (p/2) * int( (G_trial*alpha)^T C0 (G_trial*alpha) dV )
    // i.e. the enhanced modes now have to pay to switch on. That buys a genuine
    // energetic barrier against them localizing into one gauss point layer at
    // near-zero cost -- which operator regularization alone cannot provide -- and
    // pays for it by reintroducing part of the locking the EAS exists to remove.
    //
    // The patch test survives at any p (C3 gives alpha = 0 for a constant stress
    // state, and the penalty is proportional to alpha, so alpha = 0 stays an exact
    // solution). What degrades is BENDING.
    //
    // A constant is a parameter, not state: nothing to commit or revert. An
    // ADAPTIVE p would be different -- it would enter the residual, so the answer
    // would depend on the iterate history, and it would have to be frozen per step
    // with committed state. That is deliberately not offered.
    double m_eas_penalty = 0.0;

    // Whether the AUTOMATIC A-op arming is allowed to raise s on its own when a
    // gauss point's material tangent loses definiteness. Default true -- that is
    // the safety net, and turning it off gives back the pre-existing behaviour.
    //
    // It exists so the two stabilizations can be told apart experimentally: with
    // the automatic path on, A-op has already fixed a softening element by the
    // time -easPenalty gets a chance to do anything, so the penalty's own effect
    // is unobservable. -noEasAuto is what isolates it.
    //
    // Stored as a DOUBLE in the serialized option block rather than as a flag in
    // the ID: a new int would force ID(31) -> ID(32) in both sendSelf and
    // recvSelf, and recvSelf sizes its ID before recvID, so getting that pair out
    // of step is silent garbage rather than an error. Not worth the risk for one
    // bool.
    bool m_eas_auto = true;

#ifdef ASDHEX_EAS_IMPERF
    // EXPERIMENT (see the ASDHEX_EAS_IMPERF block in ASDHex.cpp): amplitude of the
    // per-gauss-point imperfection, and the 8 factors derived from it. The factors
    // are cached rather than recomputed because they depend only on the element tag
    // and the gauss point index, and caching makes them printable.
    double m_imperfection = 0.0;
    double m_imperfection_f[8] = { 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0 };
#endif

    // EASData class
    EASData* m_eas = nullptr;
    // vectors for applying load (allocated only if necessary)
    Vector* m_load = nullptr;

    // flag PG-EAS
    bool m_use_corotational;

    ASDSolidHexCorotationalTransformation* m_transformation;

	bool m_initialized;

    // Initial displacement offset for the LINEAR kinematics path, in global cs
    // (8 nodes x 3 translations): the nodal displacement present when the element
    // enters the domain - or when it is activated - subtracted from the trial
    // displacements so that an element born in an already displaced mesh starts
    // strain free. In the corotational path the equivalent offset lives inside
    // m_transformation. Captured under the same m_initialized latch as the
    // transformation's, and serialized like it, so a recvSelf does not re-capture
    // it from the already displaced nodes.
    Vector m_U0 = Vector(NDOF);

    // Cache of everything metric_basis::initialize_metric() and its
    // orthogonalize() produce. Both depend on the REFERENCE geometry only --
    // calculateAll fills X from Node::getCrds() in the linear and in the
    // corotational path alike, never from deformed coordinates -- so they are
    // computed once per element instead of once per Newton iteration. Measured
    // at 30.5 % (corotational) to 33.0 % (linear) of calculateAll before this.
    // Opaque here on purpose: the type needs skew_frame, which lives in the
    // translation unit. Not serialized -- recomputed on demand after recvSelf.
    ASDSolidHexRefMetric* m_ref = nullptr;
};

#endif // ASDHex_h
