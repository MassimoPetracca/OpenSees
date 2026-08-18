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
        const double* body = nullptr
    );
    virtual ~ASDSolidHex();

    // this element owns m_material[8], m_eas and m_transformation as raw
    // pointers, so the implicitly generated copy operations would double-free
    ASDSolidHex(const ASDSolidHex&) = delete;
    ASDSolidHex& operator = (const ASDSolidHex&) = delete;

    const char* getClassType(void) const { return "ASDSolidHex"; }

    // domain
    void setDomain(Domain* theDomain);

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

    // EASData class
    EASData* m_eas = nullptr;
    // vectors for applying load (allocated only if necessary)
    Vector* m_load = nullptr;

    // flag PG-EAS
    bool m_use_corotational;

    ASDSolidHexCorotationalTransformation* m_transformation;

	bool m_initialized;

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
