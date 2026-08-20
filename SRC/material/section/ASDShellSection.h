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
// A layered shell section that improves on LayeredShellFiberSection in
// four ways that cannot be retrofitted there:
//
// - each ply carries its own Gauss-Lobatto rule (nip >= 3): the elastic
//   bending stiffness is integrated exactly at any nip (the composite
//   midpoint rule of LayeredShell is 1 - 1/nip^2 of exact), and a node
//   sits exactly on each ply face, where cracking and crushing start;
// - gaps between plies (sandwich skins, hollow cores);
// - rebar layers on a separate list, free to overlap the plies: each
//   layer is (uniaxial, equivalent thickness, angle, z) and is wrapped
//   internally in a PlateRebarMaterial, so the steel no longer needs to
//   be carved into the ply stack;
// - a section offset: every fiber integrates with lever zeta = z + offset,
//   producing the membrane-bending coupling in the section tangent.
//
// Transverse shear: either integrated from the ply fibers with a
// configurable correction factor (default 5/6; note LayeredShell applies
// none, so k = 1 reproduces it), or frozen elastic (-elasticShear), where
// no transverse shear strain reaches the fibers and the two channels
// carry a constant stiffness (from the initial ply tangents, or given).
//
// Sign conventions mirror LayeredShellFiberSection exactly (fiber strain
// eps = eps_m - zeta*kappa, moment resultants m = +int(zeta*sigma)), so
// this section is a drop-in for the same shell elements.
//
// Design and red-gate oracle:
// OpenSees-Testing/new-asd-elements/ASDShellSection/

#ifndef ASDShellSection_h
#define ASDShellSection_h

#include <Vector.h>
#include <Matrix.h>
#include <ID.h>
#include <NDMaterial.h>
#include <SectionForceDeformation.h>
#include <vector>

class UniaxialMaterial;

class ASDShellSection : public SectionForceDeformation {

public:

    // one entry of the ply/gap stack, bottom-to-top. mat == nullptr marks
    // a gap (nip is ignored there)
    struct StackItem {
        NDMaterial* mat = nullptr;
        double t = 0.0;
        int nip = 0;
    };
    // one rebar layer: equivalent thickness t_eq at stack coordinate z,
    // uniaxial oriented at angle (degrees, from the local x axis)
    struct RebarItem {
        UniaxialMaterial* mat = nullptr;
        double teq = 0.0;
        double angle = 0.0;
        double z = 0.0;
    };
    // how the two transverse shear channels are handled
    enum ShearMode {
        Shear_Integrated = 0,  // from the ply fibers, scaled by sqrt(k)
        Shear_ElasticAuto = 1, // frozen, k * sum(G_i t_i) from initial tangents
        Shear_ElasticUser = 2  // frozen, user-given section stiffnesses
    };

    ASDShellSection();
    ASDShellSection(int tag,
        const std::vector<StackItem>& stack,
        const std::vector<RebarItem>& rebars,
        double offset,
        double k_shear,
        int shear_mode,
        double S1, double S2,
        double rho_extra,
        const std::vector<double>& modifiers = std::vector<double>());
    virtual ~ASDShellSection();

    const char* getClassType(void) const { return "ASDShellSection"; }

    SectionForceDeformation* getCopy();
    double getRho();
    int getOrder() const;
    const ID& getType();

    int commitState();
    int revertToLastCommit();
    int revertToStart();

    int setTrialSectionDeformation(const Vector& strain_from_element);
    const Vector& getSectionDeformation();
    const Vector& getStressResultant();
    const Matrix& getSectionTangent();
    const Matrix& getInitialTangent();

    void Print(OPS_Stream& s, int flag);

    int sendSelf(int commitTag, Channel& theChannel);
    int recvSelf(int commitTag, Channel& theChannel, FEM_ObjectBroker& theBroker);

    Response* setResponse(const char** argv, int argc, OPS_Stream& s);
    int getResponse(int responseID, Information& info);

    int setParameter(const char** argv, int argc, Parameter& param);
    int updateParameter(int parameterID, Information& info);

private:

    // builds the fiber arrays from the stack/rebar description. Called by
    // the full constructor; recvSelf receives the fiber arrays directly.
    void buildFibers(const std::vector<StackItem>& stack,
        const std::vector<RebarItem>& rebars);
    // Shear_ElasticAuto: resolves S1/S2 from the initial ply tangents on
    // FIRST USE, not at construction: materials with fracture-energy
    // regularization need an active element (their lch) before they can be
    // queried, and at parse time there is none. getCopy resets the flag so
    // every element-owned copy resolves on its own fibers.
    void ensureAutoShear();

    // stack description, kept for Print/getCopy/serialization:
    // m_stack_t[i] thickness, m_stack_nip[i] > 0 for a ply, -1 for a gap
    std::vector<double> m_stack_t;
    std::vector<int> m_stack_nip;
    // rebar description (t_eq, angle, z), same purpose
    std::vector<double> m_reb_teq;
    std::vector<double> m_reb_angle;
    std::vector<double> m_reb_z;

    // the fibers: stack coordinate z (offset NOT included), thickness
    // weight, material (an owned PlateFiber-type copy), rebar flag
    std::vector<double> m_fib_z;
    std::vector<double> m_fib_w;
    std::vector<NDMaterial*> m_fib_mat;
    std::vector<char> m_fib_rebar;

    double m_h = 0.0;         // total stack height (gaps included)
    double m_offset = 0.0;    // stack midplane position in the element frame
    double m_k = 5.0 / 6.0;   // shear correction factor (integrated mode)
    int m_shear_mode = Shear_Integrated;
    double m_S1 = 0.0;        // frozen shear stiffnesses (elastic modes),
    double m_S2 = 0.0;        // channels (6) and (7) of the resultants
    bool m_S_computed = false; // Shear_ElasticAuto: S1/S2 resolved (lazily)
    // per-component stiffness modifiers (ETABS-style: f11 f22 f12 m11 m22
    // m12 v13 v23, the getType order). Applied as sqrt(c) on the strain
    // side and sqrt(c) on the stress side - the only placement that keeps
    // the tangent symmetric and the section variationally consistent when
    // the components are coupled (Poisson, offset). Diagonal terms scale
    // by c_i exactly, couplings by sqrt(c_i c_j).
    double m_mod[8] = { 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0 };
    double m_rho_extra = 0.0; // additional mass per unit area (e.g. rebar)

    Vector strainResultant;   // 8, the committed trial section deformation

    static Vector stressResultant;
    static Matrix tangent;
    static ID array;
};

#endif // ASDShellSection_h
