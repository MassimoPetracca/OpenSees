/* ****************************************************************** **
**    OpenSees - Open System for Earthquake Engineering Simulation    **
**          Pacific Earthquake Engineering Research Center            **
** ****************************************************************** */

// $Revision: 1.0 $
// $Date: 2026/08/23 $

// Original implementation: Massimo Petracca (ASDEA)
//
// ASDHingeCorotationalTransformation
// ==================================
//
// The two nodes of a hinge are COINCIDENT, so every lever arm is zero: the
// translational EICR used by ASDHex and ASDTet is undefined here (Compute_S
// returns zero, and the polar decomposition of the fitted deformation
// gradient divides by a null determinant).  The frame can only come from the
// nodal rotations -- which is why the element requires six dofs per node.
//
// FRAME: the geodesic MIDPOINT of the two nodal triads,
//
//     Theta = log(R_I^T R_J)              (material relative rotation)
//     R_e   = R_I exp(Theta/2)
//
// It is symmetric on purpose: in STKO PRO the two nodes of a hinge are I and J
// by node ordering, not by geometry, so anchoring the frame to one of them
// would make the answer depend on which one got the lower tag.
//
// DEFORMATIONS
//
//     e_r = R0^T Theta
//     e_t = R0^T R_e^T d,     d = u_J - u_I
//
// e_r falls out of the midpoint choice for free: R_e^T R_J = R_I^T R_e =
// exp(Theta/2), so the pull-back of the relative rotation into the frame is
// Theta itself (the two exponentials share the axis and commute).  And Theta
// is a MATERIAL quantity -- under a superposed rigid rotation Q both triads
// are premultiplied by Q and R_I^T R_J is untouched -- so e_r is exactly
// objective with no pull-back at all.
//
// VARIATIONS.  OpenSees accumulates the nodal rotational dofs additively and
// this class composes them on the LEFT (as ASDShellQ4CorotationalTransformation
// does), so the dof variation IS the spatial spin: dR_I = spin(dth_I) R_I.
// With T = J_l, the left Jacobian of SO(3):
//
//     dTheta = T^-1(Theta) R_I^T (dth_J - dth_I)      =: Psi_th (dth_J - dth_I)
//     dw_e   = (I - Phi) dth_I + Phi dth_J,   Phi = 1/2 R_I C(Theta) R_I^T
//                                             C   = T(Theta/2) T^-1(Theta)
//     de_t   = A^T [ dd + spin(d) dw_e ],     A = R_e R0
//     de_r   = R0^T dTheta
//
// The tangent is the EXACT linearization of B^T fL, geometric terms included.
// It is NOT symmetric away from equilibrium: with rotational dofs parametrised
// by left incremental rotations the second variation of the energy is
// symmetric only where the first variation vanishes.  Measured asymmetry at a
// strongly deformed state: 1.8e-2 of the symmetric part.  Do not use this
// element with a symmetric-storage system (ProfileSPD, BandSPD): they read the
// lower triangle only and would silently discard half of the geometric terms.
//
// Everything here is gated by the oracle
// OpenSees-Testing/asd-hinge-element/verify_asdhinge_corot.py: 204 checks,
// including K == FD(p) at every state to the finite-difference floor.  Six
// deliberate injections were measured to break it; dropping the whole
// geometric part is a 36% error on K, and freezing Phi in the linearization
// (the tempting simplification) is 3.9e-4 -- small, but a hundred times the
// FD floor, and it vanishes only when the translational gap does.

#ifndef ASDHingeCorotationalTransformation_h
#define ASDHingeCorotationalTransformation_h

#include "ASDHingeTransformation.h"

// --------------------------------------------------------------------------
// quaternions.  Four doubles per node instead of nine, and renormalisable:
// a rotation matrix accumulated by repeated multiplication drifts off SO(3).
// --------------------------------------------------------------------------

namespace ASDHingeUtils {

inline void quatIdentity(double* q)
{
    q[0] = 1.0; q[1] = 0.0; q[2] = 0.0; q[3] = 0.0;
}

inline void quatNormalize(double* q)
{
    double n = std::sqrt(q[0] * q[0] + q[1] * q[1] + q[2] * q[2] + q[3] * q[3]);
    if (n > 0.0) {
        n = 1.0 / n;
        q[0] *= n; q[1] *= n; q[2] *= n; q[3] *= n;
    }
}

/// q = (w, x, y, z) of the rotation vector v
inline void quatFromRotationVector(const double* v, double* q)
{
    double a = norm3(v);
    if (a < 1.0e-12) {
        q[0] = 1.0;
        q[1] = 0.5 * v[0];
        q[2] = 0.5 * v[1];
        q[3] = 0.5 * v[2];
    }
    else {
        double s = std::sin(0.5 * a) / a;
        q[0] = std::cos(0.5 * a);
        q[1] = s * v[0];
        q[2] = s * v[1];
        q[3] = s * v[2];
    }
    quatNormalize(q);
}

/// c = a * b (Hamilton, same convention as ASDMath::ASDQuaternion)
inline void quatMul(const double* a, const double* b, double* c)
{
    double w = a[0] * b[0] - a[1] * b[1] - a[2] * b[2] - a[3] * b[3];
    double x = a[0] * b[1] + a[1] * b[0] + a[2] * b[3] - a[3] * b[2];
    double y = a[0] * b[2] + a[2] * b[0] + a[3] * b[1] - a[1] * b[3];
    double z = a[0] * b[3] + a[3] * b[0] + a[1] * b[2] - a[2] * b[1];
    c[0] = w; c[1] = x; c[2] = y; c[3] = z;
}

inline void quatToMatrix(const double* q, double R[3][3])
{
    double w = q[0], x = q[1], y = q[2], z = q[3];
    R[0][0] = 1.0 - 2.0 * (y * y + z * z);
    R[0][1] = 2.0 * (x * y - w * z);
    R[0][2] = 2.0 * (x * z + w * y);
    R[1][0] = 2.0 * (x * y + w * z);
    R[1][1] = 1.0 - 2.0 * (x * x + z * z);
    R[1][2] = 2.0 * (y * z - w * x);
    R[2][0] = 2.0 * (x * z - w * y);
    R[2][1] = 2.0 * (y * z + w * x);
    R[2][2] = 1.0 - 2.0 * (x * x + y * y);
}

} // namespace ASDHingeUtils


// --------------------------------------------------------------------------
// isotropic functions of a rotation vector: f(Th) = a0 I + a1 X + a2 X^2,
// X = skew(Th).  g[i] = (d a[i] / d th) / th, finite at th -> 0.
// --------------------------------------------------------------------------

/**
The series threshold is NOT cosmetic.  g1 of T is
(th sin th - 2(1-cos th))/th^4, whose numerator cancels down to th^4/12: with
the textbook 1-cos th the absolute rounding is eps*1, and the coefficient is
wrong in the SECOND DIGIT at th = 1e-3.  Writing 1-cos th as 2 sin^2(th/2)
makes the rounding scale with th^2, and the series covers the rest.  Measured
in the oracle: the finite-difference gate failed by 2.8e-2 at th = 1e-3 with
the naive form and a 1e-4 threshold.
*/
#define ASDHINGE_TH_SMALL 0.1

struct ASDHingeIso
{
    double a[3];
    double g[3];
    double th;

    /// 1 - cos(th), without the cancellation
    static double omc(double th)
    {
        double s = std::sin(0.5 * th);
        return 2.0 * s * s;
    }

    /// J_l(Th)
    static ASDHingeIso T(double th)
    {
        ASDHingeIso r;
        r.th = th;
        r.a[0] = 1.0;
        r.g[0] = 0.0;
        double t2 = th * th;
        if (th < ASDHINGE_TH_SMALL) {
            double t4 = t2 * t2;
            double t6 = t4 * t2;
            r.a[1] = 0.5 - t2 / 24.0 + t4 / 720.0 - t6 / 40320.0;
            r.a[2] = 1.0 / 6.0 - t2 / 120.0 + t4 / 5040.0 - t6 / 362880.0;
            r.g[1] = -1.0 / 12.0 + t2 / 180.0 - t4 / 6720.0 + t6 / 453600.0;
            r.g[2] = -1.0 / 60.0 + t2 / 1260.0 - t4 / 60480.0 + t6 / 4989600.0;
        }
        else {
            double s = std::sin(th);
            double c = omc(th);
            r.a[1] = c / t2;
            r.a[2] = (th - s) / (t2 * th);
            r.g[1] = (th * s - 2.0 * c) / (t2 * t2);
            r.g[2] = (th * c - 3.0 * (th - s)) / (t2 * t2 * th);
        }
        return r;
    }

    /// J_l^-1(Th).  The cot form is deliberate: the textbook
    /// 1/th^2 - (1+cos)/(2 th sin) is 0/0 at th = pi.
    static ASDHingeIso Tinv(double th)
    {
        ASDHingeIso r;
        r.th = th;
        r.a[0] = 1.0;
        r.a[1] = -0.5;
        r.g[0] = 0.0;
        r.g[1] = 0.0;
        double t2 = th * th;
        if (th < ASDHINGE_TH_SMALL) {
            double t4 = t2 * t2;
            double t6 = t4 * t2;
            r.a[2] = 1.0 / 12.0 + t2 / 720.0 + t4 / 30240.0 + t6 / 1209600.0;
            r.g[2] = 1.0 / 360.0 + t2 / 7560.0 + t4 / 201600.0 + t6 / 5987520.0;
        }
        else {
            double u = 0.5 * th;
            double su = std::sin(u);
            double cot = std::cos(u) / su;
            r.a[2] = (1.0 - u * cot) / t2;
            // h(u) = (1 - u cot u)/(4 u^2) = a2 ;  d a2 / d th = h'(u) / 2
            double hp = (u * cot + u * u / (su * su) - 2.0) / (4.0 * u * u * u);
            r.g[2] = 0.5 * hp / th;
        }
        return r;
    }

    /// T(Th/2), as a polynomial in X = skew(Th), NOT in skew(Th/2).
    /// The coefficients are those of T at u = th/2, rescaled; du/dth = 1/2
    /// gives g1(u)/8 and g2(u)/16.  Deriving them again by hand would be
    /// three more chances to slip.
    static ASDHingeIso Thalf(double th)
    {
        ASDHingeIso t = T(0.5 * th);
        ASDHingeIso r;
        r.th = th;
        r.a[0] = 1.0;
        r.a[1] = t.a[1] / 2.0;
        r.a[2] = t.a[2] / 4.0;
        r.g[0] = 0.0;
        r.g[1] = t.g[1] / 8.0;
        r.g[2] = t.g[2] / 16.0;
        return r;
    }

    /// product of two isotropic functions, using X^3 = -th^2 X
    static ASDHingeIso mul(const ASDHingeIso& A, const ASDHingeIso& B)
    {
        ASDHingeIso r;
        r.th = A.th;
        double t2 = A.th * A.th;
        double a0 = A.a[0], a1 = A.a[1], a2 = A.a[2];
        double b0 = B.a[0], b1 = B.a[1], b2 = B.a[2];
        double ga0 = A.g[0], ga1 = A.g[1], ga2 = A.g[2];
        double gb0 = B.g[0], gb1 = B.g[1], gb2 = B.g[2];
        r.a[0] = a0 * b0;
        r.a[1] = a0 * b1 + a1 * b0 - t2 * (a1 * b2 + a2 * b1);
        r.a[2] = a0 * b2 + a1 * b1 + a2 * b0 - t2 * a2 * b2;
        r.g[0] = ga0 * b0 + a0 * gb0;
        r.g[1] = ga0 * b1 + a0 * gb1 + ga1 * b0 + a1 * gb0
               - 2.0 * (a1 * b2 + a2 * b1)
               - t2 * (ga1 * b2 + a1 * gb2 + ga2 * b1 + a2 * gb1);
        r.g[2] = ga0 * b2 + a0 * gb2 + ga1 * b1 + a1 * gb1
               + ga2 * b0 + a2 * gb0
               - 2.0 * a2 * b2
               - t2 * (ga2 * b2 + a2 * gb2);
        return r;
    }

    /// C(Th) = T(Th/2) T^-1(Th)
    static ASDHingeIso C(double th)
    {
        return mul(Thalf(th), Tinv(th));
    }

    ASDHingeIso transposed() const
    {
        ASDHingeIso r = *this;
        r.a[1] = -r.a[1];
        r.g[1] = -r.g[1];
        return r;
    }

    void matrix(const double* Th, double M[3][3]) const
    {
        double X[3][3], X2[3][3];
        ASDHingeUtils::skew3(Th, X);
        ASDHingeUtils::mul3(X, X, X2);
        for (int i = 0; i < 3; ++i)
            for (int j = 0; j < 3; ++j)
                M[i][j] = (i == j ? a[0] : 0.0) + a[1] * X[i][j] + a[2] * X2[i][j];
    }

    void apply(const double* Th, const double* v, double* r) const
    {
        double Txv[3], TxTxv[3];
        ASDHingeUtils::cross3(Th, v, Txv);
        ASDHingeUtils::cross3(Th, Txv, TxTxv);
        for (int i = 0; i < 3; ++i)
            r[i] = a[0] * v[i] + a[1] * Txv[i] + a[2] * TxTxv[i];
    }

    /// d( f(Th) v ) / d Th, v fixed
    void dapply(const double* Th, const double* v, double D[3][3]) const
    {
        double Txv[3], TxTxv[3], Sv[3][3];
        ASDHingeUtils::cross3(Th, v, Txv);
        ASDHingeUtils::cross3(Th, Txv, TxTxv);
        ASDHingeUtils::skew3(v, Sv);
        double Tdotv = ASDHingeUtils::dot3(Th, v);
        for (int i = 0; i < 3; ++i) {
            for (int j = 0; j < 3; ++j) {
                D[i][j] = g[0] * v[i] * Th[j]
                        + g[1] * Txv[i] * Th[j]
                        + g[2] * TxTxv[i] * Th[j]
                        - a[1] * Sv[i][j]
                        // d/dTh [ Th x (Th x v) ] = Th (x) v + (Th.v) I - 2 v (x) Th
                        + a[2] * (Th[i] * v[j] + (i == j ? Tdotv : 0.0)
                                  - 2.0 * v[i] * Th[j]);
            }
        }
    }
};


class ASDHingeCorotationalTransformation : public ASDHingeTransformation
{
public:
    ASDHingeCorotationalTransformation()
        : ASDHingeTransformation()
    {
        this->resetRotations();
    }

    bool isLinear() const { return false; }

    // -- life cycle --------------------------------------------------------

    void revertToStart()
    {
        ASDHingeTransformation::revertToStart();
        this->resetRotations();
    }

    void commit()
    {
        for (int k = 0; k < 2; ++k) {
            for (int i = 0; i < 4; ++i) m_qn_cvg[k][i] = m_qn[k][i];
            for (int i = 0; i < 3; ++i) m_rv_cvg[k][i] = m_rv[k][i];
        }
    }

    void revertToLastCommit()
    {
        for (int k = 0; k < 2; ++k) {
            for (int i = 0; i < 4; ++i) m_qn[k][i] = m_qn_cvg[k][i];
            for (int i = 0; i < 3; ++i) m_rv[k][i] = m_rv_cvg[k][i];
        }
    }

    void forceCaptureInitialDisp()
    {
        // the baseline moves to the current configuration, so the triads must
        // go back to the identity as well: otherwise a hinge activated in a
        // construction stage is born with a finite e_r
        ASDHingeTransformation::forceCaptureInitialDisp();
        this->resetRotations();
    }

    /// accumulate the nodal triads from the dof increments.  Must be called
    /// once per state update and only from the element update(): the
    /// composition is incremental, calling it twice applies it twice.
    void update(const Vector& UG)
    {
        for (int k = 0; k < 2; ++k) {
            int base = 6 * k + 3;
            double cur[3], incr[3], dq[4], q[4];
            for (int i = 0; i < 3; ++i) {
                cur[i] = UG(base + i);
                incr[i] = cur[i] - m_rv[k][i];
                m_rv[k][i] = cur[i];
            }
            ASDHingeUtils::quatFromRotationVector(incr, dq);
            ASDHingeUtils::quatMul(dq, m_qn[k], q);
            ASDHingeUtils::quatNormalize(q);
            for (int i = 0; i < 4; ++i) m_qn[k][i] = q[i];
        }
    }

    // -- kinematics --------------------------------------------------------

    void computeDeformations(const Vector& UG, Vector& e)
    {
        Frame f;
        this->computeFrame(UG, f);
        double et[3];
        ASDHingeUtils::matTVec3(f.A, f.d, et);
        for (int j = 0; j < 3; ++j) {
            e(j) = et[j];
            e(3 + j) = m_R0[0][j] * f.Th[0] + m_R0[1][j] * f.Th[1]
                     + m_R0[2][j] * f.Th[2];
        }
    }

    void computeB(const Vector& UG, Matrix& B)
    {
        Frame f;
        this->computeFrame(UG, f);
        this->assembleB(f, B);
    }

    void transformToGlobal(const Vector& UG, const Vector& fL,
                           const Matrix& kL, Matrix& LHS, Vector& RHS,
                           bool LHSrequired)
    {
        Frame f;
        this->computeFrame(UG, f);

        Matrix& B = ASDHingeGlobals::instance().B;
        this->assembleB(f, B);
        RHS.addMatrixTransposeVector(0.0, B, fL, 1.0);
        if (!LHSrequired)
            return;

        Matrix& BtK = ASDHingeGlobals::instance().BtK;
        BtK.addMatrixTransposeProduct(0.0, B, kL, 1.0);
        LHS.addMatrixProduct(0.0, BtK, B, 1.0);
        this->addGeometricStiffness(f, fL, LHS);
    }

    // -- serialization -----------------------------------------------------

    int internalDataSize() const { return 12 + 8 + 6; }

    void saveInternalData(Vector& v, int pos) const
    {
        ASDHingeTransformation::saveInternalData(v, pos);
        int p = pos + 12;
        for (int k = 0; k < 2; ++k) {
            for (int i = 0; i < 4; ++i) v(p++) = m_qn_cvg[k][i];
            for (int i = 0; i < 3; ++i) v(p++) = m_rv_cvg[k][i];
        }
    }

    void restoreInternalData(const Vector& v, int pos)
    {
        ASDHingeTransformation::restoreInternalData(v, pos);
        int p = pos + 12;
        for (int k = 0; k < 2; ++k) {
            for (int i = 0; i < 4; ++i) m_qn_cvg[k][i] = v(p++);
            for (int i = 0; i < 3; ++i) m_rv_cvg[k][i] = v(p++);
        }
        this->revertToLastCommit();
    }

private:

    struct Frame
    {
        double Th[3];          // log(R_I^T R_J)
        double d[3];           // u_J - u_I
        double RI[3][3];
        double Re[3][3];
        double A[3][3];        // R_e R0
        double Psi_th[3][3];   // T^-1(Th) R_I^T
        double Phi[3][3];      // 1/2 R_I C R_I^T
        ASDHingeIso C;
        double th;
    };

    void resetRotations()
    {
        for (int k = 0; k < 2; ++k) {
            ASDHingeUtils::quatIdentity(m_qn[k]);
            ASDHingeUtils::quatIdentity(m_qn_cvg[k]);
            for (int i = 0; i < 3; ++i) {
                m_rv[k][i] = 0.0;
                m_rv_cvg[k][i] = 0.0;
            }
        }
    }

    void computeFrame(const Vector& UG, Frame& f) const
    {
        using namespace ASDHingeUtils;

        double RJ[3][3], Lam[3][3], H[3][3], Tmp[3][3];
        quatToMatrix(m_qn[0], f.RI);
        quatToMatrix(m_qn[1], RJ);

        // Theta = log(R_I^T R_J)
        mulTN3(f.RI, RJ, Lam);
        logmap3(Lam, f.Th);
        f.th = norm3(f.Th);

        // d = u_J - u_I
        for (int i = 0; i < 3; ++i)
            f.d[i] = UG(6 + i) - UG(i);

        // R_e = R_I exp(Theta/2),  A = R_e R0
        double half[3] = { 0.5 * f.Th[0], 0.5 * f.Th[1], 0.5 * f.Th[2] };
        expmap3(half, H);
        mul3(f.RI, H, f.Re);
        mul3(f.Re, m_R0, f.A);

        // Psi_th = T^-1(Theta) R_I^T
        ASDHingeIso Ti = ASDHingeIso::Tinv(f.th);
        Ti.matrix(f.Th, Tmp);
        for (int i = 0; i < 3; ++i)
            for (int j = 0; j < 3; ++j) {
                double s = 0.0;
                for (int k = 0; k < 3; ++k) s += Tmp[i][k] * f.RI[j][k];
                f.Psi_th[i][j] = s;
            }

        // Phi = 1/2 R_I C R_I^T
        f.C = ASDHingeIso::C(f.th);
        f.C.matrix(f.Th, Tmp);
        double RC[3][3];
        mul3(f.RI, Tmp, RC);
        for (int i = 0; i < 3; ++i)
            for (int j = 0; j < 3; ++j) {
                double s = 0.0;
                for (int k = 0; k < 3; ++k) s += RC[i][k] * f.RI[j][k];
                f.Phi[i][j] = 0.5 * s;
            }
    }

    void assembleB(const Frame& f, Matrix& B) const
    {
        using namespace ASDHingeUtils;
        double AT[3][3], Sd[3][3], ATSd[3][3], Psi[3][3], ImPhi[3][3];
        double blockI[3][3], blockJ[3][3];

        for (int i = 0; i < 3; ++i)
            for (int j = 0; j < 3; ++j)
                AT[i][j] = f.A[j][i];

        skew3(f.d, Sd);
        mul3(AT, Sd, ATSd);

        for (int i = 0; i < 3; ++i)
            for (int j = 0; j < 3; ++j)
                ImPhi[i][j] = (i == j ? 1.0 : 0.0) - f.Phi[i][j];

        mul3(ATSd, ImPhi, blockI);
        mul3(ATSd, f.Phi, blockJ);
        mulTN3(m_R0, f.Psi_th, Psi);

        B.Zero();
        for (int i = 0; i < 3; ++i) {
            for (int j = 0; j < 3; ++j) {
                B(i, j) = -AT[i][j];
                B(i, 3 + j) = blockI[i][j];
                B(i, 6 + j) = AT[i][j];
                B(i, 9 + j) = blockJ[i][j];
                B(3 + i, 3 + j) = -Psi[i][j];
                B(3 + i, 9 + j) = Psi[i][j];
            }
        }
    }

    // --- helpers on 3x12 blocks -------------------------------------------

    static void rowsFrom(const double bI[3][3], const double bthI[3][3],
                         const double bJ[3][3], const double bthJ[3][3],
                         double R[3][12])
    {
        for (int i = 0; i < 3; ++i)
            for (int j = 0; j < 3; ++j) {
                R[i][j] = bI ? bI[i][j] : 0.0;
                R[i][3 + j] = bthI ? bthI[i][j] : 0.0;
                R[i][6 + j] = bJ ? bJ[i][j] : 0.0;
                R[i][9 + j] = bthJ ? bthJ[i][j] : 0.0;
            }
    }

    /// out = M * R
    static void mulRows(const double M[3][3], const double R[3][12],
                        double out[3][12])
    {
        for (int i = 0; i < 3; ++i)
            for (int j = 0; j < 12; ++j) {
                double s = 0.0;
                for (int k = 0; k < 3; ++k) s += M[i][k] * R[k][j];
                out[i][j] = s;
            }
    }

    static void addRows(double dst[3][12], const double src[3][12], double fac)
    {
        for (int i = 0; i < 3; ++i)
            for (int j = 0; j < 12; ++j)
                dst[i][j] += fac * src[i][j];
    }

    static void copyRows(double dst[3][12], const double src[3][12])
    {
        for (int i = 0; i < 3; ++i)
            for (int j = 0; j < 12; ++j)
                dst[i][j] = src[i][j];
    }

    /**
    The exact linearization of p = B^T fL at fixed fL.  With

        F = A f_t,  m = R0 f_r,  M = R_I T^-T(Th) m,  q = d x F

    the internal force reads

        p = [ -F,  -(I-Phi)^T q - M,  +F,  -Phi^T q + M ]

    and the four rows below are its derivative.  See the oracle for the
    derivation and for the gate that verifies it.
    */
    void addGeometricStiffness(const Frame& f, const Vector& fL, Matrix& LHS) const
    {
        using namespace ASDHingeUtils;

        double I3[3][3], Z3[3][3];
        identity3(I3);
        zero3(Z3);

        double ft[3] = { fL(0), fL(1), fL(2) };
        double fr[3] = { fL(3), fL(4), fL(5) };

        double F[3], m[3], M[3], q[3];
        matVec3(f.A, ft, F);
        matVec3(m_R0, fr, m);

        ASDHingeIso TinvT = ASDHingeIso::Tinv(f.th).transposed();
        double tmp3[3];
        TinvT.apply(f.Th, m, tmp3);
        matVec3(f.RI, tmp3, M);

        cross3(f.d, F, q);

        double SF[3][3], Sd[3][3], SM[3][3], Sq[3][3];
        skew3(F, SF);
        skew3(f.d, Sd);
        skew3(M, SM);
        skew3(q, Sq);

        double mSF[3][3], ImPhi[3][3];
        for (int i = 0; i < 3; ++i)
            for (int j = 0; j < 3; ++j) {
                mSF[i][j] = -SF[i][j];
                ImPhi[i][j] = I3[i][j] - f.Phi[i][j];
            }

        double mI3[3][3];
        for (int i = 0; i < 3; ++i)
            for (int j = 0; j < 3; ++j)
                mI3[i][j] = -I3[i][j];

        double Dd[3][12], Dwe[3][12], DTh[3][12], DthI[3][12];
        double mPsi[3][3];
        for (int i = 0; i < 3; ++i)
            for (int j = 0; j < 3; ++j)
                mPsi[i][j] = -f.Psi_th[i][j];

        rowsFrom(mI3, Z3, I3, Z3, Dd);
        rowsFrom(Z3, ImPhi, Z3, f.Phi, Dwe);
        rowsFrom(Z3, mPsi, Z3, f.Psi_th, DTh);
        rowsFrom(Z3, I3, Z3, Z3, DthI);

        // DF = -skew(F) Dwe
        double DF[3][12];
        mulRows(mSF, Dwe, DF);

        // Dq = -skew(F) Dd + skew(d) DF
        double Dq[3][12], t1[3][12];
        mulRows(mSF, Dd, Dq);
        mulRows(Sd, DF, t1);
        addRows(Dq, t1, 1.0);

        // DM = -skew(M) DthI + R_I dTinvT DTh
        double DM[3][12];
        double mSM[3][3];
        for (int i = 0; i < 3; ++i)
            for (int j = 0; j < 3; ++j)
                mSM[i][j] = -SM[i][j];
        mulRows(mSM, DthI, DM);
        double DTinvT[3][3], P[3][3];
        TinvT.dapply(f.Th, m, DTinvT);
        mul3(f.RI, DTinvT, P);
        mulRows(P, DTh, t1);
        addRows(DM, t1, 1.0);

        // D(Phi^T q) = -skew(Phi^T q) DthI + 1/2 R_I dC^T DTh
        //              + Phi^T ( skew(q) DthI + Dq )
        double s[3], PhiTq[3], PhiT[3][3];
        for (int i = 0; i < 3; ++i)
            for (int j = 0; j < 3; ++j)
                PhiT[i][j] = f.Phi[j][i];
        matTVec3(f.RI, q, s);
        matVec3(PhiT, q, PhiTq);

        double DPhiTq[3][12];
        double mSPhiTq[3][3];
        skew3(PhiTq, mSPhiTq);
        for (int i = 0; i < 3; ++i)
            for (int j = 0; j < 3; ++j)
                mSPhiTq[i][j] = -mSPhiTq[i][j];
        mulRows(mSPhiTq, DthI, DPhiTq);

        ASDHingeIso CT = f.C.transposed();
        double DCT[3][3], RDCT[3][3];
        CT.dapply(f.Th, s, DCT);
        mul3(f.RI, DCT, RDCT);
        mulRows(RDCT, DTh, t1);
        addRows(DPhiTq, t1, 0.5);

        double t2[3][12], t3[3][12];
        mulRows(Sq, DthI, t2);
        addRows(t2, Dq, 1.0);
        mulRows(PhiT, t2, t3);
        addRows(DPhiTq, t3, 1.0);

        // assemble
        for (int j = 0; j < 12; ++j) {
            for (int i = 0; i < 3; ++i) {
                LHS(i, j) += -DF[i][j];
                LHS(6 + i, j) += DF[i][j];
                LHS(3 + i, j) += -(Dq[i][j] - DPhiTq[i][j]) - DM[i][j];
                LHS(9 + i, j) += -DPhiTq[i][j] + DM[i][j];
            }
        }
    }

private:
    double m_qn[2][4];
    double m_rv[2][3];
    double m_qn_cvg[2][4];
    double m_rv_cvg[2][3];
};

#endif // ASDHingeCorotationalTransformation_h
