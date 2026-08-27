/* ****************************************************************** **
**    OpenSees - Open System for Earthquake Engineering Simulation    **
**          Pacific Earthquake Engineering Research Center            **
** ****************************************************************** */

// $Revision: 1.0 $
// $Date: 2026/08/23 $

// Original implementation: Massimo Petracca (ASDEA)
//
// ASDHingeTransformation: the LINEAR transformation of ASDHinge, and the base
// class of ASDHingeCorotationalTransformation.  Same linear/corotational
// virtual pair as the shells (ASDShellQ4Transformation); the solids use a null
// pointer for the linear case instead, which forces an if at every call site.
//
// Conventions
// -----------
// m_R0[i][j] is the i-th GLOBAL component of the j-th LOCAL axis, i.e. the
// columns of R0 are the local axes.  So (R0^T v)_j = sum_i m_R0[i][j] v_i.
//
// The global displacement vector UG (12) handed around is ALWAYS the trial
// displacement minus the captured baseline m_U0.  m_U0 is what makes a hinge
// born inside a construction stage strain free, and for the corotational
// transformation it also re-bases the nodal triads to the identity: no
// separate d0 is needed.

#ifndef ASDHingeTransformation_h
#define ASDHingeTransformation_h

#include <Vector.h>
#include <Matrix.h>
#include <Node.h>
#include <OPS_Globals.h>
#include <cmath>

// --------------------------------------------------------------------------
// small fixed-size linear algebra: no Vector/Matrix here on purpose, they are
// dynamically sized in OpenSees and every temporary is a new/delete
// --------------------------------------------------------------------------

namespace ASDHingeUtils {

inline void zero3(double A[3][3])
{
    for (int i = 0; i < 3; ++i)
        for (int j = 0; j < 3; ++j)
            A[i][j] = 0.0;
}

inline void identity3(double A[3][3])
{
    zero3(A);
    A[0][0] = A[1][1] = A[2][2] = 1.0;
}

inline void copy3(const double A[3][3], double B[3][3])
{
    for (int i = 0; i < 3; ++i)
        for (int j = 0; j < 3; ++j)
            B[i][j] = A[i][j];
}

// C = A * B
inline void mul3(const double A[3][3], const double B[3][3], double C[3][3])
{
    for (int i = 0; i < 3; ++i)
        for (int j = 0; j < 3; ++j) {
            double s = 0.0;
            for (int k = 0; k < 3; ++k) s += A[i][k] * B[k][j];
            C[i][j] = s;
        }
}

// C = A^T * B
inline void mulTN3(const double A[3][3], const double B[3][3], double C[3][3])
{
    for (int i = 0; i < 3; ++i)
        for (int j = 0; j < 3; ++j) {
            double s = 0.0;
            for (int k = 0; k < 3; ++k) s += A[k][i] * B[k][j];
            C[i][j] = s;
        }
}

// r = A * v
inline void matVec3(const double A[3][3], const double* v, double* r)
{
    for (int i = 0; i < 3; ++i)
        r[i] = A[i][0] * v[0] + A[i][1] * v[1] + A[i][2] * v[2];
}

// r = A^T * v
inline void matTVec3(const double A[3][3], const double* v, double* r)
{
    for (int i = 0; i < 3; ++i)
        r[i] = A[0][i] * v[0] + A[1][i] * v[1] + A[2][i] * v[2];
}

inline void cross3(const double* a, const double* b, double* r)
{
    r[0] = a[1] * b[2] - a[2] * b[1];
    r[1] = a[2] * b[0] - a[0] * b[2];
    r[2] = a[0] * b[1] - a[1] * b[0];
}

inline double dot3(const double* a, const double* b)
{
    return a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
}

inline double norm3(const double* a)
{
    return std::sqrt(dot3(a, a));
}

inline void skew3(const double* v, double S[3][3])
{
    S[0][0] = 0.0;   S[0][1] = -v[2]; S[0][2] = v[1];
    S[1][0] = v[2];  S[1][1] = 0.0;   S[1][2] = -v[0];
    S[2][0] = -v[1]; S[2][1] = v[0];  S[2][2] = 0.0;
}

// Rodrigues
inline void expmap3(const double* th, double R[3][3])
{
    double a = norm3(th);
    double c1, c2;
    if (a < 1.0e-8) {
        c1 = 1.0 - a * a / 6.0;
        c2 = 0.5 - a * a / 24.0;
    }
    else {
        c1 = std::sin(a) / a;
        c2 = (1.0 - std::cos(a)) / (a * a);
    }
    double X[3][3], X2[3][3];
    skew3(th, X);
    mul3(X, X, X2);
    for (int i = 0; i < 3; ++i)
        for (int j = 0; j < 3; ++j)
            R[i][j] = (i == j ? 1.0 : 0.0) + c1 * X[i][j] + c2 * X2[i][j];
}

// rotation vector of R, through the quaternion (robust up to |th| -> pi)
inline void logmap3(const double R[3][3], double* th)
{
    double tr = R[0][0] + R[1][1] + R[2][2];
    double xyz[3], w;
    if (tr > -0.99) {
        double t = 1.0 + tr;
        w = 0.5 * std::sqrt(t > 0.0 ? t : 0.0);
        double f = 0.25 / w;
        xyz[0] = (R[2][1] - R[1][2]) * f;
        xyz[1] = (R[0][2] - R[2][0]) * f;
        xyz[2] = (R[1][0] - R[0][1]) * f;
    }
    else {
        int i = 0;
        if (R[1][1] > R[i][i]) i = 1;
        if (R[2][2] > R[i][i]) i = 2;
        int j = (i + 1) % 3;
        int k = (i + 2) % 3;
        double t = 1.0 + R[i][i] - R[j][j] - R[k][k];
        double s = std::sqrt(t > 0.0 ? t : 0.0);
        xyz[i] = 0.5 * s;
        xyz[j] = (R[j][i] + R[i][j]) / (2.0 * s);
        xyz[k] = (R[k][i] + R[i][k]) / (2.0 * s);
        w = (R[k][j] - R[j][k]) / (2.0 * s);
    }
    double n = norm3(xyz);
    if (n < 1.0e-14) {
        th[0] = th[1] = th[2] = 0.0;
        return;
    }
    double a = 2.0 * std::atan2(n, std::fabs(w));
    double sg = (w >= 0.0) ? 1.0 : -1.0;
    double f = sg * a / n;
    th[0] = xyz[0] * f;
    th[1] = xyz[1] * f;
    th[2] = xyz[2] * f;
}

} // namespace ASDHingeUtils


/**
Working storage shared by every ASDHinge, on the model of ASDShellQ4Globals.
Vector and Matrix are dynamically sized in OpenSees: a local temporary is a
new/delete on every call, and getTangentStiff is called once per element per
iteration.
*/
class ASDHingeGlobals
{
private:
    ASDHingeGlobals() = default;

public:
    Vector UG = Vector(12);   // trial displacements minus the baseline
    Vector VG = Vector(12);   // trial velocities
    Vector e = Vector(6);     // local deformations
    Vector edot = Vector(6);  // local deformation rates
    Matrix B = Matrix(6, 12); // de/dUG
    Matrix BtK = Matrix(12, 6);

    static ASDHingeGlobals& instance()
    {
        static ASDHingeGlobals x;
        return x;
    }
};


class ASDHingeTransformation
{
public:
    ASDHingeTransformation()
        : m_U0(12)
    {
        ASDHingeUtils::identity3(m_R0);
        m_nodes[0] = 0;
        m_nodes[1] = 0;
        m_U0.Zero();
    }

    virtual ~ASDHingeTransformation() {}

    virtual bool isLinear() const { return true; }

    // -- reference orientation --------------------------------------------

    /// columns of R0 are the local axes, in global components
    void setOrientation(const double R0[3][3])
    {
        ASDHingeUtils::copy3(R0, m_R0);
    }

    /// local axis i (0,1,2) in global components
    void getLocalAxis(int i, double* v) const
    {
        v[0] = m_R0[0][i];
        v[1] = m_R0[1][i];
        v[2] = m_R0[2][i];
    }

    // -- life cycle --------------------------------------------------------

    virtual void setDomain(Node** nodes, bool initialized)
    {
        m_nodes[0] = nodes[0];
        m_nodes[1] = nodes[1];
        if (!initialized)
            this->forceCaptureInitialDisp();
    }

    virtual void revertToStart()
    {
        m_U0.Zero();
    }

    virtual void commit() {}
    virtual void revertToLastCommit() {}

    /// only the corotational transformation has anything to do here
    virtual void update(const Vector& UG) { (void)UG; }

    /// re-baseline at the current configuration: a hinge activated inside a
    /// construction stage must be born with zero deformation
    virtual void forceCaptureInitialDisp()
    {
        if (m_nodes[0] == 0 || m_nodes[1] == 0)
            return;
        const Vector& d0 = m_nodes[0]->getTrialDisp();
        const Vector& d1 = m_nodes[1]->getTrialDisp();
        for (int i = 0; i < 6; ++i) {
            m_U0(i) = (i < d0.Size()) ? d0(i) : 0.0;
            m_U0(6 + i) = (i < d1.Size()) ? d1(i) : 0.0;
        }
    }

    // -- kinematics --------------------------------------------------------

    void computeGlobalDisplacements(Vector& UG) const
    {
        const Vector& d0 = m_nodes[0]->getTrialDisp();
        const Vector& d1 = m_nodes[1]->getTrialDisp();
        for (int i = 0; i < 6; ++i) {
            UG(i) = d0(i) - m_U0(i);
            UG(6 + i) = d1(i) - m_U0(6 + i);
        }
    }

    void computeGlobalVelocities(Vector& VG) const
    {
        const Vector& v0 = m_nodes[0]->getTrialVel();
        const Vector& v1 = m_nodes[1]->getTrialVel();
        for (int i = 0; i < 6; ++i) {
            VG(i) = v0(i);
            VG(6 + i) = v1(i);
        }
    }

    /// e (6) from UG (12)
    virtual void computeDeformations(const Vector& UG, Vector& e)
    {
        for (int j = 0; j < 3; ++j) {
            double et = 0.0, er = 0.0;
            for (int i = 0; i < 3; ++i) {
                et += m_R0[i][j] * (UG(6 + i) - UG(i));
                er += m_R0[i][j] * (UG(9 + i) - UG(3 + i));
            }
            e(j) = et;
            e(3 + j) = er;
        }
    }

    /// B in the REFERENCE configuration.  Also what getInitialStiff uses, in
    /// the corotational case as well.
    void computeLinearB(Matrix& B) const
    {
        B.Zero();
        for (int j = 0; j < 3; ++j) {
            for (int i = 0; i < 3; ++i) {
                B(j, i) = -m_R0[i][j];
                B(j, 6 + i) = m_R0[i][j];
                B(3 + j, 3 + i) = -m_R0[i][j];
                B(3 + j, 9 + i) = m_R0[i][j];
            }
        }
    }

    /// B = de/dUG (6x12)
    virtual void computeB(const Vector& UG, Matrix& B)
    {
        (void)UG;
        this->computeLinearB(B);
    }

    /**
    RHS = B^T fL, LHS = B^T kL B (+ the geometric part, in the derived class).
    RHS is the INTERNAL force, positive.
    */
    virtual void transformToGlobal(const Vector& UG, const Vector& fL,
                                   const Matrix& kL, Matrix& LHS, Vector& RHS,
                                   bool LHSrequired)
    {
        Matrix& B = ASDHingeGlobals::instance().B;
        this->computeB(UG, B);
        RHS.addMatrixTransposeVector(0.0, B, fL, 1.0);
        if (LHSrequired) {
            Matrix& BtK = ASDHingeGlobals::instance().BtK;
            BtK.addMatrixTransposeProduct(0.0, B, kL, 1.0);
            LHS.addMatrixProduct(0.0, BtK, B, 1.0);
        }
    }

    // -- serialization -----------------------------------------------------

    virtual int internalDataSize() const { return 12; }

    virtual void saveInternalData(Vector& v, int pos) const
    {
        for (int i = 0; i < 12; ++i)
            v(pos + i) = m_U0(i);
    }

    virtual void restoreInternalData(const Vector& v, int pos)
    {
        for (int i = 0; i < 12; ++i)
            m_U0(i) = v(pos + i);
    }

protected:
    double m_R0[3][3];
    Node* m_nodes[2];
    Vector m_U0;
};

#endif // ASDHingeTransformation_h
