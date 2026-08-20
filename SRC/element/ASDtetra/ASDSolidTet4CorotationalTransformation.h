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
// Corotational (EICR, translational-only) coordinate transformation for
// the 4-node tetrahedron. This is the constant-F special case of
// ASDSolidHexCorotationalTransformation:
//
// - the deformation gradient F = sum_a (x_a - c) (x) g_a is CONSTANT over
//   the element (g_a are the constant reference gradients, sum g_a = 0),
//   so ONE polar decomposition R = polar(F) defines the frame exactly:
//   R = I under pure stretch (patch-test exact) and R equals the rigid
//   rotation under rigid motion;
// - the rotation gradient has the closed form (in the corotated frame)
//   Gl_b = (tr(U) I - U)^{-1} skew(g_b), U = R^T F;
// - the tangent is the three-term Felippa form of the hexa transformation
//   (transformToGlobal): P = Pu - S*Gl,
//   f_g = T^T P^T f_L,
//   K_g = T^T (P^T K_L P - Gl^T Fnm^T P - Fnm Gl) T,
//   which for this translational EICR with the exact G was MEASURED to be
//   the EXACT Hessian of the corotational energy (mismatch vs FD at the
//   noise floor at any state) - see the red gate in
//   OpenSees-Testing/new-asd-elements/ASDTet/verify_asdtet_corot.py.
//
// The class is stateless (the frame is recomputed from scratch at every
// call, no incremental data): nothing to commit, revert or serialize.
// The input displacement vector is ALWAYS the displacement since
// activation (U - U0, handled by the element).

#ifndef ASDSolidTet4CorotationalTransformation_h
#define ASDSolidTet4CorotationalTransformation_h

#include <Node.h>
#include <Vector.h>
#include <Matrix.h>
#include <array>

class ASDSolidTet4CorotationalTransformation
{

public:

    ASDSolidTet4CorotationalTransformation() = default;

    // stores the node pointers and the constant reference gradients
    // (rows of g: dN_a/dX, computed by the element, sum = 0)
    void setReference(const std::array<Node*, 4>& nodes, const Matrix& g)
    {
        m_nodes = nodes;
        for (int a = 0; a < 4; ++a)
            for (int i = 0; i < 3; ++i)
                m_g[a][i] = g(a, i);
        for (int i = 0; i < 3; ++i) {
            double s = 0.0;
            for (int a = 0; a < 4; ++a)
                s += m_nodes[a]->getCrds()(i);
            m_c0[i] = 0.25 * s;
        }
    }

    // deformational displacements: uL_a = R^T (x_a - c) - (X_a - c0),
    // with x = X + d (d = displacements since activation)
    void calculateLocalDisplacements(const Vector& d, Vector& UL)
    {
        double x[4][3], y[4][3], R[3][3];
        computeFrame(d, x, R, y);
        for (int a = 0; a < 4; ++a) {
            const Vector& Xa = m_nodes[a]->getCrds();
            for (int i = 0; i < 3; ++i)
                UL(3 * a + i) = y[a][i] - (Xa(i) - m_c0[i]);
        }
    }

    // Felippa consistent transformation. On input: RHS = +f_L (local
    // internal force), LHS = K_L. On output both are global.
    void transformToGlobal(const Vector& d, Matrix& LHS, Vector& RHS, bool LHSrequired)
    {
        double x[4][3], y[4][3], R[3][3];
        computeFrame(d, x, R, y);

        // Gl (3x12): dw = Gl * dv_local, Gl_b = Minv * skew(g_b)
        // (in the corotated frame the R^T factors cancel)
        double F[3][3];
        deformationGradient(x, F);
        double U[3][3];
        for (int i = 0; i < 3; ++i)
            for (int j = 0; j < 3; ++j) {
                double s = 0.0;
                for (int k = 0; k < 3; ++k)
                    s += R[k][i] * F[k][j];
                U[i][j] = s;
            }
        double trU = U[0][0] + U[1][1] + U[2][2];
        double M[3][3], Minv[3][3];
        for (int i = 0; i < 3; ++i)
            for (int j = 0; j < 3; ++j)
                M[i][j] = (i == j ? trU : 0.0) - U[i][j];
        inv3(M, Minv);
        static Matrix Gl(3, 12);
        for (int b = 0; b < 4; ++b) {
            double Sk[3][3];
            skew3(m_g[b], Sk);
            for (int i = 0; i < 3; ++i)
                for (int j = 0; j < 3; ++j) {
                    double s = 0.0;
                    for (int k = 0; k < 3; ++k)
                        s += Minv[i][k] * Sk[k][j];
                    Gl(i, 3 * b + j) = s;
                }
        }

        // P = Pu - S*Gl, with Pu the translational projector and
        // S_a = -skew(y_a)
        static Matrix P(12, 12);
        P.Zero();
        for (int a = 0; a < 4; ++a)
            for (int b = 0; b < 4; ++b) {
                double blk = (a == b ? 1.0 : 0.0) - 0.25;
                for (int i = 0; i < 3; ++i)
                    P(3 * a + i, 3 * b + i) = blk;
            }
        static Matrix S(12, 3);
        S.Zero();
        for (int a = 0; a < 4; ++a) {
            double Sy[3][3];
            skew3(y[a], Sy);
            for (int i = 0; i < 3; ++i)
                for (int j = 0; j < 3; ++j)
                    S(3 * a + i, j) = -Sy[i][j];
        }
        static Matrix SG(12, 12);
        SG.addMatrixProduct(0.0, S, Gl, 1.0);
        P.addMatrix(1.0, SG, -1.0);

        // projected local forces: pe = P^T f_L
        static Vector pe(12);
        pe.addMatrixTransposeVector(0.0, P, RHS, 1.0);

        // global RHS = T^T pe, with T = blockdiag(R^T):
        // (T^T pe)_a = R * pe_a
        for (int a = 0; a < 4; ++a)
            for (int i = 0; i < 3; ++i) {
                double s = 0.0;
                for (int k = 0; k < 3; ++k)
                    s += R[i][k] * pe(3 * a + k);
                RHS(3 * a + i) = s;
            }

        if (!LHSrequired)
            return;

        // K_local = P^T K_L P - Gl^T Fnm^T P - Fnm Gl,
        // Fnm = spins of the projected forces
        static Matrix tmp(12, 12);
        tmp.addMatrixProduct(0.0, LHS, P, 1.0);
        LHS.addMatrixTransposeProduct(0.0, P, tmp, 1.0);

        static Matrix Fnm(12, 3);
        Fnm.Zero();
        for (int a = 0; a < 4; ++a) {
            double fa[3] = { pe(3 * a), pe(3 * a + 1), pe(3 * a + 2) };
            double Sf[3][3];
            skew3(fa, Sf);
            for (int i = 0; i < 3; ++i)
                for (int j = 0; j < 3; ++j)
                    Fnm(3 * a + i, j) = Sf[i][j];
        }
        static Matrix GtFt(12, 12);
        static Matrix FnmT(3, 12);
        FnmT.addMatrixTranspose(0.0, Fnm, 1.0);
        GtFt.addMatrixTransposeProduct(0.0, Gl, FnmT, 1.0);
        LHS.addMatrixProduct(1.0, GtFt, P, -1.0);
        LHS.addMatrixProduct(1.0, Fnm, Gl, -1.0);

        // K_global = T^T K_local T: block (a,b) -> R * K_ab * R^T
        static Matrix KG(12, 12);
        for (int a = 0; a < 4; ++a)
            for (int b = 0; b < 4; ++b)
                for (int i = 0; i < 3; ++i)
                    for (int j = 0; j < 3; ++j) {
                        double s = 0.0;
                        for (int k = 0; k < 3; ++k)
                            for (int l = 0; l < 3; ++l)
                                s += R[i][k] * LHS(3 * a + k, 3 * b + l) * R[j][l];
                        KG(3 * a + i, 3 * b + j) = s;
                    }
        LHS = KG;
    }

private:

    // F = sum_a (x_a - c) (x) g_a
    void deformationGradient(const double x[4][3], double F[3][3]) const
    {
        double c[3] = { 0.0, 0.0, 0.0 };
        for (int a = 0; a < 4; ++a)
            for (int i = 0; i < 3; ++i)
                c[i] += 0.25 * x[a][i];
        for (int i = 0; i < 3; ++i)
            for (int j = 0; j < 3; ++j) {
                double s = 0.0;
                for (int a = 0; a < 4; ++a)
                    s += (x[a][i] - c[i]) * m_g[a][j];
                F[i][j] = s;
            }
    }

    // current positions, frame R = polar(F), corotated levers y_a = R^T (x_a - c)
    void computeFrame(const Vector& d, double x[4][3], double R[3][3], double y[4][3]) const
    {
        for (int a = 0; a < 4; ++a) {
            const Vector& Xa = m_nodes[a]->getCrds();
            for (int i = 0; i < 3; ++i)
                x[a][i] = Xa(i) + d(3 * a + i);
        }
        double F[3][3];
        deformationGradient(x, F);
        polar3(F, R);
        double c[3] = { 0.0, 0.0, 0.0 };
        for (int a = 0; a < 4; ++a)
            for (int i = 0; i < 3; ++i)
                c[i] += 0.25 * x[a][i];
        for (int a = 0; a < 4; ++a)
            for (int i = 0; i < 3; ++i) {
                double s = 0.0;
                for (int k = 0; k < 3; ++k)
                    s += R[k][i] * (x[a][k] - c[k]);
                y[a][i] = s;
            }
    }

    static void skew3(const double v[3], double S[3][3])
    {
        S[0][0] = 0.0;   S[0][1] = -v[2]; S[0][2] = v[1];
        S[1][0] = v[2];  S[1][1] = 0.0;   S[1][2] = -v[0];
        S[2][0] = -v[1]; S[2][1] = v[0];  S[2][2] = 0.0;
    }

    static void inv3(const double A[3][3], double B[3][3])
    {
        double det =
            A[0][0] * (A[1][1] * A[2][2] - A[1][2] * A[2][1]) -
            A[0][1] * (A[1][0] * A[2][2] - A[1][2] * A[2][0]) +
            A[0][2] * (A[1][0] * A[2][1] - A[1][1] * A[2][0]);
        double id = 1.0 / det;
        B[0][0] = (A[1][1] * A[2][2] - A[1][2] * A[2][1]) * id;
        B[0][1] = (A[0][2] * A[2][1] - A[0][1] * A[2][2]) * id;
        B[0][2] = (A[0][1] * A[1][2] - A[0][2] * A[1][1]) * id;
        B[1][0] = (A[1][2] * A[2][0] - A[1][0] * A[2][2]) * id;
        B[1][1] = (A[0][0] * A[2][2] - A[0][2] * A[2][0]) * id;
        B[1][2] = (A[0][2] * A[1][0] - A[0][0] * A[1][2]) * id;
        B[2][0] = (A[1][0] * A[2][1] - A[1][1] * A[2][0]) * id;
        B[2][1] = (A[0][1] * A[2][0] - A[0][0] * A[2][1]) * id;
        B[2][2] = (A[0][0] * A[1][1] - A[0][1] * A[1][0]) * id;
    }

    // Higham polar iteration, R <- (R + R^-T)/2. F is a small perturbation
    // of a rotation times a stretch close to identity in our use, so the
    // iteration converges in a handful of steps.
    static void polar3(const double F[3][3], double R[3][3])
    {
        for (int i = 0; i < 3; ++i)
            for (int j = 0; j < 3; ++j)
                R[i][j] = F[i][j];
        for (int iter = 0; iter < 100; ++iter) {
            double Rin[3][3], Rit[3][3];
            inv3(R, Rin);
            for (int i = 0; i < 3; ++i)
                for (int j = 0; j < 3; ++j)
                    Rit[i][j] = Rin[j][i];
            double err = 0.0;
            for (int i = 0; i < 3; ++i)
                for (int j = 0; j < 3; ++j) {
                    double next = 0.5 * (R[i][j] + Rit[i][j]);
                    err += (next - R[i][j]) * (next - R[i][j]);
                    R[i][j] = next;
                }
            if (err < 1.0e-28)
                break;
        }
    }

private:

    std::array<Node*, 4> m_nodes = { nullptr, nullptr, nullptr, nullptr };
    double m_g[4][3] = {};
    double m_c0[3] = {};
};

#endif // ASDSolidTet4CorotationalTransformation_h
