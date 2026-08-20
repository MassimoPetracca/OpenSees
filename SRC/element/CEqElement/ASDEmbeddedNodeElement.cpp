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

#include <ASDEmbeddedNodeElement.h>

#include <Domain.h>
#include <Node.h>
#include <ErrorHandler.h>
#include <Channel.h>
#include <FEM_ObjectBroker.h>
#include <elementAPI.h>
#include <Renderer.h>
#include <analysis/dof_grp/DOF_Group.h>
#include <UniaxialMaterial.h>
#include <Information.h>
#include <ElementResponse.h>

#include <ASDMath.h>

#include <stdio.h>
#include <stdlib.h>
#include <algorithm>
#include <cmath>
#include <string>
#include <limits>

// anonymous namespace for utilities
namespace
{

    double det2(const Matrix& J) {
        return J(0, 0) * J(1, 1) - J(0, 1) * J(1, 0);
    }

    double det3(const Matrix& J) {
        return  J(0, 0) * J(1, 1) * J(2, 2) - J(0, 0) * J(1, 2) * J(2, 1) - J(0, 1) * J(1, 0) * J(2, 2) +
            J(0, 1) * J(1, 2) * J(2, 0) + J(0, 2) * J(1, 0) * J(2, 1) - J(0, 2) * J(1, 1) * J(2, 0);
    }

    void cross(const Vector& a, const Vector& b, Vector& c) {
        c(0) = a(1) * b(2) - a(2) * b(1);
        c(1) = a(2) * b(0) - a(0) * b(2);
        c(2) = a(0) * b(1) - a(1) * b(0);
    }

    // 3x3 inverse (no pivoting: callers pass well-conditioned matrices)
    void inv3(const Matrix& A, Matrix& B) {
        double d = det3(A);
        B(0, 0) = (A(1, 1) * A(2, 2) - A(1, 2) * A(2, 1)) / d;
        B(0, 1) = (A(0, 2) * A(2, 1) - A(0, 1) * A(2, 2)) / d;
        B(0, 2) = (A(0, 1) * A(1, 2) - A(0, 2) * A(1, 1)) / d;
        B(1, 0) = (A(1, 2) * A(2, 0) - A(1, 0) * A(2, 2)) / d;
        B(1, 1) = (A(0, 0) * A(2, 2) - A(0, 2) * A(2, 0)) / d;
        B(1, 2) = (A(0, 2) * A(1, 0) - A(0, 0) * A(1, 2)) / d;
        B(2, 0) = (A(1, 0) * A(2, 1) - A(1, 1) * A(2, 0)) / d;
        B(2, 1) = (A(0, 1) * A(2, 0) - A(0, 0) * A(2, 1)) / d;
        B(2, 2) = (A(0, 0) * A(1, 1) - A(0, 1) * A(1, 0)) / d;
    }

    // R = polar(F) by Higham's iteration R <- (R + R^-T)/2 (same rule the
    // ASDhexa corotational frame uses)
    void polar3(const Matrix& F, Matrix& R) {
        static Matrix Ri(3, 3);
        static Matrix Rn(3, 3);
        R = F;
        for (int it = 0; it < 100; ++it) {
            inv3(R, Ri);
            double diff = 0.0;
            for (int i = 0; i < 3; ++i)
                for (int j = 0; j < 3; ++j) {
                    Rn(i, j) = 0.5 * (R(i, j) + Ri(j, i));
                    double dd = Rn(i, j) - R(i, j);
                    diff += dd * dd;
                }
            R = Rn;
            if (std::sqrt(diff) < 1.0e-14)
                break;
        }
    }

    // inverse of the left SO(3) tangent map: d(log R) = Tinv * eta for
    // dR = skew(eta) R
    void leftJacobianInv(const Vector& v, Matrix& T) {
        double t = v.Norm();
        static Matrix K(3, 3);
        K.Zero();
        K(0, 1) = -v(2); K(0, 2) = v(1);
        K(1, 0) = v(2);  K(1, 2) = -v(0);
        K(2, 0) = -v(1); K(2, 1) = v(0);
        T.Zero();
        for (int i = 0; i < 3; ++i) T(i, i) = 1.0;
        T.addMatrix(1.0, K, -0.5);
        if (t >= 1.0e-12) {
            double half = 0.5 * t;
            double coeff = (1.0 - half / std::tan(half)) / (t * t);
            static Matrix K2(3, 3);
            K2.addMatrixProduct(0.0, K, K, 1.0);
            T.addMatrix(1.0, K2, coeff);
        }
    }

    namespace tri {

        double shapeFun(double x, double y, int i) {
            if (i == 0)
                return 1.0 - x - y;
            else if (i == 1)
                return x;
            else if (i == 2)
                return y;
            return 0.0;
        }

        void shapeFunDer(Matrix& dN) {
            dN(0, 0) = -1.0; dN(0, 1) = -1.0;
            dN(1, 0) = 1.0; dN(1, 1) = 0.0;
            dN(2, 0) = 0.0; dN(2, 1) = 1.0;
        }

        void globalCoord(const Matrix& X, double lx, double ly, double& gx, double& gy) {
            gx = gy = 0.0;
            for (int i = 0; i < 3; i++) {
                double N = shapeFun(lx, ly, i);
                gx += N * X(0, i);
                gy += N * X(1, i);
            }
        }

        void globalCoord(const Matrix& X, double lx, double ly, double& gx, double& gy, double& gz) {
            gx = gy = gz = 0.0;
            for (int i = 0; i < 3; i++) {
                double N = shapeFun(lx, ly, i);
                gx += N * X(0, i);
                gy += N * X(1, i);
                gz += N * X(2, i);
            }
        }

        void localCoord(const Matrix& X, const Matrix& invJ, double gx, double gy, double& lx, double& ly) {
            lx = ly = 0.0;
            double px, py;
            globalCoord(X, lx, ly, px, py);
            Vector D(2);
            Vector DL(2);
            D(0) = gx - px;
            D(1) = gy - py;
            DL.addMatrixVector(0.0, invJ, D, 1.0);
            lx = DL(0);
            ly = DL(1);
        }

        void localCoord(const Matrix& X, const Matrix& invJ, double gx, double gy, double gz, double& lx, double& ly) {
            lx = ly = 0.0;
            double px, py, pz;
            globalCoord(X, lx, ly, px, py, pz);
            Vector D(3);
            Vector DL(3);
            D(0) = gx - px;
            D(1) = gy - py;
            D(2) = gz - pz;
            DL.addMatrixVector(0.0, invJ, D, 1.0);
            lx = DL(0);
            ly = DL(1);
        }

        // A triangle in 3D has two natural coordinates but three equations, so
        // J = [a_xi, a_eta, 0] is singular. Filling the third column with the
        // unit normal makes it square and invertible, and gives the inverse a
        // useful meaning: the rows of J^-1 are the dual basis (a^xi, a^eta, n),
        // so J^-1 * D returns the offset D resolved on the triangle plane plus
        // its normal component. Reading only the first two entries is therefore
        // already the orthogonal projection onto the plane, and it works for a
        // point off the plane too.
        //
        // This is why the third COLUMN of the inverse must not be zeroed: doing
        // so drops the contribution of D_z, which is harmless only when the
        // normal is parallel to the global z axis. On any other orientation it
        // silently returns the wrong shape functions -- 25% error on the test
        // case tri3d_u_incl, which is what a shell face in a real model looks
        // like. The zeroing was present in TRI_3D_U and TRI_3D_UP and no test
        // covered a 3D triangle, so it survived from the first commit.
        void fillVzInJacobian(Matrix& J) {
            double nx = J(1, 0) * J(2, 1) - J(1, 1) * J(2, 0);
            double ny = J(0, 1) * J(2, 0) - J(0, 0) * J(2, 1);
            double nz = J(0, 0) * J(1, 1) - J(0, 1) * J(1, 0);
            double norm = std::sqrt(nx * nx + ny * ny + nz * nz);
            if (norm > std::numeric_limits<double>::epsilon()) {
                J(0, 2) = nx / norm;
                J(1, 2) = ny / norm;
                J(2, 2) = nz / norm;
            }
        }

    }

    // ------------------------------------------------------------------
    // 4-node quadrilateral (bilinear) host.
    // Unlike the simplices, the isoparametric map is NOT affine: the natural
    // coordinate of the constrained node needs a Newton iteration.
    // Node ordering (the OpenSees / MpcCore standard, counterclockwise):
    //   1(-1,-1)  2(+1,-1)  3(+1,+1)  4(-1,+1)
    // ------------------------------------------------------------------
    namespace quad {

        static const double NX[4] = { -1.0,  1.0, 1.0, -1.0 };
        static const double NY[4] = { -1.0, -1.0, 1.0,  1.0 };

        void shapeFun(double x, double y, Vector& N) {
            for (int i = 0; i < 4; ++i)
                N(i) = 0.25 * (1.0 + NX[i] * x) * (1.0 + NY[i] * y);
        }

        void shapeFunDer(double x, double y, Matrix& dN) {
            for (int i = 0; i < 4; ++i) {
                dN(i, 0) = 0.25 * NX[i] * (1.0 + NY[i] * y);
                dN(i, 1) = 0.25 * NY[i] * (1.0 + NX[i] * x);
            }
        }

        // Newton inversion of X(xi) = P, starting from the element centre.
        // The step is measured in natural coordinates, which are dimensionless
        // and O(1), so an absolute tolerance is meaningful (a tolerance on the
        // physical residual would instead scale with the element size).
        bool localCoord(const Matrix& X, double gx, double gy, double& lx, double& ly) {
            static Vector N(4);
            static Matrix dN(4, 2);
            static Matrix J(2, 2);
            static Matrix invJ(2, 2);
            double x = 0.0, y = 0.0;
            for (int iter = 0; iter < 20; ++iter) {
                shapeFun(x, y, N);
                double px = 0.0, py = 0.0;
                for (int i = 0; i < 4; ++i) {
                    px += N(i) * X(0, i);
                    py += N(i) * X(1, i);
                }
                shapeFunDer(x, y, dN);
                J.addMatrixProduct(0.0, X, dN, 1.0);
                if (J.Invert(invJ) < 0)
                    return false;
                double rx = gx - px;
                double ry = gy - py;
                double dx = invJ(0, 0) * rx + invJ(0, 1) * ry;
                double dy = invJ(1, 0) * rx + invJ(1, 1) * ry;
                x += dx;
                y += dy;
                if (std::sqrt(dx * dx + dy * dy) < 1.0e-12) {
                    lx = x; ly = y;
                    return true;
                }
            }
            lx = x; ly = y;
            return false;
        }

        // area by 2x2 Gauss (exact: detJ is affine in the natural coordinates)
        double area(const Matrix& X) {
            static const double g = 0.577350269189625764509;
            static const double p[2] = { -g, g };
            static Matrix dN(4, 2);
            static Matrix J(2, 2);
            double A = 0.0;
            for (int a = 0; a < 2; ++a) {
                for (int b = 0; b < 2; ++b) {
                    shapeFunDer(p[a], p[b], dN);
                    J.addMatrixProduct(0.0, X, dN, 1.0);
                    A += det2(J);
                }
            }
            return A;
        }

        // Local frame of a quadrilateral face in 3D, and the face coordinates
        // resolved on it.
        //
        // The frame is the one ASDShellQ4LocalCoordinateSystem builds, on
        // purpose: the constrained node is normally tied to a shell, and the
        // shape functions that must be reproduced are the ones the shell itself
        // uses. e3 comes from the cross product of the two diagonals, so it is
        // the normal of the MEAN plane of a warped face, not of any three of
        // its nodes; e1 is side 1-2 projected onto that plane.
        //
        // The four nodes are then projected onto the mean plane. This costs
        // nothing in accuracy for the translational constraint, even on a
        // warped face: the projection is affine and the bilinear shape
        // functions sum to one, so projecting commutes with interpolating,
        // Pi(sum N_i X_i) = sum N_i Pi(X_i). Inverting the planar map at the
        // projected point therefore returns the EXACT natural coordinate of a
        // node lying on the warped bilinear surface.
        //
        // What the projection does fix is a convention: for a node off the
        // surface it is the direction along which the node is attached, and for
        // the rotational mode it is the plane the bending and drilling
        // rotations refer to. Both are the choices the shell makes, which is
        // the point of borrowing its frame. The out-of-plane offsets are
        // returned in `warp` so the caller can warn when the face is warped
        // enough for that convention to be worth knowing about.
        //
        // R maps global to local: u_local = R * u_global, with the rows of R
        // holding e1, e2, e3 -- the same convention as TRI_3D_UR.
        void frame3D(const Matrix& X, Matrix& R, Vector& center,
                     Matrix& XL, double& warp)
        {
            static Vector e1(3), e2(3), e3(3), d13(3), d24(3);
            for (int i = 0; i < 3; ++i) {
                center(i) = 0.25 * (X(i, 0) + X(i, 1) + X(i, 2) + X(i, 3));
                d13(i) = X(i, 2) - X(i, 0);
                d24(i) = X(i, 3) - X(i, 1);
            }
            cross(d13, d24, e3);
            e3.Normalize();
            for (int i = 0; i < 3; ++i)
                e1(i) = X(i, 1) - X(i, 0);
            double e1_dot_e3 = e1 ^ e3;
            for (int i = 0; i < 3; ++i)
                e1(i) -= e1_dot_e3 * e3(i);
            e1.Normalize();
            cross(e3, e1, e2);
            e2.Normalize();
            for (int i = 0; i < 3; ++i) {
                R(0, i) = e1(i);
                R(1, i) = e2(i);
                R(2, i) = e3(i);
            }
            // face coordinates on the mean plane, measured from the centre so
            // that the inverse map does not carry the model's absolute
            // coordinates into the Newton residual
            warp = 0.0;
            for (int j = 0; j < 4; ++j) {
                double dx = X(0, j) - center(0);
                double dy = X(1, j) - center(1);
                double dz = X(2, j) - center(2);
                XL(0, j) = R(0, 0) * dx + R(0, 1) * dy + R(0, 2) * dz;
                XL(1, j) = R(1, 0) * dx + R(1, 1) * dy + R(1, 2) * dz;
                double off = R(2, 0) * dx + R(2, 1) * dy + R(2, 2) * dz;
                if (std::abs(off) > warp)
                    warp = std::abs(off);
            }
        }

    }

    // ------------------------------------------------------------------
    // 8-node hexahedron (trilinear) host. Same story as the quadrilateral,
    // with three natural coordinates. Node ordering (OpenSees stdBrick):
    //   bottom face 1..4 counterclockwise, top face 5..8 above them.
    // ------------------------------------------------------------------
    namespace hexa {

        static const double NX[8] = { -1.0,  1.0, 1.0, -1.0, -1.0,  1.0, 1.0, -1.0 };
        static const double NY[8] = { -1.0, -1.0, 1.0,  1.0, -1.0, -1.0, 1.0,  1.0 };
        static const double NZ[8] = { -1.0, -1.0,-1.0, -1.0,  1.0,  1.0, 1.0,  1.0 };

        void shapeFun(double x, double y, double z, Vector& N) {
            for (int i = 0; i < 8; ++i)
                N(i) = 0.125 * (1.0 + NX[i] * x) * (1.0 + NY[i] * y) * (1.0 + NZ[i] * z);
        }

        void shapeFunDer(double x, double y, double z, Matrix& dN) {
            for (int i = 0; i < 8; ++i) {
                dN(i, 0) = 0.125 * NX[i] * (1.0 + NY[i] * y) * (1.0 + NZ[i] * z);
                dN(i, 1) = 0.125 * NY[i] * (1.0 + NX[i] * x) * (1.0 + NZ[i] * z);
                dN(i, 2) = 0.125 * NZ[i] * (1.0 + NX[i] * x) * (1.0 + NY[i] * y);
            }
        }

        bool localCoord(const Matrix& X, double gx, double gy, double gz,
                        double& lx, double& ly, double& lz) {
            static Vector N(8);
            static Matrix dN(8, 3);
            static Matrix J(3, 3);
            static Matrix invJ(3, 3);
            double x = 0.0, y = 0.0, z = 0.0;
            for (int iter = 0; iter < 20; ++iter) {
                shapeFun(x, y, z, N);
                double px = 0.0, py = 0.0, pz = 0.0;
                for (int i = 0; i < 8; ++i) {
                    px += N(i) * X(0, i);
                    py += N(i) * X(1, i);
                    pz += N(i) * X(2, i);
                }
                shapeFunDer(x, y, z, dN);
                J.addMatrixProduct(0.0, X, dN, 1.0);
                if (J.Invert(invJ) < 0)
                    return false;
                double rx = gx - px;
                double ry = gy - py;
                double rz = gz - pz;
                double dx = invJ(0, 0) * rx + invJ(0, 1) * ry + invJ(0, 2) * rz;
                double dy = invJ(1, 0) * rx + invJ(1, 1) * ry + invJ(1, 2) * rz;
                double dz = invJ(2, 0) * rx + invJ(2, 1) * ry + invJ(2, 2) * rz;
                x += dx;
                y += dy;
                z += dz;
                if (std::sqrt(dx * dx + dy * dy + dz * dz) < 1.0e-12) {
                    lx = x; ly = y; lz = z;
                    return true;
                }
            }
            lx = x; ly = y; lz = z;
            return false;
        }

        // volume by 2x2x2 Gauss. Not exact for a distorted hexahedron (detJ is
        // triquadratic), but it only sets the penalty length scale through a
        // cube root, so the approximation is irrelevant there.
        double volume(const Matrix& X) {
            static const double g = 0.577350269189625764509;
            static const double p[2] = { -g, g };
            static Matrix dN(8, 3);
            static Matrix J(3, 3);
            double V = 0.0;
            for (int a = 0; a < 2; ++a) {
                for (int b = 0; b < 2; ++b) {
                    for (int c = 0; c < 2; ++c) {
                        shapeFunDer(p[a], p[b], p[c], dN);
                        J.addMatrixProduct(0.0, X, dN, 1.0);
                        V += det3(J);
                    }
                }
            }
            return V;
        }

    }

    // ------------------------------------------------------------------
    // Generic assembler shared by the isoparametric families.
    // Builds B for the requested constraint mode and returns K = B^T C B on
    // the reduced dofset, whose ordering matches the one used by the legacy
    // simplex kernels:
    //   [ constrained dofs | retained node 1 dofs | ... | retained node n ]
    // ------------------------------------------------------------------
    const Matrix& assembleConstraint(
        int ndm, int nn, const Vector& N, const Matrix& dNdX,
        int mode, double kU, double kP)
    {
        int nrot = (ndm == 2) ? 1 : 3;
        int nc = ndm;                 // dofs of the constrained node
        int nr = ndm;                 // dofs of each retained node
        if (mode == ASDEmbeddedNodeElement::Mode_UR) {
            nc = ndm + nrot;
        }
        else if (mode == ASDEmbeddedNodeElement::Mode_UP) {
            nc = ndm + 1;
            nr = ndm + 1;
        }
        int nrows = nc;               // one row per constrained dof
        int ncols = nc + nr * nn;

        static Matrix B;
        B.resize(nrows, ncols);
        B.Zero();

        // constrained node: -I
        for (int i = 0; i < nrows; ++i)
            B(i, i) = -1.0;

        // retained nodes
        for (int i = 0; i < nn; ++i) {
            int j = nc + i * nr;
            for (int d = 0; d < ndm; ++d)
                B(d, j + d) = N(i);
            if (mode == ASDEmbeddedNodeElement::Mode_UP) {
                B(ndm, j + ndm) = N(i);
            }
            else if (mode == ASDEmbeddedNodeElement::Mode_UR) {
                if (ndm == 2) {
                    // Rz = (d_uy_dx - d_ux_dy) / 2
                    B(2, j) = -dNdX(i, 1) / 2.0;
                    B(2, j + 1) = dNdX(i, 0) / 2.0;
                }
                else {
                    // R = axial vector of the skew part of the gradient
                    B(3, j + 1) = -dNdX(i, 2) / 2.0; B(3, j + 2) = dNdX(i, 1) / 2.0;
                    B(4, j) = dNdX(i, 2) / 2.0;      B(4, j + 2) = -dNdX(i, 0) / 2.0;
                    B(5, j) = -dNdX(i, 1) / 2.0;     B(5, j + 1) = dNdX(i, 0) / 2.0;
                }
            }
        }

        static Matrix K;
        K.resize(ncols, ncols);
        K.Zero();
        if (mode == ASDEmbeddedNodeElement::Mode_UP) {
            static Matrix C;
            C.resize(nrows, nrows);
            C.Zero();
            for (int i = 0; i < ndm; ++i)
                C(i, i) = kU;
            C(ndm, ndm) = kP;
            K.addMatrixTripleProduct(0.0, B, C, 1.0);
        }
        else {
            K.addMatrixTransposeProduct(0.0, B, B, kU);
        }
        return K;
    }

    namespace tet {

        double shapeFun(double x, double y, double z, int i) {
            if (i == 0)
                return 1.0 - (x + y + z);
            else if (i == 1)
                return x;
            else if (i == 2)
                return y;
            else if (i == 3)
                return z;
            return 0.0;
        }

        void shapeFunDer(Matrix& dN) {
            dN(0, 0) = -1.0; dN(0, 1) = -1.0; dN(0, 2) = -1.0;
            dN(1, 0) = 1.0; dN(1, 1) = 0.0; dN(1, 2) = 0.0;
            dN(2, 0) = 0.0; dN(2, 1) = 1.0; dN(2, 2) = 0.0;
            dN(3, 0) = 0.0; dN(3, 1) = 0.0; dN(3, 2) = 1.0;
        }

        void globalCoord(const Matrix& X, double lx, double ly, double lz, double& gx, double& gy, double& gz) {
            gx = gy = gz = 0.0;
            for (int i = 0; i < 4; i++) {
                double N = shapeFun(lx, ly, lz, i);
                gx += N * X(0, i);
                gy += N * X(1, i);
                gz += N * X(2, i);
            }
        }

        void localCoord(const Matrix& X, const Matrix& invJ, double gx, double gy, double gz, double& lx, double& ly, double& lz) {
            lx = ly = lz = 0.0;
            double px, py, pz;
            globalCoord(X, lx, ly, lz, px, py, pz);
            Vector D(3);
            Vector DL(3);
            D(0) = gx - px;
            D(1) = gy - py;
            D(2) = gz - pz;
            DL.addMatrixVector(0.0, invJ, D, 1.0);
            lx = DL(0);
            ly = DL(1);
            lz = DL(2);
        }

    }

}

void *
OPS_ASDEmbeddedNodeElement(void)
{
    static bool first_done = false;
    if (!first_done) {
        opserr << "Using ASDEmbeddedNodeElement - Developed by: Massimo Petracca, Guido Camata, ASDEA Software Technology\n";
        first_done = true;
    }

    const char* descr = "Want: element ASDEmbeddedNodeElement $tag $Cnode $Rnode1 $Rnode2 $Rnode3 <$Rnode4 ... $Rnode8> <-rot> <-shearDeformable> <-corotational> <-p> <-K $K> <-KP $KP> <-shape $shape> <-slip $slipMatTag $realNodeTag $KS $Xx $Xy $Xz>\n"
        "   3 retained nodes = triangle (2D or 3D)\n"
        "   4 retained nodes = quadrilateral in 2D; in 3D a tetrahedron, or a\n"
        "                      quadrilateral face with -shape quad\n"
        "   8 retained nodes = hexahedron in 3D\n"
        "   -shape tri|quad|tet|hexa: only needed to tell a quadrilateral face\n"
        "                      from a tetrahedron, the one ambiguous case.\n"
        "   -shearDeformable: with -rot on a 3D surface host whose nodes carry\n"
        "                      rotations (a shell), tie the bending rotations of\n"
        "                      the constrained node to the interpolated nodal\n"
        "                      rotations of the host instead of the slope of the\n"
        "                      transverse displacement (correct for thick shells,\n"
        "                      where rotation = slope + shear deformation).\n"
        "   -corotational: with -rot on a 3D host, make the constraint exact\n"
        "                      under finite rotations of the host patch.\n"
        "   -slip: absorb the rebar-slip zeroLength: $Cnode is the AUX node\n"
        "                      embedded in the host, $realNodeTag is the real\n"
        "                      rebar node, tied to it by the uniaxial material\n"
        "                      $slipMatTag along the bar axis ($Xx $Xy $Xz, in\n"
        "                      the reference configuration) and by the stiff\n"
        "                      elastic constant $KS [F/L] on every other\n"
        "                      relative dof. With -corotational the bar axis\n"
        "                      rotates with the host frame.\n";

    int numArgs = OPS_GetNumRemainingInputArgs();
    if (numArgs < 5) {
        opserr << "ASDEmbeddedNodeElement ERROR : Few arguments:\n" << descr;
        return 0;
    }
    
    // mandatory parameters
    int iData[5];
    int numData = 5;
    if (OPS_GetInt(&numData, iData) != 0) {
        opserr << "ASDEmbeddedNodeElement ERROR: Invalid integer mandatory values: element ASDEmbeddedNodeElement wants at least 5 integer parameters\n" << descr;
        return 0;
    }

    // the retained nodes: at least 3, and further ones may follow before the
    // first keyword. Reading them greedily (instead of only accepting one extra
    // node at a fixed position) is what allows 4 and 8 node hosts.
    ID rNodes(0, 8);
    rNodes[0] = iData[2];
    rNodes[1] = iData[3];
    rNodes[2] = iData[4];

    // parse optional parameters
    bool rot = false;
    bool shear = false;
    bool corot = false;
    bool pressure = false;
    bool keywords_started = false;
    double K = 1.0e18;
    double KP = 1.0e18;
    bool KP_set = false;
    int shape = ASDEmbeddedNodeElement::Fam_Unknown;
    UniaxialMaterial* slip_mat = nullptr;
    int slip_node = 0;
    double slip_KS = 0.0;
    Vector slip_x(3);
    for (int i = 5; i < numArgs; i++) {
        const char* what = OPS_GetString();
        if (strcmp(what, "-rot") == 0) {
            rot = true;
            keywords_started = true;
        }
        else if (strcmp(what, "-shearDeformable") == 0) {
            shear = true;
            keywords_started = true;
        }
        else if (strcmp(what, "-corotational") == 0) {
            corot = true;
            keywords_started = true;
        }
        else if (strcmp(what, "-p") == 0) {
            pressure = true;
            keywords_started = true;
        }
        else if (strcmp(what, "-K") == 0) {
            keywords_started = true;
            if (i == numArgs - 1) {
                opserr << "ASDEmbeddedNodeElement ERROR: The -K keyword should be followed by a floating point number.\n" << descr;
                return 0;
            }
            ++i;
            numData = 1;
            if (OPS_GetDouble(&numData, &K) != 0) {
                opserr << "ASDEmbeddedNodeElement ERROR invalid floating point number for -K keyword.\n";
                return 0;
            }
        }
        else if (strcmp(what, "-KP") == 0) {
            keywords_started = true;
            if (i == numArgs - 1) {
                opserr << "ASDEmbeddedNodeElement ERROR: The -KP keyword should be followed by a floating point number.\n" << descr;
                return 0;
            }
            ++i;
            numData = 1;
            if (OPS_GetDouble(&numData, &KP) != 0) {
                opserr << "ASDEmbeddedNodeElement ERROR invalid floating point number for -K keyword.\n";
                return 0;
            }
            KP_set = true;
        }
        else if (strcmp(what, "-slip") == 0) {
            keywords_started = true;
            // 2 integers (slipMatTag, realNodeTag) + 4 doubles (KS, Xx, Xy, Xz)
            if (numArgs - i - 1 < 6) {
                opserr << "ASDEmbeddedNodeElement ERROR: the -slip keyword wants "
                    << "$slipMatTag $realNodeTag $KS $Xx $Xy $Xz.\n" << descr;
                return 0;
            }
            int slipInt[2];
            numData = 2;
            if (OPS_GetInt(&numData, slipInt) != 0) {
                opserr << "ASDEmbeddedNodeElement ERROR: invalid integer values for the -slip "
                    << "keyword: it wants $slipMatTag $realNodeTag $KS $Xx $Xy $Xz.\n" << descr;
                return 0;
            }
            double slipDouble[4];
            numData = 4;
            if (OPS_GetDouble(&numData, slipDouble) != 0) {
                opserr << "ASDEmbeddedNodeElement ERROR: invalid floating point values for the "
                    << "-slip keyword: it wants $slipMatTag $realNodeTag $KS $Xx $Xy $Xz.\n" << descr;
                return 0;
            }
            i += 6;
            slip_mat = OPS_getUniaxialMaterial(slipInt[0]);
            if (slip_mat == nullptr) {
                opserr << "ASDEmbeddedNodeElement ERROR: -slip refers to uniaxialMaterial "
                    << slipInt[0] << ", which does not exist. Define the tau-slip law "
                    << "before the element.\n";
                return 0;
            }
            slip_node = slipInt[1];
            slip_KS = slipDouble[0];
            if (slip_KS <= 0.0) {
                opserr << "ASDEmbeddedNodeElement ERROR: -slip wants a positive rigid-tie "
                    << "stiffness $KS [F/L], got " << slip_KS << ". Use the same value the "
                    << "rigid Elastic material of the zeroLength assembly would take.\n";
                return 0;
            }
            slip_x(0) = slipDouble[1];
            slip_x(1) = slipDouble[2];
            slip_x(2) = slipDouble[3];
            if (slip_x.Norm() < 1.0e-12) {
                opserr << "ASDEmbeddedNodeElement ERROR: -slip wants a non-zero bar axis "
                    << "($Xx $Xy $Xz).\n";
                return 0;
            }
            slip_x.Normalize();
        }
        else if (strcmp(what, "-shape") == 0) {
            keywords_started = true;
            if (i == numArgs - 1) {
                opserr << "ASDEmbeddedNodeElement ERROR: The -shape keyword should be followed by "
                    << "one of tri, quad, tet, hexa.\n" << descr;
                return 0;
            }
            ++i;
            const char* sname = OPS_GetString();
            if (strcmp(sname, "tri") == 0)
                shape = ASDEmbeddedNodeElement::Fam_Tri;
            else if (strcmp(sname, "quad") == 0)
                shape = ASDEmbeddedNodeElement::Fam_Quad;
            else if (strcmp(sname, "tet") == 0)
                shape = ASDEmbeddedNodeElement::Fam_Tet;
            else if (strcmp(sname, "hexa") == 0)
                shape = ASDEmbeddedNodeElement::Fam_Hexa;
            else {
                opserr << "ASDEmbeddedNodeElement ERROR: unknown -shape \"" << sname
                    << "\". Use one of tri, quad, tet, hexa.\n" << descr;
                return 0;
            }
        }
        else {
            // an extra retained node: only accepted before any keyword, so that
            // a stray token later on is reported instead of being swallowed
            if (keywords_started) {
                opserr << "ASDEmbeddedNodeElement ERROR: unexpected argument \"" << what
                    << "\" after the keywords.\n" << descr;
                return 0;
            }
            char* endptr = 0;
            long value = strtol(what, &endptr, 10);
            if (endptr == what || *endptr != '\0') {
                opserr << "ASDEmbeddedNodeElement ERROR: expected an integer node tag or a "
                    << "keyword, got \"" << what << "\".\n" << descr;
                return 0;
            }
            if (rNodes.Size() >= 8) {
                opserr << "ASDEmbeddedNodeElement ERROR: at most 8 retained nodes are allowed.\n" << descr;
                return 0;
            }
            rNodes[rNodes.Size()] = static_cast<int>(value);
        }
    }

    // only 3, 4 and 8 retained nodes correspond to a supported host
    int nret = rNodes.Size();
    if (nret != 3 && nret != 4 && nret != 8) {
        opserr << "ASDEmbeddedNodeElement ERROR: " << nret << " retained nodes given. "
            << "Only 3 (triangle), 4 (quadrilateral in 2D / tetrahedron in 3D) "
            << "and 8 (hexahedron) are supported.\n" << descr;
        return 0;
    }
    if (!KP_set)
        KP = K;

    // a -shape that contradicts the node count is a modelling mistake, not
    // something to silently resolve in favour of one of the two
    if (shape != ASDEmbeddedNodeElement::Fam_Unknown) {
        int want = (shape == ASDEmbeddedNodeElement::Fam_Tri) ? 3 :
                   (shape == ASDEmbeddedNodeElement::Fam_Hexa) ? 8 : 4;
        if (nret != want) {
            opserr << "ASDEmbeddedNodeElement ERROR: -shape asks for a host with " << want
                << " nodes but " << nret << " retained nodes were given.\n" << descr;
            return 0;
        }
    }

    // check
    if (pressure && rot) {
        opserr << "ASDEmbeddedNodeElement ERROR: Cannot use both -rot and -p flags.\n" << descr;
        return 0;
    }
    if (slip_mat && pressure) {
        opserr << "ASDEmbeddedNodeElement ERROR: -slip cannot be combined with -p: the "
            << "rebar-slip assembly is not defined on u-p nodes.\n" << descr;
        return 0;
    }
    if (slip_mat && (slip_node == iData[1] || rNodes.getLocation(slip_node) >= 0)) {
        opserr << "ASDEmbeddedNodeElement ERROR: -slip $realNodeTag (" << slip_node
            << ") must be a node OUTSIDE the element: not the constrained (AUX) node "
            << "and not a retained host node.\n";
        return 0;
    }
    if (shear && !rot) {
        opserr << "ASDEmbeddedNodeElement ERROR: -shearDeformable only modifies the rotational "
            << "constraint, so it requires -rot.\n" << descr;
        return 0;
    }
    if (corot) {
        if (!rot || pressure) {
            opserr << "ASDEmbeddedNodeElement ERROR: -corotational requires -rot and cannot be "
                << "combined with -p.\n" << descr;
            return 0;
        }
        if (nret != 8 && nret != 4 && nret != 3) {
            opserr << "ASDEmbeddedNodeElement ERROR: -corotational needs a 3D host "
                << "(3/4/8 retained nodes).\n" << descr;
            return 0;
        }
    }

    // done
    return new ASDEmbeddedNodeElement(iData[0], iData[1], rNodes, rot, pressure, K, KP, shape, shear, corot,
        slip_mat, slip_node, slip_KS, slip_mat ? &slip_x : nullptr);
}

ASDEmbeddedNodeElement::ASDEmbeddedNodeElement() 
    : Element(0, ELE_TAG_ASDEmbeddedNodeElement)
{
}

ASDEmbeddedNodeElement::ASDEmbeddedNodeElement(int tag, int cNode, const ID& rNodes, bool rot_flag, bool p_flag, double K, double KP, int shape_request, bool shear_flag, bool corot_flag,
    UniaxialMaterial* slip_mat, int slip_node, double KS, const Vector* slip_x)
    : Element(tag, ELE_TAG_ASDEmbeddedNodeElement)
    , m_shape_request(shape_request)
    , m_rot_c_flag(rot_flag)
    , m_p_flag(p_flag)
    , m_shear_flag(shear_flag)
    , m_corot_flag(corot_flag)
    , m_K(K)
    , m_KP(KP)
    , m_KS(KS)
{
    int nn = rNodes.Size();
    m_slip = (slip_mat != nullptr);
    int extra = m_slip ? 2 : 1; // constrained node (+ the real rebar node)
    m_node_ids.resize(nn + extra);
    m_node_ids(0) = cNode;
    for (int i = 0; i < nn; ++i)
        m_node_ids(i + 1) = rNodes(i);
    if (m_slip) {
        m_node_ids(nn + 1) = slip_node;
        m_slip_mat = slip_mat->getCopy();
        if (m_slip_mat == nullptr) {
            opserr << "ASDEmbeddedNodeElement ERROR: failed to copy the -slip uniaxial material "
                << slip_mat->getTag() << "\n";
            exit(-1);
        }
        m_slip_x0 = *slip_x;
    }
    m_nodes.resize(static_cast<std::size_t>(nn + extra), nullptr);
}

ASDEmbeddedNodeElement::~ASDEmbeddedNodeElement( )
{
    if (m_slip_mat)
        delete m_slip_mat;
}

const char* ASDEmbeddedNodeElement::getClassType(void) const
{
    return "ASDEmbeddedNodeElement";
}

int ASDEmbeddedNodeElement::numRetained() const
{
    // the first node is the constrained one and, with -slip, the last one is
    // the real rebar node: neither is a host node
    return static_cast<int>(m_nodes.size()) - 1 - (m_slip ? 1 : 0);
}

int ASDEmbeddedNodeElement::resolveFamily() const
{
    // number of retained nodes (the first one is the constrained node)
    int nn = numRetained();
    if (nn == 3)
        return Fam_Tri;                                  // 2D or 3D (shell face)
    if (nn == 4) {
        if (m_ndm == 2)
            return Fam_Quad;
        // 4 retained nodes in 3D are the one genuinely ambiguous case: a
        // tetrahedron or a quadrilateral face. Only -shape can tell them apart,
        // and its absence keeps the historical meaning so that every model
        // written before the face host existed still reads the same way.
        if (m_shape_request == Fam_Quad || m_shape_request == Fam_Quad3D)
            return Fam_Quad3D;
        return Fam_Tet;
    }
    if (nn == 8)
        return (m_ndm == 3) ? Fam_Hexa : Fam_Unknown;
    return Fam_Unknown;
}

void ASDEmbeddedNodeElement::setDomain(Domain* theDomain)
{
    // check nodes
    m_num_dofs = 0;
    int local_dof_counter = 0;
    int local_pos = 0;
    int aux_ndf = 0; // dof signature of the constrained (AUX) node
    // 1 + the retained nodes: the -slip real node (last) is not one of them
    std::size_t host_end = static_cast<std::size_t>(numRetained()) + 1;
    std::vector<ID> aux_mapping(m_nodes.size());
    for (std::size_t i = 0; i < m_nodes.size(); ++i) {

        // the -slip real rebar node travels last and is NOT a host node
        bool is_slip_node = m_slip && (i + 1 == m_nodes.size());

        // check node
        int node_id = m_node_ids(static_cast<int>(i));
        Node* node = theDomain->getNode(node_id);
        if (node == nullptr) {
            opserr << "ASDEmbeddedNodeElement Error in setDomain: node " << node_id << " does not exit in the domain\n";
            exit(-1);
        }

        // store node
        m_nodes[i] = node;

        // check NDM
        int ndm = node->getCrds().Size();
        if (ndm != 2 && ndm != 3) {
            opserr << "ASDEmbeddedNodeElement Error in setDomain: Nodes should have either 2 or 3 dimensions, not " << ndm << "\n";
            exit(-1);
        }
        if (i == 0) {
            // save ndm at first node
            m_ndm = ndm;
        }
        else {
            if (m_ndm != ndm) {
                opserr << "ASDEmbeddedNodeElement Error in setDomain: Nodes should have the same dimension (2 or 3)\n";
                exit(-1);
            }
        }

        // check NDF
        int ndf = node->getNumberDOF();
        if (i == 0)
            aux_ndf = ndf;
        // -slip: the real rebar node must carry the SAME signature as the AUX
        // node (STKO generates the AUX as a copy of it); the rotational tie is
        // active iff that shared signature has rotations, independently of
        // -rot, exactly as the zeroLength ties -dir 1..ndf
        if (is_slip_node) {
            if (ndf != aux_ndf) {
                opserr << "ASDEmbeddedNodeElement Error in setDomain: element " << getTag()
                    << " - the -slip real node " << node_id << " has " << ndf
                    << " dofs but the constrained (AUX) node has " << aux_ndf
                    << ". Generate the AUX node with the same dof signature as the rebar node.\n";
                exit(-1);
            }
            m_slip_rot = (m_ndm == 3) ? (ndf == 6) : (ndf == 3);
        }
        if (m_ndm == 2) {
            if (ndf != 2 && ndf != 3) {
                opserr << "ASDEmbeddedNodeElement Error in setDomain: In 2D only 2 or 3 DOFs are allowed, not " << ndf << "\n";
                exit(-1);
            }
            if (i == 0) {
                m_rot_c = (m_rot_c_flag && ndf == 3);
                m_shear = false;
                if (m_shear_flag) {
                    opserr << "ASDEmbeddedNodeElement Error in setDomain: element " << getTag()
                        << " - -shearDeformable requires a surface host in 3D: a triangle, or a "
                        << "quadrilateral face with -shape quad.\n";
                    exit(-1);
                }
                if (m_p_flag && ndf == 3) {
                    // all others should have same ndf (u-p)
                    m_up = true;
                    for (std::size_t other_i = 1; other_i < host_end; ++other_i) {
                        int other_node_id = m_node_ids(static_cast<int>(other_i));
                        Node* other_node = theDomain->getNode(other_node_id);
                        if (other_node && (other_node->getNumberDOF() != ndf)) {
                            m_up = false;
                            break;
                        }
                    }
                }
            }
        }
        else {
            if (ndf != 3 && ndf != 4 && ndf != 6) {
                opserr << "ASDEmbeddedNodeElement Error in setDomain: In 3D only 3, 4 or 6 DOFs are allowed, not " << ndf << "\n";
                exit(-1);
            }
            if (i == 0) {
                m_rot_c = (m_rot_c_flag && ndf == 6);
                // -shearDeformable: decide it here, before the dof mapping of the
                // retained nodes is built, because accepting it makes them expose
                // their rotations to the constraint. Misuse is a modelling mistake
                // and stops the analysis; the one silent downgrade is the same one
                // -rot has: a constrained node without rotational dofs turns the
                // whole rotational constraint off, and this refinement with it.
                m_shear = false;
                if (m_shear_flag && m_rot_c) {
                    int nret = numRetained();
                    bool surface = (nret == 3) ||
                        (nret == 4 && (m_shape_request == Fam_Quad || m_shape_request == Fam_Quad3D));
                    if (!surface) {
                        opserr << "ASDEmbeddedNodeElement Error in setDomain: element " << getTag()
                            << " - -shearDeformable requires a surface host in 3D: a triangle, or a "
                            << "quadrilateral face with -shape quad.\n";
                        exit(-1);
                    }
                    for (std::size_t other_i = 1; other_i < host_end; ++other_i) {
                        Node* other_node = theDomain->getNode(m_node_ids(static_cast<int>(other_i)));
                        if (other_node && other_node->getNumberDOF() != 6) {
                            opserr << "ASDEmbeddedNodeElement Error in setDomain: element " << getTag()
                                << " - -shearDeformable ties the constrained rotations to the nodal "
                                << "rotations of the host, so every retained node needs 6 dofs; node "
                                << m_node_ids(static_cast<int>(other_i)) << " has "
                                << other_node->getNumberDOF() << ".\n";
                            exit(-1);
                        }
                    }
                    m_shear = true;
                }
                // -corotational: same silent downgrade as -rot when the
                // constrained node has no rotational dofs. It combines freely
                // with -shearDeformable, whose own checks above already
                // guarantee the surface host and the 6-dof retained nodes.
                m_corot = m_corot_flag && m_rot_c;
                if (m_p_flag && ndf == 4) {
                    // all others should have same ndf (u-p)
                    m_up = true;
                    for (std::size_t other_i = 1; other_i < host_end; ++other_i) {
                        int other_node_id = m_node_ids(static_cast<int>(other_i));
                        Node* other_node = theDomain->getNode(other_node_id);
                        if (other_node && (other_node->getNumberDOF() != ndf)) {
                            m_up = false;
                            break;
                        }
                    }
                }
            }
        }

        // set up mapping
        ID& imap = aux_mapping[i];
        int imap_size = m_ndm;
        if (is_slip_node) {
            // the real rebar node exposes its translations and, when the tie
            // is rotational, its rotations; never the -shear/-p logic of the
            // retained nodes
            if (m_slip_rot)
                imap_size += (m_ndm == 2) ? 1 : 3;
        }
        else if (m_rot_c) {
            if (i == 0) {
                if (m_ndm == 2)
                    imap_size += 1;
                else
                    imap_size += 3;
            }
            else if (m_shear) {
                // -shearDeformable: the retained nodes expose their rotations too
                imap_size += 3;
            }
        }
        else if (m_up) {
            imap_size += 1;
        }
        imap.resize(imap_size);
        imap(0) = local_pos; // Ux
        imap(1) = local_pos + 1; // Uy
        if (m_ndm == 3) {
            imap(2) = local_pos + 2; // Uz
            if ((is_slip_node && m_slip_rot) ||
                (!is_slip_node && m_rot_c && (i == 0 || m_shear))) {
                imap(3) = local_pos + 3; // Rx
                imap(4) = local_pos + 4; // Ry
                imap(5) = local_pos + 5; // Rz
            }
            else if (!is_slip_node && m_up) {
                imap(3) = local_pos + 3; // P
            }
        }
        else {
            if ((is_slip_node && m_slip_rot) || (!is_slip_node && i == 0 && m_rot_c)) {
                imap(2) = local_pos + 2; // Rz
            }
            else if (!is_slip_node && m_up) {
                imap(2) = local_pos + 2; // P
            }
        }
        local_pos += ndf;

        // update total dof counter
        m_num_dofs += ndf;
        // update local dof counter
        local_dof_counter += imap.Size();
    }

    // resolve the geometric family of the host now that ndm is known
    m_family = resolveFamily();
    if (m_family == Fam_Unknown) {
        opserr << "ASDEmbeddedNodeElement Error in setDomain: element " << getTag()
            << " has " << static_cast<int>(m_nodes.size()) - 1
            << " retained nodes in " << m_ndm << "D.\n"
            << "Supported hosts: 3 (triangle, 2D/3D), 4 (quadrilateral in 2D; in 3D a "
            << "tetrahedron, or a quadrilateral face with -shape quad), 8 (hexahedron in 3D).\n";
        exit(-1);
    }

    // A quadrilateral face is handled on its mean plane, the same plane the
    // shell it usually belongs to uses. Warn when the face is warped enough
    // that the projection stops being a detail.
    if (m_family == Fam_Quad3D) {
        static Matrix Xw(3, 4);
        for (int i = 1; i < 5; i++) {
            const Vector& c = m_nodes[static_cast<std::size_t>(i)]->getCrds();
            Xw(0, i - 1) = c(0);
            Xw(1, i - 1) = c(1);
            Xw(2, i - 1) = c(2);
        }
        static Matrix Rw(3, 3);
        static Vector cw(3);
        static Matrix XLw(2, 4);
        double warp = 0.0;
        quad::frame3D(Xw, Rw, cw, XLw, warp);
        double h = std::sqrt(std::abs(quad::area(XLw)));
        if (h > 0.0 && warp > 0.05 * h) {
            opserr << "ASDEmbeddedNodeElement WARNING: element " << getTag()
                << " - the quadrilateral face is warped: the nodes are up to "
                << warp / h * 100.0 << "% of the face size off their mean plane. "
                << "The constraint is written on that plane, as the shell formulation is.\n";
        }
    }

    // -slip: reference triad and coincidence check
    if (m_slip) {
        if (m_ndm == 2 && std::abs(m_slip_x0(2)) > 1.0e-12) {
            opserr << "ASDEmbeddedNodeElement Error in setDomain: element " << getTag()
                << " - in a 2D model the -slip bar axis must lie in the XY plane "
                << "(got Xz = " << m_slip_x0(2) << "). Set Xz to 0.\n";
            exit(-1);
        }
        // complete x0 to a triad with the same rule STKO's frame_from_x uses,
        // so the reported transverse components match the legacy zeroLength
        // assembly. Any completion is mechanically equivalent: the rigid tie
        // is isotropic in the plane orthogonal to the bar.
        static Vector ref(3), y0(3), z0(3);
        ref.Zero();
        if (std::abs(m_slip_x0(2)) > 0.99)
            ref(0) = 1.0;
        else
            ref(2) = 1.0;
        cross(ref, m_slip_x0, y0);
        y0.Normalize();
        cross(m_slip_x0, y0, z0);
        m_slip_T0.resize(3, 3);
        for (int j = 0; j < 3; ++j) {
            m_slip_T0(0, j) = m_slip_x0(j);
            m_slip_T0(1, j) = y0(j);
            m_slip_T0(2, j) = z0(j);
        }
        // the AUX node is expected to coincide with the real rebar node: a
        // gap is not an error (the tie still works on relative displacements,
        // as the zeroLength does) but it is worth knowing about
        const Vector& ca = m_nodes[0]->getCrds();
        const Vector& cr = m_nodes.back()->getCrds();
        double dist2 = 0.0, scale2 = 0.0;
        for (int j = 0; j < m_ndm; ++j) {
            double d = ca(j) - cr(j);
            dist2 += d * d;
        }
        for (int a = 1; a < 1 + numRetained(); ++a)
            for (int b = a + 1; b < 1 + numRetained(); ++b) {
                const Vector& xa = m_nodes[static_cast<std::size_t>(a)]->getCrds();
                const Vector& xb = m_nodes[static_cast<std::size_t>(b)]->getCrds();
                double s2 = 0.0;
                for (int j = 0; j < m_ndm; ++j) {
                    double d = xa(j) - xb(j);
                    s2 += d * d;
                }
                scale2 = std::max(scale2, s2);
            }
        if (scale2 > 0.0 && dist2 > 1.0e-12 * scale2) {
            opserr << "ASDEmbeddedNodeElement WARNING: element " << getTag()
                << " - the -slip real node " << m_node_ids(m_node_ids.Size() - 1)
                << " is not coincident with the constrained (AUX) node "
                << m_node_ids(0) << " (distance = " << std::sqrt(dist2) << ").\n";
        }
        // recorder buffers
        int ns = m_ndm + (m_slip_rot ? (m_ndm == 2 ? 1 : 3) : 0);
        m_slip_g.resize(ns);
        m_slip_g.Zero();
        m_slip_axis.resize(3);
        for (int j = 0; j < 3; ++j)
            m_slip_axis(j) = m_slip_x0(j);
    }

    // flatten mapping
    m_mapping.resize(local_dof_counter);
    local_pos = 0;
    for (const ID& imap : aux_mapping) {
        for (int i = 0; i < imap.Size(); ++i) {
            m_mapping(local_pos++) = imap(i);
        }
    }

    // compute initial displacement vector
    if (!m_U0_computed) {
        m_U0.resize(m_num_dofs);
        m_U0 = getGlobalDisplacements();
        m_U0_computed = true;
    }

    // -corotational + -shearDeformable: per-retained-node quaternion state.
    // Allocate to identity only when the size is wrong, so that a state
    // restored by recvSelf (already correctly sized) is preserved.
    if (m_corot && m_shear) {
        std::size_t nret = m_nodes.size() - 1;
        if (m_qa.size() != 4 * nret) {
            m_qa.assign(4 * nret, 0.0);
            m_qa_conv.assign(4 * nret, 0.0);
            m_rva.assign(3 * nret, 0.0);
            m_rva_conv.assign(3 * nret, 0.0);
            for (std::size_t a = 0; a < nret; ++a) {
                m_qa[4 * a] = 1.0;
                m_qa_conv[4 * a] = 1.0;
            }
        }
    }

    // call base class implementation
    DomainComponent::setDomain(theDomain);
}

void ASDEmbeddedNodeElement::Print(OPS_Stream& s, int flag)
{
    if (flag == -1) {
        int eleTag = this->getTag();
        s << "EL_ASDEmbeddedNodeElement\t" << eleTag << " :";
        for (int i = 0; i < m_node_ids.Size(); ++i)
            s << "\t" << m_node_ids(i);
        s << endln;
    }

    if (flag == OPS_PRINT_PRINTMODEL_JSON) {
        s << "\t\t\t{";
        s << "\"name\": " << this->getTag() << ", ";
        s << "\"type\": \"ASDEmbeddedNodeElement\", ";
        s << "\"nodes\": [";
        for (int i = 0; i < m_node_ids.Size(); ++i) {
            if (i > 0)
                s << ", ";
            s << m_node_ids(i);
        }
        s << "]";
        if (m_slip && m_slip_mat)
            s << ", \"slipMaterial\": " << m_slip_mat->getTag();
        s << "}";
    }
}

int ASDEmbeddedNodeElement::getNumExternalNodes() const
{
    return m_node_ids.Size();
}

const ID& ASDEmbeddedNodeElement::getExternalNodes()
{
    return m_node_ids;
}

Node**
ASDEmbeddedNodeElement::getNodePtrs(void)
{
    return m_nodes.data();
}

int ASDEmbeddedNodeElement::getNumDOF()
{
    return m_num_dofs;
}

int ASDEmbeddedNodeElement::update()
{
    if (m_corot) {
        // track the slave TOTAL rotation as a quaternion, composed from the
        // increments of the additive rotation dofs -- the same scheme the ASD
        // shells use for their nodal quaternions, so that two corotational
        // elements on the same node agree
        const Vector& iu = m_nodes[0]->getTrialDisp();
        double rv[3];
        for (int i = 0; i < 3; ++i) {
            rv[i] = iu(3 + i) - (m_U0_computed ? m_U0(3 + i) : 0.0);
        }
        ASDQuaternion<double> dq = ASDQuaternion<double>::FromRotationVector(
            rv[0] - m_rv[0], rv[1] - m_rv[1], rv[2] - m_rv[2]);
        ASDQuaternion<double> q(m_qs[0], m_qs[1], m_qs[2], m_qs[3]);
        q = dq * q;
        q.normalize();
        m_qs[0] = q.w(); m_qs[1] = q.x(); m_qs[2] = q.y(); m_qs[3] = q.z();
        for (int i = 0; i < 3; ++i) m_rv[i] = rv[i];
        // -shearDeformable: the same bookkeeping for every retained node
        if (m_shear) {
            int nn = numRetained();
            int pos = m_nodes[0]->getNumberDOF();
            for (int a = 0; a < nn; ++a) {
                const Vector& au = m_nodes[static_cast<std::size_t>(a + 1)]->getTrialDisp();
                double arv[3];
                for (int i = 0; i < 3; ++i)
                    arv[i] = au(3 + i) - (m_U0_computed ? m_U0(pos + 3 + i) : 0.0);
                ASDQuaternion<double> adq = ASDQuaternion<double>::FromRotationVector(
                    arv[0] - m_rva[3 * a], arv[1] - m_rva[3 * a + 1], arv[2] - m_rva[3 * a + 2]);
                ASDQuaternion<double> aq(m_qa[4 * a], m_qa[4 * a + 1], m_qa[4 * a + 2], m_qa[4 * a + 3]);
                aq = adq * aq;
                aq.normalize();
                m_qa[4 * a] = aq.w(); m_qa[4 * a + 1] = aq.x();
                m_qa[4 * a + 2] = aq.y(); m_qa[4 * a + 3] = aq.z();
                for (int i = 0; i < 3; ++i) m_rva[3 * a + i] = arv[i];
                pos += m_nodes[static_cast<std::size_t>(a + 1)]->getNumberDOF();
            }
        }
        // -slip rotational tie: the real rebar node gets the same bookkeeping
        if (m_slip && m_slip_rot) {
            const Vector& ru = m_nodes.back()->getTrialDisp();
            int rpos = m_num_dofs - m_nodes.back()->getNumberDOF();
            double rrv[3];
            for (int i = 0; i < 3; ++i)
                rrv[i] = ru(3 + i) - (m_U0_computed ? m_U0(rpos + 3 + i) : 0.0);
            ASDQuaternion<double> rdq = ASDQuaternion<double>::FromRotationVector(
                rrv[0] - m_rvr[0], rrv[1] - m_rvr[1], rrv[2] - m_rvr[2]);
            ASDQuaternion<double> qr(m_qr[0], m_qr[1], m_qr[2], m_qr[3]);
            qr = rdq * qr;
            qr.normalize();
            m_qr[0] = qr.w(); m_qr[1] = qr.x(); m_qr[2] = qr.y(); m_qr[3] = qr.z();
            for (int i = 0; i < 3; ++i) m_rvr[i] = rrv[i];
        }
    }
    // -slip: strain the tau-slip law with the current relative displacement
    // along the (frozen or corotated) bar axis. The return code travels: an
    // IMPL-EX material signals failure from setTrialStrain to ask for a step
    // cut, and Domain::update aggregates the element codes (as ZeroLength
    // propagates it for the legacy assembly).
    if (m_slip)
        return m_slip_mat->setTrialStrain(slipComputeGap());
    return 0;
}

int ASDEmbeddedNodeElement::commitState()
{
    for (int i = 0; i < 4; ++i) m_qs_conv[i] = m_qs[i];
    for (int i = 0; i < 3; ++i) m_rv_conv[i] = m_rv[i];
    m_qa_conv = m_qa;
    m_rva_conv = m_rva;
    for (int i = 0; i < 4; ++i) m_qr_conv[i] = m_qr[i];
    for (int i = 0; i < 3; ++i) m_rvr_conv[i] = m_rvr[i];
    if (m_slip_mat)
        m_slip_mat->commitState();
    return Element::commitState();
}

int ASDEmbeddedNodeElement::revertToLastCommit()
{
    for (int i = 0; i < 4; ++i) m_qs[i] = m_qs_conv[i];
    for (int i = 0; i < 3; ++i) m_rv[i] = m_rv_conv[i];
    m_qa = m_qa_conv;
    m_rva = m_rva_conv;
    for (int i = 0; i < 4; ++i) m_qr[i] = m_qr_conv[i];
    for (int i = 0; i < 3; ++i) m_rvr[i] = m_rvr_conv[i];
    if (m_slip_mat)
        m_slip_mat->revertToLastCommit();
    return 0;
}

int ASDEmbeddedNodeElement::revertToStart()
{
    m_qs[0] = 1.0; m_qs[1] = m_qs[2] = m_qs[3] = 0.0;
    m_rv[0] = m_rv[1] = m_rv[2] = 0.0;
    for (int i = 0; i < 4; ++i) m_qs_conv[i] = m_qs[i];
    for (int i = 0; i < 3; ++i) m_rv_conv[i] = m_rv[i];
    for (std::size_t a = 0; 4 * a < m_qa.size(); ++a) {
        m_qa[4 * a] = 1.0;
        m_qa[4 * a + 1] = m_qa[4 * a + 2] = m_qa[4 * a + 3] = 0.0;
    }
    std::fill(m_rva.begin(), m_rva.end(), 0.0);
    m_qa_conv = m_qa;
    m_rva_conv = m_rva;
    m_qr[0] = 1.0; m_qr[1] = m_qr[2] = m_qr[3] = 0.0;
    m_rvr[0] = m_rvr[1] = m_rvr[2] = 0.0;
    for (int i = 0; i < 4; ++i) m_qr_conv[i] = m_qr[i];
    for (int i = 0; i < 3; ++i) m_rvr_conv[i] = m_rvr[i];
    if (m_slip_mat)
        m_slip_mat->revertToStart();
    return 0;
}

void ASDEmbeddedNodeElement::corotSetup()
{
    if (m_corot_init)
        return;
    // reference = configuration at activation (X + U0). The center gradients
    // are computed ON that configuration, so F(reference) = I and R0 = I.
    // For the 4-node tet F is constant: point and center gradients coincide,
    // and the same closed-form G applies with the tet gradients.
    int nn = numRetained();
    static Matrix X;
    X.resize(3, nn);
    int pos = m_nodes[0]->getNumberDOF();
    for (int a = 0; a < nn; ++a) {
        const Vector& crd = m_nodes[static_cast<std::size_t>(a + 1)]->getCrds();
        for (int i = 0; i < 3; ++i)
            X(i, a) = crd(i) + (m_U0_computed ? m_U0(pos + i) : 0.0);
        pos += m_nodes[static_cast<std::size_t>(a + 1)]->getNumberDOF();
    }
    static Vector Xs(3);
    for (int i = 0; i < 3; ++i)
        Xs(i) = m_nodes[0]->getCrds()(i) + (m_U0_computed ? m_U0(i) : 0.0);

    m_cN.resize(nn);
    m_cD.resize(nn, 3);
    m_cgc.resize(nn, 3);
    static Matrix J(3, 3);
    static Matrix invJ(3, 3);
    double V = 0.0;
    if (nn == 8) {
        double lx, ly, lz;
        if (!hexa::localCoord(X, Xs(0), Xs(1), Xs(2), lx, ly, lz)) {
            opserr << "ASDEmbeddedNodeElement WARNING: element " << getTag()
                << " - the inverse isoparametric map did not converge on the HEXA host (corotational setup).\n";
        }
        hexa::shapeFun(lx, ly, lz, m_cN);
        static Matrix dN(8, 3);
        hexa::shapeFunDer(lx, ly, lz, dN);
        J.addMatrixProduct(0.0, X, dN, 1.0);
        inv3(J, invJ);
        m_cD.addMatrixProduct(0.0, dN, invJ, 1.0);
        hexa::shapeFunDer(0.0, 0.0, 0.0, dN);
        J.addMatrixProduct(0.0, X, dN, 1.0);
        inv3(J, invJ);
        m_cgc.addMatrixProduct(0.0, dN, invJ, 1.0);
        V = hexa::volume(X);
    }
    else if (nn == 4 && m_family == Fam_Tet) {
        // tetrahedron: affine map, constant gradients
        static Matrix dN(4, 3);
        tet::shapeFunDer(dN);
        J.addMatrixProduct(0.0, X, dN, 1.0);
        inv3(J, invJ);
        m_cgc.addMatrixProduct(0.0, dN, invJ, 1.0);
        m_cD = m_cgc;
        double lx, ly, lz;
        tet::localCoord(X, invJ, Xs(0), Xs(1), Xs(2), lx, ly, lz);
        for (int a = 0; a < 4; ++a)
            m_cN(a) = tet::shapeFun(lx, ly, lz, a);
        V = det3(J) / 6.0;
    }
    else {
        // surface host (3-node, or 4-node face): Kabsch frame. m_cgc rows hold
        // the reference local positions Xh_a, so that F = sum x_a (x) gc_a is
        // the Kabsch covariance A and the SAME closed-form G applies (M stays
        // SPD for a flat patch: eigenvalues lam_i + lam_j with lam3 = 0).
        m_corot_surf = true;
        // reference face frame E0 (columns e1, e2, e3)
        m_cE0.resize(3, 3);
        static Vector e1(3), e2(3), e3(3);
        if (nn == 3) {
            for (int i = 0; i < 3; ++i) {
                e1(i) = X(i, 1) - X(i, 0);
                e2(i) = X(i, 2) - X(i, 0);
            }
        }
        else {
            for (int i = 0; i < 3; ++i) {
                e1(i) = X(i, 1) + X(i, 2) - X(i, 0) - X(i, 3);
                e2(i) = X(i, 2) + X(i, 3) - X(i, 0) - X(i, 1);
            }
        }
        e1.Normalize();
        cross(e1, e2, e3);
        e3.Normalize();
        cross(e3, e1, e2);
        for (int i = 0; i < 3; ++i) {
            m_cE0(i, 0) = e1(i);
            m_cE0(i, 1) = e2(i);
            m_cE0(i, 2) = e3(i);
        }
        // centroid first (needed for the projected 2D coordinates)
        static Vector cc(3);
        cc.Zero();
        for (int a = 0; a < nn; ++a)
            for (int i = 0; i < 3; ++i)
                cc(i) += X(i, a) / static_cast<double>(nn);
        static Matrix XL;
        XL.resize(2, nn);
        for (int a = 0; a < nn; ++a)
            for (int j = 0; j < 2; ++j) {
                double s = 0.0;
                for (int i = 0; i < 3; ++i)
                    s += m_cE0(i, j) * (X(i, a) - cc(i));
                XL(j, a) = s;
            }
        double sxl = 0.0, syl = 0.0;
        for (int i = 0; i < 3; ++i) {
            sxl += m_cE0(i, 0) * (Xs(i) - cc(i));
            syl += m_cE0(i, 1) * (Xs(i) - cc(i));
        }
        // shape functions and 2D cartesian gradients at the material point
        static Matrix J2(2, 2);
        if (nn == 3) {
            static Matrix dN2(3, 2);
            tri::shapeFunDer(dN2);
            J2.addMatrixProduct(0.0, XL, dN2, 1.0);
            static Matrix iJ2(2, 2);
            J2.Invert(iJ2);
            static Matrix D2(3, 2);
            D2.addMatrixProduct(0.0, dN2, iJ2, 1.0);
            double lx, ly;
            tri::localCoord(XL, iJ2, sxl, syl, lx, ly);
            for (int a = 0; a < 3; ++a) {
                m_cN(a) = tri::shapeFun(lx, ly, a);
                m_cD(a, 0) = D2(a, 0);
                m_cD(a, 1) = D2(a, 1);
                m_cD(a, 2) = 0.0;
            }
            V = det2(J2) / 2.0;
        }
        else {
            double lx, ly;
            if (!quad::localCoord(XL, sxl, syl, lx, ly)) {
                opserr << "ASDEmbeddedNodeElement WARNING: element " << getTag()
                    << " - the inverse isoparametric map did not converge on the QUAD face (corotational setup).\n";
            }
            static Vector N4(4);
            quad::shapeFun(lx, ly, N4);
            static Matrix dN2(4, 2);
            quad::shapeFunDer(lx, ly, dN2);
            J2.addMatrixProduct(0.0, XL, dN2, 1.0);
            static Matrix iJ2(2, 2);
            J2.Invert(iJ2);
            static Matrix D2(4, 2);
            D2.addMatrixProduct(0.0, dN2, iJ2, 1.0);
            for (int a = 0; a < 4; ++a) {
                m_cN(a) = N4(a);
                m_cD(a, 0) = D2(a, 0);
                m_cD(a, 1) = D2(a, 1);
                m_cD(a, 2) = 0.0;
            }
            V = quad::area(XL);
        }
        // Kabsch "gradients": reference local positions
        for (int a = 0; a < nn; ++a)
            for (int i = 0; i < 3; ++i)
                m_cgc(a, i) = X(i, a) - cc(i);
    }
    // centroid and reference local positions
    m_cc0.resize(3);
    m_cc0.Zero();
    for (int a = 0; a < nn; ++a)
        for (int i = 0; i < 3; ++i)
            m_cc0(i) += X(i, a) / static_cast<double>(nn);
    m_cY0.resize(nn, 3);
    for (int a = 0; a < nn; ++a)
        for (int i = 0; i < 3; ++i)
            m_cY0(a, i) = X(i, a) - m_cc0(i);
    m_cY0s.resize(3);
    for (int i = 0; i < 3; ++i)
        m_cY0s(i) = Xs(i) - m_cc0(i);
    // penalty scaled on the reference size
    m_ciK = m_corot_surf ? m_K * std::sqrt(V) : m_K * std::cbrt(V);
    m_corot_init = true;
}

void ASDEmbeddedNodeElement::corotComputeBg(Matrix& B, Vector& g)
{
    corotSetup();
    typedef ASDQuaternion<double> Q4;
    int nn = numRetained();
    // with -shearDeformable the retained nodes expose their rotations too, so
    // the reduced dofset interleaves [u(3) r(3)] per retained node
    int nr = m_shear ? 6 : 3;
    // -slip: the real rebar node closes the reduced dofset with its 3
    // translations and, when the tie is rotational, its 3 rotations; the slip
    // rows are appended after the 6 embedding rows
    int scols = m_slip ? (3 + (m_slip_rot ? 3 : 0)) : 0;
    int nsrows = scols;
    int ncols = 6 + nr * nn + scols;

    // current configuration (getGlobalDisplacements removes U0)
    const Vector& U = getGlobalDisplacements();
    static Matrix x;
    x.resize(3, nn);
    static Vector xs(3);
    int pos = m_nodes[0]->getNumberDOF();
    for (int a = 0; a < nn; ++a) {
        for (int i = 0; i < 3; ++i)
            x(i, a) = m_cY0(a, i) + m_cc0(i) + U(pos + i);
        pos += m_nodes[static_cast<std::size_t>(a + 1)]->getNumberDOF();
    }
    for (int i = 0; i < 3; ++i)
        xs(i) = m_cY0s(i) + m_cc0(i) + U(i);

    // frame: R = polar(F), F = sum_a x_a (x) g_a; centroid
    static Matrix F(3, 3);
    static Matrix R(3, 3);
    F.Zero();
    for (int a = 0; a < nn; ++a)
        for (int i = 0; i < 3; ++i)
            for (int j = 0; j < 3; ++j)
                F(i, j) += x(i, a) * m_cgc(a, j);
    if (m_corot_surf) {
        // F is the rank-2 Kabsch covariance of a flat patch: complete the rank
        // (A += gamma * n_cur (x) n_ref) so Higham's iteration can run; the
        // completed polar equals the SVD Kabsch rotation.
        static Vector t1(3), t2(3), nc(3);
        for (int i = 0; i < 3; ++i) {
            t1(i) = F(i, 0) * m_cE0(0, 0) + F(i, 1) * m_cE0(1, 0) + F(i, 2) * m_cE0(2, 0);
            t2(i) = F(i, 0) * m_cE0(0, 1) + F(i, 1) * m_cE0(1, 1) + F(i, 2) * m_cE0(2, 1);
        }
        cross(t1, t2, nc);
        nc.Normalize();
        double gam = 0.0;
        for (int i = 0; i < 3; ++i)
            for (int j = 0; j < 3; ++j)
                gam += F(i, j) * F(i, j);
        gam = std::sqrt(gam / 2.0);
        static Matrix Faug(3, 3);
        Faug = F;
        for (int i = 0; i < 3; ++i)
            for (int j = 0; j < 3; ++j)
                Faug(i, j) += gam * nc(i) * m_cE0(j, 2);
        polar3(Faug, R);
    }
    else {
        polar3(F, R);
    }
    static Vector c(3);
    c.Zero();
    for (int a = 0; a < nn; ++a)
        for (int i = 0; i < 3; ++i)
            c(i) += x(i, a) / static_cast<double>(nn);

    // local levers and deformational displacements
    static Matrix y;
    y.resize(nn, 3);
    static Vector ys(3);
    static Matrix v;
    v.resize(nn, 3);
    static Vector vs(3);
    for (int a = 0; a < nn; ++a)
        for (int i = 0; i < 3; ++i) {
            double s = 0.0;
            for (int k = 0; k < 3; ++k)
                s += R(k, i) * (x(k, a) - c(k));
            y(a, i) = s;
            v(a, i) = s - m_cY0(a, i);
        }
    for (int i = 0; i < 3; ++i) {
        double s = 0.0;
        for (int k = 0; k < 3; ++k)
            s += R(k, i) * (xs(k) - c(k));
        ys(i) = s;
        vs(i) = s - m_cY0s(i);
    }

    // g_u and omega (1/2 sum D_a x v_a)
    static Vector gu(3);
    static Vector om(3);
    gu.Zero(); om.Zero();
    for (int a = 0; a < nn; ++a) {
        for (int i = 0; i < 3; ++i)
            gu(i) += m_cN(a) * v(a, i);
        if (m_corot_surf) {
            // local components of the deformational displacement: w = E0^T v_a
            double w0 = 0.0, w1 = 0.0, w2 = 0.0;
            for (int i = 0; i < 3; ++i) {
                w0 += m_cE0(i, 0) * v(a, i);
                w1 += m_cE0(i, 1) * v(a, i);
                w2 += m_cE0(i, 2) * v(a, i);
            }
            // bending from the slope of the transverse deformational field
            // (with -shearDeformable it comes from the nodal rotations below),
            // drilling from the in-plane skew (the linear kernel, corotated)
            if (!m_shear) {
                om(0) += m_cD(a, 1) * w2;
                om(1) += -m_cD(a, 0) * w2;
            }
            om(2) += 0.5 * (m_cD(a, 0) * w1 - m_cD(a, 1) * w0);
        }
        else {
            om(0) += 0.5 * (m_cD(a, 1) * v(a, 2) - m_cD(a, 2) * v(a, 1));
            om(1) += 0.5 * (m_cD(a, 2) * v(a, 0) - m_cD(a, 0) * v(a, 2));
            om(2) += 0.5 * (m_cD(a, 0) * v(a, 1) - m_cD(a, 1) * v(a, 0));
        }
    }
    for (int i = 0; i < 3; ++i)
        gu(i) -= vs(i);

    // slave deformational rotation: theta = rotvec(R^T Rs), via quaternions
    // (well conditioned near pi, unlike the trace-based log)
    Q4 qR = Q4::FromRotationMatrix(R);
    Q4 qs(m_qs[0], m_qs[1], m_qs[2], m_qs[3]);
    Q4 qdef = qR.conjugate() * qs;
    qdef.normalize();
    static Vector theta(3);
    qdef.toRotationVector(theta(0), theta(1), theta(2));
    static Vector thloc(3);
    if (m_corot_surf) {
        for (int i = 0; i < 3; ++i) {
            double s = 0.0;
            for (int k = 0; k < 3; ++k)
                s += m_cE0(k, i) * theta(k);
            thloc(i) = s;
        }
    }
    else {
        thloc = theta;
    }

    // -shearDeformable: the bending rows read the N-weighted deformational
    // nodal rotations theta_a = rotvec(R^T Ra) of the host, in E0 components
    static Matrix ThA;
    if (m_shear) {
        ThA.resize(nn, 3);
        for (int a = 0; a < nn; ++a) {
            Q4 qa(m_qa[4 * a], m_qa[4 * a + 1], m_qa[4 * a + 2], m_qa[4 * a + 3]);
            Q4 qda = qR.conjugate() * qa;
            qda.normalize();
            double t0, t1, t2;
            qda.toRotationVector(t0, t1, t2);
            ThA(a, 0) = t0; ThA(a, 1) = t1; ThA(a, 2) = t2;
            for (int i = 0; i < 2; ++i)
                om(i) += m_cN(a) * (m_cE0(0, i) * t0 + m_cE0(1, i) * t1 + m_cE0(2, i) * t2);
        }
    }

    g.resize(6 + nsrows);
    g.Zero();
    for (int i = 0; i < 3; ++i) {
        g(i) = gu(i);
        g(3 + i) = om(i) - thloc(i);
    }

    // ---- exact first variation (see DESIGN.md; gated by verify_corot_embed.py)
    // dw = sum_b Gb R^T du_b, Gb = (tr(U)I - U)^-1 skew(g_b)
    static Matrix RtF(3, 3);
    RtF.addMatrixTransposeProduct(0.0, R, F, 1.0);
    static Matrix M(3, 3);
    static Matrix Minv(3, 3);
    double trU = RtF(0, 0) + RtF(1, 1) + RtF(2, 2);
    for (int i = 0; i < 3; ++i)
        for (int j = 0; j < 3; ++j)
            M(i, j) = (i == j ? trU : 0.0) - RtF(i, j);
    inv3(M, Minv);
    static Matrix DW;
    DW.resize(3, ncols);
    DW.Zero();
    for (int a = 0; a < nn; ++a) {
        // Gb = Minv * skew(g_a); the a-th DW block is Gb * R^T
        static Matrix Sk(3, 3);
        Sk.Zero();
        Sk(0, 1) = -m_cgc(a, 2); Sk(0, 2) = m_cgc(a, 1);
        Sk(1, 0) = m_cgc(a, 2);  Sk(1, 2) = -m_cgc(a, 0);
        Sk(2, 0) = -m_cgc(a, 1); Sk(2, 1) = m_cgc(a, 0);
        static Matrix Gb(3, 3);
        Gb.addMatrixProduct(0.0, Minv, Sk, 1.0);
        for (int i = 0; i < 3; ++i)
            for (int j = 0; j < 3; ++j) {
                double s = 0.0;
                for (int k = 0; k < 3; ++k)
                    s += Gb(i, k) * R(j, k);
                DW(i, 6 + nr * a + j) = s;
            }
    }

    B.resize(6 + nsrows, ncols);
    B.Zero();
    // dg_u rows
    static Vector z(3);
    z.Zero();
    for (int a = 0; a < nn; ++a)
        for (int i = 0; i < 3; ++i)
            z(i) += m_cN(a) * y(a, i);
    for (int i = 0; i < 3; ++i)
        z(i) -= ys(i);
    for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < 3; ++j)
            B(i, j) = -R(j, i);                       // -R^T on slave u
        for (int a = 0; a < nn; ++a)
            for (int j = 0; j < 3; ++j)
                B(i, 6 + nr * a + j) = m_cN(a) * R(j, i);
    }
    // + skew(z) DW
    for (int col = 0; col < ncols; ++col) {
        B(0, col) += -z(2) * DW(1, col) + z(1) * DW(2, col);
        B(1, col) += z(2) * DW(0, col) - z(0) * DW(2, col);
        B(2, col) += -z(1) * DW(0, col) + z(0) * DW(1, col);
    }
    // d omega rows. OP is the row operator on dv_a: volumetric 1/2 skew(D_a),
    // or the surface slope/drilling operator premultiplied by E0^T components.
    static Matrix LEV(3, 3);
    LEV.Zero();
    // -shearDeformable: SHW accumulates sum_a N_a E0^T Tinv_a, the -dw part of
    // d theta_a = Tinv_a (R^T dphi_a - dw), folded into the frame term below
    static Matrix SHW(3, 3);
    SHW.Zero();
    static Matrix OPa(3, 3);
    static Matrix SY(3, 3);
    for (int a = 0; a < nn; ++a) {
        OPa.Zero();
        if (m_corot_surf) {
            // om_local = OPloc * (E0^T dv_a): rows [D_y w2; -D_x w2; skew/2];
            // with -shearDeformable the bending rows do not read dv_a
            static Matrix OPloc(3, 3);
            OPloc.Zero();
            if (!m_shear) {
                OPloc(0, 2) = m_cD(a, 1);
                OPloc(1, 2) = -m_cD(a, 0);
            }
            OPloc(2, 0) = -0.5 * m_cD(a, 1);
            OPloc(2, 1) = 0.5 * m_cD(a, 0);
            // OPa = OPloc * E0^T
            for (int i = 0; i < 3; ++i)
                for (int j = 0; j < 3; ++j) {
                    double s = 0.0;
                    for (int k = 0; k < 3; ++k)
                        s += OPloc(i, k) * m_cE0(j, k);
                    OPa(i, j) = s;
                }
        }
        else {
            OPa(0, 1) = -0.5 * m_cD(a, 2); OPa(0, 2) = 0.5 * m_cD(a, 1);
            OPa(1, 0) = 0.5 * m_cD(a, 2);  OPa(1, 2) = -0.5 * m_cD(a, 0);
            OPa(2, 0) = -0.5 * m_cD(a, 1); OPa(2, 1) = 0.5 * m_cD(a, 0);
        }
        // direct part: OPa * R^T on the a-th u block
        for (int i = 0; i < 3; ++i)
            for (int j = 0; j < 3; ++j) {
                double s = 0.0;
                for (int k = 0; k < 3; ++k)
                    s += OPa(i, k) * R(j, k);
                B(3 + i, 6 + nr * a + j) += s;
            }
        // lever part: LEV += OPa * skew(y_a)
        SY.Zero();
        SY(0, 1) = -y(a, 2); SY(0, 2) = y(a, 1);
        SY(1, 0) = y(a, 2);  SY(1, 2) = -y(a, 0);
        SY(2, 0) = -y(a, 1); SY(2, 1) = y(a, 0);
        LEV.addMatrixProduct(1.0, OPa, SY, 1.0);
        // -shearDeformable bending rows: + N_a (E0^T Tinv_a) R^T on the
        // rotation columns of node a (the direct part of d theta_a), and the
        // N_a E0^T Tinv_a weight into SHW for the shared -dw part
        if (m_shear) {
            static Vector tha(3);
            for (int i = 0; i < 3; ++i)
                tha(i) = ThA(a, i);
            static Matrix Tia(3, 3);
            leftJacobianInv(tha, Tia);
            static Matrix Wa(3, 3);
            for (int i = 0; i < 3; ++i)
                for (int j = 0; j < 3; ++j) {
                    double s = 0.0;
                    for (int k = 0; k < 3; ++k)
                        s += m_cE0(k, i) * Tia(k, j);
                    Wa(i, j) = s;
                }
            for (int i = 0; i < 2; ++i) {
                for (int j = 0; j < 3; ++j) {
                    double s = 0.0;
                    for (int k = 0; k < 3; ++k)
                        s += Wa(i, k) * R(j, k);
                    B(3 + i, 6 + nr * a + 3 + j) += m_cN(a) * s;
                }
                for (int k = 0; k < 3; ++k)
                    SHW(i, k) += m_cN(a) * Wa(i, k);
            }
        }
    }
    // frame-variation terms: +LEV*DW from d(omega), +E*Tinv*DW from
    // -d(theta_local) = -E (Tinv (R^T dphi_s - dw)), with E = E0^T on
    // surfaces and I on solids
    static Matrix Tinv(3, 3);
    leftJacobianInv(theta, Tinv);
    static Matrix ETi(3, 3);
    if (m_corot_surf) {
        for (int i = 0; i < 3; ++i)
            for (int j = 0; j < 3; ++j) {
                double s = 0.0;
                for (int k = 0; k < 3; ++k)
                    s += m_cE0(k, i) * Tinv(k, j);
                ETi(i, j) = s;
            }
    }
    else {
        ETi = Tinv;
    }
    for (int i = 0; i < 3; ++i)
        for (int col = 0; col < ncols; ++col) {
            double s = 0.0;
            for (int k = 0; k < 3; ++k)
                s += (LEV(i, k) + ETi(i, k) - SHW(i, k)) * DW(k, col);
            B(3 + i, col) += s;
        }
    // -E Tinv R^T on the slave rotation columns
    for (int i = 0; i < 3; ++i)
        for (int j = 0; j < 3; ++j) {
            double s = 0.0;
            for (int k = 0; k < 3; ++k)
                s += ETi(i, k) * R(j, k);
            B(3 + i, 3 + j) -= s;
        }

    // ---- -slip rows: the zeroLength kernel corotated by the host frame ----
    if (m_slip) {
        int c_r = 6 + nr * nn; // first column of the real-node block
        // relative displacement AUX - real (deformational: the two nodes are
        // coincident in the reference configuration)
        int rpos = m_num_dofs - m_nodes.back()->getNumberDOF();
        static Vector d(3);
        for (int i = 0; i < 3; ++i)
            d(i) = U(i) - U(rpos + i);
        // translational rows: g_i = a_i . d with a_i = R * x0_i the corotated
        // triad. DW is the LOCAL (right-trivialized) spin of the frame,
        // dR = R skew(w) with w = DW dq (the same convention the skew(z) DW
        // and Tinv DW terms above rely on), so the frame variation is
        // d(a_i . d) = (R (w x x0_i)) . d = (x0_i x R^T d) . w and the host
        // columns get +(x0_i x R^T d)^T DW -- the geometric stiffness of the
        // rotating slip direction. Gated by verify_slip_corot.py.
        static Vector ai(3);
        static Vector dl(3);
        static Vector x0i(3);
        static Vector xxd(3);
        for (int k = 0; k < 3; ++k) {
            double s = 0.0;
            for (int j = 0; j < 3; ++j)
                s += R(j, k) * d(j);
            dl(k) = s;                       // R^T d
        }
        for (int i = 0; i < 3; ++i) {
            for (int k = 0; k < 3; ++k) {
                double s = 0.0;
                for (int j = 0; j < 3; ++j)
                    s += R(k, j) * m_slip_T0(i, j);
                ai(k) = s;
            }
            if (i == 0)
                for (int k = 0; k < 3; ++k)
                    m_slip_axis(k) = ai(k);
            g(6 + i) = ai ^ d;
            for (int j = 0; j < 3; ++j) {
                B(6 + i, j) += ai(j);        // AUX (slave) translations
                B(6 + i, c_r + j) -= ai(j);  // real node translations
            }
            for (int j = 0; j < 3; ++j)
                x0i(j) = m_slip_T0(i, j);
            cross(x0i, dl, xxd);
            for (int col = 0; col < ncols; ++col)
                B(6 + i, col) += xxd(0) * DW(0, col) + xxd(1) * DW(1, col) + xxd(2) * DW(2, col);
        }
        // rotational tie rows: theta = rotvec(Rs^T Rr) is invariant under a
        // superposed rigid rotation, so there is no DW coupling; the exact
        // variation is d(theta) = Tinv(theta) Rs^T (dphi_r - dphi_s)
        if (m_slip_rot) {
            Q4 qr(m_qr[0], m_qr[1], m_qr[2], m_qr[3]);
            Q4 qrel = qs.conjugate() * qr;
            qrel.normalize();
            static Vector threl(3);
            qrel.toRotationVector(threl(0), threl(1), threl(2));
            for (int i = 0; i < 3; ++i)
                g(9 + i) = threl(i);
            static Matrix Rs(3, 3);
            qs.toRotationMatrix(Rs);
            static Matrix Tir(3, 3);
            leftJacobianInv(threl, Tir);
            for (int i = 0; i < 3; ++i)
                for (int j = 0; j < 3; ++j) {
                    double s = 0.0;
                    for (int k = 0; k < 3; ++k)
                        s += Tir(i, k) * Rs(j, k); // Tinv * Rs^T
                    B(9 + i, 3 + j) -= s;          // AUX (slave) rotations
                    B(9 + i, c_r + 3 + j) += s;    // real node rotations
                }
        }
        // recorder buffers
        for (int i = 0; i < nsrows; ++i)
            m_slip_g(i) = g(6 + i);
    }
}

const Matrix& ASDEmbeddedNodeElement::embedLocalStiffness()
{
    // the constraint mode is shared by every family
    int mode = m_rot_c ? Mode_UR : (m_up ? Mode_UP : Mode_U);

    // isoparametric (non-simplex) hosts go through the generic path
    if (m_family == Fam_Quad)
        return QUAD_2D(mode);
    if (m_family == Fam_Quad3D)
        return QUAD_3D(mode);
    if (m_family == Fam_Hexa)
        return HEX_3D(mode);
    // simplex hosts keep the original, untouched kernels
    if (numRetained() == 3) {
        // support shape is a triangle ...
        if (m_ndm == 2) {
            // ... in 2D
            if (m_rot_c) {
                // ... with rotational dofs
                return TRI_2D_UR();
            }
            else if (m_up) {
                // ... with pressure dofs
                return TRI_2D_UP();
            }
            else {
                // ... without rotational dofs
                return TRI_2D_U();
            }
        }
        else {
            // ... in 3D
            if (m_rot_c) {
                // ... with rotational dofs
                return TRI_3D_UR();
            }
            else if(m_up) {
                // ... with pressure dofs
                return TRI_3D_UP();
            }
            else {
                // ... without rotational dofs
                return TRI_3D_U();
            }
        }
    }
    else {
        // support shape is a tetrahedron ...
        if (m_rot_c) {
            // ... with rotational dofs
            return TET_3D_UR();
        }
        else if (m_up) {
            // ... with pressure dofs
            return TET_3D_UP();
        }
        else {
            // ... without rotational dofs
            return TET_3D_U();
        }
    }
}

const Matrix& ASDEmbeddedNodeElement::getTangentStiff()
{
    // -corotational: K = B^T diag(w) B on the exact first variation
    // (Gauss-Newton; the k*g*d2g geometric term vanishes with the constraint
    // violation). Without -slip every row weighs iK and this is the original
    // iK * B^T B; the slip rows weigh the material tangent (row 6) and the
    // raw rigid-tie stiffness KS (the others).
    if (m_corot) {
        static Matrix B;
        static Vector g;
        corotComputeBg(B, g);
        int nred = B.noCols();
        int nrows = B.noRows();
        static Matrix WB;
        WB.resize(nrows, nred);
        for (int i = 0; i < nrows; ++i) {
            double w = (i < 6) ? m_ciK :
                (i == 6 ? m_slip_mat->getTangent() : m_KS);
            for (int j = 0; j < nred; ++j)
                WB(i, j) = w * B(i, j);
        }
        static Matrix KL;
        KL.resize(nred, nred);
        KL.addMatrixTransposeProduct(0.0, B, WB, 1.0);
        static Matrix K;
        K.resize(m_num_dofs, m_num_dofs);
        K.Zero();
        for (int i = 0; i < nred; ++i) {
            int ig = m_mapping(i);
            for (int j = 0; j < nred; ++j)
                K(ig, m_mapping(j)) = KL(i, j);
        }
        return K;
    }

    // compute stiffness matrix in reduced local dofset
    const Matrix& KL = embedLocalStiffness();

    // output matrix
    static Matrix K;
    K.resize(m_num_dofs, m_num_dofs);
    K.Zero();

    // copy in global dofset
    for (int i = 0; i < KL.noRows(); ++i) {
        int ig = m_mapping(i);
        for (int j = 0; j < KL.noCols(); ++j) {
            int jg = m_mapping(j);
            K(ig, jg) = KL(i, j);
        }
    }

    // -slip: the zeroLength-equivalent block on the frozen triad
    if (m_slip)
        slipAddLinear(&K, nullptr);

    // done
    return K;
}

const Matrix& ASDEmbeddedNodeElement::getInitialStiff()
{
    return getTangentStiff();
}

const Matrix& ASDEmbeddedNodeElement::getMass()
{
    static Matrix M;
    M.resize(m_num_dofs, m_num_dofs);
    M.Zero();
    return M;
}

const Matrix& ASDEmbeddedNodeElement::getDamp()
{
    static Matrix C;
    C.resize(m_num_dofs, m_num_dofs);
    C.Zero();
    return C;
}

int
ASDEmbeddedNodeElement::addInertiaLoadToUnbalance(const Vector& accel)
{
    return 0;
}

const Vector& ASDEmbeddedNodeElement::getResistingForce()
{
    static Vector F;
    F.resize(m_num_dofs);
    // -corotational: f = B^T q, exact at any rotation (and exactly
    // self-equilibrated: the closed-form B annihilates the rigid modes).
    // Without -slip q = iK * g and this is the original iK * B^T g; the slip
    // row carries the material stress, the other tie rows KS * g.
    if (m_corot) {
        static Matrix B;
        static Vector g;
        corotComputeBg(B, g);
        int nred = B.noCols();
        int nrows = B.noRows();
        static Vector q;
        q.resize(nrows);
        for (int i = 0; i < nrows; ++i)
            q(i) = (i < 6) ? m_ciK * g(i) :
                (i == 6 ? m_slip_mat->getStress() : m_KS * g(i));
        static Vector fr;
        fr.resize(nred);
        fr.addMatrixTransposeVector(0.0, B, q, 1.0);
        F.Zero();
        for (int i = 0; i < nred; ++i)
            F(m_mapping(i)) = fr(i);
        return F;
    }
    if (m_slip) {
        // the slip row is materially non-linear: the force is NOT K * U.
        // Embedding part (pure elastic penalty): F = K_e * U; slip block
        // appended with the material stress on the bar axis.
        const Matrix& KL = embedLocalStiffness();
        static Matrix KE;
        KE.resize(m_num_dofs, m_num_dofs);
        KE.Zero();
        for (int i = 0; i < KL.noRows(); ++i) {
            int ig = m_mapping(i);
            for (int j = 0; j < KL.noCols(); ++j)
                KE(ig, m_mapping(j)) = KL(i, j);
        }
        const Vector& U = getGlobalDisplacements();
        F.addMatrixVector(0.0, KE, U, 1.0);
        slipAddLinear(nullptr, &F);
        return F;
    }
    const Matrix& K = getTangentStiff();
    const Vector& U = getGlobalDisplacements();
    F.addMatrixVector(0.0, K, U, 1.0);
    return F;
}

const Vector& ASDEmbeddedNodeElement::getResistingForceIncInertia()
{
    return getResistingForce();
}

int ASDEmbeddedNodeElement::sendSelf(int commitTag, Channel& theChannel)
{
    // Buffers are sized at run time: with an 8-node host in u-p mode the
    // element reaches 4 + 8*4 = 36 dofs, which the previous fixed ID(35) /
    // Vector(32) could not hold. Same layout as ASDConstraintEquationElement:
    // a small header first, then the variable-length payloads.
    int res = 0;
    int dataTag = getDbTag();

    int NN = m_node_ids.Size();          // 1 constrained + n retained (+ the -slip real node)
    int NRET = numRetained();
    int NMAP = m_mapping.Size();
    int NU0 = m_U0_computed ? m_U0.Size() : 0;

    // INT data 1: header with every size needed to read the rest
    static ID idData1(21);
    idData1(0) = getTag();
    idData1(1) = NN;
    idData1(2) = NMAP;
    idData1(3) = NU0;
    idData1(4) = m_ndm;
    idData1(5) = m_num_dofs;
    idData1(6) = m_rot_c_flag ? 1 : 0;
    idData1(7) = m_rot_c ? 1 : 0;
    idData1(8) = m_p_flag ? 1 : 0;
    idData1(9) = m_up ? 1 : 0;
    idData1(10) = m_U0_computed ? 1 : 0;
    idData1(11) = m_family;
    // the -shape request must travel too: setDomain runs again on the receiving
    // side and calls resolveFamily, which cannot tell a quadrilateral face from
    // a tetrahedron without it
    idData1(12) = m_shape_request;
    idData1(13) = m_shear_flag ? 1 : 0;
    idData1(14) = m_shear ? 1 : 0;
    idData1(15) = m_corot_flag ? 1 : 0;
    idData1(16) = m_corot ? 1 : 0;
    // -slip: the flag must travel BEFORE the node list is interpreted (the
    // real node is the last id, and numRetained depends on it), and the
    // material class/db tags are needed to rebuild the law on the far side
    idData1(17) = m_slip ? 1 : 0;
    idData1(18) = m_slip_rot ? 1 : 0;
    idData1(19) = 0;
    idData1(20) = 0;
    if (m_slip) {
        int matDbTag = m_slip_mat->getDbTag();
        if (matDbTag == 0) {
            matDbTag = theChannel.getDbTag();
            if (matDbTag != 0)
                m_slip_mat->setDbTag(matDbTag);
        }
        idData1(19) = m_slip_mat->getClassTag();
        idData1(20) = matDbTag;
    }
    res = theChannel.sendID(dataTag, commitTag, idData1);
    if (res < 0) {
        opserr << "WARNING ASDEmbeddedNodeElement::sendSelf() - " << this->getTag() << " failed to send ID 1\n";
        return res;
    }

    // INT data 2: node ids followed by the dof mapping
    ID idData2(NN + NMAP);
    int pos = 0;
    for (int i = 0; i < NN; ++i)
        idData2(pos++) = m_node_ids(i);
    for (int i = 0; i < NMAP; ++i)
        idData2(pos++) = m_mapping(i);
    res = theChannel.sendID(dataTag, commitTag, idData2);
    if (res < 0) {
        opserr << "WARNING ASDEmbeddedNodeElement::sendSelf() - " << this->getTag() << " failed to send ID 2\n";
        return res;
    }

    // DOUBLE data: K, KP, the corotational slave-rotation state, the
    // per-retained-node rotation state (-corotational + -shearDeformable:
    // 14 doubles per retained node), the -slip data (KS, x0, and the real
    // node rotation state when the corotational tie is rotational), and the
    // initial displacement vector
    int NQA = (m_corot && m_shear) ? 14 * NRET : 0;
    int NSLIP = m_slip ? 4 : 0;
    int NQR = (m_slip && m_corot && m_slip_rot) ? 14 : 0;
    Vector vectData(16 + NQA + NSLIP + NQR + NU0);
    pos = 0;
    vectData(pos++) = m_K;
    vectData(pos++) = m_KP;
    for (int i = 0; i < 4; ++i) vectData(pos++) = m_qs[i];
    for (int i = 0; i < 3; ++i) vectData(pos++) = m_rv[i];
    for (int i = 0; i < 4; ++i) vectData(pos++) = m_qs_conv[i];
    for (int i = 0; i < 3; ++i) vectData(pos++) = m_rv_conv[i];
    if (NQA > 0) {
        for (int a = 0; a < NRET; ++a) {
            for (int i = 0; i < 4; ++i) vectData(pos++) = m_qa[4 * a + i];
            for (int i = 0; i < 3; ++i) vectData(pos++) = m_rva[3 * a + i];
            for (int i = 0; i < 4; ++i) vectData(pos++) = m_qa_conv[4 * a + i];
            for (int i = 0; i < 3; ++i) vectData(pos++) = m_rva_conv[3 * a + i];
        }
    }
    if (m_slip) {
        vectData(pos++) = m_KS;
        for (int i = 0; i < 3; ++i) vectData(pos++) = m_slip_x0(i);
        if (NQR > 0) {
            for (int i = 0; i < 4; ++i) vectData(pos++) = m_qr[i];
            for (int i = 0; i < 3; ++i) vectData(pos++) = m_rvr[i];
            for (int i = 0; i < 4; ++i) vectData(pos++) = m_qr_conv[i];
            for (int i = 0; i < 3; ++i) vectData(pos++) = m_rvr_conv[i];
        }
    }
    for (int i = 0; i < NU0; ++i)
        vectData(pos++) = m_U0(i);
    res = theChannel.sendVector(dataTag, commitTag, vectData);
    if (res < 0) {
        opserr << "WARNING ASDEmbeddedNodeElement::sendSelf() - " << this->getTag() << " failed to send Vector\n";
        return res;
    }

    // -slip: the material state travels with its own sendSelf
    if (m_slip) {
        res = m_slip_mat->sendSelf(commitTag, theChannel);
        if (res < 0) {
            opserr << "WARNING ASDEmbeddedNodeElement::sendSelf() - " << this->getTag() << " failed to send the -slip material\n";
            return res;
        }
    }

    // done
    return res;
}

int ASDEmbeddedNodeElement::recvSelf(int commitTag, Channel& theChannel, FEM_ObjectBroker& theBroker)
{
    int res = 0;
    int dataTag = this->getDbTag();

    // INT data 1: header
    static ID idData1(21);
    res = theChannel.recvID(dataTag, commitTag, idData1);
    if (res < 0) {
        opserr << "WARNING ASDEmbeddedNodeElement::recvSelf() - " << this->getTag() << " failed to receive ID 1\n";
        return res;
    }
    setTag(idData1(0));
    int NN = idData1(1);
    int NMAP = idData1(2);
    int NU0 = idData1(3);
    m_ndm = idData1(4);
    m_num_dofs = idData1(5);
    m_rot_c_flag = idData1(6) == 1;
    m_rot_c = idData1(7) == 1;
    m_p_flag = idData1(8) == 1;
    m_up = idData1(9) == 1;
    m_U0_computed = idData1(10) == 1;
    m_family = idData1(11);
    m_shape_request = idData1(12);
    m_shear_flag = idData1(13) == 1;
    m_shear = idData1(14) == 1;
    m_corot_flag = idData1(15) == 1;
    m_corot = idData1(16) == 1;
    m_corot_init = false; // reference data is recomputed lazily from coords + U0
    m_slip = idData1(17) == 1;
    m_slip_rot = idData1(18) == 1;
    if (m_slip) {
        int matClassTag = idData1(19);
        if (m_slip_mat == nullptr || m_slip_mat->getClassTag() != matClassTag) {
            if (m_slip_mat)
                delete m_slip_mat;
            m_slip_mat = theBroker.getNewUniaxialMaterial(matClassTag);
            if (m_slip_mat == nullptr) {
                opserr << "WARNING ASDEmbeddedNodeElement::recvSelf() - " << this->getTag()
                    << " failed to create the -slip material with class tag " << matClassTag << "\n";
                return -1;
            }
        }
        m_slip_mat->setDbTag(idData1(20));
    }
    int NRET = NN - 1 - (m_slip ? 1 : 0);

    // INT data 2: node ids and dof mapping
    ID idData2(NN + NMAP);
    res = theChannel.recvID(dataTag, commitTag, idData2);
    if (res < 0) {
        opserr << "WARNING ASDEmbeddedNodeElement::recvSelf() - " << this->getTag() << " failed to receive ID 2\n";
        return res;
    }
    m_node_ids.resize(NN);
    m_nodes.assign(static_cast<std::size_t>(NN), nullptr);
    int pos = 0;
    for (int i = 0; i < NN; ++i)
        m_node_ids(i) = idData2(pos++);
    m_mapping.resize(NMAP);
    for (int i = 0; i < NMAP; ++i)
        m_mapping(i) = idData2(pos++);

    // DOUBLE data
    int NQA = (m_corot && m_shear) ? 14 * NRET : 0;
    int NSLIP = m_slip ? 4 : 0;
    int NQR = (m_slip && m_corot && m_slip_rot) ? 14 : 0;
    Vector vectData(16 + NQA + NSLIP + NQR + NU0);
    res = theChannel.recvVector(dataTag, commitTag, vectData);
    if (res < 0) {
        opserr << "WARNING ASDEmbeddedNodeElement::recvSelf() - " << this->getTag() << " failed to receive Vector\n";
        return res;
    }
    pos = 0;
    m_K = vectData(pos++);
    m_KP = vectData(pos++);
    for (int i = 0; i < 4; ++i) m_qs[i] = vectData(pos++);
    for (int i = 0; i < 3; ++i) m_rv[i] = vectData(pos++);
    for (int i = 0; i < 4; ++i) m_qs_conv[i] = vectData(pos++);
    for (int i = 0; i < 3; ++i) m_rv_conv[i] = vectData(pos++);
    if (NQA > 0) {
        m_qa.resize(4 * static_cast<std::size_t>(NRET));
        m_rva.resize(3 * static_cast<std::size_t>(NRET));
        m_qa_conv.resize(4 * static_cast<std::size_t>(NRET));
        m_rva_conv.resize(3 * static_cast<std::size_t>(NRET));
        for (int a = 0; a < NRET; ++a) {
            for (int i = 0; i < 4; ++i) m_qa[4 * a + i] = vectData(pos++);
            for (int i = 0; i < 3; ++i) m_rva[3 * a + i] = vectData(pos++);
            for (int i = 0; i < 4; ++i) m_qa_conv[4 * a + i] = vectData(pos++);
            for (int i = 0; i < 3; ++i) m_rva_conv[3 * a + i] = vectData(pos++);
        }
    }
    if (m_slip) {
        m_KS = vectData(pos++);
        m_slip_x0.resize(3);
        for (int i = 0; i < 3; ++i) m_slip_x0(i) = vectData(pos++);
        if (NQR > 0) {
            for (int i = 0; i < 4; ++i) m_qr[i] = vectData(pos++);
            for (int i = 0; i < 3; ++i) m_rvr[i] = vectData(pos++);
            for (int i = 0; i < 4; ++i) m_qr_conv[i] = vectData(pos++);
            for (int i = 0; i < 3; ++i) m_rvr_conv[i] = vectData(pos++);
        }
    }
    if (NU0 > 0) {
        m_U0.resize(NU0);
        for (int i = 0; i < NU0; ++i)
            m_U0(i) = vectData(pos++);
    }

    // -slip: the material state travels with its own recvSelf
    if (m_slip) {
        res = m_slip_mat->recvSelf(commitTag, theChannel, theBroker);
        if (res < 0) {
            opserr << "WARNING ASDEmbeddedNodeElement::recvSelf() - " << this->getTag() << " failed to receive the -slip material\n";
            return res;
        }
    }

    // done
    return res;
}

Response* ASDEmbeddedNodeElement::setResponse(const char** argv, int argc, OPS_Stream& output)
{
    // every response belongs to the -slip machinery: without it the element
    // is a pure penalty constraint and has nothing to record
    if (!m_slip || argc < 1)
        return Element::setResponse(argv, argc, output);

    Response* theResponse = nullptr;
    output.tag("ElementOutput");
    output.attr("eleType", this->getClassType());
    output.attr("eleTag", this->getTag());

    if (strcmp(argv[0], "slip") == 0) {
        // the scalar slip: the strain of the tau-slip law
        output.tag("ResponseType", "slip");
        theResponse = new ElementResponse(this, 1, Vector(1));
    }
    else if (strcmp(argv[0], "slipForce") == 0 || strcmp(argv[0], "bondForce") == 0) {
        // the scalar bond force: the stress of the tau-slip law
        output.tag("ResponseType", "slipForce");
        theResponse = new ElementResponse(this, 2, Vector(1));
    }
    else if (strcmp(argv[0], "gap") == 0) {
        // local relative displacement AUX - real: [slip, t1, (t2)]
        output.tag("ResponseType", "slip");
        for (int i = 1; i < m_ndm; ++i)
            output.tag("ResponseType", i == 1 ? "t1" : "t2");
        theResponse = new ElementResponse(this, 3, Vector(m_ndm));
    }
    else if (strcmp(argv[0], "gapForce") == 0) {
        // conjugate local forces: [bond force, KS*t1, (KS*t2)]
        output.tag("ResponseType", "slipForce");
        for (int i = 1; i < m_ndm; ++i)
            output.tag("ResponseType", i == 1 ? "Ft1" : "Ft2");
        theResponse = new ElementResponse(this, 4, Vector(m_ndm));
    }
    else if (strcmp(argv[0], "slipAxis") == 0) {
        // current bar axis (rotated by the host frame with -corotational)
        output.tag("ResponseType", "x1");
        output.tag("ResponseType", "x2");
        output.tag("ResponseType", "x3");
        theResponse = new ElementResponse(this, 5, Vector(3));
    }
    else if (strcmp(argv[0], "slipMaterial") == 0 && argc > 1) {
        // forward to the tau-slip law (e.g. stress/strain/tangent of the
        // inner material of a Parallel wrapper)
        theResponse = m_slip_mat->setResponse(&argv[1], argc - 1, output);
    }

    output.endTag();
    return theResponse;
}

int ASDEmbeddedNodeElement::getResponse(int responseID, Information& eleInfo)
{
    if (!m_slip)
        return Element::getResponse(responseID, eleInfo);

    static Vector r1(1);
    switch (responseID) {
    case 1:
        r1(0) = m_slip_mat->getStrain();
        return eleInfo.setVector(r1);
    case 2:
        r1(0) = m_slip_mat->getStress();
        return eleInfo.setVector(r1);
    case 3: {
        slipComputeGap();
        static Vector rg;
        rg.resize(m_ndm);
        for (int i = 0; i < m_ndm; ++i)
            rg(i) = m_slip_g(i);
        return eleInfo.setVector(rg);
    }
    case 4: {
        slipComputeGap();
        static Vector rf;
        rf.resize(m_ndm);
        rf(0) = m_slip_mat->getStress();
        for (int i = 1; i < m_ndm; ++i)
            rf(i) = m_KS * m_slip_g(i);
        return eleInfo.setVector(rf);
    }
    case 5: {
        slipComputeGap();
        static Vector ra(3);
        for (int i = 0; i < 3; ++i)
            ra(i) = m_slip_axis(i);
        return eleInfo.setVector(ra);
    }
    default:
        return Element::getResponse(responseID, eleInfo);
    }
}

const Vector& ASDEmbeddedNodeElement::getGlobalDisplacements() const
{
    static Vector U;
    U.resize(m_num_dofs);
    int counter = 0;
    for (Node* node : m_nodes) {
        const Vector& iu = node->getTrialDisp();
        for (int i = 0; i < iu.Size(); ++i) {
            U(counter++) = iu(i);
        }
    }
    if (m_U0_computed) {
        U.addVector(1.0, m_U0, -1.0);
    }
    return U;
}

double ASDEmbeddedNodeElement::slipComputeGap()
{
    if (m_corot) {
        // the corotational kernel computes the slip rows (and refreshes the
        // recorder buffers) as part of B and g
        static Matrix B;
        static Vector g;
        corotComputeBg(B, g);
        return m_slip_g(0);
    }
    // frozen frame: exactly the zeroLength kinematics, gap = u_AUX - u_real
    // on the deformational displacements
    const Vector& U = getGlobalDisplacements();
    int rpos = m_num_dofs - m_nodes.back()->getNumberDOF();
    for (int i = 0; i < m_ndm; ++i) {
        double s = 0.0;
        for (int j = 0; j < m_ndm; ++j)
            s += m_slip_T0(i, j) * (U(j) - U(rpos + j));
        m_slip_g(i) = s;
    }
    if (m_slip_rot) {
        int nrot = (m_ndm == 2) ? 1 : 3;
        for (int i = 0; i < nrot; ++i)
            m_slip_g(m_ndm + i) = U(m_ndm + i) - U(rpos + m_ndm + i);
    }
    for (int j = 0; j < 3; ++j)
        m_slip_axis(j) = m_slip_x0(j);
    return m_slip_g(0);
}

void ASDEmbeddedNodeElement::slipAddLinear(Matrix* K, Vector* F)
{
    // zeroLength-equivalent block on the frozen triad: the slip law on the
    // bar axis, the raw rigid-tie stiffness KS on every other relative dof.
    // Written directly in the FULL element dofset: the AUX dofs start at 0,
    // the real node dofs at rpos, and no other dof is touched.
    int rpos = m_num_dofs - m_nodes.back()->getNumberDOF();
    int nrot = m_slip_rot ? ((m_ndm == 2) ? 1 : 3) : 0;
    if (K) {
        // C = k_mat * x (x) x + KS * (t_r (x) t_r) on the translations
        static Matrix C(3, 3);
        C.Zero();
        double kmat = m_slip_mat->getTangent();
        for (int r = 0; r < m_ndm; ++r) {
            double kr = (r == 0) ? kmat : m_KS;
            for (int i = 0; i < m_ndm; ++i)
                for (int j = 0; j < m_ndm; ++j)
                    C(i, j) += kr * m_slip_T0(r, i) * m_slip_T0(r, j);
        }
        for (int i = 0; i < m_ndm; ++i)
            for (int j = 0; j < m_ndm; ++j) {
                (*K)(i, j) += C(i, j);
                (*K)(i, rpos + j) -= C(i, j);
                (*K)(rpos + i, j) -= C(i, j);
                (*K)(rpos + i, rpos + j) += C(i, j);
            }
        for (int i = 0; i < nrot; ++i) {
            int a = m_ndm + i;
            int r = rpos + m_ndm + i;
            (*K)(a, a) += m_KS;
            (*K)(a, r) -= m_KS;
            (*K)(r, a) -= m_KS;
            (*K)(r, r) += m_KS;
        }
    }
    if (F) {
        // refresh the local gap (a recorder can ask for forces outside the
        // update sequence) and assemble f = B^T q with q = [sigma, KS * g...]
        slipComputeGap();
        static Vector fv(3);
        fv.Zero();
        for (int r = 0; r < m_ndm; ++r) {
            double qr = (r == 0) ? m_slip_mat->getStress() : m_KS * m_slip_g(r);
            for (int i = 0; i < m_ndm; ++i)
                fv(i) += qr * m_slip_T0(r, i);
        }
        for (int i = 0; i < m_ndm; ++i) {
            (*F)(i) += fv(i);
            (*F)(rpos + i) -= fv(i);
        }
        for (int i = 0; i < nrot; ++i) {
            double m = m_KS * m_slip_g(m_ndm + i);
            (*F)(m_ndm + i) += m;
            (*F)(rpos + m_ndm + i) -= m;
        }
    }
}

const Matrix& ASDEmbeddedNodeElement::TRI_2D_U()
{
    // output
    static Matrix K(8, 8);

    // collect triangle coordinates
    static Matrix X(2, 3);
    for (int i = 1; i < 4; i++) {
        const Node* node = m_nodes[static_cast<std::size_t>(i)];
        X(0, i-1) = node->getCrds()(0);
        X(1, i-1) = node->getCrds()(1);
    }

    // shape functions natural derivatives
    static Matrix dN(3, 2);
    tri::shapeFunDer(dN);

    // jacobian
    static Matrix J(2, 2);
    J.addMatrixProduct(0.0, X, dN, 1.0);
    double detJ = det2(J);
    double V = detJ / 2.0;
    static Matrix invJ(2, 2);
    J.Invert(invJ);

    // find local coordinates of constrained node
    double lx, ly;
    tri::localCoord(X, invJ, m_nodes[0]->getCrds()(0), m_nodes[0]->getCrds()(1), lx, ly);

    // compute shape functions at constrained node
    static Vector N(3);
    for(int i = 0; i < 3; ++i)
        N(i) = tri::shapeFun(lx, ly, i);

    // compute B matrix
    // UCx = sum(N*URx) -> sum(N*URx) - UCx = 0
    // UCy = sum(N*URy) -> sum(N*URy) - UCy = 0
    static Matrix B(2, 8);
    B.Zero();
    for (int i = 0; i < 2; i++)
        B(i, i) = -1.0;
    for (int i = 0; i < 3; i++) {
        int j = 2 + i * 2;
        B(0, j) = N(i);
        B(1, j + 1) = N(i);
    }

    // Penalty stiffness
    double iK = m_K * std::sqrt(V);

    // compute stiffness
    K.addMatrixTransposeProduct(0.0, B, B, iK);

    // done
    return K;
}

const Matrix& ASDEmbeddedNodeElement::TRI_2D_UR()
{
    // output
    static Matrix K(9, 9);
    
    // collect triangle coordinates
    static Matrix X(2, 3);
    for (int i = 1; i < 4; i++) {
        const Node* node = m_nodes[static_cast<std::size_t>(i)];
        X(0, i - 1) = node->getCrds()(0);
        X(1, i - 1) = node->getCrds()(1);
    }

    // shape functions natural derivatives
    static Matrix dN(3, 2);
    tri::shapeFunDer(dN);

    // jacobian
    static Matrix J(2, 2);
    J.addMatrixProduct(0.0, X, dN, 1.0);
    double detJ = det2(J);
    double V = detJ / 2.0;
    static Matrix invJ(2, 2);
    J.Invert(invJ);

    // shape functions cartesian derivatives
    static Matrix dNdX(3, 2);
    dNdX.addMatrixProduct(0.0, dN, invJ, 1.0);

    // find local coordinates of constrained node
    double lx, ly;
    tri::localCoord(X, invJ, m_nodes[0]->getCrds()(0), m_nodes[0]->getCrds()(1), lx, ly);

    // compute shape functions at constrained node
    static Vector N(3);
    for (int i = 0; i < 3; ++i)
        N(i) = tri::shapeFun(lx, ly, i);

    // compute B matrix
    // UCx = sum(N*URx) -> sum(N*URx) - UCx = 0
    // UCy = sum(N*URy) -> sum(N*URy) - UCy = 0
    // RCz = sum(d_URy_dX - d_URx_dY)/2.0 -> sum(d_URy_dX - d_URx_dY)/2.0 - RCz = 0
    static Matrix B(3, 9);
    B.Zero();
    for (int i = 0; i < 3; i++)
        B(i, i) = -1.0;
    for (int i = 0; i < 3; i++) {
        int j = 3 + i * 2;
        B(0, j) = N(i);
        B(1, j + 1) = N(i);
        B(2, j) = -dNdX(i, 1)/2.0; B(2, j + 1) = dNdX(i, 0)/2.0;
    }

    // Penalty stiffness
    double iK = m_K * std::sqrt(V);

    // compute stiffness
    K.addMatrixTransposeProduct(0.0, B, B, iK);

    // done
    return K;
}

const Matrix& ASDEmbeddedNodeElement::TRI_2D_UP()
{
    // output
    static Matrix K(12, 12);

    // collect triangle coordinates
    static Matrix X(2, 3);
    for (int i = 1; i < 4; i++) {
        const Node* node = m_nodes[static_cast<std::size_t>(i)];
        X(0, i - 1) = node->getCrds()(0);
        X(1, i - 1) = node->getCrds()(1);
    }

    // shape functions natural derivatives
    static Matrix dN(3, 2);
    tri::shapeFunDer(dN);

    // jacobian
    static Matrix J(2, 2);
    J.addMatrixProduct(0.0, X, dN, 1.0);
    double detJ = det2(J);
    double V = detJ / 2.0;
    static Matrix invJ(2, 2);
    J.Invert(invJ);

    // find local coordinates of constrained node
    double lx, ly;
    tri::localCoord(X, invJ, m_nodes[0]->getCrds()(0), m_nodes[0]->getCrds()(1), lx, ly);

    // compute shape functions at constrained node
    static Vector N(3);
    for (int i = 0; i < 3; ++i)
        N(i) = tri::shapeFun(lx, ly, i);

    // compute B matrix
    // UCx = sum(N*URx) -> sum(N*URx) - UCx = 0
    // UCy = sum(N*URy) -> sum(N*URy) - UCy = 0
    // UCp = sum(N*URp) -> sum(N*URp) - UCp = 0
    static Matrix B(3, 12);
    B.Zero();
    for (int i = 0; i < 3; i++)
        B(i, i) = -1.0;
    for (int i = 0; i < 3; i++) {
        int j = 3 + i * 3;
        B(0, j) = N(i);
        B(1, j + 1) = N(i);
        B(2, j + 2) = N(i);
    }

    // Penalty stiffness
    double iK = m_K * std::sqrt(V);
    double iKP = m_KP * std::sqrt(V);
    static Matrix C(3, 3);
    C.Zero();
    C(0, 0) = C(1, 1) = iK;
    C(2, 2) = iKP;

    // compute stiffness
    K.addMatrixTripleProduct(0.0, B, C, 1.0);

    // done
    return K;
}

const Matrix& ASDEmbeddedNodeElement::TRI_3D_U()
{
    // output
    static Matrix K(12, 12);

    // collect triangle coordinates
    static Matrix X(3, 3);
    for (int i = 1; i < 4; i++) {
        const Node* node = m_nodes[static_cast<std::size_t>(i)];
        X(0, i - 1) = node->getCrds()(0);
        X(1, i - 1) = node->getCrds()(1);
        X(2, i - 1) = node->getCrds()(2);
    }

    // shape functions natural derivatives
    static Matrix dN(3, 3);
    dN.Zero(); // note: the tri::shapeFunDer fills the first 2 columns, the 3rd should be zero!
    tri::shapeFunDer(dN);

    // jacobian
    static Matrix J(3, 3);
    J.addMatrixProduct(0.0, X, dN, 1.0);
    tri::fillVzInJacobian(J); // note: the 3rd column will be zero, so fill it with the unit normal vector
    double detJ = det3(J);
    double V = detJ / 2.0;
    static Matrix invJ(3, 3);
    J.Invert(invJ);
    // note: invJ is used as it comes out of the inversion. See the comment on
    // tri::fillVzInJacobian for why the third column must NOT be zeroed here.

    // find local coordinates of constrained node
    double lx, ly;
    tri::localCoord(X, invJ, m_nodes[0]->getCrds()(0), m_nodes[0]->getCrds()(1), m_nodes[0]->getCrds()(2), lx, ly);

    // compute shape functions at constrained node
    static Vector N(3);
    for (int i = 0; i < 3; ++i)
        N(i) = tri::shapeFun(lx, ly, i);

    // compute B matrix
    // UCx = sum(N*URx) -> sum(N*URx) - UCx = 0
    // UCy = sum(N*URy) -> sum(N*URy) - UCy = 0
    // UCz = sum(N*URz) -> sum(N*URz) - UCz = 0
    static Matrix B(3, 12);
    B.Zero();
    for (int i = 0; i < 3; i++)
        B(i, i) = -1.0;
    for (int i = 0; i < 3; i++) {
        int j = 3 + i * 3;
        B(0, j) = N(i);
        B(1, j + 1) = N(i);
        B(2, j + 2) = N(i);
    }

    // Penalty stiffness
    double iK = m_K * std::sqrt(V);

    // compute stiffness
    K.addMatrixTransposeProduct(0.0, B, B, iK);

    // done
    return K;
}

const Matrix& ASDEmbeddedNodeElement::TRI_3D_UR()
{
    // output: 6 constrained dofs + 3 retained nodes with 3 dofs each, or 6
    // each when -shearDeformable exposes their rotations
    static Matrix K;

    // collect triangle coordinates
    // in global coordinates
    static Matrix X(3, 3);
    for (int i = 1; i < 4; i++) {
        const Node* node = m_nodes[static_cast<std::size_t>(i)];
        X(0, i - 1) = node->getCrds()(0);
        X(1, i - 1) = node->getCrds()(1);
        X(2, i - 1) = node->getCrds()(2);
    }

    // compute orientation
    static Vector dx(3);
    static Vector dy(3);
    static Vector dz(3);
    for (int i = 0; i < 3; ++i) {
        dx(i) = X(i, 1) - X(i, 0);
        dy(i) = X(i, 2) - X(i, 0);
    }
    dx.Normalize();
    dy.Normalize();
    cross(dx, dy, dz);
    dz.Normalize();
    cross(dz, dx, dy);

    // assemble orientation matrix (transposed of rotation)
    static Matrix R(3, 3);
    for (int i = 0; i < 3; ++i) {
        R(0, i) = dx(i);
        R(1, i) = dy(i);
        R(2, i) = dz(i);
    }

    // triangle coordinates in local coordinates
    static Matrix XL(2, 3);
    for (int i = 0; i < 3; ++i) {
        XL(0, i) = X(0, i) * dx(0) + X(1, i) * dx(1) + X(2, i) * dx(2);
        XL(1, i) = X(0, i) * dy(0) + X(1, i) * dy(1) + X(2, i) * dy(2);
    }

    // shape functions natural derivatives
    static Matrix dN(3, 2);
    tri::shapeFunDer(dN);

    // jacobian
    static Matrix J(2, 2);
    J.addMatrixProduct(0.0, XL, dN, 1.0);
    double detJ = det2(J);
    double V = detJ / 2.0;
    static Matrix invJ(2, 2);
    J.Invert(invJ);

    // shape functions cartesian derivatives
    static Matrix dNdX(3, 2);
    dNdX.addMatrixProduct(0.0, dN, invJ, 1.0);

    // find local coordinates of constrained node
    const Vector CPos = m_nodes[0]->getCrds();
    double CPosX = CPos(0) * dx(0) + CPos(1) * dx(1) + CPos(2) * dx(2);
    double CPosY = CPos(0) * dy(0) + CPos(1) * dy(1) + CPos(2) * dy(2);
    double lx, ly;
    tri::localCoord(XL, invJ, CPosX, CPosY, lx, ly);

    // compute shape functions at constrained node
    static Vector N(3);
    for (int i = 0; i < 3; ++i)
        N(i) = tri::shapeFun(lx, ly, i);

    // compute B matrix (in local coordinates)
    // UCx = sum(N*URx) -> sum(N*URx) - UCx = 0
    // UCy = sum(N*URy) -> sum(N*URy) - UCy = 0
    // UCz = sum(N*URz) -> sum(N*URz) - UCz = 0
    // default (thin, Kirchhoff): bending rotations from the slope of the
    // transverse displacement
    // RCx = sum( d_URz_dY) -> sum( d_URz_dY) - RCx = 0 (local)
    // RCy = sum(-d_URz_dX) -> sum(-d_URz_dX) - RCy = 0 (local)
    // -shearDeformable (thick, Mindlin): bending rotations from the
    // interpolated nodal rotations of the host, theta = slope + gamma
    // RCx = sum(N*RRx) -> sum(N*RRx) - RCx = 0 (local)
    // RCy = sum(N*RRy) -> sum(N*RRy) - RCy = 0 (local)
    // the drilling rotation is not a director rotation and keeps the skew part
    // of the in-plane gradient in either case:
    // RCz = sum(d_URy_dX - d_URx_dY)/2.0 -> sum(d_URy_dX - d_URx_dY)/2.0 - RCz = 0 (local)
    int nrdof = m_shear ? 6 : 3;   // dofs each retained node exposes
    int ncols = 6 + 3 * nrdof;
    static Matrix B;
    B.resize(6, ncols);
    B.Zero();
    // fill the -identity 6x6 block (transformed to global coordinates)
    // for constrained node+
    /*for (int i = 0; i < 6; ++i)
        B(i, i) = -1.0;*/
    for (int i = 0; i < 2; ++i) {
        int j = i * 3;
        for (int row = 0; row < 3; ++row)
            for (int col = 0; col < 3; ++col)
                B(j + row, j + col) = -R(row, col);
    }
    // fill the 2 rows of 3 blocks (transformed to global coordinates)
    static Matrix BL(3, 3);
    static Matrix BG(3, 3);
    for (int i = 0; i < 3; ++i) {
        int j = 6 + i * nrdof;
        // U block
        BL.Zero();
        BL(0, 0) = N(i);
        BL(1, 1) = N(i);
        BL(2, 2) = N(i);
        BG.addMatrixProduct(0.0, BL, R, 1.0);
        for (int row = 0; row < 3; ++row)
            for (int col = 0; col < 3; ++col)
                B(row, j + col) = BG(row, col);
        // R block on the translational dofs
        BL.Zero();
        if (!m_shear) {
            BL(0, 2) = dNdX(i, 1);
            BL(1, 2) = -dNdX(i, 0);
        }
        BL(2, 0) = -dNdX(i, 1) / 2.0;
        BL(2, 1) = dNdX(i, 0) / 2.0;
        BG.addMatrixProduct(0.0, BL, R, 1.0);
        for (int row = 0; row < 3; ++row)
            for (int col = 0; col < 3; ++col)
                B(3 + row, j + col) = BG(row, col);
        // R block on the rotational dofs: the local bending rotations are the
        // interpolated nodal rotations resolved on the face frame, theta_local
        // = R * theta_global
        if (m_shear) {
            for (int col = 0; col < 3; ++col) {
                B(3, j + 3 + col) = N(i) * R(0, col);
                B(4, j + 3 + col) = N(i) * R(1, col);
            }
        }
    }

    // Penalty stiffness
    double iK = m_K * std::sqrt(V);

    // compute stiffness
    K.resize(ncols, ncols);
    K.Zero();
    K.addMatrixTransposeProduct(0.0, B, B, iK);

    // done
    return K;
}

const Matrix& ASDEmbeddedNodeElement::TRI_3D_UP()
{
    // output
    static Matrix K(16, 16);

    // collect triangle coordinates
    static Matrix X(3, 3);
    for (int i = 1; i < 4; i++) {
        const Node* node = m_nodes[static_cast<std::size_t>(i)];
        X(0, i - 1) = node->getCrds()(0);
        X(1, i - 1) = node->getCrds()(1);
        X(2, i - 1) = node->getCrds()(2);
    }

    // shape functions natural derivatives
    static Matrix dN(3, 3);
    dN.Zero(); // note: the tri::shapeFunDer fills the first 2 columns, the 3rd should be zero!
    tri::shapeFunDer(dN);

    // jacobian
    static Matrix J(3, 3);
    J.addMatrixProduct(0.0, X, dN, 1.0);
    tri::fillVzInJacobian(J); // note: the 3rd column will be zero, so fill it with the unit normal vector
    double detJ = det3(J);
    double V = detJ / 2.0;
    static Matrix invJ(3, 3);
    J.Invert(invJ);
    // note: invJ is used as it comes out of the inversion. See the comment on
    // tri::fillVzInJacobian for why the third column must NOT be zeroed here.

    // find local coordinates of constrained node
    double lx, ly;
    tri::localCoord(X, invJ, m_nodes[0]->getCrds()(0), m_nodes[0]->getCrds()(1), m_nodes[0]->getCrds()(2), lx, ly);

    // compute shape functions at constrained node
    static Vector N(3);
    for (int i = 0; i < 3; ++i)
        N(i) = tri::shapeFun(lx, ly, i);

    // compute B matrix
    // UCx = sum(N*URx) -> sum(N*URx) - UCx = 0
    // UCy = sum(N*URy) -> sum(N*URy) - UCy = 0
    // UCz = sum(N*URz) -> sum(N*URz) - UCz = 0
    // UCp = sum(N*URp) -> sum(N*URp) - UCp = 0
    static Matrix B(4, 16);
    B.Zero();
    for (int i = 0; i < 4; i++)
        B(i, i) = -1.0;
    for (int i = 0; i < 3; i++) {
        int j = 4 + i * 4;
        B(0, j) = N(i);
        B(1, j + 1) = N(i);
        B(2, j + 2) = N(i);
        B(3, j + 3) = N(i);
    }

    // Penalty stiffness
    double iK = m_K * std::sqrt(V);
    double iKP = m_KP * std::sqrt(V);
    static Matrix C(4, 4);
    C.Zero();
    C(0, 0) = C(1, 1) = C(2, 2) = iK;
    C(3, 3) = iKP;

    // compute stiffness
    K.addMatrixTripleProduct(0.0, B, C, 1.0);

    // done
    return K;
}

const Matrix& ASDEmbeddedNodeElement::QUAD_2D(int mode)
{
    // collect the coordinates of the 4 retained nodes (one column per node)
    static Matrix X(2, 4);
    for (int i = 1; i < 5; i++) {
        const Node* node = m_nodes[static_cast<std::size_t>(i)];
        X(0, i - 1) = node->getCrds()(0);
        X(1, i - 1) = node->getCrds()(1);
    }

    // natural coordinate of the constrained node: Newton, the map is not affine
    double lx, ly;
    if (!quad::localCoord(X, m_nodes[0]->getCrds()(0), m_nodes[0]->getCrds()(1), lx, ly)) {
        opserr << "ASDEmbeddedNodeElement WARNING: element " << getTag()
            << " - the inverse isoparametric map did not converge on the QUAD host. "
            << "Check that the host is convex and that the constrained node is inside it.\n";
    }

    // shape functions and their cartesian derivatives at that point
    static Vector N(4);
    quad::shapeFun(lx, ly, N);
    static Matrix dN(4, 2);
    static Matrix J(2, 2);
    static Matrix invJ(2, 2);
    static Matrix dNdX(4, 2);
    quad::shapeFunDer(lx, ly, dN);
    J.addMatrixProduct(0.0, X, dN, 1.0);
    J.Invert(invJ);
    dNdX.addMatrixProduct(0.0, dN, invJ, 1.0);

    // penalty stiffness, scaled by the element size as in the simplex kernels
    double A = quad::area(X);
    double iK = m_K * std::sqrt(A);
    double iKP = m_KP * std::sqrt(A);

    return assembleConstraint(2, 4, N, dNdX, mode, iK, iKP);
}

const Matrix& ASDEmbeddedNodeElement::QUAD_3D(int mode)
{
    // A quadrilateral host in 3D is a SURFACE: two natural coordinates, three
    // equations. It is handled the way TRI_3D_UR handles the triangle -- build
    // an orthonormal frame on the face, resolve everything on it, and work in a
    // genuine 2D problem -- rather than the way TRI_3D_U used to, which is only
    // valid for faces normal to the global z axis.
    //
    // Once the face coordinates are on the plane the map is the ordinary
    // bilinear one, so the natural coordinate comes from the same Newton
    // iteration QUAD_2D uses. That is the whole of the difference: frame, then
    // the 2D quadrilateral.

    // collect the coordinates of the 4 retained nodes (one column per node)
    static Matrix X(3, 4);
    for (int i = 1; i < 5; i++) {
        const Node* node = m_nodes[static_cast<std::size_t>(i)];
        X(0, i - 1) = node->getCrds()(0);
        X(1, i - 1) = node->getCrds()(1);
        X(2, i - 1) = node->getCrds()(2);
    }

    // local frame of the face and the face coordinates resolved on it
    static Matrix R(3, 3);
    static Vector center(3);
    static Matrix XL(2, 4);
    double warp = 0.0;
    quad::frame3D(X, R, center, XL, warp);

    // constrained node, projected onto the plane of the face
    const Vector& CPos = m_nodes[0]->getCrds();
    double cdx = CPos(0) - center(0);
    double cdy = CPos(1) - center(1);
    double cdz = CPos(2) - center(2);
    double CPosX = R(0, 0) * cdx + R(0, 1) * cdy + R(0, 2) * cdz;
    double CPosY = R(1, 0) * cdx + R(1, 1) * cdy + R(1, 2) * cdz;

    // natural coordinate of the constrained node: Newton, the map is not affine
    double lx, ly;
    if (!quad::localCoord(XL, CPosX, CPosY, lx, ly)) {
        opserr << "ASDEmbeddedNodeElement WARNING: element " << getTag()
            << " - the inverse isoparametric map did not converge on the QUAD host in 3D. "
            << "Check that the face is convex and that the constrained node projects inside it.\n";
    }

    // shape functions at that point
    static Vector N(4);
    quad::shapeFun(lx, ly, N);

    // penalty stiffness, scaled by the element size as in every other kernel
    double A = quad::area(XL);
    double iK = m_K * std::sqrt(A);
    double iKP = m_KP * std::sqrt(A);

    // The translation and pressure constraints only involve N, and N is a
    // scalar weight: they read the same in any frame, so the generic assembler
    // applies unchanged. Only the rotational mode needs the frame, because a
    // rotation is defined from the in-plane gradient and the normal direction.
    if (mode != Mode_UR) {
        static Matrix dNdX_unused(4, 3);
        dNdX_unused.Zero();
        return assembleConstraint(3, 4, N, dNdX_unused, mode, iK, iKP);
    }

    // ---- rotational mode: same construction as TRI_3D_UR, with 4 nodes ----

    // cartesian derivatives on the face
    static Matrix dN(4, 2);
    static Matrix J(2, 2);
    static Matrix invJ(2, 2);
    static Matrix dNdX(4, 2);
    quad::shapeFunDer(lx, ly, dN);
    J.addMatrixProduct(0.0, XL, dN, 1.0);
    J.Invert(invJ);
    dNdX.addMatrixProduct(0.0, dN, invJ, 1.0);

    // B in local components, then rotated to global.
    // rows: 3 translations + 3 rotations of the constrained node
    // cols: 6 constrained dofs + 4 retained nodes x 3 dofs (x 6 dofs when
    // -shearDeformable exposes their rotations)
    int nrdof = m_shear ? 6 : 3;   // dofs each retained node exposes
    int ncols = 6 + 4 * nrdof;
    static Matrix B;
    B.resize(6, ncols);
    B.Zero();
    for (int i = 0; i < 2; ++i) {
        int j = i * 3;
        for (int row = 0; row < 3; ++row)
            for (int col = 0; col < 3; ++col)
                B(j + row, j + col) = -R(row, col);
    }
    static Matrix BL(3, 3);
    static Matrix BG(3, 3);
    for (int i = 0; i < 4; ++i) {
        int j = 6 + i * nrdof;
        // U block
        BL.Zero();
        BL(0, 0) = N(i);
        BL(1, 1) = N(i);
        BL(2, 2) = N(i);
        BG.addMatrixProduct(0.0, BL, R, 1.0);
        for (int row = 0; row < 3; ++row)
            for (int col = 0; col < 3; ++col)
                B(row, j + col) = BG(row, col);
        // R block on the translational dofs: by default the two bending
        // rotations come from the gradient of the out-of-plane displacement
        // (thin, Kirchhoff); with -shearDeformable they come from the nodal
        // rotations below instead. The drilling rotation is not a director
        // rotation and keeps the skew part of the in-plane gradient in
        // either case.
        BL.Zero();
        if (!m_shear) {
            BL(0, 2) = dNdX(i, 1);
            BL(1, 2) = -dNdX(i, 0);
        }
        BL(2, 0) = -dNdX(i, 1) / 2.0;
        BL(2, 1) = dNdX(i, 0) / 2.0;
        BG.addMatrixProduct(0.0, BL, R, 1.0);
        for (int row = 0; row < 3; ++row)
            for (int col = 0; col < 3; ++col)
                B(3 + row, j + col) = BG(row, col);
        // R block on the rotational dofs: the local bending rotations are the
        // interpolated nodal rotations resolved on the face frame, theta_local
        // = R * theta_global
        if (m_shear) {
            for (int col = 0; col < 3; ++col) {
                B(3, j + 3 + col) = N(i) * R(0, col);
                B(4, j + 3 + col) = N(i) * R(1, col);
            }
        }
    }

    static Matrix K;
    K.resize(ncols, ncols);
    K.Zero();
    K.addMatrixTransposeProduct(0.0, B, B, iK);
    return K;
}

const Matrix& ASDEmbeddedNodeElement::HEX_3D(int mode)
{
    // collect the coordinates of the 8 retained nodes (one column per node)
    static Matrix X(3, 8);
    for (int i = 1; i < 9; i++) {
        const Node* node = m_nodes[static_cast<std::size_t>(i)];
        X(0, i - 1) = node->getCrds()(0);
        X(1, i - 1) = node->getCrds()(1);
        X(2, i - 1) = node->getCrds()(2);
    }

    // natural coordinate of the constrained node
    double lx, ly, lz;
    if (!hexa::localCoord(X, m_nodes[0]->getCrds()(0), m_nodes[0]->getCrds()(1),
        m_nodes[0]->getCrds()(2), lx, ly, lz)) {
        opserr << "ASDEmbeddedNodeElement WARNING: element " << getTag()
            << " - the inverse isoparametric map did not converge on the HEXA host. "
            << "Check that the host is convex and that the constrained node is inside it.\n";
    }

    // shape functions and their cartesian derivatives at that point
    static Vector N(8);
    hexa::shapeFun(lx, ly, lz, N);
    static Matrix dN(8, 3);
    static Matrix J(3, 3);
    static Matrix invJ(3, 3);
    static Matrix dNdX(8, 3);
    hexa::shapeFunDer(lx, ly, lz, dN);
    J.addMatrixProduct(0.0, X, dN, 1.0);
    J.Invert(invJ);
    dNdX.addMatrixProduct(0.0, dN, invJ, 1.0);

    // penalty stiffness
    double V = hexa::volume(X);
    double iK = m_K * std::cbrt(V);
    double iKP = m_KP * std::cbrt(V);

    return assembleConstraint(3, 8, N, dNdX, mode, iK, iKP);
}

const Matrix& ASDEmbeddedNodeElement::TET_3D_U()
{
    // output
    static Matrix K(15, 15);

    // collect tetrahedron coordinates
    static Matrix X(3, 4);
    for (int i = 1; i < 5; i++) {
        const Node* node = m_nodes[static_cast<std::size_t>(i)];
        X(0, i - 1) = node->getCrds()(0);
        X(1, i - 1) = node->getCrds()(1);
        X(2, i - 1) = node->getCrds()(2);
    }

    // shape functions natural derivatives
    static Matrix dN(4, 3);
    tet::shapeFunDer(dN);

    // jacobian
    static Matrix J(3, 3);
    J.addMatrixProduct(0.0, X, dN, 1.0);
    double detJ = det3(J);
    double V = detJ / 6.0;
    static Matrix invJ(3, 3);
    J.Invert(invJ);

    // find local coordinates of constrained node
    double lx, ly, lz;
    tet::localCoord(X, invJ, m_nodes[0]->getCrds()(0), m_nodes[0]->getCrds()(1), m_nodes[0]->getCrds()(2), lx, ly, lz);

    // compute shape functions at constrained node
    static Vector N(4);
    for (int i = 0; i < 4; ++i)
        N(i) = tet::shapeFun(lx, ly, lz, i);

    // compute B matrix
    // UCx = sum(N*URx) -> sum(N*URx) - UCx = 0
    // UCy = sum(N*URy) -> sum(N*URy) - UCy = 0
    // UCz = sum(N*URz) -> sum(N*URz) - UCz = 0
    static Matrix B(3, 15);
    B.Zero();
    for (int i = 0; i < 3; i++)
        B(i, i) = -1.0;
    for (int i = 0; i < 4; i++) {
        int j = 3 + i * 3;
        B(0, j) = N(i);
        B(1, j + 1) = N(i);
        B(2, j + 2) = N(i);
    }

    // Penalty stiffness
    double iK = m_K * std::cbrt(V);

    // compute stiffness
    K.addMatrixTransposeProduct(0.0, B, B, iK);

    // done
    return K;
}

const Matrix& ASDEmbeddedNodeElement::TET_3D_UR()
{
    // output
    static Matrix K(18, 18);

    // collect tetrahedron coordinates
    static Matrix X(3, 4);
    for (int i = 1; i < 5; i++) {
        const Node* node = m_nodes[static_cast<std::size_t>(i)];
        X(0, i - 1) = node->getCrds()(0);
        X(1, i - 1) = node->getCrds()(1);
        X(2, i - 1) = node->getCrds()(2);
    }

    // shape functions natural derivatives
    static Matrix dN(4, 3);
    tet::shapeFunDer(dN);

    // jacobian
    static Matrix J(3, 3);
    J.addMatrixProduct(0.0, X, dN, 1.0);
    double detJ = det3(J);
    double V = detJ / 6.0;
    static Matrix invJ(3, 3);
    J.Invert(invJ);

    // shape functions cartesian derivatives
    static Matrix dNdX(4, 3);
    dNdX.addMatrixProduct(0.0, dN, invJ, 1.0);

    // find local coordinates of constrained node
    double lx, ly, lz;
    tet::localCoord(X, invJ, m_nodes[0]->getCrds()(0), m_nodes[0]->getCrds()(1), m_nodes[0]->getCrds()(2), lx, ly, lz);

    // compute shape functions at constrained node
    static Vector N(4);
    for (int i = 0; i < 4; ++i)
        N(i) = tet::shapeFun(lx, ly, lz, i);

    // compute B matrix
    // UCx = sum(N*URx) -> sum(N*URx) - UCx = 0
    // UCy = sum(N*URy) -> sum(N*URy) - UCy = 0
    // UCz = sum(N*URz) -> sum(N*URz) - UCz = 0
    // RCx = sum(d_URz_dY - d_URy_dZ)/2.0 -> sum(d_URz_dY - d_URy_dZ)/2.0 - RCx = 0
    // RCy = sum(d_URx_dZ - d_URz_dX)/2.0 -> sum(d_URx_dZ - d_URz_dX)/2.0 - RCy = 0
    // RCz = sum(d_URy_dX - d_URx_dY)/2.0 -> sum(d_URy_dX - d_URx_dY)/2.0 - RCz = 0
    static Matrix B(6, 18);
    B.Zero();
    for (int i = 0; i < 6; i++)
        B(i, i) = -1.0;
    for (int i = 0; i < 4; i++) {
        int j = 6 + i * 3;
        B(0, j) = N(i);
        B(1, j + 1) = N(i);
        B(2, j + 2) = N(i);
        B(3, j + 1) = -dNdX(i, 2) / 2.0; B(3, j + 2) = dNdX(i, 1) / 2.0;
        B(4, j) = dNdX(i, 2) / 2.0; B(4, j + 2) = -dNdX(i, 0) / 2.0;
        B(5, j) = -dNdX(i, 1) / 2.0; B(5, j + 1) = dNdX(i, 0) / 2.0;
    }

    // Penalty stiffness
    double iK = m_K * std::cbrt(V);

    // compute stiffness
    K.addMatrixTransposeProduct(0.0, B, B, iK);

    // done
    return K;
}

const Matrix& ASDEmbeddedNodeElement::TET_3D_UP()
{
    // output
    static Matrix K(20, 20);

    // collect tetrahedron coordinates
    static Matrix X(3, 4);
    for (int i = 1; i < 5; i++) {
        const Node* node = m_nodes[static_cast<std::size_t>(i)];
        X(0, i - 1) = node->getCrds()(0);
        X(1, i - 1) = node->getCrds()(1);
        X(2, i - 1) = node->getCrds()(2);
    }

    // shape functions natural derivatives
    static Matrix dN(4, 3);
    tet::shapeFunDer(dN);

    // jacobian
    static Matrix J(3, 3);
    J.addMatrixProduct(0.0, X, dN, 1.0);
    double detJ = det3(J);
    double V = detJ / 6.0;
    static Matrix invJ(3, 3);
    J.Invert(invJ);

    // find local coordinates of constrained node
    double lx, ly, lz;
    tet::localCoord(X, invJ, m_nodes[0]->getCrds()(0), m_nodes[0]->getCrds()(1), m_nodes[0]->getCrds()(2), lx, ly, lz);

    // compute shape functions at constrained node
    static Vector N(4);
    for (int i = 0; i < 4; ++i)
        N(i) = tet::shapeFun(lx, ly, lz, i);

    // compute B matrix
    // UCx = sum(N*URx) -> sum(N*URx) - UCx = 0
    // UCy = sum(N*URy) -> sum(N*URy) - UCy = 0
    // UCz = sum(N*URz) -> sum(N*URz) - UCz = 0
    // UCp = sum(N*URp) -> sum(N*URp) - UCp = 0
    static Matrix B(4, 20);
    B.Zero();
    for (int i = 0; i < 4; i++)
        B(i, i) = -1.0;
    for (int i = 0; i < 4; i++) {
        int j = 4 + i * 4;
        B(0, j) = N(i);
        B(1, j + 1) = N(i);
        B(2, j + 2) = N(i);
        B(3, j + 3) = N(i);
    }

    // Penalty stiffness
    double iK = m_K * std::cbrt(V);
    double iKP = m_KP * std::cbrt(V);
    static Matrix C(4, 4);
    C.Zero();
    C(0, 0) = C(1, 1) = C(2, 2) = iK;
    C(3, 3) = iKP;

    // compute stiffness
    K.addMatrixTripleProduct(0.0, B, C, 1.0);

    // done
    return K;
}
