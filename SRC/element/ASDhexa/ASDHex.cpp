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

// 

// Original implementation:
//
// An 8-node solid hexahedral element for three-dimensional continuum analysis,
// based on an enhanced assumed strain / Petrov-Galerkin formulation (PG-EAS).
//
// It supports both linear and corotational kinematics.
//

#include <ASDHex.h>
#include <ASDSolidHexCorotationalTransformation.h>


#include <Domain.h>
#include <ErrorHandler.h>
#include <ElementResponse.h>
#include <ElementalLoad.h>
#include <Channel.h>
#include <FEM_ObjectBroker.h>
#include <elementAPI.h>
#include <Renderer.h>
#include <Damping.h>

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <cmath>
#include <cstdint>

// ASDHEX_EAS_IMPERF is declared in ASDHex.h -- it gates class members, so it has
// to be visible before the header is parsed. See the block there for what it is
// for and why it is not shipped enabled.

#include <ASDEICR3.h>

void*
OPS_ASDSolidHex(void)
{
    static bool first_done = false;
    if (!first_done) {
        opserr << "Using ASDSolidHex - Developed by:\n";
        first_done = true;
    }

    // args: tag, 8 nodes, matTag
    int numArgs = OPS_GetNumRemainingInputArgs();
    if (numArgs < 10) {
        // NOTE: the keyword registered in Tcl/Python is "ASDHex", and the only
        // option actually accepted by the loop below is -corotational.
        opserr << "Want: element ASDHex $tag $Node1 $Node2 $Node3 $Node4 $Node5 $Node6 $Node7 $Node8 $matTag "
            "<-corotational> <-damp $dampingTag> <-b $bx $by $bz> <-easStab $s> <-easPenalty $p> <-noEasAuto>\n";
        return 0;
    }

    int iData[10];
    int numData = 10;
    if (OPS_GetInt(&numData, iData) != 0) {
        opserr << "WARNING invalid integer tag: element ASDHex \n";
        return 0;
    }


    int matTag = iData[9];

    // The solid material is treated as an n-dimensional material.
    NDMaterial* mat = OPS_getNDMaterial(matTag);
    if (mat == 0) {
        opserr << "ERROR: element ASDHex " << iData[0]
            << " NDMaterial " << matTag << " not found\n";
        return 0;
    }

    bool m_use_corotational = false;
    Damping* damping = nullptr;
    double eas_stab = 0.0;
    double eas_penalty = 0.0;
    bool eas_auto = true;
#ifdef ASDHEX_EAS_IMPERF
    double imperfection = 0.0;
#endif
    double body[3] = { 0.0, 0.0, 0.0 };
    while (OPS_GetNumRemainingInputArgs() > 0) {
        const char* type = OPS_GetString();
        if (strcmp(type, "-corotational") == 0) {
            m_use_corotational = true;
        }
        else if (strcmp(type, "-damp") == 0) {
            if (OPS_GetNumRemainingInputArgs() < 1) {
                opserr << "Error: element ASDHex: -damp needs a damping tag\n";
                return 0;
            }
            int dampingTag = 0;
            int nd = 1;
            if (OPS_GetIntInput(&nd, &dampingTag) < 0) {
                opserr << "Error: element ASDHex: invalid damping tag\n";
                return 0;
            }
            damping = OPS_getDamping(dampingTag);
            if (damping == nullptr) {
                opserr << "Error: element ASDHex: damping " << dampingTag << " not found\n";
                return 0;
            }
        }
        else if (strcmp(type, "-b") == 0) {
            // body force per unit mass (same meaning as Brick's b1 b2 b3), used by
            // eleLoad -type -brickSelfWeight / -selfWeight
            if (OPS_GetNumRemainingInputArgs() < 3) {
                opserr << "Error: element ASDHex: -b needs 3 components\n";
                return 0;
            }
            int nd = 3;
            if (OPS_GetDoubleInput(&nd, body) < 0) {
                opserr << "Error: element ASDHex: invalid -b components\n";
                return 0;
            }
        }
#ifdef ASDHEX_EAS_IMPERF
        else if (strcmp(type, "-imperfection") == 0) {
            // EXPERIMENT: see the ASDHEX_EAS_IMPERF block at the top of this file
            if (OPS_GetNumRemainingInputArgs() < 1) {
                opserr << "Error: element ASDHex: -imperfection needs 1 value\n";
                return 0;
            }
            int nd = 1;
            if (OPS_GetDoubleInput(&nd, &imperfection) < 0) {
                opserr << "Error: element ASDHex: invalid -imperfection value\n";
                return 0;
            }
            if (imperfection < 0.0 || imperfection >= 1.0) {
                opserr << "Error: element ASDHex: -imperfection must be in [0, 1), got "
                    << imperfection << "\n";
                return 0;
            }
        }
#endif // ASDHEX_EAS_IMPERF
        else if (strcmp(type, "-noEasAuto") == 0) {
            // disable the AUTOMATIC A-op arming, leaving only whatever -easStab
            // and -easPenalty ask for explicitly. This is how the two
            // stabilizations are told apart: with the automatic path on, A-op has
            // already rescued a softening element before the penalty can act.
            eas_auto = false;
        }
        else if (strcmp(type, "-easPenalty") == 0) {
            // TRUE penalty on the enhanced modes. Unlike -easStab this one enters
            // the RESIDUAL, so it CHANGES THE CONVERGED ANSWER: it trades bending
            // accuracy for an energetic barrier against the enhanced modes
            // localizing. See ASDSolidHex::m_eas_penalty.
            //
            // The patch test is exact at any value; bending is what pays.
            if (OPS_GetNumRemainingInputArgs() < 1) {
                opserr << "Error: element ASDHex: -easPenalty needs 1 value\n";
                return 0;
            }
            int nd = 1;
            if (OPS_GetDoubleInput(&nd, &eas_penalty) < 0) {
                opserr << "Error: element ASDHex: invalid -easPenalty value\n";
                return 0;
            }
            if (eas_penalty < 0.0) {
                opserr << "Error: element ASDHex: -easPenalty must be non-negative, got "
                    << eas_penalty << "\n";
                return 0;
            }
        }
        else if (strcmp(type, "-easStab") == 0) {
            // floor on the regularization of the enhanced (EAS) block. See
            // ASDSolidHex::m_eas_stab: the gate in calculateAll raises the factor
            // on its own whenever the block is not positive definite, so this is
            // only needed to force a minimum amount of it.
            //
            // Because ONLY the operator is regularized -- the enhanced residual
            // stays exact, so at convergence h = 0 and the term it feeds into the
            // condensed residual vanishes -- the converged response does NOT depend
            // on this value. It trades global iterations for robustness, never
            // accuracy.
            if (OPS_GetNumRemainingInputArgs() < 1) {
                opserr << "Error: element ASDHex: -easStab needs 1 value\n";
                return 0;
            }
            int nd = 1;
            if (OPS_GetDoubleInput(&nd, &eas_stab) < 0) {
                opserr << "Error: element ASDHex: invalid -easStab value\n";
                return 0;
            }
            if (eas_stab < 0.0) {
                opserr << "Error: element ASDHex: -easStab must be non-negative, got "
                    << eas_stab << "\n";
                return 0;
            }
        }
        else {
            opserr << "Error: element ASDHex: unknown option '" << type << "'\n";
            return 0;
        }
    }

    return new ASDSolidHex(iData[0], //tag
        iData[1], iData[2], iData[3], iData[4], iData[5], iData[6], iData[7], iData[8], //8 nodes
        mat, m_use_corotational, damping, body, eas_stab, eas_penalty, eas_auto
#ifdef ASDHEX_EAS_IMPERF
        , imperfection
#endif
    );
}

// anonymous namespace for utilities
namespace
{
    // some typedefs
    typedef ASDVector3<double> Vector3Type;
    using vec3 = Vector3Type;

    constexpr int NumGP = 8;

    // calculation options
    constexpr int OPT_NONE = 0;
    constexpr int OPT_UPDATE = (1 << 0);
    constexpr int OPT_LHS = (1 << 1);
    constexpr int OPT_RHS = (1 << 2);
    constexpr int OPT_LHS_IS_INITIAL = (1 << 3);

    // ------------------------------------------------------------------
    // A-op: regularization of the enhanced (EAS) block.
    //
    // The amount of regularization is decided POINTWISE ON THE MATERIAL TANGENT,
    // never on the assembled enhanced block.
    //
    // That distinction was learned the hard way. sym(k_qq) is NOT positive
    // definite even for a purely elastic element -- k_qq = int(G_test^T C G_trial)
    // pairs two DIFFERENT interpolations, so its symmetric part carries no
    // definiteness property at all, while the operator itself stays perfectly
    // well conditioned (measured: min pivot 0 at rcond 0.5 on P2's elastic
    // bending patch). Gating on it inflated s to O(500) on healthy elements,
    // degenerated them to displacement elements and brought the parasitic shear
    // back. The material tangent, by contrast, is symmetric and its loss of
    // definiteness IS the softening the regularization exists for.
    //
    // Requiring C + s*C0 to be positive definite at every gauss point is the
    // strongest guarantee available cheaply, and it implies the assembled block is
    // the Petrov block of a definite material tangent.
    //
    // The margin is relative to the ELASTIC tangent's scale, which does not shrink
    // as the material damages, and is deliberately small: a healthy point sits at
    // O(0.1-1), a softening one goes straight through zero, so nothing in between
    // needs to be caught.
    constexpr double EAS_MAT_PD_TOL = 1.0e-6;

    // ...but CAPPED, and the cap is the honest part of this scheme.
    //
    // Measured on a single element in uniaxial tension with a brittle softening
    // ASDConcrete3D and the material's CONSISTENT tangent (-tangent): the worst
    // gauss point sits at a min pivot of -431 times the elastic diagonal, so
    // making C + s*C0 definite there needs s of order 650. At that value the
    // enhanced modes are suppressed by 99.8% -- the element has become a plain
    // displacement element, which is not a stabilized EAS element, it is a
    // different (and much worse) element. And it still did not converge.
    //
    // So s is bounded. Past the bound the operator may stay indefinite; it only
    // has to stay INVERTIBLE, which is all the condensation and the lagged alpha
    // step actually require, and the trust region in updatePG_EAS is what keeps
    // that step from running away.
    //
    // The real fix at that point is not more s. It is to stop handing the element
    // an indefinite tangent: ASDConcrete3D's DEFAULT tangent is the secant W:C0,
    // positive semi-definite whatever the backbone does, and -implex makes the
    // material exactly affine over the step with C = W:C0 as its EXACT derivative.
    // Both were measured to run the same softening branch to completion with the
    // 8 gauss points on one path; only -tangent breaks.
    constexpr double EAS_STAB_MAX = 1.0;


    // 2x2x2 gauss quadrature data, 8 gauss point (2 values for xi,2 for eta, 2 for zeta)
    // INTERNAL ordering is "corner" order, i.e. the same sequence as the element's
    // node numbering:
    //   0:(-,-,-) 1:(+,-,-) 2:(+,+,-) 3:(-,+,-) 4:(-,-,+) 5:(+,-,+) 6:(+,+,+) 7:(-,+,+)
    // See GP_REPORT_TO_INTERNAL below for why this is NOT the order used when
    // reporting results to recorders.
    constexpr double GP = 0.577350269189626;
    constexpr std::array<double, NumGP> XI = { -GP, GP, GP, -GP, -GP, GP, GP, -GP };
    constexpr std::array<double, NumGP> ETA = { -GP, -GP, GP, GP, -GP, -GP, GP, GP };
    constexpr std::array<double, NumGP> ZETA = { -GP, -GP, -GP, -GP, GP, GP, GP, GP };

    // Gauss point index mapping: REPORTING order -> INTERNAL quadrature order.
    //
    // Results are reported in the LEXICOGRAPHIC (xi,eta,zeta) order, xi outermost,
    // i.e. the order produced by the i-j-k loop of Brick, BbarBrick and BrickUP:
    //   0:(-,-,-) 1:(-,-,+) 2:(-,+,-) 3:(-,+,+) 4:(+,-,-) 5:(+,-,+) 6:(+,+,-) 7:(+,+,+)
    // This is deliberate: it makes the gauss point output of this element
    // interchangeable with that of the other 8-node bricks, so anything already
    // able to read Brick reads this element too.
    //
    // The INTERNAL order is deliberately left alone: metric_basis::orthogonalize()
    // performs a Gram-Schmidt process over the gauss points (condition C3 of the
    // Petrov-Galerkin formulation), and Gram-Schmidt is order dependent. Permuting
    // the quadrature was measured to perturb every result at round-off level
    // (~1e-14 relative; the beam-distortion sweep moved from -1.999999999999997 to
    // -2.000000000000033), which is harmless but needlessly invalidates the
    // validated reference data. Applying the permutation only where an index is
    // exposed keeps the formulation bit-identical.
    //
    // GP_REPORT_TO_INTERNAL[j] = internal index of the j-th reported gauss point.
    constexpr std::array<int, NumGP> GP_REPORT_TO_INTERNAL = { 0, 4, 3, 7, 1, 5, 2, 6 };

    // true only if every entry of m is a finite number. Used to reject an
    // inverse that LAPACK reported as successful but filled with inf/nan.
    inline bool isFinite(const Matrix& m)
    {
        for (int i = 0; i < m.noRows(); ++i)
            for (int j = 0; j < m.noCols(); ++j)
                if (!std::isfinite(m(i, j)))
                    return false;
        return true;
    }

    inline bool isFinite(const Vector& v)
    {
        for (int i = 0; i < v.Size(); ++i)
            if (!std::isfinite(v(i)))
                return false;
        return true;
    }

    // infinity norm (largest absolute row sum). Matrix offers no norms at all --
    // only Solve, Invert and a 3x3 Eigen3 -- so the two measures the enhanced
    // solve is judged by are hand rolled here.
    inline double normInf(const Matrix& m)
    {
        double n = 0.0;
        for (int i = 0; i < m.noRows(); ++i) {
            double rowsum = 0.0;
            for (int j = 0; j < m.noCols(); ++j)
                rowsum += std::abs(m(i, j));
            if (rowsum > n)
                n = rowsum;
        }
        return n;
    }

    // smallest pivot of the LDL^T factorization of the SYMMETRIC PART of m.
    // The SIGN of the result is a definiteness test: a non-positive minimum pivot
    // means sym(m) is not positive definite.
    //
    // m itself is NOT required to be symmetric, and the Petrov-Galerkin blocks of
    // this element never are: B_test != B_trial and G_test != G_trial. What decides
    // whether the enhanced solve is well posed is the definiteness of the symmetric
    // part, so that is what is factorized -- sym(m)(i,j) = (m(i,j) + m(j,i))/2,
    // whose diagonal is just m(i,i).
    //
    // No square roots and no pivoting, ~n^3/6 flops: this is called on every
    // evaluation, both as a diagnostic and as the gate of the adaptive enhanced
    // block regularization, so it has to be cheap. Without pivoting a matrix that
    // IS positive definite always survives (Cholesky-like stability), and one that
    // is not gets caught by the non-positive pivot, which is the only verdict asked
    // of it.
    inline double ldltMinPivot(const Matrix& m)
    {
        constexpr int MaxN = 12;
        const int n = m.noRows();
        if (n < 1 || n > MaxN || m.noCols() != n)
            return 0.0;
        double L[MaxN][MaxN] = { {0.0} };  // unit lower triangular
        double d[MaxN] = { 0.0 };          // the pivots
        double dmin = 0.0;
        for (int j = 0; j < n; ++j) {
            double dj = m(j, j);
            for (int k = 0; k < j; ++k)
                dj -= L[j][k] * L[j][k] * d[k];
            d[j] = dj;
            if (j == 0 || dj < dmin)
                dmin = dj;
            // a non-positive pivot has already answered the question, and the
            // elimination below would divide by it: stop and report it
            if (!(dj > 0.0))
                return dmin;
            for (int i = j + 1; i < n; ++i) {
                double lij = 0.5 * (m(i, j) + m(j, i));
                for (int k = 0; k < j; ++k)
                    lij -= L[i][k] * L[j][k] * d[k];
                L[i][j] = lij / dj;
            }
        }
        return dmin;
    }

#ifdef ASDHEX_EAS_IMPERF
    // splitmix64 on (element tag, internal gauss point index), mapped to [-1, 1).
    //
    // Deterministic on every platform, every partition count and every rerun:
    // uint64_t multiplication wraps in a way the standard defines. Deliberately
    // NOT std::hash (implementation defined) and NOT size_t arithmetic (32 vs 64
    // bit divergence), and obviously not rand() or anything seeded from a clock --
    // an MP run must produce bit-identical results to a serial one.
    //
    // The gauss point index is the INTERNAL one, never the reported one: the
    // Gram-Schmidt in orthogonalize() is order dependent, so the internal order is
    // the only one the compute path may use.
    inline double imperfectionUnit(int tag, int gp)
    {
        std::uint64_t z =
            (0x9E3779B97F4A7C15ULL * static_cast<std::uint64_t>(static_cast<std::int64_t>(tag)))
            ^ (0xBF58476D1CE4E5B9ULL * static_cast<std::uint64_t>(gp + 1));
        z += 0x9E3779B97F4A7C15ULL;
        z = (z ^ (z >> 30)) * 0xBF58476D1CE4E5B9ULL;
        z = (z ^ (z >> 27)) * 0x94D049BB133111EBULL;
        z = z ^ (z >> 31);
        // top 53 bits -> [0, 1) exactly representable, then to [-1, 1)
        return 2.0 * (static_cast<double>(z >> 11) / 9007199254740992.0) - 1.0;
    }
#endif // ASDHEX_EAS_IMPERF

    // mean diagonal of m, the intrinsic scale the minimum pivot is normalized by.
    // trace/n rather than ||m|| because a pivot is a diagonal quantity, so the
    // ratio is dimensionless and mesh-size independent.
    inline double meanDiagonal(const Matrix& m)
    {
        const int n = m.noRows() < m.noCols() ? m.noRows() : m.noCols();
        if (n < 1)
            return 0.0;
        double t = 0.0;
        for (int i = 0; i < n; ++i)
            t += m(i, i);
        return t / static_cast<double>(n);
    }

    // weights
    constexpr std::array<double, NumGP> WTS = { 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0 };

    // shape functions isoparametric : 1/8 (1+xi*xi_a)(1+eta*eta_a)(1+zeta*zeta_a)
    inline void shapeFunctions(double xi, double eta, double zeta, Vector& N)
    {
        N(0) = 0.125 * (1.0 - xi) * (1.0 - eta) * (1 - zeta);
        N(1) = 0.125 * (1.0 + xi) * (1.0 - eta) * (1 - zeta);
        N(2) = 0.125 * (1.0 + xi) * (1.0 + eta) * (1 - zeta);
        N(3) = 0.125 * (1.0 - xi) * (1.0 + eta) * (1 - zeta);
        N(4) = 0.125 * (1.0 - xi) * (1.0 - eta) * (1 + zeta);
        N(5) = 0.125 * (1.0 + xi) * (1.0 - eta) * (1 + zeta);
        N(6) = 0.125 * (1.0 + xi) * (1.0 + eta) * (1 + zeta);
        N(7) = 0.125 * (1.0 - xi) * (1.0 + eta) * (1 + zeta);
    }

    // shape function derivatives dN/d* : [8 x 3] columns (d/dxi, d/deta, d/dzeta)
    inline void dshape(double xi, double eta, double zeta, Matrix& dNdh)
    {
        dNdh.Zero();

        // d/dxi
        dNdh(0, 0) = -0.125 * (1.0 - eta) * (1 - zeta);
        dNdh(1, 0) = 0.125 * (1.0 - eta) * (1 - zeta);
        dNdh(2, 0) = 0.125 * (1.0 + eta) * (1 - zeta);
        dNdh(3, 0) = -0.125 * (1.0 + eta) * (1 - zeta);
        dNdh(4, 0) = -0.125 * (1.0 - eta) * (1 + zeta);
        dNdh(5, 0) = 0.125 * (1.0 - eta) * (1 + zeta);
        dNdh(6, 0) = 0.125 * (1.0 + eta) * (1 + zeta);
        dNdh(7, 0) = -0.125 * (1.0 + eta) * (1 + zeta);


        // d/deta
        dNdh(0, 1) = -0.125 * (1.0 - xi) * (1 - zeta);
        dNdh(1, 1) = -0.125 * (1.0 + xi) * (1 - zeta);
        dNdh(2, 1) = 0.125 * (1.0 + xi) * (1 - zeta);
        dNdh(3, 1) = 0.125 * (1.0 - xi) * (1 - zeta);
        dNdh(4, 1) = -0.125 * (1.0 - xi) * (1 + zeta);
        dNdh(5, 1) = -0.125 * (1.0 + xi) * (1 + zeta);
        dNdh(6, 1) = 0.125 * (1.0 + xi) * (1 + zeta);
        dNdh(7, 1) = 0.125 * (1.0 - xi) * (1 + zeta);

        // d/dzeta
        dNdh(0, 2) = -0.125 * (1.0 - xi) * (1 - eta);
        dNdh(1, 2) = -0.125 * (1.0 + xi) * (1 - eta);
        dNdh(2, 2) = -0.125 * (1.0 + xi) * (1 + eta);
        dNdh(3, 2) = -0.125 * (1.0 - xi) * (1 + eta);
        dNdh(4, 2) = 0.125 * (1.0 - xi) * (1 - eta);
        dNdh(5, 2) = 0.125 * (1.0 + xi) * (1 - eta);
        dNdh(6, 2) = 0.125 * (1.0 + xi) * (1 + eta);
        dNdh(7, 2) = 0.125 * (1.0 - xi) * (1 + eta);

    }

    /** \brie JacobianOperator
     *
     * This class is a utility to compute at a given integration point,
     * the Jacobian, its inverse, its determinant
     * and the derivatives of the shape functions in the local
     * cartesian coordinate system.
     */

     // ==================================================================
     // Struct: Jacobian3d hold the jacobian matrix data, J J^-1 |j|
     // ==================================================================
    struct Jacobian3d
    {
        // Jacobian matrix
        Matrix J = Matrix(3, 3);
        // Jacobian inverse
        Matrix invJ = Matrix(3, 3);
        // Determinant of the Jacobian matrix
        double detJ = 0.0;

        /**
        Returns false when the jacobian is numerically singular, i.e. when the
        inverse below would be meaningless.

        The test is SCALE-FREE on purpose. det(J) has units of length^3, so an
        absolute epsilon is meaningless: the same element in mm and in m differs
        by 1e9. Hadamard's inequality gives

            |det J| <= ||c1|| ||c2|| ||c3||

        for the columns c_i, so the ratio of the two is dimensionless and lies in
        [0, 1]: it is 1/sqrt(3)^3-ish for a well-conditioned cube and goes to zero
        as the element flattens. That ratio is what is compared against a
        tolerance.

        Note what is NOT rejected: a NEGATIVE determinant. It means the
        connectivity is reversed, which the rest of OpenSees tolerates (Brick does)
        and which the isotropy test of Pfefferkorn & Betsch section 5.2 actually
        requires -- two of its three prescribed node orderings are reversed, and
        the displacement magnitudes come out identical. Reversal is reported once
        per element by setDomain instead, where it is a modelling matter rather
        than an inner-loop one.
        */
        bool calculate(const vec3 X[NumGP], const Matrix& dNdh)
        {
            // jacobian (3x3)
            J(0, 0) = dNdh(0, 0) * X[0].x() + dNdh(1, 0) * X[1].x() + dNdh(2, 0) * X[2].x() + dNdh(3, 0) * X[3].x() + dNdh(4, 0) * X[4].x() + dNdh(5, 0) * X[5].x() + dNdh(6, 0) * X[6].x() + dNdh(7, 0) * X[7].x();
            J(1, 0) = dNdh(0, 0) * X[0].y() + dNdh(1, 0) * X[1].y() + dNdh(2, 0) * X[2].y() + dNdh(3, 0) * X[3].y() + dNdh(4, 0) * X[4].y() + dNdh(5, 0) * X[5].y() + dNdh(6, 0) * X[6].y() + dNdh(7, 0) * X[7].y();
            J(2, 0) = dNdh(0, 0) * X[0].z() + dNdh(1, 0) * X[1].z() + dNdh(2, 0) * X[2].z() + dNdh(3, 0) * X[3].z() + dNdh(4, 0) * X[4].z() + dNdh(5, 0) * X[5].z() + dNdh(6, 0) * X[6].z() + dNdh(7, 0) * X[7].z();


            J(0, 1) = dNdh(0, 1) * X[0].x() + dNdh(1, 1) * X[1].x() + dNdh(2, 1) * X[2].x() + dNdh(3, 1) * X[3].x() + dNdh(4, 1) * X[4].x() + dNdh(5, 1) * X[5].x() + dNdh(6, 1) * X[6].x() + dNdh(7, 1) * X[7].x();
            J(1, 1) = dNdh(0, 1) * X[0].y() + dNdh(1, 1) * X[1].y() + dNdh(2, 1) * X[2].y() + dNdh(3, 1) * X[3].y() + dNdh(4, 1) * X[4].y() + dNdh(5, 1) * X[5].y() + dNdh(6, 1) * X[6].y() + dNdh(7, 1) * X[7].y();
            J(2, 1) = dNdh(0, 1) * X[0].z() + dNdh(1, 1) * X[1].z() + dNdh(2, 1) * X[2].z() + dNdh(3, 1) * X[3].z() + dNdh(4, 1) * X[4].z() + dNdh(5, 1) * X[5].z() + dNdh(6, 1) * X[6].z() + dNdh(7, 1) * X[7].z();


            J(0, 2) = dNdh(0, 2) * X[0].x() + dNdh(1, 2) * X[1].x() + dNdh(2, 2) * X[2].x() + dNdh(3, 2) * X[3].x() + dNdh(4, 2) * X[4].x() + dNdh(5, 2) * X[5].x() + dNdh(6, 2) * X[6].x() + dNdh(7, 2) * X[7].x();
            J(1, 2) = dNdh(0, 2) * X[0].y() + dNdh(1, 2) * X[1].y() + dNdh(2, 2) * X[2].y() + dNdh(3, 2) * X[3].y() + dNdh(4, 2) * X[4].y() + dNdh(5, 2) * X[5].y() + dNdh(6, 2) * X[6].y() + dNdh(7, 2) * X[7].y();
            J(2, 2) = dNdh(0, 2) * X[0].z() + dNdh(1, 2) * X[1].z() + dNdh(2, 2) * X[2].z() + dNdh(3, 2) * X[3].z() + dNdh(4, 2) * X[4].z() + dNdh(5, 2) * X[5].z() + dNdh(6, 2) * X[6].z() + dNdh(7, 2) * X[7].z();

            // determinant
            detJ = J(0, 0) * (J(1, 1) * J(2, 2) - J(1, 2) * J(2, 1))
                - J(0, 1) * (J(1, 0) * J(2, 2) - J(1, 2) * J(2, 0))
                + J(0, 2) * (J(1, 0) * J(2, 1) - J(1, 1) * J(2, 0));

            // scale-free singularity test, see the note above
            {
                double n1 = std::sqrt(J(0,0)*J(0,0) + J(1,0)*J(1,0) + J(2,0)*J(2,0));
                double n2 = std::sqrt(J(0,1)*J(0,1) + J(1,1)*J(1,1) + J(2,1)*J(2,1));
                double n3 = std::sqrt(J(0,2)*J(0,2) + J(1,2)*J(1,2) + J(2,2)*J(2,2));
                double bound = n1 * n2 * n3;
                if (!(bound > 0.0) || std::abs(detJ) < 1.0e-12 * bound) {
                    opserr << "ASDSolidHex: singular jacobian at a gauss point. "
                        << "|detJ| = " << std::abs(detJ)
                        << " , |detJ|/(|c1||c2||c3|) = "
                        << (bound > 0.0 ? std::abs(detJ) / bound : 0.0)
                        << " (tolerance 1e-12). The element is flat or has "
                        << "coincident nodes.\n";
                    invJ.Zero();
                    return false;
                }
            }

            double mult = 1.0 / detJ;

            // inv(jacobian) (cofattori/adjoint)
            invJ(0, 0) = (J(1, 1) * J(2, 2) - J(1, 2) * J(2, 1)) * mult;
            invJ(0, 1) = -(J(0, 1) * J(2, 2) - J(0, 2) * J(2, 1)) * mult;
            invJ(0, 2) = (J(0, 1) * J(1, 2) - J(0, 2) * J(1, 1)) * mult;

            invJ(1, 0) = -(J(1, 0) * J(2, 2) - J(1, 2) * J(2, 0)) * mult;
            invJ(1, 1) = (J(0, 0) * J(2, 2) - J(0, 2) * J(2, 0)) * mult;
            invJ(1, 2) = -(J(0, 0) * J(1, 2) - J(0, 2) * J(1, 0)) * mult;

            invJ(2, 0) = (J(1, 0) * J(2, 1) - J(1, 1) * J(2, 0)) * mult;
            invJ(2, 1) = -(J(0, 0) * J(2, 1) - J(0, 1) * J(2, 0)) * mult;
            invJ(2, 2) = (J(0, 0) * J(1, 1) - J(0, 1) * J(1, 0)) * mult;

            return true;
        }

    };

    //==================================================================
    // Struct: skew_frame hold the centroid x0, jacobian at the center J0, and its inverse a= J0^-1
    // ==================================================================
    struct skew_frame
    {
        vec3   x0;          // centroid
        Matrix J0 = Matrix(3, 3);
        Matrix J0_inv = Matrix(3, 3);
        Matrix J0_invT = Matrix(3, 3);

        double detJ0 = 0.0;

        // X: 8 nodal (physical) coordinates, dNdh0: dN/d* at the center (0,0,0).
        // Returns false when J0 is numerically singular; same scale-free test as
        // Jacobian3d::calculate, see the note there.
        bool compute(const vec3 X[NumGP], const Matrix& dNdh0)
        {
            // centroid
            x0 = vec3(0, 0, 0);
            for (int i = 0; i < NumGP; i++) {
                x0 += X[i];
            }
            x0 /= 8.0;

            // jacobian at (0,0,0)
            J0(0, 0) = dNdh0(0, 0) * X[0].x() + dNdh0(1, 0) * X[1].x() + dNdh0(2, 0) * X[2].x() + dNdh0(3, 0) * X[3].x()
                + dNdh0(4, 0) * X[4].x() + dNdh0(5, 0) * X[5].x() + dNdh0(6, 0) * X[6].x() + dNdh0(7, 0) * X[7].x();
            J0(1, 0) = dNdh0(0, 0) * X[0].y() + dNdh0(1, 0) * X[1].y() + dNdh0(2, 0) * X[2].y() + dNdh0(3, 0) * X[3].y()
                + dNdh0(4, 0) * X[4].y() + dNdh0(5, 0) * X[5].y() + dNdh0(6, 0) * X[6].y() + dNdh0(7, 0) * X[7].y();
            J0(2, 0) = dNdh0(0, 0) * X[0].z() + dNdh0(1, 0) * X[1].z() + dNdh0(2, 0) * X[2].z() + dNdh0(3, 0) * X[3].z()
                + dNdh0(4, 0) * X[4].z() + dNdh0(5, 0) * X[5].z() + dNdh0(6, 0) * X[6].z() + dNdh0(7, 0) * X[7].z();

            J0(0, 1) = dNdh0(0, 1) * X[0].x() + dNdh0(1, 1) * X[1].x() + dNdh0(2, 1) * X[2].x() + dNdh0(3, 1) * X[3].x()
                + dNdh0(4, 1) * X[4].x() + dNdh0(5, 1) * X[5].x() + dNdh0(6, 1) * X[6].x() + dNdh0(7, 1) * X[7].x();
            J0(1, 1) = dNdh0(0, 1) * X[0].y() + dNdh0(1, 1) * X[1].y() + dNdh0(2, 1) * X[2].y() + dNdh0(3, 1) * X[3].y()
                + dNdh0(4, 1) * X[4].y() + dNdh0(5, 1) * X[5].y() + dNdh0(6, 1) * X[6].y() + dNdh0(7, 1) * X[7].y();
            J0(2, 1) = dNdh0(0, 1) * X[0].z() + dNdh0(1, 1) * X[1].z() + dNdh0(2, 1) * X[2].z() + dNdh0(3, 1) * X[3].z()
                + dNdh0(4, 1) * X[4].z() + dNdh0(5, 1) * X[5].z() + dNdh0(6, 1) * X[6].z() + dNdh0(7, 1) * X[7].z();

            J0(0, 2) = dNdh0(0, 2) * X[0].x() + dNdh0(1, 2) * X[1].x() + dNdh0(2, 2) * X[2].x() + dNdh0(3, 2) * X[3].x()
                + dNdh0(4, 2) * X[4].x() + dNdh0(5, 2) * X[5].x() + dNdh0(6, 2) * X[6].x() + dNdh0(7, 2) * X[7].x();
            J0(1, 2) = dNdh0(0, 2) * X[0].y() + dNdh0(1, 2) * X[1].y() + dNdh0(2, 2) * X[2].y() + dNdh0(3, 2) * X[3].y()
                + dNdh0(4, 2) * X[4].y() + dNdh0(5, 2) * X[5].y() + dNdh0(6, 2) * X[6].y() + dNdh0(7, 2) * X[7].y();
            J0(2, 2) = dNdh0(0, 2) * X[0].z() + dNdh0(1, 2) * X[1].z() + dNdh0(2, 2) * X[2].z() + dNdh0(3, 2) * X[3].z()
                + dNdh0(4, 2) * X[4].z() + dNdh0(5, 2) * X[5].z() + dNdh0(6, 2) * X[6].z() + dNdh0(7, 2) * X[7].z();

            // J inverse and determinant
            detJ0 = J0(0, 0) * (J0(1, 1) * J0(2, 2) - J0(1, 2) * J0(2, 1))
                - J0(0, 1) * (J0(1, 0) * J0(2, 2) - J0(1, 2) * J0(2, 0))
                + J0(0, 2) * (J0(1, 0) * J0(2, 1) - J0(1, 1) * J0(2, 0));

            {
                double n1 = std::sqrt(J0(0,0)*J0(0,0) + J0(1,0)*J0(1,0) + J0(2,0)*J0(2,0));
                double n2 = std::sqrt(J0(0,1)*J0(0,1) + J0(1,1)*J0(1,1) + J0(2,1)*J0(2,1));
                double n3 = std::sqrt(J0(0,2)*J0(0,2) + J0(1,2)*J0(1,2) + J0(2,2)*J0(2,2));
                double bound = n1 * n2 * n3;
                if (!(bound > 0.0) || std::abs(detJ0) < 1.0e-12 * bound) {
                    opserr << "ASDSolidHex: singular jacobian at the element centre. "
                        << "|detJ0| = " << std::abs(detJ0)
                        << " , |detJ0|/(|c1||c2||c3|) = "
                        << (bound > 0.0 ? std::abs(detJ0) / bound : 0.0)
                        << " (tolerance 1e-12). The skew frame, and with it the "
                        << "whole metric basis, cannot be built.\n";
                    J0_inv.Zero();
                    J0_invT.Zero();
                    return false;
                }
            }

            const double m = 1.0 / detJ0;

            J0_inv(0, 0) = (J0(1, 1) * J0(2, 2) - J0(1, 2) * J0(2, 1)) * m;
            J0_inv(0, 1) = -(J0(0, 1) * J0(2, 2) - J0(0, 2) * J0(2, 1)) * m;
            J0_inv(0, 2) = (J0(0, 1) * J0(1, 2) - J0(0, 2) * J0(1, 1)) * m;

            J0_inv(1, 0) = -(J0(1, 0) * J0(2, 2) - J0(1, 2) * J0(2, 0)) * m;
            J0_inv(1, 1) = (J0(0, 0) * J0(2, 2) - J0(0, 2) * J0(2, 0)) * m;
            J0_inv(1, 2) = -(J0(0, 0) * J0(1, 2) - J0(0, 2) * J0(1, 0)) * m;

            J0_inv(2, 0) = (J0(1, 0) * J0(2, 1) - J0(1, 1) * J0(2, 0)) * m;
            J0_inv(2, 1) = -(J0(0, 0) * J0(2, 1) - J0(0, 1) * J0(2, 0)) * m;
            J0_inv(2, 2) = (J0(0, 0) * J0(1, 1) - J0(0, 1) * J0(1, 0)) * m;

            // Jinv^T
            J0_invT.Zero();
            for (int ii = 0; ii < 3; ii++) {
                for (int jj = 0; jj < 3; jj++) {
                    J0_invT(ii, jj) = J0_inv(jj, ii);
                }
            }

            return true;
        }


    };

    // ==================================================================
    // Struct: evaluate dM/dx = J0^-T dM/dxi, from skew to physical domain
    // ==================================================================
    struct metric_basis
    {
        Matrix V_inv = Matrix(8, 8); // Vandermonde matrix evaluated in the skew frame
        Matrix T0 = Matrix(6, 6);    // transformation matrix for stresses
        Matrix T0_inv = Matrix(6, 6);
        vec3 Xi[NumGP]; // nodal coordinates in the skew frame
        vec3 X[NumGP];  // nodal coordinates in the physical frame
        Matrix B_trial = Matrix(6, 24);
        Matrix G_trial = Matrix(6, 12);
        Matrix G_test = Matrix(6, 12);
        skew_frame sf;
        Matrix G_at_gp[NumGP] = {
            Matrix(6,12), Matrix(6,12), Matrix(6,12), Matrix(6,12),
            Matrix(6,12), Matrix(6,12), Matrix(6,12), Matrix(6,12)
        };


        // member used to trasform gauss point into sjkew coordinate
        const double h1[8] = { +1, +1, -1, -1, -1, -1, +1, +1 };
        const double h2[8] = { +1, -1, -1, +1, -1, +1, +1, -1 };
        const double h3[8] = { +1, -1, +1, -1, +1, -1, +1, -1 };
        const double h4[8] = { -1, +1, -1, +1, +1, -1, +1, -1 };

        Vector c1 = Vector(3);
        Vector c2 = Vector(3);
        Vector c3 = Vector(3);
        Vector c4 = Vector(3);

        // scratch for orthogonalize(): |J| at the 8 gauss points and the
        // Gram-Schmidt denominators, both of which used to be recomputed inside
        // the innermost loops. This struct is already a member of the
        // ASDSolidHexGlobals singleton, so it costs no per-element memory.
        double jdet_gp[NumGP] = { 0.0 };
        double wq_gp[NumGP] = { 0.0 };   // C3 inner-product weight, see orthogonalize()
        double den_sigma[18] = { 0.0 };

        // Gram-Schmidt workspace. These used to be four std::array<Matrix,8>
        // constructed inside orthogonalize(), i.e. 32 heap Matrix allocations per
        // call. Same reasoning as above: this struct already lives in the singleton.
        Matrix S_orig[NumGP] = {
            Matrix(6,18), Matrix(6,18), Matrix(6,18), Matrix(6,18),
            Matrix(6,18), Matrix(6,18), Matrix(6,18), Matrix(6,18) };
        Matrix S_ortho[NumGP] = {
            Matrix(6,18), Matrix(6,18), Matrix(6,18), Matrix(6,18),
            Matrix(6,18), Matrix(6,18), Matrix(6,18), Matrix(6,18) };
        Matrix E_orig[NumGP] = {
            Matrix(6,12), Matrix(6,12), Matrix(6,12), Matrix(6,12),
            Matrix(6,12), Matrix(6,12), Matrix(6,12), Matrix(6,12) };
        Matrix E_ortho[NumGP] = {
            Matrix(6,12), Matrix(6,12), Matrix(6,12), Matrix(6,12),
            Matrix(6,12), Matrix(6,12), Matrix(6,12), Matrix(6,12) };

        inline vec3 computeXiBar(const double xi, const double eta, const double zeta) {
            double H1 = eta * zeta;
            double H2 = xi * zeta;
            double H3 = xi * eta;
            double H4 = xi * eta * zeta;

            Vector sum_cH(3);
            sum_cH.Zero();
            for (int i = 0; i < 3; i++)
                sum_cH(i) = c1(i) * H1 + c2(i) * H2 + c3(i) * H3 + c4(i) * H4;

            Vector corr(3);
            corr.Zero();
            corr.addMatrixVector(0.0, sf.J0_inv, sum_cH, 1.0);

            return vec3(xi + corr(0), eta + corr(1), zeta + corr(2));
        }

        inline double compute_jdet(const double xi, const double eta, const double zeta) {
            // dN/dxi for the 8 isoparametric shape functions
            static const double xi_n[8] = { -1,+1,+1,-1,-1,+1,+1,-1 };
            static const double eta_n[8] = { -1,-1,+1,+1,-1,-1,+1,+1 };
            static const double zet_n[8] = { -1,-1,-1,-1,+1,+1,+1,+1 };

            double J[3][3] = {};

            for (int i = 0; i < 8; i++) {
                double dNdxi = 0.125 * xi_n[i] * (1.0 + eta_n[i] * eta) * (1.0 + zet_n[i] * zeta);
                double dNdeta = 0.125 * eta_n[i] * (1.0 + xi_n[i] * xi) * (1.0 + zet_n[i] * zeta);
                double dNdzet = 0.125 * zet_n[i] * (1.0 + xi_n[i] * xi) * (1.0 + eta_n[i] * eta);

                J[0][0] += X[i].x() * dNdxi;   J[0][1] += X[i].x() * dNdeta;   J[0][2] += X[i].x() * dNdzet;
                J[1][0] += X[i].y() * dNdxi;   J[1][1] += X[i].y() * dNdeta;   J[1][2] += X[i].y() * dNdzet;
                J[2][0] += X[i].z() * dNdxi;   J[2][1] += X[i].z() * dNdeta;   J[2][2] += X[i].z() * dNdzet;
            }

            // det(J) by Sarrus' rule
            return J[0][0] * (J[1][1] * J[2][2] - J[1][2] * J[2][1])
                - J[0][1] * (J[1][0] * J[2][2] - J[1][2] * J[2][0])
                + J[0][2] * (J[1][0] * J[2][1] - J[1][1] * J[2][0]);
        }


        inline void orthogonalize()
        {
            const int n_sigma = 18;
            const int n_eps = 12;

            // |J| depends only on the gauss point. It used to be evaluated 2960
            // times per call for these 8 numbers (153 stress pairs x 8, 216 strain
            // pairs x 8, plus the 8 of the final loop).
            for (int gp = 0; gp < NumGP; gp++)
                jdet_gp[gp] = compute_jdet(XI[gp], ETA[gp], ZETA[gp]);

            // Quadrature weight of the inner product used by the C3
            // orthogonalisation below.
            //
            // C3 requires  int_Omega G_test^T sigma dV = 0  and, with (eq. 56)
            //     G_test = (1/jdet) * T0 * E_ortho ,     dV = jdet * dOmega_hat
            // it becomes
            //     int_hat E_ortho^T T0^T sigma dOmega_hat = 0
            // i.e. the jdet CANCELS and the condition lives in the natural domain
            // with NO jacobian weight. Both Gram-Schmidt passes used to carry
            // jdet_gp[gp] in numerator and denominator, which is a different inner
            // product. When jdet is constant over the gauss points -- a
            // parallelepiped, or any affine element -- that is a global factor and
            // cancels in num/den, which is why the cube, the sheared prism and
            // every planar-faced mesh passed the patch test. As soon as jdet varies,
            // i.e. any element with a NON-PLANAR face, the weight corrupts the
            // orthogonalisation: C3 is violated, alpha becomes non-zero for a
            // constant stress state, and the patch test fails (measured: 10 % on a
            // single warped element, 24 % on the MacNeal-Harder inner hexahedron,
            // 128 % of the stress on the 7-element distorted patch test).
            // Define ASDHEX_C3_WEIGHT_JDET to restore the old behaviour.
#ifdef ASDHEX_C3_WEIGHT_JDET
            for (int gp = 0; gp < NumGP; gp++) wq_gp[gp] = jdet_gp[gp] * WTS[gp];
#else
            for (int gp = 0; gp < NumGP; gp++) wq_gp[gp] = WTS[gp];
#endif

            // lambda: valuta sigma_v (6x18) in coordinate skew (eq. 51)
            auto eval_sigma_v = [](const vec3& xb, Matrix& sv) {
                sv.Zero();
                const double x = xb.x(), y = xb.y(), z = xb.z();

                // constant modes
                sv(0, 0) = 1.0;
                sv(1, 1) = 1.0;
                sv(2, 2) = 1.0;
                sv(3, 3) = 1.0;
                sv(4, 4) = 1.0;
                sv(5, 5) = 1.0;

                // linear bending / torsion-compatible modes
                sv(0, 6) = y;
                sv(0, 7) = z;

                sv(1, 8) = z;
                sv(1, 9) = x;

                sv(2, 10) = x;
                sv(2, 11) = y;

                sv(3, 12) = z;
                sv(4, 13) = x;
                sv(5, 14) = y;

                // incompatible bilinear stress modes
                sv(0, 15) = y * z;
                sv(1, 16) = z * x;
                sv(2, 17) = x * y;
                };

            // lambda: valuta ev (6x12) in coordinate skew (eq. 54)
            auto eval_ev = [](const vec3& xb, Matrix& ev) {
                ev.Zero();
                const double x = xb.x(), y = xb.y(), z = xb.z();

                // Columns 0-2: linear normal modes
                ev(0, 0) = x;   // eps_xx ~ xi
                ev(1, 1) = y;   // eps_yy ~ eta
                ev(2, 2) = z;   // eps_zz ~ zeta

                // Columns 3-8: linear shear modes
                ev(3, 3) = x;   ev(3, 4) = y;   // gamma_xy
                ev(4, 5) = y;   ev(4, 6) = z;   // gamma_yz
                ev(5, 7) = x;   ev(5, 8) = z;   // gamma_zx

                // Columns 9-11: bilinear volumetric modes
                ev(0, 9) = x * y;  ev(1, 9) = x * y;  ev(2, 9) = x * y;
                ev(0, 10) = y * z;  ev(1, 10) = y * z;  ev(2, 10) = y * z;
                ev(0, 11) = z * x;  ev(1, 11) = z * x;  ev(2, 11) = z * x;
                };

            // =========================================================
            // Build the stress transformation consistent with T0 (strain)
            //
            // Conventions used in the code:
            // strain Voigt = [exx eyy ezz gxy gyz gzx]^T  with g = 2e
            // stress Voigt = [sxx syy szz sxy syz szx]^T
            //
            // If eps_phys = T0 * eps_skew,
            // then sig_phys = W^{-1} * T0^{-T} * W * sig_skew
            // with W = diag(1,1,1,1/2,1/2,1/2)
            // =========================================================
            // ============================================
            // PRELIMINARY: 
            // ============================================
            // i) transform the GP from natural domain to skew domain
            // ii) evaluate sigma_v ed ev in skew domain
            // iii) define starting matrix s_v as sigma_v
            // TO DO for the II step -> starting check only the first step

            for (int gp = 0; gp < NumGP; gp++) {

                vec3 xi_bar = computeXiBar(XI[gp], ETA[gp], ZETA[gp]);

                eval_sigma_v(xi_bar, S_orig[gp]);
                eval_ev(xi_bar, E_orig[gp]);

                S_ortho[gp] = S_orig[gp];
                E_ortho[gp] = E_orig[gp];

            }

            // =========================================================
            // STEP 1: Gram-Schmidt on the stresses (eq. 55, first block)
            // HERE jdet DOES belong in both numerator and denominator
            // =========================================================
            for (int i = 0; i < n_sigma; i++) {
                // den depends only on k, not on i, and column k of S_ortho is final
                // once k has been processed -- which for every k < i it has. It used
                // to be recomputed for each of the 153 (i,k) pairs.
                if (i > 0) {
                    const int kf = i - 1;
                    double den_k = 0.0;
                    for (int gp = 0; gp < NumGP; gp++) {
                        double dot_kk = 0.0;
                        for (int r = 0; r < 6; r++)
                            dot_kk += S_ortho[gp](r, kf) * S_ortho[gp](r, kf);
                        den_k += dot_kk * wq_gp[gp];
                    }
                    den_sigma[kf] = den_k;
                }
                for (int k = 0; k < i; k++) {
                    double num = 0.0;
                    const double den = den_sigma[k];

                    for (int gp = 0; gp < NumGP; gp++) {
                        double dot_ik = 0.0;
                        for (int r = 0; r < 6; r++)
                            dot_ik += S_orig[gp](r, i) * S_ortho[gp](r, k);

                        num += dot_ik * wq_gp[gp];
                    }

                    if (std::abs(den) > 1.0e-14) {
                        double coeff = num / den;
                        for (int gp = 0; gp < NumGP; gp++) {
                            for (int r = 0; r < 6; r++) {
                                S_ortho[gp](r, i) -= coeff * S_ortho[gp](r, k);
                            }
                        }
                    }
                }
            }

            // =========================================================
            // STEP 2: Orthogonalize E_ortho against S_ortho by
            // enforcing C3 directly in skew space:
            //
            // ∫_hat E_ortho · S_ortho dΩ = 0
            //
            // Hence both numerator and denominator MUST be
            // formed without jdet.
            // =========================================================
            // S_ortho is final by now, so all 18 denominators are formed once here
            // instead of once per (j,k) pair (216 times).
            for (int k = 0; k < n_sigma; k++) {
                double den_k = 0.0;
                for (int gp = 0; gp < NumGP; gp++) {
                    double dot_kk = 0.0;
                    for (int r = 0; r < 6; r++)
                        dot_kk += S_ortho[gp](r, k) * S_ortho[gp](r, k);
                    den_k += dot_kk * wq_gp[gp];
                }
                den_sigma[k] = den_k;
            }
            for (int j = 0; j < n_eps; j++) {
                for (int k = 0; k < n_sigma; k++) {
                    double num = 0.0;
                    const double den = den_sigma[k];

                    for (int gp = 0; gp < NumGP; gp++) {
                        double dot_jk = 0.0;
                        for (int r = 0; r < 6; r++)
                            dot_jk += E_orig[gp](r, j) * S_ortho[gp](r, k);

                        num += dot_jk * wq_gp[gp];
                    }

                    if (std::abs(den) > 1.0e-14) {
                        double coeff = num / den;

                        for (int gp = 0; gp < NumGP; gp++) {
                            for (int r = 0; r < 6; r++) {
                                E_ortho[gp](r, j) -= coeff * S_ortho[gp](r, k);
                            }
                        }
                    }
                }
            }

            // =========================================================
            // FINAL TRANSFORMATION (eq. 56)
            // G_test = (1/jdet) * T0 * E_ortho
            // =========================================================
            for (int gp = 0; gp < NumGP; gp++) {
                const double jdet = jdet_gp[gp];

                G_at_gp[gp].Zero();
                G_at_gp[gp].addMatrixProduct(0.0, T0, E_ortho[gp], 1.0 / jdet);
            }


        }
        inline bool initialize_metric(const vec3 Xc[NumGP]) {
            // initialize the components of the metric functions ,
            // which are the same for the 8 gauss points
            // only depends on the element geometry and the skew frame

            // evaluate the nodal coordinates in the skew frame: x_skew = J0^-1 (x - x0)

            for (int ii = 0; ii < NumGP; ii++) {
                X[ii] = Xc[ii];
            }

            Matrix dNdh0 = Matrix(NumGP, 3);
            dshape(0.0, 0.0, 0.0, dNdh0);

            if (!sf.compute(X, dNdh0))
                return false;

            for (int ii = 0; ii < NumGP; ii++) {
                vec3 dx = X[ii] - sf.x0;
                Xi[ii] = vec3(sf.J0_inv(0, 0) * dx.x() + sf.J0_inv(0, 1) * dx.y() + sf.J0_inv(0, 2) * dx.z(),
                    sf.J0_inv(1, 0) * dx.x() + sf.J0_inv(1, 1) * dx.y() + sf.J0_inv(1, 2) * dx.z(),
                    sf.J0_inv(2, 0) * dx.x() + sf.J0_inv(2, 1) * dx.y() + sf.J0_inv(2, 2) * dx.z());
            }


            // evaluate the P vector of monomials at the skew frame pointsl  
            auto P = [](const vec3& xi, Vector& p) {
                const double x = xi.x();
                const double y = xi.y();
                const double z = xi.z();

                p(0) = 1.0;
                p(1) = x;
                p(2) = y;
                p(3) = z;
                p(4) = y * z;
                p(5) = x * z;
                p(6) = x * y;
                p(7) = x * y * z;
                };

            // build the Vandermonde matrix V_ij = P_j(Xi_i)
            Matrix V = Matrix(NumGP, NumGP);
            Vector p = Vector(NumGP);
            for (int ii = 0; ii < NumGP; ii++) {
                P(Xi[ii], p);
                for (int jj = 0; jj < NumGP; jj++) {
                    V(ii, jj) = p(jj);
                }
            }
            // invert the Vandermonde Matrix
            int info = V.Invert(V_inv);
            if (info != 0) {
                opserr << "ERROR: metric_basis::computeMetricShape Vandermonde matrix is singular\n";
            }


            // Transformation matrix for strains
            // Formula: epsilon_phys = J0^{-T} * epsilon_skew * J0^{-1}
            // In Voigt notation with [ε_xx, ε_yy, ε_zz, γ_xy, γ_yz, γ_zx]^T
            // where γ = 2ε for the shear components
            auto compute_T_eps = [](const Matrix& J0_inv) -> Matrix {
                Matrix T_eps(6, 6);
                T_eps.Zero();

                // Extract by COLUMNS of J0_inv (not rows!)
                // col 0
                double a = J0_inv(0, 0);
                double d = J0_inv(1, 0);
                double g = J0_inv(2, 0);
                // col 1
                double b = J0_inv(0, 1);
                double e = J0_inv(1, 1);
                double h = J0_inv(2, 1);
                // col 2
                double c = J0_inv(0, 2);
                double f = J0_inv(1, 2);
                double ii = J0_inv(2, 2);

                // Row 0: ε_xx → uses col0 × col0
                T_eps(0, 0) = a * a;
                T_eps(0, 1) = d * d;
                T_eps(0, 2) = g * g;
                T_eps(0, 3) = a * d;
                T_eps(0, 4) = d * g;
                T_eps(0, 5) = a * g;

                // Row 1: ε_yy → uses col1 × col1
                T_eps(1, 0) = b * b;
                T_eps(1, 1) = e * e;
                T_eps(1, 2) = h * h;
                T_eps(1, 3) = b * e;
                T_eps(1, 4) = e * h;
                T_eps(1, 5) = b * h;

                // Row 2: ε_zz → uses col2 × col2
                T_eps(2, 0) = c * c;
                T_eps(2, 1) = f * f;
                T_eps(2, 2) = ii * ii;
                T_eps(2, 3) = c * f;
                T_eps(2, 4) = f * ii;
                T_eps(2, 5) = c * ii;

                // Row 3: γ_xy = 2ε_xy → uses col0 × col1
                T_eps(3, 0) = 2.0 * a * b;
                T_eps(3, 1) = 2.0 * d * e;
                T_eps(3, 2) = 2.0 * g * h;
                T_eps(3, 3) = a * e + b * d;
                T_eps(3, 4) = d * h + e * g;
                T_eps(3, 5) = a * h + b * g;

                // Row 4: γ_yz = 2ε_yz → uses col1 × col2
                T_eps(4, 0) = 2.0 * b * c;
                T_eps(4, 1) = 2.0 * e * f;
                T_eps(4, 2) = 2.0 * h * ii;
                T_eps(4, 3) = b * f + c * e;
                T_eps(4, 4) = e * ii + f * h;
                T_eps(4, 5) = b * ii + c * h;

                // Row 5: γ_xz = 2ε_xz → uses col0 × col2
                T_eps(5, 0) = 2.0 * a * c;
                T_eps(5, 1) = 2.0 * d * f;
                T_eps(5, 2) = 2.0 * g * ii;
                T_eps(5, 3) = a * f + c * d;
                T_eps(5, 4) = d * ii + g * f;
                T_eps(5, 5) = a * ii + c * g;

                return T_eps;
                };

            // compute the transformation matrix for stresses

            T0 = compute_T_eps(sf.J0_inv);
            info = T0.Invert(T0_inv);
            if (info != 0) {
                opserr << "ERROR: metric_basis::T0^-1 is singular\n";
            }


            // compute usefull memebers

            c1.Zero(); c2.Zero(); c3.Zero(); c4.Zero();

            for (int ii = 0; ii < 8; ii++) {
                c1(0) += X[ii].x() * h1[ii];
                c1(1) += X[ii].y() * h1[ii];
                c1(2) += X[ii].z() * h1[ii];

                c2(0) += X[ii].x() * h2[ii];
                c2(1) += X[ii].y() * h2[ii];
                c2(2) += X[ii].z() * h2[ii];

                c3(0) += X[ii].x() * h3[ii];
                c3(1) += X[ii].y() * h3[ii];
                c3(2) += X[ii].z() * h3[ii];

                c4(0) += X[ii].x() * h4[ii];
                c4(1) += X[ii].y() * h4[ii];
                c4(2) += X[ii].z() * h4[ii];
            }
            // after the loop over the c_A, before *= 1/8
            c1 *= (1.0 / 8.0);
            c2 *= (1.0 / 8.0);
            c3 *= (1.0 / 8.0);
            c4 *= (1.0 / 8.0);

            orthogonalize();

            return true;

        }


        // define a method to compute the metric shape funcitons M
        inline void computeMetricShape(const double xiP, const double etaP, const double zitaP, const double jdet) {
            // transform gauss point from natural to skew domain
            double H1 = etaP * zitaP;
            double H2 = xiP * zitaP;
            double H3 = xiP * etaP;
            double H4 = xiP * etaP * zitaP;

            Vector sum_cH(3);
            sum_cH.Zero();

            for (int ii = 0; ii < 3; ii++) {
                sum_cH(ii) = c1(ii) * H1 + c2(ii) * H2 + c3(ii) * H3 + c4(ii) * H4;
            }

            Vector correction(3);
            correction.Zero();
            correction.addMatrixVector(0.0, sf.J0_inv, sum_cH, 1.0);

            Vector XGP_skew(3);
            XGP_skew.Zero();
            XGP_skew(0) = xiP + correction(0);
            XGP_skew(1) = etaP + correction(1);
            XGP_skew(2) = zitaP + correction(2);

            // evaluate the monomial
            auto Pmono = [](const Vector& xi, Vector& p) {
                const double x = xi(0);
                const double y = xi(1);
                const double z = xi(2);

                p(0) = 1.0;
                p(1) = x;
                p(2) = y;
                p(3) = z;
                p(4) = y * z;
                p(5) = x * z;
                p(6) = x * y;
                p(7) = x * y * z;
                };

            // evaluate the monomial vector in the gauss point in the skew frame
            Vector p = Vector(NumGP);
            Pmono(XGP_skew, p);

            // metric shape function
            Vector M = Vector(NumGP);
            M.Zero();
            for (int ii = 0; ii < NumGP; ii++) {
                for (int jj = 0; jj < NumGP; jj++) {
                    M(ii) += p(jj) * V_inv(jj, ii);
                }
            }

            // compute the derivatives of monomials wrt skew coordinates
            auto dPdxi = [](const Vector& xi, Matrix& dpdxi) {
                const double x = xi(0);
                const double y = xi(1);
                const double z = xi(2);

                dpdxi.Zero();
                dpdxi(0, 1) = 1.0;
                dpdxi(1, 2) = 1.0;
                dpdxi(2, 3) = 1.0;

                dpdxi(1, 4) = z;
                dpdxi(2, 4) = y;

                dpdxi(0, 5) = z;
                dpdxi(2, 5) = x;

                dpdxi(0, 6) = y;
                dpdxi(1, 6) = x;

                dpdxi(0, 7) = y * z;
                dpdxi(1, 7) = x * z;
                dpdxi(2, 7) = x * y;
                };

            // compute the gradient in the skew coordinates
            Matrix dpdxi = Matrix(3, NumGP);
            dpdxi.Zero();
            dPdxi(XGP_skew, dpdxi);

            // dM/dxi = dP/dxi * V^-1
            Matrix dMdxi = Matrix(3, NumGP);
            dMdxi.Zero();
            dMdxi.addMatrixProduct(0.0, dpdxi, V_inv, 1.0);

            // compute the gradient in physical coordinates: dM/dx = J0^-T dM/dxi
            Matrix dMdx = Matrix(3, NumGP);
            dMdx.Zero();
            dMdx.addMatrixTransposeProduct(0.0, sf.J0_inv, dMdxi, 1.0);


            // compute the B trial matrix in metric coordinates
            B_trial.Zero();

            for (int ii = 0; ii < NumGP; ii++) {
                const int idx1 = ii * 3;
                const int idx2 = ii * 3 + 1;
                const int idx3 = ii * 3 + 2;

                B_trial(0, idx1) = dMdx(0, ii);  // dM/dx
                B_trial(1, idx2) = dMdx(1, ii);  // dM/dy
                B_trial(2, idx3) = dMdx(2, ii);  // dM/dz

                B_trial(3, idx1) = dMdx(1, ii);  // gxy
                B_trial(3, idx2) = dMdx(0, ii);
                B_trial(4, idx2) = dMdx(2, ii);  // gyz
                B_trial(4, idx3) = dMdx(1, ii);
                B_trial(5, idx1) = dMdx(2, ii);  // gzx
                B_trial(5, idx3) = dMdx(0, ii);
            }

            // compute G_trial for the 9 quadratic + 3 volumetric modes
            auto dMTILDEDXI = [](const vec3& xiGP,
                const vec3 Xi_nodes[8],
                const Matrix& dMdxi) -> Matrix {
                    Matrix dMtildedxi(3, 3);
                    dMtildedxi.Zero();

                    double x = xiGP.x(), y = xiGP.y(), z = xiGP.z();

                    dMtildedxi(0, 0) = 2.0 * x;
                    dMtildedxi(1, 1) = 2.0 * y;
                    dMtildedxi(2, 2) = 2.0 * z;

                    for (int i = 0; i < 8; i++) {
                        double xi_n = Xi_nodes[i].x();
                        double eta_n = Xi_nodes[i].y();
                        double zet_n = Xi_nodes[i].z();

                        double m1 = xi_n * xi_n;
                        double m2 = eta_n * eta_n;
                        double m3 = zet_n * zet_n;

                        for (int r = 0; r < 3; r++) {
                            dMtildedxi(r, 0) -= dMdxi(r, i) * m1;
                            dMtildedxi(r, 1) -= dMdxi(r, i) * m2;
                            dMtildedxi(r, 2) -= dMdxi(r, i) * m3;
                        }
                    }
                    return dMtildedxi;
                };

            vec3 xiGP_vec(XGP_skew(0), XGP_skew(1), XGP_skew(2));
            Matrix dMtildedxi = dMTILDEDXI(xiGP_vec, Xi, dMdxi);

            Matrix dMtildedx(3, 3);
            dMtildedx.Zero();
            dMtildedx.addMatrixProduct(0.0, sf.J0_invT, dMtildedxi, 1.0);

            G_trial.Zero();

            // Part 1: quadratic Wilson modes (columns 0-8)
            for (int j = 0; j < 3; j++) {
                double d1 = dMtildedx(0, j);
                double d2 = dMtildedx(1, j);
                double d3 = dMtildedx(2, j);

                int c1 = 3 * j;
                int c2 = 3 * j + 1;
                int c3 = 3 * j + 2;

                G_trial(0, c1) = d1;
                G_trial(1, c2) = d2;
                G_trial(2, c3) = d3;
                G_trial(3, c1) = d2;  G_trial(3, c2) = d1;
                G_trial(4, c2) = d3;  G_trial(4, c3) = d2;
                G_trial(5, c1) = d3;  G_trial(5, c3) = d1;
            }

            // Part 2: volumetric modes (columns 9-11)
            // Eq. (53): J0^{-T} [ p(xi,eta,zeta) * I ] J0^{-1}
            // In engineering Voigt notation => use the full P = J0^{-T} J0^{-1},
            // including the shear terms 2*P12, 2*P23, 2*P13.
            double xi = XGP_skew(0);
            double eta = XGP_skew(1);
            double zeta = XGP_skew(2);

            Matrix P(3, 3);
            P.Zero();
            P.addMatrixProduct(0.0, sf.J0_invT, sf.J0_inv, 1.0);

            double prod[3] = { xi * eta, eta * zeta, zeta * xi };

            for (int m = 0; m < 3; m++) {
                int c = 9 + m;
                double pm = prod[m];

                G_trial(0, c) = pm * P(0, 0);          // eps_xx
                G_trial(1, c) = pm * P(1, 1);          // eps_yy
                G_trial(2, c) = pm * P(2, 2);          // eps_zz
                G_trial(3, c) = pm * 2.0 * P(0, 1);    // gamma_xy
                G_trial(4, c) = pm * 2.0 * P(1, 2);    // gamma_yz
                G_trial(5, c) = pm * 2.0 * P(0, 2);    // gamma_zx
            }
        }

    };

    // ==================================================================
    // Function: compute the B_test matrix for the isoparametric test functions based on dN/dx
    // ===================================================================
    // fills B_test in place. It used to return Matrix(6,24) BY VALUE, i.e. one heap
    // allocation and one copy per gauss point per calculateAll.
    inline void compute_B_test(const Matrix& dNdh, const Matrix& invJ, Matrix& B_test)
    {
        B_test.Zero();

        for (int a = 0; a < NumGP; a++) {
            const int c0 = 3 * a, c1 = c0 + 1, c2 = c0 + 2;

            // CORRECT: use the columns of invJ (= rows of invJ^T)
            const double dNa_dx = invJ(0, 0) * dNdh(a, 0) + invJ(1, 0) * dNdh(a, 1) + invJ(2, 0) * dNdh(a, 2);
            const double dNa_dy = invJ(0, 1) * dNdh(a, 0) + invJ(1, 1) * dNdh(a, 1) + invJ(2, 1) * dNdh(a, 2);
            const double dNa_dz = invJ(0, 2) * dNdh(a, 0) + invJ(1, 2) * dNdh(a, 1) + invJ(2, 2) * dNdh(a, 2);

            B_test(0, c0) = dNa_dx;
            B_test(1, c1) = dNa_dy;
            B_test(2, c2) = dNa_dz;
            B_test(3, c0) = dNa_dy;  B_test(3, c1) = dNa_dx;
            B_test(4, c1) = dNa_dz;  B_test(4, c2) = dNa_dy;
            B_test(5, c0) = dNa_dz;  B_test(5, c2) = dNa_dx;
        }
    }
    // ==================================================================
    // Struct: copmute the matrices and modes for the element
    // ==================================================================
    /** \brief ASDSolidHexGlobals
     *
     * This singleton class stores some data for the hexaedro calculations that
     * can be statically instantiated to avoid useless re-allocations
     *
     */
    class ASDSolidHexGlobals
    {
    private:
        ASDSolidHexGlobals() = default;

    public:

        static constexpr int nnode = 8;     // nodes
        static constexpr int ndofn = 3;     // DOFs per node (ux,uy,uz)
        static constexpr int ndofe = nnode * ndofn; // 24
        static constexpr int nvoigt = 6;    // (xx,yy,zz,xy,yz,zx)
        static constexpr int nq = 12;       // selected internal EAS DOFs

        Matrix J = Matrix(3, 3);
        Matrix invJ = Matrix(3, 3);
        double detJ = 0.0;

        // metric basis. The skew frame it needs is metric_basis::sf, computed
        // inside initialize_metric; the singleton no longer carries a second one.
        metric_basis mb;
        //compute_blocks_hex8 compute_blocks;

        // workspace for coordinates and shape functions
        vec3   X[nnode];         // physical nodal coordinates
        Vector N = Vector(nnode);          // shape Hex8
        Matrix dNdh = Matrix(nnode, 3);     // dN/d(hat)
        Matrix dNdx = Matrix(nnode, 3);     // dN/dx

        // degrees of freedom/displacements
        Vector UG = Vector(ndofe); // global displacements
        Vector UL = Vector(ndofe); // (if a local frame or a copy is needed)

        // B-matrices (Voigt)
        Matrix B_test = Matrix(nvoigt, ndofe);      // test: isoparametric
        Matrix B_trial = Matrix(nvoigt, ndofe);      // trial: Petrov metrics
        Matrix G_trial = Matrix(nvoigt, nq);         // EAS trial
        Matrix G_test = Matrix(nvoigt, nq);         // EAS test (orthonormalized)

        // material data and state variables (in Voigt form)
        Matrix C = Matrix(nvoigt, nvoigt); // 3D elastic tangent matrix
        Matrix C0 = Matrix(nvoigt, nvoigt); // initial (undamaged) tangent at GP
        Vector eps = Vector(nvoigt);         // strain at GP
        Vector sig = Vector(nvoigt);         // stress at GP

        // Petrov blocks (non-symmetric)
        Matrix k_uu = Matrix(ndofe, ndofe); // int b_v^T C b_u
        Matrix k_qq = Matrix(nq, nq);       // int bq_test^T C bq_trial
        Matrix k_qq_inv = Matrix(nq, nq);     // inv(k_qq)

        // the A-op regularizer: the same Petrov block built on the INITIAL material
        // tangent, and k_qq + s*k_qq0, the operator actually inverted. Scratch, like
        // everything else here: accumulated and consumed within one calculateAll.
        Matrix k_qq0 = Matrix(nq, nq);      // int bq_test^T C0 bq_trial  (A-op)
        Matrix k_qq0_sym = Matrix(nq, nq);  // int bq_trial^T C0 bq_trial (A-pen)
        Matrix k_qq_stab = Matrix(nq, nq);  // k_qq + s*k_qq0 + p*k_qq0_sym
        Matrix C0tBq = Matrix(nvoigt, nq);  // C0 G_trial, shared by both

        // global element LHS/RHS matrices and vectors
        Matrix LHS = Matrix(ndofe, ndofe); // LHS matrix (tangent stiffness)
        Matrix LHS_initial = Matrix(ndofe, ndofe); // LHS matrix (initial stiffness)
        Matrix LHS_mass = Matrix(ndofe, ndofe); // LHS matrix (mass matrix)

        // scratch for the PG-EAS blocks. These used to be constructed inside the
        // gauss loop (3 per gauss point, so 24 heap Matrix per calculateAll) and
        // after it (2 more), while the pre-allocated members below sat unused.
        // Same arrangement as ASDShellQ4Globals' B1TD / BQTD / DBQ.
        Matrix BtC = Matrix(ndofe, nvoigt);   // B_test^T C
        Matrix CtBu = Matrix(nvoigt, ndofe);  // C B_trial
        Matrix CtBq = Matrix(nvoigt, nq);     // C G_trial
        Matrix Kuq_Kqqinv = Matrix(ndofe, nq);        // K_uq * inv(K_qq)
        Matrix Kuq_Kqqinv_Kqu = Matrix(ndofe, ndofe); // ... * K_qu
        Vector RHS = Vector(ndofe); // RHS vector (residual vector)
        Vector RHS_winertia = Vector(ndofe); // RHS vector (residual vector with inertia terms)


    public:
        static ASDSolidHexGlobals& instance() {
            static ASDSolidHexGlobals _instance;
            return _instance;
        }
    };


    // build_b_v: B (6x24) for the isoparametric test functions based on dN/dx (8x3)
    // Voigt: [exx eyy ezz gxy gyz gzx]^T, shear = engineering gamma
    inline void build_b_v(const Matrix& dNdx, Matrix& B_v)
    {
        B_v.Zero();
        // column mapping: for node a, dof = (ux,uy,uz) -> col = 3*a + {0,1,2}
        for (int a = 0; a < NumGP; a++) {
            const int c0 = 3 * a;     // ux
            const int c1 = c0 + 1;  // uy
            const int c2 = c0 + 2;  // uz

            const double dNa_dx = dNdx(a, 0);
            const double dNa_dy = dNdx(a, 1);
            const double dNa_dz = dNdx(a, 2);

            // normal components
            B_v(0, c0) += dNa_dx;           // exx <- ux,x
            B_v(1, c1) += dNa_dy;           // eyy <- uy,y
            B_v(2, c2) += dNa_dz;           // ezz <- uz,z

            // shear (engineering)
            B_v(3, c0) += dNa_dy;           // gxy <- ux,y
            B_v(3, c1) += dNa_dx;           // gxy <- uy,x

            B_v(4, c1) += dNa_dz;           // gyz <- uy,z
            B_v(4, c2) += dNa_dy;           // gyz <- uz,y

            B_v(5, c2) += dNa_dx;           // gzx <- uz,x
            B_v(5, c0) += dNa_dz;           // gzx <- ux,z
        }
    }

    // build_b_u: B (6x24) for the metric trial functions based on dM/dx (3x8)
    // dMdx has shape (3x8): rows = (d/dx, d/dy, d/dz), columns = i=0..7 (nodal metric functions)

    inline void build_b_u(const Matrix& dMdx, Matrix& b_u)
    {
        // build the B_u matrix
        b_u.Zero();

        //  dMdx_T (8x3), without allocating extra memory
        for (int a = 0; a < NumGP; a++) {
            const int c0 = 3 * a;
            const int c1 = c0 + 1;
            const int c2 = c0 + 2;

            const double dMa_dx = dMdx(0, a);
            const double dMa_dy = dMdx(1, a);
            const double dMa_dz = dMdx(2, a);

            // normal components
            b_u(0, c0) = dMa_dx;           // exx
            b_u(1, c1) = dMa_dy;           // eyy
            b_u(2, c2) = dMa_dz;           // ezz

            // shear (engineering)
            b_u(3, c0) = dMa_dy;           // gxy
            b_u(3, c1) = dMa_dx;

            b_u(4, c1) = dMa_dz;           // gyz
            b_u(4, c2) = dMa_dy;

            b_u(5, c2) = dMa_dx;           // gzx
            b_u(5, c0) = dMa_dz;
        }
    }

}

/**
Per-element cache of the reference-geometry metric basis.

metric_basis::initialize_metric() builds the skew frame, the Vandermonde inverse
and the four bilinear correction vectors, then calls orthogonalize(), which runs
the two-step Gram-Schmidt of Pfefferkorn & Betsch eq. (55) over 18 stress modes
and 12 enhanced modes at 8 gauss points and leaves the result in G_at_gp. Every
one of those quantities is a function of the REFERENCE geometry alone:
calculateAll fills X from Node::getCrds() in the linear and the corotational path
alike (the corotational frame is measured from the reference configuration, so
R0 = I and the reference local coordinates are the global ones), and nothing in
the chain reads a deformed coordinate. Rebuilding it every Newton iteration cost
a measured 30.5 % of calculateAll in a corotational run and 33.0 % in a linear
one, so it is built once per element and copied back in afterwards.

Only the fields that are read AFTER initialize_metric returns are stored:
computeMetricShape() needs V_inv, sf, c1..c4 and Xi (the last one through the
dMTILDEDXI lambda, which is easy to miss), and the gauss loop needs G_at_gp. T0,
X, jdet_gp, wq_gp, den_sigma and the four Gram-Schmidt workspaces are internal to
initialize_metric/orthogonalize and stay in the shared singleton.

Cost: 707 doubles, about 5.7 kB per element -- the same order as the EAS state
the element already carries (792 doubles). The copy back in is 5.5 kB against the
several tens of thousands of flops it replaces.
*/
struct ASDSolidHexRefMetric
{
    bool valid = false;
    Matrix V_inv = Matrix(8, 8);
    vec3 Xi[NumGP];        // nodal coordinates in the skew frame
    skew_frame sf;
    Vector c1 = Vector(3);
    Vector c2 = Vector(3);
    Vector c3 = Vector(3);
    Vector c4 = Vector(3);
    Matrix G_at_gp[NumGP] = {
        Matrix(6,12), Matrix(6,12), Matrix(6,12), Matrix(6,12),
        Matrix(6,12), Matrix(6,12), Matrix(6,12), Matrix(6,12)
    };
};

namespace {

    inline void asdhex_ref_save(const metric_basis& mb, ASDSolidHexRefMetric& r)
    {
        r.V_inv = mb.V_inv;
        for (int i = 0; i < NumGP; i++) r.Xi[i] = mb.Xi[i];
        r.sf = mb.sf;
        r.c1 = mb.c1; r.c2 = mb.c2; r.c3 = mb.c3; r.c4 = mb.c4;
        for (int gp = 0; gp < NumGP; gp++)
            r.G_at_gp[gp] = mb.G_at_gp[gp];
        r.valid = true;
    }

    inline void asdhex_ref_load(const ASDSolidHexRefMetric& r, metric_basis& mb)
    {
        mb.V_inv = r.V_inv;
        for (int i = 0; i < NumGP; i++) mb.Xi[i] = r.Xi[i];
        mb.sf = r.sf;
        mb.c1 = r.c1; mb.c2 = r.c2; mb.c3 = r.c3; mb.c4 = r.c4;
        for (int gp = 0; gp < NumGP; gp++)
            mb.G_at_gp[gp] = r.G_at_gp[gp];
    }

}

ASDSolidHex::ASDSolidHex()
    : Element(0, ELE_TAG_ASDSolidHex)
    , m_use_corotational(false)
    , m_transformation(nullptr)
    , m_load(nullptr)
    , m_initialized(false)
{
    // This is the constructor used by FEM_ObjectBroker before recvSelf.
    // m_material[] MUST be nulled here: the destructor does
    // "if (m_material[gp]) delete m_material[gp]" and would otherwise delete
    // indeterminate pointers. m_eas is allocated unconditionally because
    // commitState/revertToLastCommit/revertToStart/calculateAll all dereference
    // it without a null check.
    m_node_ids = ID(NumNodes);
    for (int i = 0; i < NumNodes; i++)
        nodePtrs[i] = nullptr;
    for (int i = 0; i < NumNodes; i++)
        m_material[i] = nullptr;
    for (int i = 0; i < NumGP; i++)
        m_damping[i] = nullptr;
    m_eas = new EASData();
}

ASDSolidHex::ASDSolidHex(
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
    bool corotational,
    Damping* damping,
    const double* body,
    double eas_stab,
    double eas_penalty,
    bool eas_auto
#ifdef ASDHEX_EAS_IMPERF
    , double imperfection
#endif
    )
    : Element(tag, ELE_TAG_ASDSolidHex)
    , m_eas_stab(eas_stab)
    , m_eas_penalty(eas_penalty)
    , m_eas_auto(eas_auto)
#ifdef ASDHEX_EAS_IMPERF
    , m_imperfection(imperfection)
#endif
    , m_use_corotational(corotational)
    , m_transformation(corotational ? new ASDSolidHexCorotationalTransformation() : nullptr)
    , m_load(nullptr)
    , m_initialized(false)
{
    // save node ids
    m_node_ids = ID(8);
    m_node_ids(0) = node1;
    m_node_ids(1) = node2;
    m_node_ids(2) = node3;
    m_node_ids(3) = node4;
    m_node_ids(4) = node5;
    m_node_ids(5) = node6;
    m_node_ids(6) = node7;
    m_node_ids(7) = node8;

    // allocate EAS variables
    m_eas = new EASData();

    for (int i = 0; i < NumGP; i++) {
        nodePtrs[i] = nullptr;
    }
    // copy ND material
    if (mat == nullptr) {
        opserr << "ASDSolidHex::ASDSolidHex - NULL NDMaterial pointer passed to constructor\n";
        exit(-1);
    }

    for (int gp = 0; gp < NumGP; gp++) {
        m_material[gp] = mat->getCopy("ThreeDimensional");
        if (m_material[gp] == nullptr) {
            opserr << "ASDSolidHex::ASDSolidHex - failed to getCopy(ThreeDimensional) for GP " << gp << endln;
            exit(-1);
        }
    }

    // one Damping copy per gauss point, mirroring ASDShellQ4 (one per section)
    for (int gp = 0; gp < NumGP; gp++)
        m_damping[gp] = nullptr;
    if (damping) {
        for (int gp = 0; gp < NumGP; gp++) {
            m_damping[gp] = damping->getCopy();
            if (m_damping[gp] == nullptr) {
                opserr << "ASDSolidHex::ASDSolidHex - failed to get copy of damping for GP " << gp << endln;
                exit(-1);
            }
        }
    }

    // body force per unit mass (gravity), used by addLoad
    if (body) {
        for (int i = 0; i < 3; i++)
            m_body[i] = body[i];
    }
}

ASDSolidHex::~ASDSolidHex()
{
    // clean up material
    for (int gp = 0; gp < NumGP; gp++) {
        if (m_material[gp]) {
            delete m_material[gp];
            m_material[gp] = nullptr;
        }
    }

    // clean up damping
    for (int gp = 0; gp < NumGP; gp++) {
        if (m_damping[gp]) {
            delete m_damping[gp];
            m_damping[gp] = nullptr;
        }
    }

    // clean up coordinate transformation.
    // NOTE: gate on the POINTER, not on m_use_corotational -- the two can
    // disagree on a partially constructed element.
    if (m_transformation) {
        delete m_transformation;
        m_transformation = nullptr;
    }

    // clean up load vector
    if (m_load) {
        delete m_load;
        m_load = nullptr;
    }

    // clean up the enhanced-assumed-strain data (this used to leak)
    if (m_eas) {
        delete m_eas;
        m_eas = nullptr;
    }

    // clean up the reference metric basis cache
    if (m_ref) {
        delete m_ref;
        m_ref = nullptr;
    }
}

void  ASDSolidHex::setDomain(Domain* theDomain)
{
    // if domain is null
    if (theDomain == nullptr) {
        for (int i = 0; i < NumNodes; i++)
            nodePtrs[i] = nullptr;
        m_initialized = false;
        // the nodes are gone, so the cached reference geometry means nothing
        if (m_ref)
            m_ref->valid = false;
        // call base class implementation
        DomainComponent::setDomain(theDomain);
        return;
    }

    // node pointers.
    // NOTE: on failure we still fall through to DomainComponent::setDomain --
    // returning early used to leave the element half-built (base class never
    // informed) which crashes later instead of at the point of the error.
    // the node pointers are about to be (re)bound, so anything cached from the
    // previous reference geometry is stale. This is the only place nodePtrs can
    // change, so it is the only place the cache has to be dropped.
    if (m_ref)
        m_ref->valid = false;

    bool ok = true;
    for (int i = 0; i < NumNodes; i++) {
        Node* node = theDomain->getNode(m_node_ids(i));
        if (node == nullptr) {
            opserr << "ASDSolidHex::setDomain -- node " << m_node_ids(i)
                << " not found in the domain\n";
            ok = false;
            break;
        }
        if (node->getNumberDOF() != 3) {
            opserr << "ASDSolidHex::setDomain -- node " << m_node_ids(i)
                << " has " << node->getNumberDOF() << " DOFs, expected 3\n";
            ok = false;
            break;
        }
        nodePtrs[i] = node;
    }

    if (!ok) {
        for (int i = 0; i < NumNodes; i++)
            nodePtrs[i] = nullptr;
        DomainComponent::setDomain(theDomain);
        return;
    }

    // if using the corotational formulation set up the domain
    if (m_use_corotational)
        m_transformation->setDomain(theDomain, m_node_ids, m_initialized);

    // damping: 6 stress components for a 3D solid (the shell passes 8 because it
    // works with section resultants)
    for (int gp = 0; gp < NumGP; gp++) {
        if (m_damping[gp] && m_damping[gp]->setDomain(theDomain, 6)) {
            opserr << "ASDSolidHex::setDomain -- error initializing damping for GP " << gp << endln;
            return;
        }
    }

    // One-time geometry validation. Both jacobian routines now refuse to invert a
    // singular matrix, but they do so from inside the newton loop, where the only
    // sensible response is to abort the step. Checking here instead means a bad
    // mesh is named at model-build time.
    //
    // A NEGATIVE determinant is reported but not rejected. It means the
    // connectivity is reversed: the stiffness comes out as -K, so displacements
    // come back with the wrong sign while their magnitudes are unchanged (this is
    // exactly what two of the three node orderings of the Pfefferkorn & Betsch
    // isotropy test do, deliberately). It is a modelling error worth one line of
    // output, not a reason to refuse the element.
    {
        vec3 Xc[NumGP];
        for (int i = 0; i < NumNodes; i++) {
            const Vector& crd = nodePtrs[i]->getCrds();
            Xc[i] = vec3(crd(0), crd(1), crd(2));
        }
        Matrix dNdh(NumGP, 3);
        Jacobian3d jac;
        int n_neg = 0;
        bool singular = false;
        for (int gp = 0; gp < NumGP; gp++) {
            dshape(XI[gp], ETA[gp], ZETA[gp], dNdh);
            if (!jac.calculate(Xc, dNdh)) {
                singular = true;
                break;
            }
            if (jac.detJ < 0.0)
                n_neg++;
        }
        if (singular) {
            opserr << "ASDSolidHex::setDomain - element " << this->getTag()
                << " has a singular geometry (see above) and cannot be used\n";
        }
        else if (n_neg == NumGP) {
            opserr << "ASDSolidHex::setDomain - element " << this->getTag()
                << " has REVERSED connectivity (detJ < 0 at all 8 gauss points). "
                << "Displacement magnitudes are unaffected but signs are flipped; "
                << "swap the two faces of the node list to fix it.\n";
        }
        else if (n_neg > 0) {
            opserr << "ASDSolidHex::setDomain - element " << this->getTag()
                << " has detJ < 0 at " << n_neg << " of 8 gauss points, i.e. it is "
                << "folded rather than merely reversed. Results are meaningless.\n";
        }
    }

#ifdef ASDHEX_EAS_IMPERF
    // Filled here, not in the constructor: the FEM_ObjectBroker path builds the
    // element with tag 0 and recvSelf installs the real tag afterwards, so a
    // constructor-time hash would key off the wrong number on every process but
    // the one that owns the element. setDomain runs after recvSelf in both paths.
    for (int gp = 0; gp < NumGP; ++gp)
        m_imperfection_f[gp] = 1.0 + m_imperfection * imperfectionUnit(this->getTag(), gp);
#endif

    // only if not already initialized from recvSelf: recvSelf restores the
    // EAS state, the initial displacement offset and the transformation's
    // internal data, so re-initializing here would throw them away.
    if (!m_initialized) {
        initializePG_EAS();
        captureInitialDisp();
        m_initialized = true;
    }

    // call base class implementation
    DomainComponent::setDomain(theDomain);
}

void ASDSolidHex::captureInitialDisp()
{
    for (int i = 0; i < NumNodes; ++i) {
        const Vector& iU = nodePtrs[i]->getTrialDisp();
        for (int j = 0; j < 3; ++j)
            m_U0(3 * i + j) = iU(j);
    }
}

void ASDSolidHex::onActivate()
{
    // Re-capture the initial displacement offset at the current configuration, so a
    // staged element is born strain free: the element-side one for the linear path,
    // the transformation's own for the corotational path (the reference frame is
    // not touched: it comes from the nodal coordinates alone). setDomain() is
    // deliberately NOT re-run: it would re-initialize the damping and re-validate
    // the geometry for nothing.
    captureInitialDisp();
    if (m_use_corotational)
        m_transformation->forceCaptureInitialDisp();
    // The PG-EAS internal DOFs are baselined on the offset displacements, so they
    // have to follow the new offset or they would keep a stale baseline - same as
    // ASDShellQ4 re-running AGQIinitialize() on activation. The local displacements
    // at the activation configuration are exactly the new offset, so the correct
    // baseline is zero, which is what initializePG_EAS() seeds.
    initializePG_EAS();
    this->update();
}

void ASDSolidHex::onDeactivate()
{
}

int ASDSolidHex::setDamping(Domain* theDomain, Damping* damping)
{
    if (theDomain && damping) {
        for (int gp = 0; gp < NumGP; gp++) {
            if (m_damping[gp])
                delete m_damping[gp];
            m_damping[gp] = damping->getCopy();
            if (m_damping[gp] == nullptr) {
                opserr << "ASDSolidHex::setDamping - failed to get copy of damping\n";
                return -1;
            }
            if (m_damping[gp]->setDomain(theDomain, 6)) {
                opserr << "ASDSolidHex::setDamping - error initializing damping\n";
                return -2;
            }
        }
    }
    return 0;
}

void ASDSolidHex::Print(OPS_Stream& s, int flag)
{
    if (flag == OPS_PRINT_CURRENTSTATE) {
        s << "ASDSolidHex, element id: " << this->getTag() << endln;
        s << "   kinematics: " << (m_use_corotational ? "corotational" : "linear") << endln;
        s << "   Connected external nodes: ";
        for (int i = 0; i < NumNodes; i++)
            s << m_node_ids(i) << " ";
        s << endln;
        s << "   Material: " << (m_material[0] ? m_material[0]->getTag() : -1) << endln;
        s << "   EAS stabilization floor: " << m_eas_stab << endln;
        s << "   EAS penalty: " << m_eas_penalty << endln;
        s << "   EAS automatic arming: " << (m_eas_auto ? "on" : "off") << endln;
#ifdef ASDHEX_EAS_IMPERF
        s << "   imperfection amplitude: " << m_imperfection << endln;
        s << "   imperfection factors (internal gauss order):";
        for (int gp = 0; gp < NumGP; ++gp)
            s << " " << m_imperfection_f[gp];
        s << endln;
#endif
    }
    else if (flag == OPS_PRINT_PRINTMODEL_JSON) {
        s << "\t\t\t{";
        s << "\"name\": " << this->getTag() << ", ";
        s << "\"type\": \"ASDSolidHex\", ";
        s << "\"nodes\": [";
        for (int i = 0; i < NumNodes; i++) {
            s << m_node_ids(i);
            if (i < NumNodes - 1) s << ", ";
        }
        s << "], ";
        s << "\"corotational\": " << (m_use_corotational ? "true" : "false") << ", ";
        s << "\"easStab\": " << m_eas_stab << ", ";
        s << "\"easPenalty\": " << m_eas_penalty << ", ";
        s << "\"easAuto\": " << (m_eas_auto ? "true" : "false") << ", ";
        s << "\"material\": \"" << (m_material[0] ? m_material[0]->getTag() : -1) << "\"}";
    }
}

int  ASDSolidHex::getNumExternalNodes() const
{
    return 8;
}

const ID& ASDSolidHex::getExternalNodes()
{
    return m_node_ids;
}

Node**
ASDSolidHex::getNodePtrs(void)
{
    return nodePtrs;
}

int  ASDSolidHex::getNumDOF()
{
    return 24;
}

int  ASDSolidHex::commitState()
{
    int success = 0;

    if(m_use_corotational)
		m_transformation->commit();

    // save the enhanced state. NOTE: m_eas->U holds the last trial local
    // displacement (written by updatePG_EAS every iteration); committing means
    // snapshotting it, NOT overwriting it -- doing the latter reset U to its
    // value at initialization, so the next step's dU became a TOTAL rather than
    // an incremental displacement.
    m_eas->alpha_commit = m_eas->alpha;
    m_eas->U_converged = m_eas->U;

    // NDMaterial at gauss Points
    for (int gp = 0; gp < NumGP; ++gp) {
        if (m_material[gp])
            success += m_material[gp]->commitState();
    }

    // damping
    for (int gp = 0; gp < NumGP; ++gp) {
        if (m_damping[gp])
            success += m_damping[gp]->commitState();
    }

    // done
    return success;
}

int  ASDSolidHex::revertToLastCommit()
{
    int success = 0;

    if(m_use_corotational)
		m_transformation->revertToLastCommit();

    m_eas->alpha = m_eas->alpha_commit;
    m_eas->U = m_eas->U_converged;

    // NDMaterial at gauss Points
    for (int gp = 0; gp < NumGP; ++gp) {
        if (m_material[gp])
            success += m_material[gp]->revertToLastCommit();
    }

    // damping
    for (int gp = 0; gp < NumGP; ++gp) {
        if (m_damping[gp])
            success += m_damping[gp]->revertToLastCommit();
    }

    // done
    return success;
}

int  ASDSolidHex::revertToStart()
{
    int success = 0;

    if(m_use_corotational)
		m_transformation->revertToStart();

    initializePG_EAS();

    // NDMaterial at the integration points
    for (int gp = 0; gp < NumGP; ++gp) {
        if (m_material[gp])
            success += m_material[gp]->revertToStart();
    }

    // damping
    for (int gp = 0; gp < NumGP; ++gp) {
        if (m_damping[gp])
            success += m_damping[gp]->revertToStart();
    }

    return success;
}


int ASDSolidHex::update()
{

    auto& LHS = ASDSolidHexGlobals::instance().LHS;
    auto& RHS = ASDSolidHexGlobals::instance().RHS;
    return calculateAll(LHS, RHS, (OPT_UPDATE));
}

const Matrix& ASDSolidHex::getTangentStiff()
{

    auto& LHS = ASDSolidHexGlobals::instance().LHS;
    auto& RHS = ASDSolidHexGlobals::instance().RHS;
    calculateAll(LHS, RHS, (OPT_LHS));
    return LHS;
}

const Matrix& ASDSolidHex::getInitialStiff()
{

    auto& LHS = ASDSolidHexGlobals::instance().LHS_initial;
    auto& RHS = ASDSolidHexGlobals::instance().RHS;

    // The condensation blocks in m_eas are the ones the NEXT update() uses to
    // advance alpha, and they must correspond to the CURRENT material tangent.
    // This call deliberately uses the INITIAL tangent, so it would otherwise
    // leave initial-tangent Kqu / Kuq / Kqq_inv / alpha_residual behind and the
    // following alpha update would be computed from the wrong operator.
    // It is not hypothetical: getInitialStiff() is called every step by
    // Rayleigh damping with betaK0 or betaKc, and by "algorithm Newton -initial".
    // Harmless while alpha was frozen at zero; live now that alpha moves.
    // Save and restore around the call rather than threading a flag through
    // calculateAll.
    static Vector saved_alpha_residual(12);
    static Matrix saved_Kqq_inv(12, 12);
    static Matrix saved_Kqu(12, 24);
    static Matrix saved_Kuq(24, 12);
    saved_alpha_residual = m_eas->alpha_residual;
    saved_Kqq_inv = m_eas->Kqq_inv;
    saved_Kqu = m_eas->Kqu;
    saved_Kuq = m_eas->Kuq;
    // the enhanced-solve diagnostics belong to the MATERIAL tangent evaluation
    // too: this call would otherwise leave a perfectly healthy elastic rcond and
    // pivot behind, and the 'easState' recorder would report them instead of the
    // softening state that is the whole reason for looking.
    const double saved_rcond = m_eas->rcond;
    const double saved_min_pivot = m_eas->min_pivot;
    const double saved_stab_used = m_eas->stab_used;

    calculateAll(LHS, RHS, (OPT_LHS | OPT_LHS_IS_INITIAL));

    m_eas->alpha_residual = saved_alpha_residual;
    m_eas->Kqq_inv = saved_Kqq_inv;
    m_eas->Kqu = saved_Kqu;
    m_eas->Kuq = saved_Kuq;
    m_eas->rcond = saved_rcond;
    m_eas->min_pivot = saved_min_pivot;
    m_eas->stab_used = saved_stab_used;

    return LHS;
}

const Matrix& ASDSolidHex::getMass()
{
    // Row-sum (lumped) translational mass, same concept as ASDShellQ4::getMass.
    // M_jj = sum_gp N_j * rho * detJ * w  =  integral(N_j * rho dV), and since
    // sum_j N_j = 1 the diagonal sums exactly to rho*V. The Jacobian determinant
    // carries the mesh distortion. Rotational inertia is neglected.
    auto& LHS = ASDSolidHexGlobals::instance().LHS_mass;
    LHS.Zero();

    // nothing to do for a massless element: skip the whole Gauss loop.
    // Element::getDamp() calls getMass() whenever alphaM != 0, and a transient
    // integrator calls it every step, so this early-out is worth having.
    bool has_mass = false;
    for (int i = 0; i < NumGP; i++) {
        if (m_material[i] && m_material[i]->getRho() != 0.0) {
            has_mass = true;
            break;
        }
    }
    if (!has_mass)
        return LHS;

    // a broker-constructed element, or one whose setDomain failed, has no nodes
    for (int a = 0; a < NumNodes; a++) {
        if (nodePtrs[a] == nullptr) {
            opserr << "ASDSolidHex::getMass - element " << this->getTag()
                << " has no node pointers (setDomain not called or failed)\n";
            return LHS;
        }
    }

    // Some matrices
    auto& X = ASDSolidHexGlobals::instance().X;  //nodal coordinates
    auto& dNdh = ASDSolidHexGlobals::instance().dNdh;
    auto& N = ASDSolidHexGlobals::instance().N;

    for (int a = 0; a < NumNodes; a++) {
        const Vector& crd = nodePtrs[a]->getCrds();
        X[a] = vec3(crd(0), crd(1), crd(2));
    }

    // Gauss loop
    for (int i = 0; i < NumGP; i++)
    {
        // Current integration point data
        double xi = XI[i];
        double eta = ETA[i];
        double zeta = ZETA[i];
        double w = WTS[i];
        shapeFunctions(xi, eta, zeta, N);
        dshape(xi, eta, zeta, dNdh);

        // jacobian (Matrix's constructor already zeroes it)
        Matrix J(3, 3);
        for (int a = 0; a < NumNodes; ++a) {
            const double dNa_dxi = dNdh(a, 0);
            const double dNa_deta = dNdh(a, 1);
            const double dNa_dzeta = dNdh(a, 2);
            const vec3& Xa = X[a];

            J(0, 0) += dNa_dxi * Xa.x();
            J(1, 0) += dNa_dxi * Xa.y();
            J(2, 0) += dNa_dxi * Xa.z();

            J(0, 1) += dNa_deta * Xa.x();
            J(1, 1) += dNa_deta * Xa.y();
            J(2, 1) += dNa_deta * Xa.z();

            J(0, 2) += dNa_dzeta * Xa.x();
            J(1, 2) += dNa_dzeta * Xa.y();
            J(2, 2) += dNa_dzeta * Xa.z();
        }

        const double detJ = J(0, 0) * (J(1, 1) * J(2, 2) - J(1, 2) * J(2, 1))
            - J(0, 1) * (J(1, 0) * J(2, 2) - J(1, 2) * J(2, 0))
            + J(0, 2) * (J(1, 0) * J(2, 1) - J(1, 1) * J(2, 0));

        // NOT abs(detJ), deliberately. A reversed connectivity gives detJ < 0 at
        // every gauss point, hence a negative lumped mass -- and that is the
        // consistent answer, because the same sign flip runs through the stiffness
        // and the internal force. With -M and -K the element is a clean negation of
        // itself: statics gives u = -u_correct and a transient gives the same
        // response with the sign of the load reversed. Forcing +M while K stays -K
        // produces a mixture that is neither, and is harder to diagnose. Brick does
        // the same (dvol = wg * xsj, no abs). The reversal is reported once by
        // setDomain, which is where the user can act on it.
        const double dV = w * detJ;

        // Mass density at this integration point
        double rho = m_material[i]->getRho();

        // Add current integration point contribution
        for (int j = 0; j < NumNodes; j++)
        {
            int index = j * 3;  //ux,uy,uz

            // Translational mass contribution
            double Tmass = N(j) * rho * dV;
            for (int q = 0; q < 3; q++)
                LHS(index + q, index + q) += Tmass;

            // Rotational mass neglected
        }
    }

    // Done
    return LHS;
}

void  ASDSolidHex::zeroLoad()
{
    if (m_load)
        m_load->Zero();
}

int
ASDSolidHex::addLoad(ElementalLoad* theLoad, double loadFactor)
{
    // Body force (self weight).
    //
    // The consistent nodal load for a body force b per unit mass is
    //     f_a = integral(N_a * rho * b dV)
    // and integral(N_a * rho dV) is EXACTLY the a-th diagonal entry of the
    // row-sum lumped mass matrix built by getMass(). So f_a = M_aa * b, with no
    // approximation: the same |J| weighting carries the mesh distortion, and the
    // total load is rho*V*b because sum_a N_a = 1.
    //
    // Reusing getMass() also makes this path numerically identical to
    // addInertiaLoadToUnbalance(), so -selfWeight with b = g and a UniformExcitation
    // with ground acceleration -g produce the same nodal forces bit for bit.
    int type = 0;
    const Vector& data = theLoad->getData(type, loadFactor);

    double b[3] = { 0.0, 0.0, 0.0 };
    if (type == LOAD_TAG_BrickSelfWeight) {
        for (int i = 0; i < 3; i++)
            b[i] = loadFactor * m_body[i];
    }
    else if (type == LOAD_TAG_SelfWeight) {
        // compatibility with the -selfWeight class used by the continuum elements
        // (C. McGann, U.W.): data holds a factor per direction
        for (int i = 0; i < 3; i++)
            b[i] = loadFactor * data(i) * m_body[i];
    }
    else {
        opserr << "ASDSolidHex::addLoad - load type " << type
            << " unknown for ele with tag: " << this->getTag() << endln;
        return -1;
    }

    if (b[0] == 0.0 && b[1] == 0.0 && b[2] == 0.0)
        return 0;

    const auto& M = getMass();
    bool has_mass = false;
    for (int i = 0; i < NDOF; i++) {
        if (M(i, i) != 0.0) { has_mass = true; break; }
    }
    if (!has_mass)
        return 0;

    if (m_load == nullptr)
        m_load = new Vector(NDOF);
    auto& F = *m_load;

    for (int a = 0; a < NumNodes; a++) {
        int index = a * 3;
        for (int j = 0; j < 3; j++)
            F(index + j) += M(index + j, index + j) * b[j];
    }

    return 0;
}
int
ASDSolidHex::addInertiaLoadToUnbalance(const Vector& accel)
{
    // Get mass matrix. Note getMass() returns an all-zero matrix when the
    // material density is zero, so bail out before allocating m_load.
    const auto& M = getMass();
    bool has_mass = false;
    for (int i = 0; i < NDOF; i++) {
        if (M(i, i) != 0.0) {
            has_mass = true;
            break;
        }
    }
    if (!has_mass)
        return 0;

    // Allocate load vector if necessary
    if (m_load == nullptr)
        m_load = new Vector(NDOF);
    auto& F = *m_load;

    // Add -M*R*acc to unbalance, taking advantage of the lumped mass matrix.
    // NOTE: loop over the 8 NODES (not the gauss points), and take the nodes
    // from nodePtrs -- the corotational transformation's getNodes() is private
    // and is a null pointer in the linear case anyway.
    for (int i = 0; i < NumNodes; i++)
    {
        const auto& RV = nodePtrs[i]->getRV(accel);
        int index = i * 3;
        for (int j = 0; j < 3; j++)
            F(index + j) -= M(index + j, index + j) * RV(j);
    }

    // Done
    return 0;
}

const Vector& ASDSolidHex::getResistingForce()
{
    // calculate
    auto& LHS = ASDSolidHexGlobals::instance().LHS;
    auto& RHS = ASDSolidHexGlobals::instance().RHS;
    calculateAll(LHS, RHS, (OPT_RHS));
    return RHS;
}


const Vector& ASDSolidHex::getResistingForceIncInertia()
{
    // calculate the static resisting force
    auto& LHS = ASDSolidHexGlobals::instance().LHS;
    auto& RHS = ASDSolidHexGlobals::instance().RHS_winertia;
    calculateAll(LHS, RHS, (OPT_RHS));

    // Add damping terms.
    // NOTE: Element::getDamp() (not overridden here) already assembles
    // alphaM*M + betaK*K into its own buffer, so leaving this out made the
    // damping matrix and the damping force inconsistent with each other.
    if (alphaM != 0.0 || betaK != 0.0 || betaK0 != 0.0 || betaKc != 0.0)
        RHS.addVector(1.0, getRayleighDampingForces(), 1.0);

    // Compute mass
    const auto& M = getMass();

    // Add M*acc to unbalance, taking advantage of the lumped mass matrix
    for (int i = 0; i < NumNodes; i++)
    {
        const auto& A = nodePtrs[i]->getTrialAccel();
        int index = i * 3;
        for (int j = 0; j < 3; j++)
            RHS(index + j) += M(index + j, index + j) * A(j);
    }

    // Done
    return RHS;
}

// number of doubles needed by the EAS state:
//   alpha(12) + alpha_commit(12) + alpha_residual(12)
// + U(24) + U_converged(24)
// + Kqq_inv(12x12) + Kqu(12x24) + Kuq(24x12)
//
// EASData's rcond / min_pivot / stab_used are deliberately NOT in here. They are
// diagnostics recomputed by every condensation, so there is nothing to restore:
// the first evaluation after recvSelf refills them. Adding them would change this
// constant, and therefore the message size, for no gain.
static const int ASDSolidHex_EAS_DATA_SIZE = 12 + 12 + 12 + 24 + 24 + 144 + 288 + 288;

// element level options that are neither EAS state nor Rayleigh/body data:
//   m_eas_stab(1) + m_eas_penalty(1) + m_eas_auto(1, as a double: see the
//   comment on ASDSolidHex::m_eas_auto for why it is not a flag in the ID)
// Kept as its own named term so the two vectData sizes -- in sendSelf and in
// recvSelf -- cannot drift apart, and so that adding an option here does not
// silently redefine what ASDSolidHex_EAS_DATA_SIZE means. These are doubles
// only: a new INT would instead have to grow the ID in BOTH sendSelf and
// recvSelf, and recvSelf sizes its ID before recvID, so a mismatch there is
// silent garbage rather than an error.
#ifdef ASDHEX_EAS_IMPERF
static const int ASDSolidHex_OPTS_DATA_SIZE = 4;  // + m_imperfection
#else
static const int ASDSolidHex_OPTS_DATA_SIZE = 3;
#endif

int ASDSolidHex::sendSelf(int commitTag, Channel& theChannel)
{
    int res = 0;

    // note: we don't check for dataTag == 0 for Element objects as that is
    // taken care of in a commit by the Domain object
    int dataTag = this->getDbTag();
    int counter;

    // has load flag
    bool has_load = m_load != nullptr;

    // INT data
    // 1 tag + 8 node tags + 1 corotational flag + 1 initialization flag
    // + 1 has_load flag + 16 -> 8 pairs of (material class tag, material db tag)
    // + 2 -> damping class tag + damping db tag (0, 0 when there is no damping)
    // + 1 activation flag
    static ID idData(31);
    counter = 0;
    idData(counter++) = this->getTag();
    for (int i = 0; i < NumNodes; ++i)
        idData(counter++) = m_node_ids(i);
    idData(counter++) = static_cast<int>(m_use_corotational);
    idData(counter++) = static_cast<int>(m_initialized);
    idData(counter++) = static_cast<int>(has_load);
    for (int i = 0; i < NumGP; i++) {
        idData(counter++) = m_material[i]->getClassTag();
        int matDbTag = m_material[i]->getDbTag();
        // NOTE: we do have to ensure that the material has a database tag if we
        // are sending to a database channel.
        if (matDbTag == 0) {
            matDbTag = theChannel.getDbTag();
            if (matDbTag != 0)
                m_material[i]->setDbTag(matDbTag);
        }
        idData(counter++) = matDbTag;
    }
    // damping: one shared db tag for all 8 copies, as ASDShellQ4 does
    if (m_damping[0]) {
        idData(counter++) = m_damping[0]->getClassTag();
        int dmpDbTag = m_damping[0]->getDbTag();
        if (dmpDbTag == 0) {
            dmpDbTag = theChannel.getDbTag();
            if (dmpDbTag != 0) {
                for (int i = 0; i < NumGP; i++)
                    m_damping[i]->setDbTag(dmpDbTag);
            }
        }
        idData(counter++) = dmpDbTag;
    }
    else {
        idData(counter++) = 0;
        idData(counter++) = 0;
    }
    // activation state: an element deactivated before the transfer must come back deactivated
    idData(counter++) = is_this_element_active ? 1 : 0;

    res = theChannel.sendID(dataTag, commitTag, idData);
    if (res < 0) {
        opserr << "WARNING ASDSolidHex::sendSelf() - " << this->getTag() << " failed to send ID\n";
        return res;
    }

    // DOUBLE data
    // 4 rayleigh damping factors
    // + 3 body force
    // + 24 initial displacement offset (see m_U0), so a restore does not
    //   re-capture it from the already displaced nodes
    // + EAS state
    // + (optional) 24 load
    // + (optional) transformation internal data
    int NLoad = has_load ? NDOF : 0;
    int NT = m_use_corotational ? m_transformation->internalDataSize() : 0;
    Vector vectData(4 + 3 + NDOF + ASDSolidHex_EAS_DATA_SIZE + ASDSolidHex_OPTS_DATA_SIZE + NLoad + NT);
    counter = 0;
    vectData(counter++) = alphaM;
    vectData(counter++) = betaK;
    vectData(counter++) = betaK0;
    vectData(counter++) = betaKc;
    for (int i = 0; i < 3; ++i) vectData(counter++) = m_body[i];
    vectData(counter++) = m_eas_stab;
    vectData(counter++) = m_eas_penalty;
    vectData(counter++) = m_eas_auto ? 1.0 : 0.0;
#ifdef ASDHEX_EAS_IMPERF
    // the 8 factors are NOT sent: setDomain recomputes them from the tag, which
    // recvSelf has already installed by then
    vectData(counter++) = m_imperfection;
#endif
    for (int i = 0; i < NDOF; ++i) vectData(counter++) = m_U0(i);
    for (int i = 0; i < 12; ++i) vectData(counter++) = m_eas->alpha(i);
    for (int i = 0; i < 12; ++i) vectData(counter++) = m_eas->alpha_commit(i);
    for (int i = 0; i < 12; ++i) vectData(counter++) = m_eas->alpha_residual(i);
    for (int i = 0; i < NDOF; ++i) vectData(counter++) = m_eas->U(i);
    for (int i = 0; i < NDOF; ++i) vectData(counter++) = m_eas->U_converged(i);
    for (int i = 0; i < 12; ++i)
        for (int j = 0; j < 12; ++j) vectData(counter++) = m_eas->Kqq_inv(i, j);
    for (int i = 0; i < 12; ++i)
        for (int j = 0; j < NDOF; ++j) vectData(counter++) = m_eas->Kqu(i, j);
    for (int i = 0; i < NDOF; ++i)
        for (int j = 0; j < 12; ++j) vectData(counter++) = m_eas->Kuq(i, j);
    if (has_load) {
        for (int i = 0; i < NDOF; ++i) vectData(counter++) = (*m_load)(i);
    }
    if (m_use_corotational)
        m_transformation->saveInternalData(vectData, counter);

    res = theChannel.sendVector(dataTag, commitTag, vectData);
    if (res < 0) {
        opserr << "WARNING ASDSolidHex::sendSelf() - " << this->getTag() << " failed to send Vector\n";
        return res;
    }

    // send all materials
    for (int i = 0; i < NumGP; i++) {
        res = m_material[i]->sendSelf(commitTag, theChannel);
        if (res < 0) {
            opserr << "WARNING ASDSolidHex::sendSelf() - " << this->getTag() << " failed to send its NDMaterial\n";
            return res;
        }
    }

    // send the damping objects
    if (m_damping[0]) {
        for (int i = 0; i < NumGP; i++) {
            res = m_damping[i]->sendSelf(commitTag, theChannel);
            if (res < 0) {
                opserr << "WARNING ASDSolidHex::sendSelf() - " << this->getTag()
                    << " failed to send its Damping\n";
                return res;
            }
        }
    }

    // done
    return res;
}

int  ASDSolidHex::recvSelf(int commitTag, Channel& theChannel, FEM_ObjectBroker& theBroker)
{
    int res = 0;

    int dataTag = this->getDbTag();
    int counter;

    // INT data
    static ID idData(31);
    res = theChannel.recvID(dataTag, commitTag, idData);
    if (res < 0) {
        opserr << "WARNING ASDSolidHex::recvSelf() - failed to receive ID\n";
        return res;
    }

    counter = 0;
    this->setTag(idData(counter++));
    for (int i = 0; i < NumNodes; ++i)
        m_node_ids(i) = idData(counter++);
    m_use_corotational = static_cast<bool>(idData(counter++));
    m_initialized = static_cast<bool>(idData(counter++));
    bool has_load = static_cast<bool>(idData(counter++));

    // Allocate the transformation now that we know the kinematics.
    // NOTE: do NOT delete and recreate an existing one. restoreInternalData()
    // restores m_U0/m_Q0/m_C0 but NOT the node pointers, and the database
    // restore path calls recvSelf WITHOUT a following setDomain -- a fresh
    // transformation would then be left with null nodes and segfault on the next
    // computeGlobalDisplacements(). Reusing it keeps the pointers setDomain
    // already installed, and in the object-broker path there is nothing to reuse
    // so a new one is built and setDomain fills it in later.
    if (m_use_corotational) {
        if (m_transformation == nullptr)
            m_transformation = new ASDSolidHexCorotationalTransformation();
    }
    else if (m_transformation) {
        delete m_transformation;
        m_transformation = nullptr;
    }

    // allocate the load vector if the sender had one
    if (has_load) {
        if (m_load == nullptr)
            m_load = new Vector(NDOF);
    }
    else {
        if (m_load) {
            delete m_load;
            m_load = nullptr;
        }
    }

    // create the materials, re-using the existing ones when the class tag matches
    for (int i = 0; i < NumGP; i++) {
        int matClassTag = idData(counter++);
        int matDbTag = idData(counter++);
        if (m_material[i] == nullptr || m_material[i]->getClassTag() != matClassTag) {
            if (m_material[i])
                delete m_material[i];
            m_material[i] = theBroker.getNewNDMaterial(matClassTag);
            if (m_material[i] == nullptr) {
                opserr << "ASDSolidHex::recvSelf() - failed to get a new NDMaterial of class tag "
                    << matClassTag << endln;
                return -1;
            }
        }
        m_material[i]->setDbTag(matDbTag);
    }

    // damping tags. The objects themselves are rebuilt and received at the end,
    // after the materials, in the same order sendSelf wrote them.
    int dmpClassTag = idData(counter++);
    int dmpDbTag = idData(counter++);

    // activation state: an element deactivated before the transfer must come back deactivated
    is_this_element_active = idData(counter++) == 1;

    // DOUBLE data
    int NLoad = has_load ? NDOF : 0;
    int NT = m_use_corotational ? m_transformation->internalDataSize() : 0;
    Vector vectData(4 + 3 + NDOF + ASDSolidHex_EAS_DATA_SIZE + ASDSolidHex_OPTS_DATA_SIZE + NLoad + NT);
    res = theChannel.recvVector(dataTag, commitTag, vectData);
    if (res < 0) {
        opserr << "WARNING ASDSolidHex::recvSelf() - failed to receive Vector\n";
        return res;
    }

    counter = 0;
    alphaM = vectData(counter++);
    betaK = vectData(counter++);
    betaK0 = vectData(counter++);
    betaKc = vectData(counter++);
    for (int i = 0; i < 3; ++i) m_body[i] = vectData(counter++);
    m_eas_stab = vectData(counter++);
    m_eas_penalty = vectData(counter++);
    m_eas_auto = (vectData(counter++) != 0.0);
#ifdef ASDHEX_EAS_IMPERF
    m_imperfection = vectData(counter++);
#endif
    for (int i = 0; i < NDOF; ++i) m_U0(i) = vectData(counter++);
    if (m_eas == nullptr)
        m_eas = new EASData();
    for (int i = 0; i < 12; ++i) m_eas->alpha(i) = vectData(counter++);
    for (int i = 0; i < 12; ++i) m_eas->alpha_commit(i) = vectData(counter++);
    for (int i = 0; i < 12; ++i) m_eas->alpha_residual(i) = vectData(counter++);
    for (int i = 0; i < NDOF; ++i) m_eas->U(i) = vectData(counter++);
    for (int i = 0; i < NDOF; ++i) m_eas->U_converged(i) = vectData(counter++);
    for (int i = 0; i < 12; ++i)
        for (int j = 0; j < 12; ++j) m_eas->Kqq_inv(i, j) = vectData(counter++);
    for (int i = 0; i < 12; ++i)
        for (int j = 0; j < NDOF; ++j) m_eas->Kqu(i, j) = vectData(counter++);
    for (int i = 0; i < NDOF; ++i)
        for (int j = 0; j < 12; ++j) m_eas->Kuq(i, j) = vectData(counter++);
    if (has_load) {
        for (int i = 0; i < NDOF; ++i) (*m_load)(i) = vectData(counter++);
    }
    if (m_use_corotational)
        m_transformation->restoreInternalData(vectData, counter);

    // receive all materials
    for (int i = 0; i < NumGP; i++) {
        res = m_material[i]->recvSelf(commitTag, theChannel, theBroker);
        if (res < 0) {
            opserr << "WARNING ASDSolidHex::recvSelf() - failed to receive its NDMaterial\n";
            return res;
        }
    }

    // damping: rebuild only when absent or of the wrong class, then receive.
    // dmpClassTag == 0 means the sender had no damping, so drop ours if any.
    if (dmpClassTag == 0) {
        for (int i = 0; i < NumGP; i++) {
            if (m_damping[i]) {
                delete m_damping[i];
                m_damping[i] = nullptr;
            }
        }
    }
    else {
        for (int i = 0; i < NumGP; i++) {
            if (m_damping[i] == nullptr || m_damping[i]->getClassTag() != dmpClassTag) {
                if (m_damping[i])
                    delete m_damping[i];
                m_damping[i] = theBroker.getNewDamping(dmpClassTag);
                if (m_damping[i] == nullptr) {
                    opserr << "ASDSolidHex::recvSelf() - failed to get a new Damping of class tag "
                        << dmpClassTag << endln;
                    return -1;
                }
            }
            m_damping[i]->setDbTag(dmpDbTag);
            res = m_damping[i]->recvSelf(commitTag, theChannel, theBroker);
            if (res < 0) {
                opserr << "WARNING ASDSolidHex::recvSelf() - failed to receive its Damping\n";
                return res;
            }
        }
    }

    // done
    return res;
}


Response*
ASDSolidHex::setResponse(const char** argv, int argc, OPS_Stream& output)
{
    Response* theResponse = 0;

    char outputData[32];

    output.tag("ElementOutput");
    output.attr("eleType", "ASDSolidHex");
    output.attr("eleTag", this->getTag());

    // use the stored node ids, not nodePtrs: setResponse can be reached before
    // setDomain has installed the pointers (or after it failed)
    for (int i = 1; i <= NumNodes; i++) {
        sprintf(outputData, "node%d", i);
        output.attr(outputData, m_node_ids(i - 1));
    }

    if (argc < 1) {
        output.endTag();
        return theResponse;
    }


    if (strcmp(argv[0], "force") == 0 || strcmp(argv[0], "forces") == 0) {

        for (int i = 1; i <= 8; i++) {
            sprintf(outputData, "P1_%d", i);
            output.tag("ResponseType", outputData);
            sprintf(outputData, "P2_%d", i);
            output.tag("ResponseType", outputData);
            sprintf(outputData, "P3_%d", i);
            output.tag("ResponseType", outputData);
        }

        theResponse = new ElementResponse(this, 1, this->getResistingForce());

    }
    else if (strcmp(argv[0], "material") == 0 || strcmp(argv[0], "Material") == 0 ||
             strcmp(argv[0], "integrPoint") == 0) {

        // argv[1] is the 1-based gauss point, argv[2..] the material's own request
        if (argc < 3) {
            opserr << "ASDSolidHex::setResponse() - need to specify more data\n";
            output.endTag();
            return 0;
        }
        int pointNum = atoi(argv[1]);
        if (pointNum > 0 && pointNum <= NumGP) {
            // pointNum is in REPORTING order; translate to the internal one
            const int ig = GP_REPORT_TO_INTERNAL[pointNum - 1];
            if (m_material[ig]) {
                output.tag("GaussPoint");
                output.attr("number", pointNum);
                output.attr("xi", XI[ig]);
                output.attr("eta", ETA[ig]);
                output.attr("zeta", ZETA[ig]);

                theResponse = m_material[ig]->setResponse(&argv[2], argc - 2, output);

                output.endTag(); // GaussPoint
            }
        }


    }
    else if (strcmp(argv[0], "stresses") == 0) {

        for (int i = 0; i < NumGP; i++) {
            const int ig = GP_REPORT_TO_INTERNAL[i];
            output.tag("GaussPoint");
            output.attr("number", i + 1);
            output.attr("xi", XI[ig]);
            output.attr("eta", ETA[ig]);
            output.attr("zeta", ZETA[ig]);
            output.tag("NdMaterialOutput");
            output.attr("classType", m_material[ig]->getClassTag());
            output.attr("tag", m_material[ig]->getTag());

            output.tag("ResponseType", "sigma11");
            output.tag("ResponseType", "sigma22");
            output.tag("ResponseType", "sigma33");
            output.tag("ResponseType", "sigma12");
            output.tag("ResponseType", "sigma23");
            output.tag("ResponseType", "sigma13");

            output.endTag(); // NdMaterialOutput
            output.endTag(); // GaussPoint
        }
        theResponse = new ElementResponse(this, 3, Vector(48));

    }
    else if (strcmp(argv[0], "strains") == 0) {

        for (int i = 0; i < NumGP; i++) {
            const int ig = GP_REPORT_TO_INTERNAL[i];
            output.tag("GaussPoint");
            output.attr("number", i + 1);
            output.attr("xi", XI[ig]);
            output.attr("eta", ETA[ig]);
            output.attr("zeta", ZETA[ig]);
            output.tag("NdMaterialOutput");
            output.attr("classType", m_material[ig]->getClassTag());
            output.attr("tag", m_material[ig]->getTag());

            output.tag("ResponseType", "eps11");
            output.tag("ResponseType", "eps22");
            output.tag("ResponseType", "eps33");
            output.tag("ResponseType", "eps12");
            output.tag("ResponseType", "eps23");
            output.tag("ResponseType", "eps13");

            output.endTag(); // NdMaterialOutput
            output.endTag(); // GaussPoint
        }
        theResponse = new ElementResponse(this, 4, Vector(48));

    }
    else if (strcmp(argv[0], "dampingStresses") == 0 && m_damping[0]) {
        // only advertised when the element actually has Damping objects: before,
        // this branch handed out responseID 5 unconditionally while getResponse
        // had no case 5, so the recorder silently produced nothing.
        for (int i = 0; i < NumGP; i++) {
            const int ig = GP_REPORT_TO_INTERNAL[i];
            output.tag("GaussPoint");
            output.attr("number", i + 1);
            output.attr("xi", XI[ig]);
            output.attr("eta", ETA[ig]);
            output.attr("zeta", ZETA[ig]);
            output.tag("NdMaterialOutput");
            output.attr("classType", m_material[ig]->getClassTag());
            output.attr("tag", m_material[ig]->getTag());

            output.tag("ResponseType", "sigma11");
            output.tag("ResponseType", "sigma22");
            output.tag("ResponseType", "sigma33");
            output.tag("ResponseType", "sigma12");
            output.tag("ResponseType", "sigma23");
            output.tag("ResponseType", "sigma13");

            output.endTag(); // NdMaterialOutput
            output.endTag(); // GaussPoint
        }
        theResponse = new ElementResponse(this, 5, Vector(48));

    }
    else if (strcmp(argv[0], "stress3D6") == 0) {
        // this response if for vtkhdf recorder which averages the stresses
        // over the 8 gauss points and returns the 6 values
        output.tag("GaussPoint");
        output.attr("number", 1);
        output.tag("NdMaterialOutput");
        output.attr("classType", m_material[0]->getClassTag());
        output.attr("tag", m_material[0]->getTag());

        output.tag("ResponseType", "sigma11");
        output.tag("ResponseType", "sigma22");
        output.tag("ResponseType", "sigma33");
        output.tag("ResponseType", "sigma12");
        output.tag("ResponseType", "sigma23");
        output.tag("ResponseType", "sigma13");

        output.endTag(); // NdMaterialOutput
        output.endTag(); // GaussPoint
        theResponse = new ElementResponse(this, 6, Vector(6));

    }
    else if (strcmp(argv[0], "strain3D6") == 0) {
        // this response if for vtkhdf recorder which averages the strains
        // over the 8 gauss points and returns the 6 values
        output.tag("GaussPoint");
        output.attr("number", 1);
        output.tag("NdMaterialOutput");
        output.attr("classType", m_material[0]->getClassTag());
        output.attr("tag", m_material[0]->getTag());

        output.tag("ResponseType", "eps11");
        output.tag("ResponseType", "eps22");
        output.tag("ResponseType", "eps33");
        output.tag("ResponseType", "eps12");
        output.tag("ResponseType", "eps23");
        output.tag("ResponseType", "eps13");

        output.endTag(); // NdMaterialOutput
        output.endTag(); // GaussPoint

        theResponse = new ElementResponse(this, 7, Vector(6));
    }
    // ------------------------------------------------------------------
    // PG-EAS diagnostics. The enhanced solve is entirely internal to the
    // element -- 12 parameters condensed out, one lagged step per global
    // iteration, no residual check -- so without these there is no way to see
    // what it is doing, and a softening analysis that misbehaves gives no
    // evidence of WHY.
    //
    // No GP_REPORT_TO_INTERNAL permutation applies here: alpha and its residual
    // are element quantities, not per-gauss-point ones.
    // ------------------------------------------------------------------
    else if (strcmp(argv[0], "easAlpha") == 0) {
        for (int i = 0; i < 12; ++i)
            output.tag("ResponseType", "alpha");
        theResponse = new ElementResponse(this, 8, Vector(12));
    }
    else if (strcmp(argv[0], "easResidual") == 0) {
        // h = -int(G_test^T sigma dV). NOTE: this has units of STRESS, not force:
        // G_test carries a 1/jdet and dV carries a jdet, and they cancel (see the
        // C3 weight note in metric_basis::orthogonalize). So h is NOT work
        // conjugate to alpha, and no energy norm applies to it.
        for (int i = 0; i < 12; ++i)
            output.tag("ResponseType", "h");
        theResponse = new ElementResponse(this, 9, Vector(12));
    }
    else if (strcmp(argv[0], "easState") == 0) {
        output.tag("ResponseType", "hNorm2");
        output.tag("ResponseType", "hNormInf");
        output.tag("ResponseType", "rcond");
        output.tag("ResponseType", "minPivot");
        output.tag("ResponseType", "stab");
        output.tag("ResponseType", "penalty");
        theResponse = new ElementResponse(this, 10, Vector(6));
    }
    output.endTag(); // ElementOutput
    return theResponse;
}

int
ASDSolidHex::getResponse(int responseID, Information& eleInfo)
{
    static Vector stresses(48);

    if (responseID == 1)
        return eleInfo.setVector(this->getResistingForce());

    else if (responseID == 3) {

        // Loop over the integration points in REPORTING order
        int cnt = 0;
        for (int i = 0; i < NumGP; i++) {

            const Vector& sigma = m_material[GP_REPORT_TO_INTERNAL[i]]->getStress();
            stresses(cnt++) = sigma(0);
            stresses(cnt++) = sigma(1);
            stresses(cnt++) = sigma(2);
            stresses(cnt++) = sigma(3);
            stresses(cnt++) = sigma(4);
            stresses(cnt++) = sigma(5);
        }
        return eleInfo.setVector(stresses);

    }
    else if (responseID == 4) {

        // Loop over the integration points in REPORTING order
        int cnt = 0;
        for (int i = 0; i < NumGP; i++) {

            const Vector& sigma = m_material[GP_REPORT_TO_INTERNAL[i]]->getStrain();
            stresses(cnt++) = sigma(0);
            stresses(cnt++) = sigma(1);
            stresses(cnt++) = sigma(2);
            stresses(cnt++) = sigma(3);
            stresses(cnt++) = sigma(4);
            stresses(cnt++) = sigma(5);
        }
        return eleInfo.setVector(stresses);
    }
    else if (responseID == 5) {

        // damping stresses, in REPORTING order. setResponse only hands out this
        // id when m_damping[0] is non null, so no null check is needed here.
        int cnt = 0;
        for (int i = 0; i < NumGP; i++) {

            const Vector& sigma = m_damping[GP_REPORT_TO_INTERNAL[i]]->getDampingForce();
            stresses(cnt++) = sigma(0);
            stresses(cnt++) = sigma(1);
            stresses(cnt++) = sigma(2);
            stresses(cnt++) = sigma(3);
            stresses(cnt++) = sigma(4);
            stresses(cnt++) = sigma(5);
        }
        return eleInfo.setVector(stresses);
    }
    else if (responseID == 6) {

        // Loop over the integration points
        Vector tmpStress(6);
        for (int i = 0; i < 8; i++) {

            // Get material stress
            const Vector& sigma = m_material[i]->getStress();
            tmpStress(0) += sigma(0) * 0.125;
            tmpStress(1) += sigma(1) * 0.125;
            tmpStress(2) += sigma(2) * 0.125;
            tmpStress(3) += sigma(3) * 0.125;
            tmpStress(4) += sigma(4) * 0.125;
            tmpStress(5) += sigma(5) * 0.125;
        }

        return eleInfo.setVector(tmpStress);
    }
    else if (responseID == 7) {

        // Loop over the integration points
        int cnt = 0;
        Vector tmpStrain(6);
        for (int i = 0; i < 8; i++) {

            // Get material strain
            const Vector& sigma = m_material[i]->getStrain();
            tmpStrain(0) += sigma(0);
            tmpStrain(1) += sigma(1);
            tmpStrain(2) += sigma(2);
            tmpStrain(3) += sigma(3);
            tmpStrain(4) += sigma(4);
            tmpStrain(5) += sigma(5);
        }
        tmpStrain /= 8.0;

        return eleInfo.setVector(tmpStrain);
    }
    else if (responseID == 8) {
        return eleInfo.setVector(m_eas->alpha);
    }
    else if (responseID == 9) {
        return eleInfo.setVector(m_eas->alpha_residual);
    }
    else if (responseID == 10) {

        // rcond, min_pivot and stab_used are snapshots taken during the last
        // condensation: they are read from EASData, never from the globals
        // singleton, which the next element to be evaluated overwrites.
        // NOTE: with -easPenalty on, hNorm2/hNormInf are the norms of the
        // PENALIZED residual h - p*Kqq0_sym*alpha, which is the equation actually
        // being solved. That is the right thing to watch converge; it is not the
        // unpenalized int(G_test^T sigma dV), which no longer goes to zero.
        static Vector eas_state(6);
        const Vector& h = m_eas->alpha_residual;
        double h_inf = 0.0;
        for (int i = 0; i < h.Size(); ++i) {
            double a = std::abs(h(i));
            if (a > h_inf)
                h_inf = a;
        }
        eas_state(0) = h.Norm();
        eas_state(1) = h_inf;
        eas_state(2) = m_eas->rcond;
        eas_state(3) = m_eas->min_pivot;
        eas_state(4) = m_eas->stab_used;
        eas_state(5) = m_eas_penalty;
        return eleInfo.setVector(eas_state);
    }
    else
        return -1;
}

int ASDSolidHex::setParameter(const char** argv, int argc, Parameter& param)
{
    int res = -1;

    int matRes = res;

    for (int i = 0; i < NumGP; i++)
    {
        if (!m_material[i]) continue;
        matRes = m_material[i]->setParameter(argv, argc, param);
        if (matRes != -1)
            res = matRes;
    }
    return res;
}

double ASDSolidHex::getCharacteristicLength(void)
{
    // The base class returns the minimum distance between two nodes, which for a
    // well-shaped hexahedron is the shortest edge.
    //
    // It is then HALVED, exactly as ASDShellQ4 does when its EAS is active: the
    // enhanced modes let the element localise into a single layer of gauss points
    // rather than across its whole width, so the crack band a regularised
    // material should smear over is half the element size. Unlike the shell,
    // where the EAS is optional (ASDShellQ4::m_eas is a pointer and the
    // halving is conditional on it), H1U/E12 is EAS unconditionally -- there are
    // always 12 enhanced parameters and no option to switch them off -- so the
    // halving is unconditional too.
    //
    // This matters only for materials that regularise on lch (ASDConcrete3D and
    // relatives): with lch twice too large they dissipate twice the fracture
    // energy. It is inert for elastic and for non-regularised materials.
    double lch = Element::getCharacteristicLength();
    lch /= 2.0;
    return lch;
}

int ASDSolidHex::calculateAll(Matrix& LHS, Vector& RHS, int options)
{
    // Check options
    // The tangent LHS always requires the RHS:
    //  - corotationally, because the geometric stiffness is built from the
    //    projected internal forces;
    //  - in both cases, because alpha_residual (h) is accumulated in the RHS
    //    branch and is an input to the PG-EAS static condensation and to
    //    updatePG_EAS. Computing the tangent without it left h = 0 on
    //    tangent-only calls in the linear path.
    if (options & OPT_LHS) {
        options |= OPT_RHS;
    }
    // =========================================================
    // Preliminary operations
    // =========================================================
    int result = 0;
    if (options & OPT_RHS)
        RHS.Zero();
    if (options & OPT_LHS)
        LHS.Zero();

    // =========================================================
    // Compute the global displacement vector of nodes
    // =========================================================
    auto& UG = ASDSolidHexGlobals::instance().UG;
    auto& UL = ASDSolidHexGlobals::instance().UL;
    UG.Zero();
    UL.Zero();

    ASDSolidHexLocalCoordinateSystem local_cs;
    for (int a = 0; a < ASDSolidHexGlobals::nnode; a++) {
        Node* nodeA = nodePtrs[a];
        const Vector& ua = nodeA->getTrialDisp();  // (ux, uy, uz)
        const int index = 3 * a;
        UG(index) = ua(0);  // ux_a
        UG(index + 1) = ua(1);  // uy_a
        UG(index + 2) = ua(2);  // uz_a
        const Vector& crd = nodeA->getCrds();  // (x,y,z)
        ASDSolidHexGlobals::instance().X[a] = vec3(crd(0), crd(1), crd(2));
    }

    if (!m_use_corotational) {
		// if is linear use that displacmenets as local displacements as well,
        // net of the initial displacement offset (see m_U0). In the corotational
        // path the offset is applied by the transformation itself.
        UL = UG;
        UL.addVector(1.0, m_U0, -1.0);
    } else {
		// if is corotaional compute the deformational part of displacements as local displacements
		// Global Displacements        
		m_transformation->computeGlobalDisplacements(UG);

        // global displacement

        if (options & OPT_UPDATE)
            m_transformation->update(UG);


        // compute the local coordinate system
	    local_cs = m_transformation->createLocalCoordinateSystem(UG);

        // local displacements
        m_transformation->calculateLocalDisplacements(local_cs,UG, UL);

    }

    // =========================================================
    // Initialize structs for the current iteration
    // =========================================================
    auto& metric_basis = ASDSolidHexGlobals::instance().mb;
    {
        // The metric basis is a function of the reference geometry only, so it is
        // built on the first call and copied back in on every later one. See the
        // comment on ASDSolidHexRefMetric.
        if (m_ref == nullptr)
            m_ref = new ASDSolidHexRefMetric();
        if (m_ref->valid) {
            asdhex_ref_load(*m_ref, metric_basis);
        }
        else {
            if (!metric_basis.initialize_metric(ASDSolidHexGlobals::instance().X)) {
                opserr << "ASDSolidHex::calculateAll - element " << this->getTag()
                    << " has a singular geometry, see above\n";
                return -1;
            }
            asdhex_ref_save(metric_basis, *m_ref);
        }
    }

    // ==========================================================
    // Before the Gauss loop intialize the GP-EAS paramreters
    // ==========================================================
    // initialize the EAS internal parameters (alpha)

    // 3) pgeas (frame skew, V_inv, R^{-T})
    // m_pgeas.compute(ASDSolidHexGlobals::instance().X, dummy, XI, ETA, ZETA, WTS);
    // =========================================================
    // Initialize some variables from global storage
    // =========================================================
    auto& B_test = ASDSolidHexGlobals::instance().B_test;
    auto& B_trial = ASDSolidHexGlobals::instance().B_trial;
    auto& G_trial = ASDSolidHexGlobals::instance().G_trial;
    auto& G_test = ASDSolidHexGlobals::instance().G_test;
    // For PG-EAS: initialize the strain stress and material matrix
    auto& eps = ASDSolidHexGlobals::instance().eps;
    auto& sig = ASDSolidHexGlobals::instance().sig;
    auto& C = ASDSolidHexGlobals::instance().C;
    // For PG-EAS: initialize the 4 stiffness blocks
    auto& k_uu = ASDSolidHexGlobals::instance().k_uu;
    auto& k_qq = ASDSolidHexGlobals::instance().k_qq;
    auto& k_qq0 = ASDSolidHexGlobals::instance().k_qq0;
    auto& k_qq0_sym = ASDSolidHexGlobals::instance().k_qq0_sym;

    // A-op: worst material-tangent definiteness over the gauss points, and the
    // regularization factor that worst point requires. See the EAS_MAT_PD_TOL
    // comment: the decision is pointwise on C, never on the assembled block.
    double mat_min_pivot = 0.0;      // min over GP of pivot(C) / meanDiag(C0)
    double stab_needed = 0.0;
    bool mat_pivot_seen = false;

    // =========================================================
    // Update the PG-EAS enhanced parameters alpha
    // =========================================================
    // This MUST happen here: before the gauss loop, so that the loop feeds the
    // material the strain built with the CURRENT alpha, and before the blocks
    // below are zeroed, so that it still sees Kqu / Kqq_inv / alpha_residual
    // from the previous force or stiffness evaluation.
    //
    // It used to run at the very END of calculateAll. That was not merely a
    // lagged update, it disabled the enhanced modes completely: update() calls
    // calculateAll with OPT_UPDATE alone, the (OPT_RHS || OPT_LHS) block that
    // fills Kqu / Kqq_inv / alpha_residual is then skipped, so all three inputs
    // had just been zeroed a few lines below and the increment
    //   d_alpha = -Kqq^-1 (Kqu dU - h)
    // evaluated to exactly zero on every call. alpha stayed at its initial
    // value of zero for the whole analysis.
    //
    // Displacements were still right because the static condensation of the
    // residual and tangent is algebraically exact for a linear material, but
    // the material never saw the enhanced strain: gauss point stress and strain
    // came out equal to a plain displacement element's, bit for bit (verified
    // against stdBrick), which contradicts the pointwise exactness the
    // formulation is built for, and integrates any nonlinear material on the
    // wrong strain. Same order as ASDShellQ4::calculateAll -> AGQIupdate.
    if (options & OPT_UPDATE)
        updatePG_EAS(UL);

    // Zero the stiffness blocks -- ONLY when this call is going to refill them.
    // They are members of m_eas and are inputs to updatePG_EAS on the next
    // update() call, which never refills them (it passes OPT_UPDATE alone).
    // Zeroing them unconditionally therefore destroyed the condensation data
    // between the evaluation that produced it and the update that consumes it.
    // Kqq_inv is not zeroed at all: it is assigned wholesale after the loop.
    if ((options & OPT_RHS) || (options & OPT_LHS)) {
        k_uu.Zero();
        k_qq.Zero();
        k_qq0.Zero();
        k_qq0_sym.Zero();
        m_eas->Kqu.Zero();
        m_eas->Kuq.Zero();
        m_eas->alpha_residual.Zero();
    }

    // Gauss loop 2x2x2
    for (int igauss = 0; igauss < NumGP; igauss++)
    {
        // We're assuming that EAS is always considere and setted to true
        // Current integration point data
        double xi = XI[igauss];
        double eta = ETA[igauss];
        double zeta = ZETA[igauss];
        double w = WTS[igauss];
        double dv = 0.;

        // ==================================================
        // Compute the Jacobian and the elementary volume
        // ==================================================
        Vector N(NumGP);
        Matrix dNdh(NumGP, 3);
        shapeFunctions(xi, eta, zeta, N);
        dshape(xi, eta, zeta, dNdh);
        Jacobian3d J;
        if (!J.calculate(ASDSolidHexGlobals::instance().X, dNdh)) {
            opserr << "ASDSolidHex::calculateAll - element " << this->getTag()
                << " , gauss point " << igauss << " , see above\n";
            return -1;
        }

        dv = w * J.detJ;


        // ===================================================
        // Compute B_test in physical coordinates
        // ===================================================
        compute_B_test(dNdh, J.invJ, B_test);

        // ===================================================
        // Compute the shape functions matrices -> B_test, B_trial, G_test, G_trial
        // ===================================================

        metric_basis.computeMetricShape(xi, eta, zeta, J.detJ);
        B_trial = metric_basis.B_trial;
        G_trial = metric_basis.G_trial;

        // Recall G_test
        G_test = metric_basis.G_at_gp[igauss];


        // ===================================================
        // Update the Strain and send them to the material
        // ===================================================
        if (options & OPT_UPDATE)
        {
            // compatible strain at the current GP
            // eps = B_v * u

            eps.addMatrixVector(0.0, B_trial, UL, 1.0);
            // add incompatble strain from EAS at the current GP
            eps.addMatrixVector(1.0, G_trial, m_eas->alpha, 1.0);
            // never feed the material a non-finite strain: a trial state built
            // on inf/nan can poison its internal variables even if this
            // iteration is later discarded. Fail the update instead, so the
            // algorithm cuts the step while the committed state is still good.
            if (!isFinite(eps)) {
                opserr << "ASDSolidHex::calculateAll() - element " << this->getTag()
                    << ", gauss point " << igauss
                    << ": non-finite trial strain\n";
                return -1;
            }
            // set the trial strain to the material allocated to the Gauss point
            result += m_material[igauss]->setTrialStrain(eps);
        }


        // Integrate RHS
        if (options & OPT_RHS)
        {
            // add the contribution of the current GP to the RHS: f_int = ∫ B_v^T sigma dV
            sig = m_material[igauss]->getStress();
#ifdef ASDHEX_EAS_IMPERF
            // EXPERIMENT: scale sigma BEFORE the damping force is added -- that
            // force is a separate physical contribution and is not part of the
            // material's constitutive response. The matching tangent scaling is
            // applied below, so both the residual and its linearization see the
            // same perturbed material.
            sig *= m_imperfection_f[igauss];
#endif
            // add the damping stress, same pattern as ASDShellQ4: the Damping
            // object is driven by the total stress and returns the additional
            // (rate dependent) stress to be integrated with it
            if (m_damping[igauss]) {
                m_damping[igauss]->update(sig);
                sig += m_damping[igauss]->getDampingForce();
            }
            // the stress goes into the RHS as it is: no tangent fallback can
            // mask a non-finite value here (a division by zero inside the
            // material, typically). Fail the evaluation instead of handing
            // inf/nan to the solver, so the algorithm cuts the step.
            if (!isFinite(sig)) {
                opserr << "ASDSolidHex::calculateAll() - element " << this->getTag()
                    << ", gauss point " << igauss << ", material " << m_material[igauss]->getTag()
                    << ": non-finite stress\n";
                return -1;
            }
            // add internal force contribution to the global RHS
            RHS.addMatrixTransposeVector(1.0, B_test, sig, dv);
            // upload the alpha residual for the current GP
            m_eas->alpha_residual.addMatrixTransposeVector(1.0, G_test, sig, -dv);

        }

        // Compute the sub-matrices blocks for the PG-EAS static condensation

        if ((options & OPT_RHS) || (options & OPT_LHS)) {

            // check for the intila flag
            const Matrix* Cmat = (options & OPT_LHS_IS_INITIAL)
                ? &m_material[igauss]->getInitialTangent()
                : &m_material[igauss]->getTangent();
            // a non-finite material tangent (as opposed to a singular one,
            // which the condensation below detects and handles) would poison
            // every block it is integrated into before the inversion is even
            // attempted: swap it for the initial tangent at this gauss point.
            if (!(options & OPT_LHS_IS_INITIAL) && !isFinite(*Cmat)) {
                opserr << "ASDSolidHex::calculateAll() - element " << this->getTag()
                    << ", gauss point " << igauss << ", material " << m_material[igauss]->getTag()
                    << ": non-finite material tangent, using the initial one\n";
                Cmat = &m_material[igauss]->getInitialTangent();
            }

            C = *Cmat;

#ifdef ASDHEX_EAS_IMPERF
            // EXPERIMENT: the same factor that scaled sigma. Applied to the local
            // C after the initial/tangent choice above, so it lands consistently on
            // all four condensation blocks, on the OPT_LHS_IS_INITIAL path and on
            // getInitialStiff(). Scaling sigma without this would destroy the
            // consistency of the tangent.
            C *= m_imperfection_f[igauss];
#endif

            // the Damping object contributes a stiffness proportional term; it is
            // reported as a multiplier on the material tangent
            if (m_damping[igauss])
                C *= m_damping[igauss]->getStiffnessMultiplier();

            // compute K_uq += B_v^T C Bq_trial dV
            auto& BtC = ASDSolidHexGlobals::instance().BtC;
            BtC.addMatrixTransposeProduct(0.0, B_test, C, 1.0);
            m_eas->Kuq.addMatrixProduct(1.0, BtC, G_trial, dv);

            // compute K_qu += Bq_test^T C B_u dV
            auto& CtBu = ASDSolidHexGlobals::instance().CtBu;
            CtBu.addMatrixProduct(0.0, C, B_trial, 1.0);
            m_eas->Kqu.addMatrixTransposeProduct(1.0, G_test, CtBu, dv);

            // compute K_qq += Bq_test^T C Bq_trial dV
            auto& CtBq = ASDSolidHexGlobals::instance().CtBq;
            CtBq.addMatrixProduct(0.0, C, G_trial, 1.0);
            k_qq.addMatrixTransposeProduct(1.0, G_test, CtBq, dv);

            // compute K_uu += B_v^T C B_u dV
            k_uu.addMatrixProduct(1.0, BtC, B_trial, dv);

            // compute K_qq0 += Bq_test^T C0 Bq_trial dV, the A-op regularizer.
            //
            // Accumulated UNCONDITIONALLY, not only when the user asked for a
            // non-zero s: the gate after the loop raises s from zero on its own
            // when k_qq is not definite, and it cannot do that without this block
            // in hand. Making it lazy would mean re-entering the gauss loop, which
            // costs a whole material sweep -- 1296 mults per gauss point against
            // k_uu's 3456 is the far cheaper end of that trade.
            //
            // NOT cached across calls: getInitialTangent() is constant in time for
            // the ASD materials, but E and nu are live setParameter targets and
            // ASDSolidHex::setParameter broadcasts to all 8 of them, so a cache
            // would go silently stale. ASDSolidHexRefMetric is the wrong home for
            // it anyway -- that struct is reference GEOMETRY only.
            //
            // The damping multiplier is applied here too, so that on the
            // OPT_LHS_IS_INITIAL path, where C is already C0, k_qq0 comes out equal
            // to k_qq and the regularization degenerates to a plain scaling.
            auto& C0 = ASDSolidHexGlobals::instance().C0;
            C0 = m_material[igauss]->getInitialTangent();
#ifdef ASDHEX_EAS_IMPERF
            // the regularizer must be the elastic block of the SAME perturbed
            // material, or the gate would compare k_qq against a differently
            // scaled reference
            C0 *= m_imperfection_f[igauss];
#endif
            if (m_damping[igauss])
                C0 *= m_damping[igauss]->getStiffnessMultiplier();
            auto& C0tBq = ASDSolidHexGlobals::instance().C0tBq;
            C0tBq.addMatrixProduct(0.0, C0, G_trial, 1.0);
            k_qq0.addMatrixTransposeProduct(1.0, G_test, C0tBq, dv);

            // A-pen needs a DIFFERENT regularizer: the SYMMETRIC one, tested with
            // G_trial rather than G_test. That is not a stylistic choice -- it is
            // the matrix that Pi_stab = (p/2)*int((G_trial*alpha)^T C0
            // (G_trial*alpha) dV) actually differentiates to, so it is the only
            // one for which the penalty is the gradient of a potential and the
            // "energetic barrier" reading holds.
            //
            // It is also the only one with a definiteness guarantee, by
            // construction: a^T k_qq0_sym a = int((G_trial*a)^T C0 (G_trial*a) dV)
            // > 0 for any a != 0. The Petrov form k_qq0 has none -- measured,
            // sym(k_qq0) is not positive definite even elastically -- which does
            // not matter for A-op, where it is only a preconditioner, but would
            // make a "penalty" built on it meaningless.
            //
            // Reuses C0tBq: the product C0*G_trial is the same one.
            //
            // Gated, so the default path costs exactly what it did before.
            if (m_eas_penalty > 0.0)
                k_qq0_sym.addMatrixTransposeProduct(1.0, G_trial, C0tBq, dv);

            // How definite is the material tangent HERE, and how much of the
            // elastic tangent would have to be added to make it definite?
            //
            // For symmetric matrices lambda_min(C + s*C0) >= lambda_min(C) +
            // s*lambda_min(C0), so s >= (tol*scale - pivot(C)) / pivot(C0) is
            // sufficient. Min pivots stand in for the extreme eigenvalues: not
            // exact, but a cheap and conservative-in-practice surrogate, and C0 is
            // an isotropic elastic tangent whose smallest pivot is its shear term.
            //
            // On the OPT_LHS_IS_INITIAL path C is C0, so the pivots coincide, the
            // requirement is trivially met and no regularization is asked for.
            const double c0_scale = meanDiagonal(C0);
            const double c0_pivot = ldltMinPivot(C0);
            const double c_pivot = ldltMinPivot(C);
            const double c_pivot_rel = (c0_scale > 0.0) ? (c_pivot / c0_scale) : 0.0;
            if (!mat_pivot_seen || c_pivot_rel < mat_min_pivot) {
                mat_min_pivot = c_pivot_rel;
                mat_pivot_seen = true;
            }
            // Once ANY gauss point has lost definiteness, go straight to the cap
            // rather than to the smallest s that would fix THIS point.
            //
            // Grading s by the measured deficit was tried and is worse, for a
            // reason worth recording: the deficit is measured AFTER the tangent
            // has already degraded, so a graded s is applied one step too late and
            // stays too small while the divergence builds. Measured on the brittle
            // uniaxial test with -tangent: graded reached step 55 of 200 with the
            // gauss points 17*ft apart, while a constant s = 1 from the start
            // reached 199 of 200 with them 4e-11 apart -- and, because only the
            // OPERATOR is regularized, the second costs nothing in accuracy.
            if (m_eas_auto && c_pivot < EAS_MAT_PD_TOL * c0_scale)
                stab_needed = EAS_STAB_MAX;

        }

    } // End of Gauss Loop

    // AGQI: static condensation
    if ((options & OPT_RHS) || (options & OPT_LHS))
    {
        // ------------------------------------------------------------------
        // A-op: adaptive regularization of the enhanced block.
        //
        //     Kqq_s = k_qq + s * k_qq0
        //
        // Kqq_s is the operator inverted here, so it drives BOTH the static
        // condensation and -- through m_eas->Kqq_inv, which outlives this call --
        // the next updatePG_EAS. alpha_residual is left EXACTLY at
        // h = -int(G_test^T sigma dV), with no s term of any kind.
        //
        // THAT IS THE WHOLE POINT: at convergence h = 0, so the Kuq*Kqq_s^-1*h
        // contribution to the residual vanishes and the converged state is
        // INDEPENDENT of s. The patch tests and the bending accuracy are untouched
        // at any s, s may change freely from one iteration to the next, and there
        // is no new committed state to serialize. This is a Levenberg-Marquardt
        // regularization of the internal solve: it changes the path, never the
        // answer. What it buys is convergence -- for an enhanced mode with
        // eigenvalue lambda < 0, s*mu > 2|lambda| brings the rate back below one.
        //
        // It is NOT a penalty on the enhanced modes. A penalty would also subtract
        // s*k_qq0*alpha from h, which for an elastic material (where k_qq0 IS
        // k_qq) gives Kqq_s = (1+s)*k_qq exactly, hence alpha_s = alpha_0/(1+s)
        // and K_cond = Kuu - Kuq*Kqq^-1*Kqu/(1+s): s then interpolates smoothly to
        // the pure displacement element, costing about 3% of the bending
        // correction at s = 0.01 and 25% at s = 0.1, while a brittle softening
        // branch needs s of O(10-100) to be definite. A fixed penalty cannot be
        // both safe and accurate, which is why only the operator is touched here.
        //
        // s comes from the gauss loop, where it was sized pointwise on the
        // MATERIAL tangent (see EAS_MAT_PD_TOL). It is NOT searched for by testing
        // the assembled block: sym(k_qq) is not positive definite even elastically,
        // so no amount of s would ever satisfy such a test.
        //
        // The initial-tangent fallback below is kept as the last resort, but with
        // this in place it should essentially never be reached -- and that matters,
        // because that path installs inv(k_qq_elastic) into Kqq_inv without
        // correcting the residual, which is the divergent modified-Newton step
        // documented in updatePG_EAS.
        // ------------------------------------------------------------------
        m_eas->min_pivot = mat_min_pivot;
        if (stab_needed > EAS_STAB_MAX)
            stab_needed = EAS_STAB_MAX;
        // the user floor is NOT capped: -easStab is an explicit instruction
        const double stab = (stab_needed > m_eas_stab) ? stab_needed : m_eas_stab;
        m_eas->stab_used = stab;

        // ------------------------------------------------------------------
        // A-pen: the TRUE penalty, and the one line that separates it from A-op.
        //
        // The enhanced equation becomes int(G_test^T sigma dV) + p*Kqq0_sym*alpha
        // = 0, so with this file's convention h = -int(G_test^T sigma dV):
        //
        //     h_p   = h     - p * Kqq0_sym * alpha        <-- HERE, the residual
        //     Kqq_p = k_qq  + p * Kqq0_sym                <-- and below
        //
        // Kqq = -dh/dalpha, hence the plus on the operator against the minus on
        // the residual. A-op adds only the second of those two lines; that is the
        // entire difference between them, and everything else follows from it:
        // because p enters the RESIDUAL, at convergence int(G_test^T sigma dV) is
        // no longer zero, so the converged answer DEPENDS ON p. A-op's does not.
        //
        // Subtracted ONCE, after the gauss loop, not inside it: alpha is an
        // element quantity, not a per-gauss-point one. It has to happen before the
        // condensation feeds alpha_residual into the RHS below, and before the
        // member is left behind for the next updatePG_EAS to consume -- which is
        // what makes the lagged alpha step see the penalized residual too, as it
        // must for the two to describe the same problem.
        //
        // Applied on the OPT_LHS_IS_INITIAL path as well: the penalty is part of
        // the element's DEFINITION, not a property of the tangent, so an initial
        // stiffness computed without it would be inconsistent with the residual.
        // ------------------------------------------------------------------
        if (m_eas_penalty > 0.0)
            m_eas->alpha_residual.addMatrixVector(1.0, k_qq0_sym, m_eas->alpha, -m_eas_penalty);

        const Matrix* Kqq_op = &k_qq;
        if (stab > 0.0 || m_eas_penalty > 0.0) {
            auto& k_qq_stab = ASDSolidHexGlobals::instance().k_qq_stab;
            k_qq_stab = k_qq;
            if (stab > 0.0)
                k_qq_stab.addMatrix(1.0, k_qq0, stab);
            if (m_eas_penalty > 0.0)
                k_qq_stab.addMatrix(1.0, k_qq0_sym, m_eas_penalty);
            Kqq_op = &k_qq_stab;
        }

        // compute the inverse of the (possibly regularized) enhanced block.
        // NOTE: Invert() overwrites its argument, so there is no need to seed it.
        const double kqq_norm = normInf(*Kqq_op);
        auto& k_qq_inv = ASDSolidHexGlobals::instance().k_qq_inv;
        int info = Kqq_op->Invert(k_qq_inv);
        // A nearly singular k_qq can pass Invert with info = 0 and still produce
        // an inverse holding inf/nan, which then poisons LHS, RHS and -- through
        // Kqq_inv, which outlives this call -- the alpha update of every later
        // iteration. Catch it here, where the failure is still local.
        if (info == 0 && !isFinite(k_qq_inv))
            info = -1;
        if (info != 0 && !(options & OPT_LHS_IS_INITIAL)) {
            // k_qq is built from the MATERIAL tangent, so a damaging or perfectly
            // softening material can make it singular while the element state is
            // otherwise perfectly healthy. Instead of failing the evaluation,
            // rebuild the condensation from the INITIAL tangent, which stays
            // positive definite for the whole analysis.
            //
            // Everything is recomputed, not patched: keeping the material tangent
            // Kuq/Kqu/Kuu and pairing them with an initial tangent Kqq^-1 would mix
            // two different operators inside Kuq*Kqq^-1*Kqu. The internal force is
            // unaffected either way -- the RHS integrates B^T*sigma from the
            // material state, which no flag changes -- so the residual stays exact
            // and only the operator that eliminates the enhanced modes becomes
            // approximate: a modified Newton step on the enhanced equations. It
            // converges to the same solution, with more global iterations.
            //
            // OPT_UPDATE is cleared on the retry: the strains have already been
            // sent to the materials and updatePG_EAS has already advanced alpha for
            // this iteration, so repeating either would corrupt the state.
            opserr << "ASDSolidHex::calculateAll() - element " << this->getTag()
                << ": PG-EAS failed to invert k_qq (info = " << info
                << ") with the material tangent, retrying with the initial one\n";
            int retry = this->calculateAll(LHS, RHS, (options & ~OPT_UPDATE) | OPT_LHS_IS_INITIAL);
            // keep the material return codes collected by THIS pass: the retry
            // does not run setTrialStrain, so its own result is always zero
            return (retry < 0) ? retry : result;
        }
        if (info != 0) {
            // already on the initial tangent: nothing left to fall back to.
            // do NOT keep going with a garbage inverse: report the failure so the
            // algorithm can reduce the step instead of silently producing nonsense
            opserr << "ASDSolidHex::calculateAll() - element " << this->getTag()
                << ": PG-EAS failed to invert k_qq (info = " << info
                << ") with the initial tangent too\n";
            return -1;
        }

        m_eas->Kqq_inv = k_qq_inv;

        // reciprocal condition estimate, now that the inverse exists. Free: the
        // expensive half was already paid by Invert().
        const double kqq_inv_norm = normInf(k_qq_inv);
        m_eas->rcond = (kqq_norm > 0.0 && kqq_inv_norm > 0.0)
            ? 1.0 / (kqq_norm * kqq_inv_norm)
            : 0.0;

        // One warning per process, not per element per step: a softening material
        // point is the NORMAL state of every element in a localizing band, and
        // warning on each of them buries the analysis log. Nor does it fire merely
        // because the regularization engaged -- that is the mechanism working. It
        // fires when s has grown so large that the enhanced modes are effectively
        // gone (the element has degenerated towards a displacement element, so
        // whatever accuracy the EAS was there to provide is no longer there), or
        // when the operator is genuinely ill conditioned.
        //
        // Per-element, per-step detail is what the 'easState' recorder is for.
        //
        // This fires long BEFORE the "failed to invert k_qq" message below, which
        // by construction only reports the cases that already went singular.
        // just under the cap: the automatic path cannot exceed EAS_STAB_MAX, so a
        // threshold above it would only ever fire for an explicit -easStab
        constexpr double EAS_STAB_WARN = 0.9 * EAS_STAB_MAX;
        if (!(options & OPT_LHS_IS_INITIAL)) {
            static bool eas_conditioning_warned = false;
            if (!eas_conditioning_warned &&
                (stab > EAS_STAB_WARN || (m_eas->rcond > 0.0 && m_eas->rcond < 1.0e-10)))
            {
                eas_conditioning_warned = true;
                opserr << "ASDSolidHex::calculateAll() - element " << this->getTag()
                    << ": the PG-EAS enhanced block needed a large regularization"
                    << " (min pivot of the material tangent / elastic diagonal = "
                    << m_eas->min_pivot
                    << ", s = " << m_eas->stab_used
                    << ", rcond = " << m_eas->rcond << ").\n"
                    << "   A material tangent this far from definite is what an"
                    << " indefinite enhanced block is made of. Regularizing the"
                    << " operator bounds the damage but does not remove the cause:"
                    << " prefer the material's default secant tangent (drop"
                    << " -tangent) or -implex, both of which keep the tangent"
                    << " positive semi-definite. Record 'easState' to follow it."
                    << " Further warnings suppressed.\n";
            }
        }

        auto& K_uq_K_qq_inv = ASDSolidHexGlobals::instance().Kuq_Kqqinv;
        K_uq_K_qq_inv.addMatrixProduct(0.0, m_eas->Kuq, m_eas->Kqq_inv, 1.0);
        auto& K_uq_K_qq_inv_K_qu = ASDSolidHexGlobals::instance().Kuq_Kqqinv_Kqu;
        K_uq_K_qq_inv_K_qu.addMatrixProduct(0.0, K_uq_K_qq_inv, m_eas->Kqu, 1.0);

        if (options & OPT_RHS) {
            RHS.addMatrixVector(1.0, K_uq_K_qq_inv, m_eas->alpha_residual, 1.0);

        }

        if (options & OPT_LHS) {
            LHS.addMatrix(0.0, k_uu, 1.0);
            LHS.addMatrix(1.0, K_uq_K_qq_inv_K_qu, -1.0);
            // After the condensation, before returning LHS:

        }
    }
    // if corotational update from the local to global
    if (m_use_corotational) {
        m_transformation->transformToGlobal(local_cs, UG, UL, LHS, RHS, (options & OPT_LHS));
    }

    // Subtract external loads if any.
    // NOTE: this MUST come after transformToGlobal -- m_load is already a
    // vector of global components (see addInertiaLoadToUnbalance), so applying
    // it before the transformation would rotate it with the corotational frame.
    // Same order as ASDShellQ4::calculateAll.
    if ((options & OPT_RHS) && m_load)
        RHS.addVector(1.0, *m_load, -1.0);

    // Done
    return result;
}

void ASDSolidHex::updatePG_EAS(const Vector& U)
{
    // Static condensation of the enhanced parameters:
    //   d_alpha = -Kqq^-1 * (Kqu*dU - h)
    // where dU is the ITERATIVE increment of the local displacements and
    // h = alpha_residual is the residual of the enhanced equilibrium equation.
    // Same scheme as ASDShellQ4::AGQIupdate.
    static Vector dU(24);
    dU = U;
    dU.addVector(1.0, m_eas->U, -1.0);

    // save the current trial displacements for the next iteration
    m_eas->U = U;

    static Vector tmp(12);
    tmp.addMatrixVector(0.0, m_eas->Kqu, dU, 1.0);
    tmp.addVector(1.0, m_eas->alpha_residual, -1.0);

    static Vector dalpha(12);
    dalpha.addMatrixVector(0.0, m_eas->Kqq_inv, tmp, -1.0);

    // ------------------------------------------------------------------
    // TRUST REGION on the enhanced increment.
    //
    // This is ONE lagged linear step per global iteration, and it is only a
    // consistent Newton step while Kqq is positive definite. Under softening it
    // is not: k_qq goes indefinite, and when it goes singular the condensation
    // rebuilds Kqq_inv from the INITIAL tangent and that inverse outlives the
    // call, so the next iteration lands here as a modified Newton step with
    // convergence rate rho = |1 - C_t/C_0| -- which exceeds 1 for ANY negative
    // material tangent, not merely a steep one. A single such step can throw
    // alpha far enough that the materials are handed a meaningless strain and
    // never recover, because nothing downstream bounds it: alpha is not a
    // solver-owned unknown, so no line search or step reduction ever sees it.
    //
    // The cap does not change the converged answer: at convergence dalpha -> 0
    // and the limit is inactive. It only stops one bad iterate from destroying
    // the state, at worst costing a few extra global iterations, since alpha may
    // still grow geometrically (GROWTH-fold per iteration).
    //
    // alpha is length-like (G_trial carries 1/length, exactly as B does), so
    // both terms of the limit are lengths and the criterion is dimensionally
    // consistent. The dU term lets alpha grow from zero at the start of the
    // analysis; the U term keeps the limit finite when a commit moved the
    // material state without moving the nodes, which is where a pure ||alpha||
    // bound would leave dalpha unbounded.
    //
    // A zero limit means the element is completely undeformed: h is zero too,
    // so let the (zero) increment through rather than scaling by 0/0.
    //
    // ARMED ONLY WHEN THE ENHANCED BLOCK IS ACTUALLY INDEFINITE, which is what
    // min_pivot (measured on the raw k_qq by the last condensation) reports. While
    // the block is positive definite this IS a consistent Newton step and clamping
    // it would be wrong -- and, more practically, alpha starts every analysis at
    // zero, so the ||alpha|| term vanishes on the first iteration and a cap left
    // permanently armed could fire on a perfectly healthy elastic element, moving
    // results that are currently correct. Linear and hardening analyses therefore
    // never reach the body below.
    //
    // min_pivot is zero before the first condensation has run, but so is Kqq_inv
    // (initializePG_EAS zeroes it), hence dalpha is zero there and the clamp is a
    // no-op regardless.
    // ------------------------------------------------------------------
    constexpr double EAS_TRUST_GROWTH = 4.0;    // alpha may grow 5x per iteration
    constexpr double EAS_TRUST_REL = 1.0e-3;    // floor: 0.1% of the current motion
    const double dalpha_norm = dalpha.Norm();
    const double trust = EAS_TRUST_GROWTH * m_eas->alpha.Norm()
        + dU.Norm() + EAS_TRUST_REL * U.Norm();
    if (m_eas->min_pivot <= 0.0 && trust > 0.0 && dalpha_norm > trust) {
        dalpha *= (trust / dalpha_norm);
        static bool eas_trust_warned = false;
        if (!eas_trust_warned) {
            eas_trust_warned = true;
            opserr << "ASDSolidHex::updatePG_EAS() - element " << this->getTag()
                << ": the enhanced increment was clamped by the trust region"
                << " (|dalpha| = " << dalpha_norm << " > " << trust << ").\n"
                << "   The enhanced solve is diverging, typically because k_qq is"
                << " no longer positive definite. Record 'easState' to follow it."
                << " Further warnings suppressed.\n";
        }
    }

    m_eas->alpha += dalpha;
}

void ASDSolidHex::initializePG_EAS()
{
    // Reset the enhanced state to the undeformed configuration.
    //
    // NOTE: this used to do  m_eas->U = m_eas->U_commit = ASDSolidHexGlobals::UL,
    // i.e. it copied whatever local displacement vector the LAST element to run
    // had left behind in the shared singleton. That happens to be zero when
    // called from setDomain (nothing has run calculateAll yet) but is arbitrary
    // when called from revertToStart() in the middle of an analysis.
    //
    // Both call sites (setDomain, revertToStart) refer to the undeformed state,
    // so the correct seed is exactly zero. Deriving it from the current node
    // state instead -- as ASDShellQ4::AGQIinitialize does -- was tried and
    // measurably perturbed the converged corotational results, because
    // calculateLocalDisplacements() does not return exactly zero when called at
    // setDomain time; zero is both correct and reproducible.
    m_eas->U.Zero();
    m_eas->U_converged.Zero();

    // initialize the EAS internal parameters (alpha)
    m_eas->alpha.Zero();
    m_eas->alpha_commit.Zero();
    m_eas->alpha_residual.Zero();
    m_eas->Kqq_inv.Zero();
    m_eas->Kqu.Zero();
    m_eas->Kuq.Zero();

    // the diagnostics are not state, but a stale rcond / pivot / s surviving a
    // revertToStart would be read by the 'easState' recorder as if it described
    // the fresh element. m_eas_stab is deliberately NOT touched: that one IS a
    // user parameter, not state.
    m_eas->rcond = 0.0;
    m_eas->min_pivot = 0.0;
    m_eas->stab_used = 0.0;
}


int
ASDSolidHex::displaySelf(Renderer& theViewer, int displayMode, float fact, const char** modes, int numMode)
{
    // Nothing to draw before setDomain has bound the nodes. The broker builds
    // this element with null node pointers and a renderer can reach it first.
    for (int i = 0; i < NumNodes; i++)
        if (nodePtrs[i] == nullptr)
            return 0;

    // All of this used to be static. That is not a style point: a Renderer may
    // walk elements from more than one thread, and two elements sharing one
    // coords/values buffer draw each other's geometry. They are locals now.
    Vector v[NumNodes] = {
        Vector(3), Vector(3), Vector(3), Vector(3),
        Vector(3), Vector(3), Vector(3), Vector(3)
    };
    for (int i = 0; i < NumNodes; i++)
        nodePtrs[i]->getDisplayCrds(v[i], fact, displayMode);

    // The colour carried to each vertex. It used to be hard-wired to zero, so
    // every face came out flat whatever the display mode asked for.
    //
    // displayMode 1..6 selects a stress component, in the same spirit as
    // Brick::displaySelf, and each vertex takes the stress of the gauss point
    // nearest it. That nearest point is simply the one with the SAME INDEX:
    //
    //   node a       xi_n  = { -1,+1,+1,-1,-1,+1,+1,-1 }  (dshape)
    //   gauss pt a   XI    = { -g,+g,+g,-g,-g,+g,+g,-g }
    //
    // and likewise for eta and zeta, so gauss point a sits in the corner of
    // node a. Note that this is the INTERNAL gauss order and the mapping must
    // NOT go through GP_REPORT_TO_INTERNAL: that permutation exists only for the
    // recorder layer, which reports in the lexicographic i-j-k order of the other
    // bricks. Sending display values through it would rotate the colours around
    // the element.
    double nodeValue[NumNodes] = { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 };
    if (displayMode > 0 && displayMode <= 6) {
        const int comp = displayMode - 1;
        for (int i = 0; i < NumNodes; i++) {
            if (m_material[i] == nullptr)
                continue;
            const Vector& sig = m_material[i]->getStress();
            if (comp < sig.Size())
                nodeValue[i] = sig(comp);
        }
    }

    Matrix coords(4, 3);
    Vector values(4);

    auto drawFace = [&](int n0, int n1, int n2, int n3) -> int {
        const int n[4] = { n0, n1, n2, n3 };
        for (int c = 0; c < 4; c++) {
            for (int k = 0; k < 3; k++)
                coords(c, k) = v[n[c]](k);
            values(c) = nodeValue[n[c]];
        }
        return theViewer.drawPolygon(coords, values, this->getTag());
        };

    int res = 0;

    // Draw 6 faces

    // bottom
    res += drawFace(0, 1, 2, 3);
    // top
    res += drawFace(4, 5, 6, 7);
    // lateral faces
    res += drawFace(0, 1, 5, 4);
    res += drawFace(1, 2, 6, 5);
    res += drawFace(2, 3, 7, 6);
    res += drawFace(3, 0, 4, 7);

    return res;

}


