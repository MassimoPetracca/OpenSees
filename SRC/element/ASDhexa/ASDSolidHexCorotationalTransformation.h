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
** ****************************************************************** */

// Original implementation: Massimo Petracca (ASDEA)
//
// Implementation of a corotational coordinate transformation for
// 8-node solid hexahedral elements with 3 translational DOFs per node.
//
// Key difference from ASDShellQ4CorotationalTransformation:
// - 8 nodes instead of 4
// - 3 DOFs per node (translations only), NO rotational DOFs
// - No nodal quaternions needed (no nodal rotations to track)
// - The corotational frame carries all rigid body rotation info
//

#ifndef ASDSolidHexCorotationalTransformation_h
#define ASDSolidHexCorotationalTransformation_h


// === THESE MACROS MUST BE DEFINED BEFORE ASDSolidHexLocalCoordinateSystem.h ===
#define USE_POLAR_DECOMP_ALIGN
#define USE_OLD_FELIPPA_FRAME 0
// CR-frame: 1 = polar decomp of the volume-weighted average F̄ over the 2x2x2 GPs (first-order best-fit)
//           0 = polar decomp of F evaluated at the parametric center only (old behavior)
#define USE_FBAR_GP_AVG 0
// CR-frame: 1 = nodal Kabsch/Procrustes best-fit (= discrete Rankin on the nodes):
//               R = polar(A), A = Σ_i (current_i - C̄) ⊗ (initial_i - C̄_0)
//           0 = use what USE_FBAR_GP_AVG selects (F at the center or averaged F̄)
// When 1, it overrides USE_FBAR_GP_AVG.
#define USE_KABSCH_RIGID 0
// LCS frame composition: 1 = m_R = matrixPtr^T directly (assumes Rtilde_init = I)
//                        0 = m_R = matrixPtr^T * m_Rtilde_curr (old behavior, doubles the rotation)
// To be used when matrixPtr is the "active" best-fit rotation (e.g. polar/Kabsch).
#define USE_SIMPLE_RFRAME 1

#include <ASDEICR3.h>
#include <ASDSolidHexLocalCoordinateSystem.h>
#include <Node.h>
#include <Domain.h>
#include <ID.h>

/**
* \brief ASDSolidHexCorotationalTransformation
*
* Corotational (nonlinear) coordinate transformation for 8-node
* hexahedral solid elements with 3 translational DOFs per node
* (24 DOFs total).
*
* Its main aim is to:
* 1) Create the local corotational coordinate system from nodal positions
* 2) Transform incoming global displacements to local deformational
*    displacements, removing rigid body translations and rotations.
* 3) Transform outgoing matrices and vectors back to the global
*    coordinate system including rigid body motion effects.
*
* Since there are NO rotational DOFs at the nodes:
* - No nodal quaternions are needed
* - No incremental rotation tracking is needed
* - The update() method is a no-op
* - Only the element frame rotation R0 = T0^T * TR matters
*
* Node numbering (standard Hex8):
*
*      8 -------- 7
*     /|          /|
*    5 -------- 6  |
*    |  |       |  |
*    |  4 ------|- 3
*    | /        | /
*    1 -------- 2
*
* References:
* - C.A. Felippa, B. Haugen, "Unified formulation of small-strain
*   corotational finite elements: I. Theory",
*   CU-CAS-05-02, January 2005
*/
class ASDSolidHexCorotationalTransformation
{

public:

    typedef ASDVector3<double> Vector3Type;

    typedef ASDQuaternion<double> QuaternionType;

    typedef Vector VectorType;

    typedef Matrix MatrixType;

    typedef std::array<Node*, 8> NodeContainerType;

public:

    ASDSolidHexCorotationalTransformation()
    {
    }

    ~ASDSolidHexCorotationalTransformation()
    {
    }

public:


    bool isLinear() const
    {
        return false;
    }

    // -----------------------------------------------------------------------
    // revertToStart
    // Called once at initialization.
    // Builds T0 (reference frame) from the undeformed nodal positions and
    // saves the reference centroid C0.
    // No nodal quaternions to initialize (no rotational DOFs).
    // -----------------------------------------------------------------------
    void revertToStart()
    {
        // compute reference (undeformed) coordinate system
        // using the 8 nodes in their initial configuration
        std::array<Vector3Type, 8> X0;
        for (int i = 0; i < 8; i++)
            X0[i] = Vector3Type(m_nodes[i]->getCrds());

        // build the reference coordinate system from undeformed nodal positions
        ASDSolidHexLocalCoordinateSystem LCS0(X0[0], X0[1], X0[2], X0[3], X0[4], X0[5], X0[6], X0[7]);

        // save the reference centroid C0 from the reference coordinate system
        m_C0 = LCS0.Origin();

        // save the T0 matrix in the reference coordiante system 
        m_Q0 = QuaternionType::FromRotationMatrix(LCS0.getRotationMatrix());


    }

    // -----------------------------------------------------------------------
    // setDomain
    // -----------------------------------------------------------------------
    void setDomain(Domain* domain, const ID& node_ids, bool initialized)
    {
        // call base class to get nodes and save initial displacements
       // if domain is null
        if (domain == nullptr) {
            for (size_t i = 0; i < 8; i++) {
                m_nodes[i] = nullptr;
            }
            return;
        }

        // get nodes and save initial displacements and rotations
        for (size_t i = 0; i < 8; i++) {
            m_nodes[i] = domain->getNode(node_ids(i));

            if (m_nodes[i] == nullptr) {
                opserr << "ASDSolidHexCorotationalTransformation::setDomain - no node " << node_ids(i)
                    << " exists in the model\n";
                exit(-1);
            }
            if (!initialized) {
                const Vector& iU = m_nodes[i]->getTrialDisp();
                if (iU.Size() != 3) {
                    opserr << "ASDSolidHexCorotationalTransformation::setDomain - node " << node_ids(i)
                        << " has " << iU.Size() << " DOFs, while 3 are expected\n";
                    exit(-1);
                }
                size_t index = i * 3;
                for (size_t j = 0; j < 3; j++)
                    m_U0(index + j) = iU(j);
            }
        }

        // The reference shape function gradients depend only on the reference
        // geometry, so they are computed unconditionally -- also on the recvSelf
        // path, where 'initialized' is true and revertToStart() is skipped.
        computeCentreGradients();

        // quick return
        if (domain == nullptr || initialized)
            return;

        // init state variables
        revertToStart();
    }

    // -----------------------------------------------------------------------
    // revertToLastCommit
    // No nodal rotation state to revert for Hex8.
    // -----------------------------------------------------------------------
    void revertToLastCommit()
    {
        // nothing to do: no nodal quaternions to revert
    }

    // -----------------------------------------------------------------------
    // commit
    // No nodal rotation state to commit for Hex8.
    // -----------------------------------------------------------------------
    void commit()
    {
        // nothing to do: no nodal quaternions to commit
    }

    // -----------------------------------------------------------------------
    // update
    // No nodal rotations to track for Hex8.
    // The corotational frame is recomputed from scratch at each call to
    // createLocalCoordinateSystem(), so no incremental state is needed here.
    // -----------------------------------------------------------------------
    void update(const VectorType& globalDisplacements)
    {
        // nothing to do: no nodal rotational DOFs to update
    }

    // -----------------------------------------------------------------------
    // createLocalCoordinateSystem
    // Computes the corotational frame TR from the deformed nodal positions,
    // then extracts the rigid body rotation R0 = T0^T * TR.
    // Returns a local coordinate system object carrying TR and the
    // deformed centroid.
    // -----------------------------------------------------------------------
    ASDSolidHexLocalCoordinateSystem createReferenceCoordinateSystem() const
    {   // the reference coordinate system in the underformed configuration
        // using the default alignment to the first column of the jacobian at center

        return ASDSolidHexLocalCoordinateSystem(
            Vector3Type(m_nodes[0]->getCrds()),
            Vector3Type(m_nodes[1]->getCrds()),
            Vector3Type(m_nodes[2]->getCrds()),
            Vector3Type(m_nodes[3]->getCrds()),
            Vector3Type(m_nodes[4]->getCrds()),
            Vector3Type(m_nodes[5]->getCrds()),
            Vector3Type(m_nodes[6]->getCrds()),
			Vector3Type(m_nodes[7]->getCrds()));
     }


    ASDSolidHexLocalCoordinateSystem createLocalCoordinateSystem(
        const VectorType& globalDisplacements) const
    {
        // create a reference coordinate system 
		ASDSolidHexLocalCoordinateSystem a = createReferenceCoordinateSystem();

        // compute nodal positions at current configuration removing initial displacements if any
        std::array<Vector3Type, 8> def = {
            Vector3Type(m_nodes[0]->getCrds()),
            Vector3Type(m_nodes[1]->getCrds()),
            Vector3Type(m_nodes[2]->getCrds()),
            Vector3Type(m_nodes[3]->getCrds()),
            Vector3Type(m_nodes[4]->getCrds()),
            Vector3Type(m_nodes[5]->getCrds()),
            Vector3Type(m_nodes[6]->getCrds()),
            Vector3Type(m_nodes[7]->getCrds()) };   

        for (int i = 0; i < 8; i++) {
            int index = i * 3;
            Vector3Type& iP = def[i];
			iP(0) += globalDisplacements(index) - m_U0(index);
			iP(1) += globalDisplacements(index + 1) - m_U0(index + 1);
			iP(2) += globalDisplacements(index + 2) - m_U0(index + 2);
        }
        // build b in GLOBAL axes (no matrixPtr -> m_R = I) so that the deformation
        // gradient F = J_b * J_a^-1 is the true global F (independent of the edge
        // frame Rtilde and hence of node ordering). polar(F) is then the rigid
        // rotation: = I for pure strain (patch-test exact) and = R for rigid motion.
		ASDSolidHexLocalCoordinateSystem b(def[0], def[1], def[2], def[3], def[4], def[5], def[6], def[7]);


#if 0
        return b;
#else // !0

        // compute the deformation gradient F -> correlation matrix
        // Local coordinates in the reference frame (undeformed configuration)
        const double aX[8] = { a.X1(), a.X2(), a.X3(), a.X4(), a.X5(), a.X6(), a.X7(), a.X8() };
        const double aY[8] = { a.Y1(), a.Y2(), a.Y3(), a.Y4(), a.Y5(), a.Y6(), a.Y7(), a.Y8() };
        const double aZ[8] = { a.Z1(), a.Z2(), a.Z3(), a.Z4(), a.Z5(), a.Z6(), a.Z7(), a.Z8() };

        // Local coordinates in the current frame (deformed configuration)
        const double bX[8] = { b.X1(), b.X2(), b.X3(), b.X4(), b.X5(), b.X6(), b.X7(), b.X8() };
        const double bY[8] = { b.Y1(), b.Y2(), b.Y3(), b.Y4(), b.Y5(), b.Y6(), b.Y7(), b.Y8() };
        const double bZ[8] = { b.Z1(), b.Z2(), b.Z3(), b.Z4(), b.Z5(), b.Z6(), b.Z7(), b.Z8() };

        // Node indices of the unit cube in natural coordinates
        static const int xi_n[8]   = { -1, +1, +1, -1, -1, +1, +1, -1 };
        static const int eta_n[8]  = { -1, -1, +1, +1, -1, -1, +1, +1 };
        static const int zeta_n[8] = { -1, -1, -1, -1, +1, +1, +1, +1 };

        // Lambda: computes F(xi,eta,zeta) and |J_a(xi,eta,zeta)|
        auto compute_F_at = [&](double xi, double eta, double zeta,
                                double F_out[3][3], double& detJa_out)
        {
            double dN[8][3];
            for (int n = 0; n < 8; n++) {
                dN[n][0] = 0.125 * xi_n[n]   * (1.0 + eta_n[n] * eta) * (1.0 + zeta_n[n] * zeta);
                dN[n][1] = 0.125 * eta_n[n]  * (1.0 + xi_n[n]  * xi)  * (1.0 + zeta_n[n] * zeta);
                dN[n][2] = 0.125 * zeta_n[n] * (1.0 + xi_n[n]  * xi)  * (1.0 + eta_n[n]  * eta);
            }
            double Ja[3][3] = { {0,0,0},{0,0,0},{0,0,0} };
            double Jb[3][3] = { {0,0,0},{0,0,0},{0,0,0} };
            for (int n = 0; n < 8; n++) {
                for (int k = 0; k < 3; k++) {
                    Ja[0][k] += dN[n][k] * aX[n];
                    Ja[1][k] += dN[n][k] * aY[n];
                    Ja[2][k] += dN[n][k] * aZ[n];
                    Jb[0][k] += dN[n][k] * bX[n];
                    Jb[1][k] += dN[n][k] * bY[n];
                    Jb[2][k] += dN[n][k] * bZ[n];
                }
            }
            const double dJa = Ja[0][0] * (Ja[1][1] * Ja[2][2] - Ja[1][2] * Ja[2][1])
                             - Ja[0][1] * (Ja[1][0] * Ja[2][2] - Ja[1][2] * Ja[2][0])
                             + Ja[0][2] * (Ja[1][0] * Ja[2][1] - Ja[1][1] * Ja[2][0]);
            detJa_out = dJa;
            const double inv_dJa = 1.0 / dJa;
            double iJa[3][3];
            iJa[0][0] = (Ja[1][1] * Ja[2][2] - Ja[1][2] * Ja[2][1]) * inv_dJa;
            iJa[0][1] = (Ja[0][2] * Ja[2][1] - Ja[0][1] * Ja[2][2]) * inv_dJa;
            iJa[0][2] = (Ja[0][1] * Ja[1][2] - Ja[0][2] * Ja[1][1]) * inv_dJa;
            iJa[1][0] = (Ja[1][2] * Ja[2][0] - Ja[1][0] * Ja[2][2]) * inv_dJa;
            iJa[1][1] = (Ja[0][0] * Ja[2][2] - Ja[0][2] * Ja[2][0]) * inv_dJa;
            iJa[1][2] = (Ja[0][2] * Ja[1][0] - Ja[0][0] * Ja[1][2]) * inv_dJa;
            iJa[2][0] = (Ja[1][0] * Ja[2][1] - Ja[1][1] * Ja[2][0]) * inv_dJa;
            iJa[2][1] = (Ja[0][1] * Ja[2][0] - Ja[0][0] * Ja[2][1]) * inv_dJa;
            iJa[2][2] = (Ja[0][0] * Ja[1][1] - Ja[0][1] * Ja[1][0]) * inv_dJa;
            for (int i = 0; i < 3; i++)
                for (int j = 0; j < 3; j++) {
                    double s = 0.0;
                    for (int k = 0; k < 3; k++) s += Jb[i][k] * iJa[k][j];
                    F_out[i][j] = s;
                }
        };

        Matrix F(3, 3);
        F.Zero();

#if USE_KABSCH_RIGID
        // ===== Nodal Kabsch / Procrustes =====
        // R minimizing Σ_i ||R · (initial_i - C_0) - (current_i - C̄)||²
        // Build A = Σ_i (current_i - C̄) ⊗ (initial_i - C_0); the polar decomp
        // of A (performed by the Higham code below) then returns the best-fit
        // rotation directly. Here F holds the covariance A.
        // The aX[], aY[], aZ[] / bX[], bY[], bZ[] coordinates above are ALREADY
        // centered, because ASDSolidHexLocalCoordinateSystem.m_P[i] = m_R * (P_i - origin).
        // For the reference LCS (a) m_R = I → aX/Y/Z are in the global frame, centered.
        // For the current LCS   (b) m_R = R_polar^T * Rtilde → bX/Y/Z are already in the CR frame.
        // The nodal Kabsch best-fit needs both in the same (global) frame, so the
        // centered current coordinates are rebuilt in global below.
        double Ccurr[3] = { 0.0, 0.0, 0.0 };
        double Cinit[3] = { 0.0, 0.0, 0.0 };
        double Ycurr_g[8][3];
        double Xinit_g[8][3];
        for (int n = 0; n < 8; n++) {
            const auto& Pn = m_nodes[n]->getCrds();
            Xinit_g[n][0] = Pn(0);
            Xinit_g[n][1] = Pn(1);
            Xinit_g[n][2] = Pn(2);
            Ycurr_g[n][0] = Pn(0) + globalDisplacements(3 * n + 0) - m_U0(3 * n + 0);
            Ycurr_g[n][1] = Pn(1) + globalDisplacements(3 * n + 1) - m_U0(3 * n + 1);
            Ycurr_g[n][2] = Pn(2) + globalDisplacements(3 * n + 2) - m_U0(3 * n + 2);
            for (int k = 0; k < 3; k++) {
                Cinit[k] += Xinit_g[n][k];
                Ccurr[k] += Ycurr_g[n][k];
            }
        }
        for (int k = 0; k < 3; k++) { Cinit[k] *= 0.125; Ccurr[k] *= 0.125; }
        // Build A = Σ Y_i ⊗ X_i  (3x3)
        for (int n = 0; n < 8; n++) {
            const double yx = Ycurr_g[n][0] - Ccurr[0];
            const double yy = Ycurr_g[n][1] - Ccurr[1];
            const double yz = Ycurr_g[n][2] - Ccurr[2];
            const double xx = Xinit_g[n][0] - Cinit[0];
            const double xy = Xinit_g[n][1] - Cinit[1];
            const double xz = Xinit_g[n][2] - Cinit[2];
            F(0, 0) += yx * xx; F(0, 1) += yx * xy; F(0, 2) += yx * xz;
            F(1, 0) += yy * xx; F(1, 1) += yy * xy; F(1, 2) += yy * xz;
            F(2, 0) += yz * xx; F(2, 1) += yz * xy; F(2, 2) += yz * xz;
        }
#elif USE_FBAR_GP_AVG
        // ===== Volume-weighted average F̄ over the 2x2x2 Gauss points =====
        // F̄ = (1/V_a) Σ_g w_g · |J_a(xi_g)| · F(xi_g)   (w_g = 1 per 2x2x2 Gauss)
        const double g_loc = 1.0 / 1.7320508075688772;  // 1/sqrt(3)
        const double xi_gp[8]   = { -g_loc, +g_loc, +g_loc, -g_loc, -g_loc, +g_loc, +g_loc, -g_loc };
        const double eta_gp[8]  = { -g_loc, -g_loc, +g_loc, +g_loc, -g_loc, -g_loc, +g_loc, +g_loc };
        const double zeta_gp[8] = { -g_loc, -g_loc, -g_loc, -g_loc, +g_loc, +g_loc, +g_loc, +g_loc };
        double V_a = 0.0;
        double Fgp[3][3];
        double detJa_gp;
        for (int gp = 0; gp < 8; gp++) {
            compute_F_at(xi_gp[gp], eta_gp[gp], zeta_gp[gp], Fgp, detJa_gp);
            for (int i = 0; i < 3; i++)
                for (int j = 0; j < 3; j++)
                    F(i, j) += Fgp[i][j] * detJa_gp;
            V_a += detJa_gp;
        }
        const double inv_Va = 1.0 / V_a;
        for (int i = 0; i < 3; i++)
            for (int j = 0; j < 3; j++)
                F(i, j) *= inv_Va;
#else
        // ===== F evaluated at the parametric center only (xi=eta=zeta=0) =====
        double Fc[3][3];
        double detJa_c;
        compute_F_at(0.0, 0.0, 0.0, Fc, detJa_c);
        for (int i = 0; i < 3; i++)
            for (int j = 0; j < 3; j++)
                F(i, j) = Fc[i][j];
#endif

        // Polar decomposition by Higham's iteration
        // F = R * U, where R is the rotation and U the right stretch
        Matrix R(3, 3);
        Matrix R_old(3, 3);
        Matrix R_inv_trans(3, 3);

        // Initialize R = F
        for (int i = 0; i < 3; i++) {
            for (int j = 0; j < 3; j++) {
                R(i, j) = F(i, j);
            }
        }

        // Newton iteration for the polar decomposition (Higham's method)
        const int maxIter = 20;
        const double tol = 1e-12;
        for (int iter = 0; iter < maxIter; iter++) {
            // Salva R precedente per controllo convergenza
            for (int i = 0; i < 3; i++) {
                for (int j = 0; j < 3; j++) {
                    R_old(i, j) = R(i, j);
                }
            }

            // Compute the inverse of R
            double detR = R(0, 0) * (R(1, 1) * R(2, 2) - R(1, 2) * R(2, 1)) -
                R(0, 1) * (R(1, 0) * R(2, 2) - R(1, 2) * R(2, 0)) +
                R(0, 2) * (R(1, 0) * R(2, 1) - R(1, 1) * R(2, 0));
            double inv_detR = 1.0 / detR;

            double invR11 = (R(1, 1) * R(2, 2) - R(1, 2) * R(2, 1)) * inv_detR;
            double invR12 = (R(0, 2) * R(2, 1) - R(0, 1) * R(2, 2)) * inv_detR;
            double invR13 = (R(0, 1) * R(1, 2) - R(0, 2) * R(1, 1)) * inv_detR;
            double invR21 = (R(1, 2) * R(2, 0) - R(1, 0) * R(2, 2)) * inv_detR;
            double invR22 = (R(0, 0) * R(2, 2) - R(0, 2) * R(2, 0)) * inv_detR;
            double invR23 = (R(0, 2) * R(1, 0) - R(0, 0) * R(1, 2)) * inv_detR;
            double invR31 = (R(1, 0) * R(2, 1) - R(1, 1) * R(2, 0)) * inv_detR;
            double invR32 = (R(0, 1) * R(2, 0) - R(0, 0) * R(2, 1)) * inv_detR;
            double invR33 = (R(0, 0) * R(1, 1) - R(0, 1) * R(1, 0)) * inv_detR;

            // R_inv_trans = (R^{-1})^T
            R_inv_trans(0, 0) = invR11;
            R_inv_trans(0, 1) = invR21;
            R_inv_trans(0, 2) = invR31;
            R_inv_trans(1, 0) = invR12;
            R_inv_trans(1, 1) = invR22;
            R_inv_trans(1, 2) = invR32;
            R_inv_trans(2, 0) = invR13;
            R_inv_trans(2, 1) = invR23;
            R_inv_trans(2, 2) = invR33;

            // R = 0.5 * (R + R_inv_trans)
            for (int i = 0; i < 3; i++) {
                for (int j = 0; j < 3; j++) {
                    R(i, j) = 0.5 * (R(i, j) + R_inv_trans(i, j));
                }
            }

            // Controllo convergenza
            double norm_diff = 0.0;
            for (int i = 0; i < 3; i++) {
                for (int j = 0; j < 3; j++) {
                    double diff = R(i, j) - R_old(i, j);
                    norm_diff += diff * diff;
                }
            }
            if (sqrt(norm_diff) < tol) break;
        }

        // R now holds the rotation matrix
        // U can be computed as U = R^T * F
        Matrix U(3, 3);
        for (int i = 0; i < 3; i++) {
            for (int j = 0; j < 3; j++) {
                U(i, j) = 0.0;
                for (int k = 0; k < 3; k++) {
                    U(i, j) += R(k, i) * F(k, j);  // R^T * F
                }
            }
        }


        //R.addMatrixProduct(0.0, R, a.Orientation(), 1.0);
        ASDSolidHexLocalCoordinateSystem c(def[0], def[1], def[2], def[3], def[4], def[5], def[6], def[7], &R);

        return c;
#endif
    }


    void calculateLocalDisplacements(
        const ASDSolidHexLocalCoordinateSystem& LCS,
        const VectorType& globalDisplacements,
        VectorType& localDisplacements)
    {

        // orientation and center of the current local coordinate systen
		QuaternionType Q = QuaternionType::FromRotationMatrix(LCS.getRotationMatrix());
        const Vector3Type& C = LCS.Origin();


        for (int i = 0; i < 8; i++) {

			int index = i * 3;

            // centered underformed configuration
            Vector3Type X0 = Vector3Type(m_nodes[i]->getCrds());
            X0 -= m_C0;

            // entered deformed position
            Vector3Type X = Vector3Type(m_nodes[i]->getCrds())
                + Vector3Type(globalDisplacements, index) - C;

            // get deformatinoal displacmeents
            Q.rotateVector(X);
            m_Q0.rotateVector(X0);
			Vector3Type deformationalDisplacement = X - X0;


            localDisplacements(index) = deformationalDisplacement(0);
            localDisplacements(index + 1) = deformationalDisplacement(1);
			localDisplacements(index + 2) = deformationalDisplacement(2);

        }

    }

    void transformToGlobal(
        const ASDSolidHexLocalCoordinateSystem& LCS,
        const VectorType& globalDisplacements,
        const VectorType& localDisplacements,
        MatrixType& LHS,
        VectorType& RHS,
        bool LHSrequired)
    {

        // T: bloc-diagonal rotation matrix (global-to-local), 24x24
        // Each 3x3 diagonal block is the orientation matrix of the CR frame
        static MatrixType T(24, 24);
        LCS.ComputeTotalRotationMatrix(T);

        // P: projector, 24x24. Pu alone removes only the rigid TRANSLATIONS;
        // the rotational part comes from -S*G, which accounts for the fact that
        // the corotational frame itself rotates when the nodes move.
        // Pu depends only on the node count, so it is built once. It used to be
        // rebuilt on every call. P itself must stay a separate buffer because the
        // -S*G term is subtracted into it in place.
        static const MatrixType Pu = [] {
            MatrixType m(24, 24);
            ASDEICR3::Compute_Pt_3dof(8, m);  // 8 nodes, 3 DOF each
            return m;
        }();
        static MatrixType P(24, 24);
        P = Pu;

        // G: rotation-gradient (spin-lever) matrix, 3x24, dw = G*dv
        static MatrixType G(3, 24);
        RotationGradient(LCS, globalDisplacements, G);

        // S: 24x3, S_a = -Spin(x_bar_a), so that du_local = (Pu - S*G) dv
        static MatrixType S(24, 3);
        SpinLeverMatrix(LCS, S);

        // full projector P = Pu - S*G
        static MatrixType tmp(24, 24);
        tmp.addMatrixProduct(0.0, S, G, 1.0);
        P.addMatrix(1.0, tmp, -1.0);

        // Projected local forces: pe = P^T * f_local
        // (RHS coming in is already the negative internal force vector)
        static VectorType projectedLocalForces(24);
        projectedLocalForces.addMatrixTransposeVector(0.0, P, RHS, 1.0);

        // Global RHS: T^T * P^T * f_local
        RHS.addMatrixTransposeVector(0.0, T, projectedLocalForces, 1.0);

        if (!LHSrequired)
            return;

        // --- Step 1: Material stiffness K.M ---
        // Ke = P^T * Km * P
        // Note: no H matrix because there are no rotational DOFs
        static MatrixType temp(24, 24);
        temp.addMatrixProduct(0.0, LHS, P, 1.0);
        LHS.addMatrixTransposeProduct(0.0, P, temp, 1.0);

        // --- Step 2: Equilibrium projection geometric stiffness K.GP ---
        // Fnm: matrix of spins of nodal forces (translational only)
        // Assemble spin matrices at translational force rows only
        // (rows 0, 3, 6, 9, 12, 15, 18, 21 -> one per node, 3 DOF each)
        static MatrixType Fnm(24, 3);
        Fnm.Zero();
        for (int i = 0; i < 8; i++)
            ASDEICR3::Spin_AtRow(projectedLocalForces, Fnm, i * 3);

        static MatrixType FnmT(3, 24);
        FnmT.addMatrixTranspose(0.0, Fnm, 1.0);


        // SIGNS. The incoming RHS of this element is the POSITIVE internal force
        // (calculateAll accumulates +B^T sigma and getResistingForce returns it
        // unchanged, which is what OpenSees expects), so projectedLocalForces is
        // +P^T p_L and both geometric terms enter with a MINUS:
        //
        //     K_local = P^T K_L P  -  Fnm G  -  G^T Fnm^T P
        //
        // Do not copy the '+' from ASDShellQ4CorotationalTransformation: the
        // comment there claims the incoming RHS is already negated, which is not
        // what its own calculateAll() produces. It has never mattered for the
        // shell because that class fixes G = 0, so both terms vanish identically.

        // --- Step 2: Equilibrium projection geometric stiffness K.GP ---
        // -G^T * Fn^T * P
        temp.addMatrixTransposeProduct(0.0, G, FnmT, 1.0);
        LHS.addMatrixProduct(1.0, temp, P, -1.0);

        // --- Step 3: Rotational geometric stiffness K.GR: -Fnm * G ---
        // Counterpart of Step 3 of the shell transformation. There, Fnm is first
        // extended with the spins of the nodal MOMENTS; here the nodes carry no
        // moments, so Fnm holds the spins of the nodal forces only and is used as
        // it stands. The term was absent while G was zero (it contributed
        // nothing); it is required for the tangent to be consistent.
        LHS.addMatrixProduct(1.0, Fnm, G, -1.0);

        // --- Step 4: Transform to global ---
        // K_global = T^T * K_local * T
        temp.addMatrixProduct(0.0, LHS, T, 1.0);
        LHS.addMatrixTransposeProduct(0.0, T, temp, 1.0);
    }

    // -----------------------------------------------------------------------
    // transformToGlobal  (without explicit globalDisplacements)
    // Convenience overload: fetches global displacements and local
    // displacements internally before calling the main overload.
    // -----------------------------------------------------------------------
    void transformToGlobal(
        const ASDSolidHexLocalCoordinateSystem& LCS,
        MatrixType& LHS,
        VectorType& RHS,
        bool LHSrequired)
    {
        static VectorType globalDisplacements(24);
        static VectorType localDisplacements(24);
        computeGlobalDisplacements(globalDisplacements);
        calculateLocalDisplacements(LCS, globalDisplacements, localDisplacements);
        transformToGlobal(LCS, globalDisplacements, localDisplacements,
            LHS, RHS, LHSrequired);
    }

    // -----------------------------------------------------------------------
    // internalDataSize
    // How many doubles are needed to save/restore the internal state.
    //
    // Shell Q4 needed: 24 (U0) + 4*4 (Q0 + 4*QN) + 4*4 (4*QN_conv)
    //                + 3 (C0) + 4*3 (RV) + 4*3 (RV_conv) = 87
    //
    // Hex8 needs:     24 (U0) + 4 (Q0) + 3 (C0) = 31
    // (No nodal quaternions, no rotation vectors)
    // -----------------------------------------------------------------------
    int internalDataSize() const
    {
        // 24: initial displacements m_U0
        //  4: reference frame quaternion m_Q0 (w, x, y, z)
        //  3: reference centroid m_C0 (x, y, z)
        return 31;
    }

    // -----------------------------------------------------------------------
    // saveInternalData
    // -----------------------------------------------------------------------
    void saveInternalData(VectorType& v, int pos) const
    {
        if ((v.Size() - pos) < internalDataSize()) {
            opserr << "ASDSolidHexCorotationalTransformation - "
                << "failed to save internal data: vector too small\n";
            exit(-1);
        }

        // 24: initial displacements
        for (int i = 0; i < 24; i++)
            v(pos++) = m_U0(i);

        // 4: reference frame quaternion
        v(pos++) = m_Q0.w();
        v(pos++) = m_Q0.x();
        v(pos++) = m_Q0.y();
        v(pos++) = m_Q0.z();

        // 3: reference centroid
        v(pos++) = m_C0(0);
        v(pos++) = m_C0(1);
        v(pos++) = m_C0(2);
    }

    // -----------------------------------------------------------------------
    // restoreInternalData
    // -----------------------------------------------------------------------
    void restoreInternalData(const VectorType& v, int pos)
    {
        if ((v.Size() - pos) < internalDataSize()) {
            opserr << "ASDSolidHexCorotationalTransformation - "
                << "failed to restore internal data: vector too small\n";
            exit(-1);
        }

        // 24: initial displacements
        for (int i = 0; i < 24; i++)
            m_U0(i) = v(pos++);

        // 4: reference frame quaternion
        m_Q0 = QuaternionType(v(pos), v(pos + 1), v(pos + 2), v(pos + 3));
        pos += 4;

        // 3: reference centroid
        m_C0 = Vector3Type(v(pos), v(pos + 1), v(pos + 2));
        pos += 3;
    }

    void computeGlobalDisplacements(VectorType& globalDisplacements) const {

        for (int i = 0; i < 8; i++) {
            int index = i * 3;
            const VectorType& iU = m_nodes[i]->getTrialDisp();
            for (int j = 0; j < 3; j++) {
                globalDisplacements(index + j) = iU(j) - m_U0(index + j);
            }
        }

    }

    const MatrixType& computeTransformationMatrix(const ASDSolidHexLocalCoordinateSystem& LCS) const
    {

        static MatrixType T(24, 24);

            LCS.ComputeTotalRotationMatrix(T);

        return T;
    }


    //void calculateLocalDisplacements(
    //    const ASDSolidHexLocalCoordinateSystem& LCS,
    //    const VectorType& globalDisplacements,
    //    VectorType& localDisplacements)
    //{
    //    const MatrixType& R = computeTransformationMatrix(LCS);
    //    localDisplacements.addMatrixVector(0.0, R, globalDisplacements, 1.0);
    //}

private:

    // -----------------------------------------------------------------------
    // computeRtilde
    // Builds the orthonormal local frame [e1 | e2 | e3] from the 8 nodal
    // positions of a Hex8 element (as rows of a 3x3 rotation matrix).
    //
    // e1: mean tangent in the xi  direction (faces xi=-1 vs xi=+1)
    // e2: Gram-Schmidt orthogonalization of eta tangent w.r.t. e1
    // e3: e1 x e2_tmp  (right-hand normal)
    //
    // Node numbering (0-based):
    //   xi=-1 face: nodes 0,3,4,7
    //   xi=+1 face: nodes 1,2,5,6
    //  eta=-1 face: nodes 0,1,4,5
    //  eta=+1 face: nodes 2,3,6,7
    // -----------------------------------------------------------------------
    MatrixType computeRtilde(const std::array<Vector3Type, 8>& P) const
    {
        // e1: mean tangent along xi
        // = mean(xi=+1 face) - mean(xi=-1 face)
        Vector3Type e1 = (P[1] + P[2] + P[5] + P[6]) / 4.0
            - (P[0] + P[3] + P[4] + P[7]) / 4.0;
        e1.normalize();

        // e2_tmp: mean tangent along eta
        // = mean(eta=+1 face) - mean(eta=-1 face)
        Vector3Type e2t = (P[2] + P[3] + P[6] + P[7]) / 4.0
            - (P[0] + P[1] + P[4] + P[5]) / 4.0;

        // e3 = e1 x e2_tmp  (normal, not yet in the e1-e2 plane)
        Vector3Type e3 = e1.cross(e2t);
        e3.normalize();

        // e2 = e3 x e1  (Gram-Schmidt: guarantees orthonormality)
        Vector3Type e2 = e3.cross(e1);
        e2.normalize();

        // assemble rotation matrix: rows = local axes
        MatrixType R(3, 3);
        R(0, 0) = e1(0);  R(0, 1) = e1(1);  R(0, 2) = e1(2);
        R(1, 0) = e2(0);  R(1, 1) = e2(1);  R(1, 2) = e2(2);
        R(2, 0) = e3(0);  R(2, 1) = e3(1);  R(2, 2) = e3(2);
        return R;
    }

    // -----------------------------------------------------------------------
    // computeCentreGradients
    // g_a = dN_a/dX at the parametric centre xi = eta = zeta = 0, in reference
    // axes. A material (reference) object, so its components in the corotated
    // frame are the same numbers -- which is what makes RotationGradient() below
    // frame-consistent.
    //
    // Two identities used later, both exact by construction:
    //   sum_a g_a         = 0     (partition of unity  ->  G kills translations)
    //   sum_a X_a (x) g_a = I     (->  G returns exactly theta for a rigid rotation)
    // -----------------------------------------------------------------------
    void computeCentreGradients()
    {
        static const double s_xi[8]   = { -1.0, +1.0, +1.0, -1.0, -1.0, +1.0, +1.0, -1.0 };
        static const double s_eta[8]  = { -1.0, -1.0, +1.0, +1.0, -1.0, -1.0, +1.0, +1.0 };
        static const double s_zeta[8] = { -1.0, -1.0, -1.0, -1.0, +1.0, +1.0, +1.0, +1.0 };

        // dN_a/dxi at the centre
        double dN[8][3];
        for (int a = 0; a < 8; ++a) {
            dN[a][0] = 0.125 * s_xi[a];
            dN[a][1] = 0.125 * s_eta[a];
            dN[a][2] = 0.125 * s_zeta[a];
        }

        // J0(i,k) = dX_i/dxi_k = sum_a X_a(i) dN_a/dxi_k
        double J0[3][3] = { {0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}, {0.0, 0.0, 0.0} };
        for (int a = 0; a < 8; ++a) {
            const Vector3Type Xa(m_nodes[a]->getCrds());
            for (int i = 0; i < 3; ++i)
                for (int k = 0; k < 3; ++k)
                    J0[i][k] += Xa(i) * dN[a][k];
        }

        const double detJ0 =
            J0[0][0] * (J0[1][1] * J0[2][2] - J0[1][2] * J0[2][1]) -
            J0[0][1] * (J0[1][0] * J0[2][2] - J0[1][2] * J0[2][0]) +
            J0[0][2] * (J0[1][0] * J0[2][1] - J0[1][1] * J0[2][0]);
        if (std::abs(detJ0) < 1.0e-300) {
            // fully degenerate reference element: leave the gradients at zero, which
            // makes G = 0 and recovers the previous (inconsistent but harmless)
            // tangent instead of producing inf/nan.
            for (int a = 0; a < 8; ++a)
                m_g[a] = Vector3Type(0.0, 0.0, 0.0);
            return;
        }
        const double inv = 1.0 / detJ0;
        double iJ0[3][3];
        iJ0[0][0] = (J0[1][1] * J0[2][2] - J0[1][2] * J0[2][1]) * inv;
        iJ0[0][1] = (J0[0][2] * J0[2][1] - J0[0][1] * J0[2][2]) * inv;
        iJ0[0][2] = (J0[0][1] * J0[1][2] - J0[0][2] * J0[1][1]) * inv;
        iJ0[1][0] = (J0[1][2] * J0[2][0] - J0[1][0] * J0[2][2]) * inv;
        iJ0[1][1] = (J0[0][0] * J0[2][2] - J0[0][2] * J0[2][0]) * inv;
        iJ0[1][2] = (J0[0][2] * J0[1][0] - J0[0][0] * J0[1][2]) * inv;
        iJ0[2][0] = (J0[1][0] * J0[2][1] - J0[1][1] * J0[2][0]) * inv;
        iJ0[2][1] = (J0[0][1] * J0[2][0] - J0[0][0] * J0[2][1]) * inv;
        iJ0[2][2] = (J0[0][0] * J0[1][1] - J0[0][1] * J0[1][0]) * inv;

        // g_a(j) = sum_k dN_a/dxi_k * iJ0(k,j)
        for (int a = 0; a < 8; ++a) {
            for (int j = 0; j < 3; ++j) {
                double s = 0.0;
                for (int k = 0; k < 3; ++k)
                    s += dN[a][k] * iJ0[k][j];
                m_g[a](j) = s;
            }
        }
    }

    // -----------------------------------------------------------------------
    // RotationGradient
    // Computes the rotation-gradient (spin-lever) matrix G (3 x 24):
    //
    //     dw = G * dv
    //
    // with dw the infinitesimal spin of the corotational frame, defined by
    // R^T dR = skew(dw), and dv the nodal displacement variations expressed in
    // the corotated frame. G is therefore fixed by the frame-fitting rule, it is
    // not a modelling choice.
    //
    // This element fits the frame by the polar decomposition of the deformation
    // gradient at the parametric centre (see createLocalCoordinateSystem):
    //
    //     F = sum_a x_a (x) g_a ,     F = R U ,     R = polar(F)
    //
    // F is LINEAR in the nodal positions, which is what makes a closed form
    // possible. Differentiating F = R U and using that dU is symmetric and
    // W = R^T dR is skew:
    //
    //     R^T dF = W U + dU
    //     R^T dF - (R^T dF)^T = W U + U W
    //     axial(W U + U W) = (tr(U) I - U) w              (identity, 3x3)
    //     axial(R^T dF - ...^T) = sum_a g_a x dv_a
    //
    // hence, block by block,
    //
    //     G_a = (tr(U) I - U)^{-1} skew(g_a)
    //
    // M = tr(U) I - U is symmetric positive definite for any admissible stretch:
    // its eigenvalues are lam_j + lam_k > 0. So it is always invertible, and the
    // only way to lose it is a fully collapsed element.
    //
    // Verified against central finite differences of polar(F) on a regular and a
    // strongly distorted hexahedron, at small strain, at 63 deg + strain and at
    // 137 deg + 60% stretch: agreement to 1.4e-8 relative, which is the finite
    // difference noise floor. The two structural identities
    //     G * (rigid translation) = 0        and     G * (rigid rotation) = theta
    // hold to 6e-17.
    // -----------------------------------------------------------------------
    inline void RotationGradient(
        const ASDSolidHexLocalCoordinateSystem& LCS,
        const VectorType& globalDisplacements,
        MatrixType& G)
    {
        G.Zero();

        // Rt is the global -> corotated rotator, i.e. R^T (see the
        // USE_SIMPLE_RFRAME branch of ASDSolidHexLocalCoordinateSystem)
        const MatrixType& Rt = LCS.getRotationMatrix();

        // F = sum_a x_a (x) g_a, with x_a the current global position of node a
        // and the initial displacements removed, exactly as in
        // createLocalCoordinateSystem()
        double F[3][3] = { {0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}, {0.0, 0.0, 0.0} };
        for (int a = 0; a < 8; ++a) {
            const int index = a * 3;
            Vector3Type xa(m_nodes[a]->getCrds());
            for (int i = 0; i < 3; ++i)
                xa(i) += globalDisplacements(index + i) - m_U0(index + i);
            const Vector3Type& ga = m_g[a];
            for (int i = 0; i < 3; ++i)
                for (int j = 0; j < 3; ++j)
                    F[i][j] += xa(i) * ga(j);
        }

        // U = R^T F = Rt * F  (the right stretch: symmetric up to round-off)
        double U[3][3];
        for (int i = 0; i < 3; ++i) {
            for (int j = 0; j < 3; ++j) {
                double s = 0.0;
                for (int k = 0; k < 3; ++k)
                    s += Rt(i, k) * F[k][j];
                U[i][j] = s;
            }
        }

        // M = tr(U) I - U, then M^{-1}
        const double trU = U[0][0] + U[1][1] + U[2][2];
        double M[3][3];
        for (int i = 0; i < 3; ++i)
            for (int j = 0; j < 3; ++j)
                M[i][j] = (i == j ? trU : 0.0) - U[i][j];

        const double detM =
            M[0][0] * (M[1][1] * M[2][2] - M[1][2] * M[2][1]) -
            M[0][1] * (M[1][0] * M[2][2] - M[1][2] * M[2][0]) +
            M[0][2] * (M[1][0] * M[2][1] - M[1][1] * M[2][0]);
        if (std::abs(detM) < 1.0e-300)
            return;  // collapsed element: fall back to G = 0
        const double invM = 1.0 / detM;
        double Mi[3][3];
        Mi[0][0] = (M[1][1] * M[2][2] - M[1][2] * M[2][1]) * invM;
        Mi[0][1] = (M[0][2] * M[2][1] - M[0][1] * M[2][2]) * invM;
        Mi[0][2] = (M[0][1] * M[1][2] - M[0][2] * M[1][1]) * invM;
        Mi[1][0] = (M[1][2] * M[2][0] - M[1][0] * M[2][2]) * invM;
        Mi[1][1] = (M[0][0] * M[2][2] - M[0][2] * M[2][0]) * invM;
        Mi[1][2] = (M[0][2] * M[1][0] - M[0][0] * M[1][2]) * invM;
        Mi[2][0] = (M[1][0] * M[2][1] - M[1][1] * M[2][0]) * invM;
        Mi[2][1] = (M[0][1] * M[2][0] - M[0][0] * M[2][1]) * invM;
        Mi[2][2] = (M[0][0] * M[1][1] - M[0][1] * M[1][0]) * invM;

        // G_a = M^{-1} skew(g_a)
        for (int a = 0; a < 8; ++a) {
            const Vector3Type& g = m_g[a];
            // skew(g) = [ 0 -g3 g2 ; g3 0 -g1 ; -g2 g1 0 ]
            const double sk[3][3] = { {   0.0, -g(2),  g(1) },
                                      {  g(2),   0.0, -g(0) },
                                      { -g(1),  g(0),   0.0 } };
            const int col = a * 3;
            for (int i = 0; i < 3; ++i)
                for (int j = 0; j < 3; ++j) {
                    double s = 0.0;
                    for (int k = 0; k < 3; ++k)
                        s += Mi[i][k] * sk[k][j];
                    G(i, col + j) = s;
                }
        }
    }

    // -----------------------------------------------------------------------
// Helper: build the 3x3 skew-symmetric spin matrix from a 3-vector
//   Spin(x) * y = x cross y
// -----------------------------------------------------------------------
    inline void Spin(const Vector3Type& x, MatrixType& S3x3)
    {
        S3x3(0, 0) = 0.0;       S3x3(0, 1) = -x(2);    S3x3(0, 2) = x(1);
        S3x3(1, 0) = x(2);      S3x3(1, 1) = 0.0;     S3x3(1, 2) = -x(0);
        S3x3(2, 0) = -x(1);      S3x3(2, 1) = x(0);    S3x3(2, 2) = 0.0;
    }

    // -----------------------------------------------------------------------
    // SpinLeverMatrix
    // Computes the 24 x 3 spin-lever matrix S in the CR frame.
    // S connects the rigid body spin dw to the nodal translations:
    //     du_a^{rig,rot} = -Spin(x_bar_a) * dw
    // so that the a-th 3x3 block of S is  S_a = -Spin(x_bar_a).
    //
    // Coordinates x_bar_a are taken from the deformed configuration
    // in the CR frame (already stored in LCS.m_P, which is centroid-
    // referred and CR-rotated).
    // -----------------------------------------------------------------------
    inline void SpinLeverMatrix(
        const ASDSolidHexLocalCoordinateSystem& LCS,
        MatrixType& S)
    {
        // S is 24 x 3
        S.Zero();
        for (int a = 0; a < 8; ++a) {
            const Vector3Type& xa = LCS.P(a);  // coords in CR frame, centroid-referred
            const int row = 3 * a;

            // S_a = -Spin(xa)
            //   = [  0    xa_z  -xa_y ]
            //     [-xa_z    0    xa_x ]
            //     [ xa_y -xa_x    0   ]
            S(row + 0, 0) = 0.0;       S(row + 0, 1) = xa(2);    S(row + 0, 2) = -xa(1);
            S(row + 1, 0) = -xa(2);     S(row + 1, 1) = 0.0;      S(row + 1, 2) = xa(0);
            S(row + 2, 0) = xa(1);     S(row + 2, 1) = -xa(0);    S(row + 2, 2) = 0.0;
        }
    }

        inline const NodeContainerType& getNodes()const { return m_nodes; }
        inline NodeContainerType& getNodes() { return m_nodes; }

private:

    // Reference frame orientation (stored as quaternion, consistent
    // with the shell implementation)
    QuaternionType m_Q0;

    // Reference centroid C0 (undeformed configuration)
    Vector3Type m_C0;

    // g_a = dN_a/dX at the parametric centre, in reference axes. Constant for the
    // life of the element (it only involves the reference coordinates), so it is
    // computed once in setDomain() and never serialized. Used by RotationGradient()
    // and, through it, by the consistent corotational tangent.
    std::array<Vector3Type, 8> m_g;

    NodeContainerType m_nodes = { {nullptr, nullptr, nullptr, nullptr, nullptr, nullptr, nullptr, nullptr} };
    Vector m_U0 = Vector(24);
};

#endif // !ASDSolidHexCorotationalTransformation_h