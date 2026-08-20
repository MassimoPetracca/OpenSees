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
// $Date: 2024/03

// Original implementation: Massimo Petracca (ASDEA)
//
// Implementation of a corotational coordinate transformation 3-node shells
//

#ifndef ASDShellT3CorotationalTransformation_h
#define ASDShellT3CorotationalTransformation_h

#include <ASDEICR.h>
#include <ASDShellT3Transformation.h>

// this is experimental: it fits the corotational frame following the polar
// decomposition rather than the 1-2 side alignment as per Felippa's work
#define USE_POLAR_DECOMP_ALLIGN

// See ASDShellQ4CorotationalTransformation.h for what these select.
// 0 = central finite differences (18 frame constructions per tangent for the T3),
// 1 = closed form. For a triangle the closed form is exact and free of the two
// awkward terms of the quadrilateral: the three nodes define the plane, so their
// local z vanishes, and side 1-2 lies in the plane, so it needs no projection.
#ifndef ASDSHELL_ANALYTIC_ROTATION_GRADIENT
#define ASDSHELL_ANALYTIC_ROTATION_GRADIENT 1
#endif
#ifndef ASDSHELL_G_VERIFY
#define ASDSHELL_G_VERIFY 0
#endif

/** \brief ASDShellT3CorotationalTransformation
*
* This class represents a corotational (nonlinear) coordinate transformation
* that can be used by any element whose geometry is a Triangle 3 in 3D space,
* with 6 D.O.F.s per node.
* Its main aim is to:
* 1) Create the local coordinate system
* 2) Transform the incoming global displacements in local coordinate system
*    removing rigid body displacements and rotations.
* 3) Transform the outgoing matrices and vectors in global coordinate system
*    with rigid body displacements and rotations.
*
* - Makes use of Quaternions (Euler Parameters) to parametrize finite rotations
*   in an efficient and robust way.
*
* References:
*
* - C.A.Felippa,B.Haugen, "Unified formulation of small-strain corotational
*   finite elements: I. Theory",
*   CU-CAS-05-02, January 2005
* - C.A.Felippa, AFEM.Ch.38, "Quadrilateral Shell Elements",
*   Chapter 5 of B.Haugen's Thesis.
*   link: http://www.colorado.edu/engineering/CAS/courses.d/AFEM.d/
*/
class ASDShellT3CorotationalTransformation : public ASDShellT3Transformation
{

public:

    typedef ASDVector3<double> Vector3Type;

    typedef ASDQuaternion<double> QuaternionType;

    typedef Vector VectorType;

    typedef Matrix MatrixType;

    typedef std::array<Node*, 3> NodeContainerType;

public:

    ASDShellT3CorotationalTransformation()
        : ASDShellT3Transformation()
    {
    }

    virtual ~ASDShellT3CorotationalTransformation()
    {
    }

public:

    virtual ASDShellT3Transformation* create()const
    {
        return new ASDShellT3CorotationalTransformation();
    }

    virtual bool isLinear() const
    {
        return false;
    }

    virtual void revertToStart()
    {
        // create the reference (undeformed configuration) coordinate system
        ASDShellT3LocalCoordinateSystem LCS = createReferenceCoordinateSystem();

        // save reference orientation and center
        m_Q0 = QuaternionType::FromRotationMatrix(LCS.Orientation());
        m_C0 = LCS.Center();

        // save initial rotations, no need to take current rotation
        // since we will remove the initial ones (themselves)...
        for (int i = 0; i < 3; i++)
        {
            m_RV[i] = Vector3Type(0.0, 0.0, 0.0);
            m_QN[i] = QuaternionType::FromRotationVector(m_RV[i]);

            m_RV_converged[i] = m_RV[i];
            m_QN_converged[i] = m_QN[i];
        }
    }

    virtual void setDomain(Domain* domain, const ID& node_ids, bool initialized)
    {
        // call base class setDomain to
        // get nodes and save initial displacements and rotations
        ASDShellT3Transformation::setDomain(domain, node_ids, initialized);

        // quick return
        if (domain == nullptr || initialized)
            return;

        // init state variables
        revertToStart();
    }

    virtual void revertToLastCommit()
    {
        for (int i = 0; i < 3; i++)
        {
            m_RV[i] = m_RV_converged[i];
            m_QN[i] = m_QN_converged[i];
        }
    }

    virtual void commit()
    {
        for (int i = 0; i < 3; i++)
        {
            m_RV_converged[i] = m_RV[i];
            m_QN_converged[i] = m_QN[i];
        }
    }

    virtual void update(const VectorType& globalDisplacements)
    {
        for (int i = 0; i < 3; i++)
        {
            // compute current rotation vector removing initial rotations if any
            Vector3Type currentRotVec;
            int index = i * 6;
            currentRotVec(0) = globalDisplacements(index + 3) - m_U0(index + 3);
            currentRotVec(1) = globalDisplacements(index + 4) - m_U0(index + 4);
            currentRotVec(2) = globalDisplacements(index + 5) - m_U0(index + 5);

            // compute incremental rotation vector
            Vector3Type incrementalRotation = currentRotVec - m_RV[i];

            // save current rotation vector
            m_RV[i] = currentRotVec;

            // compute incremental quaternion from incremental rotation vector
            QuaternionType incrementalQuaternion = QuaternionType::FromRotationVector(incrementalRotation);

            // update nodal quaternion
            m_QN[i] = incrementalQuaternion * m_QN[i];
        }
    }

    virtual ASDShellT3LocalCoordinateSystem createLocalCoordinateSystem(const VectorType& globalDisplacements)const
    {
        // reference coordinate system
        ASDShellT3LocalCoordinateSystem a = createReferenceCoordinateSystem();

        // compute nodal positions at current configuration removing initial displacements if any
        std::array<Vector3Type, 3> def = {
            Vector3Type(m_nodes[0]->getCrds()),
            Vector3Type(m_nodes[1]->getCrds()),
            Vector3Type(m_nodes[2]->getCrds())
        };
        for (int i = 0; i < 3; i++) {
            int index = i * 6;
            Vector3Type& iP = def[i];
            iP(0) += globalDisplacements(index) - m_U0(index);
            iP(1) += globalDisplacements(index + 1) - m_U0(index + 1);
            iP(2) += globalDisplacements(index + 2) - m_U0(index + 2);
        }

        // current coordinate system
        ASDShellT3LocalCoordinateSystem b(def[0], def[1], def[2]);

#ifndef USE_POLAR_DECOMP_ALLIGN
        return b;
#endif // !USE_POLAR_DECOMP_ALLIGN

        double aX1 = a.X1(); double aY1 = a.Y1();
        double bX1 = b.X1(); double bY1 = b.Y1();
        double aX2 = a.X2(); double aY2 = a.Y2();
        double bX2 = b.X2(); double bY2 = b.Y2();
        double aX3 = a.X3(); double aY3 = a.Y3();
        double bX3 = b.X3(); double bY3 = b.Y3();

        // now we are in the local coordinate systems (reference and current), i.e. we are looking in the local Z direction
        // which is the same for both coordinate systems.
        // now we can compute the 2D deformation gradient between the 2 configurations, at the element center.

        double C1 = 1.0 / (aX1*aY2 - aX1*aY3 - aX2*aY1 + aX2*aY3 + aX3*aY1 - aX3*aY2);
        double f11 = C1*(-(aY1 - aY2)*(bX1 - bX3) + (aY1 - aY3)*(bX1 - bX2));
        double f12 = C1*((aX1 - aX2)*(bX1 - bX3) - (aX1 - aX3)*(bX1 - bX2));
        double f21 = C1*(-(aY1 - aY2)*(bY1 - bY3) + (aY1 - aY3)*(bY1 - bY2));
        double f22 = C1*((aX1 - aX2)*(bY1 - bY3) - (aX1 - aX3)*(bY1 - bY2));

        // now we can extrapolate the rotation angle that makes this deformation gradient symmetric.
        // F = R*U -> find R such that R'*F = U
        double alpha = std::atan2(f21 - f12, f11 + f22);

        // this final coordinate system is the one in which 
        // the deformation gradient is equal to the stretch tensor
        return ASDShellT3LocalCoordinateSystem(def[0], def[1], def[2], alpha);
    }

    virtual void calculateLocalDisplacements(
        const ASDShellT3LocalCoordinateSystem& LCS,
        const VectorType& globalDisplacements,
        VectorType& localDisplacements)
    {
        // orientation and center of current local coordinate system
        QuaternionType Q = QuaternionType::FromRotationMatrix(LCS.Orientation());
        const Vector3Type& C = LCS.Center();

        for (int i = 0; i < 3; i++)
        {
            int index = i * 6;

            // centered undeformed position
            Vector3Type X0 = Vector3Type(m_nodes[i]->getCrds());
            X0 -= m_C0;

            // centered deformed position
            Vector3Type X = X0 + Vector3Type(globalDisplacements, index);
            X -= C;

            // get deformational displacements
            Q.rotateVector(X);
            m_Q0.rotateVector(X0);
            Vector3Type deformationalDisplacements = X - X0;

            localDisplacements[index] = deformationalDisplacements[0];
            localDisplacements[index + 1] = deformationalDisplacements[1];
            localDisplacements[index + 2] = deformationalDisplacements[2];

            // get deformational rotations
            QuaternionType Qd = Q * m_QN[i] * m_Q0.conjugate();
            Qd.toRotationVector(
                localDisplacements[index + 3],
                localDisplacements[index + 4],
                localDisplacements[index + 5]);
        }
    }

    virtual void transformToGlobal(
        const ASDShellT3LocalCoordinateSystem& LCS,
        const VectorType& globalDisplacements,
        const VectorType& localDisplacements,
        MatrixType& LHS,
        VectorType& RHS,
        bool LHSrequired)
    {
        // Get the total rotation matrix (local - to - global)

        static MatrixType T(18, 18);
        LCS.ComputeTotalRotationMatrix(T);

        // Form all matrices:
        // S: Spin-Fitter matrix
        // G: Spin-Lever matrix
        // P: Projector (Translational & Rotational)
        static MatrixType P(18, 18);
        static MatrixType S(18, 3);
        static MatrixType G(3, 18);
        EICR::Compute_Pt(3, P);
        EICR::Compute_S(LCS.Nodes(), S);
        RotationGradient(LCS, globalDisplacements, G);
        P.addMatrixProduct(1.0, S, G, -1.0); // P -= S*G

        // Compute the projected local forces ( pe = P' * RHS ).
        // Note: here the RHS is already given as a residual vector -> - internalForces -> (pe = - Ke * U)
        // so projectedLocalForces = - P' * Ke * U

        static VectorType projectedLocalForces(18);
        projectedLocalForces.addMatrixTransposeVector(0.0, P, RHS, 1.0);

        // Compute the Right-Hand-Side vector in global coordinate system (- T' * P' * Km * U).
        // At this point the computation of the Right-Hand-Side is complete.

        RHS.addMatrixTransposeVector(0.0, T, projectedLocalForces, 1.0);

        // Begin the computation of the Left-Hand-Side Matrix :

        if (!LHSrequired) 
            return; // avoid useless calculations!

        // H: Axial Vector Jacobian
        static MatrixType H(18, 18);
        EICR::Compute_H(localDisplacements, H);

        // Step 1: ( K.M : Material Stiffness Matrix )
        // Apply the projector to the Material Stiffness Matrix (Ke = P' * Km * H * P)
        // At this point 'LHS' contains the 'projected' Material Stiffness matrix
        // in local corotational coordinate system

        static MatrixType temp(18, 18);
        temp.addMatrixProduct(0.0, LHS, H, 1.0);
        LHS.addMatrixProduct(0.0, temp, P, 1.0);
        temp.addMatrixTransposeProduct(0.0, P, LHS, 1.0);
        LHS = temp;

        // Step 2: ( K.GP: Equilibrium Projection Geometric Stiffness Matrix )
        // First assemble the 'Fnm' matrix with the Spins of the nodal forces.
        // Actually at this point the 'Fnm' Matrix is the 'Fn' Matrix,
        // because it only contains the spins of the 'translational' forces.
        // At this point 'LHS' contains also this term of the Geometric stiffness
        // (Ke = (P' * Km * H * P) - (G' * Fn' * P))

        static MatrixType Fnm(18, 3);
        Fnm.Zero();
        EICR::Spin_AtRow(projectedLocalForces, Fnm, 0);
        EICR::Spin_AtRow(projectedLocalForces, Fnm, 6);
        EICR::Spin_AtRow(projectedLocalForces, Fnm, 12);

        static MatrixType FnmT(3, 18);
        FnmT.addMatrixTranspose(0.0, Fnm, 1.0);

        temp.addMatrixTransposeProduct(0.0, G, FnmT, 1.0);
        // SIGN: minus, same as ASDShellQ4. ASDShellT3::calculateAll accumulates
        // +B'*S*dA, i.e. the POSITIVE internal force, so the geometric terms take
        // a minus. It never mattered while G was fixed to zero.
        LHS.addMatrixProduct(1.0, temp, P, -1.0); // -G' * Fn' * P

        // Step 3: ( K.GR: Rotational Geometric Stiffness Matrix )
        // Add the Spins of the nodal moments to 'Fnm'.
        // At this point 'LHS' contains also this term of the Geometric stiffness
        // (Ke = (P' * Km * H * P) - (G' * Fn' * P) - (Fnm * G))

        EICR::Spin_AtRow(projectedLocalForces, Fnm, 3);
        EICR::Spin_AtRow(projectedLocalForces, Fnm, 9);
        EICR::Spin_AtRow(projectedLocalForces, Fnm, 15);

        LHS.addMatrixProduct(1.0, Fnm, G, -1.0); // -Fnm * G

        // Step 4: (Global Stiffness Matrix)
        // Transform the LHS to the Global coordinate system.
        // T' * [(P' * Km * H * P) - (G' * Fn' * P) - (Fnm * G)] * T
        temp.addMatrixProduct(0.0, LHS, T, 1.0);
        LHS.addMatrixTransposeProduct(0.0, T, temp, 1.0);
    }

    virtual void transformToGlobal(
        const ASDShellT3LocalCoordinateSystem& LCS,
        MatrixType& LHS,
        VectorType& RHS,
        bool LHSrequired)
    {
        static VectorType globalDisplacements(18);
        static VectorType localDisplacements(18);
        computeGlobalDisplacements(globalDisplacements);
        calculateLocalDisplacements(LCS, globalDisplacements, localDisplacements);
        transformToGlobal(LCS, globalDisplacements, localDisplacements, LHS, RHS, LHSrequired);
    }

    virtual int internalDataSize() const
    {
        // 18 -> initial displacement +
        // 7*4 -> 7 quaternions +
        // 7*3 -> 7 3d vectors
        return 67;
    }

    virtual void saveInternalData(VectorType& v, int pos) const
    {
        if ((v.Size() - pos) < internalDataSize()) {
            opserr << "ASDShellT3CorotationalTransformation - failed to save internal data: vector too small\n";
            exit(-1);
        }

        // 18 -> initial displacement +
        for (int i = 0; i < 18; i++)
            v(pos++) = m_U0(i);

        // 7*4 -> 7 quaternions +
        auto lamq = [&v, &pos](const QuaternionType& x) {
            v(pos++) = x.w();
            v(pos++) = x.x();
            v(pos++) = x.y();
            v(pos++) = x.z();
        };
        lamq(m_Q0);
        for (int i = 0; i < 3; i++)
            lamq(m_QN[i]);
        for (int i = 0; i < 3; i++)
            lamq(m_QN_converged[i]);

        // 7*3 -> 7 3d vectors +
        auto lamv = [&v, &pos](const Vector3Type& x) {
            v(pos++) = x.x();
            v(pos++) = x.y();
            v(pos++) = x.z();
        };
        lamv(m_C0);
        for (int i = 0; i < 3; i++)
            lamv(m_RV[i]);
        for (int i = 0; i < 3; i++)
            lamv(m_RV_converged[i]);
    }

    virtual void restoreInternalData(const VectorType& v, int pos)
    {
        if ((v.Size() - pos) < internalDataSize()) {
            opserr << "ASDShellT3CorotationalTransformation - failed to restore internal data: vector too small\n";
            exit(-1);
        }
        
        // 18 -> initial displacement +
        for (int i = 0; i < 18; i++)
            m_U0(i) = v(pos++);

        // 7*4 -> 7 quaternions +
        auto lamq = [&v, &pos](QuaternionType& x) {
            x = QuaternionType(v(pos), v(pos+1), v(pos+2), v(pos+3));
            pos += 4;
        };
        lamq(m_Q0);
        for (int i = 0; i < 3; i++)
            lamq(m_QN[i]);
        for (int i = 0; i < 3; i++)
            lamq(m_QN_converged[i]);

        // 7*3 -> 7 3d vectors +
        // NOTE: assign the components one statement at a time.
        // Writing Vector3Type(v(pos++), v(pos++), v(pos++)) leaves the three
        // increments of 'pos' unsequenced: the order in which the constructor
        // arguments are evaluated is not specified by the standard, and the
        // build compiler (icx) evaluates them right-to-left, so every restored
        // 3d vector came back with its components REVERSED (x <-> z).
        auto lamv = [&v, &pos](Vector3Type& x) {
            x(0) = v(pos++);
            x(1) = v(pos++);
            x(2) = v(pos++);
        };
        lamv(m_C0);
        for (int i = 0; i < 3; i++)
            lamv(m_RV[i]);
        for (int i = 0; i < 3; i++)
            lamv(m_RV_converged[i]);
    }

private:

    /**
    * Computes the Spin Fitter Matrix, i.e. the rotation gradient G in dw = G dv,
    * with dw the spin of the corotational frame in LOCAL components.
    * This is the only matrix not included in the EICR, because it depends on how
    * the corotational frame follows the element.
    * Two implementations, selected by ASDSHELL_ANALYTIC_ROTATION_GRADIENT:
    * the closed form (default) and central finite differences (the reference).
    * @return the Spin Fitter Matrix
    */
    inline void RotationGradient(
        const ASDShellT3LocalCoordinateSystem& LCS,
        const VectorType& globalDisplacements,
        MatrixType& G)
    {
#if ASDSHELL_ANALYTIC_ROTATION_GRADIENT
        RotationGradient_Analytic(LCS, globalDisplacements, G);
#else
        RotationGradient_FD(LCS, globalDisplacements, G);
#endif

#if ASDSHELL_G_VERIFY
        static MatrixType Gother(3, 18);
#if ASDSHELL_ANALYTIC_ROTATION_GRADIENT
        RotationGradient_FD(LCS, globalDisplacements, Gother);
#else
        RotationGradient_Analytic(LCS, globalDisplacements, Gother);
#endif
        static double worstRel = -1.0;
        static double worstR2 = -1.0;
        double dmax = 0.0, gmax = 0.0;
        for (int i = 0; i < 3; i++) {
            for (int j = 0; j < 18; j++) {
                double d = std::abs(G(i, j) - Gother(i, j));
                if (d > dmax) dmax = d;
                double a = std::abs(Gother(i, j));
                if (a > gmax) gmax = a;
            }
        }
        double rel = (gmax > 0.0) ? dmax / gmax : dmax;
        double r2 = 0.0;
        for (int i = 0; i < 3; i++) {
            for (int j = 0; j < 3; j++) {
                double sum = 0.0;
                for (int a = 0; a < 3; a++) sum += G(i, a * 6 + j);
                if (std::abs(sum) > r2) r2 = std::abs(sum);
            }
        }
        // (R1): G S = I3, with S_a = -Spin(xbar_a) and xbar_a the deformed nodal
        // coordinates in the corotated frame. This is the requirement that makes
        // P = Pt - S G filter rigid rotations, and it holds for ANY correct G, so
        // it checks the derivative against the fitting rule itself rather than
        // against the finite differences.
        double r1 = 0.0;
        {
            const auto& xL = LCS.Nodes();
            double GS[3][3] = { {0.0,0.0,0.0},{0.0,0.0,0.0},{0.0,0.0,0.0} };
            for (int a = 0; a < 3; a++) {
                double x = xL[a][0], y = xL[a][1], z = xL[a][2];
                double Sa[3][3] = { {0.0,   z,  -y},
                                    { -z, 0.0,   x},
                                    {  y,  -x, 0.0} };
                for (int i = 0; i < 3; i++)
                    for (int j = 0; j < 3; j++)
                        for (int k = 0; k < 3; k++)
                            GS[i][j] += G(i, a * 6 + k) * Sa[k][j];
            }
            for (int i = 0; i < 3; i++)
                for (int j = 0; j < 3; j++) {
                    double d = std::abs(GS[i][j] - ((i == j) ? 1.0 : 0.0));
                    if (d > r1) r1 = d;
                }
        }
        static double worstR1 = -1.0;
        if (rel > worstRel || r2 > worstR2 || r1 > worstR1) {
            if (rel > worstRel) worstRel = rel;
            if (r2 > worstR2) worstR2 = r2;
            if (r1 > worstR1) worstR1 = r1;
            opserr << "ASDShellT3 G verify: |dG|/max|G| = " << worstRel
                << " , (R1) |G*S - I| = " << worstR1
                << " , (R2) residual = " << worstR2
                << " , max|G| = " << gmax << "\n";
        }
#endif
    }

    /**
    * Closed form of the rotation gradient for the triangle. See
    * SRC/element/ASDhexa/doc/claude_g_matrix.pdf section 6, in
    * particular 6.7: for a triangle the frame fitting simplifies twice, and both
    * simplifications are identities rather than approximations.
    *   - the three nodes define the plane exactly, so their local z vanishes and
    *     the warping coupling of the quadrilateral is absent;
    *   - e1 runs along side 1-2, which lies in that plane, so no projection is
    *     needed and the (v.e3) term of the quadrilateral vanishes.
    */
    inline void RotationGradient_Analytic(
        const ASDShellT3LocalCoordinateSystem& LCS,
        const VectorType& globalDisplacements,
        MatrixType& G)
    {
        G.Zero();

#ifndef USE_POLAR_DECOMP_ALLIGN
        RotationGradient_FD(LCS, globalDisplacements, G);
        return;
#else
        // deformed nodal positions, initial displacements removed
        std::array<Vector3Type, 3> p = {
            Vector3Type(m_nodes[0]->getCrds()),
            Vector3Type(m_nodes[1]->getCrds()),
            Vector3Type(m_nodes[2]->getCrds())
        };
        for (int i = 0; i < 3; i++) {
            int index = i * 6;
            for (int k = 0; k < 3; k++)
                p[i](k) += globalDisplacements(index + k) - m_U0(index + k);
        }

        // side-aligned current frame ("b"): alpha is defined from its local
        // coordinates. One construction, against 18 for the finite differences.
        ASDShellT3LocalCoordinateSystem b(p[0], p[1], p[2]);
        const MatrixType& Ob = b.Orientation();
        const auto& xb = b.Nodes();   // local coords, one Vector3 per node

        const MatrixType& Of = LCS.Orientation();
        Vector3Type e1(Of(0, 0), Of(0, 1), Of(0, 2));
        Vector3Type e2(Of(1, 0), Of(1, 1), Of(1, 2));
        Vector3Type e1t(Ob(0, 0), Ob(0, 1), Ob(0, 2));
        Vector3Type e2t(Ob(1, 0), Ob(1, 1), Ob(1, 2));

        Vector3Type r12 = p[1] - p[0];
        Vector3Type r13 = p[2] - p[0];
        double nn = 2.0 * b.Area();       // ||n|| = ||r12 x r13||
        double l12 = r12.norm();
        if (nn < 1.0e-300 || l12 < 1.0e-300)
            return;
        double inn = 1.0 / nn;

        // reference constants of the polar angle
        ASDShellT3LocalCoordinateSystem a = createReferenceCoordinateSystem();
        double aX1 = a.X1(); double aY1 = a.Y1();
        double aX2 = a.X2(); double aY2 = a.Y2();
        double aX3 = a.X3(); double aY3 = a.Y3();
        double C1 = 1.0 / (aX1 * aY2 - aX1 * aY3 - aX2 * aY1
            + aX2 * aY3 + aX3 * aY1 - aX3 * aY2);
        double k1 = -(aY1 - aY2);
        double k2 = (aY1 - aY3);
        double m1 = (aX1 - aX2);
        double m2 = -(aX1 - aX3);

        // current s = f21 - f12 and c = f11 + f22, linear in the local coords
        double bX1 = xb[0][0], bX2 = xb[1][0], bX3 = xb[2][0];
        double bY1 = xb[0][1], bY2 = xb[1][1], bY3 = xb[2][1];
        double ss = C1 * (k1 * (bY1 - bY3) + k2 * (bY1 - bY2)
            - m1 * (bX1 - bX3) - m2 * (bX1 - bX2));
        double cc = C1 * (k1 * (bX1 - bX3) + k2 * (bX1 - bX2)
            + m1 * (bY1 - bY3) + m2 * (bY1 - bY2));
        double den = ss * ss + cc * cc;
        if (den < 1.0e-300)
            return;
        double iden = 1.0 / den;

        for (int i = 0; i < 3; i++) {
            for (int j = 0; j < 3; j++) {

                Vector3Type du(Of(j, 0), Of(j, 1), Of(j, 2));

                // dn = C_i du, with C = (Spin(r13) - Spin(r12), -Spin(r13),
                // +Spin(r12))
                Vector3Type dn;
                switch (i) {
                case 0:  dn = r13.cross(du) - r12.cross(du); break;
                case 1:  dn = du.cross(r13); break;   // -(r13 x du)
                default: dn = r12.cross(du); break;
                }

                double dw0 = -e2.dot(dn) * inn;
                double dw1 = e1.dot(dn) * inn;
                double dw0t = -e2t.dot(dn) * inn;
                double dw1t = e1t.dot(dn) * inn;

                // side 1-2 is in the plane, so there is no (v.e3) term here
                double dvterm = 0.0;
                if (i == 0)      dvterm = -e2t.dot(du);
                else if (i == 1) dvterm = e2t.dot(du);
                double dw2t = dvterm / l12;

                // local coordinate variations. The local z of all three nodes is
                // identically zero, so the quadrilateral's warping terms
                // (-dw1t*z, +dw0t*z) are absent by construction.
                (void)dw0t; (void)dw1t;
                double dbX[3], dbY[3];
                for (int k = 0; k < 3; k++) {
                    double f = (k == i) ? (2.0 / 3.0) : (-1.0 / 3.0);
                    Vector3Type dk = f * du;
                    double ax = Ob(0, 0) * dk(0) + Ob(0, 1) * dk(1) + Ob(0, 2) * dk(2);
                    double ay = Ob(1, 0) * dk(0) + Ob(1, 1) * dk(1) + Ob(1, 2) * dk(2);
                    dbX[k] = ax + dw2t * xb[k][1];
                    dbY[k] = ay - dw2t * xb[k][0];
                }
                double ds = C1 * (k1 * (dbY[0] - dbY[2]) + k2 * (dbY[0] - dbY[1])
                    - m1 * (dbX[0] - dbX[2]) - m2 * (dbX[0] - dbX[1]));
                double dc = C1 * (k1 * (dbX[0] - dbX[2]) + k2 * (dbX[0] - dbX[1])
                    + m1 * (dbY[0] - dbY[2]) + m2 * (dbY[0] - dbY[1]));
                double dalpha = (cc * ds - ss * dc) * iden;

                int col = i * 6 + j;
                G(0, col) = dw0;
                G(1, col) = dw1;
                G(2, col) = dw2t + dalpha;
            }
        }
#endif // USE_POLAR_DECOMP_ALLIGN
    }

    /**
    * Central finite differences of the frame orientation: the reference
    * implementation, and the oracle for the closed form above.
    */
    inline void RotationGradient_FD(
        const ASDShellT3LocalCoordinateSystem& LCS,
        const VectorType& globalDisplacements,
        MatrixType& G)
    {
        G.Zero();

#ifdef USE_POLAR_DECOMP_ALLIGN

        /**
        Same change as ASDShellQ4CorotationalTransformation::RotationGradient --
        see the long note there. In short: the old claim that G "should be exactly
        0" for a polar-decomposition fitted frame cannot hold, because the frame
        must follow a rigid rotation of the element, which requires G * S = I3
        (equivalently P * S = 0). G = 0 leaves P = Pt and kills the geometric
        stiffness -Fnm * G - G' * Fnm' * P entirely.

        Computed by CENTRAL FINITE DIFFERENCES of the frame orientation. P, S and
        G all live in the LOCAL frame, so the perturbation is applied along the
        local axes, whose global components are the rows of LCS.Orientation().
        The frame depends on the nodal translations only, so the rotational
        columns stay zero. Spin extraction uses R0 * R1^T = I + Spin(w) with w in
        local components.
        */

        const MatrixType& R0 = LCS.Orientation();

        double pert = std::sqrt(createReferenceCoordinateSystem().Area()) * 1.0e-7;

        static VectorType UGP(18);
        UGP = globalDisplacements;

        double M[3][3];
        double wp[3], wm[3];

        for (int i = 0; i < 3; i++)
        {
            int index = i * 6;

            for (int j = 0; j < 3; j++)
            {
                for (int pass = 0; pass < 2; pass++)
                {
                    double h = (pass == 0) ? pert : -pert;
                    for (int c = 0; c < 3; c++)
                        UGP(index + c) = globalDisplacements(index + c) + h * R0(j, c);

                    // NOTE: must be a named local. A reference bound straight to
                    // createLocalCoordinateSystem(UGP).Orientation() dangles as
                    // soon as the temporary dies.
                    ASDShellT3LocalCoordinateSystem pertCS =
                        createLocalCoordinateSystem(UGP);
                    const MatrixType& R1 = pertCS.Orientation();

                    // M = R0 * R1^T
                    for (int r = 0; r < 3; r++) {
                        for (int c = 0; c < 3; c++) {
                            double sum = 0.0;
                            for (int k = 0; k < 3; k++)
                                sum += R0(r, k) * R1(c, k);
                            M[r][c] = sum;
                        }
                    }

                    double* w = (pass == 0) ? wp : wm;
                    w[0] = 0.5 * (M[2][1] - M[1][2]);
                    w[1] = 0.5 * (M[0][2] - M[2][0]);
                    w[2] = 0.5 * (M[1][0] - M[0][1]);
                }

                for (int c = 0; c < 3; c++)
                    UGP(index + c) = globalDisplacements(index + c);

                double den = 1.0 / (2.0 * pert);
                for (int r = 0; r < 3; r++)
                    G(r, index + j) = (wp[r] - wm[r]) * den;
            }
        }

#else // !USE_POLAR_DECOMP_ALLIGN

        /**
        Side-aligned frame: with USE_POLAR_DECOMP_ALLIGN undefined,
        createLocalCoordinateSystem returns b(def[0], def[1], def[2]) directly, so
        e1 runs along side 1-2 and the in-plane polar angle is absent.

        Rows 0 and 1 -- the tilt of e3. The three nodes define the plane exactly,
        so w(x,y) = sum_a N_a w_a is the CST interpolation of the local z
        displacements, and a rigid rotation theta gives w = theta_0 y - theta_1 x.
        Hence dw_0 = dw/dy and dw_1 = -dw/dx, i.e.

            G(0, w_a) = N_a,y = c_a / 2A ,   G(1, w_a) = -N_a,x = -b_a / 2A

        with the standard CST coefficients b_a = y_(a+1)(a+2), c_a = x_(a+2)(a+1):

            node 1 : c_1 = x32 , b_1 = y23   ->  x32/2A , y32/2A
            node 2 : c_2 = x13 , b_2 = y31   ->  x13/2A , y13/2A
            node 3 : c_3 = x21 , b_3 = y12   ->  x21/2A , y21/2A

        Row 2 -- the drilling rate. e1 = (p2-p1)/L12, so
        dw_2 = de1.e2 = e2.(du2-du1)/L12, which is +-1/L12 on the local y
        translations of nodes 1 and 2 and zero on node 3.

        Both structural requirements then hold identically:
          (R2) sum_a G_a = 0, because x32+x13+x21 = 0, y32+y13+y21 = 0 and
               -1/L12 + 1/L12 = 0;
          (R1) sum_a G_a S_a = I3 with S_a = -Spin(xbar_a). Entry (0,0) is
               sum_a N_a,y y_a = dy/dy = 1, entry (1,1) is sum_a N_a,x x_a = 1,
               entry (2,2) is (x2-x1)/L12 = 1, and every off-diagonal is a
               dy/dx-type derivative and vanishes.

        This is the same layout as the Q4 branch above. The two coefficients of
        node 3 used to repeat those of node 1, which broke (R2), and the drilling
        row used to carry h3 = 2A/L12 -- a LENGTH where a rotation gradient needs
        1/length, which made (R1) entry (2,2) equal 2A instead of 1.
        */

        const auto& P1 = LCS.P1();
        const auto& P2 = LCS.P2();
        const auto& P3 = LCS.P3();

        double Ap = 2.0 * LCS.Area();
        double m = 1.0 / Ap;

        double x13 = P1(0) - P3(0);
        double y13 = P1(1) - P3(1);
        double x32 = P3(0) - P2(0);
        double y32 = P3(1) - P2(1);
        double x21 = P2(0) - P1(0);
        double y21 = P2(1) - P1(1);

        double L3 = (P2 - P1).norm();
        if (L3 < 1.0e-300 || Ap < 1.0e-300)
            return;
        double iL3 = 1.0 / L3;

        // G1

        G(0, 2) = x32 * m;
        G(1, 2) = y32 * m;
        G(2, 1) = -iL3;

        // G2

        G(0, 8) = x13 * m;
        G(1, 8) = y13 * m;
        G(2, 7) = iL3;

        // G3

        G(0, 14) = x21 * m;
        G(1, 14) = y21 * m;

#endif // USE_POLAR_DECOMP_ALLIGN

    }

private:

    QuaternionType m_Q0;
    Vector3Type m_C0;
    std::array<QuaternionType, 3> m_QN;
    std::array<Vector3Type, 3> m_RV;
    std::array<QuaternionType, 3> m_QN_converged;
    std::array<Vector3Type, 3> m_RV_converged;
};

#endif // !ASDShellT3CorotationalTransformation_h
