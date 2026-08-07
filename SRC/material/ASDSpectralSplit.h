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

// Massimo Petracca - ASDEA Software, Italy
//
// The 3x3 symmetric eigen-decomposition and the fourth-order spectral split
// projectors, shared by the ASD 3D materials.
//
// WHY IT IS HERE. ASDConcrete3DMaterial.cpp carries its own copy of the Jacobi
// solver in an anonymous namespace, with the comment "Taken from Eigen3 in
// Matrix class. We need it because we have to access eigenvectors" - i.e. it is
// ALREADY a copy, made only to reach the eigenvectors that Matrix::Eigen3 does
// not return. A second 3D material needing the same thing is one copy too many:
// divergent copies of this exact kind of code are how the error metric of this
// family became non-uniform in the first place. This is the single version.
//
// ASDConcrete3DMaterial IS WIRED TO IT since 2026-08-07: its 193 private lines
// (Eigen3 = 142, computePjj = 51) are gone, replaced by an include and a
// forwarder that passes normalize=false. That keeps it BIT-IDENTICAL, measured -
// see below and the note at the forwarder. It is the one place on this
// long-running private branch where a file that also exists upstream loses lines
// rather than only gaining them, and it was a deliberate, separately verified
// step, not a side effect of adding a material.
//
// ONE CHANGE FROM THAT COPY: the input can be NORMALIZED before the sweeps, and
// that is why it is a FLAG rather than the only behaviour. The routine's
// convergence test is `sm > 1.0e-8` with sm = |a01| + |a12| + |a20|, an ABSOLUTE
// threshold in the units of whatever it is handed. Below that scale the
// while-loop never runs and the routine returns the BARE DIAGONAL, with V = I -
// no rotation at all, silently. On a strain-scale tensor (1e-3 and below) that
// always bites; at stress scale it does not, which is why it survived.
//
// WHAT NORMALIZING COSTS, MEASURED, because the first version of this comment
// claimed it "changes nothing else" and that is not true. Turning the threshold
// relative also LOOSENS it wherever the tensor is larger than 1e-8: at scale 30
// the absolute test demanded ~3e-10 relative, the normalized one demands 1e-8.
// On 4000 random tensors per scale, normalized against the un-normalized copy:
// the eigenvalues still agree to ~1e-15 relative, but the eigenVECTORS move by
// up to 5.8e-8 at scale 30 and 5.2e-8 at 1e+6. So it is a fix at small scale (at
// 1e-8 and 1e-12 the un-normalized routine puts 339 and 1832 eigenvalues on the
// wrong SIDE of zero out of 12000 - it is not computing anything there) and a
// mild loosening at large scale. A caller that wants both should tighten `tol`,
// which is now a relative quantity; that is a separate change with its own
// verification and nobody has needed it yet.
//
// CONVENTIONS, stated once:
//
// * eigenvalues come out DESCENDING, so d(0) is the maximum. That is what the
//   copy in ASDConcrete3DMaterial does, and every caller depends on it;
// * the Voigt order is (11, 22, 33, 12, 23, 13), the OpenSees one;
// * the 6x6 projectors are STRESS-IN / STRESS-OUT: they act on a 6-vector
//   holding the TENSOR components of a symmetric tensor (sigma_12, not a
//   doubled one), and the doubling that the full contraction over the
//   off-diagonal pairs implies is carried by their shear COLUMNS. So the same
//   matrix applies correctly to a stress and to a tensor-component strain, and
//   the product of two of them composes as the fourth-order composition does.

#ifndef ASDSpectralSplit_h
#define ASDSpectralSplit_h

#include <Vector.h>
#include <Matrix.h>
#include <cmath>
#include <algorithm>

namespace ASDSpectralSplit {

	/**
	The eigen-decomposition of a 3x3 symmetric matrix, by Jacobi rotations.

	On input v holds the matrix; on output d holds the eigenvalues DESCENDING
	and the columns of v the corresponding eigenvectors. Returns 0, or -1 on a
	size mismatch.

	The body is the one in ASDConcrete3DMaterial.cpp, itself taken from
	Matrix::Eigen3. Two differences, both deliberate: the normalization
	described in the file header, and the three function-level `static Vector`
	scratch arrays replaced by plain locals - three doubles each, so the static
	bought nothing and cost re-entrancy.

	WHY normalize IS A FLAG AND NOT A CONSTANT. Normalizing is a fix, but it is
	not arithmetically free: dividing by the scale and multiplying back is not
	bit-exact, so switching an existing material onto this routine WITH it would
	move that material's answers in the last digits. ASDConcrete3DMaterial passes
	false, which makes its use of this shared code bit-identical to the private
	copy it replaced - see the measurement in project_opensees_cdp3d_port. Turning
	it on there is a separate decision with a separate verification, not a side
	effect of a de-duplication.
	*/
	inline int eigen3(Vector& d, Matrix& v, bool normalize = true)
	{
		if (v.noRows() != 3 || v.noCols() != 3 || d.Size() != 3)
			return -1;

		// NORMALIZE. The sweep threshold below is absolute, so without this the
		// routine is a no-op on anything smaller than 1e-8 - see the file
		// header. max(|entry|) rather than a norm because it is what makes the
		// largest entry exactly 1 and cannot itself overflow.
		double scale = 0.0;
		bool scaled = false;
		if (normalize) {
			for (int i = 0; i < 3; ++i)
				for (int j = 0; j < 3; ++j)
					scale = std::max(scale, std::abs(v(i, j)));
			scaled = (scale > 0.0) && std::isfinite(scale) && (scale != 1.0);
			if (scaled) {
				for (int i = 0; i < 3; ++i)
					for (int j = 0; j < 3; ++j)
						v(i, j) /= scale;
			}
			if (!(scale > 0.0) || !std::isfinite(scale)) {
				// the zero tensor (or a broken one): the answer is exact and the
				// sweeps have nothing to do. Returning the identity basis here
				// also keeps a caller from reading an uninitialized V
				d.Zero();
				v.Zero();
				for (int i = 0; i < 3; ++i)
					v(i, i) = 1.0;
				if (!std::isfinite(scale))
					return -1;
				return 0;
			}
		}

		int     rot, its, i, j, k;
		double  g, h, aij, sm, thresh, t, c, s, tau;

		double a[3];
		double b[3];
		double z[3];

		static const double tol = 1.0e-08;

		//.... move array into one-d arrays
		a[0] = v(0, 1);
		a[1] = v(1, 2);
		a[2] = v(2, 0);

		for (i = 0; i < 3; i++) {
			d(i) = v(i, i);
			b[i] = v(i, i);
			z[i] = 0.0;

			for (j = 0; j < 3; j++)
				v(i, j) = 0.0;

			v(i, i) = 1.0;

		} //end for i

		rot = 0;
		its = 0;

		sm = fabs(a[0]) + fabs(a[1]) + fabs(a[2]);

		while (sm > tol) {
			//.... set convergence test and threshold
			if (its < 3)
				thresh = 0.011 * sm;
			else
				thresh = 0.0;

			//.... perform sweeps for rotations
			for (i = 0; i < 3; i++) {

				j = (i + 1) % 3;
				k = (j + 1) % 3;

				aij = a[i];

				g = 100.0 * fabs(aij);

				if (fabs(d(i)) + g != fabs(d(i)) ||
					fabs(d(j)) + g != fabs(d(j))) {

					if (fabs(aij) > thresh) {

						a[i] = 0.0;
						h = d(j) - d(i);

						if (fabs(h) + g == fabs(h))
							t = aij / h;
						else {
							double hDIVaij = h / aij;
							if (hDIVaij > 0.0)
								t = 2.0 / (hDIVaij + sqrt(4.0 + (hDIVaij * hDIVaij)));
							else
								t = -2.0 / (-hDIVaij + sqrt(4.0 + (hDIVaij * hDIVaij)));
						}

						//.... set rotation parameters

						c = 1.0 / sqrt(1.0 + t * t);
						s = t * c;
						tau = s / (1.0 + c);

						//.... rotate diagonal terms

						h = t * aij;
						z[i] = z[i] - h;
						z[j] = z[j] + h;
						d(i) = d(i) - h;
						d(j) = d(j) + h;

						//.... rotate off-diagonal terms

						h = a[j];
						g = a[k];
						a[j] = h + s * (g - h * tau);
						a[k] = g - s * (h + g * tau);

						//.... rotate eigenvectors

						for (k = 0; k < 3; k++) {
							g = v(k, i);
							h = v(k, j);
							v(k, i) = g - s * (h + g * tau);
							v(k, j) = h + s * (g - h * tau);
						} // end for k

						rot = rot + 1;

					} // end if fabs > thresh
				} //else
				else
					a[i] = 0.0;

			}  // end for i

			//.... update the diagonal terms
			for (i = 0; i < 3; i++) {
				b[i] = b[i] + z[i];
				d(i) = b[i];
				z[i] = 0.0;
			} // end for i

			its += 1;

			sm = fabs(a[0]) + fabs(a[1]) + fabs(a[2]);

		} //end while sm

		// sort in descending order (unrolled bubble sort)
		auto sortij = [&d, &v](int i, int j) {
			if (d(i) < d(j)) {
				std::swap(d(i), d(j));
				for (int k = 0; k < 3; ++k)
					std::swap(v(k, i), v(k, j));
			}
		};
		sortij(0, 1);
		sortij(1, 2);
		sortij(0, 1);

		// undo the normalization. The eigenVECTORS are untouched by it
		if (scaled) {
			for (i = 0; i < 3; ++i)
				d(i) *= scale;
		}

		// done
		return 0;
	}

	/**
	The eigenvalues alone, DESCENDING, of a symmetric tensor given in Voigt with
	its TENSOR shear components. Returns 0 or -1.
	*/
	inline int eigenvalues(const Vector& t, Vector& d)
	{
		static Matrix V(3, 3);
		V(0, 0) = t(0);
		V(1, 1) = t(1);
		V(2, 2) = t(2);
		V(0, 1) = V(1, 0) = t(3);
		V(1, 2) = V(2, 1) = t(4);
		V(0, 2) = V(2, 0) = t(5);
		return eigen3(d, V);
	}

	/**
	The eigenvalues and eigenvectors of a symmetric tensor given in Voigt with
	its TENSOR shear components. d DESCENDING, V holding the directions in
	columns. Returns 0 or -1.
	*/
	inline int spectral(const Vector& t, Vector& d, Matrix& V)
	{
		V(0, 0) = t(0);
		V(1, 1) = t(1);
		V(2, 2) = t(2);
		V(0, 1) = V(1, 0) = t(3);
		V(1, 2) = V(2, 1) = t(4);
		V(0, 2) = V(2, 0) = t(5);
		return eigen3(d, V);
	}

	/**
	The 6x6 Voigt form of vj (x) vj (x) vj (x) vj, with vj the j-th column of V.

	Verbatim from ASDConcrete3DMaterial.cpp, expression by expression. Its shear
	COLUMNS carry the factor two of the contraction over the off-diagonal pairs,
	which is what makes it applicable to a 6-vector holding tensor components -
	see the file header.
	*/
	inline void computePjj(const Matrix& V, int j, Matrix& pjj)
	{
		double A1 = V(2, j) * V(2, j);
		double A2 = V(1, j) * V(1, j);
		double A3 = V(0, j) * V(0, j);
		double A4 = V(2, j);
		double A5 = V(0, j);
		double A6 = V(1, j);
		double A7 = 2.0 * A1 * A5 * A6;
		double A8 = 2.0 * A3 * A4 * A6;
		double A9 = 2.0 * A2 * A4 * A5;
		double A10 = A1 * A2;
		pjj(0, 0) = A3 * A3;
		pjj(0, 1) = A2 * A3;
		pjj(0, 2) = A1 * A3;
		pjj(0, 3) = 2.0 * A3 * A5 * A6;
		pjj(0, 4) = A8;
		pjj(0, 5) = 2.0 * A3 * A4 * A5;
		pjj(1, 0) = A2 * A3;
		pjj(1, 1) = A2 * A2;
		pjj(1, 2) = A10;
		pjj(1, 3) = 2.0 * A2 * A5 * A6;
		pjj(1, 4) = 2.0 * A2 * A4 * A6;
		pjj(1, 5) = A9;
		pjj(2, 0) = A1 * A3;
		pjj(2, 1) = A10;
		pjj(2, 2) = A1 * A1;
		pjj(2, 3) = A7;
		pjj(2, 4) = 2.0 * A1 * A4 * A6;
		pjj(2, 5) = 2.0 * A1 * A4 * A5;
		pjj(3, 0) = A3 * A5 * A6;
		pjj(3, 1) = A2 * A5 * A6;
		pjj(3, 2) = A1 * A5 * A6;
		pjj(3, 3) = 2.0 * A2 * A3;
		pjj(3, 4) = A9;
		pjj(3, 5) = A8;
		pjj(4, 0) = A3 * A4 * A6;
		pjj(4, 1) = A2 * A4 * A6;
		pjj(4, 2) = A1 * A4 * A6;
		pjj(4, 3) = A9;
		pjj(4, 4) = 2.0 * A1 * A2;
		pjj(4, 5) = A7;
		pjj(5, 0) = A3 * A4 * A5;
		pjj(5, 1) = A2 * A4 * A5;
		pjj(5, 2) = A1 * A4 * A5;
		pjj(5, 3) = A8;
		pjj(5, 4) = A7;
		pjj(5, 5) = 2.0 * A1 * A3;
	}

	/**
	Relative half-width of the band around zero inside which an eigenvalue is
	treated as neither positive nor negative by splitFromSpectral. See there.
	*/
	static const double SplitTolerance = 1.0e-12;

	/**
	(PT, PC): the fourth-order split of a symmetric tensor by the SIGN of its
	eigenvalues, from an already computed decomposition. PT:s is the positive
	part, PC:s the negative one, and PT + PC == I.

	THE REMAINDER IS SHARED, NOT DROPPED, and that is not a detail. PO = I - PT
	- PC holds the zero eigenvalues AND the whole off-diagonal subspace;
	leaving it out makes PT + PC != I, which zeroes the shear response of the
	operator - measured on the model this was written for, the tangent came out
	with three null eigenvalues. On the VALUE the remainder changes nothing
	(both subspaces carry no stress in the principal basis), so half and half is
	free; on the TANGENT it is the symmetric element of the subdifferential.

	THE SIGN TEST IS ON A BAND AND NOT ON ZERO. A strict w > 0 puts an
	eigenvalue of +1e-18 into PT and one of -1e-18 into PC. Those are the same
	state - the difference is round-off in an eigen-solver - but the
	eigenvectors spanning a numerically degenerate subspace are ARBITRARY, so
	which direction lands where is arbitrary too, and PT then differs by a whole
	n (x) n between two evaluations of one state. The VALUE never notices, which
	is why this kind of bug is invisible: a near-zero eigenvalue carries a
	near-zero stress whichever reduction multiplies it. The DERIVATIVE notices,
	and a Newton corrector that lands on the resulting kink chatters. Putting
	the whole degenerate group into PO fixes it because PO is formed as a
	REMAINDER: the arbitrary directions never appear in it individually, only
	the group's own projector does, and that is basis-independent.

	The band is relative to the largest eigenvalue, so it is a statement about
	significant digits and not about stress units, and it is not a delicate
	number: anything from 1e-14 to 1e-10 removes the chatter and stays far below
	the last digit of any eigenvalue a real state has.

	tol_rel = 0 recovers the strict sign test, up to rounding: it is what
	ASDConcrete3DMaterial's StressDecomposition does through Heavyside(), whose
	0.5 at an exact zero puts half of that projector in each side - the same
	number this reaches through PO, by a different sum.
	*/
	inline void splitFromSpectral(const Vector& d, const Matrix& V,
		Matrix& PT, Matrix& PC, double tol_rel = SplitTolerance)
	{
		double amax = 0.0;
		for (int j = 0; j < 3; ++j)
			amax = std::max(amax, std::abs(d(j)));
		double tol = tol_rel * amax;

		static Matrix pjj(6, 6);
		PT.Zero();
		PC.Zero();
		for (int j = 0; j < 3; j++) {
			if (d(j) > tol) {
				computePjj(V, j, pjj);
				PT.addMatrix(1.0, pjj, 1.0);
			}
			else if (d(j) < -tol) {
				computePjj(V, j, pjj);
				PC.addMatrix(1.0, pjj, 1.0);
			}
		}
		// PO = I - PT - PC, then half to each
		static Matrix PO(6, 6);
		PO.addMatrix(0.0, PT, -1.0);
		PO.addMatrix(1.0, PC, -1.0);
		for (int i = 0; i < 6; ++i)
			PO(i, i) += 1.0;
		PT.addMatrix(1.0, PO, 0.5);
		PC.addMatrix(1.0, PO, 0.5);
	}

	/**
	splitFromSpectral straight off a tensor given in Voigt with its TENSOR shear
	components. Returns 0, or -1 if the decomposition failed.
	*/
	inline int split(const Vector& t, Matrix& PT, Matrix& PC,
		double tol_rel = SplitTolerance)
	{
		static Vector d(3);
		static Matrix V(3, 3);
		if (spectral(t, d, V) < 0)
			return -1;
		splitFromSpectral(d, V, PT, PC, tol_rel);
		return 0;
	}

}

#endif // ASDSpectralSplit_h
