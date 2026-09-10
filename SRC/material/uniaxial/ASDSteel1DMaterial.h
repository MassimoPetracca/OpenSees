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

// $Revision: 1.0 $
// $Date: 2025-01-03 11:29:01 $
// $Source: /usr/local/cvs/OpenSees/SRC/material/uniaxial/ASDSteel1DMaterial.cpp,v $

// Alessia Casalucci, Massimo Petracca, Guido Camata - ASDEA Software, Italy
//
// A unified and efficient plastic-damage material model for steel bars including fracture, bond-slip, and buckling via multiscale homogenization
//

#ifndef ASDSteel1DMaterial_h
#define ASDSteel1DMaterial_h

#include <UniaxialMaterial.h>
#include <Vector.h>
#include <Matrix.h>
#include <IMPLEXManager.h>
#include <cmath>
#include <memory>
#include <vector>
#include <map>

class ASDSteel1DMaterialPIMPL;

// This material runs an IMPL-EX scheme, so it takes part in the IMPL-EX error
// control: it registers itself, says when it took part in a step, and measures
// its own extrapolation error when asked. See IMPLEXManager.h
class ASDSteel1DMaterial : public UniaxialMaterial, public IMPLEXObject
{
public:
	class InputParameters {
	public:
		// Young's modulus
		double E = 0.0;
		// Yield stress
		double sy = 0.0;
		// ultimate strain for damage initialization
		double eu = 0.0;
		// Chaboche kinematic hardening parameters
		double H1 = 0.0;
		double H2 = 0.0;
		double gamma1 = 0.0;
		double gamma2 = 0.0;
		// THE BAUSCHINGER DIAL, in [0, 1]: 0 keeps the original calibration
		// (bit-identical response), 1 is the recipe tuned on Menegotto-Pinto
		// (OpenSees-Testing/asd-steel-3d, ex04/ex05). Each identity pair is a
		// Voce softening term (amplitude pair_Q >= 0, rate pair_b) that ships
		// with a kinematic twin at the SAME rate (H = Q*b, gamma = b), derived
		// - never stored - so an unpaired softening cannot exist. On the
		// monotonic backbone the twin cancels its softening identically, so
		// the backbone is the original one for EVERY dial value; on reversal
		// the twin must travel 2*Q to change sign while the softening only
		// follows the accumulated p - the released difference IS the
		// Bauschinger effect. The twin also protects the return mapping: it
		// steepens the local Newton residual two to one where the softening
		// flattens it, so -F'/E never drops below 1.
		double bauschinger = 0.0;
		static constexpr int NPAIRS = 3;
		static constexpr int NKIN = 2 + NPAIRS;
		double pair_Q[NPAIRS] = { 0.0, 0.0, 0.0 };
		double pair_b[NPAIRS] = { 0.0, 0.0, 0.0 };
		inline bool hasBauschinger() const { return bauschinger > 0.0; }
		// the isotropic part of the yield radius, R(p) <= 0, and its
		// derivative. Closed-form in the accumulated plastic multiplier: the
		// pairs own no state of their own beyond their twin backstresses
		inline double isoR(double p) const {
			double r = 0.0;
			for (int j = 0; j < NPAIRS; ++j)
				r -= pair_Q[j] * (1.0 - std::exp(-pair_b[j] * p));
			return r;
		}
		inline double isoRprime(double p) const {
			double r = 0.0;
			for (int j = 0; j < NPAIRS; ++j)
				r -= pair_Q[j] * pair_b[j] * std::exp(-pair_b[j] * p);
			return r;
		}
		// the kinematic terms, originals first, twins after. Fills H and g
		// (sized NKIN) and returns how many are active: 2 when the dial is at
		// zero, so the default pays nothing for the machinery
		inline int kinTerms(double* H, double* g) const {
			H[0] = H1; g[0] = gamma1;
			H[1] = H2; g[1] = gamma2;
			if (!hasBauschinger())
				return 2;
			for (int j = 0; j < NPAIRS; ++j) {
				H[2 + j] = pair_Q[j] * pair_b[j];
				g[2 + j] = pair_b[j];
			}
			return NKIN;
		}
		// misc
		bool implex = false;
		bool implex_control = false;
		// LEGACY: let the material itself fail the step when its own error is
		// out of tolerance. Off by default, and it should stay off: a material
		// that aborts owns a policy it cannot see the analysis to choose. The
		// rejection belongs to the convergence test wrapper - see IMPLEXManager.h
		bool implex_abort_on_error = false;
		double implex_error_tolerance = 0.0;
		double implex_time_redution_limit = 0.0;
		bool auto_regularization = true;
		bool buckling = false;
		bool fracture = false;
		bool slip = false;
		// buckling
		double radius = 0.0;
		double length = 0.0;
		double lch_element = 0.0;

		//convergence
		double K_alpha = 0.0;
		double max_iter = 0.0;
		double tolU = 0.0;
		double tolR = 0.0;
		// counter: MUST match the number of fields sendSelf/recvSelf write and
		// read here. It said 19 while 22 were being written, so the vector was
		// undersized and the tail of the parameters was landing past its end.
		// 23 original + the dial + NPAIRS*(Q, b)
		static constexpr int NDATA = 24 + 2 * NPAIRS;
	};

public:
	// life-cycle
	ASDSteel1DMaterial(
		int _tag,
		const InputParameters& _params,
		UniaxialMaterial* slip_material);
	ASDSteel1DMaterial();
	ASDSteel1DMaterial(const ASDSteel1DMaterial& other);
	~ASDSteel1DMaterial();

	// info
	const char* getClassType(void) const { return "ASDSteel1DMaterial"; }

	// set strain
	int setTrialStrain(double v, double r = 0.0);

	// get state
	double getStrain(void);
	double getStress(void);
	double getTangent(void);
	double getInitialTangent(void);

	// handle state
	int commitState(void);
	int revertToLastCommit(void);
	int revertToStart(void);

	// copy and others...
	UniaxialMaterial* getCopy(void);
	void Print(OPS_Stream& s, int flag = 0);

	// send/recv self
	int sendSelf(int commitTag, Channel& theChannel);
	int recvSelf(int commitTag, Channel& theChannel, FEM_ObjectBroker& theBroker);

	// parameters and responses
	int setParameter(const char** argv, int argc, Parameter& param);
	int updateParameter(int parameterID, Information& info);
	Response* setResponse(const char** argv, int argc, OPS_Stream& output);
	int getResponse(int responseID, Information& matInformation);
	double getEnergy(void);

	// IMPL-EX error control (see IMPLEXManager.h)
	double computeImplexErrorMetric(void);
	double implexTimeRatio(void) const;

private:
	// THE response, in one place: the explicit (extrapolated) answer when
	// do_implex is true, the implicit one when it is false. setTrialStrain, the
	// commit correction and the error peek all go through here, so the two
	// passes cannot drift apart - which is how the same step used to be written
	// four times over, once per branch of buckling x implex
	int computeResponse(bool do_implex);
	int homogenize(bool do_implex);
	// the largest stress this material can carry, the denominator of the metric
	double stressReference(void) const;
	// THE metric: a normalized stress gap
	double implexStressGap(double delivered, double stress_implicit) const;
	// a DIAGNOSTIC, not the metric: how far apart the two passes leave the RVE
	// deformed shape. See the note on it in the .cpp
	double implexDisplacementGap(const Vector& delivered, const Vector& implicit_u, bool elastic_layout) const;
	const Vector& getBucklingIndicator() const;
	const Vector& getDamage() const;
	const Vector& getEqPlStrain() const;
	const Vector& getAccumulatedPlasticStrain() const;
	const Vector& getSlipResponse() const;
	const Vector& getSteelResponse() const;
	const Vector& getTimeIncrements() const;
	const Vector& getImplexError() const;
	const Vector& getImplexErrorU() const;
	const Vector& getImplexStress() const;


 private:	
	 // common input parameters
	 InputParameters params;
	 // state variables - implex
	 double dtime_n = 0.0;
	 double dtime_n_commit = 0.0;
	 double dtime_0 = 0.0;
	 bool dtime_is_user_defined = false;
	 bool commit_done = false;
	 double implex_error = 0.0;
	 // a DIAGNOSTIC kept next to the metric and never mixed into it: see
	 // implexDisplacementGap()
	 double implex_error_u = 0.0;
	 // strain, stress and tangent (homogenized)
	 double strain = 0.0;
	 double strain_commit = 0.0;
	 double stress = 0.0;
	 double stress_commit = 0.0;
	 // what the step DELIVERED: the extrapolated stress, recorded before the
	 // commit correction overwrites 'stress' with the implicit one
	 double stress_implex = 0.0;
	 double C = 0.0;
	 double stress_rve = 0.0;
	 double stress_rve_commit = 0.0;
	 double C_rve = 0.0;
	 double N_rve_last = 0.0;
	 
	 // other variables for output purposes
	 double energy = 0.0;

	 // private implementation
	 ASDSteel1DMaterialPIMPL* pdata = nullptr;
};

#endif
