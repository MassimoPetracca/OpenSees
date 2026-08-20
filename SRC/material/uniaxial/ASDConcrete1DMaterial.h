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
// $Date: 2042-06-14 11:29:01 $
// $Source: /usr/local/cvs/OpenSees/SRC/material/uniaxial/ASDConcrete1DMaterial.cpp,v $

// Massimo Petracca - ASDEA Software, Italy
//
// A Simple and robust plastic-damage model for concrete and masonry
//

#ifndef ASDConcrete1DMaterial_h
#define ASDConcrete1DMaterial_h

#include <UniaxialMaterial.h>
#include <ASDHardeningLaw.h>
#include <IMPLEXManager.h>
#include <Vector.h>
#include <Matrix.h>
#include <cmath>
#include <memory>
#include <vector>
#include <map>

class ASDConcrete1DMaterial : public UniaxialMaterial, public IMPLEXObject
{
public:
	// sub-classes

	// The hardening law, its points and its enumerations are SHARED with the
	// other ASD materials, see ASDHardeningLaw.h. These aliases keep the
	// nested names this class has always exposed.
	using HardeningLawPoint = ASDHardeningLawPoint;
	using HardeningLawType = ASDHardeningLawType;
	using HardeningLawPointComponent = ASDHardeningLawPointComponent;
	using HardeningLaw = ASDHardeningLaw;

	// Everything the implicit pass of compute() can write. Saved and restored
	// around the non-destructive measurement of the IMPL-EX error, so that a
	// step that is measured and then rejected is a step that never happened.
	// It has to be COMPLETE: revertToLastCommit() does not restore PT_commit,
	// which is a committed quantity that the implicit pass overwrites.
	struct TrialState {
		double xt = 0.0;
		double xc = 0.0;
		double stress = 0.0;
		double stress_eff = 0.0;
		double C = 0.0;
		double dt_bar = 0.0;
		double dc_bar = 0.0;
		double PT_commit = 0.5;
	};

public:
	// life-cycle
	ASDConcrete1DMaterial(
		int _tag,
		double _E,
		double _eta,
		bool _implex,
		bool _implex_control,
		bool _implex_abort_on_error,
		double _implex_error_tolerance,
		double _implex_time_reduction_limit,
		double _implex_alpha,
		bool _tangent,
		bool _auto_regularize,
		double _lch_ref,
		const HardeningLaw& _ht,
		const HardeningLaw& _hc);
	ASDConcrete1DMaterial();
	~ASDConcrete1DMaterial();

	// info
	const char* getClassType(void) const { return "ASDConcrete1DMaterial"; }

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

	// IMPL-EX error control
	double computeImplexErrorMetric(void);
	double implexTimeRatio(void) const;

private:
	// internal computation
	int compute(bool do_implex, bool do_tangent);
	// the IMPL-EX error metric: the gap between the stress the extrapolated
	// step DELIVERS and the one the implicit solution would carry at the same
	// strain, over the largest stress this material can reach.
	//
	// It is NOT a difference of damage, which is what this material used to
	// measure. A damage metric is blind - exactly, not approximately - to any
	// error that moves the stress at FIXED damage, and the frozen sign switch
	// of the split is exactly such an error: at the worst step of a measured
	// run the effective stress, the plastic strain, both strain measures and
	// both damage values are identical between the two passes, so the damage
	// metric reads 0.000e+00, while the delivered stress is wrong by 11.6 MPa
	// out of a 30 MPa strength. The legacy quantity survives as implex_gap,
	// not used by the model, to price the difference.
	// 'delivered', not 'stress_implex': that is the name of the member now, and a
	// parameter shadowing it would read as the member to anyone skimming
	double implexStressGap(double delivered, double stress_implicit) const;
	// the legacy metric, kept for output only
	double implexDamageGap(double dt_implex, double dc_implex,
		double dt_implicit, double dc_implicit) const;
	// save/restore everything the implicit pass writes
	void saveTrialState(TrialState& x) const;
	void restoreTrialState(const TrialState& x);
	Vector getHardeningLawVector(HardeningLawType ltype, HardeningLawPointComponent c) const;
	const Vector& getStrainMeasure() const;
	const Vector& getDamage() const;
	const Vector& getEquivalentPlasticStrain() const;
	const Vector& getCrackWidth() const;
	const Vector& getCrushWidth() const;
	const Vector& getImplexError() const;
	const Vector& getImplexStress() const;
	const Vector& getImplexGap() const;
	const Vector& getTimeIncrements() const;

 private:
	 // Young's modulus
	 double E = 0.0;
	 // Viscosity for rate-dependent damage
	 double eta = 0.0;
	 // True = use the IMPL-EX algorithm
	 bool implex = false;
	 // True = measure the IMPL-EX error at every setTrialStrain and publish it.
	 // This is the LEGACY path: it pays one implicit solve per Newton
	 // iteration, at strains that are not equilibrated yet. The measurement
	 // the error control actually runs on is the one the convergence test
	 // wrapper asks for, once per step, through computeImplexErrorMetric()
	 bool implex_control = false;
	 // True = let this material fail the step on its own when the IMPL-EX
	 // error is over tolerance. OFF by default, and it should stay off: a
	 // material that aborts owns a policy it cannot see the analysis to
	 // choose, and the code it returns reaches the element as a material
	 // failure rather than as a controlled rejection. Kept for the inputs
	 // that were written against the old behaviour
	 bool implex_abort_on_error = false;
	 // Maximum allowed IMPL-EX error (default = 5%)
	 double implex_error_tolerance = 0.05;
	 // Minimum allowed time step reduction factor under which IMPL-EX error is not controlled anymore (default = 1%)
	 double implex_time_redution_limit = 0.01;
	 // Scale factor for implex extrapolation
	 double implex_alpha = 1.0;
	 // True = use the tangent matrix, False (default) = use the secant matrix
	 bool tangent = false;
	 // True = automatically regularize the fracture energy using the element's characteristic length
	 bool auto_regularize = true;
	 bool regularization_done = false;
	 double lch = 1.0; // the parent-element's characteristic length
	 double lch_ref = 1.0; // the reference characteristic length (i.e. the size the specific-fracture-energy in the hardening-law is referred to)
	 // The hardening law for the tensile response
	 HardeningLaw ht;
	 // The hardening law for the compressive response
	 HardeningLaw hc;
	 // state variables - tension
	 double xt = 0.0;
	 double xt_commit = 0.0;
	 double xt_commit_old = 0.0;
	 // state variables - compression
	 double xc = 0.0;
	 double xc_commit = 0.0;
	 double xc_commit_old = 0.0;
	 // state variables - implex
	 double dtime_n = 0.0;
	 double dtime_n_commit = 0.0;
	 double dtime_0 = 0.0;
	 bool dtime_is_user_defined = false;
	 bool commit_done = false;
	 double implex_error = 0.0;
	 // the legacy damage-based metric. Computed for output, never used by the
	 // model and never compared against a tolerance
	 double implex_gap = 0.0;
	 double PT_commit = 0.5;
	 // strain, stress and tangent
	 double strain = 0.0;
	 double strain_commit = 0.0;
	 double stress = 0.0;
	 double stress_commit = 0.0;
	 // WHAT THE STEP DELIVERED TO THE ELEMENT, kept because it is otherwise
	 // unobservable: under IMPL-EX commitState re-solves implicitly and installs
	 // that answer over `stress`, so from outside - a recorder, a test - the
	 // extrapolated stress is gone by the time anyone can ask. It is what the
	 // error metric measures the distance from, and it used to live as a local
	 // in the two places that measure. On the implicit path there is no
	 // extrapolation and it equals `stress`
	 double stress_implex = 0.0;
	 double stress_eff = 0.0;
	 double stress_eff_commit = 0.0;
	 double C = 0.0;
	 // other variables for output purposes
	 double dt_bar = 0.0;
	 double dc_bar = 0.0;
	 double energy = 0.0;
};

#endif
