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

#include <ASDConcrete1DMaterial.h>
#include <Channel.h>
#include <OPS_Globals.h>
#include <Information.h>
#include <Parameter.h>
#include <elementAPI.h>
#include <Element.h>
#include <MaterialResponse.h>
#include <cmath>
#include <algorithm>
#include <limits>
#include <string>
#include <sstream>
#include <iomanip>

// anonymous namespace for utilities
namespace {

	enum ErrorCodes {
		EC_Generic = -1,
		EC_IMPLEX_Error_Control = -10
	};

	/**
	Converts a string into a vector of doubles using whitespace as delimiter
	*/
	bool string_to_double(const std::string& text, double& num) {
		num = 0.0;
		try {
			num = std::stod(text);
			return true;
		}
		catch (...) {
			return false;
		}
	}
	bool string_to_list_of_doubles(const std::string& text, char sep, std::vector<double>& out) {
		if (out.size() > 0) out.clear();
		std::size_t start = 0, end = 0;
		double value;
		while (true) {
			end = text.find(sep, start);
			if (end == std::string::npos) {
				if (start < text.size()) {
					if (!string_to_double(text.substr(start), value))
						return false;
					out.push_back(value);
				}
				break;
			}
			std::string subs = text.substr(start, end - start);
			if (subs.size() > 0) {
				if (!string_to_double(subs, value))
					return false;
				out.push_back(value);
			}
			start = end + 1;
		}
		return true;
	}

	// Heavyside function
	inline double Heavyside(double X) { return X > 0.0 ? 1.0 : (X < 0.0 ? 0.0 : 0.5); }

	// Macauley function
	inline double Macauley(double X) { return X > 0.0 ? X : 0.0; }

	/**
	global parameters storage
	*/
	class GlobalParameters {
	private:
		double max_error = 0.0;
		double avg_error = 0.0;
		int avg_counter = 0;
	private:
		GlobalParameters() = default;
		GlobalParameters(const GlobalParameters&) = delete;
		GlobalParameters& operator = (const GlobalParameters&) = delete;
	public:
		static GlobalParameters& instance() {
			static GlobalParameters _instance;
			return _instance;
		}
		inline double getMaxError() const {
			return max_error;
		}
		inline void setMaxError(double x) {
			max_error = x;
		}
		inline double getAverageError() {
			if (avg_counter > 0) {
				avg_error /= static_cast<double>(avg_counter);
				avg_counter = 0;
			}
			return avg_error;
		}
		inline void accumulateAverageError(double x) {
			avg_error += x;
			++avg_counter;
		}
		inline void setAverageError(double x) {
			avg_error = x;
			avg_counter = 0;
		}
	};

	double bezier3(double xi,
		double x0, double x1, double x2,
		double y0, double y1, double y2)
	{
		double A = x0 - 2.0 * x1 + x2;
		double B = 2.0 * (x1 - x0);
		double C = x0 - xi;
		if (fabs(A) < 1.0e-12) {
			x1 = x1 + 1.0E-6 * (x2 - x0);
			A = x0 - 2.0 * x1 + x2;
			B = 2.0 * (x1 - x0);
			C = x0 - xi;
		}
		if (A == 0.0)
			return 0.0;

		double D = B * B - 4.0 * A * C;
		double t = (sqrt(D) - B) / (2.0 * A);

		return (y0 - 2.0 * y1 + y2) * t * t + 2.0 * (y1 - y0) * t + y0;
	}

}

void* OPS_ASDConcrete1DMaterial()
{
	// some kudos
	static bool first_done = false;
	if (!first_done) {
		opserr << "Using ASDConcrete1D - Developed by: Massimo Petracca, Guido Camata, ASDEA Software Technology\n";
		first_done = true;
	}

	// check arguments
	int numArgs = OPS_GetNumRemainingInputArgs();
	if (numArgs < 2) {
		opserr <<
			"nDMaterial ASDConcrete1D Error: Few arguments (< 2).\n"
			"nDMaterial ASDConcrete1D $tag $E "
			"<-fc $fc> <-ft $ft> "
			"<-Te $Te -Ts $Ts <-Td $Td>> <-Ce $Ce -Cs $Cs <-Cd $Cd>> "
			"<-implex> <-implexControl $implexErrorTolerance $implexTimeReductionLimit> <-implexAbort> <-implexAlpha $alpha>"
			"<-eta $eta> <-tangent> <-autoRegularization $lch_ref>\n";
		return nullptr;
	}

	// numData
	int numData = 1;

	// data
	int tag;
	double E;
	bool implex = false;
	bool implex_control = false;
	bool implex_abort_on_error = false;
	double implex_error_tolerance = 0.05;
	double implex_time_redution_limit = 0.01;
	double implex_alpha = 1.0;
	double eta = 0.0;
	bool tangent = false;
	bool auto_regularization = false;
	double lch_ref = 1.0;
	std::vector<double> Te, Ts, Td, Ce, Cs, Cd;

	// get tag
	if (OPS_GetInt(&numData, &tag) != 0) {
		opserr << "nDMaterial ASDConcrete1D Error: invalid 'tag'.\n";
		return nullptr;
	}

	// get Elasticity arguments
	if (OPS_GetDouble(&numData, &E) != 0) {
		opserr << "nDMaterial ASDConcrete1D Error: invalid 'E'.\n";
		return nullptr;
	}
	if (E <= 0.0) {
		opserr << "nDMaterial ASDConcrete1D Error: invalid value for 'E' (" << E << "). It should be strictly positive.\n";
		return nullptr;
	}

	// utilities (code re-use)
	auto lam_optional_int = [&numData](const char* variable, int& value) -> bool {
		if (OPS_GetNumRemainingInputArgs() > 0) {
			if (OPS_GetInt(&numData, &value) < 0) {
				opserr << "nDMaterial ASDConcrete1D Error: failed to get '" << variable << "'.\n";
				return false;
			}
		}
		else {
			opserr << "nDMaterial ASDConcrete1D Error: '" << variable << "' requested but not provided.\n";
			return false;
		}
		return true;
	};
	auto lam_optional_double = [&numData](const char* variable, double& value) -> bool {
		if (OPS_GetNumRemainingInputArgs() > 0) {
			if (OPS_GetDouble(&numData, &value) < 0) {
				opserr << "nDMaterial ASDConcrete1D Error: failed to get '" << variable << "'.\n";
				return false;
			}
		}
		else {
			opserr << "nDMaterial ASDConcrete1D Error: '" << variable << "' requested but not provided.\n";
			return false;
		}
		return true;
	};
	auto lam_optional_list = [&numData](const char* variable, std::vector<double>& value) -> bool {
		// first try expanded list like {*}$the_list,
		// also used in python like *the_list
		value.clear();
		while (OPS_GetNumRemainingInputArgs() > 0) {
			double item;
			auto old_num_rem = OPS_GetNumRemainingInputArgs();
			if (OPS_GetDoubleInput(&numData, &item) < 0) {
				auto new_num_rem = OPS_GetNumRemainingInputArgs();
				if (new_num_rem < old_num_rem)
					OPS_ResetCurrentInputArg(-1);
				break;
			}
			value.push_back(item);
		}
		// try Tcl list (it's a string after all...)
		if (value.size() == 0 && OPS_GetNumRemainingInputArgs() > 0) {
			std::string list_string = OPS_GetString();
			if (!string_to_list_of_doubles(list_string, ' ', value)) {
				opserr << "nDMaterial ASDConcrete1D Error: cannot parse the '" << variable << "' list.\n";
				return false;
			}
		}
		return true;
	};

	double fc;
	double ft;
	bool have_fc = false;
	bool have_ft = false;
	bool have_lch_ref = false;

	// optional parameters
	while (OPS_GetNumRemainingInputArgs() > 0) {
		const char* value = OPS_GetString();
		if (strcmp(value, "-fc") == 0) {
			if (!lam_optional_double("fc", fc))
				return nullptr;
			have_fc = true;
		}
		else if (strcmp(value, "-ft") == 0) {
			if (!lam_optional_double("ft", ft))
				return nullptr;
			have_ft = true;
		}
		else if (strcmp(value, "-implex") == 0) {
			implex = true;
		}
		else if (strcmp(value, "-implexControl") == 0) {
			implex_control = true;
			if (OPS_GetNumRemainingInputArgs() < 2) {
				opserr << "nDMaterial ASDConcrete1D Error: '-implexControl' given without the next 2 arguments $implexErrorTolerance $implexTimeReductionLimit.\n";
				return nullptr;
			}
			if (!lam_optional_double("implexErrorTolerance", implex_error_tolerance))
				return nullptr;
			if (!lam_optional_double("implexTimeReductionLimit", implex_time_redution_limit))
				return nullptr;
		}
		else if (strcmp(value, "-implexAbort") == 0) {
			// legacy: let the material fail the step by itself. See the note
			// on implex_abort_on_error in the header
			implex_abort_on_error = true;
		}
		else if (strcmp(value, "-implexAlpha") == 0) {
			if (!lam_optional_double("alpha", implex_alpha))
				return nullptr;
		}
		else if (strcmp(value, "-eta") == 0) {
			if (!lam_optional_double("eta", eta))
				return nullptr;
		}
		else if (strcmp(value, "-tangent") == 0) {
			tangent = true;
		}
		else if (strcmp(value, "-autoRegularization") == 0) {
			auto_regularization = true;
			if (OPS_GetNumRemainingInputArgs() < 1) {
				opserr << "nDMaterial ASDConcrete1D Error: '-autoRegularization' given without the next 1 argument $lch_ref.\n";
				return nullptr;
			}
			if (!lam_optional_double("lch_ref", lch_ref))
				return nullptr;
			have_lch_ref = true;
		}
		else if (strcmp(value, "-Te") == 0) {
			if (!lam_optional_list("Te", Te))
				return nullptr;
		}
		else if (strcmp(value, "-Ts") == 0) {
			if (!lam_optional_list("Ts", Ts))
				return nullptr;
		}
		else if (strcmp(value, "-Td") == 0) {
			if (!lam_optional_list("Td", Td))
				return nullptr;
		}
		else if (strcmp(value, "-Ce") == 0) {
			if (!lam_optional_list("Ce", Ce))
				return nullptr;
		}
		else if (strcmp(value, "-Cs") == 0) {
			if (!lam_optional_list("Cs", Cs))
				return nullptr;
		}
		else if (strcmp(value, "-Cd") == 0) {
			if (!lam_optional_list("Cd", Cd))
				return nullptr;
		}
	}

	// Set a default value of tension strength if none specified
	if (have_fc && !have_ft)
		ft = 0.1 * fc;

	if (have_fc) {
		double ec = 2 * fc / E;
		double Gt = 0.073 * pow(fc, 0.18);
		double Gc = 2 * Gt * (fc * fc) / (ft * ft);


		if (!have_lch_ref) {
			//
			// _get_lch_ref from ASDConcrete1D_MakeLaws.py
			//

			// min lch for tension
			double et_el = ft / E;
			double Gt_min = 0.5 * ft * et_el;
			double hmin_t = 0.01 * Gt / Gt_min;

			// min lch for compression
			double ec1 = fc / E;
			double ec_pl = (ec - ec1) * 0.4 + ec1;
			double Gc_min = 0.5 * fc * (ec - ec_pl);
			double hmin_c = 0.01 * Gc / Gc_min;

			lch_ref = std::min(hmin_c, hmin_t);
		}

		//
		// _make_tension from ASDConcrete1D_MakeLaws.py
		//

		Gt = Gt / lch_ref;

		double f0 = 0.9 * ft;
		double f1 = ft;
		double e0 = f0 / E;
		double e1 = 1.5 * f1 / E;
		double ep = e1 - f1 / E;
		double f2 = 0.2 * ft;
		double f3 = 1.0e-3 * ft;
		double w2 = Gt / ft;
		double w3 = 5.0 * w2;
		double e2 = w2 + f2 / E + ep;
		if (e2 <= e1)
			e2 = 1.001 * e1;
		double e3 = w3 + f3 / E + ep;
		if (e3 <= e2)
			e3 = 1.001 * e2;
		double e4 = 10.0 * e3;
		Te.resize(6); Te = { 0.0, e0, e1, e2, e3, e4 };
		Ts.resize(6); Ts = { 0.0, f0, f1, f2, f3, f3 };
		Td.resize(6); Td = { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 };
		double Tpl[6] = { 0.0, 0.0, ep, 0.9 * e2, 0.8 * e3, 0.8 * e3 };

		for (int i = 2; i < 6; i++) {
			double xi = Te[i];
			double si = Ts[i];
			double xipl = Tpl[i];
			double xipl_max = xi - si / E;
			xipl = std::min(xipl, xipl_max);
			double qi = (xi - xipl) * E;
			Td[i] = 1.0 - si / qi;
		}

		//
		// _make_compression from ASDConcrete1D_MakeLaws.py
		//

		Gc = Gc / lch_ref;

		double fc0 = 0.5 * fc;
		double ec0 = fc0 / E;
		double ec1 = fc / E;
		double fcr = 0.1 * fc;
		double ec_pl = (ec - ec1) * 0.4 + ec1;
		double Gc1 = 0.5 * fc * (ec - ec_pl);
		double Gc2 = std::max(0.01 * Gc1, Gc - Gc1);
		double ecr = ec + 2.0 * Gc2 / (fc + fcr);
		const int nc = 10;
		Ce.resize(nc + 3); Ce[0] = 0.0; Ce[1] = ec0;
		Cs.resize(nc + 3); Cs[0] = 0.0; Cs[1] = fc0;
		double Cpl[nc + 3]; Cpl[0] = 0.0; Cpl[1] = 0.0;
		double dec = (ec - ec0) / (nc - 1);
		for (int i = 0; i < nc - 1; i++) {
			double iec = ec0 + (i + 1) * dec;
			Ce[i + 2] = iec;
			Cs[i + 2] = bezier3(iec, ec0, ec1, ec, fc0, fc, fc);
			Cpl[i + 2] = Cpl[i + 1] + 0.7 * (iec - Cpl[i + 1]);
		}
		Ce[nc + 1] = ecr;
		Cs[nc + 1] = fcr;
		Cpl[nc + 1] = Cpl[nc] + 0.7 * (ecr - Cpl[nc]);
		Ce[nc + 2] = ecr + ec0;
		Cs[nc + 2] = fcr;
		Cpl[nc + 2] = Cpl[nc + 1];
		Cd.resize(nc + 3); Cd[0] = 0.0; Cd[1] = 0.0;
		for (int i = 2; i < nc + 3; i++) {
			double xi = Ce[i];
			double si = Cs[i];
			double xipl = Cpl[i];
			double xipl_max = xi - si / E;
			xipl = std::min(xipl, xipl_max);
			double qi = (xi - xipl) * E;
			Cd[i] = 1.0 - si / qi;
		}
	}

	// check lists
	if (Te.size() < 1) {
		opserr << "nDMaterial ASDConcrete1D Error: 'Te' list is empty. At least 1 non-zero value should be provided.\n";
		return nullptr;
	}
	if (Ts.size() != Te.size()) {
		opserr << "nDMaterial ASDConcrete1D Error: 'Te' (size = " <<
			static_cast<int>(Te.size()) << ") and 'Ts' (size = " <<
			static_cast<int>(Ts.size()) << ") lists should have the same size.\n";
		return nullptr;
	}
	if (Td.size() == 0) {
		Td.resize(Te.size(), 0.0);
	}
	else if (Td.size() != Te.size()) {
		opserr << "nDMaterial ASDConcrete1D Error: 'Te' (size = " <<
			static_cast<int>(Te.size()) << ") and 'Td' (size = " <<
			static_cast<int>(Td.size()) << ") lists should have the same size.\n";
		return nullptr;
	}
	if (Ce.size() < 1) {
		opserr << "nDMaterial ASDConcrete1D Error: 'Tc' list is empty. At least 1 non-zero value should be provided.\n";
		return nullptr;
	}
	if (Cs.size() != Ce.size()) {
		opserr << "nDMaterial ASDConcrete1D Error: 'Ce' (size = " <<
			static_cast<int>(Ce.size()) << ") and 'Cs' (size = " <<
			static_cast<int>(Cs.size()) << ") lists should have the same size.\n";
		return nullptr;
	}
	if (Cd.size() == 0) {
		Cd.resize(Ce.size(), 0.0);
	}
	else if (Cd.size() != Ce.size()) {
		opserr << "nDMaterial ASDConcrete1D Error: 'Ce' (size = " <<
			static_cast<int>(Ce.size()) << ") and 'Cd' (size = " <<
			static_cast<int>(Cd.size()) << ") lists should have the same size.\n";
		return nullptr;
	}

	// build the hardening laws
	ASDConcrete1DMaterial::HardeningLaw HT(tag, ASDConcrete1DMaterial::HardeningLawType::Tension, E, Te, Ts, Td);
	if (!HT.isValid()) {
		opserr << "nDMaterial ASDConcrete1D Error: Tensile hardening law is not valid.\n";
		return nullptr;
	}
	ASDConcrete1DMaterial::HardeningLaw HC(tag, ASDConcrete1DMaterial::HardeningLawType::Compression, E, Ce, Cs, Cd);
	if (!HC.isValid()) {
		opserr << "nDMaterial ASDConcrete1D Error: Compressive hardening law is not valid.\n";
		return nullptr;
	}

	// create the material
	UniaxialMaterial* instance = new ASDConcrete1DMaterial(
		tag,
		E, eta,
		implex, implex_control, implex_abort_on_error, implex_error_tolerance, implex_time_redution_limit, implex_alpha,
		tangent, auto_regularization, lch_ref,
		HT, HC);
	if (instance == nullptr) {
		opserr << "UniaxialMaterial ASDConcrete1D Error: failed to allocate a new material.\n";
		return nullptr;
	}
	return instance;
}

ASDConcrete1DMaterial::ASDConcrete1DMaterial(
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
	const HardeningLaw& _hc)
	: UniaxialMaterial(_tag, MAT_TAG_ASDConcrete1DMaterial)
	, E(_E)
	, eta(_eta)
	, implex(_implex)
	, implex_control(_implex_control)
	, implex_abort_on_error(_implex_abort_on_error)
	, implex_error_tolerance(_implex_error_tolerance)
	, implex_time_redution_limit(_implex_time_reduction_limit)
	, implex_alpha(_implex_alpha)
	, tangent(_tangent)
	, auto_regularize(_auto_regularize)
	, lch_ref(_lch_ref)
	, ht(_ht)
	, hc(_hc)
{
	// intialize C as C0
	C = getInitialTangent();

	// initialize PT_commit as eye(6)*0.5
	PT_commit = 0.5;
}

ASDConcrete1DMaterial::ASDConcrete1DMaterial()
	: UniaxialMaterial(0, MAT_TAG_ASDConcrete1DMaterial)
{
}

ASDConcrete1DMaterial::~ASDConcrete1DMaterial()
{
}

int ASDConcrete1DMaterial::setTrialStrain(double v, double r)
{
	// return value
	int retval = 0;

	// this material point takes part in the current step, so it takes part in
	// the IMPL-EX error aggregate (see IMPLEXManager.h)
	implexTouch();

	// get characteristic length and perform regularization
	if (!regularization_done) {
		if (ops_TheActiveElement)
			lch = ops_TheActiveElement->getCharacteristicLength();
		regularization_done = true;
		if (auto_regularize) {
			ht.regularize(lch, lch_ref);
			hc.regularize(lch, lch_ref);
		}
	}

	// save dT
	if (!dtime_is_user_defined) {
		dtime_n = ops_Dt;
		if (!commit_done) {
			dtime_0 = dtime_n;
			dtime_n_commit = dtime_n;
		}
	}

	// if the user requested the tangent matrix
	// and not the IMPL-EX (in IMPL-EX the tangent coincides with the secant) ...
	if (tangent && !implex) {
		// numerical tangent tensor
		double Cnum = 0.0;
		// strain perturbation parameter
		double PERT = (ht.strainTolerance() + hc.strainTolerance()) / 2.0;
		// compute the forward perturbed solution and store in Cnum
		strain = v + PERT;
		retval = compute(true, false);
		if (retval < 0) return retval;
		Cnum = stress;
		// compute unperturbed solution
		strain = v;
		retval = compute(true, false);
		if (retval < 0) return retval;
		Cnum = (Cnum - stress) / PERT;
		// save tangent
		C = Cnum;
	}
	else {
		strain = v;
		if (implex) {
			if (implex_control) {
				// LEGACY in-material measurement: an implicit solution first,
				// to have something to compare against, then the extrapolated
				// one, which is the answer that must reach the element.
				// The error control does not need this: the convergence test
				// wrapper measures once per step, through
				// computeImplexErrorMetric(), instead of once per iteration
				double aux = PT_commit;
				retval = compute(false, false);
				if (retval < 0) return retval;
				double stress_implicit = stress;
				double dt_implicit = dt_bar;
				double dc_implicit = dc_bar;
				// the implicit pass has overwritten the frozen switch: put it
				// back before extrapolating, or the extrapolation is not the
				// one this step would have taken. Nothing else needs restoring
				// here, because compute() always restarts from the committed
				// state
				PT_commit = aux;
				// standard call
				retval = compute(true, true);
				if (retval < 0) return retval;
				implex_error = implexStressGap(stress, stress_implicit);
				implex_gap = implexDamageGap(dt_bar, dc_bar, dt_implicit, dc_implicit);
				// and only if the user asked for the old behaviour, fail here
				if (implex_abort_on_error && implex_error > implex_error_tolerance) {
					if (dtime_n >= implex_time_redution_limit * dtime_0) {
						retval = EC_IMPLEX_Error_Control;
					}
				}
			}
			else {
				// standard call 
				retval = compute(true, true);
			}
		}
		else {
			// standard call
			retval = compute(true, true);
		}
	}

	// RECORD WHAT THIS STEP DELIVERS, in ONE place covering all four branches
	// above and not once per branch: `stress` is what getStress() is about to
	// hand the element on every path, extrapolated or implicit, and commitState()
	// is about to overwrite it with the implicit re-solve. Writing this per
	// branch is how the same member came to be missed in the numeric-tangent
	// branch of ASDHysteretic1D - the one branch that is easy to forget
	stress_implex = stress;

	// done
	return retval;
}

double ASDConcrete1DMaterial::getStress(void)
{
	return stress;
}

double ASDConcrete1DMaterial::getTangent(void)
{
	return C;
}

double ASDConcrete1DMaterial::getInitialTangent(void)
{
	return E;
}

double ASDConcrete1DMaterial::getStrain(void)
{
	return strain;
}

int ASDConcrete1DMaterial::commitState(void)
{
	// implicit stage
	if (implex) {
		// what the extrapolated step delivered is already in stress_implex,
		// recorded by setTrialStrain, and it is now also readable from outside as
		// the 'implexStress' response
		double dt = dt_bar;
		double dc = dc_bar;
		// implicit solution. Note that this is not a second answer that gets
		// thrown away: the implicit one IS what gets committed, so this same
		// call measures the error and does the first half of the commit
		compute(false, false);
		// compute implex error here always
		implex_error = implexStressGap(stress_implex, stress);
		implex_gap = implexDamageGap(dt, dc, dt_bar, dc_bar);
		GlobalParameters::instance().setMaxError(std::max(implex_error, GlobalParameters::instance().getMaxError()));
		GlobalParameters::instance().accumulateAverageError(implex_error);
	}
	// compute energy
	energy += 0.5 * (stress_commit + stress) * (strain - strain_commit);
	// store the previously committed variables for next move from n to n - 1
	xt_commit_old = xt_commit;
	xc_commit_old = xc_commit;
	// store committed variables
	xt_commit = xt;
	xc_commit = xc;
	strain_commit = strain;
	stress_commit = stress;
	stress_eff_commit = stress_eff;
	dtime_n_commit = dtime_n;
	// done
	commit_done = true;
	return 0;
}

int ASDConcrete1DMaterial::revertToLastCommit(void)
{
	// restore converged values
	xt = xt_commit;
	xc = xc_commit;
	strain = strain_commit;
	stress = stress_commit;
	stress_eff = stress_eff_commit;
	dtime_n = dtime_n_commit;
	// done
	return 0;
}

int ASDConcrete1DMaterial::revertToStart(void)
{
	// State variables
	xt = 0.0;
	xt_commit = 0.0;
	xt_commit_old = 0.0;
	xc = 0.0;
	xc_commit = 0.0;
	xc_commit_old = 0.0;

	// Time step
	dtime_n = 0.0;
	dtime_n_commit = 0.0;
	dtime_0 = 0.0;
	dtime_is_user_defined = false;

	// Commit flag
	commit_done = false;

	// IMPL-EX error
	implex_error = 0.0;
	implex_gap = 0.0;

	// Strain, Stress and Tangent
	strain = 0.0;
	strain_commit = 0.0;
	stress = 0.0;
	stress_commit = 0.0;
	stress_implex = 0.0;
	stress_eff = 0.0;
	stress_eff_commit = 0.0;
	C = getInitialTangent();
	PT_commit = 0.5;

	// Output variables
	dt_bar = 0.0;
	dc_bar = 0.0;
	energy = 0.0;

	// Done
	return 0;
}

UniaxialMaterial* ASDConcrete1DMaterial::getCopy(void)
{
	// we can safely use the default copy-constructor according to the member variables we're using
	return new ASDConcrete1DMaterial(*this);
}

void ASDConcrete1DMaterial::Print(OPS_Stream& s, int flag)
{
	s << "ASDConcrete1D Material, tag: " << this->getTag() << "\n";
}

int ASDConcrete1DMaterial::sendSelf(int commitTag, Channel &theChannel)
{
	// aux
	int counter;

	// variable DBL data size (30: stress_implex was added to the state)
	int nv_dbl = 30 +
		ht.serializationDataSize() +
		hc.serializationDataSize();

	// send INT data
	static ID idata(10);
	counter = 0;
	idata(counter++) = getTag();
	idata(counter++) = static_cast<int>(implex);
	idata(counter++) = static_cast<int>(implex_control);
	idata(counter++) = static_cast<int>(implex_abort_on_error);
	idata(counter++) = static_cast<int>(tangent);
	idata(counter++) = static_cast<int>(auto_regularize);
	idata(counter++) = static_cast<int>(regularization_done);
	idata(counter++) = static_cast<int>(dtime_is_user_defined);
	idata(counter++) = static_cast<int>(commit_done);
	idata(counter++) = nv_dbl;
	if (theChannel.sendID(getDbTag(), commitTag, idata) < 0) {
		opserr << "ASDConcrete1DMaterial::sendSelf() - failed to send INT data\n";
		return -1;
	}

	// send DBL data
	Vector ddata(nv_dbl);
	counter = 0;
	ddata(counter++) = E;
	ddata(counter++) = eta;
	ddata(counter++) = implex_error_tolerance;
	ddata(counter++) = implex_time_redution_limit;
	ddata(counter++) = implex_alpha;
	ddata(counter++) = lch;
	ddata(counter++) = lch_ref;
	ddata(counter++) = xt;
	ddata(counter++) = xt_commit;
	ddata(counter++) = xt_commit_old;
	ddata(counter++) = xc;
	ddata(counter++) = xc_commit;
	ddata(counter++) = xc_commit_old;
	ddata(counter++) = dtime_n;
	ddata(counter++) = dtime_n_commit;
	ddata(counter++) = dtime_0;
	ddata(counter++) = implex_error;
	ddata(counter++) = implex_gap;
	ddata(counter++) = PT_commit;
	ddata(counter++) = strain;
	ddata(counter++) = strain_commit;
	ddata(counter++) = stress;
	ddata(counter++) = stress_commit;
	ddata(counter++) = stress_implex;
	ddata(counter++) = stress_eff;
	ddata(counter++) = stress_eff_commit;
	ddata(counter++) = C;
	ddata(counter++) = dt_bar;
	ddata(counter++) = dc_bar;
	ddata(counter++) = energy;
	ht.serialize(ddata, counter);
	hc.serialize(ddata, counter);
	if (theChannel.sendVector(getDbTag(), commitTag, ddata) < 0) {
		opserr << "ASDConcrete1DMaterial::sendSelf() - failed to send DBL data\n";
		return -1;
	}

	// done
	return 0;
}

int ASDConcrete1DMaterial::recvSelf(int commitTag, Channel& theChannel, FEM_ObjectBroker& theBroker)
{
	// aux
	int counter;

	// recv INT data
	static ID idata(10);
	if (theChannel.recvID(getDbTag(), commitTag, idata) < 0) {
		opserr << "ASDConcrete1DMaterial::recvSelf() - failed to receive INT data\n";
		return -1;
	}
	counter = 0;
	setTag(idata(counter++));
	implex = static_cast<bool>(idata(counter++));
	implex_control = static_cast<bool>(idata(counter++));
	implex_abort_on_error = static_cast<bool>(idata(counter++));
	tangent = static_cast<bool>(idata(counter++));
	auto_regularize = static_cast<bool>(idata(counter++));
	regularization_done = static_cast<bool>(idata(counter++));
	dtime_is_user_defined = static_cast<bool>(idata(counter++));
	commit_done = static_cast<bool>(idata(counter++));
	int nv_dbl = idata(counter++);

	// recv DBL data
	Vector ddata(nv_dbl);
	if (theChannel.recvVector(getDbTag(), commitTag, ddata) < 0) {
		opserr << "ASDConcrete1DMaterial::recvSelf() - failed to receive DBL data\n";
		return -1;
	}
	counter = 0;
	E = ddata(counter++);
	eta = ddata(counter++);
	implex_error_tolerance = ddata(counter++);
	implex_time_redution_limit = ddata(counter++);
	implex_alpha = ddata(counter++);
	lch = ddata(counter++);
	lch_ref = ddata(counter++);
	xt = ddata(counter++);
	xt_commit = ddata(counter++);
	xt_commit_old = ddata(counter++);
	xc = ddata(counter++);
	xc_commit = ddata(counter++);
	xc_commit_old = ddata(counter++);
	dtime_n = ddata(counter++);
	dtime_n_commit = ddata(counter++);
	dtime_0 = ddata(counter++);
	implex_error = ddata(counter++);
	implex_gap = ddata(counter++);
	PT_commit = ddata(counter++);
	strain = ddata(counter++);
	strain_commit = ddata(counter++);
	stress = ddata(counter++);
	stress_commit = ddata(counter++);
	stress_implex = ddata(counter++);
	stress_eff = ddata(counter++);
	stress_eff_commit = ddata(counter++);
	C = ddata(counter++);
	dt_bar = ddata(counter++);
	dc_bar = ddata(counter++);
	energy = ddata(counter++);
	ht.deserialize(ddata, counter);
	hc.deserialize(ddata, counter);

	// done
	return 0;
}

int ASDConcrete1DMaterial::setParameter(const char** argv, int argc, Parameter& param)
{
	
	// 1000 - elasticity & mass & length
	if (strcmp(argv[0], "E") == 0) {
		param.setValue(E);
		return param.addObject(1000, this);
	}
	if (strcmp(argv[0], "lch_ref") == 0) {
		param.setValue(lch_ref);
		return param.addObject(1001, this);
	}

	// 2000 - time
	if (strcmp(argv[0], "dTime") == 0) {
		param.setValue(dtime_n);
		return param.addObject(2000, this);
	}
	if (strcmp(argv[0], "dTimeCommit") == 0) {
		param.setValue(dtime_n_commit);
		return param.addObject(2001, this);
	}
	if (strcmp(argv[0], "dTimeInitial") == 0) {
		param.setValue(dtime_0);
		return param.addObject(2002, this);
	}

	// 3000 - globals
	if (strcmp(argv[0], "implexError") == 0 || strcmp(argv[0], "ImplexError") == 0) {
		param.setValue(GlobalParameters::instance().getMaxError());
		return param.addObject(3000, this);
	}
	if (strcmp(argv[0], "avgImplexError") == 0 || strcmp(argv[0], "AvgImplexError") == 0) {
		param.setValue(GlobalParameters::instance().getAverageError());
		return param.addObject(3001, this);
	}

	// default
	return -1;
}

int ASDConcrete1DMaterial::updateParameter(int parameterID, Information& info)
{
	switch (parameterID) {
		// 1000 - elasticity & mass
	case 1000:
		E = info.theDouble;
		return 0;
	case 1001:
		lch_ref = info.theDouble;
		auto_regularize = true;
		regularization_done = false;
		return 0;

		// 2000 - time
	case 2000:
		dtime_n = info.theDouble;
		dtime_is_user_defined = true;
		return 0;
	case 2001:
		dtime_n_commit = info.theDouble;
		dtime_is_user_defined = true;
		return 0;
	case 2002:
		dtime_0 = info.theDouble;
		dtime_is_user_defined = true;
		return 0;

	case 3000:
		GlobalParameters::instance().setMaxError(info.theDouble);
		return 0;
	case 3001:
		GlobalParameters::instance().setAverageError(info.theDouble);
		return 0;

		// default
	default:
		return -1;
	}
}

Response* ASDConcrete1DMaterial::setResponse(const char** argv, int argc, OPS_Stream& output)
{
	// utils
	auto make_resp = [&output, this](int rid, const Vector& v, const std::vector<std::string>* labels = nullptr) -> MaterialResponse* {
		output.tag("UniaxialMaterialOutput");
		output.attr("matType", getClassType());
		output.attr("matTag", getTag());
		if (labels) {
			for (const auto& item : (*labels))
				output.tag("ResponseType", item.c_str());
		}
		MaterialResponse* resp = new MaterialResponse(this, rid, v);
		output.endTag();
		return resp;
	};

	// labels
	static std::vector<std::string> lb_damage = { "d+", "d-" };
	static std::vector<std::string> lb_eqpl_strain = { "PLE+", "PLE-" };
	static std::vector<std::string> lb_tot_strain = { "TE+", "TE-" };
	static std::vector<std::string> lb_cw = { "cw" };
	static std::vector<std::string> lb_crackpattern = { "C1x", "C1y", "C1z",   "C2x", "C2y", "C2z",   "C3x", "C3y", "C3z" };
	static std::vector<std::string> lb_implex_error = { "Error" };
	static std::vector<std::string> lb_implex_stress = { "S" };
	static std::vector<std::string> lb_implex_gap = { "Gap" };
	static std::vector<std::string> lb_time = { "dTime", "dTimeCommit", "dTimeInitial" };
	static std::vector<std::string> lb_crack_strain = { "CS+", "LchRef" };
	static std::vector<std::string> lb_crush_strain = { "CS-", "LchRef" };
	static Vector Cinfo(2);

	// check specific responses
	if (argc > 0) {
		// 1000 - compressive hardening variables
		if (strcmp(argv[0], "Ce") == 0)
			return make_resp(1000, getHardeningLawVector(HardeningLawType::Compression, HardeningLawPointComponent::TotalStrain));
		if (strcmp(argv[0], "Cs") == 0)
			return make_resp(1001, getHardeningLawVector(HardeningLawType::Compression, HardeningLawPointComponent::NominalStress));
		if (strcmp(argv[0], "Cq") == 0)
			return make_resp(1002, getHardeningLawVector(HardeningLawType::Compression, HardeningLawPointComponent::EffectiveStress));
		// 1100 - tensile hardening variables
		if (strcmp(argv[0], "Te") == 0)
			return make_resp(1100, getHardeningLawVector(HardeningLawType::Tension, HardeningLawPointComponent::TotalStrain));
		if (strcmp(argv[0], "Ts") == 0)
			return make_resp(1101, getHardeningLawVector(HardeningLawType::Tension, HardeningLawPointComponent::NominalStress));
		if (strcmp(argv[0], "Tq") == 0)
			return make_resp(1102, getHardeningLawVector(HardeningLawType::Tension, HardeningLawPointComponent::EffectiveStress));
		// 2000 - damage variables
		if (strcmp(argv[0], "damage") == 0 || strcmp(argv[0], "Damage") == 0)
			return make_resp(2000, getDamage(), &lb_damage);
		if (strcmp(argv[0], "equivalentPlasticStrain") == 0 || strcmp(argv[0], "EquivalentPlasticStrain") == 0)
			return make_resp(2001, getEquivalentPlasticStrain(), &lb_eqpl_strain);
		if (strcmp(argv[0], "equivalentTotalStrain") == 0 || strcmp(argv[0], "EquivalentTotalStrain") == 0)
			return make_resp(2002, getStrainMeasure(), &lb_tot_strain);
		if (strcmp(argv[0], "cw") == 0 || strcmp(argv[0], "crackWidth") == 0 || strcmp(argv[0], "CrackWidth") == 0)
			return make_resp(2003, getCrackWidth(), &lb_cw);
		if (strcmp(argv[0], "crackStrain") == 0 || strcmp(argv[0], "CrackStrain") == 0) {
			double user_lch_ref = this->lch_ref; // by default the input one
			if (argc > 2 && strcmp(argv[1], "-lchRef") == 0) {
				double trial_lch_ref = 0.0;
				if (string_to_double(argv[2], trial_lch_ref)) {
					user_lch_ref = trial_lch_ref;
				}
			}
			Cinfo(1) = user_lch_ref;
			Cinfo(0) = getCrackWidth()(0) / user_lch_ref;
			return make_resp(2004, Cinfo, &lb_crack_strain);
		}
		if (strcmp(argv[0], "crushStrain") == 0 || strcmp(argv[0], "CrushStrain") == 0) {
			double user_lch_ref = this->lch_ref; // by default the input one
			if (argc > 2 && strcmp(argv[1], "-lchRef") == 0) {
				double trial_lch_ref = 0.0;
				if (string_to_double(argv[2], trial_lch_ref)) {
					user_lch_ref = trial_lch_ref;
				}
			}
			Cinfo(1) = user_lch_ref;
			Cinfo(0) = getCrushWidth()(0) / user_lch_ref;
			return make_resp(2005, Cinfo, &lb_crush_strain);
		}
		// 3000 - implex error
		if (strcmp(argv[0], "implexError") == 0 || strcmp(argv[0], "ImplexError") == 0) {
			return make_resp(3000, getImplexError(), &lb_implex_error);
		}
		// 3001 - implex legacy damage gap (output only, the model does not
		// steer on it: see implexDamageGap)
		if (strcmp(argv[0], "implexGap") == 0 || strcmp(argv[0], "ImplexGap") == 0) {
			return make_resp(3001, getImplexGap(), &lb_implex_gap);
		}
		// 3003 - THE STRESS THIS STEP DELIVERED, which under IMPL-EX is not the
		// one 'stress' reports afterwards: commitState installs the implicit
		// solution over it. Same id and same name in every material of the family
		if (strcmp(argv[0], "implexStress") == 0 || strcmp(argv[0], "ImplexStress") == 0) {
			return make_resp(3003, getImplexStress(), &lb_implex_stress);
		}
		// 4000 - internal time
		if (strcmp(argv[0], "time") == 0 || strcmp(argv[0], "Time") == 0) {
			return make_resp(4000, getTimeIncrements(), &lb_time);
		}
	}

	// otherwise return base-class response
	return UniaxialMaterial::setResponse(argv, argc, output);
}

int ASDConcrete1DMaterial::getResponse(int responseID, Information& matInformation)
{
	switch (responseID) {
		// 1000 - compressive hardening variables
	case 1000: return matInformation.setVector(getHardeningLawVector(HardeningLawType::Compression, HardeningLawPointComponent::TotalStrain));
	case 1001: return matInformation.setVector(getHardeningLawVector(HardeningLawType::Compression, HardeningLawPointComponent::NominalStress));
	case 1002: return matInformation.setVector(getHardeningLawVector(HardeningLawType::Compression, HardeningLawPointComponent::EffectiveStress));
		// 1100 - tensile hardening variables
	case 1100: return matInformation.setVector(getHardeningLawVector(HardeningLawType::Tension, HardeningLawPointComponent::TotalStrain));
	case 1101: return matInformation.setVector(getHardeningLawVector(HardeningLawType::Tension, HardeningLawPointComponent::NominalStress));
	case 1102: return matInformation.setVector(getHardeningLawVector(HardeningLawType::Tension, HardeningLawPointComponent::EffectiveStress));
		// 2000 - damage variables
	case 2000: return matInformation.setVector(getDamage());
	case 2001: return matInformation.setVector(getEquivalentPlasticStrain());
	case 2002: return matInformation.setVector(getStrainMeasure());
	case 2003: return matInformation.setVector(getCrackWidth());
	case 2004:
		if (matInformation.theVector && matInformation.theVector->Size() == 2) {
			double user_lch_ref = matInformation.theVector->operator()(1);
			matInformation.theVector->operator()(0) = getCrackWidth()(0) / user_lch_ref;
			return 0;
		}
		break;
	case 2005:
		if (matInformation.theVector && matInformation.theVector->Size() == 2) {
			double user_lch_ref = matInformation.theVector->operator()(1);
			matInformation.theVector->operator()(0) = getCrushWidth()(0) / user_lch_ref;
			return 0;
		}
		break;
		// 3000 - implex error
	case 3000: return matInformation.setVector(getImplexError());
		// 3001 - implex legacy damage gap
	case 3001: return matInformation.setVector(getImplexGap());
		// 3003 - the stress this step delivered
	case 3003: return matInformation.setVector(getImplexStress());
		// 4000 - internal time
	case 4000: return matInformation.setVector(getTimeIncrements());
	default:
		break;
	}
	return UniaxialMaterial::getResponse(responseID, matInformation);
}

double ASDConcrete1DMaterial::getEnergy(void)
{
	return energy;
}

int ASDConcrete1DMaterial::compute(bool do_implex, bool do_tangent)
{
	// get committed variables
	xt = xt_commit;
	xc = xc_commit;
	stress = stress_commit;
	stress_eff = stress_eff_commit;

	// time factor for explicit extrapolation
	double time_factor = 1.0;
	if (implex && do_implex && (dtime_n_commit > 0.0))
		time_factor = dtime_n / dtime_n_commit * implex_alpha;

	// compute rate coefficients
	double rate_coeff_1 = 0.0;
	double rate_coeff_2 = 1.0;
	if ((dtime_n > 0.0) && (eta > 0.0)) {
		rate_coeff_1 = eta / (eta + dtime_n);
		rate_coeff_2 = dtime_n / (eta + dtime_n);
	}

	// compute elastic effective stress: SEFFn = C0 : (En - En-1)
	double dStrain = strain - strain_commit;
	stress_eff += E * dStrain;

	// compute stress split
	double PT = (implex && do_implex) ? PT_commit : Heavyside(stress_eff);
	double PC = 1.0 - PT;
	double ST = PT * stress_eff;
	double SC = PC * stress_eff;

	// compute committed hardening variables
	HardeningLawPoint pt = ht.evaluateAt(xt);
	HardeningLawPoint pc = hc.evaluateAt(xc);

	// temporary clone of old equivalent plastic strains
	double xt_pl = pt.plasticStrain(E);
	double xc_pl = pc.plasticStrain(E);

	// compute new trial equivalent strain measures
	if (implex && do_implex) {
		// extrapolated equivalent strain measures (explicit)
		xt = xt_commit + time_factor * (xt_commit - xt_commit_old);
		xc = xc_commit + time_factor * (xc_commit - xc_commit_old);
	}
	else {
		// compute trial strain measures (implicit)
		double xt_trial = ST / E + xt_pl;
		double xc_trial = -SC / E + xc_pl;
		// update hardening variables
		if (xt_trial > pt.x) 
			xt = rate_coeff_1 * pt.x + rate_coeff_2 * xt_trial;
		if (xc_trial > pc.x) 
			xc = rate_coeff_1 * pc.x + rate_coeff_2 * xc_trial;
	}
	pt = ht.evaluateAt(xt);
	pc = hc.evaluateAt(xc);

	// compute plastic damage
	double seff_eq_t = (pt.x - xt_pl) * E;
	double dt_plastic = seff_eq_t > 0.0 ? 1.0 - pt.q / seff_eq_t : 0.0;
	double seff_eq_c = (pc.x - xc_pl) * E;
	double dc_plastic = seff_eq_c > 0.0 ? 1.0 - pc.q / seff_eq_c : 0.0;

	// update effective stress
	stress_eff = (1.0 - dt_plastic) * ST + (1 - dc_plastic) * SC;

	// update nominal stress
	dt_bar = pt.d + dt_plastic - pt.d * dt_plastic;
	dc_bar = pc.d + dc_plastic - pc.d * dc_plastic;
	stress = (1.0 - dt_bar) * ST + (1.0 - dc_bar) * SC;

	// secant matrix
	if (do_tangent) {
		double W = 1.0 - dt_bar * PT - dc_bar * PC;
		C = W * E;
	}

	// save real PT and R, if mp.implex and !do_implex -> called from commit
	if (implex && !do_implex) {
		// save it in implex mode during implicit phase
		PT_commit = PT;
	}

	// done
	return 0;
}

Vector ASDConcrete1DMaterial::getHardeningLawVector(HardeningLawType ltype, HardeningLawPointComponent c) const
{
	Vector r;
	const HardeningLaw& h = ltype == HardeningLawType::Tension ? ht : hc;
	r.resize(static_cast<int>(h.points().size()));
	for (std::size_t i = 0; i < h.points().size(); ++i) {
		const HardeningLawPoint& p = h.points()[i];
		switch (c)
		{
		case ASDConcrete1DMaterial::HardeningLawPointComponent::TotalStrain:
			r(static_cast<int>(i)) = p.totalStrain();
			break;
		case ASDConcrete1DMaterial::HardeningLawPointComponent::EffectiveStress:
			r(static_cast<int>(i)) = p.effectiveStress();
			break;
		case ASDConcrete1DMaterial::HardeningLawPointComponent::NominalStress:
			r(static_cast<int>(i)) = p.stress();
			break;
		default:
			break;
		}
	}
	return r;
}

const Vector& ASDConcrete1DMaterial::getStrainMeasure() const
{
	static Vector d(2);
	d(0) = xt;
	d(1) = xc;
	return d;
}

const Vector& ASDConcrete1DMaterial::getDamage() const
{
	static Vector d(2);
	const Vector& x = getStrainMeasure();
	d(0) = ht.evaluateAt(x(0)).crackingDamage();
	d(1) = hc.evaluateAt(x(1)).crackingDamage();
	return d;
}

const Vector& ASDConcrete1DMaterial::getEquivalentPlasticStrain() const
{
	static Vector d(2);
	const Vector& x = getStrainMeasure();
	d(0) = ht.evaluateAt(x(0)).plasticStrain(E);
	d(1) = hc.evaluateAt(x(1)).plasticStrain(E);
	return d;
}

const Vector& ASDConcrete1DMaterial::getCrackWidth() const
{
	static Vector d(1);
	d.Zero();
	if (ht.hasStrainSoftening()) {
		double e0 = ht.strainAtOnsetOfCrack();
		const Vector& x = getStrainMeasure();
		d(0) = std::max(x(0) - e0, 0.0) * lch;
	}
	return d;
}

const Vector& ASDConcrete1DMaterial::getCrushWidth() const
{
	static Vector d(1);
	d.Zero();
	if (hc.hasStrainSoftening()) {
		double e0 = hc.strainAtOnsetOfCrack();
		const Vector& x = getStrainMeasure();
		d(0) = std::max(x(1) - e0, 0.0) * lch;
	}
	return d;
}

const Vector& ASDConcrete1DMaterial::getImplexError() const
{
	static Vector d(1);
	d(0) = implex_error;
	return d;
}

const Vector& ASDConcrete1DMaterial::getImplexStress() const
{
	static Vector d(1);
	d(0) = stress_implex;
	return d;
}

const Vector& ASDConcrete1DMaterial::getImplexGap() const
{
	static Vector d(1);
	d(0) = implex_gap;
	return d;
}

double ASDConcrete1DMaterial::computeImplexErrorMetric(void)
{
	// no extrapolation, no error. This is what makes it safe for the
	// aggregation to call it on everything it has
	if (!implex)
		return 0.0;
	// the current state holds the EXPLICIT answer: the stress this step
	// delivered to the element, and the state the recorders must keep seeing
	TrialState delivered;
	saveTrialState(delivered);
	// stress_implex already holds it - setTrialStrain records it - and compute()
	// below does not write that member, so it survives the throw-away solve.
	// READ here, never written: measuring is not allowed to move the state
	// the implicit answer, at the same trial strain
	if (compute(false, false) < 0) {
		restoreTrialState(delivered);
		// no metric is not a small metric: whoever reads this must not accept
		// the step
		return std::numeric_limits<double>::quiet_NaN();
	}
	double err = implexStressGap(stress_implex, stress);
	implex_gap = implexDamageGap(delivered.dt_bar, delivered.dc_bar, dt_bar, dc_bar);
	// undo. Measuring is not allowed to move the state: the step may still be
	// rejected, and revertToLastCommit() would not put PT_commit back
	restoreTrialState(delivered);
	implex_error = err;
	return err;
}

double ASDConcrete1DMaterial::implexTimeRatio(void) const
{
	return dtime_0 > 0.0 ? dtime_n / dtime_0 : 1.0;
}

double ASDConcrete1DMaterial::implexStressGap(double delivered, double stress_implicit) const
{
	// the reference is the largest stress this material can carry. The NORM
	// of the gap may differ from model to model (an absolute value here, the
	// largest Voigt component in 3D); the DENOMINATOR must not, or the same
	// tolerance means different things in 1D and in 3D and a stepper ends up
	// comparing numbers that are not comparable
	double ref = std::max(ht.computeMaxStress(), hc.computeMaxStress());
	if (ref <= 0.0)
		ref = 1.0;
	return std::abs(delivered - stress_implicit) / ref;
}

double ASDConcrete1DMaterial::implexDamageGap(double dt_implex, double dc_implex,
	double dt_implicit, double dc_implicit) const
{
	// LEGACY, for output only: see the note in the header on why the model
	// does not steer on this
	return std::max(std::abs(dt_implex - dt_implicit), std::abs(dc_implex - dc_implicit));
}

void ASDConcrete1DMaterial::saveTrialState(TrialState& x) const
{
	x.xt = xt;
	x.xc = xc;
	x.stress = stress;
	x.stress_eff = stress_eff;
	x.C = C;
	x.dt_bar = dt_bar;
	x.dc_bar = dc_bar;
	x.PT_commit = PT_commit;
}

void ASDConcrete1DMaterial::restoreTrialState(const TrialState& x)
{
	xt = x.xt;
	xc = x.xc;
	stress = x.stress;
	stress_eff = x.stress_eff;
	C = x.C;
	dt_bar = x.dt_bar;
	dc_bar = x.dc_bar;
	PT_commit = x.PT_commit;
}

const Vector& ASDConcrete1DMaterial::getTimeIncrements() const
{
	static Vector d(3);
	d(0) = dtime_n;
	d(1) = dtime_n_commit;
	d(2) = dtime_0;
	return d;
}

