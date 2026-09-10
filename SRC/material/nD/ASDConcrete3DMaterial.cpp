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
                                                                        
// $Revision: 1.4 $
// $Date: 2020-04-19 23:01:25 $
// $Source: /usr/local/cvs/OpenSees/SRC/material/nD/ASDConcrete3DMaterial.h,v $

// Massimo Petracca, Guido Camata - ASDEA Software, Italy
//
// A Simple and robust plastic-damage model for concrete and masonry
//

#include <ASDConcrete3DMaterial.h>
#include <ASDSpectralSplit.h>
#include <Channel.h>
#include <FEM_ObjectBroker.h>
#include <OPS_Globals.h>
#include <elementAPI.h>
#include <Element.h>
#include <MaterialResponse.h>
#include <Parameter.h>
#include <algorithm>
#include <limits>
#include <string>
#include <sstream>
#include <iomanip>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

#if defined(_WIN32)
#ifndef NOMINMAX
#define NOMINMAX
#endif
#endif

#define _FFMT(X) std::setw(10) << std::setprecision(3) << X

// anonymous namespace for utilities
namespace {

	enum ErrorCodes {
		EC_Generic = -1,
		EC_IMPLEX_Error_Control = -10,
		EC_Eigen_Error = -1000
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

	/*
	The 3x3 symmetric eigen-decomposition and the vj x vj x vj x vj projector used
	to live here, 193 lines of them, copied out of Matrix::Eigen3 to reach the
	eigenvectors that Matrix does not return. ASDPlasticDamageConcrete3D needs the
	same two things, and a second copy is one too many: they now live once, in
	ASDSpectralSplit.h. Read that file for the conventions (eigenvalues
	DESCENDING, the Voigt order, what the shear columns of pjj carry).

	normalize=false IS THE POINT OF THIS FORWARDER. The shared routine can
	normalize its input, which fixes a real defect - the sweep threshold
	`sm > 1.0e-8` is ABSOLUTE, so on a tensor smaller than that the loop never
	runs and the bare diagonal comes back with V = I. This material has always run
	without it and passing false keeps it bit-identical: measured on 48000 tensors
	over six scales from 1e+6 down to 1e-12, random and near-degenerate,
	eigenvalues eigenvectors and projectors all agree to 0.000e+00, every case.

	Turning it on here would NOT be free, and not only near zero: normalizing also
	makes that threshold relative, which at stress scale is LOOSER than the
	absolute one it replaces. Same measurement, normalize=true: the eigenvalues
	still agree to ~1e-15 relative but the eigenvectors move by up to 5e-8 at
	scale 1e+6 and 30. So it is a separate decision with its own verification.
	*/
	inline int Eigen3(Vector& d, Matrix& v)
	{
		return ASDSpectralSplit::eigen3(d, v, false);
	}

	using ASDSpectralSplit::computePjj;

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
    
    double D = B*B - 4.0*A*C;
    double t = (sqrt(D) - B) / (2.0*A);
    
    return (y0 - 2.0*y1 + y2)*t*t + 2.0*(y1 - y0)*t + y0;
  }

  
}

void *OPS_ASDConcrete3DMaterial(void)
{
	// some kudos
	static bool first_done = false;
	if (!first_done) {
		opserr << "Using ASDConcrete3D - Developed by: Massimo Petracca, Guido Camata, ASDEA Software Technology\n";
		first_done = true;
	}

	// check arguments
	int numArgs = OPS_GetNumRemainingInputArgs();
	if (numArgs < 3) {
		opserr << 
			"nDMaterial ASDConcrete3D Error: Few arguments (< 3).\n"
			"nDMaterial ASDConcrete3D $tag $E $v "
			"-Te $Te -Ts $Ts <-Td $Td> -Ce $Ce -Cs $Cs <-Cd $Cd> "
			"<-rho $rho> <-Kc $Kc>"
			"<-implex> <-implexControl $implexErrorTolerance $implexTimeReductionLimit> <-implexAbort> <-implexAlpha $alpha>"
			"<-cdf $cdf>"
			"<-crackPlanes $nct $ncc $smoothingAngle>"
			"<-eta $eta> <-tangent> <-autoRegularization $lch_ref>\n";
		return nullptr;
	}

	// numData
	int numData = 1;

	// data
	int tag;
	double E;
	double v;
	double rho = 0.0;
	bool implex = false;
	bool implex_control = false;
	bool implex_abort_on_error = false;
	double implex_error_tolerance = 0.05;
	double implex_time_redution_limit = 0.01;
	double implex_alpha = 1.0;
	double eta = 0.0;
	double Kc = 2.0 / 3.0; // default suggested by Lubliner et al.
	bool tangent = false;
	bool auto_regularization = false;
	double lch_ref = 1.0;
	std::vector<double> Te, Ts, Td, Ce, Cs, Cd;
	double cdf = 0.0;
	int nct = 0;
	int ncc = 0;
	double smoothing_angle = 45.0;
	
	// get tag
	if (OPS_GetInt(&numData, &tag) != 0)  {
		opserr << "nDMaterial ASDConcrete3D Error: invalid 'tag'.\n";
		return nullptr;
	}

	// get Elasticity arguments
	if (OPS_GetDouble(&numData, &E) != 0) {
		opserr << "nDMaterial ASDConcrete3D Error: invalid 'E'.\n";
		return nullptr;
	}
	if (E <= 0.0) {
		opserr << "nDMaterial ASDConcrete3D Error: invalid value for 'E' (" << E << "). It should be strictly positive.\n";
		return nullptr;
	}
	if (OPS_GetDouble(&numData, &v) != 0) {
		opserr << "nDMaterial ASDConcrete3D Error: invalid 'v'.\n";
		return nullptr;
	}

	// utilities (code re-use)
	auto lam_optional_int = [&numData](const char* variable, int& value) -> bool {
		if (OPS_GetNumRemainingInputArgs() > 0) {
			if (OPS_GetInt(&numData, &value) < 0) {
				opserr << "nDMaterial ASDConcrete3D Error: failed to get '" << variable << "'.\n";
				return false;
			}
		}
		else {
			opserr << "nDMaterial ASDConcrete3D Error: '" << variable << "' requested but not provided.\n";
			return false;
		}
		return true;
	};
	auto lam_optional_double = [&numData](const char* variable, double& value) -> bool {
		if (OPS_GetNumRemainingInputArgs() > 0) {
			if (OPS_GetDouble(&numData, &value) < 0) {
				opserr << "nDMaterial ASDConcrete3D Error: failed to get '" << variable << "'.\n";
				return false;
			}
		}
		else {
			opserr << "nDMaterial ASDConcrete3D Error: '" << variable << "' requested but not provided.\n";
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
				opserr << "nDMaterial ASDConcrete3D Error: cannot parse the '" << variable << "' list.\n";
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
		if (strcmp(value, "-rho") == 0) {
			if (!lam_optional_double("rho", rho))
				return nullptr;
		}
		else if (strcmp(value, "-fc") == 0) {
			if (!lam_optional_double("fc", fc))
				return nullptr;
			have_fc = true;
		}
		else if (strcmp(value, "-ft") == 0) {
			if (!lam_optional_double("ft", ft))
				return nullptr;
			have_ft = true;
		}		
		else if (strcmp(value, "-Kc") == 0) {
			if (!lam_optional_double("Kc", Kc))
				return nullptr;
			if (Kc < 2.0 / 3.0 || Kc > 1.0) {
				opserr << "nDMaterial ASDConcrete3D Error: 'Kc' (" << Kc << ") double be >= 2/3 and <= 1.\n";
				return nullptr;
			}
		}
		else if (strcmp(value, "-implex") == 0) {
			implex = true;
		}
		else if (strcmp(value, "-implexControl") == 0) {
			implex_control = true;
			if (OPS_GetNumRemainingInputArgs() < 2) {
				opserr << "nDMaterial ASDConcrete3D Error: '-implexControl' given without the next 2 arguments $implexErrorTolerance $implexTimeReductionLimit.\n";
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
				opserr << "nDMaterial ASDConcrete3D Error: '-autoRegularization' given without the next 1 argument $lch_ref.\n";
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
		else if (strcmp(value, "-crackPlanes") == 0) {
			if (OPS_GetNumRemainingInputArgs() < 3) {
				opserr << "nDMaterial ASDConcrete3D Error: '-crackPlanes' given without the next 3 arguments $nct $ncc and $smoothingAngle.\n";
				return nullptr;
			}
			if (!lam_optional_int("nct", nct))
				return nullptr;
			if (!lam_optional_int("ncc", ncc))
				return nullptr;
			if (!lam_optional_double("smoothingAngle", smoothing_angle))
				return nullptr;
		}
		else if (strcmp(value, "-cdf") == 0) {
			if (OPS_GetNumRemainingInputArgs() < 1) {
				opserr << "nDMaterial ASDConcrete3D Error: '-cdf' given without the next 1 argument $cdf.\n";
				return nullptr;
			}
			if (!lam_optional_double("cdf", cdf))
				return nullptr;
		}
	}

	// Set a default value of tension strength if none specified
	if (have_fc && !have_ft)
	  ft = 0.1*fc;

	if (have_fc) {
	  double ec = 2*fc/E;
	  double Gt = 0.073*pow(fc,0.18);
	  double Gc = 2*Gt*(fc*fc)/(ft*ft);
	  
	  
	  if (!have_lch_ref) {
	    //
	    // _get_lch_ref from ASDConcrete3D_MakeLaws.py
	    //
	    
	    // min lch for tension
	    double et_el = ft/E;
	    double Gt_min = 0.5*ft*et_el;
	    double hmin_t = 0.01*Gt/Gt_min;
	    
	    // min lch for compression
	    double ec1 = fc/E;
	    double ec_pl = (ec-ec1)*0.4 + ec1;
	    double Gc_min = 0.5*fc*(ec-ec_pl);
	    double hmin_c = 0.01*Gc/Gc_min;
	    
	    lch_ref = std::min(hmin_c,hmin_t);
	  }
	  
	  //
	  // _make_tension from ASDConcrete3D_MakeLaws.py
	  //

	  Gt = Gt/lch_ref;
	  
	  double f0 = 0.9*ft;
	  double f1 = ft;
	  double e0 = f0/E;
	  double e1 = 1.5*f1/E;
	  double ep = e1 - f1/E;
	  double f2 = 0.2*ft;
	  double f3 = 1.0e-3*ft;
	  double w2 = Gt/ft;
	  double w3 = 5.0*w2;
	  double e2 = w2 + f2/E + ep;
	  if (e2 <= e1)
	    e2 = 1.001*e1;
	  double e3 = w3 + f3/E + ep;
	  if (e3 <= e2)
	    e3 = 1.001*e2;  
	  double e4 = 10.0*e3;
	  Te.resize(6); Te = {0.0, e0, e1, e2, e3, e4};
	  Ts.resize(6); Ts = {0.0, f0, f1, f2, f3, f3};
	  Td.resize(6); Td = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
	  double Tpl[6] = {0.0, 0.0, ep, 0.9*e2, 0.8*e3, 0.8*e3};
	  
	  for (int i = 2; i < 6; i++) {
	    double xi = Te[i];
	    double si = Ts[i];
	    double xipl = Tpl[i];
	    double xipl_max = xi-si/E;
	    xipl = std::min(xipl, xipl_max);
	    double qi = (xi-xipl)*E;
	    Td[i] = 1.0 - si/qi;
	  }
	  
	  
	  
	  //
	  // _make_compression from ASDConcrete3D_MakeLaws.py
	  //

	  Gc = Gc/lch_ref;
	  
	  double fc0 = 0.5*fc;
	  double ec0 = fc0/E;
	  double ec1 = fc/E;
	  double fcr = 0.1*fc;
	  double ec_pl = (ec-ec1)*0.4 + ec1;
	  double Gc1 = 0.5*fc*(ec-ec_pl);
	  double Gc2 = std::max(0.01*Gc1,Gc-Gc1);
	  double ecr = ec + 2.0*Gc2/(fc+fcr);
	  const int nc = 10;
	  Ce.resize(nc+3); Ce[0] = 0.0; Ce[1] = ec0;
	  Cs.resize(nc+3); Cs[0] = 0.0; Cs[1] = fc0;
	  double Cpl[nc+3]; Cpl[0] = 0.0; Cpl[1] = 0.0;
	  double dec = (ec-ec0)/(nc-1);
	  for (int i = 0; i < nc-1; i++) {
	    double iec = ec0 + (i+1)*dec;
	    Ce[i+2] = iec;
	    Cs[i+2] = bezier3(iec,  ec0, ec1, ec,  fc0, fc, fc);
	    Cpl[i+2] = Cpl[i+1] + 0.7*(iec-Cpl[i+1]);
	  }
	  Ce[nc+1] = ecr;
	  Cs[nc+1] = fcr;
	  Cpl[nc+1] = Cpl[nc] + 0.7*(ecr-Cpl[nc]);
	  Ce[nc+2] = ecr + ec0;
	  Cs[nc+2] = fcr;
	  Cpl[nc+2] = Cpl[nc+1];
	  Cd.resize(nc+3); Cd[0] = 0.0; Cd[1] = 0.0;
	  for (int i = 2; i < nc+3; i++) {
	    double xi = Ce[i];
	    double si = Cs[i];
	    double xipl = Cpl[i];
	    double xipl_max = xi-si/E;
	    xipl = std::min(xipl, xipl_max);
	    double qi = (xi-xipl)*E;
	    Cd[i] = 1.0 - si/qi;
	  }
	}
	
	// check lists
	if (Te.size() < 1) {
		opserr << "nDMaterial ASDConcrete3D Error: 'Te' list is empty. At least 1 non-zero value should be provided.\n";
		return nullptr;
	}
	if (Ts.size() != Te.size()) {
		opserr << "nDMaterial ASDConcrete3D Error: 'Te' (size = " <<
			static_cast<int>(Te.size()) << ") and 'Ts' (size = " <<
			static_cast<int>(Ts.size()) << ") lists should have the same size.\n";
		return nullptr;
	}
	if (Td.size() == 0) {
		Td.resize(Te.size(), 0.0);
	}
	else if (Td.size() != Te.size()) {
		opserr << "nDMaterial ASDConcrete3D Error: 'Te' (size = " <<
			static_cast<int>(Te.size()) << ") and 'Td' (size = " <<
			static_cast<int>(Td.size()) << ") lists should have the same size.\n";
		return nullptr;
	}
	if (Ce.size() < 1) {
		opserr << "nDMaterial ASDConcrete3D Error: 'Tc' list is empty. At least 1 non-zero value should be provided.\n";
		return nullptr;
	}
	if (Cs.size() != Ce.size()) {
		opserr << "nDMaterial ASDConcrete3D Error: 'Ce' (size = " <<
			static_cast<int>(Ce.size()) << ") and 'Cs' (size = " <<
			static_cast<int>(Cs.size()) << ") lists should have the same size.\n";
		return nullptr;
	}
	if (Cd.size() == 0) {
		Cd.resize(Ce.size(), 0.0);
	}
	else if (Cd.size() != Ce.size()) {
		opserr << "nDMaterial ASDConcrete3D Error: 'Ce' (size = " <<
			static_cast<int>(Ce.size()) << ") and 'Cd' (size = " <<
			static_cast<int>(Cd.size()) << ") lists should have the same size.\n";
		return nullptr;
	}

	// build the hardening laws
	ASDConcrete3DMaterial::HardeningLaw HT(tag, ASDConcrete3DMaterial::HardeningLawType::Tension, E, Te, Ts, Td);
	if (!HT.isValid()) {
		opserr << "nDMaterial ASDConcrete3D Error: Tensile hardening law is not valid.\n";
		return nullptr;
	}
	ASDConcrete3DMaterial::HardeningLaw HC(tag, ASDConcrete3DMaterial::HardeningLawType::Compression, E, Ce, Cs, Cd);
	if (!HC.isValid()) {
		opserr << "nDMaterial ASDConcrete3D Error: Compressive hardening law is not valid.\n";
		return nullptr;
	}

	// create the material
	NDMaterial* instance = new ASDConcrete3DMaterial(
		tag, 
		E, v, rho, eta, Kc,
		implex, implex_control, implex_abort_on_error, implex_error_tolerance, implex_time_redution_limit, implex_alpha,
		tangent, auto_regularization, lch_ref,
		HT, HC,
		cdf, nct, ncc, smoothing_angle);
	if (instance == nullptr) {
		opserr << "nDMaterial ASDConcrete3D Error: failed to allocate a new material.\n";
		return nullptr;
	}
	return instance;
}

int ASDConcrete3DMaterial::StressDecomposition::compute(const Vector& S, double cdf)
{
	// copy S in V
	V(0, 0) = S(0);
	V(1, 1) = S(1);
	V(2, 2) = S(2);
	V(0, 1) = V(1, 0) = S(3);
	V(1, 2) = V(2, 1) = S(4);
	V(0, 2) = V(2, 0) = S(5);

	// solve
	if (Eigen3(Si, V) < 0)
		return EC_Eigen_Error;

	// construct matrices PT and PC
	static Matrix pjj(6, 6);
	PT.Zero();
	PC.Zero();

	// method 1: compute PT = sum(H(sj)*pjj) --- PC = I - PT
	//for (int j = 0; j < 3; j++) {
	//	double hsj = Heavyside(Si(j));
	//	if (hsj > 0.0) {
	//		computePjj(V, j, pjj);
	//		PT.addMatrix(1.0, pjj, hsj);
	//	}
	//}
	//for (int i = 0; i < 6; ++i)
	//	PC(i, i) = 1.0;
	//PC.addMatrix(1.0, PT, -1.0);

	// method 2: compute PT = sum(H(sj)*pjj) --- PC = sum(H(-sj)*pjj)
	// then PO = I - PT - PC
	// finally PT += PO/2 --- PC += PO/2
	for (int j = 0; j < 3; j++) {
		computePjj(V, j, pjj);
		double hsj = Heavyside(Si(j));
		if(hsj > 0.0)
			PT.addMatrix(1.0, pjj, hsj);
		hsj = Heavyside(-Si(j));
		if (hsj > 0.0)
			PC.addMatrix(1.0, pjj, hsj);
	}

	static Matrix PO(6, 6); // PO = I - PT - PC
	PO.addMatrix(0.0, PT, -1.0);
	PO.addMatrix(1.0, PC, -1.0);
	for (int i = 0; i < 6; ++i)
		PO(i, i) += 1.0;
	PT.addMatrix(1.0, PO, 0.5);
	PC.addMatrix(1.0, PO, 0.5);

	// compute positive and negative parts of the stress
	ST.addMatrixVector(0.0, PT, S, 1.0);
	SC.addMatrixVector(0.0, PC, S, 1.0);

	// R factor
	R = 0.0;
	double Rnum = 0.0;
	double Rden = 0.0;
	for (int j = 0; j < 3; ++j) {
		double Sj = Si(j);
		if (Sj > 0.0) Sj *= cdf;
		if (Sj > 0.0) Rnum += Sj;
		Rden += std::abs(Sj);
	}
	if (Rden > 0.0)
		R = Rnum / Rden;

	// done
	return 0;
}

void ASDConcrete3DMaterial::StressDecomposition::recompose(const Vector& S, Vector& Sv) const
{
	Sv(0) = S(0) * std::pow(V(0, 0), 2) + S(1) * std::pow(V(0, 1), 2) + S(2) * std::pow(V(0, 2), 2);
	Sv(1) = S(0) * std::pow(V(1, 0), 2) + S(1) * std::pow(V(1, 1), 2) + S(2) * std::pow(V(1, 2), 2);
	Sv(2) = S(0) * std::pow(V(2, 0), 2) + S(1) * std::pow(V(2, 1), 2) + S(2) * std::pow(V(2, 2), 2);
	Sv(3) = S(0) * V(0, 0) * V(1, 0) + S(1) * V(0, 1) * V(1, 1) + S(2) * V(0, 2) * V(1, 2);
	Sv(4) = S(0) * V(1, 0) * V(2, 0) + S(1) * V(1, 1) * V(2, 1) + S(2) * V(1, 2) * V(2, 2);
	Sv(5) = S(0) * V(0, 0) * V(2, 0) + S(1) * V(0, 1) * V(2, 1) + S(2) * V(0, 2) * V(2, 2);
}

void ASDConcrete3DMaterial::StressDecomposition::recompose(Vector& Sv) const
{
	recompose(Si, Sv);
}

ASDConcrete3DMaterial::CrackPlanesStorage& ASDConcrete3DMaterial::CrackPlanesStorage::instance()
{
	static CrackPlanesStorage _instance;
	return _instance;
}

ASDConcrete3DMaterial::CrackPlanesStorage::Vector3ListPointer ASDConcrete3DMaterial::CrackPlanesStorage::get(int n90)
{
	// return nullptr for 0
	if (n90 < 1)
		return nullptr;
	// get it or make it
	auto& item = m_map[n90];
	if (item == nullptr) {
		// make it
		item = std::make_shared<Vector3List>();
		Vector3List& normals = *item;
		double dangle = M_PI / 2.0 / static_cast<double>(n90); // angle increment
		std::size_t N = static_cast<std::size_t>(n90) * 2; // number of subdivision per linear range [0 - 180[
		std::size_t M = N * (N - 1) + 1; // number of total unique samples
		normals.resize(M);
		std::size_t counter = 0;
		for (std::size_t j = 0; j < N; ++j) {
			double ay = static_cast<double>(j)* dangle;
			std::size_t Ni = j == 0 ? 1 : N;
			for (std::size_t i = 0; i < Ni; ++i) {
				double ax = static_cast<double>(i)* dangle;
				Vector3& nij = normals[counter++];
				nij.x = std::cos(ax) * std::sin(ay);
				nij.y = std::sin(ax) * std::sin(ay);
				nij.z = std::cos(ay);
			}
		}
	}
	// done
	return item;
}

ASDConcrete3DMaterial::CrackPlanes::CrackPlanes(int n90)
	: m_n90(std::max(0, n90))
	, m_normals(CrackPlanesStorage::instance().get(m_n90))
{
	std::size_t num = m_normals ? m_normals->size() : 1;
	m_equivalent_strain.resize(num, 0.0);
}

void ASDConcrete3DMaterial::CrackPlanes::setCurrentNormal(const Vector3& ni)
{
	// always reset to 0
	m_closest_normal_loc = 0;
	// search only if necessary
	if (m_normals) {
		m_current_normal = ni;
		double max_dot = -1.0;
		const auto& normals = *m_normals;
		for (std::size_t i = 0; i < normals.size(); ++i) {
			double adot = std::abs(m_current_normal.dot(normals[i]));
			if (adot > max_dot) {
				max_dot = adot;
				m_closest_normal_loc = i;
			}
		}
	}
}

double ASDConcrete3DMaterial::CrackPlanes::getCurrentEquivalentStrain() const
{
	if (m_closest_normal_loc < m_equivalent_strain.size())
		return m_equivalent_strain[m_closest_normal_loc];
	return 0.0;
}

void ASDConcrete3DMaterial::CrackPlanes::updateCurrentEquivalentStrain(double x, double smooth_angle)
{
	if (m_closest_normal_loc < m_equivalent_strain.size()) {
		// smooth only if necessary
		if (m_normals) {
			double sig = std::max(1.0e-6, smooth_angle) / 2.3546;
			double sden = 2.0 * sig * sig;
			const auto& normals = *m_normals;
			for (std::size_t i = 0; i < normals.size(); ++i) {
				double adot = std::abs(m_current_normal.dot(normals[i]));
				double angle = std::acos(adot);
				double X1 = std::exp(-std::pow(angle, 2) / sden);
				double X2 = std::exp(-std::pow(angle + M_PI, 2) / sden);
				double X3 = std::exp(-std::pow(angle - M_PI, 2) / sden);
				double XX = std::max(X1, std::max(X2, X3));
				double EEQ = m_equivalent_strain[i];
				EEQ = std::max(EEQ, x * XX);
				m_equivalent_strain[i] = EEQ;
			}
		}
		// impose the closest always
		double EEQ_closest = m_equivalent_strain[m_closest_normal_loc];
		EEQ_closest = std::max(EEQ_closest, x);
		m_equivalent_strain[m_closest_normal_loc] = EEQ_closest;
	}
}

void ASDConcrete3DMaterial::CrackPlanes::reset()
{
	for (std::size_t i = 0; i < m_equivalent_strain.size(); ++i)
		m_equivalent_strain[i] = 0.0;
}

double ASDConcrete3DMaterial::CrackPlanes::getEquivalentStrainAtNormal(std::size_t i) const
{
	if (i < m_equivalent_strain.size())
		return m_equivalent_strain[i];
	return 0.0;
}

void ASDConcrete3DMaterial::CrackPlanes::setEquivalentStrainAtNormal(std::size_t i, double x)
{
	if (i < m_equivalent_strain.size())
		m_equivalent_strain[i] = x;
}

const ASDConcrete3DMaterial::Vector3& ASDConcrete3DMaterial::CrackPlanes::getNormal(std::size_t i) const
{
	static Vector3 dummy;
	if (m_normals && i < m_normals->size())
		return m_normals->operator[](i);
	return dummy;
}

std::size_t ASDConcrete3DMaterial::CrackPlanes::getClosestNormal(const Vector3& N) const
{
	double dot_max = 0.0;
	std::size_t loc = 0;
	if (m_normals) {
		const auto& normals = *m_normals;
		for (std::size_t i = 0; i < normals.size(); ++i) {
			const Vector3& Ntrial = normals[i];
			double aid = std::abs(Ntrial.dot(N));
			if (aid > dot_max) {
				dot_max = aid;
				loc = i;
			}
		}
	}
	return loc;
}

std::vector<int> ASDConcrete3DMaterial::CrackPlanes::getMax3Normals(double smooth_angle) const
{
	std::vector<int> out;
	if (m_normals) {
		const auto& normals = *m_normals;
		double tol = smooth_angle;
		// find 1
		std::size_t p1 = 0;
		double v1 = 0.0;
		for (std::size_t i = 0; i < normals.size(); ++i) {
			double vtrial = m_equivalent_strain[i];
			if (vtrial > v1) {
				v1 = vtrial;
				p1 = i;
			}
		}
		if (v1 > 0.0) {
			out.push_back(static_cast<int>(p1));
			const Vector3& N1 = normals[p1];

			// find 2
			std::size_t p2 = 0;
			double v2 = 0.0;
			for (std::size_t i = 0; i < normals.size(); ++i) {
				double vtrial = m_equivalent_strain[i];
				const Vector3& Ntrial = normals[i];
				double A1 = std::acos(std::abs(Ntrial.dot(N1)));
				if (A1 < tol) continue;
				if (vtrial > v2) {
					v2 = vtrial;
					p2 = i;
				}
			}
			if (v2 > 0.0) {
				out.push_back(static_cast<int>(p2));
				const Vector3& N2 = normals[p2];

				// find 3
				std::size_t p3 = 0;
				double v3 = 0.0;
				for (std::size_t i = 0; i < normals.size(); ++i) {
					double vtrial = m_equivalent_strain[i];
					const Vector3& Ntrial = normals[i];
					double A1 = std::acos(std::abs(Ntrial.dot(N1)));
					if (A1 < tol) continue;
					double A2 = std::acos(std::abs(Ntrial.dot(N2)));
					if (A2 < tol) continue;
					if (vtrial > v3) {
						v3 = vtrial;
						p3 = i;
					}
				}
				if (v3 > 0.0) {
					out.push_back(static_cast<int>(p3));
				}
			}
		}
	}
	return out;
}

int ASDConcrete3DMaterial::CrackPlanes::serializationDataSize() const
{
	return 6 + static_cast<int>(m_equivalent_strain.size());
}

void ASDConcrete3DMaterial::CrackPlanes::serialize(Vector& data, int& pos)
{
	data(pos++) = static_cast<double>(m_n90);
	data(pos++) = static_cast<double>(m_equivalent_strain.size());
	data(pos++) = m_current_normal.x;
	data(pos++) = m_current_normal.y;
	data(pos++) = m_current_normal.z;
	data(pos++) = static_cast<double>(m_closest_normal_loc);
	for (double i : m_equivalent_strain)
		data(pos++) = i;
}

void ASDConcrete3DMaterial::CrackPlanes::deserialize(Vector& data, int& pos)
{
	m_n90 = static_cast<int>(data(pos++));
	m_normals = CrackPlanesStorage::instance().get(m_n90);
	m_equivalent_strain.resize(static_cast<std::size_t>(data(pos++)));
	m_current_normal.x = data(pos++);
	m_current_normal.y = data(pos++);
	m_current_normal.z = data(pos++);
	m_closest_normal_loc = static_cast<std::size_t>(data(pos++));
	for (std::size_t i = 0; i < m_equivalent_strain.size(); ++i)
		m_equivalent_strain[i] = data(pos++);
}

ASDConcrete3DMaterial::ASDConcrete3DMaterial(
	int _tag,
	double _E,
	double _v,
	double _rho,
	double _eta,
	double _Kc,
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
	const HardeningLaw& _hc,
	double _cdf,
	int _nct,
	int _ncc,
	double _smoothing_angle)
	: NDMaterial(_tag, ND_TAG_ASDConcrete3DMaterial)
	, E(_E)
	, v(_v)
	, rho(_rho)
	, eta(_eta)
	, Kc(_Kc)
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
	, cdf(_cdf)
	, nct(std::max(0, _nct))
	, ncc(std::max(0, _ncc))
	, smoothing_angle(std::abs(_smoothing_angle) * M_PI / 180.0)
	, svt(nct)
	, svt_commit(nct)
	, svt_commit_old(nct)
	, svc(ncc)
	, svc_commit(ncc)
	, svc_commit_old(ncc)
{
	// intialize C as C0
	C = getInitialTangent();

	// initialize PT_commit as eye(6)*0.5
	for (int i = 0; i < 6; ++i)
		PT_commit(i, i) = 0.5;

	// compute fc/ft ratio.
	// Note that the Lubliner surface, being developed for concrete-like materials, assumes fc > ft!
	// do not allow fc/ft < 5
	fcft_ratio = 1.0; // start with the default
	double fcmax = hc.computeMaxStress();
	double ftmax = ht.computeMaxStress();
	if (ftmax > 0.0) {
		fcft_ratio = std::max(5.0, fcmax / ftmax);
	}
}

ASDConcrete3DMaterial::ASDConcrete3DMaterial()
	: NDMaterial(0, ND_TAG_ASDConcrete3DMaterial)
{
}

ASDConcrete3DMaterial::~ASDConcrete3DMaterial()
{ 
}

double ASDConcrete3DMaterial::getRho(void)
{
	return rho;
}

int ASDConcrete3DMaterial::setTrialStrain(const Vector& v)
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
		static Matrix Cnum(6, 6);
		// strain perturbation parameter
		double PERT = (ht.strainTolerance() + hc.strainTolerance()) / 2.0;
		// compute the forward perturbed solution and store in Cnum columns
		for (int j = 0; j < 6; ++j) {
			strain = v;
			strain(j) += PERT;
			retval = compute(true, false);
			if (retval < 0) return retval;
			for (int i = 0; i < 6; ++i)
				Cnum(i, j) = stress(i);
		}
		// compute unperturbed solution
		strain = v;
		retval = compute(true, false);
		if (retval < 0) return retval;
		for (int j = 0; j < 6; ++j) {
			for (int i = 0; i < 6; ++i)
				Cnum(i, j) = (Cnum(i, j) - stress(i)) / PERT;
		}
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
				static Matrix aux = Matrix(6, 6);
				static Vector stress_implicit(6);
				aux = PT_commit;
				double R_aux = R_commit;
				retval = compute(false, false);
				if (retval < 0) return retval;
				stress_implicit = stress;
				double dt_implicit = dt_bar;
				double dc_implicit = dc_bar;
				// the implicit pass has overwritten the frozen split: put it
				// back before extrapolating, or the extrapolation is not the
				// one this step would have taken. Nothing else needs restoring
				// here, because compute() always restarts from the committed
				// state
				PT_commit = aux;
				R_commit = R_aux;
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

int ASDConcrete3DMaterial::setTrialStrain(const Vector& v, const Vector& /*r*/)
{
	return setTrialStrain(v);
}

int ASDConcrete3DMaterial::setTrialStrainIncr(const Vector& v)
{
	static Vector aux(6);
	aux = strain;
	aux.addVector(1.0, v, 1.0);
	return setTrialStrain(aux);
}

int ASDConcrete3DMaterial::setTrialStrainIncr(const Vector& v, const Vector& /*r*/)
{
	return setTrialStrainIncr(v);
}

const Vector &ASDConcrete3DMaterial::getStrain(void)
{
	return strain;
}

const Vector &ASDConcrete3DMaterial::getStress(void)
{
	return stress;
}

const Matrix &ASDConcrete3DMaterial::getTangent(void)
{
	return C;
}

const Matrix &ASDConcrete3DMaterial::getInitialTangent(void)
{
	static Matrix D(6, 6);
	D.Zero();
	double mu2 = E / (1.0 + v);
	double lam = v * mu2 / (1.0 - 2.0 * v);
	double mu = 0.50 * mu2;
	mu2 += lam;
	D(0, 0) = D(1, 1) = D(2, 2) = mu2;
	D(0, 1) = D(1, 0) = D(0, 2) = D(2, 0) = D(1, 2) = D(2, 1) = lam;
	D(3, 3) = D(4, 4) = D(5, 5) = mu;
	return D;
}

int ASDConcrete3DMaterial::commitState(void)
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
	// store the previously committed variables for next move from n to n - 1
	svt_commit_old = svt_commit;
	svc_commit_old = svc_commit;
	// store committed variables
	svt_commit = svt;
	svc_commit = svc;
	strain_commit = strain;
	stress_eff_commit = stress_eff;
	dtime_n_commit = dtime_n;
	// for output purposes
	xt_max_commit = xt_max;
	xc_max_commit = xc_max;
	// done
	commit_done = true;
	return 0;
}

int ASDConcrete3DMaterial::revertToLastCommit(void)
{
	// restore converged values
	svt = svt_commit;
	svc = svc_commit;
	strain = strain_commit;
	stress_eff = stress_eff_commit;
	dtime_n = dtime_n_commit;
	// for output purposes
	xt_max = xt_max_commit;
	xc_max = xc_max_commit;
	// done
	return 0;
}

int ASDConcrete3DMaterial::revertToStart(void)
{
	// State variables
	svt.reset();
	svt_commit.reset();
	svt_commit_old.reset();
	svc.reset();
	svc_commit.reset();
	svc_commit_old.reset();

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
	strain.Zero();
	strain_commit.Zero();
	stress.Zero();
	stress_implex.Zero();
	stress_eff.Zero();
	stress_eff_commit.Zero();
	C = getInitialTangent();
	PT_commit.Zero();
	for (int i = 0; i < 6; ++i)
		PT_commit(i, i) = 0.5;
	R_commit = 0.0;

	// Output variables
	dt_bar = 0.0;
	dc_bar = 0.0;
	xt_max = 0.0;
	xt_max_commit = 0.0;
	xc_max = 0.0;
	xc_max_commit = 0.0;
	iso_crack_normal.Zero();
	iso_crush_normal.Zero();

	// Done
	return 0;
}

NDMaterial * ASDConcrete3DMaterial::getCopy(void)
{
	// we can safely use the default copy-constructor according to the member variables we're using
	return new ASDConcrete3DMaterial(*this);
}

NDMaterial* ASDConcrete3DMaterial::getCopy(const char* code)
{
	if (strcmp(code, "ThreeDimensional") == 0)
		return getCopy();
	return NDMaterial::getCopy(code);
}

const char* ASDConcrete3DMaterial::getType(void) const
{
	return "ThreeDimensional";
}

int ASDConcrete3DMaterial::getOrder(void) const
{
	return 6;
}

void ASDConcrete3DMaterial::Print(OPS_Stream &s, int flag)
{
	s << "ASDConcrete3D Material, tag: " << this->getTag() << "\n";
}

int ASDConcrete3DMaterial::sendSelf(int commitTag, Channel &theChannel)
{
	// aux
	int counter;

	// variable DBL data size (138: the six components of stress_implex were
	// added to the state)
	int nv_dbl = 138 +
		ht.serializationDataSize() +
		hc.serializationDataSize() +
		svt.serializationDataSize() +
		svt_commit.serializationDataSize() +
		svt_commit_old.serializationDataSize() +
		svc.serializationDataSize() +
		svc_commit.serializationDataSize() +
		svc_commit_old.serializationDataSize();

	// send INT data
	static ID idata(12);
	counter = 0;
	idata(counter++) = getTag();
	idata(counter++) = static_cast<int>(implex);
	idata(counter++) = static_cast<int>(implex_control);
	idata(counter++) = static_cast<int>(implex_abort_on_error);
	idata(counter++) = static_cast<int>(tangent);
	idata(counter++) = static_cast<int>(auto_regularize);
	idata(counter++) = static_cast<int>(regularization_done);
	idata(counter++) = nct;
	idata(counter++) = ncc;
	idata(counter++) = static_cast<int>(dtime_is_user_defined);
	idata(counter++) = static_cast<int>(commit_done);
	idata(counter++) = nv_dbl;
	if (theChannel.sendID(getDbTag(), commitTag, idata) < 0) {
		opserr << "ASDConcrete3DMaterial::sendSelf() - failed to send INT data\n";
		return -1;
	}

	// send DBL data
	Vector ddata(nv_dbl);
	counter = 0;
	ddata(counter++) = E;
	ddata(counter++) = v;
	ddata(counter++) = rho;
	ddata(counter++) = eta;
	ddata(counter++) = Kc;
	ddata(counter++) = implex_error_tolerance;
	ddata(counter++) = implex_time_redution_limit;
	ddata(counter++) = implex_alpha;
	ddata(counter++) = lch;
	ddata(counter++) = lch_ref;
	ddata(counter++) = smoothing_angle;
	ddata(counter++) = dtime_n;
	ddata(counter++) = dtime_n_commit;
	ddata(counter++) = dtime_0;
	ddata(counter++) = implex_error;
	ddata(counter++) = implex_gap;
	for (int i = 0; i < 6; ++i)
		for (int j = 0; j < 6; ++j)
			ddata(counter++) = PT_commit(i, j);
	ddata(counter++) = R_commit;
	ddata(counter++) = cdf;
	for (int i = 0; i < 6; ++i) ddata(counter++) = strain(i);
	for (int i = 0; i < 6; ++i) ddata(counter++) = strain_commit(i);
	for (int i = 0; i < 6; ++i) ddata(counter++) = stress(i);
	for (int i = 0; i < 6; ++i) ddata(counter++) = stress_implex(i);
	for (int i = 0; i < 6; ++i) ddata(counter++) = stress_eff(i);
	for (int i = 0; i < 6; ++i) ddata(counter++) = stress_eff_commit(i);
	for (int i = 0; i < 6; ++i)
		for (int j = 0; j < 6; ++j)
			ddata(counter++) = C(i, j);
	ddata(counter++) = dt_bar;
	ddata(counter++) = dc_bar;
	for (int i = 0; i < 3; ++i) ddata(counter++) = iso_crack_normal(i);
	for (int i = 0; i < 3; ++i) ddata(counter++) = iso_crush_normal(i);
	ddata(counter++) = xt_max;
	ddata(counter++) = xt_max_commit;
	ddata(counter++) = xc_max;
	ddata(counter++) = xc_max_commit;
	ht.serialize(ddata, counter);
	hc.serialize(ddata, counter);
	svt.serialize(ddata, counter);
	svt_commit.serialize(ddata, counter);
	svt_commit_old.serialize(ddata, counter);
	svc.serialize(ddata, counter);
	svc_commit.serialize(ddata, counter);
	svc_commit_old.serialize(ddata, counter);
	if (theChannel.sendVector(getDbTag(), commitTag, ddata) < 0) {
		opserr << "ASDConcrete3DMaterial::sendSelf() - failed to send DBL data\n";
		return -1;
	}

	// done
	return 0;
}

int ASDConcrete3DMaterial::recvSelf(int commitTag, Channel & theChannel, FEM_ObjectBroker & theBroker)
{
	// aux
	int counter;

	// recv INT data
	static ID idata(12);
	if (theChannel.recvID(getDbTag(), commitTag, idata) < 0) {
		opserr << "ASDConcrete3DMaterial::recvSelf() - failed to receive INT data\n";
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
	nct = idata(counter++);
	ncc = idata(counter++);
	dtime_is_user_defined = static_cast<bool>(idata(counter++));
	commit_done = static_cast<bool>(idata(counter++));
	int nv_dbl = idata(counter++);
	
	// recv DBL data
	Vector ddata(nv_dbl);
	if (theChannel.recvVector(getDbTag(), commitTag, ddata) < 0) {
		opserr << "ASDConcrete3DMaterial::recvSelf() - failed to receive DBL data\n";
		return -1;
	}
	counter = 0;
	E = ddata(counter++);
	v = ddata(counter++);
	rho = ddata(counter++);
	eta = ddata(counter++);
	Kc = ddata(counter++);
	implex_error_tolerance = ddata(counter++);
	implex_time_redution_limit = ddata(counter++);
	implex_alpha = ddata(counter++);
	lch = ddata(counter++);
	lch_ref = ddata(counter++);
	smoothing_angle = ddata(counter++);
	dtime_n = ddata(counter++);
	dtime_n_commit = ddata(counter++);
	dtime_0 = ddata(counter++);
	implex_error = ddata(counter++);
	implex_gap = ddata(counter++);
	for (int i = 0; i < 6; ++i)
		for (int j = 0; j < 6; ++j)
			PT_commit(i, j) = ddata(counter++);
	R_commit = ddata(counter++);
	cdf = ddata(counter++);
	for (int i = 0; i < 6; ++i) strain(i) = ddata(counter++);
	for (int i = 0; i < 6; ++i) strain_commit(i) = ddata(counter++);
	for (int i = 0; i < 6; ++i) stress(i) = ddata(counter++);
	for (int i = 0; i < 6; ++i) stress_implex(i) = ddata(counter++);
	for (int i = 0; i < 6; ++i) stress_eff(i) = ddata(counter++);
	for (int i = 0; i < 6; ++i) stress_eff_commit(i) = ddata(counter++);
	for (int i = 0; i < 6; ++i)
		for (int j = 0; j < 6; ++j)
			C(i, j) = ddata(counter++);
	dt_bar = ddata(counter++);
	dc_bar = ddata(counter++);
	for (int i = 0; i < 3; ++i) iso_crack_normal(i) = ddata(counter++);
	for (int i = 0; i < 3; ++i) iso_crush_normal(i) = ddata(counter++);
	xt_max = ddata(counter++);
	xt_max_commit = ddata(counter++);
	xc_max = ddata(counter++);
	xc_max_commit = ddata(counter++);
	ht.deserialize(ddata, counter);
	hc.deserialize(ddata, counter);
	svt.deserialize(ddata, counter);
	svt_commit.deserialize(ddata, counter);
	svt_commit_old.deserialize(ddata, counter);
	svc.deserialize(ddata, counter);
	svc_commit.deserialize(ddata, counter);
	svc_commit_old.deserialize(ddata, counter);

	// done
	return 0;
}

int ASDConcrete3DMaterial::setParameter(const char** argv, int argc, Parameter& param)
{
	// 1000 - elasticity & mass
	if (strcmp(argv[0], "E") == 0) {
		param.setValue(E);
		return param.addObject(1000, this);
	}
	if (strcmp(argv[0], "v") == 0) {
		param.setValue(v);
		return param.addObject(1001, this);
	}
	if (strcmp(argv[0], "rho") == 0) {
		param.setValue(rho);
		return param.addObject(1002, this);
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

int ASDConcrete3DMaterial::updateParameter(int parameterID, Information& info)
{
	switch (parameterID) {
		// 1000 - elasticity & mass
	case 1000:
		E = info.theDouble;
		return 0;
	case 1001:
		v = info.theDouble;
		return 0;
	case 1002:
		rho = info.theDouble;
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

Response* ASDConcrete3DMaterial::setResponse(const char** argv, int argc, OPS_Stream& output)
{
	// utils
	auto make_resp = [&output, this](int rid, const Vector& v, const std::vector<std::string>* labels = nullptr) -> MaterialResponse* {
		output.tag("NdMaterialOutput");
		output.attr("matType", getClassType());
		output.attr("matTag", getTag());
		if (labels) {
			for(const auto& item : (*labels))
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
	// same six labels as the rest of the family, see ASDPlasticDamageConcrete3D
	static std::vector<std::string> lb_tensor = { "11", "22", "33", "12", "23", "13" };
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
		if (strcmp(argv[0], "damage") == 0 || strcmp(argv[0], "Damage") == 0) {
			if(argc > 1 && strcmp(argv[1], "-avg") == 0)
				return make_resp(2000, getAvgDamage(), &lb_damage);
			else
				return make_resp(2001, getMaxDamage(), &lb_damage);
		}
		if (strcmp(argv[0], "equivalentPlasticStrain") == 0 || strcmp(argv[0], "EquivalentPlasticStrain") == 0) {
			if (argc > 1 && strcmp(argv[1], "-avg") == 0)
				return make_resp(2002, getAvgEquivalentPlasticStrain(), &lb_eqpl_strain);
			else
				return make_resp(2003, getMaxEquivalentPlasticStrain(), &lb_eqpl_strain);
		}
		if (strcmp(argv[0], "equivalentTotalStrain") == 0 || strcmp(argv[0], "EquivalentTotalStrain") == 0) {
			if (argc > 1 && strcmp(argv[1], "-avg") == 0)
				return make_resp(2004, getAvgStrainMeasure(), &lb_tot_strain);
			else
				return make_resp(2005, getMaxStrainMeasure(), &lb_tot_strain);
		}
		if (strcmp(argv[0], "cw") == 0 || strcmp(argv[0], "crackWidth") == 0 || strcmp(argv[0], "CrackWidth") == 0) {
			if (argc > 1 && strcmp(argv[1], "-avg") == 0)
				return make_resp(2006, getAvgCrackWidth(), &lb_cw);
			else
				return make_resp(2007, getMaxCrackWidth(), &lb_cw);
		}
		if (strcmp(argv[0], "crackPattern") == 0 || strcmp(argv[0], "CrackPattern") == 0) {
			return make_resp(2008, getCrackPattern(), &lb_crackpattern);
		}
		if (strcmp(argv[0], "crushPattern") == 0 || strcmp(argv[0], "CrushPattern") == 0) {
			return make_resp(2009, getCrushPattern(), &lb_crackpattern);
		}
		if (strcmp(argv[0], "crackInfo") == 0 || strcmp(argv[0], "CrackInfo") == 0) {
			if (argc > 3) {
				Vector3 N;
				if (string_to_double(argv[1], N.x) && 
					string_to_double(argv[2], N.y) && 
					string_to_double(argv[3], N.z)) {
					std::size_t Npos = svt.getClosestNormal(N);
					Cinfo(0) = static_cast<double>(Npos);
					Cinfo(1) = svt.getEquivalentStrainAtNormal(Npos);
					return make_resp(2010, Cinfo);
				}
			}
		}
		if (strcmp(argv[0], "crushInfo") == 0 || strcmp(argv[0], "CrushInfo") == 0) {
			if (argc > 3) {
				Vector3 N;
				if (string_to_double(argv[1], N.x) &&
					string_to_double(argv[2], N.y) &&
					string_to_double(argv[3], N.z)) {
					std::size_t Npos = svc.getClosestNormal(N);
					Cinfo(0) = static_cast<double>(Npos);
					Cinfo(1) = svc.getEquivalentStrainAtNormal(Npos);
					return make_resp(2011, Cinfo);
				}
			}
		}
		if (strcmp(argv[0], "crackStrain") == 0 || strcmp(argv[0], "CrackStrain") == 0) {
			double user_lch_ref = this->lch_ref; // by default the input one
			if (argc > 2 && strcmp(argv[1], "-lchRef") == 0) {
				double trial_lch_ref = 0.0;
				if (string_to_double(argv[2], trial_lch_ref)) {
					user_lch_ref = trial_lch_ref;
				}
			}
			Cinfo(1) = user_lch_ref;
			Cinfo(0) = getMaxCrackWidth()(0) / user_lch_ref;
			return make_resp(2012, Cinfo, &lb_crack_strain);
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
			Cinfo(0) = getMaxCrushWidth()(0) / user_lch_ref;
			return make_resp(2013, Cinfo, &lb_crush_strain);
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
			return make_resp(3003, getImplexStress(), &lb_tensor);
		}
		// 4000 - internal time
		if (strcmp(argv[0], "time") == 0 || strcmp(argv[0], "Time") == 0) {
			return make_resp(4000, getTimeIncrements(), &lb_time);
		}
	}

	// otherwise return base-class response
	return NDMaterial::setResponse(argv, argc, output);
}

int ASDConcrete3DMaterial::getResponse(int responseID, Information& matInformation)
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
	case 2000: return matInformation.setVector(getAvgDamage());
	case 2001: return matInformation.setVector(getMaxDamage());
	case 2002: return matInformation.setVector(getAvgEquivalentPlasticStrain());
	case 2003: return matInformation.setVector(getMaxEquivalentPlasticStrain());
	case 2004: return matInformation.setVector(getAvgStrainMeasure());
	case 2005: return matInformation.setVector(getMaxStrainMeasure());
	case 2006: return matInformation.setVector(getAvgCrackWidth());
	case 2007: return matInformation.setVector(getMaxCrackWidth());
	case 2008: return matInformation.setVector(getCrackPattern());
	case 2009: return matInformation.setVector(getCrushPattern());
	case 2010:
		if (matInformation.theVector && matInformation.theVector->Size() == 2) {
			std::size_t Npos = static_cast<std::size_t>(matInformation.theVector->operator()(0));
			matInformation.theVector->operator()(1) = svt.getEquivalentStrainAtNormal(Npos);
			return 0;
		}
		break;
	case 2011:
		if (matInformation.theVector && matInformation.theVector->Size() == 2) {
			std::size_t Npos = static_cast<std::size_t>(matInformation.theVector->operator()(0));
			matInformation.theVector->operator()(1) = svc.getEquivalentStrainAtNormal(Npos);
			return 0;
		}
		break;
	case 2012:
		if (matInformation.theVector && matInformation.theVector->Size() == 2) {
			double user_lch_ref = matInformation.theVector->operator()(1);
			matInformation.theVector->operator()(0) = getMaxCrackWidth()(0) / user_lch_ref;
			return 0;
		}
		break;
	case 2013:
		if (matInformation.theVector && matInformation.theVector->Size() == 2) {
			double user_lch_ref = matInformation.theVector->operator()(1);
			matInformation.theVector->operator()(0) = getMaxCrushWidth()(0) / user_lch_ref;
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
	return NDMaterial::getResponse(responseID, matInformation);
}

int ASDConcrete3DMaterial::compute(bool do_implex, bool do_tangent)
{
	// get committed variables
	svt = svt_commit;
	svc = svc_commit;
	stress_eff = stress_eff_commit;
	xt_max = xt_max_commit;
	xc_max = xc_max_commit;

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
	static Vector dStrain(6);
	dStrain = strain;
	dStrain.addVector(1.0, strain_commit, -1.0);
	stress_eff.addMatrixVector(1.0, getInitialTangent(), dStrain, 1.0);

	// compute stress split
	static StressDecomposition D;
	if (implex && do_implex) {
		// explicit: PT and R from committed
		D.R = R_commit;
		D.PT = PT_commit;
		D.PC.Zero();
		for (int i = 0; i < 6; ++i)
			D.PC(i, i) = 1.0;
		D.PC.addMatrix(1.0, D.PT, -1.0);
		// update the ST and SC with the committed PT and PC
		D.ST.addMatrixVector(0.0, D.PT, stress_eff, 1.0);
		D.SC.addMatrixVector(0.0, D.PC, stress_eff, 1.0);
	}
	else {
		// implicit
		D.compute(stress_eff, cdf);
		// update normals
		Vector3 Tnormal(D.V(0, 0), D.V(1, 0), D.V(2, 0));
		Vector3 Cnormal(D.V(0, 2), D.V(1, 2), D.V(2, 2));
		svt.setCurrentNormal(Tnormal);
		svc.setCurrentNormal(Cnormal);
	}

	// compute committed hardening variables
	HardeningLawPoint pt = ht.evaluateAt(svt.getCurrentEquivalentStrain());
	HardeningLawPoint pc = hc.evaluateAt(svc.getCurrentEquivalentStrain());

	// temporary clone of old equivalent plastic strains
	double xt_pl = pt.plasticStrain(E);
	double xc_pl = pc.plasticStrain(E);

	// compute new trial equivalent strain measures
	if (implex && do_implex) {
		// extrapolated equivalent strain measures (explicit)
		for (std::size_t i = 0; i < svt.count(); ++i) {
			double xn = svt_commit.getEquivalentStrainAtNormal(i);
			double xnn = svt_commit_old.getEquivalentStrainAtNormal(i);
			double x_new = xn + time_factor * (xn - xnn);
			svt.setEquivalentStrainAtNormal(i, x_new);
		}
		for (std::size_t i = 0; i < svc.count(); ++i) {
			double xn = svc_commit.getEquivalentStrainAtNormal(i);
			double xnn = svc_commit_old.getEquivalentStrainAtNormal(i);
			double x_new = xn + time_factor * (xn - xnn);
			svc.setEquivalentStrainAtNormal(i, x_new);
		}
		pt = ht.evaluateAt(svt.getCurrentEquivalentStrain());
		pc = hc.evaluateAt(svc.getCurrentEquivalentStrain());
	}
	else {
		// compute trial strain measures (implicit)
		double xt_trial = equivalentTensileStrainMeasure(D.Si(0), D.Si(1), D.Si(2)) + xt_pl;
		double xc_trial = equivalentCompressiveStrainMeasure(D.Si(0), D.Si(1), D.Si(2)) + xc_pl;
		// update hardening variables
		if (xt_trial > pt.x) {
			double xt_new = rate_coeff_1 * pt.x + rate_coeff_2 * xt_trial;
			pt = ht.evaluateAt(xt_new);
			svt.updateCurrentEquivalentStrain(pt.x, smoothing_angle);
			if (xt_trial > xt_max) {
				xt_max = xt_trial;
				iso_crack_normal(0) = D.V(0, 0);
				iso_crack_normal(1) = D.V(1, 0);
				iso_crack_normal(2) = D.V(2, 0);
			}
		}
		if (xc_trial > pc.x) {
			double xc_new = rate_coeff_1 * pc.x + rate_coeff_2 * xc_trial;
			pc = hc.evaluateAt(xc_new);
			svc.updateCurrentEquivalentStrain(pc.x, smoothing_angle);
			if (xc_trial > xc_max) {
				xc_max = xc_trial;
				iso_crush_normal(0) = D.V(0, 2);
				iso_crush_normal(1) = D.V(1, 2);
				iso_crush_normal(2) = D.V(2, 2);
			}
		}
	}

	// compute plastic damage
	double seff_eq_t = (pt.x - xt_pl) * E;
	double dt_plastic = seff_eq_t > 0.0 ? 1.0 - pt.q / seff_eq_t : 0.0;
	double seff_eq_c = (pc.x - xc_pl) * E;
	double dc_plastic = seff_eq_c > 0.0 ? 1.0 - pc.q / seff_eq_c : 0.0;

	// mix damages
	auto mix_dam = [](double dt, double dc) -> double {
		double alpha = 1.0;
		double W = (1.0 - dc) * (1.0 - alpha * D.R * dt);
		double D = 1.0 - W;
		return D;
	};
	{
		auto auxp = hc.evaluateAt(pt.x * fcft_ratio);
		double aux_pl = auxp.plasticStrain(E);
		double aux_seff_eq = (auxp.x - aux_pl) * E;
		double aux_d_plastic = aux_seff_eq > 0.0 ? 1.0 - auxp.q / aux_seff_eq : 0.0;
		dc_plastic = mix_dam(aux_d_plastic, dc_plastic);
		pc.d = mix_dam(auxp.d, pc.d);
	}

	// update effective stress
	stress_eff.addVector(0.0, D.ST, 1.0 - dt_plastic);
	stress_eff.addVector(1.0, D.SC, 1.0 - dc_plastic);

	// update nominal stress
	dt_bar = pt.d + dt_plastic - pt.d * dt_plastic;
	dc_bar = pc.d + dc_plastic - pc.d * dc_plastic;
	stress.addVector(0.0, D.ST, 1.0 - dt_bar);
	stress.addVector(1.0, D.SC, 1.0 - dc_bar);

	// tangent matrix
	if (do_tangent) {
		static Matrix W(6, 6);
		W.Zero();
		for (int i = 0; i < 6; ++i)
			W(i, i) = 1.0;
		W.addMatrix(1.0, D.PT, -dt_bar);
		W.addMatrix(1.0, D.PC, -dc_bar);
		C.addMatrixProduct(0.0, W, getInitialTangent(), 1.0);
	}

	// save real PT and R, if mp.implex and !do_implex -> called from commit
	if (implex && !do_implex) {
		// save it in implex mode during implicit phase
		PT_commit = D.PT;
		R_commit = D.R;
	}

	// done
	return 0;
}

double ASDConcrete3DMaterial::lublinerCriterion(double s1, double s2, double s3, double ft, double fc, double k1, double scale) const
{
	double fb = 1.16 * fc;
	double gamma = 3.0 * (1.0 - Kc) / (2.0 * Kc - 1.0);
	double alpha = (fb - fc) / (2.0 * fb - fc);
	double I1 = s1 + s2 + s3;
	double p = I1 / 3.0;
	Vector3 Sdev(s1 - p, s2 - p, s3 - p);
	double J2 = 0.5 * Sdev.dot(Sdev);
	double beta = fc / ft * (1.0 - alpha) - (1.0 + alpha);
	double smax_alg = s1;
	double smax = Macauley(smax_alg);
	double smin = Macauley(-smax_alg);
	return (1.0 / (1.0 - alpha) * (alpha * I1 + std::sqrt(3.0 * J2) + k1 * beta * smax - gamma * smin)) * scale;
}

double ASDConcrete3DMaterial::equivalentTensileStrainMeasure(double s1, double s2, double s3) const
{
	// skip if maximum principal stress is not strictly positive
	if (s1 < ht.stressTolerance())
		return 0.0;

	// compute equivalent strain measure
	// - Lubliner criterion
	// - normalized to fc
	// - assume ratio ft/fc = 0.1 (don't take the ratio from the hardening laws!)
	return lublinerCriterion(s1, s2, s3, 1.0/fcft_ratio, 1.0, 1.0, 1.0/fcft_ratio) / E;
}

double ASDConcrete3DMaterial::equivalentCompressiveStrainMeasure(double s1, double s2, double s3) const
{
	// consider only the negative part of the stress
	s1 = std::min(s1, 0.0);
	s2 = std::min(s2, 0.0);
	s3 = std::min(s3, 0.0);

	// compute equivalent strain measure
	// - Lubliner criterion
	// - normalized to fc
	return lublinerCriterion(s1, s2, s3, 1.0/fcft_ratio, 1.0, 1.0, 1.0) / E;
}

Vector ASDConcrete3DMaterial::getHardeningLawVector(HardeningLawType ltype, HardeningLawPointComponent c) const
{
	Vector r;
	const HardeningLaw& h = ltype == HardeningLawType::Tension ? ht : hc;
	r.resize(static_cast<int>(h.points().size()));
	for (std::size_t i = 0; i < h.points().size(); ++i) {
		const HardeningLawPoint& p = h.points()[i];
		switch (c)
		{
		case ASDConcrete3DMaterial::HardeningLawPointComponent::TotalStrain:
			r(static_cast<int>(i)) = p.totalStrain();
			break;
		case ASDConcrete3DMaterial::HardeningLawPointComponent::EffectiveStress:
			r(static_cast<int>(i)) = p.effectiveStress();
			break;
		case ASDConcrete3DMaterial::HardeningLawPointComponent::NominalStress:
			r(static_cast<int>(i)) = p.stress();
			break;
		default:
			break;
		}
	}
	return r;
}

const Vector& ASDConcrete3DMaterial::getMaxStrainMeasure() const
{
	static Vector d(2);
	double xt_max = 0.0;
	double xc_max = 0.0;
	for (std::size_t i = 0; i < svt.count(); ++i)
		xt_max = std::max(xt_max, svt.getEquivalentStrainAtNormal(i));
	for (std::size_t i = 0; i < svc.count(); ++i)
		xc_max = std::max(xc_max, svc.getEquivalentStrainAtNormal(i));
	d(0) = xt_max;
	d(1) = xc_max;
	return d;
}

const Vector& ASDConcrete3DMaterial::getAvgStrainMeasure() const
{
	static Vector d(2);
	double xt = 0.0;
	double xc = 0.0;
	if (svt.count() > 0) {
		for (std::size_t i = 0; i < svt.count(); ++i)
			xt += svt.getEquivalentStrainAtNormal(i);
		xt /= static_cast<double>(svt.count());
	}
	if (svc.count() > 0) {
		for (std::size_t i = 0; i < svc.count(); ++i)
			xc += svc.getEquivalentStrainAtNormal(i);
		xc /= static_cast<double>(svc.count());
	}
	d(0) = xt;
	d(1) = xc;
	return d;
}

const Vector& ASDConcrete3DMaterial::getMaxDamage() const
{
	static Vector d(2);
	const Vector& x = getMaxStrainMeasure();
	d(0) = ht.evaluateAt(x(0)).crackingDamage();
	d(1) = hc.evaluateAt(x(1)).crackingDamage();
	return d;
}

const Vector& ASDConcrete3DMaterial::getAvgDamage() const
{
	static Vector d(2);
	const Vector& x = getAvgStrainMeasure();
	d(0) = ht.evaluateAt(x(0)).crackingDamage();
	d(1) = hc.evaluateAt(x(1)).crackingDamage();
	return d;
}

const Vector& ASDConcrete3DMaterial::getMaxEquivalentPlasticStrain() const
{
	static Vector d(2);
	const Vector& x = getMaxStrainMeasure();
	d(0) = ht.evaluateAt(x(0)).plasticStrain(E);
	d(1) = hc.evaluateAt(x(1)).plasticStrain(E);
	return d;
}

const Vector& ASDConcrete3DMaterial::getAvgEquivalentPlasticStrain() const
{
	static Vector d(2);
	const Vector& x = getAvgStrainMeasure();
	d(0) = ht.evaluateAt(x(0)).plasticStrain(E);
	d(1) = hc.evaluateAt(x(1)).plasticStrain(E);
	return d;
}

const Vector& ASDConcrete3DMaterial::getMaxCrackWidth() const
{
	static Vector d(1);
	d.Zero();
	if (ht.hasStrainSoftening()) {
		double e0 = ht.strainAtOnsetOfCrack();
		const Vector& x = getMaxStrainMeasure();
		d(0) = std::max(x(0) - e0, 0.0) * lch;
	}
	return d;
}

const Vector& ASDConcrete3DMaterial::getAvgCrackWidth() const
{
	static Vector d(1);
	d.Zero();
	if (ht.hasStrainSoftening()) {
		double e0 = ht.strainAtOnsetOfCrack();
		const Vector& x = getAvgStrainMeasure();
		d(0) = std::max(x(0) - e0, 0.0) * lch;
	}
	return d;
}

const Vector& ASDConcrete3DMaterial::getMaxCrushWidth() const
{
	static Vector d(1);
	d.Zero();
	if (hc.hasStrainSoftening()) {
		double e0 = hc.strainAtOnsetOfCrack();
		const Vector& x = getMaxStrainMeasure();
		d(0) = std::max(x(1) - e0, 0.0) * lch;
	}
	return d;
}

const Vector& ASDConcrete3DMaterial::getAvgCrushWidth() const
{
	static Vector d(1);
	d.Zero();
	if (hc.hasStrainSoftening()) {
		double e0 = hc.strainAtOnsetOfCrack();
		const Vector& x = getAvgStrainMeasure();
		d(0) = std::max(x(1) - e0, 0.0) * lch;
	}
	return d;
}

const Vector& ASDConcrete3DMaterial::getCrackPattern() const
{
	static Vector d(9);
	d.Zero();
	if (ht.hasStrainSoftening()) {
		double e0 = ht.strainAtOnsetOfCrack();
		if (svt.count() > 1) {
			std::vector<int> normals = svt.getMax3Normals(smoothing_angle);
			int pos = 0;
			for (int nid : normals) {
				double crstrain = std::max(svt.getEquivalentStrainAtNormal(static_cast<std::size_t>(nid)) - e0, 0.0);
				double crdisp = crstrain * lch;
				const Vector3& N = svt.getNormal(static_cast<std::size_t>(nid));
				d(pos + 0) = N.x * crdisp;
				d(pos + 1) = N.y * crdisp;
				d(pos + 2) = N.z * crdisp;
				pos += 3;
			}
		}
		else {
			double crstrain = std::max(xt_max - e0, 0.0);
			double crdisp = crstrain * lch;
			d(0) = iso_crack_normal(0) * crdisp;
			d(1) = iso_crack_normal(1) * crdisp;
			d(2) = iso_crack_normal(2) * crdisp;
		}
	}
	return d;
}

const Vector& ASDConcrete3DMaterial::getCrushPattern() const
{
	static Vector d(9);
	d.Zero();
	if (hc.hasStrainSoftening()) {
		double e0 = hc.strainAtOnsetOfCrack();
		if (svc.count() > 1) {
			std::vector<int> normals = svc.getMax3Normals(smoothing_angle);
			int pos = 0;
			for (int nid : normals) {
				double crstrain = std::max(svc.getEquivalentStrainAtNormal(static_cast<std::size_t>(nid)) - e0, 0.0);
				double crdisp = crstrain * lch;
				const Vector3& N = svc.getNormal(static_cast<std::size_t>(nid));
				d(pos + 0) = N.x * crdisp;
				d(pos + 1) = N.y * crdisp;
				d(pos + 2) = N.z * crdisp;
				pos += 3;
			}
		}
		else {
			double crstrain = std::max(xc_max - e0, 0.0);
			double crdisp = crstrain * lch;
			d(0) = iso_crush_normal(0) * crdisp;
			d(1) = iso_crush_normal(1) * crdisp;
			d(2) = iso_crush_normal(2) * crdisp;
		}
	}
	return d;
}

const Vector& ASDConcrete3DMaterial::getImplexError() const
{
	static Vector d(1);
	d(0) = implex_error;
	return d;
}

const Vector& ASDConcrete3DMaterial::getImplexStress() const
{
	return stress_implex;
}

const Vector& ASDConcrete3DMaterial::getImplexGap() const
{
	static Vector d(1);
	d(0) = implex_gap;
	return d;
}

double ASDConcrete3DMaterial::computeImplexErrorMetric(void)
{
	// no extrapolation, no error. This is what makes it safe for the
	// aggregation to call it on everything it has
	if (!implex)
		return 0.0;
	// the current state holds the EXPLICIT answer: the stress this step
	// delivered to the element, and the state the recorders must keep seeing
	static TrialState delivered;
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
	// rejected, and revertToLastCommit() would not put PT_commit and R_commit
	// back
	restoreTrialState(delivered);
	implex_error = err;
	return err;
}

double ASDConcrete3DMaterial::implexTimeRatio(void) const
{
	return dtime_0 > 0.0 ? dtime_n / dtime_0 : 1.0;
}

double ASDConcrete3DMaterial::implexStressGap(const Vector& delivered, const Vector& stress_implicit) const
{
	// the NORM is the largest Voigt component - the 1D twin takes an absolute
	// value - but the DENOMINATOR is the same in both, the largest stress the
	// material can carry, or the same tolerance would mean different things in
	// 1D and in 3D and a stepper would be comparing numbers that are not
	// comparable
	double gap = 0.0;
	for (int i = 0; i < 6; ++i)
		gap = std::max(gap, std::abs(delivered(i) - stress_implicit(i)));
	double ref = std::max(ht.computeMaxStress(), hc.computeMaxStress());
	if (ref <= 0.0)
		ref = 1.0;
	return gap / ref;
}

double ASDConcrete3DMaterial::implexDamageGap(double dt_implex, double dc_implex,
	double dt_implicit, double dc_implicit) const
{
	// LEGACY, for output only: see the note in the header on why the model
	// does not steer on this
	return std::max(std::abs(dt_implex - dt_implicit), std::abs(dc_implex - dc_implicit));
}

void ASDConcrete3DMaterial::saveTrialState(TrialState& x) const
{
	x.svt = svt;
	x.svc = svc;
	x.stress = stress;
	x.stress_eff = stress_eff;
	x.C = C;
	x.dt_bar = dt_bar;
	x.dc_bar = dc_bar;
	x.iso_crack_normal = iso_crack_normal;
	x.iso_crush_normal = iso_crush_normal;
	x.xt_max = xt_max;
	x.xc_max = xc_max;
	x.PT_commit = PT_commit;
	x.R_commit = R_commit;
}

void ASDConcrete3DMaterial::restoreTrialState(const TrialState& x)
{
	svt = x.svt;
	svc = x.svc;
	stress = x.stress;
	stress_eff = x.stress_eff;
	C = x.C;
	dt_bar = x.dt_bar;
	dc_bar = x.dc_bar;
	iso_crack_normal = x.iso_crack_normal;
	iso_crush_normal = x.iso_crush_normal;
	xt_max = x.xt_max;
	xc_max = x.xc_max;
	PT_commit = x.PT_commit;
	R_commit = x.R_commit;
}

const Vector& ASDConcrete3DMaterial::getTimeIncrements() const
{
	static Vector d(3);
	d(0) = dtime_n;
	d(1) = dtime_n_commit;
	d(2) = dtime_0;
	return d;
}

