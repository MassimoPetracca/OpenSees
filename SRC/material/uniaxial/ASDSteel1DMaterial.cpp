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
// $Source: /usr/local/cvs/OpenSees/SRC/material/uniaxial/ASDSteel1DMaterial.cpp,v $

// Massimo Petracca - ASDEA Software, Italy
//
// A Simple model for reinforcing steel
//

#include <ASDSteel1DMaterial.h>
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
#include <FEM_ObjectBroker.h>
#include <ASDConcrete1DMaterial.h>
#include <ElasticMaterial.h>
#include <ParallelMaterial.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

void* OPS_ASDSteel1DMaterial()
{
	// some kudos
	static bool first_done = false;
	if (!first_done) {
		opserr << "Using ASDSteel1D - Developed by: Massimo Petracca, Guido Camata, ASDEA Software Technology\n";
		first_done = true;
	}

	// the command description
	static const char* cmd_descr =
		"uniaxialMaterial ASDSteel1D $tag $Es $fy $Hs $fc $diam $Leff "
		"<-good_bond> "
		"<-mm $mm> <-N $N> "
		"<-tangent> "
		"<-implex> "
		"<-eta $eta> "
		"<-autoRegularization>";

	// a utility to print error message
	auto print_err = [](const char* msg) {
		opserr << "uniaxialMaterial ASDSteel1D Error:\n" << msg << "\n" << cmd_descr << "\n";
	};

	// check arguments
	int numArgs = OPS_GetNumRemainingInputArgs();
	if (numArgs < 7) {
		print_err("Few arguments (< 7)");
		return nullptr;
	}

	// numData
	int numData = 1;

	// allocate input arguments and initialize optional ones
	int tag;
	double Es = 0.0;
	double fy = 0.0;
	double Hs = 0.0;
	double fc = 0.0;
	double diam = 0.0;
	double Leff = 0.0;
	bool good_bond = false;
	double mm = 1.0;
	double N = 1.0;
	UniaxialMaterial* steel_mat = nullptr;
	UniaxialMaterial* bond_mat = nullptr;
	// optional arguments to be used to construct sub-materials
	bool tangent = false;
	bool implex = false;
	double eta = 0.0;
	bool auto_reg = false;

	// parse mandatory arguments
	if (OPS_GetInt(&numData, &tag) != 0) {
		print_err("invalid 'tag'");
		return nullptr;
	}
	else if (OPS_GetDouble(&numData, &Es) != 0) {
		print_err("invalid 'Es'");
		return nullptr;
	}
	else if (OPS_GetDouble(&numData, &fy) != 0) {
		print_err("invalid 'fy'");
		return nullptr;
	}
	else if (OPS_GetDouble(&numData, &Hs) != 0) {
		print_err("invalid 'Hs'");
		return nullptr;
	}
	else if (OPS_GetDouble(&numData, &fc) != 0) {
		print_err("invalid 'fc'");
		return nullptr;
	}
	else if (OPS_GetDouble(&numData, &diam) != 0) {
		print_err("invalid 'diam'");
		return nullptr;
	}
	else if (OPS_GetDouble(&numData, &Leff) != 0) {
		print_err("invalid 'Leff'");
		return nullptr;
	}

	// parse optional parameters
	while (OPS_GetNumRemainingInputArgs() > 0) {
		const char* value = OPS_GetString();
		if (strcmp(value, "-good_bond") == 0) {
			good_bond = true;
		}
		else if (strcmp(value, "-mm") == 0) {
			if (OPS_GetNumRemainingInputArgs() > 0) {
				if (OPS_GetDouble(&numData, &mm) < 0) {
					print_err("Failed to get optional -mm argument");
					return nullptr;
				}
			}
			else {
				print_err("-mm argument requested but not provided");
				return nullptr;
			}
		}
		else if (strcmp(value, "-N") == 0) {
			if (OPS_GetNumRemainingInputArgs() > 0) {
				if (OPS_GetDouble(&numData, &N) < 0) {
					print_err("Failed to get optional -N argument");
					return nullptr;
				}
			}
			else {
				print_err("-N argument requested but not provided");
				return nullptr;
			}
		}
		else if (strcmp(value, "-tangent") == 0) {
			tangent = true;
		}
		else if (strcmp(value, "-implex") == 0) {
			implex = true;
		}
		else if (strcmp(value, "-eta") == 0) {
			if (OPS_GetNumRemainingInputArgs() > 0) {
				if (OPS_GetDouble(&numData, &eta) < 0) {
					print_err("Failed to get optional -eta argument");
					return nullptr;
				}
			}
			else {
				print_err("-eta argument requested but not provided");
				return nullptr;
			}
		}
		else if (strcmp(value, "-autoRegularization") == 0) {
			auto_reg = true;
		}
	}

	// do some checks of the input parameters
	if (Es <= 0.0) {
		print_err("Es must be strictly positive");
		return nullptr;
	}
	if (fy <= 0.0) {
		print_err("fy must be strictly positive");
		return nullptr;
	}
	if (Hs < 0.0) {
		print_err("Hs must be positive or zero");
		return nullptr;
	}
	if (fc <= 0.0) {
		print_err("fc must be strictly positive");
		return nullptr;
	}
	if (diam <= 0.0) {
		print_err("diam must be strictly positive");
		return nullptr;
	}
	if (Leff <= 0.0) {
		print_err("Leff must be strictly positive");
		return nullptr;
	}

	// create sub materials
	int base_tag = -1000001;

	// create the perfectly plastic material
	{
		ASDConcrete1DMaterial::HardeningLaw ht_pp(base_tag--, ASDConcrete1DMaterial::HardeningLawType::Tension, Es,
			{ 0.0, fy / Es, 2.0 * fy / Es }, { 0.0, fy, fy }, { 0.0, 0.0, 0.0 });
		ASDConcrete1DMaterial::HardeningLaw hc_pp(base_tag--, ASDConcrete1DMaterial::HardeningLawType::Compression, Es,
			{ 0.0, fy / Es, 2.0 * fy / Es }, { 0.0, fy, fy }, { 0.0, 0.0, 0.0 });
		UniaxialMaterial* steel_pp = new ASDConcrete1DMaterial(
			base_tag--, Es, eta,
			implex, false, 0.0, 0.0, 1.0,
			tangent, auto_reg, Leff, ht_pp, hc_pp);
		// create the linear material for kinematic hardening
		UniaxialMaterial* steel_hk = new ElasticMaterial(base_tag--, Es);
		// combine them into a parallel material
		UniaxialMaterial** parallel_mats = new UniaxialMaterial * [2];
		Vector* parallel_factors = new Vector(2);
		parallel_mats[0] = steel_pp;
		parallel_mats[1] = steel_hk;
		(*parallel_factors)[0] = 1.0 - Hs;
		(*parallel_factors)[1] = Hs;
		steel_mat = new ParallelMaterial(base_tag--, 2, parallel_mats, parallel_factors);
		// clean memory used to construct parallel material
		delete steel_pp;
		delete steel_hk;
		delete[] parallel_mats;
		delete parallel_factors;
	}

	// create the bond material
	// using the Bond-Slip law as per Model Code 2010 Volume 1 @ 6.1.1
	{
		// some parameters are required in N-mm
		double MPa = N / mm / mm;
		double _fc = fc / MPa;
		double _diam = diam / mm;
		// some settings for Bond-Slip (tau-slip) response in N-mm
		double alpha = 0.4;
		double s3 = std::max(25.4, _diam);
		double beta = 1.25;
		double s1 = 1.8;
		double s2 = 3.6;
		if (good_bond) {
			beta = 2.5;
			s1 = 1.0;
			s2 = 2.0;
		}
		double tau_max = beta * std::sqrt(_fc);
		double tau_f = 0.4 * tau_max;
		// define discrete slip values in mm (20 points in the first parabolic parts)
		std::vector<double> slip(22);
		std::vector<double> bond(22);
		std::vector<double> damg(22);
		for (std::size_t i = 1; i < 20; ++i)
			slip[i] = static_cast<double>(i) / 19.0 * s1;
		slip[20] = s2;
		slip[21] = s3;
		// fill bond (tau) stresses in MPa
		for (std::size_t i = 0; i < slip.size(); ++i) {
			double s = slip[i];
			double tau = tau_f;
			if (s <= s1)
				tau = tau_max * std::pow(s / s1, alpha);
			else if (s <= s2)
				tau = tau_max;
			else if (s <= s3)
				tau = tau_max - (tau_max - tau_f) * (s - s2) / (s3 - s2);
			bond[i] = tau;
		}
		// define the slip scale (from mm to user length, and from disp to strain)
		double slip_scale = mm / Leff;
		// define the bond scale (from MPa to user pressure, and from tau to equivalent stress)
		double radius = diam / 2.0;
		double Asteel = M_PI * radius * radius;
		double Asurf = 2.0 * M_PI * radius * Leff;
		double bond_scale = MPa * Asurf / Asteel;
		// convert
		for (std::size_t i = 0; i < slip.size(); ++i) {
			slip[i] *= slip_scale;
			bond[i] *= bond_scale;
			damg[i] = 0.0;
		}
		// compute the initial stiffness of the equivalent bond material
		double bond_stiffness = bond[1] / slip[1];
		// create dummy vectors for the zero-stress zone
		double zero_tau = tau_f * 1.0e-4;
		double zero_tau_strain = zero_tau / bond_stiffness;
		std::vector<double> Ne = { 0.0, zero_tau_strain, zero_tau_strain * 2.0 };
		std::vector<double> Ns = { 0.0, zero_tau, zero_tau };
		std::vector<double> Nd = { 1.0, 1.0, 1.0 };
		// the bond-slip material in positive direction
		UniaxialMaterial* bond_pos = new ASDConcrete1DMaterial(
			base_tag--, bond_stiffness, eta,
			implex, false, 0.0, 0.0, 1.0,
			tangent, auto_reg, Leff,
			ASDConcrete1DMaterial::HardeningLaw(base_tag--, ASDConcrete1DMaterial::HardeningLawType::Tension, bond_stiffness, slip, bond, damg),
			ASDConcrete1DMaterial::HardeningLaw(base_tag--, ASDConcrete1DMaterial::HardeningLawType::Compression, bond_stiffness, Ne, Ns, Nd));
		// the bond-slip material in negative direction
		UniaxialMaterial* bond_neg = new ASDConcrete1DMaterial(
			base_tag--, bond_stiffness, eta,
			implex, false, 0.0, 0.0, 1.0,
			tangent, auto_reg, Leff,
			ASDConcrete1DMaterial::HardeningLaw(base_tag--, ASDConcrete1DMaterial::HardeningLawType::Tension, bond_stiffness, Ne, Ns, Nd),
			ASDConcrete1DMaterial::HardeningLaw(base_tag--, ASDConcrete1DMaterial::HardeningLawType::Compression, bond_stiffness, slip, bond, damg));
		// the bond-slip material in both direction for the frictional part of the cycle
		UniaxialMaterial* bond_res = new ASDConcrete1DMaterial(
			base_tag--, bond_stiffness, eta,
			implex, false, 0.0, 0.0, 1.0,
			tangent, auto_reg, Leff,
			ASDConcrete1DMaterial::HardeningLaw(base_tag--, ASDConcrete1DMaterial::HardeningLawType::Tension, bond_stiffness, slip, bond, damg),
			ASDConcrete1DMaterial::HardeningLaw(base_tag--, ASDConcrete1DMaterial::HardeningLawType::Compression, bond_stiffness, slip, bond, damg));
		// combine them into a parallel material
		UniaxialMaterial** parallel_mats = new UniaxialMaterial * [3];
		Vector* parallel_factors = new Vector(3);
		parallel_mats[0] = bond_pos;
		parallel_mats[1] = bond_neg;
		parallel_mats[2] = bond_res;
		(*parallel_factors)[0] = 0.5;
		(*parallel_factors)[1] = 0.5;
		(*parallel_factors)[2] = 0.5;
		bond_mat = new ParallelMaterial(base_tag--, 3, parallel_mats, parallel_factors);
		// clean memory used to construct parallel material
		delete bond_pos;
		delete bond_neg;
		delete bond_res;
		delete[] parallel_mats;
		delete parallel_factors;
	}

	// create the material
	UniaxialMaterial* instance = new ASDSteel1DMaterial(
		tag, Es, fy, Hs, fc, diam, Leff, steel_mat, bond_mat
	);

	// clean memory used to construct the new material (steel_mat and bond_mat are copyied by this material!)
	delete steel_mat;
	delete bond_mat;

	// final checks
	if (instance == nullptr) {
		opserr << "UniaxialMaterial ASDSteel1D Error: failed to allocate a new material.\n";
		return nullptr;
	}
	return instance;
}

ASDSteel1DMaterial::ASDSteel1DMaterial(
	int _tag,
	double _Es,
	double _fy,
	double _Hs,
	double _fc,
	double _diam,
	double _Leff,
	UniaxialMaterial* _steel_mat,
	UniaxialMaterial* _bond_mat)
	: UniaxialMaterial(_tag, MAT_TAG_ASDSteel1DMaterial)
	, Es(_Es)
	, fy(_fy)
	, Hs(_Hs)
	, fc(_fc)
	, diam(_diam)
	, Leff(_Leff)
	, steel_mat(_steel_mat->getCopy())
	, bond_mat(_bond_mat->getCopy())
{
	// you can do some extra initialization here
}

ASDSteel1DMaterial::ASDSteel1DMaterial()
	: UniaxialMaterial(0, MAT_TAG_ASDSteel1DMaterial)
{
	// this is an empty constructor.
	// no need to do anything here, all variables are initialized in the class definition.
}

ASDSteel1DMaterial::~ASDSteel1DMaterial()
{
	// clean memory if needed.
	if (steel_mat) {
		delete steel_mat;
		steel_mat = nullptr;
	}
	if (bond_mat) {
		delete bond_mat;
		bond_mat = nullptr;
	}
}

int ASDSteel1DMaterial::setTrialStrain(double v, double r)
{
	// return value
	int retval = 0;

	// do your computation here and return 0 if no error occourred.

	// set current total strain
	strain = v;

	// obj function
	auto iso_stress_cond = [this](double S, double& F, double& dFdS) {
		double steel_strain = strain - S;
		steel_mat->setTrialStrain(steel_strain);
		bond_mat->setTrialStrain(S);
		double steel_stress = steel_mat->getStress();
		double bond_stress = bond_mat->getStress();
		F = bond_stress - steel_stress;
		double Es = steel_mat->getTangent();
		double Eb = bond_mat->getTangent();
		dFdS = Eb + Es;
	};

	// initial guess for the slip strain
	double s = bond_mat->getStrain();

	// iterative procedure to impose the iso-stress condition
	int max_iter = 100;
	double stress_tol = std::max(1.0e-6, 1.0e-4 * fy);
	double strain_tol = 1.0e-6;
	bool converged = false;
	double F, dFdS;
	for (int iter = 0; iter < max_iter; ++iter) {
		iso_stress_cond(s, F, dFdS);
		double ds = F / dFdS;
		s -= ds;
		/*opserr 
			<< "iter: " << iter
			<< " | F: " << std::abs(F) 
			<< " | ds: " << std::abs(ds) 
			<< " | dF = " << dFdS << "\n";*/
		if (std::abs(F) < stress_tol || std::abs(ds) < strain_tol) {
			converged = true;
			break;
		}
	}
	if (!converged)
		return -1;

	// done
	return retval;
}

double ASDSteel1DMaterial::getStress(void)
{
	// return the stress
	// note: the two materials are in series (same stress)
	// so we could pick any of them.
	// anyway we do an avg...
	return (steel_mat->getStress() + bond_mat->getStress()) / 2.0;
}

double ASDSteel1DMaterial::getTangent(void)
{
	// return the current tangent K = 1/(1/Ks + 1/Kb)

	// guard against div-by-zero
	constexpr double rel_tol = 1.0e-8;
	constexpr double abs_tol = 1.0e-12;
	double Ks_zero = std::max(abs_tol, rel_tol * steel_mat->getInitialTangent());
	double Kb_zero = std::max(abs_tol, rel_tol * bond_mat->getInitialTangent());

	double Ks = steel_mat->getTangent();
	double Kb = bond_mat->getTangent();
	if (std::abs(Ks) == 0.0)
		Ks = Ks_zero;
	if (std::abs(Kb) == 0.0)
		Kb = Kb_zero;

	return 1.0 / (1.0 / Ks + 1.0 / Kb);
}

double ASDSteel1DMaterial::getInitialTangent(void)
{
	// return the initial tangent  K = 1/(1/Ks + 1/Kb)
	
	// guard against div-by-zero
	constexpr double abs_tol = 1.0e-12;
	double Ks = steel_mat->getInitialTangent();
	double Kb = bond_mat->getInitialTangent();
	if (std::abs(Ks) == 0.0)
		Ks = abs_tol;
	if (std::abs(Kb) == 0.0)
		Kb = abs_tol;

	return 1.0 / (1.0 / Ks + 1.0 / Kb);
}

double ASDSteel1DMaterial::getStrain(void)
{
	// return the current strain
	return strain;
}

int ASDSteel1DMaterial::commitState(void)
{
	// save current state
	strain_commit = strain;
	int res = 0;
	res += steel_mat->commitState();
	res += bond_mat->commitState();
	return res;
}

int ASDSteel1DMaterial::revertToLastCommit(void)
{
	// restore converged (previous) state
	strain = strain_commit;
	int res = 0;
	res += steel_mat->revertToLastCommit();
	res += bond_mat->revertToLastCommit();
	return res;
}

int ASDSteel1DMaterial::revertToStart(void)
{
	// restart
	strain = 0.0;
	strain_commit = 0.0;
	energy = 0.0;
	int res = 0;
	res += steel_mat->revertToStart();
	res += bond_mat->revertToStart();
	return res;;
}

UniaxialMaterial* ASDSteel1DMaterial::getCopy(void)
{
	// use the default copy-constructor to copy all POD variables
	ASDSteel1DMaterial* copy = new ASDSteel1DMaterial(*this);

	// the pointers to the 2 uniaxial materials are not POD! copy them
	copy->steel_mat = steel_mat->getCopy();
	copy->bond_mat = bond_mat->getCopy();

	return copy;
}

void ASDSteel1DMaterial::Print(OPS_Stream& s, int flag)
{
	// print info
	s << "ASDSteel1D Material, tag: " << this->getTag() << "\n";
}

int ASDSteel1DMaterial::sendSelf(int commitTag, Channel& theChannel)
{
	// we can pack:
	// * this tag (1)
	// * all variables (9) +
	// * 2 tags for each sub-material (4)
	// in a vector
	// with a total size = 1+9+4=14
	static Vector D(14);
	int counter = 0;

	// put the tag first
	D(counter++) = static_cast<double>(getTag());

	// then all variables
	D(counter++) = Es;
	D(counter++) = fy;
	D(counter++) = Hs;
	D(counter++) = fc;
	D(counter++) = diam;
	D(counter++) = Leff;
	D(counter++) = strain;
	D(counter++) = strain_commit;
	D(counter++) = energy;

	// now put tags of sub-materials, we need them to restore those materials later on
	auto lambda_mat = [&theChannel, &counter](UniaxialMaterial* imat) {
		D(counter++) = static_cast<double>(imat->getClassTag());
		int mat_db_tag = imat->getDbTag();
		if (mat_db_tag == 0) {
			mat_db_tag = theChannel.getDbTag();
			if (mat_db_tag != 0)
				imat->setDbTag(mat_db_tag);
		}
		D(counter++) = static_cast<double>(mat_db_tag);
	};
	lambda_mat(steel_mat);
	lambda_mat(bond_mat);

	// now we can send the vector for this material
	if (theChannel.sendVector(getDbTag(), commitTag, D) < 0) {
		opserr << "ASDSteel1DMaterial::sendSelf() - failed to send DBL data\n";
		return -1;
	}

	// now we can send the 2 sub-materials
	if (steel_mat->sendSelf(commitTag, theChannel) < 0) {
		opserr << "ASDSteel1DMaterial::sendSelf() - failed to send steel material\n";
		return -1;
	}
	if (bond_mat->sendSelf(commitTag, theChannel) < 0) {
		opserr << "ASDSteel1DMaterial::sendSelf() - failed to send bond material\n";
		return -1;
	}

	// done
	return 0;
}

int ASDSteel1DMaterial::recvSelf(int commitTag, Channel& theChannel, FEM_ObjectBroker& theBroker)
{
	// receive the 14-component DBL vector
	static Vector D(14);
	if (theChannel.recvVector(getDbTag(), commitTag, D) < 0) {
		opserr << "ASDSteel1DMaterial::recvSelf() - failed to receive DBL data\n";
		return -1;
	}

	// parse
	int counter = 0;

	// get the tag first
	setTag(static_cast<int>(D(counter++)));

	// then all variables
	Es = D(counter++);
	fy = D(counter++);
	Hs = D(counter++);
	fc = D(counter++);
	diam = D(counter++);
	Leff = D(counter++);
	strain = D(counter++);
	strain_commit = D(counter++);
	energy = D(counter++);

	// receive sub-materials
	auto lambda_mat = [&theChannel, &theBroker, &commitTag, &counter]() -> UniaxialMaterial* {
		int mat_class_tag = D(counter++);
		int mat_db_tag = D(counter++);
		UniaxialMaterial* the_material = theBroker.getNewUniaxialMaterial(mat_class_tag);
		if (the_material == nullptr) {
			opserr << "ASDSteel1DMaterial::recvSelf() - could not get a new UniaxialMaterial\n";
			return nullptr;
		}
		the_material->setDbTag(mat_db_tag);
		if (the_material->recvSelf(commitTag, theChannel, theBroker) < 0) {
			opserr << "ASDSteel1DMaterial::recvSelf() - failed to receive material\n";
			return nullptr;
		}
		return the_material;
	};
	if (steel_mat) delete steel_mat;
	if (bond_mat) delete bond_mat;
	steel_mat = lambda_mat();
	if (steel_mat == nullptr) return -1;
	bond_mat = lambda_mat();
	if (bond_mat == nullptr) return -1;

	// done
	return 0;
}

int ASDSteel1DMaterial::setParameter(const char** argv, int argc, Parameter& param)
{
	// first try to handle parameters for this material
	// like this:

	// 1000 - elasticity
	// if (strcmp(argv[0], "E") == 0) {
		// param.setValue(E);
		// return param.addObject(1000, this);
	// }

	// if the input parameter is not handled by this material, forward it to the
	// sub-materials
	int res = -1;
	if (steel_mat->setParameter(argv, argc, param) == 0)
		res = 0;
	if (bond_mat->setParameter(argv, argc, param) == 0)
		res = 0;
	return res;
}

int ASDSteel1DMaterial::updateParameter(int parameterID, Information& info)
{
	// first try to handle parameters for this material
	// like this:
	switch (parameterID) {
		// 1000 - elasticity & mass
	//case 1000:
	//	E = info.theDouble;
	//	return 0;

		// default
	default:
		return -1;
	}
}

Response* ASDSteel1DMaterial::setResponse(const char** argv, int argc, OPS_Stream& output)
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
	static std::vector<std::string> lb_steel_failure = { "Fracture", "Buckling" };
	static std::vector<std::string> lb_bond_resp = { "Slip", "Tau" };
	static std::vector<std::string> lb_strain_stress = { "strain", "stress" };
	static std::vector<std::string> lb_implex_error = { "Error" };

	// check specific responses
	if (argc > 0) {
		// 1000 - sub-material responses
		if (strcmp(argv[0], "SteelFailure") == 0)
			return make_resp(1000, getSteelFailure(), &lb_steel_failure);
		if (strcmp(argv[0], "BondResponse") == 0)
			return make_resp(1001, getBondResponse(), &lb_bond_resp);
		if (strcmp(argv[0], "EquivalentBondResponse") == 0)
			return make_resp(1002, getEquivalentBondResponse(), &lb_strain_stress);
		if (strcmp(argv[0], "SteelResponse") == 0)
			return make_resp(1003, getSteelResponse(), &lb_strain_stress);
		// 3000 - implex error
		if (strcmp(argv[0], "implexError") == 0 || strcmp(argv[0], "ImplexError") == 0) {
			return make_resp(3000, getImplexError(), &lb_implex_error);
		}
	}

	// otherwise return base-class response
	return UniaxialMaterial::setResponse(argv, argc, output);
}

int ASDSteel1DMaterial::getResponse(int responseID, Information& matInformation)
{
	switch (responseID) {
		// 1000 - sub-material responses
	case 1000: return matInformation.setVector(getSteelFailure());
	case 1001: return matInformation.setVector(getBondResponse());
	case 1002: return matInformation.setVector(getEquivalentBondResponse());
	case 1003: return matInformation.setVector(getSteelResponse());
		// 3000 - implex error
	case 3000: return matInformation.setVector(getImplexError());
	default:
		break;
	}
	return UniaxialMaterial::getResponse(responseID, matInformation);
}

double ASDSteel1DMaterial::getEnergy(void)
{
	return energy;
}

const Vector& ASDSteel1DMaterial::getSteelFailure() const
{
	static Vector d(2);
	d.Zero(); // todo: not now...
	return d;
}

const Vector& ASDSteel1DMaterial::getBondResponse() const
{
	static Vector d(2);
	double radius = diam / 2.0;
	double Asteel = M_PI * radius * radius;
	double Asurf = 2.0 * M_PI * radius * Leff;
	d(0) = bond_mat->getStrain() * Leff;
	d(1) = bond_mat->getStress() * Asteel / Asurf;
	return d;
}

const Vector& ASDSteel1DMaterial::getEquivalentBondResponse() const
{
	static Vector d(2);
	d(0) = bond_mat->getStrain();
	d(1) = bond_mat->getStress();
	return d;
}

const Vector& ASDSteel1DMaterial::getSteelResponse() const
{
	static Vector d(2);
	d(0) = steel_mat->getStrain();
	d(1) = steel_mat->getStress();
	return d;
}

const Vector& ASDSteel1DMaterial::getImplexError() const
{
	static Vector d(1);
	d.Zero(); // todo: not now...
	return d;
}

