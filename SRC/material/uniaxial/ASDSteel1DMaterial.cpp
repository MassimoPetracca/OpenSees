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
	// ...
	int tag;
	
	// parse mandatory arguments
	
	// this is how you get 1 integer
	//if (OPS_GetInt(&numData, &tag) != 0) {
	//	print_err("invalid 'tag'");
	//	return nullptr;
	//}
	
	// this is how you get 1 integer
	//if (OPS_GetDouble(&numData, &tag) != 0) {
	//	print_err("invalid 'tag'");
	//	return nullptr;
	//}
	
	// parse optional parameters
	while (OPS_GetNumRemainingInputArgs() > 0) {
		const char* value = OPS_GetString();
		
		// this is how you get an optional string argument
		//if (strcmp(value, "-option") == 0) {
		//	the_option = true;
		//}
		
		//else if (strcmp(value, "-option") == 0) {
		//	if (OPS_GetNumRemainingInputArgs() > 0) {
		//		if (OPS_GetDouble(&numData, &the_option_value) < 0) {
		//			print_err("Failed to get optional -option argument");
		//			return nullptr;
		//		}
		//	}
		//	else {
		//		print_err("-option argument requested but not provided");
		//		return nullptr;
		//	}
		//}
		
	}

	// do some checks of the input parameters
	// ...
	
	// create the material
	UniaxialMaterial* instance = new ASDSteel1DMaterial(
		// place here the input arguments
		// ...
	);
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
	bool _good_bond,
	double _mm,
	double _N,
	UniaxialMaterial* _steel_mat,
	UniaxialMaterial* _bond_mat)
	: UniaxialMaterial(_tag, MAT_TAG_ASDSteel1DMaterial)
	, Es(_Es)
	, fy(_fy)
	, Hs(_Hs)
	, fc(_fc)
	, diam(_diam)
	, Leff(_Leff)
	, good_bond(_good_bond)
	, mm(_mm)
	, N(_N)
	, steel_mat(_steel_mat)
	, bond_mat(_bond_mat)
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
	if(steel_mat) {
		delete steel_mat;
		steel_mat = nullptr;
	}
	if(bond_mat) {
		delete bond_mat;
		bond_mat = nullptr;
	}
}

int ASDSteel1DMaterial::setTrialStrain(double v, double r)
{
	// return value
	int retval = 0;
	
	// do your computation here and return 0 if no error occourred.
	
	// done
	return retval;
}

double ASDSteel1DMaterial::getStress(void)
{
	// return the stress
	return 0.0;
}

double ASDSteel1DMaterial::getTangent(void)
{
	// return the current tangent
	return 0.0;
}

double ASDSteel1DMaterial::getInitialTangent(void)
{
	// return the initial tangent
	return 0.0;
}

double ASDSteel1DMaterial::getStrain(void)
{
	// return the current strain
	return strain;
}

int ASDSteel1DMaterial::commitState(void)
{
	// save current state
	return 0;
}

int ASDSteel1DMaterial::revertToLastCommit(void)
{
	// restore converged (previous) state
	return 0;
}

int ASDSteel1DMaterial::revertToStart(void)
{
	// restart
	return 0;
}

UniaxialMaterial* ASDSteel1DMaterial::getCopy(void)
{
	// use the default copy-constructor to copy all POD variables
	ASDSteel1DMaterial* copy = new ASDSteel1DMaterial(*this);
	
	// the pointers to the 2 uniaxial materials are not POD! copy them
	// ...
	
	return copy;
}

void ASDSteel1DMaterial::Print(OPS_Stream& s, int flag)
{
	// print info
	s << "ASDSteel1D Material, tag: " << this->getTag() << "\n";
}

int ASDSteel1DMaterial::sendSelf(int commitTag, Channel &theChannel)
{
	// we can pack:
	// * this tag (1)
	// * all variables (12) +
	// * 2 tags for each sub-material (4)
	// in a vector
	// with a total size = 1+12+4=17
	static Vector D(17);
	int counter = 0;
	
	// put the tag first
	D(counter++) = static_cast<double>(getTag());
	
	// then all variables
	D(counter++) = Es;
	// continue here
	// ...
	
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
	// receive the 17-component DBL vector
	static Vector D(17);
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
	// continue here
	// ...
	
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
	// ...

	// default
	return -1;
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
	static Vector d;
	// resize d and put proper values in it
	return d;
}

const Vector& ASDSteel1DMaterial::getBondResponse() const
{
	static Vector d;
	// resize d and put proper values in it
	return d;
}

const Vector& ASDSteel1DMaterial::getEquivalentBondResponse() const
{
	static Vector d;
	// resize d and put proper values in it
	return d;
}

const Vector& ASDSteel1DMaterial::getSteelResponse() const
{
	static Vector d;
	// resize d and put proper values in it
	return d;
}

const Vector& ASDSteel1DMaterial::getImplexError() const
{
	static Vector d;
	// resize d and put proper values in it
	return d;
}

