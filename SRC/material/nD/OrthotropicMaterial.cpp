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
// $Source: /usr/local/cvs/OpenSees/SRC/material/nD/OrthotropicMaterial.h,v $

// Davide Raino, Massimo Petracca - ASDEA Software, Italy
//
// A Generic Orthotropic Material Wrapper that can convert any
// nonlinear isotropic material into an orthotropic one by means of tensor
// mapping
//

#include <OrthotropicMaterial.h>
#include <Channel.h>
#include <FEM_ObjectBroker.h>
#include <OPS_Globals.h>
#include <elementAPI.h>

void *OPS_OrthotropicMaterial(void)
{
	NDMaterial *theOrthotropicMaterial = 0;

	int numArgs= OPS_GetNumRemainingInputArgs();
	if (numArgs < 17)
	{
		opserr << 
			"Want: nDMaterial Orthotropic $tag $theIsoMat $Ex $Ey $Ez $Gxy $Gyz $Gzx $vxy $vyz $vzx $Asigmaxx $Asigmayy $Asigmazz $Asigmaxyxy $Asigmayzyz $Asigmaxzxz.\n";
		return 0;
	}
	
	int iData[2];
	int numData = 2;

	if (OPS_GetInt(&numData, iData) != 0) 
	{
		opserr << "WARNING: invalid nDMaterial Orthotropic tags.\n";
		return 0;
	}

	double dData[15];
	numData = 15;

	if (OPS_GetDouble(&numData, dData) != 0) 
	{
		opserr << "WARNING: invalid data nDMaterial Orthotropic: " << iData[0] << ".\n";
		return 0;
	}

	NDMaterial *theIsoMaterial = OPS_getNDMaterial(iData[1]);
	if (theIsoMaterial == 0)
	{
		opserr << "WARNING: nDMaterial does not exist.\n";
		opserr << "nDMaterial: " << iData[1] << "\n";
		opserr << "nDMaterial Orthotropic: " << iData[0] << "\n";
		return 0;
	}

	theOrthotropicMaterial = new OrthotropicMaterial(iData[0], *theIsoMaterial,
		dData[0], dData[1], dData[2], dData[3],
		dData[4], dData[5], dData[6], dData[7],
		dData[8], dData[9], dData[10], dData[11],
		dData[12], dData[13], dData[14]);

	if (theOrthotropicMaterial == 0)
	{
		opserr << "WARNING: failed to create OrthotropicMaterial.\n";
		return 0;
	}

	return theOrthotropicMaterial;
}

OrthotropicMaterial::OrthotropicMaterial(int tag, NDMaterial &theIsoMat,
	double Ex, double Ey, double Ez, double Gxy, double Gyz, double Gzx,
	double vxy, double vyz, double vzx,
	double Asigmaxx, double Asigmayy, double Asigmazz, double Asigmaxyxy, double Asigmayzyz, double Asigmaxzxz)
	: NDMaterial(tag, ND_TAG_OrthotropicMaterial)
	, C(6,6)
	, C0(6,6)
	, epsilon(6)
	, sigma(6)
	, Asigma(6,6)
	, Aepsilon(6,6)
	, commit_C(6,6)
	, commit_epsilon(6)
	, commit_sigma(6)
{
	theIsotropicMaterial = theIsoMat.getCopy("ThreeDimensional");
	if (theIsotropicMaterial == 0)
	{
		opserr << "Error: OrthotropicMaterial - failed to get a copy of the isotropic material\n";
		exit(-1);
	}

	double vyx = vxy * Ey / Ex;
	double vzy = vyz * Ez / Ey;
	double vxz = vzx * Ex / Ez;

	double d = (1.0 - vxy * vyx - vyz * vzy - vzx * vxz - 2.0*vxy*vyz*vzx) / (Ex*Ey*Ez);

	C0(0, 0) = (1.0 - vyz * vzy) / (Ey*Ez*d);
	C0(1, 1) = (1.0 - vzx * vxz) / (Ez*Ex*d);
	C0(2, 2) = (1.0 - vxy * vyx) / (Ex*Ey*d);

	C0(1, 0) = (vxy + vxz * vzy) / (Ez*Ex*d);
	C0(0, 1) = C0(1, 0);

	C0(2, 0) = (vxz + vxy * vyz) / (Ex*Ey*d);
	C0(0, 2) = C0(2, 0);

	C0(2, 1) = (vyz + vxz * vyx) / (Ex*Ey*d);
	C0(1, 2) = C0(2, 1);

	C0(3, 3) = Gxy;
	C0(4, 4) = Gyz;
	C0(5, 5) = Gzx;

	if (Asigmaxx <= 0 || Asigmayy <= 0 || Asigmazz <= 0 || Asigmaxyxy <= 0 || Asigmayzyz <= 0 || Asigmaxzxz <= 0)
	{
		opserr << "Error: OrthotropicMaterial - Asigma11, Asigma22, Asigma33, Asigma12, Asigma23, Asigma13 must be greater than 0.\n";
		exit(-1);
	}

	Asigma(0, 0) = Asigmaxx;
	Asigma(1, 1) = Asigmayy;
	Asigma(2, 2) = Asigmazz;
	Asigma(3, 3) = Asigmaxyxy;
	Asigma(4, 4) = Asigmayzyz;
	Asigma(5, 5) = Asigmaxzxz;

	Matrix C0iso(6, 6);
	Matrix C0iso_inv(6, 6);

	C0iso = theIsotropicMaterial->getInitialTangent();
	if (C0iso.noRows() != 6 || C0iso.noCols() != 6) 
	{
		opserr << "Error: OrthotropicMaterial - the isotropic material must be a 3D material\n";
		exit(-1);
	}
	int res = C0iso.Invert(C0iso_inv);
	if (res < 0)
	{
		opserr << "Error: OrthotropicMaterial - the isotropic material gave a singular initial tangent\n";
		exit(-1);
	}

	Aepsilon = C0iso_inv*Asigma*C0;
}

OrthotropicMaterial::OrthotropicMaterial()
	: NDMaterial(0, ND_TAG_OrthotropicMaterial)
	, C(6, 6)
	, C0(6, 6)
	, epsilon(6)
	, sigma(6)
	, Asigma(6, 6)
	, Aepsilon(6, 6)
	, theIsotropicMaterial(0)
{ }

OrthotropicMaterial::OrthotropicMaterial(int tag_copy, NDMaterial &theIsoMat_copy,
	Matrix C0_copy, Matrix Aepsilon_copy, Matrix Asigma_copy)
	: NDMaterial(tag_copy, ND_TAG_OrthotropicMaterial)
	, C(6, 6)
	, C0(6, 6)
	, epsilon(6)
	, sigma(6)
	, Asigma(6, 6)
	, Aepsilon(6, 6)
	, commit_C(6, 6)
	, commit_epsilon(6)
	, commit_sigma(6)
{
	theIsotropicMaterial = theIsoMat_copy.getCopy();
	C0 = C0_copy;
	Aepsilon = Aepsilon_copy;
	Asigma = Asigma_copy;
}

OrthotropicMaterial::~OrthotropicMaterial()
{ 
	if (theIsotropicMaterial)
		delete theIsotropicMaterial;
}

int OrthotropicMaterial::setTrialStrain(const Vector & strain)
{
	epsilon = strain;

	// move to isotropic space
	static Vector eps_iso(6);
	eps_iso = Aepsilon * epsilon;

	// call isotropic material
	int iso_res = theIsotropicMaterial->setTrialStrain(eps_iso);
	if (iso_res != 0)
	{
		opserr << "Error: OrthotropicMaterial - the isotropic material failed in setTrialStrain.\n";
		return iso_res;
	}

	// move to orthotropic space
	isoToOrtho();

	// done
	return 0;
}

Vector & OrthotropicMaterial::getStrain(void)
{
	return epsilon;
}

Vector & OrthotropicMaterial::getStress(void)
{
	return sigma;
}

Matrix & OrthotropicMaterial::getTangent(void)
{
	return C;
}

Matrix & OrthotropicMaterial::getInitialTangent(void)
{
	return C0;
}

int OrthotropicMaterial::commitState(void)
{
	int iso_comState = theIsotropicMaterial->commitState();

	if (iso_comState != 0)
	{
		opserr << " Error: OrthotropicMaterial - the isotropic material failed in commitState.\n";
		return iso_comState;
	}

	/* move to orthotropic space before saving sigma and C.
	Why? some classes of material (IMPLEX) may change them in the commit state! */
	isoToOrtho();

	commit_C = C;
	commit_epsilon = epsilon;
	commit_sigma = sigma;

	return 0;
}

int OrthotropicMaterial::revertToLastCommit(void)
{
	int iso_revToLastCommit = theIsotropicMaterial->revertToLastCommit();
	if (iso_revToLastCommit != 0)
	{
		opserr << "Error: OrthotropicMaterial - the isotropic material failed in revertToLastCommit.\n";
		return iso_revToLastCommit;
	}

	C = commit_C;
	epsilon = commit_epsilon;
	sigma = commit_sigma;

	return 0;
}

int OrthotropicMaterial::revertToStart(void)
{
	epsilon.Zero();
	sigma.Zero();
	C.Zero();

	int iso_tevToStart = theIsotropicMaterial->revertToStart();
	if (iso_tevToStart!=0)
	{
		opserr << "Error: OrthotropicMaterial - the isotropic material failed in revertToStart.\n";
		return iso_tevToStart;
	}

	return 0;
}

NDMaterial * OrthotropicMaterial::getCopy(void)
{
	OrthotropicMaterial *theCopy = new OrthotropicMaterial(this->getTag(), *theIsotropicMaterial, C0, Aepsilon, Asigma);

	return theCopy;
}

NDMaterial * OrthotropicMaterial::getCopy(const char * type)
{
	if (strcmp(type, "ThreeDimensional") == 0 || strcmp(type, "3D") == 0) 
		return getCopy();
	return nullptr;
}

const char* OrthotropicMaterial::getType(void) const
{
	return "ThreeDimensional";
}

void OrthotropicMaterial::Print(OPS_Stream &s, int flag)
{
	s << "OrthotropicMaterial, tag: " << this->getTag() << "\n";
	s << "                 epsilon: " << epsilon << "\n";
	s << "                      C0: " << C0 << "\n";
	s << "                       C: " << C << "\n";
	s << "                   sigma: " << sigma << "\n";

	return;
}

int OrthotropicMaterial::sendSelf(int commitTag, Channel &theChannel)
{
	return 0;
}

int OrthotropicMaterial::recvSelf(int commitTag, Channel & theChannel, FEM_ObjectBroker & theBroker)
{
	return 0;
}

int OrthotropicMaterial::setParameter(const char** argv, int argc, Parameter& param)
{
	return theIsotropicMaterial->setParameter(argv, argc, param);
}

Response* OrthotropicMaterial::setResponse(const char** argv, int argc, OPS_Stream& s)
{
	// for strain, stress and tangent use the base class implementation
		// so that the output will be that of the adapter
	if (strcmp(argv[0], "Tangent") == 0 ||
		strcmp(argv[0], "tangent") == 0 ||
		strcmp(argv[0], "stress") == 0 ||
		strcmp(argv[0], "stresses") == 0 ||
		strcmp(argv[0], "strain") == 0 ||
		strcmp(argv[0], "strains") == 0
		) {
		return NDMaterial::setResponse(argv, argc, s);
	}
	// otherwise, for other custom results, forward the call to the adaptee
	return theIsotropicMaterial->setResponse(argv, argc, s);
}

void OrthotropicMaterial::isoToOrtho()
{
	// move back to orthotropic space
	const Vector &sig_iso = theIsotropicMaterial->getStress();
	const Matrix &C_iso = theIsotropicMaterial->getTangent();

	// invert Asigma, note: Asigma is a diagonal matrix!
	static Matrix Asigma_inv(6, 6);
	for (int i = 0; i < 6; i++)
		Asigma_inv(i, i) = 1.0 / Asigma(i, i);

	// tensor calculation of orthotropic material stresses
	sigma.addMatrixVector(0.0, Asigma_inv, sig_iso, 1.0);

	// compute orthotripic tangent
	static Matrix temp(6, 6);
	temp.addMatrixProduct(0.0, C_iso, Aepsilon, 1.0);   // temp = C_iso*Aepsilon
	C.addMatrixProduct(0.0, Asigma_inv, temp, 1.0);     //C = Asigma-1 *  C_iso * Aepsilon
}