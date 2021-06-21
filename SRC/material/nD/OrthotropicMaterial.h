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

#ifndef OrthotropicMaterial_h
#define OrthotropicMaterial_h

#include <NDMaterial.h>
#include <Vector.h>
#include <Matrix.h>

#include<ID.h>

class OrthotropicMaterial : public NDMaterial 
{
public:
	OrthotropicMaterial(int tag, NDMaterial &theIsoMat,
		double Ex, double Ey, double Ez, double Gxy, double Gyz, double Gzx,
		double vxy, double vyz, double vzx,
		double Asigmaxx, double Asigmayy, double Asigmazz, double Asigmaxyxy, double Asigmayzyz, double Asigmaxzxz);

	OrthotropicMaterial(int tag_copy, NDMaterial &theIsoMat_copy, Matrix C0_copy, Matrix Aepsilon_copy, Matrix Asigma_copy);
	OrthotropicMaterial();
	virtual ~OrthotropicMaterial();

	virtual int setTrialStrain(const Vector &strain);

	virtual Vector &getStrain(void);
	virtual Vector &getStress(void);

	virtual Matrix &getTangent(void);
	virtual Matrix &getInitialTangent(void);

	virtual int commitState(void);
	virtual int revertToLastCommit(void);
	virtual int revertToStart(void);

	virtual NDMaterial *getCopy(void);
	virtual NDMaterial *getCopy(const char *type);
	virtual const char *getType(void) const;

	virtual void Print(OPS_Stream &s, int flag=0);

	virtual int sendSelf(int commitTag, Channel &theChannel);
	int recvSelf(int commitTag, Channel &theChannel, FEM_ObjectBroker &theBroker);

	int setParameter(const char** argv, int argc, Parameter& param);
	virtual Response* setResponse(const char** argv, int argc, OPS_Stream& s);

private:
	void isoToOrtho();

private:
	NDMaterial *theIsotropicMaterial;

	Matrix C;
	Matrix C0;
	Vector epsilon;
	Vector sigma;

	Matrix Aepsilon;
	Matrix Asigma;

	Matrix commit_C;
	Vector commit_epsilon;
	Vector commit_sigma;
};
#endif