/*
Developed by: Davide Noe Gorini (davideno.gorini@uniroma1.it)				  
Date: 2023
Description: This file contains the implementation of the TIM simulating the triaxial force-deformation relationship at the deck-abutment contact for soil-abutment systems
Reference: Gorini, D.N., Callisto, L., Whittle A.J., Sessa, S. (2023): "A multiaxial inertial macroelement for bridge abutments", International Journal for Numerical and Analytical Methods in Geomechanics, vol. 47(5), pp. 793-816, doi: 10.1002/nag.3493
*/

#ifndef TIMSoilAbutment3D_h
#define TIMSoilAbutment3D_h

#include <NDMaterial.h>

#include <T2Vector.h>
#include <Matrix.h>
#include <Vector.h>


class TIMSoilAbutment3D : public NDMaterial
{
public:

	TIMSoilAbutment3D(int tag,
		int nYield, double H11el, double H22el, double H33el, double aMult, double aiult, double alult, double delta, double Mc1, double Mc3, double tol, int niter);

	TIMSoilAbutment3D();
	virtual ~TIMSoilAbutment3D();

	const char* getClassType(void) const { return "TIMSoilAbutment3D"; };
	const char* getType(void) const { return "ThreeDimensional"; };
	int setTrialStrain(const Vector& strain);
	int setTrialStrain(const Vector& v, const Vector& r);
	int setTrialStrainIncr(const Vector& v);
	int setTrialStrainIncr(const Vector& v, const Vector& r);
	int getOrder() const; // reduction to 3 DoFs 20200219

						  // Calculates current tangent stiffness.
	const Matrix& getTangent(void);
	const Matrix& getInitialTangent(void);

	// Calculates the corresponding stress increment (rate), for a given strain increment. 
	const Vector& getStress(void);
	const Vector& getStrain(void);
	const Vector& getCommittedStress(void);
	const Vector& getCommittedStrain(void);


	int commitState(void);
	int revertToLastCommit(void);
	int revertToStart(void);

	NDMaterial* getCopy(void);

	NDMaterial* getCopy(const char* code);


	int sendSelf(int commitTag, Channel& theChannel);
	int recvSelf(int commitTag, Channel& theChannel,
		FEM_ObjectBroker& theBroker);

	Response* setResponse(const char** argv, int argc, OPS_Stream& s);
	int getResponse(int responseID, Information& matInformation);
	void Print(OPS_Stream& s, int flag = 0);

	int setParameter(const char** argv, int argc, Parameter& param);
	int updateParameter(int responseID, Information& eleInformation);

protected:

private:

	int ComputeResponse();

	int nYield;
	int niter;
	double H11el;
	double H22el;
	double H33el;

	Vector a3;
	Vector a2;
	Vector a1;

	double m1;
	double m2;
	double m3;

	Vector H11;
	Vector H22;
	Vector H33;

	double mh1;
	double mh2;
	double mh3;

	double tol;
	double delta;

	double q1el;
	double q2el;
	double q3el;

	Vector t13;
	Vector t23;
	Vector t12;
	Vector t13_2;
	Vector t23_2;
	Vector t12_2;

	// --- trial variables ----------------
	Vector Tstress;

	Vector Tstrain;

	int setFlows;

	int TrialStressElasticFlag;

	Vector Tc1;
	Vector Tc2;
	Vector Tc3;
	Vector Tc1Comm;
	Vector Tc2Comm;
	Vector Tc3Comm;

	Vector Qy1_iter;
	Vector Qy2_iter;
	Vector Qy3_iter;
	Vector Qy1_iter_2;
	Vector Qy2_iter_2;
	Vector Qy3_iter_2;

	int checkUnload;
	int unload;
	int numUnload;
	// ------------------------------------


	// ---- committed variables -----------
	Vector Cstress;

	Vector Cstrain;

	int CnFlows;
	int PrnFlows;

	double PStr1;
	double PStr2;
	double PStr3;
	double PStr1_2;
	double PStr2_2;
	double PStr3_2;
	double TStr1_2;
	double TStr2_2;
	double TStr3_2;
	double CStr1_2;
	double CStr2_2;
	double CStr3_2;
	double Peps1;
	double Peps2;
	double Peps3;

	Vector CQy1;
	Vector CQy2;
	Vector CQy3;
	Vector CQy1_2;
	Vector CQy2_2;
	Vector CQy3_2;

	int PcheckUnload;
	int Punload;
	int PnumUnload;

	Vector PCQy1;
	Vector PCQy2;
	Vector PCQy3;
	Vector PCQy1_2;
	Vector PCQy2_2;
	Vector PCQy3_2;

	Vector Cc1;
	Vector Cc2;
	Vector Cc3;

	Matrix PKt;

	int CommittedStressElasticFlag;

	// ------------------------------------

	Matrix theTangent;

};

#endif
