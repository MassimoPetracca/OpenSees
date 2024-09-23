/*
Developed by: Davide Noe Gorini (davideno.gorini@uniroma1.it)
Date: 2023
Description: This file contains the implementation of the TIM simulating the triaxial force-deformation relationship at the deck-abutment contact for soil-abutment systems
Reference: Gorini, D.N., Callisto, L., Whittle A.J., Sessa, S. (2023): "A multiaxial inertial macroelement for bridge abutments", International Journal for Numerical and Analytical Methods in Geomechanics, vol. 47(5), pp. 793-816, doi: 10.1002/nag.3493
*/

#include <math.h>
#include <stdlib.h>
#include <TIMSoilAbutment3D.h>
#include <Information.h>
#include <ID.h>
#include <MaterialResponse.h>
#include <Parameter.h>

#include <elementAPI.h>


void* OPS_TIMSoilAbutment3D(void) {

	int tag, nYield, niter;
	double H11el, H22el, H33el, aMult, aiult, alult, delta, Mc1, Mc3, tol;

	int numArgs = OPS_GetNumRemainingInputArgs();
	if (numArgs != 13) {
		opserr << "ndMaterial TIMSoilAbutment3D incorrect num args: want tag nYield H11 H22 H33 aM ai al delta Mc1 Mc3 tol niter\n";
		return 0;
	}

	int iData[2];
	double dData[10];

	int numData = 2;
	if (OPS_GetInt(&numData, iData) != 0) {
		opserr << "WARNING invalid integer values: nDMaterial TIMSoilAbutment3D \n";
		return 0;
	}
	tag = iData[0];
	nYield = iData[1];

	if (nYield > 20) {
		opserr << "WARNING too many yield surfaces: nDMaterial TIMSoilAbutment3D; nYield must be <=20 \n";
		return 0;
	}


	numData = 10;
	if (OPS_GetDouble(&numData, dData) != 0) {
		opserr << "WARNING invalid double values: nDMaterial TIMSoilAbutment3D " << tag << endln;
		return 0;
	}
	H11el = dData[0];
	H22el = dData[1];
	H33el = dData[2];
	aMult = dData[3];
	aiult = dData[4];
	alult = dData[5];
	delta = dData[6];
	Mc1 = dData[7];
	Mc3 = dData[8];
	tol = dData[9];


	numData = 1;
	int iterData[1];
	if (OPS_GetInt(&numData, iterData) != 0) {
		opserr << "WARNING invalid integer values: nDMaterial TIMSoilAbutment3D \n";
		return 0;
	}
	niter = iterData[0];

	NDMaterial* theMaterial = new TIMSoilAbutment3D(tag,
		nYield, H11el, H22el, H33el, aMult, aiult, alult, delta, Mc1, Mc3, tol, niter);
	
	return theMaterial;
}

TIMSoilAbutment3D::TIMSoilAbutment3D(int ptag,
	int pnYield, double pH11el, double pH22el, double pH33el, double paMult, double paiult, double palult, double pdelta, double pMc1, 
	double pMc3, double ptol, int pniter)
	: NDMaterial(ptag, ND_TAG_TIMSoilAbutment3D), Tstress(6),
	Tstrain(6), Cstress(6), Cstrain(6), theTangent(6, 6), PKt(6, 6), a3(pnYield), a2(pnYield),
	a1(pnYield), H11(pnYield), H22(pnYield), H33(pnYield),
	Tc1(pnYield), Tc2(pnYield), Tc3(pnYield), Tc1Comm(pnYield), Tc2Comm(pnYield), Tc3Comm(pnYield), Cc1(pnYield), Cc2(pnYield), Cc3(pnYield),
	CQy1(pnYield), CQy2(pnYield), CQy3(pnYield), CQy1_2(pnYield), CQy2_2(pnYield), CQy3_2(pnYield),
	PCQy1(pnYield), PCQy2(pnYield), PCQy3(pnYield),
	PCQy1_2(pnYield), PCQy2_2(pnYield), PCQy3_2(pnYield),
	Qy1_iter(pnYield), Qy2_iter(pnYield), Qy3_iter(pnYield), Qy1_iter_2(pnYield), Qy2_iter_2(pnYield), Qy3_iter_2(pnYield),
	t13(3600), t23(360), t12(360), t13_2(7199), t23_2(719), t12_2(719)
{
	this->nYield = pnYield;
	this->niter = pniter;
	this->tol = ptol;
	this->H11el = pH11el;
	this->H22el = pH22el;
	this->H33el = pH33el;

	this->CnFlows = 0;
	this->PrnFlows = 0;
	this->PStr1 = 0.0;
	this->PStr2 = 0.0;
	this->PStr3 = 0.0;
	this->Peps1 = 0.0;
	this->Peps2 = 0.0;
	this->Peps3 = 0.0;
	this->PStr1_2 = 0.0;
	this->PStr2_2 = 0.0;
	this->PStr3_2 = 0.0;
	this->TStr1_2 = 0.0;
	this->TStr2_2 = 0.0;
	this->TStr3_2 = 0.0;
	this->CStr1_2 = 0.0;
	this->CStr2_2 = 0.0;
	this->CStr3_2 = 0.0;

	CQy1.Zero();
	CQy2.Zero();
	CQy3.Zero();
	CQy1_2.Zero();
	CQy2_2.Zero();
	CQy3_2.Zero();

	PcheckUnload = 0;
	Punload = 0;
	PnumUnload = 0;

	checkUnload = 0;
	numUnload = 0;
	q1el = 0.0;
	q2el = 0.0;
	q3el = 0.0;
	unload = 0;

	Qy1_iter.Zero();
	Qy2_iter.Zero();
	Qy3_iter.Zero();
	Qy1_iter_2.Zero();
	Qy2_iter_2.Zero();
	Qy3_iter_2.Zero();
	this->setFlows = 0;

	PCQy1.Zero();
	PCQy2.Zero();
	PCQy3.Zero();
	PCQy1_2.Zero();
	PCQy2_2.Zero();
	PCQy3_2.Zero();

	Tstress.Zero();
	Tstrain.Zero();

	Cstress.Zero();
	Cstrain.Zero();

	theTangent.Zero();
	theTangent(0, 0) = H11el;
	theTangent(1, 1) = H22el;
	theTangent(2, 2) = H33el;

	PKt.Zero();
	PKt(0, 0) = H11el;
	PKt(1, 1) = H22el;
	PKt(2, 2) = H33el;

	this->delta = pdelta / 180 * 4 * atan(1.0);

	this->a3(nYield - 1) = paMult;
	this->a2(nYield - 1) = paiult;
	this->a1(nYield - 1) = palult;

	double Mc1 = pMc1;
	double Mc3 = pMc3;

	this->Tc3(nYield - 1) = Mc3 * a3(nYield - 1);
	this->Tc2(nYield - 1) = 0.0;
	this->Tc1(nYield - 1) = Tc3(nYield - 1) / Mc1;

	this->a3(0) = a3(nYield - 1) * .1;
	this->a2(0) = a2(nYield - 1) * .1;
	this->a1(0) = a1(nYield - 1) * .1;

	this->Tc3(0) = Tc3(nYield - 1) * .1;
	this->Tc2(0) = Tc2(nYield - 1) * .1;
	this->Tc1(0) = Tc1(nYield - 1) * .1;

	// increasing size of the yield surfaces, from the innermost to the outermost
	m1 = (a1(nYield - 1) - a1(0)) / (nYield - 1);
	m2 = (a2(nYield - 1) - a2(0)) / (nYield - 1);
	m3 = (a3(nYield - 1) - a3(0)) / (nYield - 1);

	double mc1 = (Tc1(nYield - 1) - Tc1(0)) / (nYield - 1);
	double mc2 = (Tc2(nYield - 1) - Tc2(0)) / (nYield - 1);
	double mc3 = (Tc3(nYield - 1) - Tc3(0)) / (nYield - 1);

	for (int j = 1; j < nYield; j++)
	{
		this->a1(j) = a1(j - 1) + m1;
		this->a2(j) = a2(j - 1) + m2;
		this->a3(j) = a3(j - 1) + m3;
		this->Tc1(j) = Tc1(j - 1) + mc1;
		this->Tc2(j) = Tc2(j - 1) + mc2;
		this->Tc3(j) = Tc3(j - 1) + mc3;
	}

	// translation of the centers of the yield surfaces
	for (int j = 0; j < nYield; j++)
	{
		this->Tc1(j) = Tc1(j);
		this->Tc2(j) = Tc2(j);
		this->Tc3(j) = Tc3(j);
		this->Cc1(j) = Tc1(j);
		this->Cc2(j) = Tc2(j);
		this->Cc3(j) = Tc3(j);
	}
	for (int j = 0; j < nYield; j++)
	{
		this->Tc1Comm(j) = Tc1(j);
		this->Tc2Comm(j) = Tc2(j);
		this->Tc3Comm(j) = Tc3(j);
	}

	this->H11(0) = H11el;
	this->H22(0) = H22el;
	this->H33(0) = H33el;

	this->H11(nYield - 1) = 0.001 * H11(0);
	this->H22(nYield - 1) = 0.001 * H22(0);
	this->H33(nYield - 1) = 0.001 * H33(0);

	mh1 = (H11(nYield - 1) - H11(0)) / (nYield - 1);
	mh2 = (H22(nYield - 1) - H22(0)) / (nYield - 1);
	mh3 = (H33(nYield - 1) - H33(0)) / (nYield - 1);

	for (int j = 1; j < nYield; j++)
	{
		this->H11(j) = H11(j - 1) + mh1;
		this->H22(j) = H22(j - 1) + mh2;
		this->H33(j) = H33(j - 1) + mh3;
	}

	this->TrialStressElasticFlag = 0;
	this->CommittedStressElasticFlag = 0;

}

TIMSoilAbutment3D::TIMSoilAbutment3D() : NDMaterial(0, ND_TAG_TIMSoilAbutment3D), Tstress(6),
Tstrain(6), Cstress(6), Cstrain(6), theTangent(6, 6), PKt(6, 6), a3(nYield), a2(nYield),
a1(nYield), H11(nYield), H22(nYield), H33(nYield),
Tc1(nYield), Tc2(nYield), Tc3(nYield), Tc1Comm(nYield), Tc2Comm(nYield), Tc3Comm(nYield), Cc1(nYield), Cc2(nYield), Cc3(nYield),
CQy1(nYield), CQy2(nYield), CQy3(nYield), CQy1_2(nYield), CQy2_2(nYield), CQy3_2(nYield),
PCQy1(nYield), PCQy2(nYield), PCQy3(nYield),
PCQy1_2(nYield), PCQy2_2(nYield), PCQy3_2(nYield),
Qy1_iter(nYield), Qy2_iter(nYield), Qy3_iter(nYield), Qy1_iter_2(nYield), Qy2_iter_2(nYield), Qy3_iter_2(nYield),
t13(3600), t23(360), t12(360), t13_2(7199), t23_2(719), t12_2(719)
{
	this->nYield = 0.0;
	this->niter = 0.0;
	this->tol = 0.0;

	this->H11el = 0.0;
	this->H22el = 0.0;
	this->H33el = 0.0;

	Tstress.Zero();
	Tstrain.Zero();

	Cstress.Zero();
	Cstrain.Zero();

	PcheckUnload = 0;
	Punload = 0;
	PnumUnload = 0;

	this->CnFlows = 0;
	this->PrnFlows = 0;
	this->PStr1 = 0.0;
	this->PStr2 = 0.0;
	this->PStr3 = 0.0;
	this->Peps1 = 0.0;
	this->Peps2 = 0.0;
	this->Peps3 = 0.0;
	this->PStr1_2 = 0.0;
	this->PStr2_2 = 0.0;
	this->PStr3_2 = 0.0;
	this->TStr1_2 = 0.0;
	this->TStr2_2 = 0.0;
	this->TStr3_2 = 0.0;
	this->CStr1_2 = 0.0;
	this->CStr2_2 = 0.0;
	this->CStr3_2 = 0.0;
	CQy1.Zero();
	CQy2.Zero();
	CQy3.Zero();
	CQy1_2.Zero();
	CQy2_2.Zero();
	CQy3_2.Zero();
	PCQy1.Zero();
	PCQy2.Zero();
	PCQy3.Zero();
	PCQy1_2.Zero();
	PCQy2_2.Zero();
	PCQy3_2.Zero();

	Qy1_iter.Zero();
	Qy2_iter.Zero();
	Qy3_iter.Zero();
	Qy1_iter_2.Zero();
	Qy2_iter_2.Zero();
	Qy3_iter_2.Zero();
	this->setFlows = 0;

	theTangent.Zero();

	PKt.Zero();

	this->delta = 0.0;

	for (int j = 0; j < nYield; j++)
	{
		this->a1(j) = 0.0;
		this->a2(j) = 0.0;
		this->a3(j) = 0.0;
	}

	// translation of the centers of the yield surfaces
	for (int j = 0; j < nYield; j++)
	{
		this->Tc1(j) = 0.0;
		this->Tc2(j) = 0.0;
		this->Tc3(j) = 0.0;
		this->Cc1(j) = Tc1(j);
		this->Cc2(j) = Tc2(j);
		this->Cc3(j) = Tc3(j);
		this->Tc1Comm(j) = 0.0;
		this->Tc2Comm(j) = 0.0;
		this->Tc3Comm(j) = 0.0;
	}

	for (int j = 0; j < nYield; j++)
	{
		this->H11(j) = 0.0;
		this->H22(j) = 0.0;
		this->H33(j) = 0.0;
	}

	this->TrialStressElasticFlag = 0;
	this->CommittedStressElasticFlag = 0;

}


TIMSoilAbutment3D::~TIMSoilAbutment3D() {

	return;
};


int TIMSoilAbutment3D::setTrialStrain(const Vector& pStrain) {


	Tstrain = pStrain;
	double checkpStrain0 = pStrain(0);
	double checkpStrain1 = pStrain(1);
	double checkpStrain2 = pStrain(2);
	double checkpStrain3 = pStrain(3);
	double checkpStrain4 = pStrain(4);
	double checkpStrain5 = pStrain(5);

	Tc1 = Cc1;
	Tc2 = Cc2;
	Tc3 = Cc3;

	this->ComputeResponse();

	return 0;

};

int TIMSoilAbutment3D::setTrialStrain(const Vector& v, const Vector& r) {

	return this->setTrialStrain(v);

};

int TIMSoilAbutment3D::setTrialStrainIncr(const Vector& v) {

	Tstrain = Cstrain + v;

	Tc1 = Cc1;
	Tc2 = Cc2;
	Tc3 = Cc3;

	this->ComputeResponse();
	return 0;

};

int TIMSoilAbutment3D::getOrder() const
{
	return 3;
}

int TIMSoilAbutment3D::setTrialStrainIncr(const Vector& v, const Vector& r) {

	return this->setTrialStrainIncr(v);


};


const Matrix& TIMSoilAbutment3D::getTangent() {
	
	return this->theTangent;
};


const Matrix& TIMSoilAbutment3D::getInitialTangent() {

	static Matrix initTangent(6, 6);
	initTangent.Zero();
	double initTangent0 = initTangent(0, 0);
	initTangent(0, 0) = H11el;
	initTangent(1, 1) = H22el;
	initTangent(2, 2) = H33el;
	return initTangent;
};


const Vector& TIMSoilAbutment3D::getStress(void) {

	return Tstress;

};

const Vector& TIMSoilAbutment3D::getStrain(void) {

	return Tstrain;
};

const Vector& TIMSoilAbutment3D::getCommittedStress(void) {

	return Cstress;
};

const Vector& TIMSoilAbutment3D::getCommittedStrain(void) {

	return Cstrain;
	double checkgetCommittedStrain = Cstrain(0);

};


int TIMSoilAbutment3D::commitState(void) {

	CnFlows = setFlows;
	PrnFlows = CnFlows;

	Cstress = Tstress;

	Cstrain = Tstrain;
	PStr1 = Cstress(0);
	PStr2 = Cstress(1);
	PStr3 = Cstress(2);
	Peps1 = Cstrain(0);
	Peps2 = Cstrain(1);
	Peps3 = Cstrain(2);

	CStr1_2 = TStr1_2;
	CStr2_2 = TStr2_2;
	CStr3_2 = TStr3_2;

	PStr1_2 = CStr1_2;
	PStr2_2 = CStr2_2;
	PStr3_2 = CStr3_2;

	Cc1 = Tc1Comm;
	Cc2 = Tc2Comm;
	Cc3 = Tc3Comm;

	PcheckUnload = checkUnload;
	Punload = unload;
	PnumUnload = numUnload;

	PCQy1 = CQy1;
	PCQy2 = CQy2;
	PCQy3 = CQy3;
	PCQy1_2 = CQy1_2;
	PCQy2_2 = CQy2_2;
	PCQy3_2 = CQy3_2;

	double A1y0 = PCQy1(0);
	double A1yN = PCQy1(nYield - 1);

	PKt = theTangent;
	
	CommittedStressElasticFlag = TrialStressElasticFlag;

	double checkcommitStateCstress = Cstress(0);
	double checkcommitStateCstrain = Cstrain(0);
	double checkcommitStateCc1 = Cc1(0);
	double checkcommitStateCc2 = Cc2(0);
	double checkcommitStateCc3 = Cc3(0);

	return 0;

};

int TIMSoilAbutment3D::revertToLastCommit(void) {

	Tstress = Cstress;

	Tstrain = Cstrain;

	Tc1 = Cc1;
	Tc2 = Cc2;
	Tc3 = Cc3;

	TrialStressElasticFlag = CommittedStressElasticFlag;

	return 0;
};

int TIMSoilAbutment3D::revertToStart(void) {
	Tstress.Zero();

	Tstrain.Zero();

	Tc1.Zero();
	Tc2.Zero();
	Tc3.Zero();

	TrialStressElasticFlag = 0;

	return 0;
}



NDMaterial* TIMSoilAbutment3D::getCopy(void) {

	TIMSoilAbutment3D* theSA = new TIMSoilAbutment3D(*this);
	return theSA;


};



NDMaterial* TIMSoilAbutment3D::getCopy(const char* code)
{
	if (strcmp(code, "ThreeDimensional") == 0) {
		TIMSoilAbutment3D* theSA = new TIMSoilAbutment3D(*this);
		return theSA;
	}

	return 0;
}


int TIMSoilAbutment3D::sendSelf(int commitTag, Channel& theChannel) {
	
	int res = 0;
	
	static Vector data(123 + 24 * nYield + 9);

	data(0) = this->getTag();
	data(1) = nYield;
	data(2) = niter;
	data(3) = tol;
	data(4) = H11el;
	data(5) = H22el;
	data(6) = H33el;
	data(7) = delta;
	data(8) = m1;
	data(9) = m2;
	data(10) = m3;
	data(11) = mh1;
	data(12) = mh2;
	data(13) = mh3;
	data(14) = TrialStressElasticFlag;
	data(15) = CommittedStressElasticFlag;

	data(16) = Tstress(0);
	data(17) = Tstress(1);
	data(18) = Tstress(2);
	data(19) = Tstress(3);
	data(20) = Tstress(4);
	data(21) = Tstress(5);

	data(22) = Tstrain(0);
	data(23) = Tstrain(1);
	data(24) = Tstrain(2);
	data(25) = Tstrain(3);
	data(26) = Tstrain(4);
	data(27) = Tstrain(5);

	data(28) = Cstress(0);
	data(29) = Cstress(1);
	data(30) = Cstress(2);
	data(31) = Cstress(3);
	data(32) = Cstress(4);
	data(33) = Cstress(5);

	data(34) = Cstrain(0);
	data(35) = Cstrain(1);
	data(36) = Cstrain(2);
	data(37) = Cstrain(3);
	data(38) = Cstrain(4);
	data(39) = Cstrain(5);

	data(40) = CnFlows;
	data(41) = PrnFlows;
	data(42) = PStr1;
	data(43) = PStr2;
	data(44) = PStr3;
	data(45) = Peps1;
	data(46) = Peps2;
	data(47) = Peps3;
	data(48) = PcheckUnload;
	data(49) = Punload;
	data(50) = PnumUnload;

	int idata = 50;
	for (int i = 0; i < 6; i++) {
		for (int j = 0; j < 6; j++) {
			idata = idata + 1;
			data(idata) = theTangent(i, j);
		}
	}
	for (int i = 0; i < 6; i++) {
		for (int j = 0; j < 6; j++) {
			idata = idata + 1;
			data(idata) = PKt(i, j);
		}
	}
	for (int i = 0; i < nYield; i++) {
		idata = idata + 1;
		data(idata) = a3(i);
	}
	for (int i = 0; i < nYield; i++) {
		idata = idata + 1;
		data(idata) = a2(i);
	}
	for (int i = 0; i < nYield; i++) {
		idata = idata + 1;
		data(idata) = a1(i);
	}
	for (int i = 0; i < nYield; i++) {
		idata = idata + 1;
		data(idata) = H11(i);
	}
	for (int i = 0; i < nYield; i++) {
		idata = idata + 1;
		data(idata) = H22(i);
	}
	for (int i = 0; i < nYield; i++) {
		idata = idata + 1;
		data(idata) = H33(i);
	}
	for (int i = 0; i < nYield; i++) {
		idata = idata + 1;
		data(idata) = Tc1(i);
	}
	for (int i = 0; i < nYield; i++) {
		idata = idata + 1;
		data(idata) = Tc2(i);
	}
	for (int i = 0; i < nYield; i++) {
		idata = idata + 1;
		data(idata) = Tc3(i);
	}
	for (int i = 0; i < nYield; i++) {
		idata = idata + 1;
		data(idata) = CQy1(i);
	}
	for (int i = 0; i < nYield; i++) {
		idata = idata + 1;
		data(idata) = CQy2(i);
	}
	for (int i = 0; i < nYield; i++) {
		idata = idata + 1;
		data(idata) = CQy3(i);
	}
	for (int i = 0; i < nYield; i++) {
		idata = idata + 1;
		data(idata) = CQy1_2(i);
	}
	for (int i = 0; i < nYield; i++) {
		idata = idata + 1;
		data(idata) = CQy2_2(i);
	}
	for (int i = 0; i < nYield; i++) {
		idata = idata + 1;
		data(idata) = CQy3_2(i);
	}
	for (int i = 0; i < nYield; i++) {
		idata = idata + 1;
		data(idata) = PCQy1(i);
	}
	for (int i = 0; i < nYield; i++) {
		idata = idata + 1;
		data(idata) = PCQy2(i);
	}
	for (int i = 0; i < nYield; i++) {
		idata = idata + 1;
		data(idata) = PCQy3(i);
	}
	for (int i = 0; i < nYield; i++) {
		idata = idata + 1;
		data(idata) = PCQy1_2(i);
	}
	for (int i = 0; i < nYield; i++) {
		idata = idata + 1;
		data(idata) = PCQy2_2(i);
	}
	for (int i = 0; i < nYield; i++) {
		idata = idata + 1;
		data(idata) = PCQy3_2(i);
	}

	for (int i = 0; i < nYield; i++) {
		idata = idata + 1;
		data(idata) = Tc1Comm(i);
	}
	for (int i = 0; i < nYield; i++) {
		idata = idata + 1;
		data(idata) = Tc2Comm(i);
	}
	for (int i = 0; i < nYield; i++) {
		idata = idata + 1;
		data(idata) = Tc3Comm(i);
	}

	idata = idata + 1;
	data(idata) = PStr1_2;
	idata = idata + 1;
	data(idata) = PStr2_2;
	idata = idata + 1;
	data(idata) = PStr3_2;

	idata = idata + 1;
	data(idata) = TStr1_2;
	idata = idata + 1;
	data(idata) = TStr2_2;
	idata = idata + 1;
	data(idata) = TStr3_2;

	idata = idata + 1;
	data(idata) = CStr1_2;
	idata = idata + 1;
	data(idata) = CStr2_2;
	idata = idata + 1;
	data(idata) = CStr3_2;

	res = theChannel.sendVector(this->getDbTag(), commitTag, data);

	return 0;
};

int TIMSoilAbutment3D::recvSelf(int commitTag, Channel& theChannel, FEM_ObjectBroker& theBroker) {
	
	static Vector data(123 + 24 * nYield + 9);
	int res = theChannel.recvVector(0, commitTag, data);
	this->setTag(data(0));
	nYield = data(1);
	niter = data(2);
	tol = data(3);
	H11el = data(4);
	H22el = data(5);
	H33el = data(6);
	delta = data(7);
	m1 = data(8);
	m2 = data(9);
	m3 = data(10);
	mh1 = data(11);
	mh2 = data(12);
	mh3 = data(13);
	TrialStressElasticFlag = data(14);
	CommittedStressElasticFlag = data(15);

	Tstress(0) = data(16);
	Tstress(1) = data(17);
	Tstress(2) = data(18);
	Tstress(3) = data(19);
	Tstress(4) = data(20);
	Tstress(5) = data(21);

	Tstrain(0) = data(22);
	Tstrain(1) = data(23);
	Tstrain(2) = data(24);
	Tstrain(3) = data(25);
	Tstrain(4) = data(26);
	Tstrain(5) = data(27);

	Cstress(0) = data(28);
	Cstress(1) = data(29);
	Cstress(2) = data(30);
	Cstress(3) = data(31);
	Cstress(4) = data(32);
	Cstress(5) = data(33);

	Cstrain(0) = data(34);
	Cstrain(1) = data(35);
	Cstrain(2) = data(36);
	Cstrain(3) = data(37);
	Cstrain(4) = data(38);
	Cstrain(5) = data(39);

	CnFlows = data(40);
	PrnFlows = data(41);
	PStr1 = data(42);
	PStr2 = data(43);
	PStr3 = data(44);
	Peps1 = data(45);
	Peps2 = data(46);
	Peps3 = data(47);
	PcheckUnload = data(48);
	Punload = data(49);
	PnumUnload = data(50);

	int idata = 50;
	for (int i = 0; i < 6; i++) {
		for (int j = 0; j < 6; j++) {
			idata = idata + 1;
			theTangent(i, j) = data(idata);
		}
	};
	for (int i = 0; i < 6; i++) {
		for (int j = 0; j < 6; j++) {
			idata = idata + 1;
			PKt(i, j) = data(idata);
		}
	};
	for (int i = 0; i < nYield; i++) {
		idata = idata + 1;
		a3(i) = data(idata);
	}
	for (int i = 0; i < nYield; i++) {
		idata = idata + 1;
		a2(i) = data(idata);
	}
	for (int i = 0; i < nYield; i++) {
		idata = idata + 1;
		a1(i) = data(idata);
	}
	for (int i = 0; i < nYield; i++) {
		idata = idata + 1;
		H11(i) = data(idata);
	}
	for (int i = 0; i < nYield; i++) {
		idata = idata + 1;
		H22(i) = data(idata);
	}
	for (int i = 0; i < nYield; i++) {
		idata = idata + 1;
		H33(i) = data(idata);
	}
	for (int i = 0; i < nYield; i++) {
		idata = idata + 1;
		Tc1(i) = data(idata);
	}
	for (int i = 0; i < nYield; i++) {
		idata = idata + 1;
		Tc2(i) = data(idata);
	}
	for (int i = 0; i < nYield; i++) {
		idata = idata + 1;
		Tc3(i) = data(idata);
	}
	for (int i = 0; i < nYield; i++) {
		idata = idata + 1;
		CQy1(i) = data(idata);
	}
	for (int i = 0; i < nYield; i++) {
		idata = idata + 1;
		CQy2(i) = data(idata);
	}
	for (int i = 0; i < nYield; i++) {
		idata = idata + 1;
		CQy3(i) = data(idata);
	}
	for (int i = 0; i < nYield; i++) {
		idata = idata + 1;
		CQy1_2(i) = data(idata);
	}
	for (int i = 0; i < nYield; i++) {
		idata = idata + 1;
		CQy2_2(i) = data(idata);
	}
	for (int i = 0; i < nYield; i++) {
		idata = idata + 1;
		CQy3_2(i) = data(idata);
	}
	for (int i = 0; i < nYield; i++) {
		idata = idata + 1;
		PCQy1(i) = data(idata);
	}
	for (int i = 0; i < nYield; i++) {
		idata = idata + 1;
		PCQy2(i) = data(idata);
	}
	for (int i = 0; i < nYield; i++) {
		idata = idata + 1;
		PCQy3(i) = data(idata);
	}
	for (int i = 0; i < nYield; i++) {
		idata = idata + 1;
		PCQy1_2(i) = data(idata);
	}
	for (int i = 0; i < nYield; i++) {
		idata = idata + 1;
		PCQy2_2(i) = data(idata);
	}
	for (int i = 0; i < nYield; i++) {
		idata = idata + 1;
		PCQy3_2(i) = data(idata);
	};

	for (int i = 0; i < nYield; i++) {
		idata = idata + 1;
		Tc1Comm(i) = data(idata);
	}
	for (int i = 0; i < nYield; i++) {
		idata = idata + 1;
		Tc2Comm(i) = data(idata);
	}
	for (int i = 0; i < nYield; i++) {
		idata = idata + 1;
		Tc3Comm(i) = data(idata);
	}

	idata = idata + 1;
	PStr1_2 = data(idata);
	idata = idata + 1;
	PStr2_2 = data(idata);
	idata = idata + 1;
	PStr3_2 = data(idata);

	idata = idata + 1;
	TStr1_2 = data(idata);
	idata = idata + 1;
	TStr2_2 = data(idata);
	idata = idata + 1;
	TStr3_2 = data(idata);

	idata = idata + 1;
	CStr1_2 = data(idata);
	idata = idata + 1;
	CStr2_2 = data(idata);
	idata = idata + 1;
	CStr3_2 = data(idata);

	return 0;
};

Response*
TIMSoilAbutment3D::setResponse(const char** argv, int argc, OPS_Stream& s)
{

	if (strcmp(argv[0], "stress") == 0 || strcmp(argv[0], "stresses") == 0)
		return new MaterialResponse(this, 1, this->getStress());

	else if (strcmp(argv[0], "strain") == 0 || strcmp(argv[0], "strains") == 0)
		return new MaterialResponse(this, 2, this->getStrain());

	else if (strcmp(argv[0], "tangent") == 0 || strcmp(argv[0], "Tangent") == 0)
		return new MaterialResponse(this, 3, this->getTangent());
	else
		return 0;

}



int TIMSoilAbutment3D::getResponse(int responseID, Information& matInfo) {

	switch (responseID) {
	case -1:
		return -1;
	case 1:
		if (matInfo.theVector != 0)
			*(matInfo.theVector) = getStress();
		return 0;

	case 2:
		if (matInfo.theVector != 0)
			*(matInfo.theVector) = getStrain();
		return 0;

	case 3:
		if (matInfo.theMatrix != 0)
			*(matInfo.theMatrix) = getTangent();
		return 0;
	}

	return 0;
};

void TIMSoilAbutment3D::Print(OPS_Stream& s, int flag) {
	s << "TIMSoilAbutment3D tag: " << this->getTag() << endln;
	s << "  nYield: " << nYield << endln;
	s << "  H11el: " << H11el << endln;
	s << "  H22el: " << H22el << endln;
	s << "  H33el: " << H33el << endln;
	s << " Yield Centers 1:" << Cc1 << endln;
	s << " Yield Centers 2:" << Cc2 << endln;
	s << " Yield Centers 3:" << Cc3 << endln;
	s << " Active Yield Surfaces:" << CnFlows << endln;
	s << "  stress: " << Tstress << " tangent: " << theTangent << endln;
	return;
};


int TIMSoilAbutment3D::setParameter(const char** argv, int argc, Parameter& param) {
	// -- to be implemented if necessary

	return 0;
};

int TIMSoilAbutment3D::updateParameter(int responseID, Information& eleInformation) {
	// -- to be implemented if necessary
	return 0;
};


int TIMSoilAbutment3D::ComputeResponse() {
	Vector csi1_iter(nYield), csi2_iter(nYield), csi3_iter(nYield), Qy1(nYield),
		Qy2(nYield), Qy3(nYield), Qy1_iter(nYield), Qy2_iter(nYield), Qy3_iter(nYield),
		Qy1_iter_2(nYield), Qy2_iter_2(nYield), Qy3_iter_2(nYield), Rn(nYield),
		Rn_iter(nYield), Rn_iter_2(nYield),
		Fn1(nYield), Fn2(nYield), Fn3(nYield), Fn4(nYield), Fy1(nYield), Fy2(nYield),
		Fy3(nYield), Fy4(nYield), Fd(nYield), DQ1(nYield), DQ2(nYield), DQ3(nYield),
		csi1(nYield), csi2(nYield), csi3(nYield), nFlows(niter),
		Qy1_iterS(nYield), Qy2_iterS(nYield), Qy3_iterS(nYield), Qy1_iter_2S(nYield),
		Qy2_iter_2S(nYield), Qy3_iter_2S(nYield);

	csi1.Zero();
	csi2.Zero();
	csi3.Zero();
	csi1_iter.Zero();
	csi2_iter.Zero();
	csi3_iter.Zero();
	Qy1.Zero();
	Qy2.Zero();
	Qy3.Zero();
	Qy1_iter.Zero();
	Qy2_iter.Zero();
	Qy3_iter.Zero();
	Qy1_iter_2.Zero();
	Qy2_iter_2.Zero();
	Qy3_iter_2.Zero();
	Qy1_iterS.Zero();
	Qy2_iterS.Zero();
	Qy3_iterS.Zero();
	Qy1_iter_2S.Zero();
	Qy2_iter_2S.Zero();
	Qy3_iter_2S.Zero();
	Rn.Zero();
	Rn_iter.Zero();
	Rn_iter_2.Zero();
	Fn1.Zero();
	Fn2.Zero();
	Fn3.Zero();
	Fn4.Zero();
	Fy1.Zero();
	Fy2.Zero();
	Fy3.Zero();
	Fy4.Zero();
	Fd.Zero();
	DQ1.Zero();
	DQ2.Zero();
	DQ3.Zero();
	nFlows.Zero();

	double t1_iterS = 0.0;
	double t2_iterS = 0.0;
	double t3_iterS = 0.0;

	double signtgang31 = 0.0;
	double signtgang31_2 = 0.0;
	double signRQ_2Rt = 0.0;
	double signRQ_2Rtm = 0.0;
	double tgang31 = 0.0;
	double tgang31_2 = 0.0;
	double signQ3_eq13tim13 = 0.0;
	double signQ1_eqtim13 = 0.0;
	double signQ2_eqtim13 = 0.0;
	double signQ2_eq_2 = 0.0;
	double signRQRt = 0.0;
	double signRQRtm = 0.0;

	Vector signQy1S(nYield), signQy2S(nYield), signQy3S(nYield), signQy1_2S(nYield), signQy2_2S(nYield), signQy3_2S(nYield);
	signQy1S.Zero();
	signQy2S.Zero();
	signQy3S.Zero();
	signQy1_2S.Zero();
	signQy2_2S.Zero();
	signQy3_2S.Zero();

	Vector signQy1(nYield), signQy2(nYield), signQy3(nYield), signQy1_2(nYield), signQy2_2(nYield), signQy3_2(nYield);
	signQy1.Zero();
	signQy2.Zero();
	signQy3.Zero();
	signQy1_2.Zero();
	signQy2_2.Zero();
	signQy3_2.Zero();

	double H12el = 0.0;
	double H13el = 0.0;
	double H14el = 0.0;
	double H15el = 0.0;
	double H21el = 0.0;
	double H23el = 0.0;
	double H24el = 0.0;
	double H25el = 0.0;
	double H31el = 0.0;
	double H32el = 0.0;
	double H34el = 0.0;
	double H35el = 0.0;
	double H41el = 0.0;
	double H42el = 0.0;
	double H43el = 0.0;
	double H44el = 0.0;
	double H45el = 0.0;
	double H51el = 0.0;
	double H52el = 0.0;
	double H53el = 0.0;
	double H54el = 0.0;
	double H55el = 0.0;

	Vector H12(nYield), H13(nYield), H14(nYield), H15(nYield), H21(nYield), H23(nYield), H24(nYield), H25(nYield),
		H31(nYield), H32(nYield), H34(nYield), H35(nYield), H41(nYield), H42(nYield), H43(nYield), H45(nYield),
		H51(nYield), H52(nYield), H53(nYield), H54(nYield);
	H12.Zero();
	H13.Zero();
	H14.Zero();
	H15.Zero();
	H21.Zero();
	H23.Zero();
	H24.Zero();
	H25.Zero();
	H31.Zero();
	H32.Zero();
	H34.Zero();
	H35.Zero();
	H41.Zero();
	H42.Zero();
	H43.Zero();
	H45.Zero();
	H51.Zero();
	H52.Zero();
	H53.Zero();
	H54.Zero();

	Vector Fd_1M(nYield), Fd_2M(nYield), Fd_3M(nYield);
	Fd_1M.Zero();
	Fd_2M.Zero();
	Fd_3M.Zero();

	Vector N1y_n(nYield), N2y_n(nYield), N3y_n(nYield), N4y_n(nYield), N5y_n(nYield), Fd_1(nYield), Fd_2(nYield), Fd_3(nYield), Fd_4(nYield), Fd_5(nYield), G11_1(nYield),
		G11_2(nYield), G11_3(nYield), G11_4(nYield), G11_5(nYield), G12_1(nYield), G12_2(nYield), G12_3(nYield), G12_4(nYield), G12_5(nYield), G13_1(nYield), G13_2(nYield),
		G13_3(nYield), G13_4(nYield), G13_5(nYield), G14_1(nYield), G14_2(nYield), G14_3(nYield), G14_4(nYield), G14_5(nYield), G15_1(nYield), G15_2(nYield), G15_3(nYield), G15_4(nYield),
		G15_5(nYield), G16_1(nYield), G16_2(nYield), G16_3(nYield), G16_4(nYield), G16_5(nYield), G21_1(nYield), G21_2(nYield), G21_3(nYield), G21_4(nYield), G21_5(nYield), G22_1(nYield),
		G22_2(nYield), G22_3(nYield), G22_4(nYield), G22_5(nYield), G23_1(nYield), G23_2(nYield), G23_3(nYield), G23_4(nYield), G23_5(nYield), G24_1(nYield), G24_2(nYield),
		G24_3(nYield), G24_4(nYield), G24_5(nYield), G25_1(nYield), G25_2(nYield), G25_3(nYield), G25_4(nYield), G25_5(nYield), G26_1(nYield), G26_2(nYield), G26_3(nYield),
		G26_4(nYield), G26_5(nYield), G31_1(nYield), G31_2(nYield), G31_3(nYield), G31_4(nYield), G31_5(nYield), G32_1(nYield), G32_2(nYield), G32_3(nYield), G32_4(nYield),
		G32_5(nYield), G33_1(nYield), G33_2(nYield), G33_3(nYield), G33_4(nYield), G33_5(nYield), G34_1(nYield), G34_2(nYield), G34_3(nYield), G34_4(nYield), G34_5(nYield),
		G35_1(nYield), G35_2(nYield), G35_3(nYield), G35_4(nYield), G35_5(nYield), G36_1(nYield), G36_2(nYield), G36_3(nYield), G36_4(nYield), G36_5(nYield), G41_1(nYield),
		G41_2(nYield), G41_3(nYield), G41_4(nYield), G41_5(nYield), G42_1(nYield), G42_2(nYield), G42_3(nYield), G42_4(nYield), G42_5(nYield), G43_1(nYield), G43_2(nYield),
		G43_3(nYield), G43_4(nYield), G43_5(nYield), G44_1(nYield), G44_2(nYield), G44_3(nYield), G44_4(nYield), G44_5(nYield), G45_1(nYield), G45_2(nYield), G45_3(nYield),
		G45_4(nYield), G45_5(nYield), G46_1(nYield), G46_2(nYield), G46_3(nYield), G46_4(nYield), G46_5(nYield), G51_1(nYield), G51_2(nYield), G51_3(nYield), G51_4(nYield),
		G51_5(nYield), G52_1(nYield), G52_2(nYield), G52_3(nYield), G52_4(nYield), G52_5(nYield), G53_1(nYield), G53_2(nYield), G53_3(nYield), G53_4(nYield), G53_5(nYield),
		G54_1(nYield), G54_2(nYield), G54_3(nYield), G54_4(nYield), G54_5(nYield), G55_1(nYield), G55_2(nYield), G55_3(nYield), G55_4(nYield), G55_5(nYield), G56_1(nYield),
		G56_2(nYield), G56_3(nYield), G56_4(nYield), G56_5(nYield), G11(nYield), G12(nYield), G13(nYield), G14(nYield), G15(nYield), G21(nYield), G22(nYield), G23(nYield),
		G24(nYield), G25(nYield), G31(nYield), G32(nYield), G33(nYield), G34(nYield), G35(nYield), G41(nYield), G42(nYield), G43(nYield), G44(nYield), G45(nYield), G51(nYield),
		G52(nYield), G53(nYield), G54(nYield), G55(nYield);

	double Q1 = 0.0;
	double Q2 = 0.0;
	double Q3 = 0.0;
	double Q4 = 0.0;
	double Q5 = 0.0;
	double Q6 = 0.0;

	double pi = 3.141592653589793;

	double Q01 = 0.0;
	double Q02 = 0.0;
	double Q03 = 0.0;
	double Q04 = 0.0;
	double Q05 = 0.0;
	double Q06 = 0.0;

	double q01 = 0.0;
	double q02 = 0.0;
	double q03 = 0.0;
	double q04 = 0.0;
	double q05 = 0.0;
	double q06 = 0.0;

	double Q = 0.0;

	Matrix Ct(6, 6);
	Ct.Zero();
	Matrix Kt(6, 6);
	Kt.Zero();

	double q1 = 0.0;
	double q2 = 0.0;
	double q3 = 0.0;
	double q4 = 0.0;
	double q5 = 0.0;
	double q6 = 0.0;

	double t1 = 0.0;
	double t2 = 0.0;
	double t3 = 0.0;
	double t4 = 0.0;
	double t5 = 0.0;
	double t6 = 0.0;

	q1 = Tstrain(0);
	q2 = Tstrain(1);
	q3 = Tstrain(2);
	q4 = Tstrain(3);
	q5 = Tstrain(4);
	q6 = Tstrain(5);

	if (Peps1 != q1 || Peps2 != q2 || Peps3 != q3) {

		// sign functions --------------------------------------
		double signq1 = 0.0;
		double signq2 = 0.0;
		double signq3 = 0.0;
		double signF2 = 0.0;
		double signPeps1 = 0.0;
		double signPeps2 = 0.0;
		double signPeps3 = 0.0;
		// --------------------------------------------------
		if (abs(q1) > 0) {
			signq1 = q1 / abs(q1);
		}
		if (abs(q1) == 0) {
			signq1 = 0.0;
		}
		if (abs(q2) > 0) {
			signq2 = q2 / abs(q2);
		}
		if (abs(q2) == 0) {
			signq2 = 0.0;
		}
		if (abs(q3) > 0) {
			signq3 = q3 / abs(q3);
		}
		if (abs(q3) == 0) {
			signq3 = 0.0;
		}
		if (abs(Peps1) > 0) {
			signPeps1 = Peps1 / abs(Peps1);
		}
		if (abs(Peps1) == 0) {
			signPeps1 = 0.0;
		}
		if (abs(Peps2) > 0) {
			signPeps2 = Peps2 / abs(Peps2);
		}
		if (abs(Peps2) == 0) {
			signPeps2 = 0.0;
		}
		if (abs(Peps3) > 0) {
			signPeps3 = Peps3 / abs(Peps3);
		}
		if (abs(Peps3) == 0) {
			signPeps3 = 0.0;
		}
		// -----------------------------------------------------

		// Trial forces
		t1 = H11el * q1;
		t2 = H22el * q2;
		t3 = H33el * q3;
		t4 = 0.0;
		t5 = 0.0;
		t6 = 0.0;

		double Qprev = pow((pow(t1, 2.0) + pow(t2, 2.0) + pow(t3, 2.0)), 0.5);
		double QprevH= pow((pow(t1, 2.0) + pow(t2, 2.0)), 0.5);

		int cont = 0;
		int convergence = 0;

		double t1_iter = t1;
		double t2_iter = t2;
		double t3_iter = t3;
		double t4_iter = t4;
		double t5_iter = t5;
		double t6_iter = t6;

		if (abs(t1_iter) < 0.0001) {
			t1_iter = 0.0;
		}
		if (abs(t2_iter) < 0.0001) {
			t2_iter = 0.0;
		}
		if (abs(t3_iter) < 0.0001) {
			t3_iter = 0.0;
		}

		double RatiotH = 0.0;
		double RatiotZ = 0.0;

		double checkRH = 0.0;
		double checkRZ = 0.0;

		double Rt3h = 0.0;
		double Rt21 = 0.0;

		// internal iterative procedure for the elastic-plastic response
		while ((cont <= niter) & (convergence == 0)) {

			cont = cont + 1;

			t1_iter = Cstress(0) +
				PKt(0, 0) * (q1 - Cstrain(0)) + PKt(0, 1) * (q2 - Cstrain(1)) + PKt(0, 2) * (q3 - Cstrain(2));
			t2_iter = Cstress(1) +
				PKt(1, 0) * (q1 - Cstrain(0)) + PKt(1, 1) * (q2 - Cstrain(1)) + PKt(1, 2) * (q3 - Cstrain(2));
			t3_iter = Cstress(2) +
				PKt(2, 0) * (q1 - Cstrain(0)) + PKt(2, 1) * (q2 - Cstrain(1)) + PKt(2, 2) * (q3 - Cstrain(2));

			double tnorm_iter = pow(pow(t1_iter, 2.0) + pow(t2_iter, 2.0) + pow(t3_iter, 2.0), 0.5);

			double signt3_iter = 0.0;
			double signt2_iter = 0.0;
			double signt1_iter = 0.0;

			if (abs(t3_iter) > 0) {
				signt3_iter = t3_iter / abs(t3_iter);
			}
			if (abs(t3_iter) == 0) {
				signt3_iter = 0.0;
			}

			if (abs(t2_iter) > 0) {
				signt2_iter = t2_iter / abs(t2_iter);
			}
			if (abs(t2_iter) == 0) {
				signt2_iter = 0.0;
			}

			if (abs(t1_iter) > 0) {
				signt1_iter = t1_iter / abs(t1_iter);
			}
			if (abs(t1_iter) == 0) {
				signt1_iter = 0.0;
			}

			double tHnorm_iter = pow((pow(t1_iter, 2) + pow(t2_iter, 2)), 0.5);
			double tH3norm_iter = pow(pow(t1, 2.0) + pow(t2, 2.0) + pow(t3, 2.0), 0.5);

			double t31N = pow(pow(t1_iter, 2.0) + pow(t3_iter, 2.0), 0.5);
			double t3hN = pow(pow(t1_iter, 2.0) + pow(t2_iter, 2.0) + pow(t3_iter, 2.0), 0.5);

			RatiotH = t2_iter / t1_iter;
			RatiotZ = tHnorm_iter / t3_iter;

			if ((abs(t1_iter) < 0.0001) && (abs(t2_iter) != 0.0)) {
				RatiotH = 1000000.0;
			}
			if ((abs(t1_iter) < 0.0001) && (abs(t2_iter) < 0.0001)) {
				RatiotH = 0.0;
			}

			if ((abs(t3_iter) < 0.0001) && (tHnorm_iter != 0.0)) {
				RatiotZ = 1000000.0;
			}
			if ((abs(t3_iter) < 0.0001) && (tHnorm_iter < 0.0001)) {
				RatiotZ = 0.0;
			}
			// ---------------------------------------------------------------------


			Rt3h = abs(t3_iter) / pow(pow(t1_iter, 2.), 0.5);
			Rt21 = abs(t2_iter) / abs(t1_iter);

			if (pow(pow(t1_iter, 2), 0.5) < 0.0001) {
				Rt3h = 1.0e18;
			}
			if ((abs(t1_iter) < 0.0001) && (abs(t2_iter) > 0.0001)) {
				Rt21 = 1.0e18;
				t1_iter = 0.;
			}
			if ((abs(t1_iter) < 0.0001) && (abs(t2_iter) < 0.0001)) {
				Rt21 = 0.0;
			}

			
			double signRt21 = 0.0;

			if (abs(Rt21) > 0) {
				signRt21 = Rt21 / abs(Rt21);
			}
			if (abs(Rt21) == 0) {
				signRt21 = 0.0;
			}

			double signRt3h = 0.0;

			if (abs(Rt3h) > 0) {
				signRt3h = Rt3h / abs(Rt3h);
			}
			if (abs(Rt3h) == 0) {
				signRt3h = 0.0;
			}

			
			// tolerances -----------------------------------------
			double tolin = 0.;
			double tolForce = tolin;
			double tolEps = 0.;
			double TolY = 0.;
			double tolRatios = 0.2;


			// sign functions --------------------------------------
			double signCstress3 = 0.0;
			double signCstress2 = 0.0;
			double signCstress1 = 0.0;

			if (abs(t3_iter) > 0) {
				signt3_iter = t3_iter / abs(t3_iter);
			}
			if (abs(t3_iter) == 0) {
				signt3_iter = 0.0;
			}

			if (abs(t2_iter) > 0) {
				signt2_iter = t2_iter / abs(t2_iter);
			}
			if (abs(t2_iter) == 0) {
				signt2_iter = 0.0;
			}

			if (abs(t1_iter) > 0) {
				signt1_iter = t1_iter / abs(t1_iter);
			}
			if (abs(t1_iter) == 0) {
				signt1_iter = 0.0;
			}

			if (abs(Cstress(2)) > 0) {
				signCstress3 = Cstress(2) / abs(Cstress(2));
			}
			if (abs(Cstress(2)) == 0) {
				signCstress3 = 0.0;
			}

			if (abs(Cstress(1)) > 0) {
				signCstress2 = Cstress(1) / abs(Cstress(1));
			}
			if (abs(Cstress(1)) == 0) {
				signCstress2 = 0.0;
			}

			if (abs(Cstress(0)) > 0) {
				signCstress1 = Cstress(0) / abs(Cstress(0));
			}
			if (abs(Cstress(0)) == 0) {
				signCstress1 = 0.0;
			}
			// ---------------------------------------------------------------------


			double Q1_eq = 0.0;
			double Q2_eq = 0.0;
			double Q12_eq = 0.0;
			double Q3_eq13 = 0.0;
			double Q3_eq23 = 0.0;
			double Q13_eq = 0.0;
			double Qh3_eq = 0.0;
			double Q23_eq = 0.0;
			double RQ31_eq = 0.0;
			double RQ21_eq = 0.0;
			double RQ3h_eq = 0.0;
			double RQ32_eq = 0.0;
			double DQ13_eq = 0.0;
			double DQh3_eq = 0.0;
			double DQ23_eq = 0.0;
			double Q1_eq_2 = 0.0;
			double Q2_eq_2 = 0.0;
			double Q12_eq_2 = 0.0;
			double Q3_eq13_2 = 0.0;
			double Q3_eq23_2 = 0.0;
			double Q13_eq_2 = 0.0;
			double Qh3_eq_2 = 0.0;
			double Q23_eq_2 = 0.0;
			double RQ31_eq_2 = 0.0;
			double RQ21_eq_2 = 0.0;
			double RQ3h_eq_2 = 0.0;
			double RQ32_eq_2 = 0.0;
			double DQ13_eq_2 = 0.0;
			double DQh3_eq_2 = 0.0;
			double DQ23_eq_2 = 0.0;
			
			double Q1_eqP = 0.0;
			double Q2_eqP = 0.0;
			double Q12_eqP = 0.0;
			double Q3_eq13P = 0.0;
			double Q3_eq23P = 0.0;
			double Q13_eqP = 0.0;
			double Qh3_eqP = 0.0;
			double Q23_eqP = 0.0;
			double RQ31_eqP = 0.0;
			double RQ21_eqP = 0.0;
			double RQ3h_eqP = 0.0;
			double RQ32_eqP = 0.0;
			double DQ13_eqP = 0.0;
			double DQh3_eqP = 0.0;
			double DQ23_eqP = 0.0;
			double Q1_eq_2P = 0.0;
			double Q2_eq_2P = 0.0;
			double Q12_eq_2P = 0.0;
			double Q3_eq13_2P = 0.0;
			double Q3_eq23_2P = 0.0;
			double Q13_eq_2P = 0.0;
			double Qh3_eq_2P = 0.0;
			double Q23_eq_2P = 0.0;
			double RQ31_eq_2P = 0.0;
			double RQ21_eq_2P = 0.0;
			double RQ3h_eq_2P = 0.0;
			double RQ32_eq_2P = 0.0;
			double DQ13_eq_2P = 0.0;
			double DQh3_eq_2P = 0.0;
			double DQ23_eq_2P = 0.0;
			int ActiveFlows = 0;

			Vector	Qy3_iter13(nYield), Qy3_iter13_2(nYield), Yield_1(nYield), Yield_2(nYield), Yield_3(nYield), Yield(nYield);

			// identification of the (image) yield points ------------------------------------------------
			for (int j = 0; j < nYield; j++) {

				int flag13 = 0;
				int flag23 = 0;

				int flag13_2 = flag13;
				int flag23_2 = flag23;

				int der = j;

				int trick1 = 0;

				if ((abs(t1_iter) < tolForce) && (abs(t2_iter) < tolForce) && (abs(t3_iter) < tolForce)) {
					if ((abs(PCQy1(j)) > 0.0) || (abs(PCQy2(j)) > 0.0) || (abs(PCQy3(j)) > 0.0))
					{
						Qy1_iter(j) = PCQy1(j);
						Qy2_iter(j) = PCQy2(j);
						Qy3_iter(j) = PCQy3(j);
						Qy1_iter_2(j) = PCQy1_2(j);
						Qy2_iter_2(j) = PCQy2_2(j);
						Qy3_iter_2(j) = PCQy3_2(j);
						trick1 = 1;
					}
					else {
						Qy1_iter(j) = 100000000000.0;
						Qy2_iter(j) = 100000000000.0;
						Qy3_iter(j) = 100000000000.0;
						Qy1_iter_2(j) = -Qy1_iter(j);
						Qy2_iter_2(j) = -Qy2_iter(j);
						Qy3_iter_2(j) = -Qy3_iter(j);
						trick1 = 1;
					}
					if (abs(t1_iter) < tolForce) {
						t1_iter = 0.0;
					}
					if (abs(t2_iter) < tolForce) {
						t2_iter = 0.0;
					}
					if (abs(t3_iter) < tolForce) {
						t3_iter = 0.0;
					}

					if (abs(t3_iter) > 0) {
						signt3_iter = t3_iter / abs(t3_iter);
					}
					if (abs(t3_iter) == 0) {
						signt3_iter = 0.0;
					}

					if (abs(t2_iter) > 0) {
						signt2_iter = t2_iter / abs(t2_iter);
					}
					if (abs(t2_iter) == 0) {
						signt2_iter = 0.0;
					}

					if (abs(t1_iter) > 0) {
						signt1_iter = t1_iter / abs(t1_iter);
					}
					if (abs(t1_iter) == 0) {
						signt1_iter = 0.0;
					}
				}


				if ((abs(q1) < tolEps) && (abs(q2) < tolEps) && (abs(q3) < tolEps)) {
					if ((abs(PCQy1(j)) > 0.0) || (abs(PCQy2(j)) > 0.0) || (abs(PCQy3(j)) > 0.0))
					{
						Qy1_iter(j) = PCQy1(j);
						Qy2_iter(j) = PCQy2(j);
						Qy3_iter(j) = PCQy3(j);
						Qy1_iter_2(j) = PCQy1_2(j);
						Qy2_iter_2(j) = PCQy2_2(j);
						Qy3_iter_2(j) = PCQy3_2(j);
						trick1 = 1;
					}
					else {
						Qy1_iter(j) = 100000000000.0;
						Qy2_iter(j) = 100000000000.0;
						Qy3_iter(j) = 100000000000.0;
						Qy1_iter_2(j) = -Qy1_iter(j);
						Qy2_iter_2(j) = -Qy2_iter(j);
						Qy3_iter_2(j) = -Qy3_iter(j);
						trick1 = 1;
					}
					if (abs(t1_iter) < tolForce) {
						t1_iter = 0.0;
					}
					if (abs(t2_iter) < tolForce) {
						t2_iter = 0.0;
					}
					if (abs(t3_iter) < tolForce) {
						t3_iter = 0.0;
					}

					if (abs(t3_iter) > 0) {
						signt3_iter = t3_iter / abs(t3_iter);
					}
					if (abs(t3_iter) == 0) {
						signt3_iter = 0.0;
					}

					if (abs(t2_iter) > 0) {
						signt2_iter = t2_iter / abs(t2_iter);
					}
					if (abs(t2_iter) == 0) {
						signt2_iter = 0.0;
					}

					if (abs(t1_iter) > 0) {
						signt1_iter = t1_iter / abs(t1_iter);
					}
					if (abs(t1_iter) == 0) {
						signt1_iter = 0.0;
					}
				}

				if (PnumUnload > 0) {
					Qy1_iter(j) = PCQy1(j);
					Qy2_iter(j) = PCQy2(j);
					Qy3_iter(j) = PCQy3(j);
					Qy1_iter_2(j) = PCQy1_2(j);
					Qy2_iter_2(j) = PCQy2_2(j);
					Qy3_iter_2(j) = PCQy3_2(j);
					trick1 = 1;

					if (abs(t3_iter) > 0) {
						signt3_iter = t3_iter / abs(t3_iter);
					}
					if (abs(t3_iter) == 0) {
						signt3_iter = 0.0;
					}

					if (abs(t2_iter) > 0) {
						signt2_iter = t2_iter / abs(t2_iter);
					}
					if (abs(t2_iter) == 0) {
						signt2_iter = 0.0;
					}

					if (abs(t1_iter) > 0) {
						signt1_iter = t1_iter / abs(t1_iter);
					}
					if (abs(t1_iter) == 0) {
						signt1_iter = 0.0;
					}
				}

				if (trick1 == 0)
				{
					
					Yield_1(der) = pow((t1_iter- Tc1(der)) * sin(delta) + (t3_iter - Tc3(der)) * cos(delta), 2.) / pow(a3(der),2.);
					Yield_2(der) = pow(t2_iter - Tc2(der), 2.) / pow(a2(der), 2.);
					Yield_3(der) = pow((t1_iter - Tc1(der)) * cos(delta) - (t3_iter - Tc3(der)) * sin(delta), 2.) / pow(a1(der), 2.);
					Yield(der)	 = Yield_1(der) + Yield_2(der) + Yield_3(der) - 1.;

					double mm1 = Tc1(der);
					double mm3 = Tc3(der);
					double aa1 = a1(der);
					double aa2 = a2(der);
					double aa3 = a3(der);
					double HereTheYield_1 = Yield_1(der);
					double HereTheYield_2 = Yield_2(der);
					double HereTheYield_3 = Yield_3(der);
					double HereTheYield = Yield(der);
					double HereTheYield0 = Yield(0);


					Q1_eq = 0.0;
					Q2_eq = 0.0;
					Q12_eq = 0.0;
					Q3_eq13 = 0.0;
					Q3_eq23 = 0.0;
					Q13_eq = 0.0;
					Qh3_eq = 0.0;
					Q23_eq = 0.0;
					RQ31_eq = 0.0;
					RQ21_eq = .0;
					RQ3h_eq = 0.0;
					RQ32_eq = 0.0;
					DQ13_eq = 0.0;
					DQh3_eq = 0.0;
					DQ23_eq = 0.0;

					Q1_eqP = 0.0;
					Q2_eqP = 0.0;
					Q12_eqP = 0.0;
					Q3_eq13P = 0.0;
					Q3_eq23P = 0.0;
					Q13_eqP = 0.0;
					Qh3_eqP = 0.0;
					Q23_eqP = 0.0;
					RQ31_eqP = 0.0;
					RQ21_eqP = 0.0;
					RQ3h_eqP = 0.0;
					RQ32_eqP = 0.0;
					DQ13_eqP = 0.0;
					DQh3_eqP = 0.0;
					DQ23_eqP = 0.0;

					double t1_iter_L = cos(delta) * (t1_iter - Tc1(der)) - sin(delta) * (t3_iter - Tc3(der));
					double t2_iter_L = t2_iter - Tc2(der);
					double t3_iter_L = sin(delta) * (t1_iter - Tc1(der)) + cos(delta) * (t3_iter - Tc3(der));

					double Rt3h_L = abs(t3_iter_L / t1_iter_L);
					double Rt21_L = abs(t2_iter_L / t1_iter_L);

					if (t1_iter>6900.) {
						double eee=1.;
					}

					if (pow(pow(t1_iter_L, 2), 0.5) < tolForce) {
						Rt3h_L = 1.0e18;
					}
					if ((abs(t1_iter_L) < tolForce) && (abs(t2_iter_L) > tolForce)) {
						Rt21_L = 1.0e18;
						t1_iter_L = 0.;
					}
					if ((abs(t1_iter_L) < tolForce) && (abs(t2_iter_L) < tolForce)) {
						Rt21_L = 0.0;
					}

					double t12_L = abs(atan(Rt21_L));
					if ((t1_iter_L < 0) && (t2_iter_L > 0)) {
						t12_L = t12_L + pi;
					}
					if ((t1_iter_L < 0) && (t2_iter_L < 0)) {
						t12_L = t12_L + pi;
					}
					if ((t1_iter_L > 0) && (t2_iter_L < 0)) {
						t12_L = t12_L + 2. * pi;
					}

					double t13_L = abs(atan(Rt3h_L));
					if ((t1_iter_L < 0) && (t3_iter_L > 0)) {
						t13_L = t13_L + pi;
					}
					if ((t1_iter_L < 0) && (t3_iter_L < 0)) {
						t13_L = t13_L + pi;
					}
					if ((t1_iter_L > 0) && (t3_iter_L < 0)) {
						t13_L = t13_L + 2. * pi;
					}

					t13_L = pi / 2. - t13_L + 2. * pi;

					Q1_eq = a3(der) * cos(t13_L) * sin(delta) + a1(der) * sin(t13_L) * cos(t12_L) * cos(delta) + Tc1(der);
					Q2_eq = a2(der) * sin(t12_L) * sin(t13_L) + Tc2(der);
					Q3_eq13 = a3(der) * cos(t13_L) * cos(delta) - a1(der) * sin(t13_L) * cos(t12_L) * sin(delta) + Tc3(der);

					double tc11 = Tc1(der);
					double tc22 = Tc2(der);
					double tc33 = Tc3(der);
					Qy1_iter(j) = Q1_eq;
					Qy2_iter(j) = Q2_eq;
					Qy3_iter(j) = Q3_eq13;

					if (Yield(der) > TolY) {
						ActiveFlows = der + 1;
					}
				
				}
				

				// sign functions
				if (abs(Qy1_iter(j)) > 0) {
					signQy1(j) = Qy1_iter(j) / abs(Qy1_iter(j));
				}
				if (abs(Qy1_iter(j)) == 0) {
					signQy1(j) = 0.0;
				}
				if (abs(Qy2_iter(j)) > 0) {
					signQy2(j) = Qy2_iter(j) / abs(Qy2_iter(j));
				}
				if (abs(Qy2_iter(j)) == 0) {
					signQy2(j) = 0.0;
				}
				if (abs(Qy3_iter(j)) > 0) {
					signQy3(j) = Qy3_iter(j) / abs(Qy3_iter(j));
				}
				if (abs(Qy3_iter(j)) == 0) {
					signQy3(j) = 0.0;
				}
				

				// norm of the image vectors
				Rn_iter(j) = pow(pow(Qy1_iter(j), 2.0) + pow(Qy2_iter(j), 2.0) + pow(Qy3_iter(j), 2.0), 0.5);
				
				if ((Qy1_iter(j) == 0) && (abs(t1_iter) > 0)) {
					Qy1_iter(j) = PCQy1(j);
					Qy2_iter(j) = PCQy2(j);
					Qy3_iter(j) = PCQy3(j);
				}
				if ((Qy2_iter(j) == 0) && (abs(t2_iter) > 0)) {
					Qy1_iter(j) = PCQy1(j);
					Qy2_iter(j) = PCQy2(j);
					Qy3_iter(j) = PCQy3(j);
				}
				if ((Qy3_iter(j) == 0) && (abs(t3_iter) > 0)) {
					Qy1_iter(j) = PCQy1(j);
					Qy2_iter(j) = PCQy2(j);
					Qy3_iter(j) = PCQy3(j);
				}

				if ((abs(Qy1_iter(j)) > tolForce) && (abs(t1_iter) < tolForce)) {
					Qy1_iter(j) = 0.0;
				}
				if ((abs(Qy2_iter(j)) > tolForce) && (abs(t2_iter) < tolForce)) {
					Qy2_iter(j) = 0.0;
				}
				if ((abs(Qy3_iter(j)) > tolForce) && (abs(t3_iter) < tolForce)) {
					Qy3_iter(j) = 0.0;
				}
				
				if ((abs(t1_iter) < tolForce) && (abs(t2_iter) < tolForce) && (abs(t3_iter) < tolForce)) {
					if ((abs(PCQy1(j)) > 0.0) || (abs(PCQy2(j)) > 0.0) || (abs(PCQy3(j)) > 0.0)) {
						Qy1_iter(j) = PCQy1(j);
						Qy2_iter(j) = PCQy2(j);
						Qy3_iter(j) = PCQy3(j);
					}
					else {
						Qy1_iter(j) = pow(10.,11);
						Qy2_iter(j) = pow(10., 11);
						Qy3_iter(j) = pow(10., 11);
					}
					if (abs(t1_iter) < tolForce) {
						t1_iter = 0.0;
					}
					if (abs(t2_iter) < tolForce) {
						t2_iter = 0.0;
					}
					if (abs(t3_iter) < tolForce) {
						t3_iter = 0.0;
					}
				}
				if ((abs(q1) < tolEps) && (abs(q2) < tolEps) && (abs(q3) < tolEps)) {
					if ((abs(PCQy1(j)) > 0.0) || (abs(PCQy2(j)) > 0.0) || (abs(PCQy3(j)) > 0.0)) {
						Qy1_iter(j) = PCQy1(j);
						Qy2_iter(j) = PCQy2(j);
						Qy3_iter(j) = PCQy3(j);
					}
					else {
						Qy1_iter(j) = pow(10., 11);
						Qy2_iter(j) = pow(10., 11);
						Qy3_iter(j) = pow(10., 11);
					}
					if (abs(t1_iter) < tolForce) {
						t1_iter = 0.0;
					}
					if (abs(t2_iter) < tolForce) {
						t2_iter = 0.0;
					}
					if (abs(t3_iter) < tolForce) {
						t3_iter = 0.0;
					}
				}
				if ((abs(Qy1_iter(j)) < tolForce) && (abs(Qy2_iter(j)) < tolForce) && (abs(Qy3_iter(j)) < tolForce)) {
					if ((abs(PCQy1(j)) > 0.0) || (abs(PCQy2(j)) > 0.0) || (abs(PCQy3(j)) > 0.0)) {
						Qy1_iter(j) = PCQy1(j);
						Qy2_iter(j) = PCQy2(j);
						Qy3_iter(j) = PCQy3(j);
					}
					else {
						Qy1_iter(j) = 100000000000.0;
						Qy2_iter(j) = 100000000000.0;
						Qy3_iter(j) = 100000000000.0;
					}
				}

				if ((abs(Qy1_iter(j)) == 1000000000000000.0)) {
					if ((abs(PCQy1(j)) > 0.0) || (abs(PCQy1_2(j)) > 0.0)) {
						Qy1_iter(j) = PCQy1(j);
						Qy1_iter_2(j) = PCQy1_2(j);
					}
				}
				if ((abs(Qy2_iter(j)) == 1000000000000000.0)) {
					if ((abs(PCQy2(j)) > 0.0)) {
						Qy2_iter(j) = PCQy2(j);
					}
				}
				if (abs(Qy3_iter(j)) == 1000000000000000.0) {
					if ((abs(PCQy3(j)) > 0.0)) {
						Qy3_iter(j) = PCQy3(j);
					}
				}

				Rn_iter(j) = pow(pow(Qy1_iter(j), 2.0) + pow(Qy2_iter(j), 2.0) + pow(Qy3_iter(j), 2.0), 0.5);

				double Qy1_iterCd = Qy1_iter(j);
				double Qy2_iterCd = Qy2_iter(j);
				double Qy3_iterCd = Qy3_iter(j);

				// dissipative forces
				csi1_iter(j) = Qy1_iter(j) - Tc1(j);
				csi2_iter(j) = Qy2_iter(j) - Tc2(j);
				csi3_iter(j) = Qy3_iter(j) - Tc3(j);

			} // loop on the yield surfaces finished

			double PStrH = pow(pow(PStr1, 2.) + pow(PStr2, 2.), .5);
			double RatiotHP = PStr2 / PStr1;
			double RatiotZP = PStrH / PStr3;

			if ((abs(PStr1) < 0.0001) && (abs(PStr2) != 0.0)) {
				RatiotHP = 1000000.0;
			}
			if ((abs(PStr1) < 0.0001) && (abs(PStr2) < 0.0001)) {
				RatiotHP = 0.0;
			}

			if ((abs(PStr3) < 0.0001) && (PStrH != 0.0)) {
				RatiotZP = 1000000.0;
			}
			if ((abs(PStr3) < 0.0001) && (PStrH < 0.0001)) {
				RatiotZP = 0.0;
			}

			double IncrQ1 = t1_iter - PStr1;
			double IncrQ2 = t2_iter - PStr2;
			double IncrQ3 = t3_iter - PStr3;

			double IncrQ1_2 = PStr1 - PStr1_2;
			double IncrQ2_2 = PStr2 - PStr2_2;
			double IncrQ3_2 = PStr3 - PStr3_2;

			double signIncrQ1 = 0.0;
			double signIncrQ2 = 0.0;
			double signIncrQ3 = 0.0;
			double signIncrQ1_2 = 0.0;
			double signIncrQ2_2 = 0.0;
			double signIncrQ3_2 = 0.0;

			if (abs(IncrQ1) > 0) {
				signIncrQ1 = IncrQ1 / abs(IncrQ1);
			}
			if (abs(IncrQ2) > 0) {
				signIncrQ2 = IncrQ2 / abs(IncrQ2);
			}
			if (abs(IncrQ3) > 0) {
				signIncrQ3 = IncrQ3 / abs(IncrQ3);
			}

			if (abs(IncrQ1_2) > 0) {
				signIncrQ1_2 = IncrQ1_2 / abs(IncrQ1_2);
			}
			if (abs(IncrQ2_2) > 0) {
				signIncrQ2_2 = IncrQ2_2 / abs(IncrQ2_2);
			}
			if (abs(IncrQ3_2) > 0) {
				signIncrQ3_2 = IncrQ3_2 / abs(IncrQ3_2);
			}

			// enforce initial stiffness upon load reversal
			if ((pow(pow(PStr1, 2.) + pow(PStr2, 2.) + pow(PStr3, 2.), 0.5) >
				pow(pow(t1_iter, 2.) + pow(t2_iter, 2.) + pow(t3_iter, 2.), 0.5)) &&
				abs(RatiotHP-RatiotH) < tolRatios && abs(RatiotZP-RatiotZ) < tolRatios
				&& PrnFlows == ActiveFlows && ActiveFlows > 0) {
				if ((signIncrQ1 != signIncrQ1_2) || (signIncrQ2 != signIncrQ2_2) || (signIncrQ3 != signIncrQ3_2)) {
					ActiveFlows = 0;
				}
			}
			// enforce ultimate capacity
			if ((pow(pow(PStr1, 2.) + pow(PStr2, 2.) + pow(PStr3, 2.), 0.5) <=
				pow(pow(t1_iter, 2.) + pow(t2_iter, 2.) + pow(t3_iter, 2.), 0.5)) &&
				abs(RatiotHP-RatiotH) < tolRatios && abs(RatiotZP-RatiotZ) < tolRatios && PrnFlows == ActiveFlows-1 && ActiveFlows == nYield) {
				ActiveFlows = nYield;
			}
			// avoid jumps between too many surfaces
			if (ActiveFlows > PrnFlows+2) {
				ActiveFlows = PrnFlows;
			}
			if (PrnFlows == 0 && ActiveFlows > 1) {
				ActiveFlows = PrnFlows;
			}


			nFlows(cont) = ActiveFlows;

			double Qy1_iterC0 = Qy1_iter(0);
			double Qy2_iterC0 = Qy2_iter(0);
			double Qy3_iterC0 = Qy3_iter(0);

			double Qy1_iterCN = Qy1_iter(nYield - 1);
			double Qy2_iterCN = Qy2_iter(nYield - 1);
			double Qy3_iterCN = Qy3_iter(nYield - 1);

			for (int j = 0; j < nYield; j++) {
				if (abs(Qy1_iter(j)) > 1000000000.0) {
					Qy1_iter(j) = PCQy1(j);
				}
				if (abs(Qy2_iter(j)) > 1000000000.0) {
					Qy2_iter(j) = PCQy2(j);
				}
				if (abs(Qy3_iter(j)) > 1000000000.0) {
					Qy3_iter(j) = PCQy3(j);
				}
			}

			
			int numpl = 0;

			if ((cont > 1) && (nFlows(cont - 1) == nYield) && (nFlows(cont) == nYield) && (PrnFlows < nYield - 2)) {
				nFlows(cont) = PrnFlows;
			}
			if ((cont > 2) && (nFlows(cont) == nYield) && (PrnFlows < nYield - 2)) {
				nFlows(cont) = PrnFlows;
			}

			if ((unload == 0) && (PnumUnload == 0) && (abs(q3 / pow(pow(q1, 2.) + pow(q2, 2.), .5) - Peps3 / pow(pow(Peps1, 2.) + pow(Peps2, 2.), .5)) < tolin) && (nFlows(cont) == 0) && (PrnFlows > 0)) {
				nFlows(cont) = PrnFlows;
			}
			if ((unload == 0) && (PnumUnload == 0) && (abs(q3) > tolEps) && (pow(pow(q1, 2.) + pow(q2, 2.), .5) < tolEps) && (nFlows(cont) == 0) && (PrnFlows > 0)) {
				nFlows(cont) = PrnFlows;
			}


			if ((nFlows(cont) == nYield) && (abs(t1_iter) < tolEps) && (abs(q1) > tolEps)) {
				nFlows(cont) = PrnFlows;
			}
			if ((nFlows(cont) == nYield) && (abs(t2_iter) < tolEps) && (abs(q2) > tolEps)) {
				nFlows(cont) = PrnFlows;
			}
			if ((nFlows(cont) == nYield) && (abs(t3_iter) < tolEps) && (abs(q3) > tolEps)) {
				nFlows(cont) = PrnFlows;
			}
			


			// elastic displacements of the iteration procedure
			double q1el_iter = 1 / H11(0) * t1_iter;
			double q2el_iter = 1 / H22(0) * t2_iter;
			double q3el_iter = 1 / H33(0) * t3_iter;
			
			// elastic-plastic stiffness matrix
			N1y_n.Zero();
			N2y_n.Zero();
			N3y_n.Zero();
			N4y_n.Zero();
			N5y_n.Zero();

			Fd_1.Zero();
			Fd_2.Zero();
			Fd_3.Zero();
			Fd_4.Zero();
			Fd_5.Zero();

			G11_1.Zero();
			G11_2.Zero();
			G11_3.Zero();
			G11_4.Zero();
			G11_5.Zero();

			G12_1.Zero();
			G12_2.Zero();
			G12_3.Zero();
			G12_4.Zero();
			G12_5.Zero();

			G13_1.Zero();
			G13_2.Zero();
			G13_3.Zero();
			G13_4.Zero();
			G13_5.Zero();

			G14_1.Zero();
			G14_2.Zero();
			G14_3.Zero();
			G14_4.Zero();
			G14_5.Zero();

			G15_1.Zero();
			G15_2.Zero();
			G15_3.Zero();
			G15_4.Zero();
			G15_5.Zero();

			G16_1.Zero();
			G16_2.Zero();
			G16_3.Zero();
			G16_4.Zero();
			G16_5.Zero();

			G21_1.Zero();
			G21_2.Zero();
			G21_3.Zero();
			G21_4.Zero();
			G21_5.Zero();

			G22_1.Zero();
			G22_2.Zero();
			G22_3.Zero();
			G22_4.Zero();
			G22_5.Zero();

			G23_1.Zero();
			G23_2.Zero();
			G23_3.Zero();
			G23_4.Zero();
			G23_5.Zero();

			G24_1.Zero();
			G24_2.Zero();
			G24_3.Zero();
			G24_4.Zero();
			G24_5.Zero();

			G25_1.Zero();
			G25_2.Zero();
			G25_3.Zero();
			G25_4.Zero();
			G25_5.Zero();

			G26_1.Zero();
			G26_2.Zero();
			G26_3.Zero();
			G26_4.Zero();
			G26_5.Zero();

			G31_1.Zero();
			G31_2.Zero();
			G31_3.Zero();
			G31_4.Zero();
			G31_5.Zero();

			G32_1.Zero();
			G32_2.Zero();
			G32_3.Zero();
			G32_4.Zero();
			G32_5.Zero();

			G33_1.Zero();
			G33_2.Zero();
			G33_3.Zero();
			G33_4.Zero();
			G33_5.Zero();

			G34_1.Zero();
			G34_2.Zero();
			G34_3.Zero();
			G34_4.Zero();
			G34_5.Zero();

			G35_1.Zero();
			G35_2.Zero();
			G35_3.Zero();
			G35_4.Zero();
			G35_5.Zero();

			G36_1.Zero();
			G36_2.Zero();
			G36_3.Zero();
			G36_4.Zero();
			G36_5.Zero();

			G41_1.Zero();
			G41_2.Zero();
			G41_3.Zero();
			G41_4.Zero();
			G41_5.Zero();

			G42_1.Zero();
			G42_2.Zero();
			G42_3.Zero();
			G42_4.Zero();
			G42_5.Zero();

			G43_1.Zero();
			G43_2.Zero();
			G43_3.Zero();
			G43_4.Zero();
			G43_5.Zero();

			G44_1.Zero();
			G44_2.Zero();
			G44_3.Zero();
			G44_4.Zero();
			G44_5.Zero();

			G45_1.Zero();
			G45_2.Zero();
			G45_3.Zero();
			G45_4.Zero();
			G45_5.Zero();

			G46_1.Zero();
			G46_2.Zero();
			G46_3.Zero();
			G46_4.Zero();
			G46_5.Zero();

			G51_1.Zero();
			G51_2.Zero();
			G51_3.Zero();
			G51_4.Zero();
			G51_5.Zero();

			G52_1.Zero();
			G52_2.Zero();
			G52_3.Zero();
			G52_4.Zero();
			G52_5.Zero();

			G53_1.Zero();
			G53_2.Zero();
			G53_3.Zero();
			G53_4.Zero();
			G53_5.Zero();

			G54_1.Zero();
			G54_2.Zero();
			G54_3.Zero();
			G54_4.Zero();
			G54_5.Zero();

			G55_1.Zero();
			G55_2.Zero();
			G55_3.Zero();
			G55_4.Zero();
			G55_5.Zero();

			G56_1.Zero();
			G56_2.Zero();
			G56_3.Zero();
			G56_4.Zero();
			G56_5.Zero();


			G11.Zero();
			G12.Zero();
			G13.Zero();
			G14.Zero();
			G15.Zero();

			G21.Zero();
			G22.Zero();
			G23.Zero();
			G24.Zero();
			G25.Zero();

			G31.Zero();
			G32.Zero();
			G33.Zero();
			G34.Zero();
			G35.Zero();

			G41.Zero();
			G42.Zero();
			G43.Zero();
			G44.Zero();
			G45.Zero();

			G51.Zero();
			G52.Zero();
			G53.Zero();
			G54.Zero();
			G55.Zero();


			double G11_f = 0.0;
			double G12_f = 0.0;
			double G13_f = 0.0;
			double G14_f = 0.0;
			double G15_f = 0.0;

			double G21_f = 0.0;
			double G22_f = 0.0;
			double G23_f = 0.0;
			double G24_f = 0.0;
			double G25_f = 0.0;

			double G31_f = 0.0;
			double G32_f = 0.0;
			double G33_f = 0.0;
			double G34_f = 0.0;
			double G35_f = 0.0;

			double G41_f = 0.0;
			double G42_f = 0.0;
			double G43_f = 0.0;
			double G44_f = 0.0;
			double G45_f = 0.0;

			double G51_f = 0.0;
			double G52_f = 0.0;
			double G53_f = 0.0;
			double G54_f = 0.0;
			double G55_f = 0.0;
			double Nas = 0.;
			double AssFl = 0.;

			// elastic-plastic stiffness
			if (nFlows(cont) > 0.0) {
				for (int j = 0; j < nFlows(cont); j++) {

					N1y_n(j) = 2.0 * cos(delta) / pow(a1(j), 2.0) *
						(cos(delta) * (csi1_iter(j) - Tc1(j)) -
							sin(delta) * (csi3_iter(j) - Tc3(j))) +
						2.0 * sin(delta) / pow(a3(j), 2.0) *
						(sin(delta) * (csi1_iter(j) - Tc1(j)) +
							cos(delta) * (csi3_iter(j) - Tc3(j)));

					N2y_n(j) = 2.0 / pow(a2(j), 2.0) * (csi2_iter(j) - Tc2(j));

					N3y_n(j) = 2.0 * sin(delta) / pow(a1(j), 2.0) *
						(cos(delta) * (csi1_iter(j) - Tc1(j)) -
							sin(delta) * (csi3_iter(j) - Tc3(j))) +
						2.0 * cos(delta) / pow(a3(j), 2.0) *
						(sin(delta) * (csi1_iter(j) - Tc1(j)) +
							cos(delta) * (csi3_iter(j) - Tc3(j)));

					N4y_n(j) = 0.0;

					N5y_n(j) = 0.0;

					Fd_1(j) = abs(pow(N1y_n(j), 2.0) * H11(j) + pow(N2y_n(j), 2.0) * H12(j) +
						pow(N3y_n(j), 2.0) * H13(j));
					Fd_2(j) = abs(pow(N1y_n(j), 2.0) * H21(j) + pow(N2y_n(j), 2.0) * H22(j) +
						pow(N3y_n(j), 2.0) * H23(j));
					Fd_3(j) = abs(pow(N1y_n(j), 2.0) * H31(j) + pow(N2y_n(j), 2.0) * H32(j) +
						pow(N3y_n(j), 2.0) * H33(j));
					Fd_4(j) = 10000000000000.0;
					Fd_5(j) = 10000000000000.0;

					Fd_1M(j) = abs(pow(N1y_n(j), 2.0) * H11(j) + pow(N2y_n(j), 2.0) * H12(j) +
						pow(N3y_n(j), 2.0) * H13(j) + pow(N4y_n(j), 2.0) * H14(j) +
						pow(N5y_n(j), 2.0) * H15(j));
					Fd_2M(j) = abs(pow(N1y_n(j), 2.0) * H21(j) + pow(N2y_n(j), 2.0) * H22(j) +
						pow(N3y_n(j), 2.0) * H23(j) + pow(N4y_n(j), 2.0) * H24(j) +
						pow(N5y_n(j), 2.0) * H25(j));
					Fd_3M(j) = abs(pow(N1y_n(j), 2.0) * H31(j) + pow(N2y_n(j), 2.0) * H32(j) +
						pow(N3y_n(j), 2.0) * H33(j) + pow(N4y_n(j), 2.0) * H34(j) +
						pow(N5y_n(j), 2.0) * H35(j));


					if ((abs(H11el * N1y_n(j) * N1y_n(j)) < 0.0000001) && (abs(Fd_1(j)) < 0.0000001)) {
						G11_1(j) = 0.0;
					}
					else {
						G11_1(j) = H11el * N1y_n(j) * N1y_n(j) / Fd_1(j);
					}
					if ((abs(H12el * N2y_n(j) * N1y_n(j)) < 0.0000001) && (abs(Fd_2(j)) < 0.0000001)) {
						G11_2(j) = 0.0;
					}
					else {
						G11_2(j) = H12el * N2y_n(j) * N1y_n(j) / Fd_2(j);
					}
					if ((abs(H13el * N3y_n(j) * N1y_n(j)) < 0.0000001) && (abs(Fd_3(j)) < 0.0000001)) {
						G11_3(j) = 0.0;
					}
					else {
						G11_3(j) = H13el * N3y_n(j) * N1y_n(j) / Fd_3(j);
					}
					if ((abs(H14el * N4y_n(j) * N1y_n(j)) < 0.0000001) && (abs(Fd_4(j)) < 0.0000001)) {
						G11_4(j) = 0.0;
					}
					else {
						G11_4(j) = H14el * N4y_n(j) * N1y_n(j) / Fd_4(j);
					}
					if ((abs(H15el * N5y_n(j) * N1y_n(j)) < 0.0000001) && (abs(Fd_5(j)) < 0.0000001)) {
						G11_5(j) = 0.0;
					}
					else {
						G11_5(j) = H15el * N5y_n(j) * N1y_n(j) / Fd_5(j);
					}
					G11(j) = G11_1(j) + G11_2(j) + G11_3(j) + G11_4(j) +
						G11_5(j);
					G11_f = G11_f + G11(j);
					// -------------------------------------------------------------

					if ((abs(H21el * N1y_n(j) * N1y_n(j)) < 0.0000001) && (abs(Fd_1(j)) < 0.0000001)) {
						G12_1(j) = 0.0;
					}
					else {
						G12_1(j) = H21el * N1y_n(j) * N1y_n(j) / Fd_1(j);
					}
					if ((abs(H22el * N2y_n(j) * N1y_n(j)) < 0.0000001) && (abs(Fd_2(j)) < 0.0000001)) {
						G12_2(j) = 0.0;
					}
					else {
						G12_2(j) = H22el * N2y_n(j) * N1y_n(j) / Fd_2(j);
					}
					if ((abs(H23el * N3y_n(j) * N1y_n(j)) < 0.0000001) && (abs(Fd_3(j)) < 0.0000001)) {
						G12_3(j) = 0.0;
					}
					else {
						G12_3(j) = H23el * N3y_n(j) * N1y_n(j) / Fd_3(j);
					}
					if ((abs(H24el * N4y_n(j) * N1y_n(j)) < 0.0000001) && (abs(Fd_4(j)) < 0.0000001)) {
						G12_4(j) = 0.0;
					}
					else {
						G12_4(j) = H24el * N4y_n(j) * N1y_n(j) / Fd_4(j);
					}
					if ((abs(H25el * N5y_n(j) * N1y_n(j)) < 0.0000001) && (abs(Fd_5(j)) < 0.0000001)) {
						G12_5(j) = 0.0;
					}
					else {
						G12_5(j) = H25el * N5y_n(j) * N1y_n(j) / Fd_5(j);
					}
					G12(j) = G12_1(j) + G12_2(j) + G12_3(j) + G12_4(j) +
						G12_5(j);
					if (j >= int(nYield*AssFl)) {
						G12(j) = 0.;
					}
					G12_f = G12_f + G12(j);
					// -------------------------------------------------------------

					if ((abs(H31el * N1y_n(j) * N1y_n(j)) < 0.0000001) && (abs(Fd_1(j)) < 0.0000001)) {
						G13_1(j) = 0.0;
					}
					else {
						G13_1(j) = H31el * N1y_n(j) * N1y_n(j) / Fd_1(j);
					}
					if ((abs(H32el * N2y_n(j) * N1y_n(j)) < 0.0000001) && (abs(Fd_2(j)) < 0.0000001)) {
						G13_2(j) = 0.0;
					}
					else {
						G13_2(j) = H32el * N2y_n(j) * N1y_n(j) / Fd_2(j);
					}
					if ((abs(H33el * N3y_n(j) * N1y_n(j)) < 0.0000001) && (abs(Fd_3(j)) < 0.0000001)) {
						G13_3(j) = 0.0;
					}
					else {
						G13_3(j) = H33el * N3y_n(j) * N1y_n(j) / Fd_3(j);
					}
					if ((abs(H34el * N4y_n(j) * N1y_n(j)) < 0.0000001) && (abs(Fd_4(j)) < 0.0000001)) {
						G13_4(j) = 0.0;
					}
					else {
						G13_4(j) = H34el * N4y_n(j) * N1y_n(j) / Fd_4(j);
					}
					if ((abs(H35el * N5y_n(j) * N1y_n(j)) < 0.0000001) && (abs(Fd_5(j)) < 0.0000001)) {
						G13_5(j) = 0.0;
					}
					else {
						G13_5(j) = H35el * N5y_n(j) * N1y_n(j) / Fd_5(j);
					}
					G13(j) = G13_1(j) + G13_2(j) + G13_3(j) + G13_4(j) +
						G13_5(j);
					if (j >= int(nYield*AssFl)) {
						G13(j) = 0.;
					}
					G13_f = G13_f + G13(j);
					// -------------------------------------------------------------

					if ((abs(H41el * N1y_n(j) * N1y_n(j)) < 0.0000001) && (abs(Fd_1(j)) < 0.0000001)) {
						G14_1(j) = 0.0;
					}
					else {
						G14_1(j) = H41el * N1y_n(j) * N1y_n(j) / Fd_1(j);
					}
					if ((abs(H42el * N2y_n(j) * N1y_n(j)) < 0.0000001) && (abs(Fd_2(j)) < 0.0000001)) {
						G14_2(j) = 0.0;
					}
					else {
						G14_2(j) = H42el * N2y_n(j) * N1y_n(j) / Fd_2(j);
					}
					if ((abs(H43el * N3y_n(j) * N1y_n(j)) < 0.0000001) && (abs(Fd_3(j)) < 0.0000001)) {
						G14_3(j) = 0.0;
					}
					else {
						G14_3(j) = H43el * N3y_n(j) * N1y_n(j) / Fd_3(j);
					}
					if ((abs(H44el * N4y_n(j) * N1y_n(j)) < 0.0000001) && (abs(Fd_4(j)) < 0.0000001)) {
						G14_4(j) = 0.0;
					}
					else {
						G14_4(j) = H44el * N4y_n(j) * N1y_n(j) / Fd_4(j);
					}
					if ((abs(H45el * N5y_n(j) * N1y_n(j)) < 0.0000001) && (abs(Fd_5(j)) < 0.0000001)) {
						G14_5(j) = 0.0;
					}
					else {
						G14_5(j) = H45el * N5y_n(j) * N1y_n(j) / Fd_5(j);
					}
					G14(j) = G14_1(j) + G14_2(j) + G14_3(j) + G14_4(j) +
						G14_5(j);
					if (j >= int(nYield*AssFl)) {
						G14(j) = 0.;
					}
					G14_f = G14_f + G14(j);
					// -------------------------------------------------------------

					if ((abs(H51el * N1y_n(j) * N1y_n(j)) < 0.0000001) && (abs(Fd_1(j)) < 0.0000001)) {
						G15_1(j) = 0.0;
					}
					else {
						G15_1(j) = H51el * N1y_n(j) * N1y_n(j) / Fd_1(j);
					}
					if ((abs(H52el * N2y_n(j) * N1y_n(j)) < 0.0000001) && (abs(Fd_2(j)) < 0.0000001)) {
						G15_2(j) = 0.0;
					}
					else {
						G15_2(j) = H52el * N2y_n(j) * N1y_n(j) / Fd_2(j);
					}
					if ((abs(H53el * N3y_n(j) * N1y_n(j)) < 0.0000001) && (abs(Fd_3(j)) < 0.0000001)) {
						G15_3(j) = 0.0;
					}
					else {
						G15_3(j) = H53el * N3y_n(j) * N1y_n(j) / Fd_3(j);
					}
					if ((abs(H54el * N4y_n(j) * N1y_n(j)) < 0.0000001) && (abs(Fd_4(j)) < 0.0000001)) {
						G15_4(j) = 0.0;
					}
					else {
						G15_4(j) = H54el * N4y_n(j) * N1y_n(j) / Fd_4(j);
					}
					if ((abs(H55el * N5y_n(j) * N1y_n(j)) < 0.0000001) && (abs(Fd_5(j)) < 0.0000001)) {
						G15_5(j) = 0.0;
					}
					else {
						G15_5(j) = H55el * N5y_n(j) * N1y_n(j) / Fd_5(j);
					}
					G15(j) = G15_1(j) + G15_2(j) + G15_3(j) + G15_4(j) +
						G15_5(j);
					if (j >= int(nYield*AssFl)) {
						G15(j) = 0.;
					}
					G15_f = G15_f + G15(j);
					// -------------------------------------------------------------
					// -------------------------------------------------------------

					if ((abs(H11el * N1y_n(j) * N2y_n(j)) < 0.0000001) && (abs(Fd_1(j)) < 0.0000001)) {
						G21_1(j) = 0.0;
					}
					else {
						G21_1(j) = H11el * N1y_n(j) * N2y_n(j) / Fd_1(j);
					}
					if ((abs(H12el * N2y_n(j) * N2y_n(j)) < 0.0000001) && (abs(Fd_2(j)) < 0.0000001)) {
						G21_2(j) = 0.0;
					}
					else {
						G21_2(j) = H12el * N2y_n(j) * N2y_n(j) / Fd_2(j);
					}
					if ((abs(H13el * N3y_n(j) * N2y_n(j)) < 0.0000001) && (abs(Fd_3(j)) < 0.0000001)) {
						G21_3(j) = 0.0;
					}
					else {
						G21_3(j) = H13el * N3y_n(j) * N2y_n(j) / Fd_3(j);
					}
					if ((abs(H14el * N4y_n(j) * N2y_n(j)) < 0.0000001) && (abs(Fd_4(j)) < 0.0000001)) {
						G21_4(j) = 0.0;
					}
					else {
						G21_4(j) = H14el * N4y_n(j) * N2y_n(j) / Fd_4(j);
					}
					if ((abs(H15el * N5y_n(j) * N2y_n(j)) < 0.0000001) && (abs(Fd_5(j)) < 0.0000001)) {
						G21_5(j) = 0.0;
					}
					else {
						G21_5(j) = H15el * N5y_n(j) * N2y_n(j) / Fd_5(j);
					}
					G21(j) = G21_1(j) + G21_2(j) + G21_3(j) + G21_4(j) +
						G21_5(j);
					if (j >= int(nYield*AssFl)) {
						G21(j) = 0.;
					}
					G21_f = G21_f + G21(j);
					// -------------------------------------------------------------

					if ((abs(H21el * N1y_n(j) * N2y_n(j)) < 0.0000001) && (abs(Fd_1(j)) < 0.0000001)) {
						G22_1(j) = 0.0;
					}
					else {
						G22_1(j) = H21el * N1y_n(j) * N2y_n(j) / Fd_1(j);
					}
					if ((abs(H22el * N2y_n(j) * N2y_n(j)) < 0.0000001) && (abs(Fd_2(j)) < 0.0000001)) {
						G22_2(j) = 0.0;
					}
					else {
						G22_2(j) = H22el * N2y_n(j) * N2y_n(j) / Fd_2(j);
					}
					if ((abs(H23el * N3y_n(j) * N2y_n(j)) < 0.0000001) && (abs(Fd_3(j)) < 0.0000001)) {
						G22_3(j) = 0.0;
					}
					else {
						G22_3(j) = H23el * N3y_n(j) * N2y_n(j) / Fd_3(j);
					}
					if ((abs(H24el * N4y_n(j) * N2y_n(j)) < 0.0000001) && (abs(Fd_4(j)) < 0.0000001)) {
						G22_4(j) = 0.0;
					}
					else {
						G22_4(j) = H24el * N4y_n(j) * N2y_n(j) / Fd_4(j);
					}
					if ((abs(H25el * N5y_n(j) * N2y_n(j)) < 0.0000001) && (abs(Fd_5(j)) < 0.0000001)) {
						G22_5(j) = 0.0;
					}
					else {
						G22_5(j) = H25el * N5y_n(j) * N2y_n(j) / Fd_5(j);
					}
					G22(j) = G22_1(j) + G22_2(j) + G22_3(j) + G22_4(j) +
						G22_5(j);
					G22_f = G22_f + G22(j);
					// -------------------------------------------------------------

					if ((abs(H31el * N1y_n(j) * N2y_n(j)) < 0.0000001) && (abs(Fd_1(j)) < 0.0000001)) {
						G23_1(j) = 0.0;
					}
					else {
						G23_1(j) = H31el * N1y_n(j) * N2y_n(j) / Fd_1(j);
					}
					if ((abs(H32el * N2y_n(j) * N2y_n(j)) < 0.0000001) && (abs(Fd_2(j)) < 0.0000001)) {
						G23_2(j) = 0.0;
					}
					else {
						G23_2(j) = H32el * N2y_n(j) * N2y_n(j) / Fd_2(j);
					}
					if ((abs(H33el * N3y_n(j) * N2y_n(j)) < 0.0000001) && (abs(Fd_3(j)) < 0.0000001)) {
						G23_3(j) = 0.0;
					}
					else {
						G23_3(j) = H33el * N3y_n(j) * N2y_n(j) / Fd_3(j);
					}
					if ((abs(H34el * N4y_n(j) * N2y_n(j)) < 0.0000001) && (abs(Fd_4(j)) < 0.0000001)) {
						G23_4(j) = 0.0;
					}
					else {
						G23_4(j) = H34el * N4y_n(j) * N2y_n(j) / Fd_4(j);
					}
					if ((abs(H35el * N5y_n(j) * N2y_n(j)) < 0.0000001) && (abs(Fd_5(j)) < 0.0000001)) {
						G23_5(j) = 0.0;
					}
					else {
						G23_5(j) = H35el * N5y_n(j) * N2y_n(j) / Fd_5(j);
					}
					G23(j) = G23_1(j) + G23_2(j) + G23_3(j) + G23_4(j) +
						G23_5(j);
					if (j >= int(nYield*AssFl)) {
						G23(j) = 0.;
					}
					G23_f = G23_f + G23(j);
					// -------------------------------------------------------------

					if ((abs(H41el * N1y_n(j) * N2y_n(j)) < 0.0000001) && (abs(Fd_1(j)) < 0.0000001)) {
						G24_1(j) = 0.0;
					}
					else {
						G24_1(j) = H41el * N1y_n(j) * N2y_n(j) / Fd_1(j);
					}
					if ((abs(H42el * N2y_n(j) * N2y_n(j)) < 0.0000001) && (abs(Fd_2(j)) < 0.0000001)) {
						G24_2(j) = 0.0;
					}
					else {
						G24_2(j) = H42el * N2y_n(j) * N2y_n(j) / Fd_2(j);
					}
					if ((abs(H43el * N3y_n(j) * N2y_n(j)) < 0.0000001) && (abs(Fd_3(j)) < 0.0000001)) {
						G24_3(j) = 0.0;
					}
					else {
						G24_3(j) = H43el * N3y_n(j) * N2y_n(j) / Fd_3(j);
					}
					if ((abs(H44el * N4y_n(j) * N2y_n(j)) < 0.0000001) && (abs(Fd_4(j)) < 0.0000001)) {
						G24_4(j) = 0.0;
					}
					else {
						G24_4(j) = H44el * N4y_n(j) * N2y_n(j) / Fd_4(j);
					}
					if ((abs(H45el * N5y_n(j) * N2y_n(j)) < 0.0000001) && (abs(Fd_5(j)) < 0.0000001)) {
						G24_5(j) = 0.0;
					}
					else {
						G24_5(j) = H45el * N5y_n(j) * N2y_n(j) / Fd_5(j);
					}
					G24(j) = G24_1(j) + G24_2(j) + G24_3(j) + G24_4(j) +
						G24_5(j);
					if (j > int(nYield*AssFl)) {
						G24(j) = 0.;
					}
					G24_f = G24_f + G24(j);
					// -------------------------------------------------------------

					if ((abs(H51el * N1y_n(j) * N2y_n(j)) < 0.0000001) && (abs(Fd_1(j)) < 0.0000001)) {
						G25_1(j) = 0.0;
					}
					else {
						G25_1(j) = H51el * N1y_n(j) * N2y_n(j) / Fd_1(j);
					}
					if ((abs(H52el * N2y_n(j) * N2y_n(j)) < 0.0000001) && (abs(Fd_2(j)) < 0.0000001)) {
						G25_2(j) = 0.0;
					}
					else {
						G25_2(j) = H52el * N2y_n(j) * N2y_n(j) / Fd_2(j);
					}
					if ((abs(H53el * N3y_n(j) * N2y_n(j)) < 0.0000001) && (abs(Fd_3(j)) < 0.0000001)) {
						G25_3(j) = 0.0;
					}
					else {
						G25_3(j) = H53el * N3y_n(j) * N2y_n(j) / Fd_3(j);
					}
					if ((abs(H54el * N4y_n(j) * N2y_n(j)) < 0.0000001) && (abs(Fd_4(j)) < 0.0000001)) {
						G25_4(j) = 0.0;
					}
					else {
						G25_4(j) = H54el * N4y_n(j) * N2y_n(j) / Fd_4(j);
					}
					if ((abs(H55el * N5y_n(j) * N2y_n(j)) < 0.0000001) && (abs(Fd_5(j)) < 0.0000001)) {
						G25_5(j) = 0.0;
					}
					else {
						G25_5(j) = H55el * N5y_n(j) * N2y_n(j) / Fd_5(j);
					}
					G25(j) = G25_1(j) + G25_2(j) + G25_3(j) + G25_4(j) +
						G25_5(j);
					if (j >= int(nYield*AssFl)) {
						G25(j) = 0.;
					}
					G25_f = G25_f + G25(j);
					// -------------------------------------------------------------
					// -------------------------------------------------------------

					if ((abs(H11el * N1y_n(j) * N3y_n(j)) < 0.0000001) && (abs(Fd_1(j)) < 0.0000001)) {
						G31_1(j) = 0.0;
					}
					else {
						G31_1(j) = H11el * N1y_n(j) * N3y_n(j) / Fd_1(j);
					}
					if ((abs(H12el * N2y_n(j) * N3y_n(j)) < 0.0000001) && (abs(Fd_2(j)) < 0.0000001)) {
						G31_2(j) = 0.0;
					}
					else {
						G31_2(j) = H12el * N2y_n(j) * N3y_n(j) / Fd_2(j);
					}
					if ((abs(H13el * N3y_n(j) * N3y_n(j)) < 0.0000001) && (abs(Fd_3(j)) < 0.0000001)) {
						G31_3(j) = 0.0;
					}
					else {
						G31_3(j) = H13el * N3y_n(j) * N3y_n(j) / Fd_3(j);
					}
					if ((abs(H14el * N4y_n(j) * N3y_n(j)) < 0.0000001) && (abs(Fd_4(j)) < 0.0000001)) {
						G31_4(j) = 0.0;
					}
					else {
						G31_4(j) = H14el * N4y_n(j) * N3y_n(j) / Fd_4(j);
					}
					if ((abs(H15el * N5y_n(j) * N3y_n(j)) < 0.0000001) && (abs(Fd_5(j)) < 0.0000001)) {
						G31_5(j) = 0.0;
					}
					else {
						G31_5(j) = H15el * N5y_n(j) * N3y_n(j) / Fd_5(j);
					}
					G31(j) = G31_1(j) + G31_2(j) + G31_3(j) + G31_4(j) +
						G31_5(j);
					if (j >= int(nYield*AssFl)) {
						G31(j) = 0.;
					}
					G31_f = G31_f + G31(j);
					// -------------------------------------------------------------

					if ((abs(H21el * N1y_n(j) * N3y_n(j)) < 0.0000001) && (abs(Fd_1(j)) < 0.0000001)) {
						G32_1(j) = 0.0;
					}
					else {
						G32_1(j) = H21el * N1y_n(j) * N3y_n(j) / Fd_1(j);
					}
					if ((abs(H22el * N2y_n(j) * N3y_n(j)) < 0.0000001) && (abs(Fd_2(j)) < 0.0000001)) {
						G32_2(j) = 0.0;
					}
					else {
						G32_2(j) = H22el * N2y_n(j) * N3y_n(j) / Fd_2(j);
					}
					if ((abs(H23el * N3y_n(j) * N3y_n(j)) < 0.0000001) && (abs(Fd_3(j)) < 0.0000001)) {
						G32_3(j) = 0.0;
					}
					else {
						G32_3(j) = H23el * N3y_n(j) * N3y_n(j) / Fd_3(j);
					}
					if ((abs(H24el * N4y_n(j) * N3y_n(j)) < 0.0000001) && (abs(Fd_4(j)) < 0.0000001)) {
						G32_4(j) = 0.0;
					}
					else {
						G32_4(j) = H24el * N4y_n(j) * N3y_n(j) / Fd_4(j);
					}
					if ((abs(H25el * N5y_n(j) * N3y_n(j)) < 0.0000001) && (abs(Fd_5(j)) < 0.0000001)) {
						G32_5(j) = 0.0;
					}
					else {
						G32_5(j) = H25el * N5y_n(j) * N3y_n(j) / Fd_5(j);
					}
					G32(j) = G32_1(j) + G32_2(j) + G32_3(j) + G32_4(j) +
						G32_5(j);
					if (j >= int(nYield*AssFl)) {
						G32(j) = 0.;
					}
					G32_f = G32_f + G32(j);
					// -------------------------------------------------------------

					if ((abs(H31el * N1y_n(j) * N3y_n(j)) < 0.0000001) && (abs(Fd_1(j)) < 0.0000001)) {
						G33_1(j) = 0.0;
					}
					else {
						G33_1(j) = H31el * N1y_n(j) * N3y_n(j) / Fd_1(j);
					}
					if ((abs(H32el * N2y_n(j) * N3y_n(j)) < 0.0000001) && (abs(Fd_2(j)) < 0.0000001)) {
						G33_2(j) = 0.0;
					}
					else {
						G33_2(j) = H32el * N2y_n(j) * N3y_n(j) / Fd_2(j);
					}
					if ((abs(H33el * N3y_n(j) * N3y_n(j)) < 0.0000001) && (abs(Fd_3(j)) < 0.0000001)) {
						G33_3(j) = 0.0;
					}
					else {
						G33_3(j) = H33el * N3y_n(j) * N3y_n(j) / Fd_3(j);
					}
					if ((abs(H34el * N4y_n(j) * N3y_n(j)) < 0.0000001) && (abs(Fd_4(j)) < 0.0000001)) {
						G33_4(j) = 0.0;
					}
					else {
						G33_4(j) = H34el * N4y_n(j) * N3y_n(j) / Fd_4(j);
					}
					if ((abs(H35el * N5y_n(j) * N3y_n(j)) < 0.0000001) && (abs(Fd_5(j)) < 0.0000001)) {
						G33_5(j) = 0.0;
					}
					else {
						G33_5(j) = H35el * N5y_n(j) * N3y_n(j) / Fd_5(j);
					}
					G33(j) = G33_1(j) + G33_2(j) + G33_3(j) + G33_4(j) +
						G33_5(j);
					G33_f = G33_f + G33(j);
					// -------------------------------------------------------------

					if ((abs(H41el * N1y_n(j) * N3y_n(j)) < 0.0000001) && (abs(Fd_1(j)) < 0.0000001)) {
						G34_1(j) = 0.0;
					}
					else {
						G34_1(j) = H41el * N1y_n(j) * N3y_n(j) / Fd_1M(j);
					}
					if ((abs(H42el * N2y_n(j) * N3y_n(j)) < 0.0000001) && (abs(Fd_2(j)) < 0.0000001)) {
						G34_2(j) = 0.0;
					}
					else {
						G34_2(j) = H42el * N2y_n(j) * N3y_n(j) / Fd_2(j);
					}
					if ((abs(H43el * N3y_n(j) * N3y_n(j)) < 0.0000001) && (abs(Fd_3(j)) < 0.0000001)) {
						G34_3(j) = 0.0;
					}
					else {
						G34_3(j) = H43el * N3y_n(j) * N3y_n(j) / Fd_3(j);
					}
					if ((abs(H44el * N4y_n(j) * N3y_n(j)) < 0.0000001) && (abs(Fd_4(j)) < 0.0000001)) {
						G34_4(j) = 0.0;
					}
					else {
						G34_4(j) = H44el * N4y_n(j) * N3y_n(j) / Fd_4(j);
					}
					if ((abs(H45el * N5y_n(j) * N3y_n(j)) < 0.0000001) && (abs(Fd_5(j)) < 0.0000001)) {
						G34_5(j) = 0.0;
					}
					else {
						G34_5(j) = H45el * N5y_n(j) * N3y_n(j) / Fd_5(j);
					}
					G34(j) = G34_1(j) + G34_2(j) + G34_3(j) + G34_4(j) +
						G34_5(j);
					if (j >= int(nYield*AssFl)) {
						G34(j) = 0.;
					}
					G34_f = G34_f + G34(j);
					// -------------------------------------------------------------

					if ((abs(H51el * N1y_n(j) * N3y_n(j)) < 0.0000001) && (abs(Fd_1(j)) < 0.0000001)) {
						G35_1(j) = 0.0;
					}
					else {
						G35_1(j) = H51el * N1y_n(j) * N3y_n(j) / Fd_1(j);
					}
					if ((abs(H52el * N2y_n(j) * N3y_n(j)) < 0.0000001) && (abs(Fd_2(j)) < 0.0000001)) {
						G35_2(j) = 0.0;
					}
					else {
						G35_2(j) = H52el * N2y_n(j) * N3y_n(j) / Fd_2(j);
					}
					if ((abs(H53el * N3y_n(j) * N3y_n(j)) < 0.0000001) && (abs(Fd_3(j)) < 0.0000001)) {
						G35_3(j) = 0.0;
					}
					else {
						G35_3(j) = H53el * N3y_n(j) * N3y_n(j) / Fd_3(j);
					}
					if ((abs(H54el * N4y_n(j) * N3y_n(j)) < 0.0000001) && (abs(Fd_4(j)) < 0.0000001)) {
						G35_4(j) = 0.0;
					}
					else {
						G35_4(j) = H54el * N4y_n(j) * N3y_n(j) / Fd_4(j);
					}
					if ((abs(H55el * N5y_n(j) * N3y_n(j)) < 0.0000001) && (abs(Fd_5(j)) < 0.0000001)) {
						G35_5(j) = 0.0;
					}
					else {
						G35_5(j) = H55el * N5y_n(j) * N3y_n(j) / Fd_5(j);
					}
					G35(j) = G35_1(j) + G35_2(j) + G35_3(j) + G35_4(j) +
						G35_5(j);
					if (j >= int(nYield*AssFl)) {
						G35(j) = 0.;
					}
					G35_f = G35_f + G35(j);
					// -------------------------------------------------------------
					// -------------------------------------------------------------

					if ((abs(H11el * N1y_n(j) * N4y_n(j)) < 0.0000001) && (abs(Fd_1M(j)) < 0.0000001)) {
						G41_1(j) = 0.0;
					}
					else {
						G41_1(j) = H11el * N1y_n(j) * N4y_n(j) / Fd_1M(j);
					}
					if ((abs(H12el * N2y_n(j) * N4y_n(j)) < 0.0000001) && (abs(Fd_2M(j)) < 0.0000001)) {
						G41_2(j) = 0.0;
					}
					else {
						G41_2(j) = H12el * N2y_n(j) * N4y_n(j) / Fd_2M(j);
					}
					if ((abs(H13el * N3y_n(j) * N4y_n(j)) < 0.0000001) && (abs(Fd_3M(j)) < 0.0000001)) {
						G41_3(j) = 0.0;
					}
					else {
						G41_3(j) = H13el * N3y_n(j) * N4y_n(j) / Fd_3M(j);
					}
					if ((abs(H14el * N4y_n(j) * N4y_n(j)) < 0.0000001) && (abs(Fd_4(j)) < 0.0000001)) {
						G41_4(j) = 0.0;
					}
					else {
						G41_4(j) = H14el * N4y_n(j) * N4y_n(j) / Fd_4(j);
					}
					if ((abs(H15el * N5y_n(j) * N4y_n(j)) < 0.0000001) && (abs(Fd_5(j)) < 0.0000001)) {
						G41_5(j) = 0.0;
					}
					else {
						G41_5(j) = H15el * N5y_n(j) * N4y_n(j) / Fd_5(j);
					}
					G41(j) = G41_1(j) + G41_2(j) + G41_3(j) + G41_4(j) +
						G41_5(j);
					if (j >= int(nYield*AssFl)) {
						G41(j) = 0.;
					}
					G41_f = G41_f + G41(j);
					// -------------------------------------------------------------

					if ((abs(H21el * N1y_n(j) * N4y_n(j)) < 0.0000001) && (abs(Fd_1M(j)) < 0.0000001)) {
						G42_1(j) = 0.0;
					}
					else {
						G42_1(j) = H21el * N1y_n(j) * N4y_n(j) / Fd_1M(j);
					}
					if ((abs(H22el * N2y_n(j) * N4y_n(j)) < 0.0000001) && (abs(Fd_2M(j)) < 0.0000001)) {
						G42_2(j) = 0.0;
					}
					else {
						G42_2(j) = H22el * N2y_n(j) * N4y_n(j) / Fd_2M(j);
					}
					if ((abs(H23el * N3y_n(j) * N4y_n(j)) < 0.0000001) && (abs(Fd_3M(j)) < 0.0000001)) {
						G42_3(j) = 0.0;
					}
					else {
						G42_3(j) = H23el * N3y_n(j) * N4y_n(j) / Fd_3M(j);
					}
					if ((abs(H24el * N4y_n(j) * N4y_n(j)) < 0.0000001) && (abs(Fd_4(j)) < 0.0000001)) {
						G42_4(j) = 0.0;
					}
					else {
						G42_4(j) = H24el * N4y_n(j) * N4y_n(j) / Fd_4(j);
					}
					if ((abs(H25el * N5y_n(j) * N4y_n(j)) < 0.0000001) && (abs(Fd_5(j)) < 0.0000001)) {
						G42_5(j) = 0.0;
					}
					else {
						G42_5(j) = H25el * N5y_n(j) * N4y_n(j) / Fd_5(j);
					}
					G42(j) = G42_1(j) + G42_2(j) + G42_3(j) + G42_4(j) +
						G42_5(j);
					if (j >= int(nYield*AssFl)) {
						G42(j) = 0.;
					}
					G42_f = G42_f + G42(j);
					// -------------------------------------------------------------

					if ((abs(H31el * N1y_n(j) * N4y_n(j)) < 0.0000001) && (abs(Fd_1M(j)) < 0.0000001)) {
						G43_1(j) = 0.0;
					}
					else {
						G43_1(j) = H31el * N1y_n(j) * N4y_n(j) / Fd_1M(j);
					}
					if ((abs(H32el * N2y_n(j) * N4y_n(j)) < 0.0000001) && (abs(Fd_2M(j)) < 0.0000001)) {
						G43_2(j) = 0.0;
					}
					else {
						G43_2(j) = H32el * N2y_n(j) * N4y_n(j) / Fd_2M(j);
					}
					if ((abs(H33el * N3y_n(j) * N4y_n(j)) < 0.0000001) && (abs(Fd_3M(j)) < 0.0000001)) {
						G43_3(j) = 0.0;
					}
					else {
						G43_3(j) = H33el * N3y_n(j) * N4y_n(j) / Fd_3M(j);
					}
					if ((abs(H34el * N4y_n(j) * N4y_n(j)) < 0.0000001) && (abs(Fd_4(j)) < 0.0000001)) {
						G43_4(j) = 0.0;
					}
					else {
						G43_4(j) = H34el * N4y_n(j) * N4y_n(j) / Fd_4(j);
					}
					if ((abs(H35el * N5y_n(j) * N4y_n(j)) < 0.0000001) && (abs(Fd_5(j)) < 0.0000001)) {
						G43_5(j) = 0.0;
					}
					else {
						G43_5(j) = H35el * N5y_n(j) * N4y_n(j) / Fd_5(j);
					}
					G43(j) = G43_1(j) + G43_2(j) + G43_3(j) + G43_4(j) +
						G43_5(j);
					if (j >= int(nYield*AssFl)) {
						G43(j) = 0.;
					}
					G43_f = G43_f + G43(j);
					// -------------------------------------------------------------

					if ((abs(H41el * N1y_n(j) * N4y_n(j)) < 0.0000001) && (abs(Fd_1M(j)) < 0.0000001)) {
						G44_1(j) = 0.0;
					}
					else {
						G44_1(j) = H41el * N1y_n(j) * N4y_n(j) / Fd_1M(j);
					}
					if ((abs(H42el * N2y_n(j) * N4y_n(j)) < 0.0000001) && (abs(Fd_2M(j)) < 0.0000001)) {
						G44_2(j) = 0.0;
					}
					else {
						G44_2(j) = H42el * N2y_n(j) * N4y_n(j) / Fd_2M(j);
					}
					if ((abs(H43el * N3y_n(j) * N4y_n(j)) < 0.0000001) && (abs(Fd_3M(j)) < 0.0000001)) {
						G44_3(j) = 0.0;
					}
					else {
						G44_3(j) = H43el * N3y_n(j) * N4y_n(j) / Fd_3M(j);
					}
					if ((abs(H44el * N4y_n(j) * N4y_n(j)) < 0.0000001) && (abs(Fd_4(j)) < 0.0000001)) {
						G44_4(j) = 0.0;
					}
					else {
						G44_4(j) = H44el * N4y_n(j) * N4y_n(j) / Fd_4(j);
					}
					if ((abs(H45el * N5y_n(j) * N4y_n(j)) < 0.0000001) && (abs(Fd_5(j)) < 0.0000001)) {
						G44_5(j) = 0.0;
					}
					else {
						G44_5(j) = H45el * N5y_n(j) * N4y_n(j) / Fd_5(j);
					}
					G44(j) = G44_1(j) + G44_2(j) + G44_3(j) + G44_4(j) +
						G44_5(j);
					G44_f = G44_f + G44(j);
					// -------------------------------------------------------------

					if ((abs(H51el * N1y_n(j) * N4y_n(j)) < 0.0000001) && (abs(Fd_1M(j)) < 0.0000001)) {
						G45_1(j) = 0.0;
					}
					else {
						G45_1(j) = H51el * N1y_n(j) * N4y_n(j) / Fd_1M(j);
					}
					if ((abs(H52el * N2y_n(j) * N4y_n(j)) < 0.0000001) && (abs(Fd_2M(j)) < 0.0000001)) {
						G45_2(j) = 0.0;
					}
					else {
						G45_2(j) = H52el * N2y_n(j) * N4y_n(j) / Fd_2M(j);
					}
					if ((abs(H53el * N3y_n(j) * N4y_n(j)) < 0.0000001) && (abs(Fd_3M(j)) < 0.0000001)) {
						G45_3(j) = 0.0;
					}
					else {
						G45_3(j) = H53el * N3y_n(j) * N4y_n(j) / Fd_3M(j);
					}
					if ((abs(H54el * N4y_n(j) * N4y_n(j)) < 0.0000001) && (abs(Fd_4(j)) < 0.0000001)) {
						G45_4(j) = 0.0;
					}
					else {
						G45_4(j) = H54el * N4y_n(j) * N4y_n(j) / Fd_4(j);
					}
					if ((abs(H55el * N5y_n(j) * N4y_n(j)) < 0.0000001) && (abs(Fd_5(j)) < 0.0000001)) {
						G45_5(j) = 0.0;
					}
					else {
						G45_5(j) = H55el * N5y_n(j) * N4y_n(j) / Fd_5(j);
					}
					G45(j) = G45_1(j) + G45_2(j) + G45_3(j) + G45_4(j) +
						G45_5(j);
					if (j >= int(nYield*AssFl)) {
						G45(j) = 0.;
					}
					G45_f = G45_f + G45(j);
					// -------------------------------------------------------------
					// -------------------------------------------------------------

					if ((abs(H11el * N1y_n(j) * N5y_n(j)) < 0.0000001) && (abs(Fd_1M(j)) < 0.0000001)) {
						G51_1(j) = 0.0;
					}
					else {
						G51_1(j) = H11el * N1y_n(j) * N5y_n(j) / Fd_1M(j);
					}
					if ((abs(H12el * N2y_n(j) * N5y_n(j)) < 0.0000001) && (abs(Fd_2M(j)) < 0.0000001)) {
						G51_2(j) = 0.0;
					}
					else {
						G51_2(j) = H12el * N2y_n(j) * N5y_n(j) / Fd_2M(j);
					}
					if ((abs(H13el * N3y_n(j) * N5y_n(j)) < 0.0000001) && (abs(Fd_3M(j)) < 0.0000001)) {
						G51_3(j) = 0.0;
					}
					else {
						G51_3(j) = H13el * N3y_n(j) * N5y_n(j) / Fd_3M(j);
					}
					if ((abs(H14el * N4y_n(j) * N5y_n(j)) < 0.0000001) && (abs(Fd_4(j)) < 0.0000001)) {
						G51_4(j) = 0.0;
					}
					else {
						G51_4(j) = H14el * N4y_n(j) * N5y_n(j) / Fd_4(j);
					}
					if ((abs(H15el * N5y_n(j) * N5y_n(j)) < 0.0000001) && (abs(Fd_5(j)) < 0.0000001)) {
						G51_5(j) = 0.0;
					}
					else {
						G51_5(j) = H15el * N5y_n(j) * N5y_n(j) / Fd_5(j);
					}
					G51(j) = G51_1(j) + G51_2(j) + G51_3(j) + G51_4(j) +
						G51_5(j);
					if (j >= int(nYield*AssFl)) {
						G51(j) = 0.;
					}
					G51_f = G51_f + G51(j);
					// -------------------------------------------------------------

					if ((abs(H21el * N1y_n(j) * N5y_n(j)) < 0.0000001) && (abs(Fd_1M(j)) < 0.0000001)) {
						G52_1(j) = 0.0;
					}
					else {
						G52_1(j) = H21el * N1y_n(j) * N5y_n(j) / Fd_1M(j);
					}
					if ((abs(H22el * N2y_n(j) * N5y_n(j)) < 0.0000001) && (abs(Fd_2M(j)) < 0.0000001)) {
						G52_2(j) = 0.0;
					}
					else {
						G52_2(j) = H22el * N2y_n(j) * N5y_n(j) / Fd_2M(j);
					}
					if ((abs(H23el * N3y_n(j) * N5y_n(j)) < 0.0000001) && (abs(Fd_3M(j)) < 0.0000001)) {
						G52_3(j) = 0.0;
					}
					else {
						G52_3(j) = H23el * N3y_n(j) * N5y_n(j) / Fd_3M(j);
					}
					if ((abs(H24el * N4y_n(j) * N5y_n(j)) < 0.0000001) && (abs(Fd_4(j)) < 0.0000001)) {
						G52_4(j) = 0.0;
					}
					else {
						G52_4(j) = H24el * N4y_n(j) * N5y_n(j) / Fd_4(j);
					}
					if ((abs(H25el * N5y_n(j) * N5y_n(j)) < 0.0000001) && (abs(Fd_5(j)) < 0.0000001)) {
						G52_5(j) = 0.0;
					}
					else {
						G52_5(j) = H25el * N5y_n(j) * N5y_n(j) / Fd_5(j);
					}
					G52(j) = G52_1(j) + G52_2(j) + G52_3(j) + G52_4(j) +
						G52_5(j);
					if (j >= int(nYield*AssFl)) {
						G52(j) = 0.;
					}
					G52_f = G52_f + G52(j);
					// -------------------------------------------------------------

					if ((abs(H31el * N1y_n(j) * N5y_n(j)) < 0.0000001) && (abs(Fd_1M(j)) < 0.0000001)) {
						G53_1(j) = 0.0;
					}
					else {
						G53_1(j) = H31el * N1y_n(j) * N5y_n(j) / Fd_1M(j);
					}
					if ((abs(H32el * N2y_n(j) * N5y_n(j)) < 0.0000001) && (abs(Fd_2M(j)) < 0.0000001)) {
						G53_2(j) = 0.0;
					}
					else {
						G53_2(j) = H32el * N2y_n(j) * N5y_n(j) / Fd_2M(j);
					}
					if ((abs(H33el * N3y_n(j) * N5y_n(j)) < 0.0000001) && (abs(Fd_3M(j)) < 0.0000001)) {
						G53_3(j) = 0.0;
					}
					else {
						G53_3(j) = H33el * N3y_n(j) * N5y_n(j) / Fd_3M(j);
					}
					if ((abs(H34el * N4y_n(j) * N5y_n(j)) < 0.0000001) && (abs(Fd_4(j)) < 0.0000001)) {
						G53_4(j) = 0.0;
					}
					else {
						G53_4(j) = H34el * N4y_n(j) * N5y_n(j) / Fd_4(j);
					}
					if ((abs(H35el * N5y_n(j) * N5y_n(j)) < 0.0000001) && (abs(Fd_5(j)) < 0.0000001)) {
						G53_5(j) = 0.0;
					}
					else {
						G53_5(j) = H35el * N5y_n(j) * N5y_n(j) / Fd_5(j);
					}
					G53(j) = G53_1(j) + G53_2(j) + G53_3(j) + G53_4(j) +
						G53_5(j);
					G53_f = G53_f + G53(j);
					if (j >= int(nYield*AssFl)) {
						G53(j) = 0.;
					}
					// -------------------------------------------------------------

					if ((abs(H41el * N1y_n(j) * N5y_n(j)) < 0.0000001) && (abs(Fd_1M(j)) < 0.0000001)) {
						G54_1(j) = 0.0;
					}
					else {
						G54_1(j) = H41el * N1y_n(j) * N5y_n(j) / Fd_1M(j);
					}
					if ((abs(H42el * N2y_n(j) * N5y_n(j)) < 0.0000001) && (abs(Fd_2M(j)) < 0.0000001)) {
						G54_2(j) = 0.0;
					}
					else {
						G54_2(j) = H42el * N2y_n(j) * N5y_n(j) / Fd_2M(j);
					}
					if ((abs(H43el * N3y_n(j) * N5y_n(j)) < 0.0000001) && (abs(Fd_3M(j)) < 0.0000001)) {
						G54_3(j) = 0.0;
					}
					else {
						G54_3(j) = H43el * N3y_n(j) * N5y_n(j) / Fd_3M(j);
					}
					if ((abs(H44el * N4y_n(j) * N5y_n(j)) < 0.0000001) && (abs(Fd_4(j)) < 0.0000001)) {
						G54_4(j) = 0.0;
					}
					else {
						G54_4(j) = H44el * N4y_n(j) * N5y_n(j) / Fd_4(j);
					}
					if ((abs(H45el * N5y_n(j) * N5y_n(j)) < 0.0000001) && (abs(Fd_5(j)) < 0.0000001)) {
						G54_5(j) = 0.0;
					}
					else {
						G54_5(j) = H45el * N5y_n(j) * N5y_n(j) / Fd_5(j);
					}
					G54(j) = G54_1(j) + G54_2(j) + G54_3(j) + G54_4(j) +
						G54_5(j);
					G54_f = G54_f + G54(j);
					if (j >= int(nYield*AssFl)) {
						G54(j) = 0.;
					}
					// -------------------------------------------------------------

					if ((abs(H51el * N1y_n(j) * N5y_n(j)) < 0.0000001) && (abs(Fd_1M(j)) < 0.0000001)) {
						G55_1(j) = 0.0;
					}
					else {
						G55_1(j) = H51el * N1y_n(j) * N5y_n(j) / Fd_1M(j);
					}
					if ((abs(H52el * N2y_n(j) * N5y_n(j)) < 0.0000001) && (abs(Fd_2M(j)) < 0.0000001)) {
						G55_2(j) = 0.0;
					}
					else {
						G55_2(j) = H52el * N2y_n(j) * N5y_n(j) / Fd_2M(j);
					}
					if ((abs(H53el * N3y_n(j) * N5y_n(j)) < 0.0000001) && (abs(Fd_3M(j)) < 0.0000001)) {
						G55_3(j) = 0.0;
					}
					else {
						G55_3(j) = H53el * N3y_n(j) * N5y_n(j) / Fd_3M(j);
					}
					if ((abs(H54el * N4y_n(j) * N5y_n(j)) < 0.0000001) && (abs(Fd_4(j)) < 0.0000001)) {
						G55_4(j) = 0.0;
					}
					else {
						G55_4(j) = H54el * N4y_n(j) * N5y_n(j) / Fd_4(j);
					}
					if ((abs(H55el * N5y_n(j) * N5y_n(j)) < 0.0000001) && (abs(Fd_5(j)) < 0.0000001)) {
						G55_5(j) = 0.0;
					}
					else {
						G55_5(j) = H55el * N5y_n(j) * N5y_n(j) / Fd_5(j);
					}
					G55(j) = G55_1(j) + G55_2(j) + G55_3(j) + G55_4(j) +
						G55_5(j);
					G55_f = G55_f + G55(j);


				} // for j = 1: 1: nFlows(cont)


				// Tangent compliance matrix ---------------------------------------
				Ct(0, 0) = 1.0 / H11el * (1.0 + G11_f);
				Ct(0, 1) = 1.0 / H11el * G12_f * Nas;
				Ct(0, 2) = 1.0 / H11el * G13_f * Nas;
				Ct(0, 3) = 0.0;
				Ct(0, 4) = 0.0;
				Ct(0, 5) = 0.0;

				Ct(1, 0) = 1.0 / H22el * G21_f * Nas;
				Ct(1, 1) = 1.0 / H22el * (1.0 + G22_f);
				Ct(1, 2) = 1.0 / H22el * G23_f * Nas;
				Ct(1, 3) = 0.0;
				Ct(1, 4) = 0.0;
				Ct(1, 5) = 0.0;

				Ct(2, 0) = 1.0 / H33el * G31_f * Nas;
				Ct(2, 1) = 1.0 / H33el * G32_f * Nas;
				Ct(2, 2) = 1.0 / H22el * (1.0 + G33_f);
				Ct(2, 3) = 0.0;
				Ct(2, 4) = 0.0;
				Ct(2, 5) = 0.0;

				Ct(3, 0) = 1000000.0;
				Ct(3, 1) = 1000000.0;
				Ct(3, 2) = 1000000.0;
				Ct(3, 3) = 1000000.0;
				Ct(3, 4) = 1000000.0;
				Ct(3, 5) = 1000000.0;

				Ct(4, 0) = 1000000.0;
				Ct(4, 1) = 1000000.0;
				Ct(4, 2) = 1000000.0;
				Ct(4, 3) = 1000000.0;
				Ct(4, 4) = 1000000.0;
				Ct(4, 5) = 1000000.0;

				Ct(5, 0) = 1000000.0;
				Ct(5, 1) = 1000000.0;
				Ct(5, 2) = 1000000.0;
				Ct(5, 3) = 1000000.0;
				Ct(5, 4) = 1000000.0;
				Ct(5, 5) = 1000000.0;

				for (int ii = 3; ii < 6; ii++) {
					for (int jj = 3; jj < 6; jj++) {
						if (ii != jj) {
							Ct(ii, jj) = 0.0;
						}
					}
				}

				// tangent stiffness matrix ----------------------------------------
				Ct.Invert(Kt);

				Kt(0, 3) = 0.0;
				Kt(0, 4) = 0.0;
				Kt(0, 5) = 0.0;

				Kt(1, 3) = 0.0;
				Kt(1, 4) = 0.0;
				Kt(1, 5) = 0.0;

				Kt(2, 3) = 0.0;
				Kt(2, 4) = 0.0;
				Kt(2, 5) = 0.0;

				Kt(3, 0) = 0.0;
				Kt(3, 1) = 0.0;
				Kt(3, 2) = 0.0;
				Kt(3, 3) = 0.0;
				Kt(3, 4) = 0.0;
				Kt(3, 5) = 0.0;

				Kt(4, 0) = 0.0;
				Kt(4, 1) = 0.0;
				Kt(4, 2) = 0.0;
				Kt(4, 3) = 0.0;
				Kt(4, 4) = 0.0;
				Kt(4, 5) = 0.0;

				Kt(5, 0) = 0.0;
				Kt(5, 1) = 0.0;
				Kt(5, 2) = 0.0;
				Kt(5, 3) = 0.0;
				Kt(5, 4) = 0.0;
				Kt(5, 5) = 0.0;

				Kt(0, 1) = Kt(0, 1) * Nas;
				Kt(0, 2) = Kt(0, 2) * Nas;
				Kt(1, 0) = Kt(1, 0) * Nas;
				Kt(1, 2) = Kt(1, 2) * Nas;
				Kt(3, 0) = Kt(3, 0) * Nas;
				Kt(3, 1) = Kt(3, 1) * Nas;

				// incremental response
				Q01 = Cstress(0);
				Q02 = Cstress(1);
				Q03 = Cstress(2);
				Q04 = Cstress(3);
				Q05 = Cstress(4);
				Q06 = Cstress(5);

				q01 = Cstrain(0);
				q02 = Cstrain(1);
				q03 = Cstrain(2);
				q04 = Cstrain(3);
				q05 = Cstrain(4);
				q06 = Cstrain(5);

				Q1 = Q01 +
					Kt(0, 0) * (q1 - q01) + Kt(0, 1) * (q2 - q02) + Kt(0, 2) * (q3 - q03) +
					Kt(0, 1) * pow(H22el, -1.0) * H21el * (q1 - q01) +
					Kt(1, 2) * pow(H33el, -1.0) * H31el * (q1 - q01);
				Q2 = Q02 +
					Kt(1, 0) * (q1 - q01) + Kt(1, 1) * (q2 - q02) + Kt(1, 2) * (q3 - q03) +
					Kt(1, 0) * pow(H11el, -1.0) * H12el * (q2 - q02) +
					Kt(1, 2) * pow(H33el, -1.0) * H32el * (q2 - q02);
				Q3 = Q03 +
					Kt(2, 0) * (q1 - q01) + Kt(2, 1) * (q2 - q02) + Kt(2, 2) * (q3 - q03) +
					Kt(2, 0) * pow(H11el, -1.0) * H13el * (q3 - q03) +
					Kt(2, 1) * pow(H22el, -1.0) * H23el * (q3 - q03);
				Q4 = 0.0;
				Q5 = 0.0;
				Q6 = 0.0;

				TrialStressElasticFlag = 1;

			}

			else { // purely elastic response

				Q01 = Cstress(0);
				Q02 = Cstress(1);
				Q03 = Cstress(2);
				Q04 = Cstress(3);
				Q05 = Cstress(4);
				Q06 = Cstress(5);

				q01 = Cstrain(0);
				q02 = Cstrain(1);
				q03 = Cstrain(2);
				q04 = Cstrain(3);
				q05 = Cstrain(4);
				q06 = Cstrain(5);

				Kt(0, 0) = H11el;
				Kt(0, 1) = 0.0;
				Kt(0, 2) = 0.0;
				Kt(0, 3) = 0.0;
				Kt(0, 4) = 0.0;
				Kt(0, 5) = 0.0;
				Kt(1, 0) = 0.0;
				Kt(1, 1) = H22el;
				Kt(1, 2) = 0.0;
				Kt(1, 3) = 0.0;
				Kt(1, 4) = 0.0;
				Kt(1, 5) = 0.0;
				Kt(2, 0) = 0.0;
				Kt(2, 1) = 0.0;
				Kt(2, 2) = H33el;
				Kt(2, 3) = 0.0;
				Kt(2, 4) = 0.0;
				Kt(2, 5) = 0.0;
				Kt(3, 0) = 0.0;
				Kt(3, 1) = 0.0;
				Kt(3, 2) = 0.0;
				Kt(3, 3) = 0.0;
				Kt(3, 4) = 0.0;
				Kt(3, 5) = 0.0;
				Kt(4, 0) = 0.0;
				Kt(4, 1) = 0.0;
				Kt(4, 2) = 0.0;
				Kt(4, 3) = 0.0;
				Kt(4, 4) = 0.0;
				Kt(4, 5) = 0.0;
				Kt(5, 0) = 0.0;
				Kt(5, 1) = 0.0;
				Kt(5, 2) = 0.0;
				Kt(5, 3) = 0.0;
				Kt(5, 4) = 0.0;
				Kt(5, 5) = 0.0;

				Q1 = Q01 +
					Kt(0, 0) * (q1 - q01) + Kt(0, 1) * (q2 - q02) + Kt(0, 2) * (q3 - q03);
				Q2 = Q02 +
					Kt(1, 0) * (q1 - q01) + Kt(1, 1) * (q2 - q02) + Kt(1, 2) * (q3 - q03);
				Q3 = Q03 +
					Kt(2, 0) * (q1 - q01) + Kt(2, 1) * (q2 - q02) + Kt(2, 2) * (q3 - q03);
				Q4 = 0.0;
				Q5 = 0.0;
				Q6 = 0.0;

				TrialStressElasticFlag = 0;
			}

			if ((abs(Q1) < tolForce) && (abs(Q2) < tolForce) && (abs(Q3) < tolForce)) {
				if ((abs(PCQy1(nYield - 1)) > 0.0) || (abs(PCQy2(nYield - 1)) > 0.0) || (abs(PCQy3(nYield - 1)) > 0.0)) {
					Q1 = 0;
					Q2 = 0;
					Q3 = 0;
				}
			}
			if ((abs(Qy1_iter(nYield - 1)) < tolForce) && (abs(Qy2_iter(nYield - 1)) < tolForce) && (abs(Qy3_iter(nYield - 1)) < tolForce)) {
				if ((abs(PCQy1(nYield - 1)) > 0.0) || (abs(PCQy2(nYield - 1)) > 0.0) || (abs(PCQy3(nYield - 1)) > 0.0)) {
					Q1 = Q01;
					Q2 = Q02;
					Q3 = Q03;
				}
			}
			if (abs(Q1) < tolForce) {
				Q1 = 0.0;
			}
			if (abs(Q2) < tolForce) {
				Q2 = 0.0;
			}
			if (abs(Q3) < tolForce) {
				Q3 = 0.0;
			}
			Q = pow(pow(Q1, 2.0) + pow(Q2, 2.0) + pow(Q3, 2.0), 0.5);


			// convergence check -----------------------------------------------
			if (abs(Q - Qprev) <= tol) {
				convergence = 1;
			}
			else {
				Qprev = Q;
				t1 = Q1;
				t2 = Q2;
				t3 = Q3;
				t4 = 0.0;
				t5 = 0.0;
				t6 = 0.0;
			}

			int ll = 1;
			// conditions to get convergence when the trial force oscillates around a plastic threshold
			double Conver = 0.0;
			for (int ll = 0; ll < nYield; ll++) {
				if (cont > 4) {
					if ((nFlows(cont) == ll) &&
						(nFlows(cont - 1) == ll - 1) &&
						(nFlows(cont - 2) == ll) &&
						(nFlows(cont - 3) == ll - 1) &&
						(nFlows(cont - 4) == ll)) {
						Conver = 1.0;
						break;
					}
				}
			}

			for (int ll = 0; ll < nYield; ll++) {
				if (cont > 5) {
					if ((nFlows(cont) > 0) && (nFlows(cont) == nFlows(cont - 2)) &&
						(nFlows(cont - 1) == nFlows(cont - 3))) {
						Conver = 1.0;
						break;
					}
				}
			}

			if (Conver == 1.0) {
				break;
			}

			if (cont > 2) {
				if ((nFlows(cont) == nYield) && (PrnFlows == nYield)) {
					Q1 = PStr1;
					Q2 = PStr2;
					Q3 = PStr3;
					break;
				}
			}

			if (cont == niter) {
				int ll = 1;
				Conver = 0.0;
				for (int ll = 0; ll < nYield; ll++) {
					if ((nFlows(cont) == ll - 1) &&
						(nFlows(cont - 1) == ll) &&
						(nFlows(cont - 2) == ll - 1)) {
						Conver = 1.0;
						break;
					}
					if ((nFlows(cont) == ll) &&
						(nFlows(cont - 1) == ll - 1) &&
						(nFlows(cont - 2) == ll)) {
						Conver = 1.0;
						break;
					}
				}
				if (Conver == 1.0) {
					break;
				}
				// ---------------------------------------------------------------------
				if (abs(Q - Qprev) > tol) {
					"disp(internal loop failed to converge)\n";
					return 0;
				}
			}

			t1_iter = Q1;
			t2_iter = Q2;
			t3_iter = Q3;

		}
		// INTERNAL LOOP FINISHED
		// -------------------------------------------------------------------------

		setFlows = nFlows(cont);

		double Qnorm = pow(pow(Q1, 2.0) + pow(Q2, 2.0) + pow(Q3, 2.0), 0.5);

		for (int j = 0; j < nYield - 1; j++) {
			csi1(j) = csi1_iter(j);
			csi2(j) = csi2_iter(j);
			csi3(j) = csi3_iter(j);
		}

		double DeltaQ1 = 0.0;
		double DeltaQ2 = 0.0;
		double DeltaQ3 = 0.0;
		
		Tstress(0) = Q1;
		Tstress(1) = Q2;
		Tstress(2) = Q3;
		Tstress(3) = 0.0;
		Tstress(4) = 0.0;
		Tstress(5) = 0.0;

		TStr1_2 = PStr1;
		TStr2_2 = PStr2;
		TStr3_2 = PStr3;

		// ---------------------------------------------------------------------
		// consistency condition and plastic flow

		double tolCons = 1.0;

		if (setFlows > 0) {
			for (int j = 0; j < setFlows; j++) {
				DQ1(j) = (Q1 - PStr1) * tolCons;
				DQ2(j) = (Q2 - PStr2) * tolCons;
				DQ3(j) = (Q3 - PStr3) * tolCons;
			}
		}
		// ---------------------------------------------------------------------


		for (int j = 0; j < nYield - 2; j++) {
			if (j < setFlows) {
				Tc1Comm(j) = Cc1(j) + DQ1(j);
				Tc2Comm(j) = Cc2(j) + DQ2(j);
				Tc3Comm(j) = Cc3(j) + DQ3(j);
			}
			else {
				Tc1Comm(j) = Cc1(j);
				Tc2Comm(j) = Cc2(j);
				Tc3Comm(j) = Cc3(j);
			}
		}

		int j = nYield - 1;
		DQ1(j) = 0.0;
		DQ2(j) = 0.0;
		DQ3(j) = 0.0;
		Tc1Comm(j) = Cc1(j) + DQ1(j);
		Tc2Comm(j) = Cc2(j) + DQ2(j);
		Tc3Comm(j) = Cc3(j) + DQ3(j);

		for (int i = 0; i < 6; i++) {
			for (int j = 0; j < 6; j++) {
				theTangent(i, j) = Kt(i, j);
			}
		}

		for (int j = 0; j < nYield; j++) {
			CQy1(j) = Qy1_iter(j);
			CQy2(j) = Qy2_iter(j);
			CQy3(j) = Qy3_iter(j);
		}

	}
	else {

		setFlows = CnFlows;

		Q1 = PStr1;
		Q2 = PStr2;
		Q3 = PStr3;

		Tstress(0) = PStr1;
		Tstress(1) = PStr2;
		Tstress(2) = PStr3;
		Tstress(3) = 0;
		Tstress(4) = 0;
		Tstress(5) = 0;

		TStr1_2 = PStr1_2;
		TStr2_2 = PStr2_2;
		TStr3_2 = PStr3_2;

		for (int j = 0; j < nYield - 2; j++) {
			if (j < setFlows) {
				Tc1Comm(j) = Cc1(j);
				Tc2Comm(j) = Cc2(j);
				Tc3Comm(j) = Cc3(j);
			}
			else {
				Tc1Comm(j) = Cc1(j);
				Tc2Comm(j) = Cc2(j);
				Tc3Comm(j) = Cc3(j);
			}
		}

		int j = nYield - 1;
		Tc1Comm(j) = Cc1(j);
		Tc2Comm(j) = Cc2(j);
		Tc3Comm(j) = Cc3(j);

		for (int i = 0; i < 6; i++) {
			for (int j = 0; j < 6; j++) {
				theTangent(i, j) = PKt(i, j);
			}
		}


		for (int j = 0; j < nYield; j++) {
			CQy1(j) = Qy1_iter(j);
			CQy2(j) = Qy2_iter(j);
			CQy3(j) = Qy3_iter(j);
		}

	}

	return 0;

};