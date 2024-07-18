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
// $Date: 2024-06-14 11:29:01 $
// $Source: /usr/local/cvs/OpenSees/SRC/material/uniaxial/ASDSteel1DMaterial.cpp,v $

// Massimo Petracca - ASDEA Software, Italy
//
// A Simple model for reinforcing steel
//

#ifndef ASDSteel1DMaterial_h
#define ASDSteel1DMaterial_h

#include <UniaxialMaterial.h>
#include <Vector.h>
#include <Matrix.h>
#include <cmath>

class ASDSteel1DMaterial : public UniaxialMaterial
{
public:
	// life-cycle
	ASDSteel1DMaterial(
		int _tag, /** the material tag */
		double _Es, /** steel Young's modulus */
		double _fy, /** steel yield strength */
		double _Hs, /** steel kinematic hardening ratio */
		double _fc, /** concrete compressive strength */
		double _diam, /** rebar diameter */
		double _Leff, /** effective bond length */
		bool _good_bond, /** good bond flag */
		double _mm, /** 1 millimeter in your unit system */
		double _N, /** 1 Newton in your unit system */
		UniaxialMaterial* _steel_mat, /** the steel material */
		UniaxialMaterial* _bond_mat /** the bond material */
		);
	ASDSteel1DMaterial();
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

private:
	// internal computation
	const Vector& getSteelFailure() const;
	const Vector& getBondResponse() const;
	const Vector& getEquivalentBondResponse() const;
	const Vector& getSteelResponse() const;
	const Vector& getImplexError() const;

private:
	// Steel Young's modulus
	double Es = 0.0;
	// Steel yield strength
	double fy = 0.0;
	// Steel kinematic hardening ratio
	double Hs = 0.0;
	// Concrete compressive strength
	double fc = 0.0;
	// Rebar diameter
	double diam = 0.0;
	// Effective bond length
	double Leff = 0.0;
	// good bond flag
	bool good_bond = false;
	// 1 millimeter in your unit system
	double mm = 1.0;
	// 1 Newton in your unit system
	double N = 1.0;
	
	// strain, stress and tangent
	double strain = 0.0;
	double strain_commit = 0.0;
	
	// other variables for output purposes
	double energy = 0.0;
	
	// the materials
	UniaxialMaterial* steel_mat = nullptr;
	UniaxialMaterial* bond_mat = nullptr;
};

#endif
