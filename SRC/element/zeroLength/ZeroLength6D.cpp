/*
Developed by: Davide Noè Gorini (davideno.gorini@uniroma1.it)
Date: 2023
Description: This file contains the implementation for the ZeroLength6D class
Reference: Gorini, D.N. and Callisto, L. (2023): "A multiaxial inertial macroelement for deep foundations", Computers and Geotechnics, vol. 155, doi: doi.org/10.1016/j.compgeo.2022.105222
*/

#include <ZeroLength6D.h>
#include <Information.h>

#include <Domain.h>
#include <Node.h>
#include <Channel.h>
#include <FEM_ObjectBroker.h>
#include <NDMaterial.h>
#include <Renderer.h>
#include <ElementResponse.h>

#include <G3Globals.h>

#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <elementAPI.h>

Matrix ZeroLength6D::K12(12, 12);

Vector ZeroLength6D::P12(12);

Vector ZeroLength6D::v3(6);

//  Constructor:
//  responsible for allocating the necessary space needed by each object
//  and storing the tags of the ZeroLength6D end nodes.

static Node *theNodes[2];

void* OPS_ZeroLength6D()
{
	int ndm = OPS_GetNDM();
	int numdata = OPS_GetNumRemainingInputArgs();
	if (numdata < 4) {
		opserr << "WARNING too few arguments " <<
			"want - element zeroLength6D eleTag? iNode? jNode? " <<
			"NDTag? <1DTag?>" <<
			"<-orient x1? x2? x3? y1? y2? y3?>\n";

		return 0;
	}

	int idata[4];
	numdata = 4;
	if (OPS_GetIntInput(&numdata, idata) < 0) {
		opserr << "WARNING: failed to get integer data\n";
		return 0;
	}
	NDMaterial* nmat = OPS_getNDMaterial(idata[3]);
	if (nmat == 0) {
		opserr << "WARNING: NDMaterial " << idata[3] << " is not defined\n";
		return 0;
	}

	const char* type = OPS_GetString();
	Vector x(3); x(0) = 1.0; x(1) = 0.0; x(2) = 0.0;
	Vector y(3); y(0) = 0.0; y(1) = 1.0; y(2) = 0.0;
	if (strcmp(type, "orient") == 0) {
		if (OPS_GetNumRemainingInputArgs() < 6) {
			opserr << "WARNING: insufficient orient values\n";
			return 0;
		}
		numdata = 3;
		if (OPS_GetDoubleInput(&numdata, &x(0)) < 0) {
			opserr << "WARNING: invalid double input\n";
			return 0;
		}
		if (OPS_GetDoubleInput(&numdata, &y(0)) < 0) {
			opserr << "WARNING: invalid double input\n";
			return 0;
		}
	}

	return new ZeroLength6D(idata[0], ndm, idata[1], idata[2], x, y, *nmat);
}

ZeroLength6D::ZeroLength6D(int tag, int dim, int Nd1, int Nd2,
	const Vector& x, const Vector& yprime,
	NDMaterial &theNDmat) :
	Element(tag, ELE_TAG_ZeroLength6D),
	connectedExternalNodes(2),
	dimension(dim), numDOF(0),
	transformation(6, 6), A(0), v(0), K(0), P(0),
	end1Ptr(0), end2Ptr(0), theNDMaterial(0), order(0)
{
	// Obtain copy of Nd material model
	theNDMaterial = theNDmat.getCopy();

	if (theNDMaterial == 0) {
		opserr << "ZeroLength6D::zeroLength6D-- failed to get copy of NDMaterial\n";
		exit(-1);
	}
	// Get the material order
	order = theNDMaterial->getOrder();

	// Set up the transformation matrix of direction cosines
	this->setUp(Nd1, Nd2, x, yprime);
}

ZeroLength6D::ZeroLength6D() :
	Element(0, ELE_TAG_ZeroLength6D),
	connectedExternalNodes(2),
	dimension(0), numDOF(0),
	transformation(6, 6), A(0), v(0), K(0), P(0),
	end1Ptr(0), end2Ptr(0), theNDMaterial(0), order(0)
{

}

ZeroLength6D::~ZeroLength6D()
{
	// invoke the destructor on any objects created by the object
	// that the object still holds a pointer to

	if (theNDMaterial != 0)
		delete theNDMaterial;
	if (A != 0)
		delete A;
}

int
ZeroLength6D::getNumExternalNodes(void) const
{
	return 2;
}

Node **
ZeroLength6D::getNodePtrs(void)
{
	theNodes[0] = end1Ptr;
	theNodes[1] = end2Ptr;
	return theNodes;
}

const ID &
ZeroLength6D::getExternalNodes(void)
{
	return connectedExternalNodes;
}



int
ZeroLength6D::getNumDOF(void)
{
	return numDOF;
}

void
ZeroLength6D::setDomain(Domain *theDomain)
{
	// check Domain is not null - invoked when object removed from a domain
	if (theDomain == 0) {
		end1Ptr = 0;
		end2Ptr = 0;
		return;
	}

	// first set the node pointers
	int Nd1 = connectedExternalNodes(0);
	int Nd2 = connectedExternalNodes(1);
	end1Ptr = theDomain->getNode(Nd1);
	end2Ptr = theDomain->getNode(Nd2);

	// if can't find both - send a warning message
	if (end1Ptr == 0 || end2Ptr == 0) {
		if (end1Ptr == 0)
			opserr << "ZeroLength6D::setDomain()-- Nd1 does not exist in model\n";
		else
			opserr << "ZeroLength6D::setDomain -- Nd2 does not exist in model\n";

		end1Ptr = 0;
		end2Ptr = 0;
		return;
	}

	// now determine the number of dof and the dimension    
	int dofNd1 = end1Ptr->getNumberDOF();
	int dofNd2 = end2Ptr->getNumberDOF();

	// if differing dof at the ends - print a warning message
	if (dofNd1 != dofNd2) {
		opserr << "ZeroLength6D::setDomain -- nodes have differing dof's at end\n";
		end1Ptr = 0;
		end2Ptr = 0;
		return;
	}

	numDOF = 2 * dofNd1;

	if (numDOF != 6 && numDOF != 12) {
		opserr << "ZeroLength6D::setDomain -  element only works for 3 (2d) or 6 (3d) dof per node\n";
		end1Ptr = 0;
		end2Ptr = 0;
		return;
	}

	// Check that length is zero within tolerance
	const Vector &end1Crd = end1Ptr->getCrds();
	const Vector &end2Crd = end2Ptr->getCrds();
	const Vector     diff = end1Crd - end2Crd;
	double L = diff.Norm();
	double v1 = end1Crd.Norm();
	double v2 = end2Crd.Norm();
	double vm;
	double end1Crd0 = end1Crd(0);
	double end1Crd2 = end1Crd(2);
	vm = (v1<v2) ? v2 : v1;

	if (L > LENTOL*vm)
		opserr << "ZeroLength6D::setDomain -- Element has L=, which is greater than the tolerance\n";

	// call the base class method
	this->DomainComponent::setDomain(theDomain);

	// Set up the A matrix
	this->setTransformation();
}

int
ZeroLength6D::commitState()
{
	int err = 0;

	// Commit the NDMaterial
	err += theNDMaterial->commitState();

	return err;
}

int
ZeroLength6D::revertToLastCommit()
{
	int err = 0;

	// Revert the NDMaterial
	err += theNDMaterial->revertToLastCommit();

	return err;
}

int
ZeroLength6D::revertToStart()
{
	int err = 0;

	// Revert the NDMaterial to start
	err += theNDMaterial->revertToStart();

	return err;
}

const Matrix &
ZeroLength6D::getTangentStiff(void)
{
	// Compute material strains
	this->computeStrain();

	// Set trial strain for NDMaterial
	theNDMaterial->setTrialStrain(*v);

	// Get NDMaterial tangent, the element basic stiffness
	const Matrix &kb = theNDMaterial->getTangent();

	// Set some references to make the syntax nicer
	Matrix &stiff = *K;
	const Matrix &tran = *A; // A = transformation matrix
	double tran00e = tran(0, 0);
	double tran01e = tran(0, 1);
	double tran02e = tran(0, 2);
	double tran03e = tran(0, 3);
	double tran04e = tran(0, 4);
	double tran05e = tran(0, 5);
	double tran10e = tran(1, 0);
	double tran11e = tran(1, 1);
	double tran12e = tran(1, 2);
	double tran13e = tran(1, 3);
	double tran14e = tran(1, 4);
	double tran15e = tran(1, 5);
	double tran20e = tran(2, 0);
	double tran21e = tran(2, 1);
	double tran22e = tran(2, 2);
	double tran23e = tran(2, 3);
	double tran24e = tran(2, 4);
	double tran25e = tran(2, 5);
	double tran30e = tran(3, 0);
	double tran31e = tran(3, 1);
	double tran32e = tran(3, 2);
	double tran33e = tran(3, 3);
	double tran34e = tran(3, 4);
	double tran35e = tran(3, 5);
	double tran40e = tran(4, 0);
	double tran41e = tran(4, 1);
	double tran42e = tran(4, 2);
	double tran43e = tran(4, 3);
	double tran44e = tran(4, 4);
	double tran45e = tran(4, 5);
	double tran50e = tran(5, 0);
	double tran51e = tran(5, 1);
	double tran52e = tran(5, 2);
	double tran53e = tran(5, 3);
	double tran54e = tran(5, 4);
	double tran55e = tran(5, 5);

	stiff.Zero();

	double E;

	// Compute element stiffness ... K = A^*kb*A
	int k = 0;
	int l = 0;
	int i = 0;
	int j = 0;
	for (k = 0; k < order * 2; k++) {
		for (l = 0; l < order * 2; l++) {
			E = kb(k, l);
			for (i = 0; i < numDOF; i++)
				for (j = 0; j < i + 1; j++)
					stiff(i, j) += tran(k, i) * E * tran(l, j);
			double icc = i;
			double jcc = j;
			double kcc = k;
			double lcc = l;
			if (i == 0 && j == 0)
			{
				double stiff00 = stiff(i, j);
				double trank00 = tran(k, i);
				double tranl00 = tran(l, j);
				double E00 = E;
				double pau00 = 1;
			}
			if (i == 1 && j == 1)
			{
				double stiff11 = stiff(i, j);
				double trank11 = tran(k, i);
				double tranl11 = tran(l, j);
				double E11 = E;
				double pau11 = 1;
			}
			if (i == 2 && j == 2)
			{
				double stiff22 = stiff(i, j);
				double trank22 = tran(k, i);
				double tranl22 = tran(l, j);
				double E22 = E;
				double pau22 = 1;
			}
			if (i == 3 && j == 3)
			{
				double stiff33 = stiff(i, j);
				double trank33 = tran(k, i);
				double tranl33 = tran(l, j);
				double E33 = E;
				double pau33 = 1;
			}
			if (i == 4 && j == 4)
			{
				double stiff44 = stiff(i, j);
				double trank44 = tran(k, i);
				double tranl44 = tran(l, j);
				double E44 = E;
				double pau44 = 1;
			}
			if (i == 5 && j == 5)
			{
				double stiff55 = stiff(i, j);
				double trank55 = tran(k, i);
				double tranl55 = tran(l, j);
				double E55 = E;
				double pau55 = 1;
			}
		}
	}

	double tran00 = tran(0, 0);
	double tran11 = tran(1, 1);
	double tran22 = tran(2, 2);
	double tran33 = tran(3, 3);
	double tran44 = tran(4, 4);
	double tran55 = tran(5, 5);

	double stiff00 = stiff(0, 0);
	double stiff01 = stiff(0, 1);
	double stiff02 = stiff(0, 2);
	double stiff03 = stiff(0, 3);
	double stiff04 = stiff(0, 4);
	double stiff05 = stiff(0, 5);
	double stiff10 = stiff(1, 0);
	double stiff11 = stiff(1, 1);
	double stiff12 = stiff(1, 2);
	double stiff13 = stiff(1, 3);
	double stiff14 = stiff(1, 4);
	double stiff15 = stiff(1, 5);
	double stiff20 = stiff(2, 0);
	double stiff21 = stiff(2, 1);
	double stiff22 = stiff(2, 2);
	double stiff23 = stiff(2, 3);
	double stiff24 = stiff(2, 4);
	double stiff25 = stiff(2, 5);
	double stiff30 = stiff(3, 0);
	double stiff31 = stiff(3, 1);
	double stiff32 = stiff(3, 2);
	double stiff33 = stiff(3, 3);
	double stiff34 = stiff(3, 4);
	double stiff35 = stiff(3, 5);
	double stiff40 = stiff(4, 0);
	double stiff41 = stiff(4, 1);
	double stiff42 = stiff(4, 2);
	double stiff43 = stiff(4, 3);
	double stiff44 = stiff(4, 4);
	double stiff45 = stiff(4, 5);
	double stiff50 = stiff(5, 0);
	double stiff51 = stiff(5, 1);
	double stiff52 = stiff(5, 2);
	double stiff53 = stiff(5, 3);
	double stiff54 = stiff(5, 4);
	double stiff55 = stiff(5, 5);

	// Complete symmetric stiffness matrix
	for (int i = 0; i < numDOF; i++)
		for (int j = 0; j < i; j++)
			stiff(j, i) = stiff(i, j);

	return stiff;
}

const Matrix &
ZeroLength6D::getInitialStiff(void)
{
	// Get NDMaterial tangent, the element basic stiffness
	const Matrix &kb = theNDMaterial->getInitialTangent();

	// Set some references to make the syntax nicer
	Matrix &stiff = *K;
	const Matrix &tran = *A;

	stiff.Zero();

	double E;

	// Compute element stiffness ... K = A^*kb*A
	for (int k = 0; k < order * 2; k++) {
		for (int l = 0; l < order * 2; l++) {
			E = kb(k, l);
			for (int i = 0; i < numDOF; i++)
				for (int j = 0; j < i + 1; j++)
					stiff(i, j) += tran(k, i) * E * tran(l, j);
		}
	}

	// Complete symmetric stiffness matrix
	for (int i = 0; i < numDOF; i++)
		for (int j = 0; j < i; j++)
			stiff(j, i) = stiff(i, j);

	return stiff;
}

const Matrix &
ZeroLength6D::getDamp(void)
{
	// Return zero damping
	K->Zero();

	return *K;
}

const Matrix &
ZeroLength6D::getMass(void)
{
	// Return zero mass
	K->Zero();

	return *K;
}

void
ZeroLength6D::zeroLoad(void)
{
	// does nothing now
}

int
ZeroLength6D::addLoad(ElementalLoad *theLoad, double loadFactor)
{
	opserr << "ZeroLength::addLoad - load type unknown for ZeroLength6D\n";

	return -1;
}

int
ZeroLength6D::addInertiaLoadToUnbalance(const Vector &accel)
{
	// does nothing as element has no mass yet!
	return 0;
}

const Vector &
ZeroLength6D::getResistingForce()
{
	// Compute material strains
	this->computeStrain();

	// Set trial strain for NDMaterial
	theNDMaterial->setTrialStrain(*v);

	// Get NDMaterial stress, the element basic force
	const Vector &q = theNDMaterial->getStress();

	// Set some references to make the syntax nicer
	Vector &force = *P;
	const Matrix &tran = *A;

	force.Zero();

	double s;

	// Compute element resisting force ... P = A^*q
	for (int k = 0; k < order * 2; k++) {
		s = q(k);
		for (int i = 0; i < numDOF; i++)
			force(i) += tran(k, i) * s;
	}

	return force;
}

const Vector &
ZeroLength6D::getResistingForceIncInertia()
{
	// There is no mass, so return
	return this->getResistingForce();
}

int
ZeroLength6D::sendSelf(int commitTag, Channel &theChannel)
{
	int res = 0;

	// note: we don't check for dataTag == 0 for Element
	// objects as that is taken care of in a commit by the Domain
	// object - don't want to have to do the check if sending data
	int dataTag = this->getDbTag();

	// ZeroLength6D packs its data into an ID and sends this to theChannel
	// along with its dbTag and the commitTag passed in the arguments

	static ID idData(8);

	idData(0) = this->getTag();
	idData(1) = dimension;
	idData(2) = numDOF;
	idData(3) = order;
	idData(4) = connectedExternalNodes(0);
	idData(5) = connectedExternalNodes(1);
	idData(6) = theNDMaterial->getClassTag();

	int dbTag = theNDMaterial->getDbTag();
	if (dbTag == 0) {
		dbTag = theChannel.getDbTag();
		if (dbTag != 0)
			theNDMaterial->setDbTag(dbTag);
	}
	idData(7) = dbTag;

	res += theChannel.sendID(dataTag, commitTag, idData);
	if (res < 0) {
		opserr << "ZeroLength6D::sendSelf() -- failed to send ID data\n";
		return res;
	}

	res += theChannel.sendMatrix(dataTag, commitTag, transformation);
	if (res < 0) {
		opserr << "ZeroLength6D::sendSelf -- failed to send transformation Matrix\n";
		return res;
	}

	// Send the NDMaterial
	res += theNDMaterial->sendSelf(commitTag, theChannel);
	if (res < 0) {
		opserr << "ZeroLength6D::  -- failed to send NDMaterial\n";
		return res;
	}

	return res;
}

int
ZeroLength6D::recvSelf(int commitTag, Channel &theChannel, FEM_ObjectBroker &theBroker)
{
	int res = 0;

	int dataTag = this->getDbTag();

	// ZeroLength6D creates an ID, receives the ID and then sets the 
	// internal data with the data in the ID

	static ID idData(8);

	res += theChannel.recvID(dataTag, commitTag, idData);
	if (res < 0) {
		opserr << "ZeroLength6D::recvSelf -- failed to receive ID data\n";
		return res;
	}

	res += theChannel.recvMatrix(dataTag, commitTag, transformation);
	if (res < 0) {
		opserr << "zeroLength6D::revbSelf -- failed to receive transformation Matrix\n";
		return res;
	}

	this->setTag(idData(0));
	dimension = idData(1);
	numDOF = idData(2);
	connectedExternalNodes(0) = idData(4);
	connectedExternalNodes(1) = idData(5);

	if (order != idData(3)) {

		order = idData(3);

		// Allocate transformation matrix
		if (A != 0)
			delete A;

		A = new Matrix(order * 2, numDOF);

		K = &K12;
		P = &P12;
		v = &v3;
	}

	int classTag = idData(6);

	// If null, get a new one from the broker
	if (theNDMaterial == 0)
		theNDMaterial = theBroker.getNewNDMaterial(classTag);

	// If wrong type, get a new one from the broker
	if (theNDMaterial->getClassTag() != classTag) {
		delete theNDMaterial;
		theNDMaterial = theBroker.getNewNDMaterial(classTag);
	}

	// Check if either allocation failed from broker
	if (theNDMaterial == 0) {
		opserr << "ZeroLength6D::  -- failed to allocate new NDMaterial\n";
		return -1;
	}

	// Receive the NDMaterial
	theNDMaterial->setDbTag(idData(7));
	res += theNDMaterial->recvSelf(commitTag, theChannel, theBroker);
	if (res < 0) {
		opserr << "ZeroLength6D::  -- failed to receive NDMaterial\n";
		return res;
	}

	return res;
}

int
ZeroLength6D::displaySelf(Renderer &theViewer, int displayMode, float fact, const char **modes, int numMode)
{
	// ensure setDomain() worked
	if (end1Ptr == 0 || end2Ptr == 0)
		return 0;

	// first determine the two end points of the ZeroLength6D based on
	// the display factor (a measure of the distorted image)
	// store this information in 2 3d vectors v1 and v2
	const Vector &end1Crd = end1Ptr->getCrds();
	const Vector &end2Crd = end2Ptr->getCrds();
	const Vector &end1Disp = end1Ptr->getDisp();
	const Vector &end2Disp = end2Ptr->getDisp();

	if (displayMode == 1 || displayMode == 2) {
		static Vector v1(6);
		static Vector v2(6);
		for (int i = 0; i < dimension * 2; i++) {
			v1(i) = end1Crd(i) + end1Disp(i)*fact;
			v2(i) = end2Crd(i) + end2Disp(i)*fact;
		}

		return theViewer.drawLine(v1, v2, 0.0, 0.0);
	}

	return 0;
}

void
ZeroLength6D::Print(OPS_Stream &s, int flag)
{
	if (flag == OPS_PRINT_CURRENTSTATE) {
		s << "ZeroLength6D, tag: " << this->getTag() << endln;
		s << "\tConnected Nodes: " << connectedExternalNodes << endln;
		s << "\tNDMaterial, tag: " << theNDMaterial->getTag() << endln;
	}

	if (flag == OPS_PRINT_PRINTMODEL_JSON) {
		s << "\t\t\t{";
		s << "\"name\": \"" << this->getTag() << "\", ";
		s << "\"type\": \"ZeroLength6D\", ";
		s << "\"nodes\": [\"" << connectedExternalNodes(0) << "\", \"" << connectedExternalNodes(1) << "\"], ";
		s << "\"ndMaterial\": \"" << theNDMaterial->getTag() << "\", ";
		s << "\"transMatrix\": [[";
		for (int i = 0; i < 6; i++) {
			for (int j = 0; j < 6; j++) {
				if (j < 5)
					s << transformation(i, j) << ", ";
				else if (j == 5 && i < 5)
					s << transformation(i, j) << "], [";
				else if (j == 5 && i == 5)
					s << transformation(i, j) << "]]}";
			}
		}
	}
}

Response*
ZeroLength6D::setResponse(const char **argv, int argc, OPS_Stream &output)
{
	Response *theResponse = 0;

	output.tag("ElementOutput");
	output.attr("eleType", "ZeroLength");
	output.attr("eleTag", this->getTag());
	output.attr("node1", connectedExternalNodes[0]);
	output.attr("node2", connectedExternalNodes[1]);

	char outputData[10];

	if ((strcmp(argv[0], "force") == 0) || (strcmp(argv[0], "forces") == 0)
		|| (strcmp(argv[0], "globalForces") == 0) || (strcmp(argv[0], "globalforces") == 0)) {

		char outputData[10];
		int numDOFperNode = numDOF / 2;
		for (int i = 0; i<numDOFperNode; i++) {
			sprintf(outputData, "P1_%d", i + 1);
			output.tag("ResponseType", outputData);
		}
		for (int j = 0; j<numDOFperNode; j++) {
			sprintf(outputData, "P2_%d", j + 1);
			output.tag("ResponseType", outputData);
		}
		theResponse = new ElementResponse(this, 1, *P);

	}
	else if (strcmp(argv[0], "basicForce") == 0 || strcmp(argv[0], "basicForces") == 0) {

		for (int i = 0; i<order * 2; i++) {
			sprintf(outputData, "P%d", i + 1);
			output.tag("ResponseType", outputData);
		}
		theResponse = new ElementResponse(this, 2, Vector(order * 2));

	}
	else if (strcmp(argv[0], "defo") == 0 || strcmp(argv[0], "deformations") == 0 ||
		strcmp(argv[0], "deformation") == 0) {

		for (int i = 0; i<order * 2; i++) {
			sprintf(outputData, "e%d", i + 1);
			output.tag("ResponseType", outputData);
		}
		theResponse = new ElementResponse(this, 3, Vector(order * 2));

		// a material quantity
	}
	else if (strcmp(argv[0], "material") == 0) {
		// See if NDMaterial can handle request ...
		theResponse = theNDMaterial->setResponse(&argv[1], argc - 1, output);
	}

	output.endTag();
	return theResponse;
}

int
ZeroLength6D::getResponse(int responseID, Information &eleInfo)
{
	switch (responseID) {
	case 1:
		return eleInfo.setVector(this->getResistingForce());

	case 2:
		if (eleInfo.theVector != 0) {
			const Vector &tmp = theNDMaterial->getStress();
			Vector &force = *(eleInfo.theVector);
			for (int i = 0; i < order * 2; i++)
				force(i) = tmp(i);
		}
		return 0;

	case 3:
		if (eleInfo.theVector != 0) {
			this->computeStrain();
			const Vector &tmp = *v;	// NDMaterial strains
			Vector &def = *(eleInfo.theVector);
			for (int i = 0; i < order * 2; i++)
				def(i) = tmp(i);
		}
		return 0;

	default:
		return -1;
	}
}

// Private methods


// Establish the external nodes and set up the transformation matrix
// for orientation
void
ZeroLength6D::setUp(int Nd1, int Nd2, const Vector &x, const Vector &yp)
{
	// ensure the connectedExternalNode ID is of correct size & set values
	if (connectedExternalNodes.Size() != 2) {
		opserr << "ZeroLength6D::setUp -- failed to create an ID of correct size\n";
		exit(-1);
	}

	connectedExternalNodes(0) = Nd1;
	connectedExternalNodes(1) = Nd2;

	// check that vectors for orientation are correct size
	if (x.Size() != 3 || yp.Size() != 3) {
		opserr << "ZeroLength6D -- incorrect dimension of orientation vectors\n";
		exit(-1);
	}
	// establish orientation of element for the tranformation matrix
	// z = x cross yp
	static Vector z(3);
	z(0) = x(1)*yp(2) - x(2)*yp(1);
	z(1) = x(2)*yp(0) - x(0)*yp(2);
	z(2) = x(0)*yp(1) - x(1)*yp(0);

	// y = z cross x
	static Vector y(3);
	y(0) = z(1)*x(2) - z(2)*x(1);
	y(1) = z(2)*x(0) - z(0)*x(2);
	y(2) = z(0)*x(1) - z(1)*x(0);

	// compute length(norm) of vectors
	double xn = x.Norm();
	double yn = y.Norm();
	double zn = z.Norm();

	// check valid x and y vectors, i.e. not parallel and of zero length
	if (xn == 0 || yn == 0 || zn == 0) {
		opserr << "ZeroLength6D::setUP -- invalid vectors to constructor\n";
		exit(-1);
	}

	// create transformation matrix of direction cosines
	for (int i = 0; i < 6; i++) {
		transformation(0, i) = x(i) / xn;
		transformation(1, i) = y(i) / yn;
		transformation(2, i) = z(i) / zn;
	}
	transformation(3, 3) = 1.0;
	transformation(4, 4) = 1.0;
	transformation(5, 5) = 1.0;
}

// Set basic deformation-displacement transformation matrix for the materials
void
ZeroLength6D::setTransformation(void)
{
	// Allocate transformation matrix
	if (A != 0)
		delete A;

	A = (1 == 0) ? new Matrix(order * 2, numDOF) : new Matrix(order * 2 + 1, numDOF);

	if (A == 0) {
		opserr << "ZeroLength6D::setTransformation -- failed to allocate transformation Matrix\n";
		exit(-1);
	}
	K = &K12;
	P = &P12;
	v = &v3;

	// Set a reference to make the syntax nicer
	Matrix &tran = *A;

	// Loop over the NDMaterial order
	for (int i = 0; i < order * 2; i++) {

		tran(i, 6) = transformation(i, 0);
		tran(i, 7) = transformation(i, 1);
		tran(i, 8) = transformation(i, 2);
		tran(i, 9) = transformation(i, 3);
		tran(i, 10) = transformation(i, 4);
		tran(i, 11) = transformation(i, 5);

		// Fill in first half of transformation matrix with negative sign
		for (int j = 0; j < numDOF / 2; j++)
			tran(i, j) = -tran(i, j + numDOF / 2);
	}
}

void
ZeroLength6D::computeStrain(void)
{
	// Get nodal displacements
	const Vector &u1 = end1Ptr->getTrialDisp();
	const Vector &u2 = end2Ptr->getTrialDisp();

	// Compute differential displacements
	const Vector diff = u2 - u1;
	// Set some references to make the syntax nicer
	Vector &def = *v;
	const Matrix &tran = *A;

	def.Zero();

	// Compute element basic deformations ... v = A*(u2-u1)
	for (int i = 0; i < order * 2; i++)
		for (int j = 0; j < numDOF / 2; j++)
			def(i) += -diff(j)*tran(i, j);

}
