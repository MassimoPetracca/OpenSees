/* 
Developed by: Davide Noè Gorini (davideno.gorini@uniroma1.it)
Date: 2023
Description: This file contains the implementation for the ZeroLength6D class
Reference: Gorini, D.N. and Callisto, L. (2023): "A multiaxial inertial macroelement for deep foundations", Computers and Geotechnics, vol. 155, doi: doi.org/10.1016/j.compgeo.2022.105222
*/

#ifndef ZeroLength6D_h
#define ZeroLength6D_h

#include <Element.h>
#include <Matrix.h>

// Tolerance for zero length of element
#define	LENTOL 1.0e-6

class Node;
class Channel;
class NDMaterial;
class Response;

class ZeroLength6D : public Element
{
public:

	// Constructor for a single Nd material model of order 2 or 3
	ZeroLength6D(int tag,
		int dimension,
		int Nd1, int Nd2,
		const Vector& x,
		const Vector& yprime,
		NDMaterial& theNDMaterial);

	ZeroLength6D();
	~ZeroLength6D();

	// public methods to obtain inforrmation about dof & connectivity    
	int getNumExternalNodes(void) const;
	const ID &getExternalNodes(void);
	Node **getNodePtrs(void);

	int getNumDOF(void);
	void setDomain(Domain *theDomain);

	// public methods to set the state of the element    
	int commitState(void);
	int revertToLastCommit(void);
	int revertToStart(void);

	// public methods to obtain stiffness, mass, damping and residual information    
	const Matrix &getTangentStiff(void);
	const Matrix &getInitialStiff(void);
	const Matrix &getDamp(void);
	const Matrix &getMass(void);

	void zeroLoad(void);
	int addLoad(ElementalLoad *theLoad, double loadFactor);
	int addInertiaLoadToUnbalance(const Vector &accel);

	const Vector &getResistingForce(void);
	const Vector &getResistingForceIncInertia(void);

	// public methods for element output
	int sendSelf(int commitTag, Channel &theChannel);
	int recvSelf(int commitTag, Channel &theChannel, FEM_ObjectBroker &theBroker);
	int displaySelf(Renderer &, int mode, float fact, const char **displayModes = 0, int numModes = 0);
	void Print(OPS_Stream &s, int flag = 0);

	Response *setResponse(const char **argv, int argc, OPS_Stream &theHandler);
	int getResponse(int responseID, Information &eleInformation);

protected:

private:
	// private methods
	void setUp(int Nd1, int Nd2, const Vector& x, const Vector& y);
	void setTransformation(void);
	void computeStrain(void);

	// private attributes - a copy for each object of the class
	ID  connectedExternalNodes;         // contains the tags of the end nodes
	int dimension;                      // 3 dimensions
	int numDOF;	                        // number of dof for ZeroLength6D
	Matrix transformation;		// transformation matrix for orientation

	Matrix *A;	// Transformation matrix ... e = A*(u2-u1)
	Vector *v;	// NDMaterial strain vector, the element basic deformations

	Matrix *K;	// Element stiffness matrix
	Vector *P;	// Element force vector

	Node *end1Ptr;      		// pointer to the end1 node object
	Node *end2Ptr;      		// pointer to the end2 node object	

	NDMaterial *theNDMaterial;	// Pointer to NDMaterial object
	int order;		// Order of the NDMaterial (2 or 3)

					// Class wide matrices for return
	static Matrix K12;

	// Class wide vectors for return
	static Vector P12;

	// Class wide vectors for storing NDMaterial strains
	static Vector v3;
};

#endif
