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

// $Revision: 1.00 $
// $Date: 2026/04/15 22:51:21 $

// Original implementation: Stefano Ercolessi & Massimo Petracca (ASDEA)

// Implementation of a 8 node solid hexa element local coordinate system 

#ifndef ASDSolidHexLocalCoordinateSystem_h
#define ASDSolidHexLocalCoordinateSystem_h

//This Class Represent the local coordinate system of any 8 node solid element in 3D space

#include <ASDMath.h>
#include <array>
#include <vector>

class ASDSolidHexLocalCoordinateSystem
{
public:
	typedef ASDVector3<double> Vector3Type;
	typedef ASDQuaternion<double> QuaternionType;
	typedef std::vector<Vector3Type> Vector3ContainerType;
	typedef Matrix MatrixType;

public:
	// default constructor
	ASDSolidHexLocalCoordinateSystem()
		: m_P(8),
		  m_origin(Vector3Type(0, 0, 0)),
		m_Rtilde(3, 3),
		m_R(3, 3)
	{
		m_Rtilde.Zero();
		m_R.Zero();
		for (int i = 0; i < 8; i++) {
			m_P[i] = Vector3Type(0, 0, 0);
		}
	}
	// constructor 
	ASDSolidHexLocalCoordinateSystem(const Vector3Type& P1global,
		const Vector3Type& P2global,
		const Vector3Type& P3global,
		const Vector3Type& P4global,
		const Vector3Type& P5global,
		const Vector3Type& P6global,
		const Vector3Type& P7global,
		const Vector3Type& P8global,
		const MatrixType* matrixPtr = nullptr) :
		m_P(8),
		m_Rtilde(3, 3),
		m_R(3,3)
		{
		// compute the origin as the average of the 8 nodes coordinates
		m_origin = (P1global + P2global + P3global + P4global + P5global + P6global + P7global + P8global) * 0.125;
		// compute the e1 e2 and e3 unit vectors of the local coordinate system
		Vector3Type e1 = 0.25 * (P2global + P3global + P6global + P7global - P1global - P4global - P5global - P8global);
		//opserr << " e1: " << e1 << endln;
		e1.normalize();
		//opserr << " e1 normalized: " << e1 << endln;
		Vector3Type e2tmp = 0.25 * (P3global + P4global + P7global + P8global - P1global - P2global - P5global - P6global);
		//opserr << " e2tmp: " << e2tmp;
		e2tmp.normalize();
		//opserr << " e2tmp normalized: " << e2tmp << endln;
		// compute e3 as cross product of e1 and e2tmp
		Vector3Type e3 = e1.cross(e2tmp);
		//opserr << " e3: " << e3 << endln;
		e3.normalize();
		//opserr << " e3 normalized: " << e3 << endln;
		// recompute e2 as cross product of e3 and e1 to ensure orthogonality
		Vector3Type e2 = e3.cross(e1);
		//opserr << " e2: " << e2 << endln;
		e2.normalize();
		//opserr << " e2 normalized: " << e2 << endln;

		//opserr << "versors system: " << " e1 = " << e1 << "e2 = " << e2 << "e3 = " << e3;
		// set the orientation matrix
		for (int ii = 0; ii < 3; ii++) {
			m_Rtilde(0, ii) = e1(ii);
			m_Rtilde(1, ii) = e2(ii);
			m_Rtilde(2, ii) = e3(ii);
		}

		//opserr << "R_tilde " << m_Rtilde;

		if (matrixPtr != nullptr) {
#if defined(USE_SIMPLE_RFRAME) && (USE_SIMPLE_RFRAME == 1)
			// m_R = matrixPtr^T  (assumendo Rtilde_init = I).
			// matrixPtr e' la rotazione "active" (polar/Kabsch di best-fit) che
			// porta initial -> current. Il rotatore global->CR-local e' la sua
			// trasposta. Si IGNORA m_Rtilde nella composizione finale per
			// evitare il "raddoppio" della rotazione (vedi diagnosi 2026-04-24).
			for (int i = 0; i < 3; i++)
				for (int j = 0; j < 3; j++)
					m_R(i, j) = (*matrixPtr)(j, i);
#else
			// if a ptr to a matrix is provided compute the R = Rrtilde * Rtilde^T
			// temporary metrix to store the result
			m_R.addMatrixTransposeProduct(0.0, *matrixPtr, m_Rtilde, 1.0);
#endif

		} else {
			// otherwise R = I
			m_R(0, 0) = 1.0;
			m_R(1, 1) = 1.0;
			m_R(2, 2) = 1.0;
		}

		//opserr << "Rcr" << m_orientation;
		//opserr << "Origin" << m_origin;
		// transform global coordinates to the local coordinate system

		for (int ii = 0; ii < 3; ii++) {
			m_P[0](ii) = m_R(ii, 0) * (P1global(0) - m_origin(0)) + m_R(ii, 1) * (P1global(1) - m_origin(1)) + m_R(ii, 2) * (P1global(2) - m_origin(2));
			m_P[1](ii) = m_R(ii, 0) * (P2global(0) - m_origin(0)) + m_R(ii, 1) * (P2global(1) - m_origin(1)) + m_R(ii, 2) * (P2global(2) - m_origin(2));
			m_P[2](ii) = m_R(ii, 0) * (P3global(0) - m_origin(0)) + m_R(ii, 1) * (P3global(1) - m_origin(1)) + m_R(ii, 2) * (P3global(2) - m_origin(2));
			m_P[3](ii) = m_R(ii, 0) * (P4global(0) - m_origin(0)) + m_R(ii, 1) * (P4global(1) - m_origin(1)) + m_R(ii, 2) * (P4global(2) - m_origin(2));
			m_P[4](ii) = m_R(ii, 0) * (P5global(0) - m_origin(0)) + m_R(ii, 1) * (P5global(1) - m_origin(1)) + m_R(ii, 2) * (P5global(2) - m_origin(2));
			m_P[5](ii) = m_R(ii, 0) * (P6global(0) - m_origin(0)) + m_R(ii, 1) * (P6global(1) - m_origin(1)) + m_R(ii, 2) * (P6global(2) - m_origin(2));
			m_P[6](ii) = m_R(ii, 0) * (P7global(0) - m_origin(0)) + m_R(ii, 1) * (P7global(1) - m_origin(1)) + m_R(ii, 2) * (P7global(2) - m_origin(2));
			m_P[7](ii) = m_R(ii, 0) * (P8global(0) - m_origin(0)) + m_R(ii, 1) * (P8global(1) - m_origin(1)) + m_R(ii, 2) * (P8global(2) - m_origin(2));

		}

		//opserr << "Rcr" << m_orientation << "\n";
		//opserr << "Origin" << m_origin << "\n";
		//opserr << "P1local" << m_P[0] << "\n";
		//opserr << "P2local" << m_P[1] << "\n";
		//opserr << "P3local" << m_P[2] << "\n";
		//opserr << "P4local" << m_P[3] << "\n";
		//opserr << "P5local" << m_P[4] << "\n";
		//opserr << "P6local" << m_P[5] << "\n";
		//opserr << "P7local" << m_P[6] << "\n";
		//opserr << "P8local" << m_P[7] << "\n";

	}

public:

	inline const Vector3ContainerType& Nodes() const{ return m_P;}

	inline const Vector3Type& P1()const { return m_P[0]; }
	inline const Vector3Type& P2()const { return m_P[1]; }
	inline const Vector3Type& P3()const { return m_P[2]; }
	inline const Vector3Type& P4()const { return m_P[3]; }
	inline const Vector3Type& P5()const { return m_P[4]; }
	inline const Vector3Type& P6()const { return m_P[5]; }
	inline const Vector3Type& P7()const { return m_P[6]; }
	inline const Vector3Type& P8()const { return m_P[7]; }

	inline const Vector3Type& Origin()const { return m_origin; }

	inline double X1()const { return m_P[0][0]; }
	inline double X2()const { return m_P[1][0]; }
	inline double X3()const { return m_P[2][0]; }
	inline double X4()const { return m_P[3][0]; }
	inline double X5()const { return m_P[4][0]; }
	inline double X6()const { return m_P[5][0]; }
	inline double X7()const { return m_P[6][0]; }
	inline double X8()const { return m_P[7][0]; }

	inline double Y1()const { return m_P[0][1]; }
	inline double Y2()const { return m_P[1][1]; }
	inline double Y3()const { return m_P[2][1]; }
	inline double Y4()const { return m_P[3][1]; }
	inline double Y5()const { return m_P[4][1]; }
	inline double Y6()const { return m_P[5][1]; }
	inline double Y7()const { return m_P[6][1]; }
	inline double Y8()const { return m_P[7][1]; }

	inline double Z1()const { return m_P[0][2]; }
	inline double Z2()const { return m_P[1][2]; }
	inline double Z3()const { return m_P[2][2]; }
	inline double Z4()const { return m_P[3][2]; }
	inline double Z5()const { return m_P[4][2]; }
	inline double Z6()const { return m_P[5][2]; }
	inline double Z7()const { return m_P[6][2]; }
	inline double Z8()const { return m_P[7][2]; }

	inline double X(size_t i)const { return m_P[i][0]; }
	inline double Y(size_t i)const { return m_P[i][1]; }
	inline double Z(size_t i)const { return m_P[i][2]; }

	inline const MatrixType& getRotationMatrix()const { return m_R; }
	inline const MatrixType& Orientation()const { return m_Rtilde; }

	inline Vector3Type Vx()const { return Vector3Type(m_R(0, 0), m_R(0, 1), m_R(0, 2)); }
	inline Vector3Type Vy()const { return Vector3Type(m_R(1, 0), m_R(1, 1), m_R(1, 2)); }
	inline Vector3Type Vz()const { return Vector3Type(m_R(2, 0), m_R(2, 1), m_R(2, 2)); }

	inline const Vector3Type& P(std::size_t i) const { return m_P[i]; }

	inline void ComputeTotalRotationMatrix(MatrixType& R)const {
		constexpr size_t mat_size = 24;
		if (R.noRows() != mat_size || R.noCols() != mat_size) {
			R.resize(mat_size, mat_size);
		}
		R.Zero();
		for (size_t k = 0; k < 8; k++) {
			size_t i = k * 3;
			R(i, i) = m_R(0, 0);   R(i, i + 1) = m_R(0, 1);   R(i, i + 2) = m_R(0, 2);
			R(i + 1, i) = m_R(1, 0);   R(i + 1, i + 1) = m_R(1, 1);   R(i + 1, i + 2) = m_R(1, 2);
			R(i + 2, i) = m_R(2, 0);   R(i + 2, i + 1) = m_R(2, 1);   R(i + 2, i + 2) = m_R(2, 2);
		}

	}


private:
	// define m_P a container for the node coordinates
	Vector3ContainerType m_P;
	// define m_origin the origin of the local coordinate system
	Vector3Type m_origin;
	// define m_orientation the orientation of the local coordinate system with respect to the global one
	// each row representas the vesor of the local coordinate system in the global coordinates
	Matrix m_Rtilde;
	Matrix m_R;


};

#endif