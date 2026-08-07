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

// Massimo Petracca - ASDEA Software, Italy
//
// The piecewise-linear hardening law shared by the ASD material family
// (ASDConcrete1D, ASDConcrete3D and the models built on top of them).
//
// It used to be duplicated, verbatim, as a nested class of each material.
// It lives here so that a fix, a new accessor or a change of convention
// cannot apply to one model and not to the other: divergent copies of this
// exact code are how the error metric of the family became non-uniform in
// the first place.
//
// THE PRISTINE LAW IS CARRIED BY THE LAW ITSELF. Regularization rescales the
// softening branch to the characteristic length of the parent element, so
// every law needs the un-regularized curve to go back to before rescaling
// again. That backup used to live in a singleton keyed by the MATERIAL TAG,
// which is wrong in two ways: the tags of two different material classes
// collide, and the singleton outlives a 'wipe', so
//
//     build a model (tag 1, law A) -> wipe -> build another (tag 1, law B)
//
// silently recovered law A inside the second model, at the first
// regularization. Here each law holds a shared_ptr to its own immutable
// pristine copy, so the backup cannot be mixed up with anyone else's, and it
// is shared by every copy of the material - one array of points for the
// millions of gauss points that use the same law, which is what the storage
// was meant to do in the first place and no longer did.

#ifndef ASDHardeningLaw_h
#define ASDHardeningLaw_h

#include <Vector.h>
#include <vector>
#include <memory>
#include <map>
#include <cstddef>

/**
The hardening law type (which side of the response it describes)
*/
enum class ASDHardeningLawType {
	Tension = 0,
	Compression
};

/**
The component of a hardening law point (for output)
*/
enum class ASDHardeningLawPointComponent {
	TotalStrain = 0,
	EffectiveStress,
	NominalStress
};

/**
A point in the hardening law
*/
struct ASDHardeningLawPoint {
	double x = 0.0; // total strain
	double y = 0.0; // final backbone stress (plasticity & damage)
	double d = 0.0; // damage variable
	double q = 0.0; // backbone stress in the effective-stress space (plasticity only)
	ASDHardeningLawPoint() = default;
	ASDHardeningLawPoint(double _x, double _y, double _d, double _q)
		: x(_x), y(_y), d(_d), q(_q) {}
	inline double totalStrain()const { return x; }
	inline double plasticStrain(double E)const { return x - q / E; }
	inline double stress()const { return y; }
	inline double effectiveStress()const { return q; }
	inline double elasticStress(double E) const { return E * x; }
	inline double crackingDamage()const { return d; }
};

/**
The hardening law
*/
class ASDHardeningLaw {
public:
	// default constructor
	ASDHardeningLaw() = default;
	// full constructor
	ASDHardeningLaw(
		int tag, ASDHardeningLawType type,
		double E,
		const std::vector<double>& x,
		const std::vector<double>& y,
		const std::vector<double>& d);

	// regularizes the hardening curve according to the 'lch'
	// characteristic length of the parent element, and the 'lch_ref'
	// characteristic length of the input hardening curve
	void regularize(double lch, double lch_ref);
	// nullify previous regularization
	void deRegularize();
	// evaluate the hardening law at a certain strain
	ASDHardeningLawPoint evaluateAt(double x) const;
	// get max stress value
	double computeMaxStress() const;
	// serialization
	int serializationDataSize() const;
	void serialize(Vector& data, int& pos);
	void deserialize(Vector& data, int& pos);
	// properties
	inline bool isValid() const { return m_valid; }
	inline int tag()const { return m_tag; }
	inline ASDHardeningLawType type()const { return m_type; }
	inline std::size_t uid()const { return m_uid; }
	inline const std::vector<ASDHardeningLawPoint>& points()const { return m_points; }
	inline double strainTolerance()const { return m_xtolerance; }
	inline double stressTolerance()const { return m_ytolerance; }
	inline bool hasStrainSoftening()const { return m_fracture_energy_is_bounded; }
	inline double strainAtOnsetOfCrack()const {
		if (m_fracture_energy_is_bounded && m_softening_begin < m_points.size())
			return m_points[m_softening_begin].totalStrain();
		return 0.0;
	}

private:
	// adjusts the points
	void adjust();
	// computes the fracture energy of the hardening curve
	void computeFractureEnergy();

private:
	// tag (same as the parent material's tag). Kept for output and for the
	// serialized layout only: it is NOT an identity, see the note on top.
	int m_tag = 0;
	// type
	ASDHardeningLawType m_type = ASDHardeningLawType::Tension;
	// the hardening points of the backbone curve in total-strain
	std::vector<ASDHardeningLawPoint> m_points;
	// the fracture energy (computed only if there is softening)
	double m_fracture_energy = 0.0;
	// true if there is softening, false otherwise
	bool m_fracture_energy_is_bounded = false;
	// the location of the point in m_points where the softening begins
	std::size_t m_softening_begin = 0;
	// the location of the point in m_points where the softening ends
	std::size_t m_softening_end = 0;
	// false if the initialization (in constructor) fails, true otherwise
	bool m_valid = false;
	// tolerances
	double m_xtolerance = 1.0e-12;
	double m_ytolerance = 1.0e-12;
	// process-wide unique identity of this law, assigned once by the full
	// constructor and NEVER reused, so that a law defined after a 'wipe'
	// cannot be mistaken for one defined before it. 0 = no identity
	// (default-constructed or invalid law).
	std::size_t m_uid = 0;
	// the un-regularized law, shared by every copy of this one. Null only for
	// invalid laws. See deRegularize().
	std::shared_ptr<const ASDHardeningLaw> m_pristine;
};

/**
Storage for the pristine laws that arrive over a Channel.

On the build path no storage is needed: the full constructor makes the
pristine copy and every material copy shares it through the shared_ptr. On the
RECEIVE path, though, each material copy is deserialized on its own - under
OpenSeesSP one per gauss point - and without a way to recognize that they carry
the same law each copy would allocate its own array of points. The sender's uid
travels with the law and is the key here: it is unique in the sending process
and never reused, so sharing by it is always correct.

Values are weak, so entries die with the last material that holds the law and
this singleton never needs to be cleared - which is what made the old one
outlive a 'wipe' and hand back a stale law.
*/
class ASDHardeningLawStorage {
public:
	using PointerType = std::shared_ptr<const ASDHardeningLaw>;
	using MapType = std::map<std::size_t, std::weak_ptr<const ASDHardeningLaw>>;

public:
	ASDHardeningLawStorage() = default;
	ASDHardeningLawStorage(const ASDHardeningLawStorage&) = delete;
	ASDHardeningLawStorage& operator = (const ASDHardeningLawStorage&) = delete;

public:
	// access to this singleton
	static ASDHardeningLawStorage& instance();
	// a new identity. Monotonic, never reset, never reused
	std::size_t generateUID();
	// share (or create and share) the pristine law received under the
	// sender's uid
	PointerType internReceived(std::size_t uid, const ASDHardeningLaw& law);

private:
	std::size_t m_counter = 0;
	MapType m_received;
};

#endif // ASDHardeningLaw_h
