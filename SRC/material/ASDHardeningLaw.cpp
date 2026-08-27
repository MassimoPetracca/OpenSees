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
// Implementation of the hardening law shared by the ASD material family.
// See ASDHardeningLaw.h for why it is here and not nested in a material.

#include <ASDHardeningLaw.h>
#include <OPS_Globals.h>
#include <cmath>
#include <algorithm>
#include <limits>

ASDHardeningLaw::ASDHardeningLaw(
	int tag, ASDHardeningLawType type,
	double E,
	const std::vector<double>& x, const std::vector<double>& y, const std::vector<double>& d)
	: m_tag(tag)
	, m_type(type)
{
	// initial checks
	if (!(x.size() > 0 && x.size() == y.size() && x.size() == d.size())) {
		opserr << "ASDHardeningLaw Fatal Error: c-tor (material " << tag << ") - found incompatible sizes.\n";
		return;
	}
	// fill it
	m_points.resize(x.size());
	// make sure they are positive values
	double xmax = 0.0;
	double ymax = 0.0;
	for (std::size_t i = 0; i < x.size(); ++i) {
		auto& pi = m_points[i];
		pi.x = std::abs(x[i]);
		pi.y = std::abs(y[i]);
		pi.d = std::min(1.0, std::max(0.0, d[i]));
		xmax = std::max(xmax, pi.x);
		ymax = std::max(ymax, pi.y);
	}
	if (xmax == 0.0) {
		opserr << "ASDHardeningLaw Fatal Error: c-tor (material " << tag << ") - max(X) == 0 " << xmax << "\n";
		return;
	}
	if (ymax == 0.0) {
		opserr << "ASDHardeningLaw Fatal Error: c-tor (material " << tag << ") - max(Y) == 0 " << ymax << "\n";
		return;
	}
	// check first values and, if needed, append (0,0,0)
	if (m_points[0].x > 0.0) {
		if (m_points[0].y == 0.0) {
			// first x is > 0, y cannot be 0, so we set it to E*x
			m_points[0].y = E * m_points[0].x;
		}
		else {
			// ok, the user did not define the first point at (0,0), we make it
			m_points.insert(m_points.begin(), ASDHardeningLawPoint());
		}
	}
	else {
		// first x is 0, set y to zero as well
		m_points[0].y = 0.0;
	}
	// make sure the first damage is 0!
	m_points[0].d = 0.0;
	// make sure the slope of the first 2 points is E (note p1 is now 0,0)
	m_points[1].y = E * m_points[1].x;
	// define tolerances: find the min(dx or dy) > 0
	double dxmin = xmax;
	double dymin = ymax;
	for (std::size_t i = 1; i < m_points.size(); ++i) {
		const auto& pi = m_points[i];
		const auto& pold = m_points[i - 1];
		double dx = std::abs(pi.x - pold.x);
		double dy = std::abs(pi.y - pold.y);
		if (dx > 0.0)
			dxmin = std::min(dx, dxmin);
		if (dy > 0.0)
			dymin = std::min(dy, dymin);
	}
	m_xtolerance = 1.0e-6 * dxmin;
	m_ytolerance = 1.0e-6 * dymin;
	// make valid
	m_valid = true;
	// make a first adjustment
	adjust();
	// compute fracture energy. If the user selects the auto-regularization,
	// this will store the real (non-regularized) fracture energy
	computeFractureEnergy();
	// take an identity and keep an immutable copy of the law as it is now:
	// this is what deRegularize() goes back to, and it is shared by every
	// copy of the parent material
	m_uid = ASDHardeningLawStorage::instance().generateUID();
	m_pristine = std::make_shared<const ASDHardeningLaw>(*this);
}

void ASDHardeningLaw::regularize(double lch, double lch_ref)
{
	// quick return if not valid
	if (!m_valid)
		return;
	// quick return if un-bounded (inf fracture energy), or invalid lch
	double lch_scale = lch > 0.0 ? lch_ref / lch : 0.0;
	if (!m_fracture_energy_is_bounded || lch_scale <= 0.0 || lch_scale == 1.0)
		return;
	// back to original
	deRegularize();
	// the initial fracture energy has been computed in the full constructor, and we
	// are back to it after deRegularize...
	// compute the required specific fracture energy
	double gnew = m_fracture_energy * lch_scale;
	// compute the minimum fracture energy (in case lch is too large)
	const auto& peak = m_points[m_softening_begin];
	double gmin = (peak.y * peak.x / 2.0) * 1.01; // make it 1% larger to ensure a monotonically increasing abscissa
	gnew = std::max(gnew, gmin);
	// iteratively scale the curve until the gnew
	double tol = 1.0e-3 * gnew;
	double dscale = gnew / m_fracture_energy;
	double E = m_points[1].y / m_points[1].x;
	double x0 = peak.x;
	double scale = dscale;
	static constexpr int max_iter = 10;
	for (int iter = 0; iter < max_iter; ++iter) {
		// scale points after the peak
		for (std::size_t i = m_softening_begin + 1; i < m_points.size(); ++i) {
			auto& pi = m_points[i];
			// save the plastic-to-inelastic ratio before
			double xi_inel = std::max(pi.x - pi.y / E, 0.0);
			double xi_pl = pi.x - pi.q / E;
			double xi_ratio = xi_inel > 0.0 ? xi_pl / xi_inel : 0.0;
			// update the total strain
			pi.x = x0 + (pi.x - x0) * dscale;
			// update damage to keep the same ratio
			xi_inel = std::max(pi.x - pi.y / E, 0.0);
			xi_pl = xi_inel * xi_ratio;
			pi.q = E * (pi.x - xi_pl);
			if (pi.q > 0)
				pi.d = 1.0 - pi.y / pi.q;
		}
		// update g and check
		adjust();
		computeFractureEnergy();
		if (std::abs(m_fracture_energy - gnew) < tol)
			break;
		// update scale
		dscale = gnew / m_fracture_energy;
		scale *= dscale;
	}
	// re-adjust.
	adjust();
}

void ASDHardeningLaw::deRegularize()
{
	// the pristine law is the one this law was born as, before any
	// regularization. Note that the identity and the pointer to the pristine
	// copy are properties of THIS law and must survive the restore.
	if (m_pristine) {
		auto keep = m_pristine;
		std::size_t uid = m_uid;
		*this = *m_pristine;
		m_pristine = keep;
		m_uid = uid;
	}
}

ASDHardeningLawPoint ASDHardeningLaw::evaluateAt(double x) const
{
	// quick return
	if (!m_valid)
		return ASDHardeningLawPoint();
	// search for x
	double x1, x2, y1, y2, q1, q2;
	bool found = false;
	for (std::size_t i = 1; i < m_points.size(); ++i) {
		const auto& p1 = m_points[i - 1];
		const auto& p2 = m_points[i];
		if (x <= p2.x + m_xtolerance) {
			x1 = p1.x;
			x2 = p2.x;
			y1 = p1.y;
			y2 = p2.y;
			q1 = p1.q;
			q2 = p2.q;
			found = true;
			break;
		}
	}
	if (!found) {
		// we're beyond last point
		x1 = m_points.back().x;
		x2 = x;
		double span = x - x1;
		// interpolate last tangent if positive, otherwise keep constant
		y1 = m_points.back().y;
		double tangent = (y1 - m_points[m_points.size() - 2].y) / (x1 - m_points[m_points.size() - 2].x);
		y2 = tangent > 0.0 ? y1 + span * tangent : y1;
		// interpolate last tangent if positive, otherwise keep constant
		q1 = m_points.back().q;
		tangent = (q1 - m_points[m_points.size() - 2].q) / (x1 - m_points[m_points.size() - 2].x);
		q2 = tangent > 0.0 ? q1 + span * tangent : q1;
	}
	// interpolate
	double xspan = x2 - x1;
	double xratio = xspan > 0.0 ? (x - x1) / xspan : 0.0;
	double y = std::max(m_ytolerance, y1 + (y2 - y1) * xratio);
	double q = std::max(m_ytolerance, q1 + (q2 - q1) * xratio);
	double d = 1.0 - y / q;
	// done
	return ASDHardeningLawPoint(x, y, d, q);
}

double ASDHardeningLaw::slopeAt(double x) const
{
	// the segment search of evaluateAt(), with the slope instead of the
	// interpolated value. Kept here rather than differenced by the caller so
	// that the two can never disagree about which segment x falls on
	if (!m_valid)
		return 0.0;
	for (std::size_t i = 1; i < m_points.size(); ++i) {
		const auto& p1 = m_points[i - 1];
		const auto& p2 = m_points[i];
		if (x <= p2.x + m_xtolerance) {
			double xspan = p2.x - p1.x;
			return xspan > 0.0 ? (p2.y - p1.y) / xspan : 0.0;
		}
	}
	// beyond the last point: same rule as evaluateAt, the last tangent if it
	// is positive, otherwise perfectly plastic
	if (m_points.size() < 2)
		return 0.0;
	const auto& pl = m_points.back();
	const auto& pp = m_points[m_points.size() - 2];
	double xspan = pl.x - pp.x;
	if (xspan <= 0.0)
		return 0.0;
	double tangent = (pl.y - pp.y) / xspan;
	return tangent > 0.0 ? tangent : 0.0;
}

double ASDHardeningLaw::computeMaxStress() const
{
	double smax = 0.0;
	for (const auto& p : m_points) {
		smax = std::max(smax, p.y);
	}
	return smax;
}

int ASDHardeningLaw::serializationDataSize() const
{
	// number of points (variable, 4 components each)
	int np = static_cast<int>(m_points.size());
	// number of fixed data
	int nn = 10;
	// we need to save 2 copies (the current and the pristine one), plus the
	// identity, which is sent once and is what lets the receiver share a
	// single pristine law among all the copies that carry it
	return (nn + np * 4) * 2 + 1;
}

void ASDHardeningLaw::serialize(Vector& data, int& pos)
{
	// internal
	auto lam = [&data, &pos](const ASDHardeningLaw& x) {
		data(pos++) = static_cast<double>(x.m_tag);
		data(pos++) = static_cast<double>(static_cast<int>(x.m_type));
		data(pos++) = static_cast<double>(x.m_points.size());
		data(pos++) = x.m_fracture_energy;
		data(pos++) = static_cast<double>(x.m_fracture_energy_is_bounded);
		data(pos++) = static_cast<double>(x.m_softening_begin);
		data(pos++) = static_cast<double>(x.m_softening_end);
		data(pos++) = static_cast<double>(x.m_valid);
		data(pos++) = x.m_xtolerance;
		data(pos++) = x.m_ytolerance;
		for (const auto& p : x.m_points) {
			data(pos++) = p.x;
			data(pos++) = p.y;
			data(pos++) = p.d;
			data(pos++) = p.q;
		}
	};
	// save the identity
	data(pos++) = static_cast<double>(m_uid);
	// save the current
	lam(*this);
	// save the original
	lam(m_pristine ? *m_pristine : *this);
}

void ASDHardeningLaw::deserialize(Vector& data, int& pos)
{
	// internal
	auto lam = [&data, &pos](ASDHardeningLaw& x) {
		x.m_tag = static_cast<int>(data(pos++));
		x.m_type = static_cast<ASDHardeningLawType>(static_cast<int>(data(pos++)));
		x.m_points.resize(static_cast<std::size_t>(data(pos++)));
		x.m_fracture_energy = data(pos++);
		x.m_fracture_energy_is_bounded = static_cast<bool>(data(pos++));
		x.m_softening_begin = static_cast<std::size_t>(data(pos++));
		x.m_softening_end = static_cast<std::size_t>(data(pos++));
		x.m_valid = static_cast<bool>(data(pos++));
		x.m_xtolerance = data(pos++);
		x.m_ytolerance = data(pos++);
		for (auto& p : x.m_points) {
			p.x = data(pos++);
			p.y = data(pos++);
			p.d = data(pos++);
			p.q = data(pos++);
		}
	};
	// recover the identity
	std::size_t uid = static_cast<std::size_t>(data(pos++));
	// recover the current
	lam(*this);
	// recover the original
	ASDHardeningLaw original;
	lam(original);
	// re-attach the identity and share the pristine law with every other copy
	// that came from the same law in the sending process. An invalid law has
	// no identity (uid 0) and no pristine copy, exactly as on the build path
	m_uid = uid;
	if (uid != 0) {
		original.m_uid = uid;
		m_pristine = ASDHardeningLawStorage::instance().internReceived(uid, original);
	}
}

void ASDHardeningLaw::adjust()
{
	// quick return
	if (!m_valid)
		return;
	// get initial tangent
	double E = m_points[1].y / m_points[1].x;
	// check all points
	for (std::size_t i = 1; i < m_points.size(); ++i) {
		const auto& pold = m_points[i - 1];
		auto& pi = m_points[i];
		double xi = pi.x;
		double xold = pold.x;
		double yi = pi.y;
		double yold = pold.y;
		double di = pi.d;
		double dold = pold.d;
		// check: strictly monotonic x
		if (xi <= xold)
			xi += m_xtolerance;
		// check yi is not exactly 0.0
		if (yi < m_ytolerance)
			yi = m_ytolerance;
		// check tangent stifness is not higher than the initial
		double Ei = (yi - yold) / (xi - xold); // denom > 0 < --strictly monotonic x
		if (Ei > E)
			yi = yold + (xi - xold) * E;
		// check damage
		// current plastic strain cannot be lower then the old one
		double eepd_old = dold < 1.0 ? xold - yold / ((1.0 - dold) * E) : xold;
		double eepd = di < 1.0 ? xi - yi / ((1.0 - di) * E) : -1.0;
		if (eepd < eepd_old)
			di = std::min(1.0, std::max(0.0, 1.0 - yi / E / (xi - eepd_old)));
		// damage cannot decrease
		if (di < dold)
			di = dold;
		// last check: make sure (1-d)*E is not < tangent (Ei)
		Ei = (yi - yold) / (xi - xold);
		double Ed = (1.0 - di) * E;
		if (Ei > Ed)
			yi = Ed * (xi - eepd_old);
		// update
		pi.x = xi;
		pi.y = yi;
		pi.d = di;
		pi.q = yi / (1.0 - di);
	}
}

void ASDHardeningLaw::computeFractureEnergy()
{
	// initialize as un-bounded
	m_fracture_energy = 0.0;
	m_fracture_energy_is_bounded = false;
	m_softening_begin = 0;
	m_softening_end = 0;
	// quick return
	if (!m_valid)
		return;
	// find the first point where the slope is negative
	std::size_t pos1 = 0;
	bool pos1_found = false;
	for (std::size_t i = 1; i < m_points.size(); ++i) {
		const auto& p1 = m_points[i - 1];
		const auto& p2 = m_points[i];
		double k = (p2.y - p1.y) / (p2.x - p1.x);
		if (k < 0.0) {
			pos1 = i - 1;
			pos1_found = true;
			break;
		}
	}
	// exit if not found -> infinite energy
	if (!pos1_found)
		return;
	// find the last point where the slope is negative
	std::size_t pos2 = 0;
	bool pos2_found = false;
	for (std::size_t i = pos1 + 1; i < m_points.size(); ++i) {
		const auto& p1 = m_points[i - 1];
		const auto& p2 = m_points[i];
		double k = (p2.y - p1.y) / (p2.x - p1.x);
		if (k >= 0.0) {
			pos2 = i - 1;
			pos2_found = true;
			break;
		}
	}
	// if the last point is not found, it means
	// that we need to extend the last portion
	double g_add = 0.0;
	if (!pos2_found) {
		const auto& p2 = m_points.back();
		if (p2.y > 0.0) {
			const auto& p1 = m_points[m_points.size() - 2];
			double k = (p2.y - p1.y) / (p2.x - p1.x);
			double x3 = p2.x - p2.y / k;
			g_add = p2.y * (x3 - p2.x) / 2.0;
		}
		pos2 = m_points.size() - 1;
	}
	// now compute g
	double g = 0.0;
	// first add the initial triangle based on the unloading
	// stiffness at the pos1 point
	double E = m_points[1].y / m_points[1].x;
	double Ed = (1.0 - m_points[pos1].d) * E;
	g += (std::pow(m_points[pos1].y, 2) / Ed / 2.0);
	// then add other components
	for (std::size_t i = pos1 + 1; i <= pos2; ++i) {
		const auto& p1 = m_points[i - 1];
		const auto& p2 = m_points[i];
		g += (p2.x - p1.x) * (p1.y + p2.y) / 2.0;
	}
	// finally add the final one
	g += g_add;
	// done. store values now
	m_fracture_energy = g;
	m_fracture_energy_is_bounded = true;
	m_softening_begin = pos1;
	m_softening_end = pos2;
}

ASDHardeningLawStorage& ASDHardeningLawStorage::instance()
{
	static ASDHardeningLawStorage _instance;
	return _instance;
}

std::size_t ASDHardeningLawStorage::generateUID()
{
	// never reset, not even by a 'wipe': that is the whole point
	return ++m_counter;
}

ASDHardeningLawStorage::PointerType ASDHardeningLawStorage::internReceived(std::size_t uid, const ASDHardeningLaw& law)
{
	// share an already received law with the same identity, if still alive
	auto it = m_received.find(uid);
	if (it != m_received.end()) {
		PointerType alive = it->second.lock();
		if (alive)
			return alive;
		m_received.erase(it);
	}
	// drop the entries whose law is gone. The map is as long as the number of
	// distinct laws ever received, so this stays cheap
	for (auto i = m_received.begin(); i != m_received.end(); ) {
		if (i->second.expired())
			i = m_received.erase(i);
		else
			++i;
	}
	// take a copy and share it
	PointerType item = std::make_shared<const ASDHardeningLaw>(law);
	m_received[uid] = item;
	return item;
}
