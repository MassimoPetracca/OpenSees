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

#include <IMPLEXManager.h>
#include <algorithm>
#include <cmath>

IMPLEXObject::IMPLEXObject()
{
	m_implex_uid = IMPLEXManager::instance().add(this);
}

IMPLEXObject::IMPLEXObject(const IMPLEXObject&)
{
	// a copy is a new object: it gets its own slot, and it has not taken part
	// in anything yet
	m_implex_uid = IMPLEXManager::instance().add(this);
}

IMPLEXObject::~IMPLEXObject()
{
	IMPLEXManager::instance().remove(this);
}

IMPLEXManager& IMPLEXManager::instance()
{
	// deliberately never destroyed. Every IMPLEXObject reaches back here from
	// its destructor, and objects held by static containers elsewhere are
	// destroyed in an order nobody controls: a registry that dies first would
	// be touched afterwards
	static IMPLEXManager* _instance = new IMPLEXManager();
	return *_instance;
}

std::size_t IMPLEXManager::add(IMPLEXObject* obj)
{
	m_objects.push_back(obj);
	return m_objects.size() - 1;
}

void IMPLEXManager::remove(IMPLEXObject* obj)
{
	std::size_t i = obj->m_implex_uid;
	// defensive: an object that is not where it says it is would corrupt the
	// registry, and the cost of the check is one comparison
	if (i >= m_objects.size() || m_objects[i] != obj) {
		auto it = std::find(m_objects.begin(), m_objects.end(), obj);
		if (it == m_objects.end())
			return;
		i = static_cast<std::size_t>(it - m_objects.begin());
	}
	// swap with the last one and tell it where it went
	IMPLEXObject* last = m_objects.back();
	m_objects[i] = last;
	last->m_implex_uid = i;
	m_objects.pop_back();
}

void IMPLEXManager::clearTouched()
{
	for (IMPLEXObject* obj : m_objects)
		obj->m_implex_touched = false;
	m_aggregate = Aggregate();
}

const IMPLEXManager::Aggregate& IMPLEXManager::aggregate(double overThreshold)
{
	for (IMPLEXObject* obj : m_objects) {
		if (!obj->m_implex_touched)
			continue;
		obj->m_implex_touched = false;
		double e = obj->computeImplexErrorMetric();
		if (std::isnan(e)) {
			m_aggregate.any_nan = true;
			continue;
		}
		m_aggregate.max = std::max(m_aggregate.max, e);
		m_aggregate.sum += e;
		++m_aggregate.count;
		// strictly over, so that a threshold of zero counts the active ones and
		// a threshold equal to the tolerance counts the violations - a point
		// exactly AT tolerance is within it, which is what the maximum-only
		// criterion has always said
		if (e > overThreshold)
			++m_aggregate.count_over;
		m_aggregate.min_time_ratio = std::min(m_aggregate.min_time_ratio, obj->implexTimeRatio());
	}
	return m_aggregate;
}
