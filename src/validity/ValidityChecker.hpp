#pragma once
#include "utils/combinatorics.hpp"
#include "utils/Timer.hpp"
#include "transMatrices.hpp"
#include "lagrangeVector.hpp"
#include "cornerIndices.hpp"

#include <queue>
#include <span>
#include <memory>

namespace element_validity {
struct Subdomain {
	public:
	std::vector<Interval> B;
	std::vector<uint> qSequence;
	Interval time;
	Interval incl;

	// Priority function
	bool operator<(const Subdomain &o) const;
	// Debug print
	friend std::ostream& operator<<(std::ostream& ost, const Subdomain &s);
	inline uint depth() const { return qSequence.size(); }
	// Copy the path I took to get here
	void copySequence(std::vector<uint> &dst) const;
};


struct CheckerInfo {
	uint spaceDepth = 0;
	enum class Status {
		unknown, completed, reachedTarget, globalCondition, maxDepth, noSplit
	};
	Status status;

	bool success() const;
	std::string description() const;
};


template<uint n, uint s, uint p>
class ContinuousValidityChecker {
	private:
	static constexpr uint subdomains = 2 << n;
	Matrix<Interval> matL2B;
	std::array<Matrix<Interval>, 2> matT;
	std::array<Matrix<Interval>, subdomains> matQ;
	std::vector<uint> interpIndices;

	const uint nThreads;

	fp_t precision = .1;
	fp_t epsilon = 0.;
	uint maxSubdiv = 0;

	fp_t maxTimeStepElement(
		std::span<const fp_t> cp,
		std::vector<uint> *adaptiveHierarchy = nullptr,
		fp_t *timeOfInversion = nullptr,
		fp_t earlyStop = 1,
		CheckerInfo *info = nullptr
	) const;

	fp_t maxTimeStepMesh(
		std::span<const fp_t> cp,
		std::vector<uint> *adaptiveHierarchy = nullptr,
		uint *invalidElemID = nullptr,
		fp_t *timeOfInversion = nullptr
	) const;

	public:
	ContinuousValidityChecker(uint nThreads = 1) : nThreads(nThreads) {
		initMatricesT<n, s, p>(matL2B, matT, matQ);
		cornerIndicesT<n, s, p>(interpIndices);
	}
	fp_t maxTimeStep(
		std::span<const fp_t> cp,
		std::vector<uint> *adaptiveHierarchy = nullptr,
		uint *invalidElemID = nullptr,
		fp_t *timeOfInversion = nullptr,
		CheckerInfo *info = nullptr
	) const;

	void setPrecisionTarget(fp_t t) { precision = t; }
	void setEpsilon(fp_t t) { epsilon = t; }
	void setMaxSubdiv(uint v) { maxSubdiv = v; }

	private:
	Interval minclusion(const std::vector<Interval> &B) const;
	Subdomain split(const Subdomain &src, uint q) const;
	Subdomain splitTime(const Subdomain &src, bool t) const;
};

template<uint n, uint s, uint p>
class ValidityChecker {
	private:
	static constexpr uint subdomains = 1 << n;
	Matrix<Interval> matL2B;
	std::array<Matrix<Interval>, subdomains> matQ;
	std::vector<uint> interpIndices;

	const uint nThreads;

	fp_t epsilon = 0.;
	uint maxSubdiv = 0;

	Validity isValidElement(
		std::span<const fp_t> cp,
		std::vector<uint> *adaptiveHierarchy = nullptr,
		CheckerInfo *info = nullptr
	) const;
	Validity isValidMesh(
		std::span<const fp_t> cp,
		std::vector<uint> *adaptiveHierarchy = nullptr,
		uint *invalidElemID = nullptr
	) const;

	public:
	ValidityChecker(uint nThreads = 1) : nThreads(nThreads) {
		initMatrices<n, s, p>(matL2B, matQ);
		cornerIndices<n, s, p>(interpIndices);
	}
	
	Validity isValid(
		std::span<const fp_t> cp,
		std::vector<uint> *adaptiveHierarchy = nullptr,
		uint *invalidElemID = nullptr,
		CheckerInfo *info = nullptr
	) const;

	void setEpsilon(fp_t t) { epsilon = t; }
	void setMaxSubdiv(uint v) { maxSubdiv = v; }

	private:
	Interval minclusion(const std::vector<Interval> &B) const;
	Subdomain split(const Subdomain &src, uint q) const;
};

}


namespace element_validity {
template<uint n, uint s, uint p>
Interval ContinuousValidityChecker<n, s, p>::minclusion(
	const std::vector<Interval>& B
) const {
	Interval lo(std::numeric_limits<fp_t>::max());
	for (const Interval &b : B) lo = min(lo, b);
	Interval hi(std::numeric_limits<fp_t>::max());
	for (const uint c : interpIndices) hi = min(hi, B.at(c));
	return {lo.lower(), hi.upper()};
}

template<uint n, uint s, uint p>
Subdomain ContinuousValidityChecker<n, s, p>::split(
	const Subdomain &src, uint q
) const {
	Subdomain res;
	matQ[q].mult(src.B, res.B);
	res.time = (q >> n) ?
		Interval(src.time.middle(), src.time.upper()) :
		Interval(src.time.lower(), src.time.middle());
	res.incl = minclusion(res.B);
	src.copySequence(res.qSequence);
	res.qSequence.push_back(q % (1U << (n-1)));
	return res;
}

template<uint n, uint s, uint p>
Subdomain ContinuousValidityChecker<n, s, p>::splitTime(
	const Subdomain &src, bool t
) const {
	Subdomain res;
	(t ? matT[1] : matT[0]).mult(src.B, res.B);
	res.time = t ?
		Interval(src.time.middle(), src.time.upper()) :
		Interval(src.time.lower(), src.time.middle());
	res.incl = minclusion(res.B);
	src.copySequence(res.qSequence);
	return res;
}


template<uint n, uint s, uint p>
Subdomain ValidityChecker<n, s, p>::split(
	const Subdomain &src, uint q
) const {
	Subdomain res;
	matQ[q].mult(src.B, res.B);
	res.time = (q >> n) ?
		Interval(src.time.middle(), src.time.upper()) :
		Interval(src.time.lower(), src.time.middle());
	res.incl = minclusion(res.B);
	src.copySequence(res.qSequence);
	res.qSequence.push_back(q % (1U << (n-1)));
	return res;
}


template<uint n, uint s, uint p>
Interval ValidityChecker<n, s, p>::minclusion(
	const std::vector<Interval>& B
) const {
	Interval lo(std::numeric_limits<fp_t>::max());
	for (const Interval &b : B) lo = min(lo, b);
	Interval hi(std::numeric_limits<fp_t>::max());
	for (const uint c : interpIndices) hi = min(hi, B.at(c));
	return {lo.lower(), hi.upper()};
}

//------------------------------------------------------------------------------

template<uint n, uint s, uint p>
fp_t ContinuousValidityChecker<n,s,p>::maxTimeStep(
	std::span<const fp_t> cp,
	std::vector<uint> *adaptiveHierarchy,
	uint *invalidElemID,
	fp_t *timeOfInversion,
	CheckerInfo *info
) const {
	const uint numCoordsPerElem = nControlGeoMap(n,s,p) * 2 * n;
	const uint numEl = cp.size() / (numCoordsPerElem);
	if (numEl == 1) return maxTimeStepElement(
		cp, adaptiveHierarchy, timeOfInversion, 2, info);
	else return maxTimeStepMesh(
		cp, adaptiveHierarchy, invalidElemID, timeOfInversion);
}

//------------------------------------------------------------------------------

template<uint n, uint s, uint p>
Validity ValidityChecker<n,s,p>::isValid(
	std::span<const fp_t> cp,
	std::vector<uint> *adaptiveHierarchy,
	uint *invalidElemID,
	CheckerInfo *info
) const {
	const uint numCoordsPerElem = nControlGeoMap(n,s,p) * n;
	const uint numEl = cp.size() / (numCoordsPerElem);
	if (numEl == 1) return isValidElement(cp, adaptiveHierarchy, info);
	else return isValidMesh(cp, adaptiveHierarchy, invalidElemID);
}

//------------------------------------------------------------------------------

template<uint n, uint s, uint p>
fp_t ContinuousValidityChecker<n, s, p>::maxTimeStepElement(
	std::span<const fp_t> cp,
	std::vector<uint> *hierarchy,
	fp_t *timeOfInversion,
	fp_t earlyStop,
	CheckerInfo *info
) const {
	assert(precision > 0 && precision <= 1);

	// Compute Lagrange coefficients
	std::vector<Interval> vL(nControlJacobian(n,s,p,true));
	lagrangeVectorT<n, s, p>(cp, vL);

	// Create initial subdomain
	Subdomain sd0;
	matL2B.mult(vL, sd0.B);
	sd0.time = Interval(0,1);
	sd0.incl = minclusion(sd0.B);
	
	// Initialize queue
	std::priority_queue<Subdomain> queue;
	queue.push(sd0);

	// Initialize auxiliary variables
	uint reachedDepthS = 0;
	const bool maxIterCheck = (maxSubdiv > 0);
	bool foundInvalid = false;
	fp_t tmin = 0.;
	fp_t tmax = 1.;

	if (hierarchy) hierarchy->clear();

	// Loop
	while(true) {
		if (queue.empty()) {
			tmin = 1;
			if (info) info->status = CheckerInfo::Status::completed;
			break;
		}
		// Check whether we reached precision
		if (tmax - tmin < precision && tmin > 0.) {
			if (info) info->status = CheckerInfo::Status::reachedTarget;
			break;
		}

		// Check whether we satisfy early termination
		if (tmin >= earlyStop) {
			if (info) info->status = CheckerInfo::Status::globalCondition;
			break;
		}

		// Get box from top of the queue
		const Subdomain dom = queue.top();
		queue.pop();

		assert(dom.time.lower() < tmax);
		assert(dom.time.lower() >= tmin);

		reachedDepthS = std::max(reachedDepthS, dom.depth());

		// Check whether we need to give up
		if (maxIterCheck && (reachedDepthS >= maxSubdiv)) {
			if (hierarchy && !foundInvalid) dom.copySequence(*hierarchy);
			if (info) info->status = CheckerInfo::Status::maxDepth;
			break;
		}

		// Update tmin
		tmin = dom.time.lower();
		
		// Subdomain contains invalidity
		if (dom.incl <= epsilon) {
			foundInvalid = true;
			if (tmax > dom.time.upper()) {
				tmax = dom.time.upper();				// update tmax
				if (hierarchy) dom.copySequence(*hierarchy);
			}
			// Split on time only and push to queue
			const auto mid = dom.time.middle();
			if (mid == dom.time.upper() || mid == dom.time.lower()) {
				if (info) info->status = CheckerInfo::Status::noSplit;
				break;
			}
			queue.push(splitTime(dom, false));
			queue.push(splitTime(dom, true));
		}
		// Subdomain is undetermined and needs refinement
		else if (!(dom.incl > epsilon)) {
			// Split on all axes and push to queue
			for (uint q=0; q<subdomains; ++q) queue.push(split(dom, q));
		}
	}

	if (info) info->spaceDepth = reachedDepthS;
	if (timeOfInversion) *timeOfInversion = foundInvalid ? tmax : -1.;
	return tmin;
}

//------------------------------------------------------------------------------

template<uint n, uint s, uint p>
fp_t ContinuousValidityChecker<n, s, p>::maxTimeStepMesh(
	std::span<const fp_t> cp,
	std::vector<uint> *adaptiveHierarchy,
	uint *invalidElemID,
	fp_t *toi
) const {
	const uint numCoordsPerElem = nControlGeoMap(n,s,p) * 2 * n;
	const uint numEl = cp.size() / (numCoordsPerElem);
	fp_t minT = 1;
	std::vector<fp_t> timeOfInversion(numEl);
	std::vector<fp_t> timings(numEl);
	std::vector<std::vector<uint>> hierarchies(numEl);
	bool foundInvalid = false;

	#pragma omp parallel for \
		reduction(min : minT) num_threads(nThreads)
	for (uint e=0; e<numEl; ++e) {
		std::span<const fp_t> element(
			cp.data() + numCoordsPerElem * e, numCoordsPerElem);
		Timer timer;
		fp_t invT;
		timer.start();
		fp_t t = maxTimeStepElement(element, &hierarchies.at(e), &invT, minT);
		timer.stop();
		timeOfInversion.at(e) = invT;
		if (invT >= 0) foundInvalid = true;
		timings.at(e) = timer.read<std::chrono::microseconds>();
		minT = std::min(minT, t);
	}
	if (foundInvalid) {
		const auto m = std::min_element(timeOfInversion.cbegin(), timeOfInversion.cend()); 
		const uint i = std::distance(timeOfInversion.cbegin(), m);
		if (toi) *toi = *m;
		if (invalidElemID) *invalidElemID = i;
		*adaptiveHierarchy = std::move(hierarchies.at(i));
	}
	else {
		std::vector<fp_t> depthOfSequence(numEl);
		for (uint j = 0; j < numEl; ++j)
			depthOfSequence.at(j) = hierarchies.at(j).size();
		const uint i = std::distance(depthOfSequence.cbegin(),
			std::max_element(
				depthOfSequence.cbegin(),
				depthOfSequence.cend()
			)
		);
		if (invalidElemID) *invalidElemID = i;
		*adaptiveHierarchy = std::move(hierarchies.at(i));
	}
	return minT;
}


//------------------------------------------------------------------------------

template<uint n, uint s, uint p>
Validity ValidityChecker<n, s, p>::isValidElement(
	std::span<const fp_t> cp,
	std::vector<uint> *hierarchy,
	CheckerInfo *info
) const {
	// Compute Lagrange coefficients
	std::vector<Interval> vL(nControlJacobian(n,s,p,false));
	lagrangeVector<n, s, p>(cp, vL);

	// Check for early termination
	for (const Interval &l : vL) {
		if (l < 0) {
			if (info) info->status = CheckerInfo::Status::reachedTarget;
			return Validity::invalid;
		}
	}


	// Create initial subdomain
	Subdomain sd0;
	matL2B.mult(vL, sd0.B);
	sd0.incl = minclusion(sd0.B);
	
	// Initialize queue
	std::priority_queue<Subdomain> queue;
	queue.push(sd0);

	// Initialize auxiliary variables
	uint reachedDepth = 0;
	const bool maxIterCheck = (maxSubdiv > 0);
	bool foundInvalid = false;
	bool gaveUp = false;

	if (hierarchy) hierarchy->clear();

	// Loop
	while(true) {
		if (queue.empty()) {
			if (info) info->status = CheckerInfo::Status::completed;
			break;
		}
	
		// Get box from top of the queue
		const Subdomain dom = queue.top();
		queue.pop();

		reachedDepth = std::max(reachedDepth, dom.depth());

		// Check whether we need to give up
		if (maxIterCheck && (reachedDepth >= maxSubdiv)) {
			gaveUp = true;
			if (hierarchy && !foundInvalid) dom.copySequence(*hierarchy);
			if (info) info->status = CheckerInfo::Status::maxDepth;
			break;
		}

		// Subdomain contains invalidity
		if (dom.incl <= epsilon) {
			foundInvalid = true;
			if (info) info->status = CheckerInfo::Status::reachedTarget;
		}
		// Subdomain is undetermined and needs refinement
		else if (!(dom.incl > epsilon)) {
			// Split on all axes and push to queue
			for (uint q=0; q<subdomains; ++q) queue.push(split(dom, q));
		}
	}

	if (info) info->spaceDepth = reachedDepth;
	return gaveUp ? Validity::uncertain : (foundInvalid ? Validity::invalid : Validity::valid);
}

//------------------------------------------------------------------------------

template<uint n, uint s, uint p>
Validity ValidityChecker<n, s, p>::isValidMesh(
	std::span<const fp_t> cp,
	std::vector<uint> *adaptiveHierarchy,
	uint *invalidElemID
) const {
	const uint numCoordsPerElem = nControlGeoMap(n,s,p) * n;
	const uint numEl = cp.size() / (numCoordsPerElem);
	std::cerr << "number of polys: " << numEl << std::endl;
	std::vector<Validity> results(numEl);
	std::vector<fp_t> timings(numEl);
	std::vector<std::vector<uint>> hierarchies(numEl);
	bool gaveUp = false;
	bool foundInvalid = false;

	#pragma omp parallel for num_threads(nThreads)
	for (uint e=0; e<numEl; ++e) {
		if (foundInvalid) {
			timings.at(e) = 0;
			continue;
		}
		std::span<const fp_t> element(
			cp.data() + numCoordsPerElem * e, numCoordsPerElem);
		Timer timer;
		timer.start();
		Validity v = isValidElement(element, &hierarchies.at(e));
		timer.stop();
		results.at(e) = v;
		if (v == Validity::invalid) foundInvalid = true;
		if (v == Validity::uncertain) gaveUp = true;
		timings.at(e) = timer.read<std::chrono::microseconds>();
	}
	if (foundInvalid) {
		uint invalidIndex = 0;
		while (results.at(invalidIndex) != Validity::invalid) ++invalidIndex;
		if (invalidElemID) *invalidElemID = invalidIndex;
		*adaptiveHierarchy = std::move(hierarchies.at(invalidIndex));
		return Validity::invalid;
	}
	else if (gaveUp) {
		std::vector<int> depthOfSequence(numEl);
		for (uint j = 0; j < numEl; ++j)
			depthOfSequence.at(j) = (results.at(j) == Validity::valid) ?
				-1 : hierarchies.at(j).size();
		const uint i = std::distance(depthOfSequence.cbegin(),
			std::max_element(
				depthOfSequence.cbegin(),
				depthOfSequence.cend()
			)
		);
		if (invalidElemID) *invalidElemID = i;
		*adaptiveHierarchy = std::move(hierarchies.at(i));
		return Validity::uncertain;
	}
	return Validity::valid;
}

}
