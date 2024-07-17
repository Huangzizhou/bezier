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
class ValidityChecker {
	private:
	static constexpr uint subdomains = 2 << n;
	Matrix<Interval> matL2B;
	std::pair<Matrix<Interval>, Matrix<Interval>> matT;
	std::array<Matrix<Interval>, subdomains> matQ;
	std::vector<uint> interpIndices;

	const uint nThreads;

	fp_t precision = .1;
	fp_t epsilon = 0.;
	uint maxSubdiv = 0;

	public:
	ValidityChecker(uint nThreads = 1) : nThreads(nThreads) {
		initMatricesT<n, s, p>(matL2B, matT, matQ);
		cornerIndicesT<n, s, p>(interpIndices);
	}
	fp_t maxTimeStep(
		std::span<const fp_t> cp,
		std::vector<uint> *adaptiveHierarchy = nullptr,
		fp_t earlyStop = 1,
		fp_t *timeOfInversion = nullptr,
		CheckerInfo *info = nullptr
	) const;

	fp_t maxTimeStepVec(
		std::span<const fp_t> cp,
		uint *invalidElemID = nullptr,
		std::vector<uint> *adaptiveHierarchy = nullptr
	) const;

	void setPrecisionTarget(fp_t t) { precision = t; }
	void setEpsilon(fp_t t) { epsilon = t; }
	void setMaxSubdiv(uint v) { maxSubdiv = v; }

	private:
	Interval minclusion(const std::vector<Interval> &B) const;
	Subdomain split(const Subdomain &src, uint q) const;
	Subdomain splitTime(const Subdomain &src, bool t) const;
};

}


namespace element_validity {
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
Subdomain ValidityChecker<n, s, p>::splitTime(
	const Subdomain &src, bool t
) const {
	Subdomain res;
	(t ? matT.second : matT.first).mult(src.B, res.B);
	res.time = t ?
		Interval(src.time.middle(), src.time.upper()) :
		Interval(src.time.lower(), src.time.middle());
	res.incl = minclusion(res.B);
	src.copySequence(res.qSequence);
	return res;
}


//------------------------------------------------------------------------------

template<uint n, uint s, uint p>
fp_t ValidityChecker<n, s, p>::maxTimeStep(
	std::span<const fp_t> cp,
	std::vector<uint> *hierarchy,
	fp_t earlyStop,
	fp_t *timeOfInversion,
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
fp_t ValidityChecker<n, s, p>::maxTimeStepVec(
	std::span<const fp_t> cp,
	uint *invalidElemID,
	std::vector<uint> *adaptiveHierarchy
) const {
	const uint numCoordsPerElem = nControlGeoMap(n,s,p) * 2 * n;
	const uint numEl = cp.size() / (numCoordsPerElem);
	std::cerr << "number of polys: " << numEl << std::endl;
	fp_t minT = 1;
	std::vector<fp_t> timeOfInversion(numEl);
	std::vector<fp_t> depthOfSequence(numEl);
	std::vector<fp_t> timings(numEl);
	std::vector<std::vector<fp_t>> hierarchies(numEl);
	bool foundInvalid = false;

	#pragma omp parallel for \
		reduction(min : minT) num_threads(nThreads)
	for (uint e=0; e<numEl; ++e) {
		std::span<const fp_t> element(
			cp.data() + numCoordsPerElem * e, numCoordsPerElem);
		Timer timer;
		fp_t invT;
		timer.start();
		fp_t t = maxTimeStep(element, hierarchies.at(e), minT, &invT);
		timer.stop();
		timeOfInversion.at(e) = invT;
		if (invT >= 0) foundInvalid = true;
		timings.at(e) = timer.read<std::chrono::microseconds>();
		minT = std::min(minT, t);
	}
	if (foundInvalid) {
		const uint i = std::distance(timeOfInversion.cbegin(),
			std::min_element(
				timeOfInversion.cbegin(),
				timeOfInversion.cend()
			)
		);
		if (invalidElemID) *invalidElemID = i;
		*adaptiveHierarchy = std::move(hierarchies.at(i));
	}
	else {
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

}