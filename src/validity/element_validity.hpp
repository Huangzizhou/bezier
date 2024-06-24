#pragma once
#include "utils/Timer.hpp"
#include "transMatrices.hpp"
#include "lagrangeVector.hpp"
#include "cornerIndices.hpp"

#include <queue>

namespace element_validity {
struct Subdomain {
	public:
	std::vector<Interval> B;
	std::vector<uint> qSequence;
	Interval time;
	Interval incl;

	// Priority function
	bool operator<(const Subdomain &o) const {
		if (time.lower() != o.time.lower())
			return time.lower() > o.time.lower();
		else return incl.lower() > o.incl.lower();
	}


	// Debug print
	friend std::ostream& operator<<(std::ostream& ost, const Subdomain &s) {
		ost << "t: " << s.time << std::endl;
		ost << "I: " << s.incl << std::endl;
		ost << "Q: ";
		for (const uint x : s.qSequence) ost << x << ", ";
		ost << std::endl;
		ost << "B: ";
		for (const Interval &b : s.B) ost << '\t' << b << std::endl;

		return ost;
	}

	// Copy the path I took to get here
	void copySequence(std::vector<uint> &dst) const {
		dst.clear();
		dst.reserve(qSequence.size());
		std::copy(
			qSequence.cbegin(), qSequence.cend(), std::back_inserter(dst)
		);
	}
};


template<uint n, uint s, uint p>
class ValidityChecker {
	private:
	static constexpr uint subdomains = 2 << n;
	Matrix<Interval> matL2B;
	std::pair<Matrix<Interval>, Matrix<Interval>> matT;
	std::array<Matrix<Interval>, subdomains> matQ;
	std::vector<uint> interpIndices;

	fp_t precision = .1;
	fp_t threshold = 0;
	uint maxSubdiv = 0;

	public:
	ValidityChecker() {
		initMatricesT<n, s, p>(matL2B, matT, matQ);
		cornerIndicesT<n, s, p>(interpIndices);
	}
	fp_t maxTimeStep(
		const std::vector<fp_t> &cp,
		std::array<uint, 3> *info = nullptr
	);

	void setPrecisionTarget(fp_t t) { precision = t; }
	void setThreshold(fp_t t) { threshold = t; }
	void setMaxSubdiv(uint v) { maxSubdiv = v; }

	private:
	Interval minclusion(const std::vector<Interval> &B);
	Subdomain split(const Subdomain &src, uint q);
	Subdomain splitTime(const Subdomain &src, bool t);
};

//------------------------------------------------------------------------------

template<uint n, uint s, uint p>
Interval ValidityChecker<n, s, p>::minclusion(const std::vector<Interval>& B) {
	Interval lo(std::numeric_limits<fp_t>::max());
	for (const Interval &b : B) lo = min(lo, b);
	Interval hi(std::numeric_limits<fp_t>::max());
	for (const uint c : interpIndices) hi = min(hi, B.at(c));
	return {lo.lower(), hi.upper()};
}

template<uint n, uint s, uint p>
Subdomain ValidityChecker<n, s, p>::split(
	const Subdomain &src, uint q
) {
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
) {
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
	const std::vector<fp_t> &cp,
	std::array<uint, 3> *info
) {
	assert(precision <= 1);
	assert(precision > 0);

	// Compute Lagrange coefficients
	std::vector<Interval> vL;
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
	uint reachedDepthT = 0;
	uint maxQueueSize = 0;
	const bool maxIterCheck = (maxSubdiv > 0);
	// bool foundInvalid = false;
	std::vector<uint> invalidSubdivSequence;
	std::vector<uint> deepestSubdivSequence;
	fp_t tmin = 0;
	fp_t tmax = 1;

	// Loop
	while(true) {
		if (queue.empty()) {
			tmin = 1;
			break;
		}
		// Check whether we reached precision
		if (tmax - tmin < precision && tmin > 0) break;

		// Check whether we satisfy early termination
		// if (tmin >= earlyStop) break;

		// Get box from top of the queue
		maxQueueSize = std::max(maxQueueSize, static_cast<uint>(queue.size()));
		const Subdomain dom = queue.top();
		queue.pop();

		assert(dom.time.lower() < tmax);
		assert(dom.time.lower() >= tmin);

		if (dom.qSequence.size() > reachedDepthS) {
			dom.copySequence(deepestSubdivSequence);
			reachedDepthS = deepestSubdivSequence.size();
		}

		uint td = 0;
		for (double d=dom.time.width(); d<1.; d*=2) { ++td; }
		reachedDepthT = std::max(reachedDepthT, td);

		// Check whether we need to give up
		if (maxIterCheck && (reachedDepthS >= maxSubdiv)) {
			std::cerr << "Reached max subdivisions: giving up" << std::endl;
			break;
		}

		// Update tmin
		tmin = dom.time.lower();
		
		// Subdomain contains invalidity
		if (dom.incl <= 0) {
			// foundInvalid = true;
			if (tmax > dom.time.upper()) {
				tmax = dom.time.upper();				// update tmax
				dom.copySequence(invalidSubdivSequence);	// save trace
			}
			// Split on time only and push to queue
			const auto mid = dom.time.middle();
			if (mid == dom.time.upper() || mid == dom.time.lower()) {
				// Cannot split
				std::cerr << "Cannot split in time: giving up" << std::endl;
				break;
			}
			queue.push(splitTime(dom, false));
			queue.push(splitTime(dom, true));
		}
		// Subdomain is undetermined and needs refinement
		else if (!(dom.incl > 0)) {
			// Split on all axes and push to queue
			for (uint q=0; q<subdomains; ++q) queue.push(split(dom, q));
		}
	}

	if (info) {
		info->at(0) = reachedDepthS;
		info->at(1) = reachedDepthT;
		info->at(2) = maxQueueSize;
	}

	return tmin;
}

//------------------------------------------------------------------------------

// template<uint n, uint s, uint p>
// fp_t ValidityChecker<n, s, p>::isValid(
// 	const std::vector<fp_t> &cp
// ) {

// }

}
