#pragma once
#include "Validator.hpp"

namespace element_validity {
template<uint n, uint s, uint p>
class StaticValidator : public Validator {
	private:
	static constexpr uint subdomains = 1U << n;
	std::array<Matrix<Interval>, subdomains> matQ;

	Validity isValidElement(
		span<const fp_t> cp,
		std::vector<uint> *adaptiveHierarchy = nullptr,
		Info *info = nullptr
	) const;
	Validity isValidMesh(
		span<const fp_t> cp,
		std::vector<uint> *adaptiveHierarchy = nullptr,
		uint *invalidElemID = nullptr,
		std::vector<uint> *invalidList = nullptr
	) const;

	Subdomain split(const Subdomain &src, uint q) const;

	public:
	StaticValidator(uint nThreads = 1) : Validator(nThreads) {
		initMatrices<n, s, p>(matL2B, matQ);
		cornerIndices<n, s, p>(interpIndices);
	}
	
	Validity isValid(
		span<const fp_t> cp,
		std::vector<uint> *adaptiveHierarchy = nullptr,
		uint *invalidElemID = nullptr,
		std::vector<uint> *invalidList = nullptr,
		Info *info = nullptr
	) const;

	Validity isValidAtTime(
		span<const fp_t> cp,
		fp_t time,
		std::vector<uint> *adaptiveHierarchy = nullptr,
		uint *invalidElemID = nullptr,
		std::vector<uint> *invalidList = nullptr,
		Info *info = nullptr
	) const;

	#ifdef EIGEN_INTERFACE
	Validity isValid(
		const Eigen::MatrixXd& cp,
		std::vector<uint> *adaptiveHierarchy = nullptr,
		uint *invalidElemID = nullptr,
		std::vector<uint> *invalidList = nullptr,
		Info *info = nullptr
	) const {
		return isValid(
			convertEigenMatrix(cp),
			adaptiveHierarchy,
			invalidElemID,
			invalidList,
			info
		);
	}
	#endif
};

//------------------------------------------------------------------------------

template<uint n, uint s, uint p>
Validator::Subdomain StaticValidator<n,s,p>::split(
	const Validator::Subdomain &src, uint q
) const {
	Subdomain res;
	matQ[q].mult(src.B, res.B);
	res.incl = minclusion(res.B);
	src.copySequence(res.qSequence);
	res.qSequence.push_back(q);
	return res;
}

//------------------------------------------------------------------------------

template<uint n, uint s, uint p>
Validity StaticValidator<n,s,p>::isValid(
	span<const fp_t> cp,
	std::vector<uint> *adaptiveHierarchy,
	uint *invalidElemID,
	std::vector<uint> *invalidList,
	Info *info
) const {
	constexpr uint numCoordsPerElem = nControlGeoMap(n,s,p) * n;
	const uint numEl = cp.size() / numCoordsPerElem;
	if (numEl == 1) return isValidElement(cp, adaptiveHierarchy, info);
	else return isValidMesh(cp, adaptiveHierarchy, invalidElemID, invalidList);
}

template<uint n, uint s, uint p>
Validity StaticValidator<n,s,p>::isValidAtTime(
	span<const fp_t> cp,
	fp_t time,
	std::vector<uint> *adaptiveHierarchy,
	uint *invalidElemID,
	std::vector<uint> *invalidList,
	Info *info
) const {
	std::vector<fp_t> frame;
	assert(cp.size() % 2 == 0);
	frame.reserve(cp.size() / 2);
	const fp_t time1m = 1.-time;
	for (uint i = 0; i < cp.size(); i+=2)
		frame.push_back(time1m*cp[i] + time*cp[i+1]);
	return isValid(frame, adaptiveHierarchy, invalidElemID, invalidList, info);
}

//------------------------------------------------------------------------------

template<uint n, uint s, uint p>
Validity StaticValidator<n, s, p>::isValidElement(
	span<const fp_t> cp,
	std::vector<uint> *hierarchy,
	Info *info
) const {
	// Compute Lagrange coefficients
	std::vector<Interval> vL(nControlJacobian(n,s,p,false));
	lagrangeVector<n, s, p>(cp, vL);

	// Check for early termination
	for (const Interval &l : vL) {
		if (l < 0) {
			if (info) info->status = Info::Status::reachedTarget;
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
			if (info) info->status = Info::Status::completed;
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
			if (info) info->status = Info::Status::maxDepth;
			break;
		}

		// Subdomain contains invalidity
		if (dom.incl <= epsilon) {
			foundInvalid = true;
			if (info) info->status = Info::Status::reachedTarget;
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
Validity StaticValidator<n, s, p>::isValidMesh(
	span<const fp_t> cp,
	std::vector<uint> *adaptiveHierarchy,
	uint *invalidElemID,
	std::vector<uint> *invalidList
) const {
	const uint numCoordsPerElem = nControlGeoMap(n,s,p) * n;
	const uint numEl = cp.size() / (numCoordsPerElem);
	std::cerr << "number of polys: " << numEl << std::endl;
	std::vector<Validity> results(numEl);
	std::vector<fp_t> timings(numEl);
	std::vector<std::vector<uint>> hierarchies(numEl);
	bool gaveUp = false;
	bool foundInvalid = false;

	if (invalidList) invalidList->clear();

	#pragma omp parallel for num_threads(nThreads)
	for (uint e=0; e<numEl; ++e) {
		if (foundInvalid) {
			timings.at(e) = 0;
			continue;
		}
		span<const fp_t> element(
			cp.data() + numCoordsPerElem * e, numCoordsPerElem);
		Timer timer;
		timer.start();
		Validity v = isValidElement(element, &hierarchies.at(e));
		timer.stop();
		results.at(e) = v;
		if (v == Validity::invalid) foundInvalid = true;
		else if (v == Validity::uncertain) gaveUp = true;
		if (invalidList &&
			(v == Validity::invalid || v == Validity::uncertain)) {
			invalidList->push_back(e);
		}
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
