#pragma once
#include "utils/combinatorics.hpp"
#include "lagrangeVector.hpp"
#include "lagrangeEvaluate.hpp"

template<uint n, uint s, uint p>
class JacobianEvaluator {
	private:
	std::vector<Interval> coeffs;
	constexpr uint numCoordsPerElem = nControlGeoMap(n,s,p) * n;
	const uint numEl = 0;

	public:
	JacobianEvaluator(const std::span<const fp_t>);
}

template<uint n, uint s, uint p>
JacobianEvaluator<n, s, p>(const std::span<const fp_t> cp) :
	numEl(cp.size() / (numCoordsPerElem)
) {
	const uint totCoord = numEl * numCoordsPerElem;
	assert(cp.size() == totCoord);
	const uint numLagCoefPerElem = nControlJacobian(n,s,p,false);
	coeffs.resize(totCoord);
	for (uint i=0; i < numEl; ++i) {
		std::span<const fp_t> in =
			cp.subspan(i * numCoordsPerElem, numCoordsPerElem);
		// Save Lagrange coefficients to coeffs vector
		std::span<const fp_t> out(
			coeffs.data() + i * numLagCoefPerElem, numLagCoefPerElem);
		lagrangeVector<n, s, p>(in, out);
	}
}

template<uint n, uint s, uint p>
fp_t JacobianEvaluator<n, s, p>::eval(uint i, std::span<const fp_t> x) const {
	assert(x.size() == n);
	std::span<const fp_t> lag(
		coeffs.data() + i * numLagCoefPerElem, numLagCoefPerElem);
	Interval res = lagrangeEvaluate(x, lag);
	return res.middle();
}
