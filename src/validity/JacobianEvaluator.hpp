#pragma once
#include "numbers/Interval.hpp"
#include "utils/combinatorics.hpp"
#include "lagrangeVector.hpp"
#include "lagrangeEvaluate.hpp"

namespace element_validity {
template<uint n, uint s, uint p>
class JacobianEvaluator {
	private:
	std::vector<Interval> coeffs;
	static constexpr uint numCoordsPerElem = nControlGeoMap(n,s,p) * n;
	static constexpr uint numLagCoefPerElem = nControlJacobian(n,s,p,false);
	const uint numEl = 0;

	public:
	JacobianEvaluator(const span<const fp_t>);
	fp_t eval(span<const fp_t>, uint element=0) const;
	fp_t eval(std::initializer_list<fp_t> il, uint element=0) const {
		std::vector<fp_t> v(il);
		return eval(v, element);
	}

	#ifdef EIGEN_INTERFACE
	JacobianEvaluator(const Eigen::MatrixXd& cp) :
		JacobianEvaluator(convertEigenMatrix(cp)) {}
	Eigen::VectorXd eval(const Eigen::MatrixXd& uv, uint element=0) const {
		std::vector<double> uvs = convertEigenMatrix(uv);
		const uint numPts = uvs.size() / n;
		Eigen::VectorXd res(numPts);
		for (uint i=0; i<numPts; ++i) {
			span<fp_t> pt(uvs.data() + i*n, n);
			res[i] = eval(pt);
		}
		return res;
	}
	#endif
};

template<uint n, uint s, uint p>
JacobianEvaluator<n, s, p>::JacobianEvaluator(const span<const fp_t> cp) :
	numEl(cp.size() / (numCoordsPerElem)
) {
	const uint totCoord = numEl * numCoordsPerElem;
	assert(cp.size() == totCoord);
	coeffs.resize(totCoord);
	for (uint i=0; i < numEl; ++i) {
		span<const fp_t> in =
			cp.subspan(i * numCoordsPerElem, numCoordsPerElem);
		// Save Lagrange coefficients to coeffs vector
		span<Interval> out(
			coeffs.data() + i * numLagCoefPerElem, numLagCoefPerElem);
		lagrangeVector<n, s, p>(in, out);
	}
}

template<uint n, uint s, uint p>
fp_t JacobianEvaluator<n, s, p>::eval(span<const fp_t> x, uint i) const {
	assert(x.size() == n);
	span<const Interval> lag(
		coeffs.data() + i * numLagCoefPerElem, numLagCoefPerElem);
	Interval res = lagrangeEvaluate<n,s,p>(x, lag);
	return res.middle();
}
}