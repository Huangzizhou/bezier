#include "../lagrangeEvaluate.hpp"

#define R(p, q) (Interval(p) / q)
#define I const Interval 

namespace element_validity {
template<>
Interval lagrangeEvaluate<3, 3, 1>(
	const std::span<const fp_t> xFP,
	const std::span<const Interval> lagVec
) {
	assert(xFP.size() == 3);
	assert(lagVec.size() == 4);
	std::array<Interval, 3> x;
	for (uint i = 0; i < 3; ++i) x[i] = xFP[i];
	Interval acc = 0.;
	acc += lagVec[0] * -x[0] - x[1] - x[2] + 1;
	acc += lagVec[1] * x[0];
	acc += lagVec[2] * x[1];
	acc += lagVec[3] * x[2];
	return acc;
}}
#undef R
#undef I