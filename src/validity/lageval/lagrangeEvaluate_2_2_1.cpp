#include "../lagrangeEvaluate.hpp"

#define R(p, q) (Interval(p) / q)
#define I const Interval 

namespace element_validity {
template<>
Interval lagrangeEvaluate<2, 2, 1>(
	const std::span<const fp_t> xFP,
	const std::span<const Interval> lagVec
) {
	assert(xFP.size() == 2);
	assert(lagVec.size() == 3);
	std::array<Interval, 2> x;
	for (uint i = 0; i < 2; ++i) x[i] = xFP[i];
	Interval acc = 0.;
	acc += lagVec[0] * -x[0] - x[1] + 1;
	acc += lagVec[1] * x[0];
	acc += lagVec[2] * x[1];
	return acc;
}}
#undef R
#undef I