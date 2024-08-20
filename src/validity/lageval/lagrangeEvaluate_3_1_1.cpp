#include "../lagrangeEvaluate.hpp"

#define R(p, q) (Interval(p) / q)
#define I const Interval 

namespace element_validity {
template<>
Interval lagrangeEvaluate<3, 1, 1>(
	const std::span<const fp_t> xFP,
	const std::span<const Interval> lagVec
) {
	assert(xFP.size() == 3);
	assert(lagVec.size() == 8);
	std::array<Interval, 3> x;
	for (uint i = 0; i < 3; ++i) x[i] = xFP[i];
	Interval acc = 0.;
	I tmp_0 = -x[0]*x[2];
	I tmp_1 = x[0]*x[1];
	I tmp_2 = tmp_1*x[2];
	I tmp_3 = -tmp_1 + tmp_2;
	I tmp_4 = tmp_0 + tmp_3 + x[0];
	I tmp_5 = -x[1]*x[2];
	I tmp_6 = tmp_5 + x[1];
	I tmp_7 = tmp_0 + tmp_2;
	acc += lagVec[0] * -tmp_4 - tmp_6 - x[2] + 1;
	acc += lagVec[1] * tmp_4;
	acc += lagVec[2] * -tmp_3;
	acc += lagVec[3] * tmp_3 + tmp_6;
	acc += lagVec[4] * tmp_5 + tmp_7 + x[2];
	acc += lagVec[5] * -tmp_7;
	acc += lagVec[6] * tmp_2;
	acc += lagVec[7] * -tmp_2 - tmp_5;
	return acc;
}}
#undef R
#undef I