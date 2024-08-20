#include "../lagrangeEvaluate.hpp"

#define R(p, q) (Interval(p) / q)
#define I const Interval 

namespace element_validity {
template<>
Interval lagrangeEvaluate<3, 3, 2>(
	const std::span<const fp_t> xFP,
	const std::span<const Interval> lagVec
) {
	assert(xFP.size() == 3);
	assert(lagVec.size() == 10);
	std::array<Interval, 3> x;
	for (uint i = 0; i < 3; ++i) x[i] = xFP[i];
	Interval acc = 0.;
	I tmp_0 = powi(x[0], 2);
	I tmp_1 = 2*tmp_0;
	I tmp_2 = powi(x[1], 2);
	I tmp_3 = 2*tmp_2;
	I tmp_4 = powi(x[2], 2);
	I tmp_5 = 2*tmp_4;
	I tmp_6 = 4*x[1]*x[2];
	I tmp_7 = 4*x[0];
	I tmp_8 = tmp_7*x[1];
	I tmp_9 = tmp_7*x[2];
	I tmp_10 = tmp_8 + tmp_9;
	acc += lagVec[0] * tmp_1 + tmp_10 + tmp_3 + tmp_5 + tmp_6 - 3*x[0] - 3*x[1] - 3*x[2] + 1;
	acc += lagVec[1] * tmp_1 - x[0];
	acc += lagVec[2] * tmp_3 - x[1];
	acc += lagVec[3] * tmp_5 - x[2];
	acc += lagVec[4] * -4*tmp_0 - tmp_10 + 4*x[0];
	acc += lagVec[5] * tmp_8;
	acc += lagVec[6] * -4*tmp_2 - tmp_6 - tmp_8 + 4*x[1];
	acc += lagVec[7] * -4*tmp_4 - tmp_6 - tmp_9 + 4*x[2];
	acc += lagVec[8] * tmp_9;
	acc += lagVec[9] * tmp_6;
	return acc;
}}
#undef R
#undef I