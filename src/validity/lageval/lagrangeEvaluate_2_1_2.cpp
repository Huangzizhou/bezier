#include "../lagrangeEvaluate.hpp"

#define R(p, q) (Interval(p) / q)
#define I const Interval 

namespace element_validity {
template<>
Interval lagrangeEvaluate<2, 1, 2>(
	const std::span<const fp_t> xFP,
	const std::span<const Interval> lagVec
) {
	assert(xFP.size() == 2);
	assert(lagVec.size() == 9);
	std::array<Interval, 2> x;
	for (uint i = 0; i < 2; ++i) x[i] = xFP[i];
	Interval acc = 0.;
	I tmp_0 = 3*x[0];
	I tmp_1 = x[0]*x[1];
	I tmp_2 = powi(x[0], 2);
	I tmp_3 = 2*tmp_2;
	I tmp_4 = tmp_2*x[1];
	I tmp_5 = tmp_3 - 6*tmp_4;
	I tmp_6 = powi(x[1], 2);
	I tmp_7 = 4*tmp_2;
	I tmp_8 = tmp_6*tmp_7;
	I tmp_9 = 2*tmp_6;
	I tmp_10 = tmp_6*x[0];
	I tmp_11 = -6*tmp_10 + tmp_8 + tmp_9;
	I tmp_12 = tmp_0*x[1];
	I tmp_13 = tmp_8 - tmp_9*x[0];
	I tmp_14 = -tmp_3*x[1];
	I tmp_15 = -8*tmp_6*x[0];
	I tmp_16 = tmp_2*tmp_6;
	I tmp_17 = 8*tmp_16;
	I tmp_18 = 12*tmp_1 + tmp_17;
	I tmp_19 = -8*tmp_2*x[1];
	I tmp_20 = 4*tmp_1 + tmp_17;
	acc += lagVec[0] * -tmp_0 + 9*tmp_1 + tmp_11 + tmp_5 - 3*x[1] + 1;
	acc += lagVec[1] * tmp_12 + tmp_13 + tmp_5 - x[0];
	acc += lagVec[2] * tmp_1 + tmp_13 + tmp_14;
	acc += lagVec[3] * tmp_11 + tmp_12 + tmp_14 - x[1];
	acc += lagVec[4] * -tmp_15 - tmp_18 + 12*tmp_2*x[1] - tmp_7 + 4*x[0];
	acc += lagVec[5] * -tmp_19 - tmp_20 + 4*tmp_6*x[0];
	acc += lagVec[6] * -tmp_15 + 4*tmp_2*x[1] - tmp_20;
	acc += lagVec[7] * -tmp_18 - tmp_19 + 12*tmp_6*x[0] - 4*tmp_6 + 4*x[1];
	acc += lagVec[8] * 16*tmp_1 - 16*tmp_10 + 16*tmp_16 - 16*tmp_4;
	return acc;
}}
#undef R
#undef I