#include "../lagrangeEvaluate.hpp"

#define R(p, q) (Interval(p) / q)
#define I const Interval 

namespace element_validity {
template<>
Interval lagrangeEvaluate<2, 2, 3>(
	const std::span<const fp_t> xFP,
	const std::span<const Interval> lagVec
) {
	assert(xFP.size() == 2);
	assert(lagVec.size() == 10);
	std::array<Interval, 2> x;
	for (uint i = 0; i < 2; ++i) x[i] = xFP[i];
	Interval acc = 0.;
	I tmp_0 = powi(x[0], 2);
	I tmp_1 = powi(x[0], 3);
	I tmp_2 = (R(9, 2))*tmp_1;
	I tmp_3 = powi(x[1], 2);
	I tmp_4 = powi(x[1], 3);
	I tmp_5 = (R(9, 2))*tmp_4;
	I tmp_6 = tmp_3*x[0];
	I tmp_7 = (R(27, 2))*tmp_6;
	I tmp_8 = tmp_0*x[1];
	I tmp_9 = (R(27, 2))*tmp_8;
	I tmp_10 = -R(45, 2)*x[0]*x[1];
	I tmp_11 = (R(27, 2))*tmp_1;
	I tmp_12 = 27*tmp_8;
	I tmp_13 = (R(9, 2))*x[0];
	I tmp_14 = -tmp_13*x[1];
	I tmp_15 = tmp_14 + tmp_9;
	I tmp_16 = tmp_14 + tmp_7;
	I tmp_17 = (R(27, 2))*tmp_4;
	I tmp_18 = 27*tmp_6;
	acc += lagVec[0] * 9*tmp_0 - tmp_2 + 9*tmp_3 - tmp_5 - tmp_7 - tmp_9 + 18*x[0]*x[1] - R(11, 2)*x[0] - R(11, 2)*x[1] + 1;
	acc += lagVec[1] * -R(9, 2)*tmp_0 + tmp_2 + x[0];
	acc += lagVec[2] * -R(9, 2)*tmp_3 + tmp_5 + x[1];
	acc += lagVec[3] * -R(45, 2)*tmp_0 + tmp_10 + tmp_11 + tmp_12 + tmp_7 + 9*x[0];
	acc += lagVec[4] * 18*tmp_0 - tmp_11 - tmp_13 - tmp_15;
	acc += lagVec[5] * tmp_15;
	acc += lagVec[6] * tmp_16;
	acc += lagVec[7] * -tmp_16 - tmp_17 + 18*tmp_3 - R(9, 2)*x[1];
	acc += lagVec[8] * tmp_10 + tmp_17 + tmp_18 - R(45, 2)*tmp_3 + tmp_9 + 9*x[1];
	acc += lagVec[9] * -tmp_12 - tmp_18 + 27*x[0]*x[1];
	return acc;
}}
#undef R
#undef I