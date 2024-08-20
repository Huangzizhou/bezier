#include "../lagrangeEvaluate.hpp"

#define R(p, q) (Interval(p) / q)
#define I const Interval 

namespace element_validity {
template<>
Interval lagrangeEvaluate<2, 2, 4>(
	const std::span<const fp_t> xFP,
	const std::span<const Interval> lagVec
) {
	assert(xFP.size() == 2);
	assert(lagVec.size() == 15);
	std::array<Interval, 2> x;
	for (uint i = 0; i < 2; ++i) x[i] = xFP[i];
	Interval acc = 0.;
	I tmp_0 = x[0]*x[1];
	I tmp_1 = powi(x[0], 2);
	I tmp_2 = powi(x[0], 3);
	I tmp_3 = powi(x[0], 4);
	I tmp_4 = (R(32, 3))*tmp_3;
	I tmp_5 = powi(x[1], 2);
	I tmp_6 = powi(x[1], 3);
	I tmp_7 = powi(x[1], 4);
	I tmp_8 = (R(32, 3))*tmp_7;
	I tmp_9 = tmp_5*x[0];
	I tmp_10 = (R(128, 3))*tmp_6*x[0];
	I tmp_11 = tmp_1*x[1];
	I tmp_12 = (R(128, 3))*tmp_2*x[1];
	I tmp_13 = tmp_1*tmp_5;
	I tmp_14 = 64*tmp_13;
	I tmp_15 = 128*tmp_2;
	I tmp_16 = tmp_15*x[1];
	I tmp_17 = (R(128, 3))*tmp_3;
	I tmp_18 = 128*tmp_13;
	I tmp_19 = (R(208, 3))*tmp_0 + tmp_18;
	I tmp_20 = 28*tmp_0;
	I tmp_21 = tmp_14 - 16*tmp_5*x[0];
	I tmp_22 = (R(16, 3))*x[0]*x[1];
	I tmp_23 = -32*tmp_11;
	I tmp_24 = tmp_12 + tmp_22 + tmp_23;
	I tmp_25 = -16*tmp_1*x[1];
	I tmp_26 = -32*tmp_9;
	I tmp_27 = tmp_10 + tmp_22 + tmp_26;
	I tmp_28 = (R(128, 3))*tmp_7;
	I tmp_29 = 128*tmp_6;
	I tmp_30 = tmp_29*x[0];
	I tmp_31 = 32*tmp_0 + tmp_18;
	acc += lagVec[0] * (R(140, 3))*tmp_0 + (R(70, 3))*tmp_1 + tmp_10 - 80*tmp_11 + tmp_12 + tmp_14 - R(80, 3)*tmp_2 + tmp_4 + (R(70, 3))*tmp_5 - R(80, 3)*tmp_6 + tmp_8 - 80*tmp_9 - R(25, 3)*x[0] - R(25, 3)*x[1] + 1;
	acc += lagVec[1] * (R(22, 3))*tmp_1 - 16*tmp_2 + tmp_4 - x[0];
	acc += lagVec[2] * (R(22, 3))*tmp_5 - 16*tmp_6 + tmp_8 - x[1];
	acc += lagVec[3] * 192*tmp_1*x[1] - R(208, 3)*tmp_1 - tmp_10 - tmp_16 - tmp_17 - tmp_19 + 96*tmp_2 + 96*tmp_5*x[0] + 16*x[0];
	acc += lagVec[4] * 76*tmp_1 - 144*tmp_11 - tmp_15 + tmp_16 + tmp_20 + tmp_21 + 64*tmp_3 - 12*x[0];
	acc += lagVec[5] * -R(112, 3)*tmp_1 - tmp_17 + (R(224, 3))*tmp_2 - tmp_24 + (R(16, 3))*x[0];
	acc += lagVec[6] * tmp_24;
	acc += lagVec[7] * 4*tmp_0 + tmp_21 + tmp_25;
	acc += lagVec[8] * tmp_27;
	acc += lagVec[9] * -tmp_27 - tmp_28 - R(112, 3)*tmp_5 + (R(224, 3))*tmp_6 + (R(16, 3))*x[1];
	acc += lagVec[10] * tmp_14 + tmp_20 + tmp_25 - tmp_29 + tmp_30 + 76*tmp_5 + 64*tmp_7 - 144*tmp_9 - 12*x[1];
	acc += lagVec[11] * 96*tmp_1*x[1] - tmp_12 - tmp_19 - tmp_28 - tmp_30 + 192*tmp_5*x[0] - R(208, 3)*tmp_5 + 96*tmp_6 + 16*x[1];
	acc += lagVec[12] * 96*tmp_0 - 224*tmp_11 + 256*tmp_13 + tmp_16 + tmp_30 - 224*tmp_9;
	acc += lagVec[13] * -tmp_23 - tmp_30 - tmp_31 + 160*tmp_5*x[0];
	acc += lagVec[14] * 160*tmp_1*x[1] - tmp_16 - tmp_26 - tmp_31;
	return acc;
}}
#undef R
#undef I