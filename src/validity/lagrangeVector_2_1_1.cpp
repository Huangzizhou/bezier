#include "lagrangeVector.hpp"

#define R(p, q) (Interval(p) / q)

namespace element_validity {
template<>
void lagrangeVectorT<2, 1, 1>(const std::vector<fp_t> &cpFP, std::vector<Interval> &out) {
	out.resize(12);
	const uint S = cpFP.size();
	std::vector<Interval> cp(S);
	for (uint i = 0; i < S; ++i) cp[i] = cpFP[i];
	const Interval tmp_0 = cp[0]*cp[14];
	const Interval tmp_1 = cp[12]*cp[6];
	const Interval tmp_2 = cp[14]*cp[4];
	const Interval tmp_3 = -cp[0]*cp[6] + cp[2]*cp[4];
	const Interval tmp_4 = -cp[12]*cp[2] + tmp_0 + tmp_1 - tmp_2 + tmp_3;
	const Interval tmp_5 = cp[0]*cp[15];
	const Interval tmp_6 = cp[12]*cp[7];
	const Interval tmp_7 = cp[13]*cp[6];
	const Interval tmp_8 = cp[13]*cp[7];
	const Interval tmp_9 = cp[14]*cp[1];
	const Interval tmp_10 = cp[14]*cp[5];
	const Interval tmp_11 = cp[15]*cp[4];
	const Interval tmp_12 = cp[15]*cp[5];
	const Interval tmp_13 = cp[15]*cp[1];
	const Interval tmp_14 = -cp[13]*cp[3] + tmp_13;
	const Interval tmp_15 = cp[3]*cp[5];
	const Interval tmp_16 = -cp[1]*cp[7] + tmp_15;
	const Interval tmp_17 = -cp[0]*cp[7] - cp[1]*cp[6] + cp[2]*cp[5] + cp[3]*cp[4];
	const Interval tmp_18 = tmp_12 - tmp_8;
	const Interval tmp_19 = cp[13]*cp[3] - tmp_13;
	const Interval tmp_20 = cp[10]*cp[12];
	const Interval tmp_21 = cp[0]*cp[10] - cp[2]*cp[8];
	const Interval tmp_22 = cp[12]*cp[2] + cp[14]*cp[8] - tmp_0 - tmp_20 + tmp_21;
	const Interval tmp_23 = cp[10]*cp[13];
	const Interval tmp_24 = cp[11]*cp[12];
	const Interval tmp_25 = cp[11]*cp[13];
	const Interval tmp_26 = cp[11]*cp[1];
	const Interval tmp_27 = cp[3]*cp[9];
	const Interval tmp_28 = tmp_26 - tmp_27;
	const Interval tmp_29 = cp[0]*cp[11] + cp[10]*cp[1] - cp[2]*cp[9] - cp[3]*cp[8];
	const Interval tmp_30 = -cp[15]*cp[9] + tmp_25;
	const Interval tmp_31 = -cp[10]*cp[4] + cp[6]*cp[8];
	const Interval tmp_32 = tmp_21 + tmp_3 + tmp_31;
	const Interval tmp_33 = -cp[11]*cp[5] + cp[7]*cp[9];
	const Interval tmp_34 = tmp_16 + tmp_28 + tmp_33;
	const Interval tmp_35 = -cp[10]*cp[5] - cp[11]*cp[4] + cp[6]*cp[9] + cp[7]*cp[8];
	const Interval tmp_36 = -cp[14]*cp[8] - tmp_1 + tmp_2 + tmp_20 + tmp_31;
	const Interval tmp_37 = tmp_18 + tmp_30 + tmp_33;
	out[0] = -tmp_4;
	out[1] = (R(1, 4))*cp[12]*cp[3] + (R(1, 4))*cp[13]*cp[2] + (R(1, 4))*tmp_10 + (R(1, 4))*tmp_11 + (R(1, 4))*tmp_12 - R(1, 4)*tmp_14 - R(1, 4)*tmp_16 - R(1, 4)*tmp_17 - R(1, 4)*tmp_4 - R(1, 4)*tmp_5 - R(1, 4)*tmp_6 - R(1, 4)*tmp_7 - R(1, 4)*tmp_8 - R(1, 4)*tmp_9;
	out[2] = cp[1]*cp[7] - tmp_15 + tmp_18 + tmp_19;
	out[3] = tmp_22;
	out[4] = (R(1, 4))*cp[12]*cp[3] + (R(1, 4))*cp[13]*cp[2] + (R(1, 4))*cp[14]*cp[9] + (R(1, 4))*cp[15]*cp[8] + (R(1, 4))*cp[15]*cp[9] + (R(1, 4))*tmp_19 + (R(1, 4))*tmp_22 - R(1, 4)*tmp_23 - R(1, 4)*tmp_24 - R(1, 4)*tmp_25 + (R(1, 4))*tmp_28 + (R(1, 4))*tmp_29 - R(1, 4)*tmp_5 - R(1, 4)*tmp_9;
	out[5] = -tmp_14 + tmp_26 - tmp_27 - tmp_30;
	out[6] = -tmp_32;
	out[7] = -R(1, 4)*tmp_17 - R(1, 4)*tmp_29 - R(1, 4)*tmp_32 - R(1, 4)*tmp_34 - R(1, 4)*tmp_35;
	out[8] = -tmp_34;
	out[9] = -tmp_36;
	out[10] = (R(1, 4))*cp[14]*cp[9] + (R(1, 4))*cp[15]*cp[8] - R(1, 4)*tmp_10 - R(1, 4)*tmp_11 - R(1, 4)*tmp_23 - R(1, 4)*tmp_24 - R(1, 4)*tmp_35 - R(1, 4)*tmp_36 - R(1, 4)*tmp_37 + (R(1, 4))*tmp_6 + (R(1, 4))*tmp_7;
	out[11] = -tmp_37;
}}
