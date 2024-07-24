#include "lagrangeVector.hpp"

#define R(p, q) (Interval(p) / q)
#define I const Interval 

namespace element_validity {
template<>
void lagrangeVectorT<2, 1, 1>(const std::span<const fp_t> cpFP, const std::span<Interval> out) {
	assert(out.size() == 12);
	const uint S = cpFP.size();
	std::vector<Interval> cp(S);
	for (uint i = 0; i < S; ++i) cp[i] = cpFP[i];
	I tmp_0 = cp[0]*cp[14];
	I tmp_1 = cp[12]*cp[6];
	I tmp_2 = cp[14]*cp[4];
	I tmp_3 = cp[2]*cp[4];
	I tmp_4 = -cp[0]*cp[6] + tmp_3;
	I tmp_5 = (R(1, 4))*tmp_2;
	I tmp_6 = (R(1, 4))*cp[14];
	I tmp_7 = cp[5]*tmp_6;
	I tmp_8 = (R(1, 4))*cp[4];
	I tmp_9 = cp[15]*tmp_8;
	I tmp_10 = cp[15]*cp[5];
	I tmp_11 = (R(1, 4))*tmp_10;
	I tmp_12 = (R(1, 4))*tmp_0;
	I tmp_13 = (R(1, 4))*cp[0];
	I tmp_14 = cp[15]*tmp_13;
	I tmp_15 = (R(1, 4))*tmp_1;
	I tmp_16 = (R(1, 4))*cp[12];
	I tmp_17 = cp[7]*tmp_16;
	I tmp_18 = (R(1, 4))*cp[13];
	I tmp_19 = cp[6]*tmp_18;
	I tmp_20 = cp[13]*cp[7];
	I tmp_21 = (R(1, 4))*tmp_20;
	I tmp_22 = cp[1]*tmp_6;
	I tmp_23 = cp[15]*cp[1];
	I tmp_24 = (R(1, 4))*tmp_23;
	I tmp_25 = (R(1, 4))*cp[2];
	I tmp_26 = cp[3]*cp[5];
	I tmp_27 = -R(1, 4)*cp[0]*cp[6] - R(1, 4)*cp[0]*cp[7] - R(1, 4)*cp[1]*cp[6] - R(1, 4)*cp[1]*cp[7] + cp[3]*tmp_8 + cp[5]*tmp_25 + (R(1, 4))*tmp_26 + (R(1, 4))*tmp_3;
	I tmp_28 = cp[13]*cp[3];
	I tmp_29 = tmp_10 - tmp_20;
	I tmp_30 = cp[12]*cp[2];
	I tmp_31 = cp[14]*cp[8];
	I tmp_32 = cp[10]*cp[12];
	I tmp_33 = cp[0]*cp[10];
	I tmp_34 = cp[2]*cp[8];
	I tmp_35 = tmp_33 - tmp_34;
	I tmp_36 = (R(1, 4))*tmp_32;
	I tmp_37 = cp[10]*tmp_18;
	I tmp_38 = cp[11]*tmp_16;
	I tmp_39 = cp[11]*cp[13];
	I tmp_40 = (R(1, 4))*tmp_39;
	I tmp_41 = (R(1, 4))*cp[8];
	I tmp_42 = cp[3]*cp[9];
	I tmp_43 = cp[11]*cp[1];
	I tmp_44 = (R(1, 4))*cp[10]*cp[1] + cp[11]*tmp_13 - cp[3]*tmp_41 - cp[9]*tmp_25 + (R(1, 4))*tmp_33 - R(1, 4)*tmp_34 - R(1, 4)*tmp_42 + (R(1, 4))*tmp_43;
	I tmp_45 = -cp[15]*cp[9] + tmp_39;
	I tmp_46 = cp[6]*cp[8];
	I tmp_47 = -cp[10]*cp[4] + tmp_46;
	I tmp_48 = cp[7]*cp[9];
	I tmp_49 = -R(1, 4)*cp[10]*cp[4] - R(1, 4)*cp[10]*cp[5] - R(1, 4)*cp[11]*cp[4] - R(1, 4)*cp[11]*cp[5] + (R(1, 4))*cp[6]*cp[9] + cp[7]*tmp_41 + (R(1, 4))*tmp_46 + (R(1, 4))*tmp_48;
	I tmp_50 = -cp[11]*cp[5] + tmp_48;
	out[0] = cp[12]*cp[2] - tmp_0 - tmp_1 + tmp_2 - tmp_4;
	out[1] = (R(1, 4))*cp[12]*cp[2] + (R(1, 4))*cp[12]*cp[3] + (R(1, 4))*cp[13]*cp[2] + (R(1, 4))*cp[13]*cp[3] + tmp_11 - tmp_12 - tmp_14 - tmp_15 - tmp_17 - tmp_19 - tmp_21 - tmp_22 - tmp_24 - tmp_27 + tmp_5 + tmp_7 + tmp_9;
	out[2] = cp[1]*cp[7] - tmp_23 - tmp_26 + tmp_28 + tmp_29;
	out[3] = -tmp_0 + tmp_30 + tmp_31 - tmp_32 + tmp_35;
	out[4] = (R(1, 4))*cp[15]*cp[9] + cp[15]*tmp_41 + cp[2]*tmp_18 + cp[3]*tmp_16 + cp[9]*tmp_6 - tmp_12 - tmp_14 - tmp_22 - tmp_24 + (R(1, 4))*tmp_28 + (R(1, 4))*tmp_30 + (R(1, 4))*tmp_31 - tmp_36 - tmp_37 - tmp_38 - tmp_40 + tmp_44;
	out[5] = cp[13]*cp[3] - tmp_23 - tmp_42 + tmp_43 - tmp_45;
	out[6] = -tmp_35 - tmp_4 - tmp_47;
	out[7] = -tmp_27 - tmp_44 - tmp_49;
	out[8] = cp[1]*cp[7] - tmp_26 + tmp_42 - tmp_43 - tmp_50;
	out[9] = cp[14]*cp[8] + tmp_1 - tmp_2 - tmp_32 - tmp_47;
	out[10] = (R(1, 4))*cp[14]*cp[8] + (R(1, 4))*cp[14]*cp[9] + (R(1, 4))*cp[15]*cp[8] + (R(1, 4))*cp[15]*cp[9] - tmp_11 + tmp_15 + tmp_17 + tmp_19 + tmp_21 - tmp_36 - tmp_37 - tmp_38 - tmp_40 - tmp_49 - tmp_5 - tmp_7 - tmp_9;
	out[11] = -tmp_29 - tmp_45 - tmp_50;
}}
#undef R
#undef I