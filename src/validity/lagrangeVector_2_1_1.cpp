#include "lagrangeVector.hpp"

#define R(p, q) (Interval(p) / q)

namespace element_validity {
#ifdef LAGVEC_GCC_O0
#pragma GCC push_options
#pragma GCC optimize ("-O0")
#endif
template<>
void lagrangeVectorT<2, 1, 1>(const std::vector<fp_t> &cpFP, std::vector<Interval> &out) {
	out.resize(12);
	const uint S = cpFP.size();
	std::vector<Interval> cp(S);
	for (uint i = 0; i < S; ++i) cp[i] = cpFP[i];
	const Interval tmp_0 = cp[0]*cp[14];
	const Interval tmp_1 = cp[12]*cp[6];
	const Interval tmp_2 = cp[14]*cp[4];
	const Interval tmp_3 = cp[2]*cp[4];
	const Interval tmp_4 = -cp[0]*cp[6] + tmp_3;
	const Interval tmp_5 = (R(1, 4))*tmp_2;
	const Interval tmp_6 = (R(1, 4))*cp[14];
	const Interval tmp_7 = cp[5]*tmp_6;
	const Interval tmp_8 = (R(1, 4))*cp[4];
	const Interval tmp_9 = cp[15]*tmp_8;
	const Interval tmp_10 = cp[15]*cp[5];
	const Interval tmp_11 = (R(1, 4))*tmp_10;
	const Interval tmp_12 = (R(1, 4))*tmp_0;
	const Interval tmp_13 = (R(1, 4))*cp[0];
	const Interval tmp_14 = cp[15]*tmp_13;
	const Interval tmp_15 = (R(1, 4))*tmp_1;
	const Interval tmp_16 = (R(1, 4))*cp[12];
	const Interval tmp_17 = cp[7]*tmp_16;
	const Interval tmp_18 = (R(1, 4))*cp[13];
	const Interval tmp_19 = cp[6]*tmp_18;
	const Interval tmp_20 = cp[13]*cp[7];
	const Interval tmp_21 = (R(1, 4))*tmp_20;
	const Interval tmp_22 = cp[1]*tmp_6;
	const Interval tmp_23 = cp[15]*cp[1];
	const Interval tmp_24 = (R(1, 4))*tmp_23;
	const Interval tmp_25 = (R(1, 4))*cp[2];
	const Interval tmp_26 = cp[3]*cp[5];
	const Interval tmp_27 = -R(1, 4)*cp[0]*cp[6] - R(1, 4)*cp[0]*cp[7] - R(1, 4)*cp[1]*cp[6] - R(1, 4)*cp[1]*cp[7] + cp[3]*tmp_8 + cp[5]*tmp_25 + (R(1, 4))*tmp_26 + (R(1, 4))*tmp_3;
	const Interval tmp_28 = cp[13]*cp[3];
	const Interval tmp_29 = tmp_10 - tmp_20;
	const Interval tmp_30 = cp[12]*cp[2];
	const Interval tmp_31 = cp[14]*cp[8];
	const Interval tmp_32 = cp[10]*cp[12];
	const Interval tmp_33 = cp[0]*cp[10];
	const Interval tmp_34 = cp[2]*cp[8];
	const Interval tmp_35 = tmp_33 - tmp_34;
	const Interval tmp_36 = (R(1, 4))*tmp_32;
	const Interval tmp_37 = cp[10]*tmp_18;
	const Interval tmp_38 = cp[11]*tmp_16;
	const Interval tmp_39 = cp[11]*cp[13];
	const Interval tmp_40 = (R(1, 4))*tmp_39;
	const Interval tmp_41 = (R(1, 4))*cp[8];
	const Interval tmp_42 = cp[3]*cp[9];
	const Interval tmp_43 = cp[11]*cp[1];
	const Interval tmp_44 = (R(1, 4))*cp[10]*cp[1] + cp[11]*tmp_13 - cp[3]*tmp_41 - cp[9]*tmp_25 + (R(1, 4))*tmp_33 - R(1, 4)*tmp_34 - R(1, 4)*tmp_42 + (R(1, 4))*tmp_43;
	const Interval tmp_45 = -cp[15]*cp[9] + tmp_39;
	const Interval tmp_46 = cp[6]*cp[8];
	const Interval tmp_47 = -cp[10]*cp[4] + tmp_46;
	const Interval tmp_48 = cp[7]*cp[9];
	const Interval tmp_49 = -R(1, 4)*cp[10]*cp[4] - R(1, 4)*cp[10]*cp[5] - R(1, 4)*cp[11]*cp[4] - R(1, 4)*cp[11]*cp[5] + (R(1, 4))*cp[6]*cp[9] + cp[7]*tmp_41 + (R(1, 4))*tmp_46 + (R(1, 4))*tmp_48;
	const Interval tmp_50 = -cp[11]*cp[5] + tmp_48;
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
#ifdef LAGVEC_GCC_O0
#pragma GCC pop_options
#endif
#undef R