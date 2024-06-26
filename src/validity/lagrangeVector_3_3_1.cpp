#include "lagrangeVector.hpp"

#define R(p, q) (Interval(p) / q)

namespace element_validity {
#ifdef LAGVEC_GCC_O0
#pragma GCC push_options
#pragma GCC optimize ("-O0")
#endif
template<>
void lagrangeVectorT<3, 3, 1>(const std::vector<fp_t> &cpFP, std::vector<Interval> &out) {
	out.resize(4);
	const uint S = cpFP.size();
	std::vector<Interval> cp(S);
	for (uint i = 0; i < S; ++i) cp[i] = cpFP[i];
	const Interval tmp_0 = cp[0] - cp[12];
	const Interval tmp_1 = -cp[4];
	const Interval tmp_2 = cp[10] + tmp_1;
	const Interval tmp_3 = -cp[2];
	const Interval tmp_4 = cp[20] + tmp_3;
	const Interval tmp_5 = cp[14] + tmp_3;
	const Interval tmp_6 = -cp[2] + cp[8];
	const Interval tmp_7 = cp[0] - cp[18];
	const Interval tmp_8 = cp[22] + tmp_1;
	const Interval tmp_9 = cp[16] + tmp_1;
	const Interval tmp_10 = (R(2, 3))*cp[0] + (R(1, 3))*cp[1];
	const Interval tmp_11 = -R(2, 3)*cp[12] - R(1, 3)*cp[13] + tmp_10;
	const Interval tmp_12 = -R(2, 3)*cp[4] - R(1, 3)*cp[5];
	const Interval tmp_13 = (R(2, 3))*cp[10] + (R(1, 3))*cp[11] + tmp_12;
	const Interval tmp_14 = (R(2, 3))*cp[2];
	const Interval tmp_15 = (R(1, 3))*cp[3];
	const Interval tmp_16 = -tmp_14 - tmp_15;
	const Interval tmp_17 = (R(2, 3))*cp[20] + (R(1, 3))*cp[21] + tmp_16;
	const Interval tmp_18 = (R(2, 3))*cp[14] + (R(1, 3))*cp[15] + tmp_16;
	const Interval tmp_19 = (R(2, 3))*cp[8] + (R(1, 3))*cp[9] - tmp_14 - tmp_15;
	const Interval tmp_20 = -R(2, 3)*cp[18] - R(1, 3)*cp[19] + tmp_10;
	const Interval tmp_21 = (R(2, 3))*cp[22] + (R(1, 3))*cp[23] + tmp_12;
	const Interval tmp_22 = (R(2, 3))*cp[16] + (R(1, 3))*cp[17] + tmp_12;
	const Interval tmp_23 = (R(1, 3))*cp[0] + (R(2, 3))*cp[1];
	const Interval tmp_24 = -R(1, 3)*cp[12] - R(2, 3)*cp[13] + tmp_23;
	const Interval tmp_25 = -R(1, 3)*cp[4] - R(2, 3)*cp[5];
	const Interval tmp_26 = (R(1, 3))*cp[10] + (R(2, 3))*cp[11] + tmp_25;
	const Interval tmp_27 = (R(2, 3))*cp[3];
	const Interval tmp_28 = (R(1, 3))*cp[2];
	const Interval tmp_29 = -tmp_27 - tmp_28;
	const Interval tmp_30 = (R(1, 3))*cp[20] + (R(2, 3))*cp[21] + tmp_29;
	const Interval tmp_31 = (R(1, 3))*cp[14] + (R(2, 3))*cp[15] + tmp_29;
	const Interval tmp_32 = (R(1, 3))*cp[8] + (R(2, 3))*cp[9] - tmp_27 - tmp_28;
	const Interval tmp_33 = -R(1, 3)*cp[18] - R(2, 3)*cp[19] + tmp_23;
	const Interval tmp_34 = (R(1, 3))*cp[22] + (R(2, 3))*cp[23] + tmp_25;
	const Interval tmp_35 = (R(1, 3))*cp[16] + (R(2, 3))*cp[17] + tmp_25;
	const Interval tmp_36 = -cp[1];
	const Interval tmp_37 = -cp[13] - tmp_36;
	const Interval tmp_38 = -cp[5];
	const Interval tmp_39 = cp[11] + tmp_38;
	const Interval tmp_40 = -cp[3];
	const Interval tmp_41 = cp[21] + tmp_40;
	const Interval tmp_42 = cp[15] + tmp_40;
	const Interval tmp_43 = -cp[3] + cp[9];
	const Interval tmp_44 = -cp[19] - tmp_36;
	const Interval tmp_45 = cp[23] + tmp_38;
	const Interval tmp_46 = cp[17] + tmp_38;
	out[0] = -tmp_0*(tmp_2*tmp_4 + tmp_5*tmp_6) - tmp_7*(tmp_2*tmp_8 + tmp_6*tmp_9) - (cp[0] - cp[6])*(-tmp_4*tmp_9 + tmp_5*tmp_8) - (-tmp_5 - tmp_8)*(tmp_0*tmp_6 + tmp_2*tmp_7);
	out[1] = -tmp_11*(tmp_13*tmp_17 + tmp_18*tmp_19) - tmp_20*(tmp_13*tmp_21 + tmp_19*tmp_22) - (-tmp_18 - tmp_21)*(tmp_11*tmp_19 + tmp_13*tmp_20) - (-tmp_17*tmp_22 + tmp_18*tmp_21)*(-R(2, 3)*cp[6] - R(1, 3)*cp[7] + tmp_10);
	out[2] = -tmp_24*(tmp_26*tmp_30 + tmp_31*tmp_32) - tmp_33*(tmp_26*tmp_34 + tmp_32*tmp_35) - (-tmp_31 - tmp_34)*(tmp_24*tmp_32 + tmp_26*tmp_33) - (-tmp_30*tmp_35 + tmp_31*tmp_34)*(-R(1, 3)*cp[6] - R(2, 3)*cp[7] + tmp_23);
	out[3] = -tmp_37*(tmp_39*tmp_41 + tmp_42*tmp_43) - tmp_44*(tmp_39*tmp_45 + tmp_43*tmp_46) - (cp[1] - cp[7])*(-tmp_41*tmp_46 + tmp_42*tmp_45) - (-tmp_42 - tmp_45)*(tmp_37*tmp_43 + tmp_39*tmp_44);
}}
#ifdef LAGVEC_GCC_O0
#pragma GCC pop_options
#endif
#undef R