#include "lagrangeVector.hpp"

#define R(p, q) (Interval(p) / q)
#define I const Interval 

namespace element_validity {
template<>
void lagrangeVectorT<3, 3, 1>(const std::span<const fp_t> cpFP, const std::span<Interval> out) {
	assert(out.size() == 4);
	const uint S = cpFP.size();
	std::vector<Interval> cp(S);
	for (uint i = 0; i < S; ++i) cp[i] = cpFP[i];
	I tmp_0 = cp[0] - cp[12];
	I tmp_1 = -cp[4];
	I tmp_2 = cp[10] + tmp_1;
	I tmp_3 = -cp[2];
	I tmp_4 = cp[20] + tmp_3;
	I tmp_5 = cp[14] + tmp_3;
	I tmp_6 = -cp[2] + cp[8];
	I tmp_7 = cp[0] - cp[18];
	I tmp_8 = cp[22] + tmp_1;
	I tmp_9 = cp[16] + tmp_1;
	I tmp_10 = (R(2, 3))*cp[0] + (R(1, 3))*cp[1];
	I tmp_11 = -R(2, 3)*cp[12] - R(1, 3)*cp[13] + tmp_10;
	I tmp_12 = -R(2, 3)*cp[4] - R(1, 3)*cp[5];
	I tmp_13 = (R(2, 3))*cp[10] + (R(1, 3))*cp[11] + tmp_12;
	I tmp_14 = (R(2, 3))*cp[2];
	I tmp_15 = (R(1, 3))*cp[3];
	I tmp_16 = -tmp_14 - tmp_15;
	I tmp_17 = (R(2, 3))*cp[20] + (R(1, 3))*cp[21] + tmp_16;
	I tmp_18 = (R(2, 3))*cp[14] + (R(1, 3))*cp[15] + tmp_16;
	I tmp_19 = (R(2, 3))*cp[8] + (R(1, 3))*cp[9] - tmp_14 - tmp_15;
	I tmp_20 = -R(2, 3)*cp[18] - R(1, 3)*cp[19] + tmp_10;
	I tmp_21 = (R(2, 3))*cp[22] + (R(1, 3))*cp[23] + tmp_12;
	I tmp_22 = (R(2, 3))*cp[16] + (R(1, 3))*cp[17] + tmp_12;
	I tmp_23 = (R(1, 3))*cp[0] + (R(2, 3))*cp[1];
	I tmp_24 = -R(1, 3)*cp[12] - R(2, 3)*cp[13] + tmp_23;
	I tmp_25 = -R(1, 3)*cp[4] - R(2, 3)*cp[5];
	I tmp_26 = (R(1, 3))*cp[10] + (R(2, 3))*cp[11] + tmp_25;
	I tmp_27 = (R(2, 3))*cp[3];
	I tmp_28 = (R(1, 3))*cp[2];
	I tmp_29 = -tmp_27 - tmp_28;
	I tmp_30 = (R(1, 3))*cp[20] + (R(2, 3))*cp[21] + tmp_29;
	I tmp_31 = (R(1, 3))*cp[14] + (R(2, 3))*cp[15] + tmp_29;
	I tmp_32 = (R(1, 3))*cp[8] + (R(2, 3))*cp[9] - tmp_27 - tmp_28;
	I tmp_33 = -R(1, 3)*cp[18] - R(2, 3)*cp[19] + tmp_23;
	I tmp_34 = (R(1, 3))*cp[22] + (R(2, 3))*cp[23] + tmp_25;
	I tmp_35 = (R(1, 3))*cp[16] + (R(2, 3))*cp[17] + tmp_25;
	I tmp_36 = -cp[1];
	I tmp_37 = -cp[13] - tmp_36;
	I tmp_38 = -cp[5];
	I tmp_39 = cp[11] + tmp_38;
	I tmp_40 = -cp[3];
	I tmp_41 = cp[21] + tmp_40;
	I tmp_42 = cp[15] + tmp_40;
	I tmp_43 = -cp[3] + cp[9];
	I tmp_44 = -cp[19] - tmp_36;
	I tmp_45 = cp[23] + tmp_38;
	I tmp_46 = cp[17] + tmp_38;
	out[0] = -tmp_0*(tmp_2*tmp_4 + tmp_5*tmp_6) - tmp_7*(tmp_2*tmp_8 + tmp_6*tmp_9) - (cp[0] - cp[6])*(-tmp_4*tmp_9 + tmp_5*tmp_8) - (-tmp_5 - tmp_8)*(tmp_0*tmp_6 + tmp_2*tmp_7);
	out[1] = -tmp_11*(tmp_13*tmp_17 + tmp_18*tmp_19) - tmp_20*(tmp_13*tmp_21 + tmp_19*tmp_22) - (-tmp_18 - tmp_21)*(tmp_11*tmp_19 + tmp_13*tmp_20) - (-tmp_17*tmp_22 + tmp_18*tmp_21)*(-R(2, 3)*cp[6] - R(1, 3)*cp[7] + tmp_10);
	out[2] = -tmp_24*(tmp_26*tmp_30 + tmp_31*tmp_32) - tmp_33*(tmp_26*tmp_34 + tmp_32*tmp_35) - (-tmp_31 - tmp_34)*(tmp_24*tmp_32 + tmp_26*tmp_33) - (-tmp_30*tmp_35 + tmp_31*tmp_34)*(-R(1, 3)*cp[6] - R(2, 3)*cp[7] + tmp_23);
	out[3] = -tmp_37*(tmp_39*tmp_41 + tmp_42*tmp_43) - tmp_44*(tmp_39*tmp_45 + tmp_43*tmp_46) - (cp[1] - cp[7])*(-tmp_41*tmp_46 + tmp_42*tmp_45) - (-tmp_42 - tmp_45)*(tmp_37*tmp_43 + tmp_39*tmp_44);
}}
#undef R
#undef I