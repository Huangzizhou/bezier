#include "lagrangeVector.hpp"

#define R(p, q) (Interval(p) / q)

namespace element_validity {
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
	const Interval tmp_6 = cp[2] - cp[8];
	const Interval tmp_7 = cp[0] - cp[18];
	const Interval tmp_8 = cp[22] + tmp_1;
	const Interval tmp_9 = cp[16] + tmp_1;
	const Interval tmp_10 = cp[0] - cp[6];
	const Interval tmp_11 = tmp_5 + tmp_8;
	const Interval tmp_12 = 2*cp[0];
	const Interval tmp_13 = cp[1] + tmp_12;
	const Interval tmp_14 = -2*cp[12] - cp[13] + tmp_13;
	const Interval tmp_15 = -2*cp[4];
	const Interval tmp_16 = -cp[5];
	const Interval tmp_17 = cp[11] + tmp_16;
	const Interval tmp_18 = 2*cp[10] + tmp_15 + tmp_17;
	const Interval tmp_19 = 2*cp[2];
	const Interval tmp_20 = -tmp_19;
	const Interval tmp_21 = -cp[3];
	const Interval tmp_22 = cp[21] + tmp_21;
	const Interval tmp_23 = 2*cp[20] + tmp_20 + tmp_22;
	const Interval tmp_24 = cp[15] + tmp_21;
	const Interval tmp_25 = 2*cp[14] + tmp_20 + tmp_24;
	const Interval tmp_26 = cp[3] - cp[9];
	const Interval tmp_27 = -2*cp[8] + tmp_19 + tmp_26;
	const Interval tmp_28 = -2*cp[18] - cp[19] + tmp_13;
	const Interval tmp_29 = cp[23] + tmp_16;
	const Interval tmp_30 = 2*cp[22] + tmp_15 + tmp_29;
	const Interval tmp_31 = cp[17] + tmp_16;
	const Interval tmp_32 = 2*cp[16] + tmp_15 + tmp_31;
	const Interval tmp_33 = cp[1] - cp[7];
	const Interval tmp_34 = 2*cp[1];
	const Interval tmp_35 = -2*cp[13] + tmp_0 + tmp_34;
	const Interval tmp_36 = -2*cp[5];
	const Interval tmp_37 = 2*cp[11] + tmp_2 + tmp_36;
	const Interval tmp_38 = 2*cp[3];
	const Interval tmp_39 = -tmp_38;
	const Interval tmp_40 = 2*cp[21] + tmp_39 + tmp_4;
	const Interval tmp_41 = 2*cp[15] + tmp_39;
	const Interval tmp_42 = tmp_41 + tmp_5;
	const Interval tmp_43 = -2*cp[9] + tmp_38 + tmp_6;
	const Interval tmp_44 = -2*cp[19] + tmp_34 + tmp_7;
	const Interval tmp_45 = 2*cp[23] + tmp_36;
	const Interval tmp_46 = tmp_45 + tmp_8;
	const Interval tmp_47 = 2*cp[17] + tmp_36 + tmp_9;
	const Interval tmp_48 = -cp[1];
	const Interval tmp_49 = cp[13] + tmp_48;
	const Interval tmp_50 = cp[19] + tmp_48;
	out[0] = -tmp_0*(tmp_2*tmp_4 - tmp_5*tmp_6) - tmp_10*(-tmp_4*tmp_9 + tmp_5*tmp_8) + tmp_11*(-tmp_0*tmp_6 + tmp_2*tmp_7) - tmp_7*(tmp_2*tmp_8 - tmp_6*tmp_9);
	out[1] = -R(1, 27)*tmp_14*(tmp_18*tmp_23 - tmp_25*tmp_27) - R(1, 27)*tmp_28*(tmp_18*tmp_30 - tmp_27*tmp_32) + (R(1, 27))*(tmp_25 + tmp_30)*(-tmp_14*tmp_27 + tmp_18*tmp_28) - R(1, 27)*(-tmp_23*tmp_32 + tmp_25*tmp_30)*(-2*cp[6] + tmp_12 + tmp_33);
	out[2] = -R(1, 27)*tmp_35*(tmp_37*tmp_40 - tmp_42*tmp_43) - R(1, 27)*tmp_44*(tmp_37*tmp_46 - tmp_43*tmp_47) + (R(1, 27))*(-tmp_35*tmp_43 + tmp_37*tmp_44)*(tmp_11 + tmp_41 + tmp_45) - R(1, 27)*(-tmp_40*tmp_47 + tmp_42*tmp_46)*(-2*cp[7] + tmp_10 + tmp_34);
	out[3] = -tmp_33*(-tmp_22*tmp_31 + tmp_24*tmp_29) + tmp_49*(tmp_17*tmp_22 - tmp_24*tmp_26) + tmp_50*(tmp_17*tmp_29 - tmp_26*tmp_31) + (tmp_24 + tmp_29)*(-tmp_17*tmp_50 + tmp_26*tmp_49);
}}
