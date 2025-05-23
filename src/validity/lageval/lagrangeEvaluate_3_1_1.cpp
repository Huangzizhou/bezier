#include "../lagrangeEvaluate.hpp"

#define R(p, q) (Interval(p) / q)
#define I const Interval 

namespace element_validity {
template<>
Interval lagrangeEvaluate<3, 1, 1>(
	const span<const fp_t> xFP,
	const span<const Interval> lagVec
) {
	assert(xFP.size() == 3);
	assert(lagVec.size() == 27);
	std::array<Interval, 3> x;
	for (int i = 0; i < 3; ++i) x[i] = xFP[i];
	Interval acc = 0.;
	I tmp_0 = 3*x[0];
	I tmp_1 = 3*x[1];
	I tmp_2 = x[0]*x[1];
	I tmp_3 = x[0]*x[2];
	I tmp_4 = x[1]*x[2];
	I tmp_5 = 9*tmp_4;
	I tmp_6 = tmp_4*x[0];
	I tmp_7 = powi(x[2], 2);
	I tmp_8 = 2*tmp_7;
	I tmp_9 = 6*x[0];
	I tmp_10 = 6*x[1];
	I tmp_11 = tmp_2*tmp_7;
	I tmp_12 = powi(x[0], 2);
	I tmp_13 = 4*tmp_7;
	I tmp_14 = tmp_12*tmp_13;
	I tmp_15 = 12*tmp_7;
	I tmp_16 = -tmp_12*tmp_15*x[1] + tmp_14;
	I tmp_17 = -tmp_10*tmp_7 + 18*tmp_11 + tmp_16 - tmp_7*tmp_9 + tmp_8;
	I tmp_18 = powi(x[1], 2);
	I tmp_19 = 2*tmp_18;
	I tmp_20 = 6*x[2];
	I tmp_21 = tmp_18*tmp_3;
	I tmp_22 = tmp_13*tmp_18;
	I tmp_23 = -tmp_15*tmp_18*x[0] + tmp_22;
	I tmp_24 = -tmp_18*tmp_20 - tmp_18*tmp_9 + tmp_19 + 18*tmp_21 + tmp_23;
	I tmp_25 = 2*tmp_12;
	I tmp_26 = tmp_12*tmp_4;
	I tmp_27 = 8*tmp_7;
	I tmp_28 = tmp_12*tmp_27;
	I tmp_29 = tmp_18*tmp_28;
	I tmp_30 = 4*tmp_18;
	I tmp_31 = -12*tmp_12*tmp_18*x[2] + tmp_12*tmp_30 + tmp_29;
	I tmp_32 = -tmp_10*tmp_12 - tmp_12*tmp_20 + tmp_25 + 18*tmp_26 + tmp_31;
	I tmp_33 = -8*tmp_18*x[2];
	I tmp_34 = 16*x[2];
	I tmp_35 = tmp_18*tmp_34;
	I tmp_36 = -tmp_12*tmp_35;
	I tmp_37 = 16*tmp_7;
	I tmp_38 = tmp_18*tmp_37;
	I tmp_39 = tmp_12*tmp_38;
	I tmp_40 = tmp_36 + tmp_39;
	I tmp_41 = 24*tmp_26;
	I tmp_42 = -8*tmp_12*x[2] + tmp_41;
	I tmp_43 = 12*tmp_4;
	I tmp_44 = -36*x[0]*x[1]*x[2];
	I tmp_45 = tmp_18*tmp_27 - 24*tmp_18*tmp_7*x[0];
	I tmp_46 = tmp_43 + tmp_44 + tmp_45;
	I tmp_47 = 24*tmp_21;
	I tmp_48 = 12*tmp_3;
	I tmp_49 = -24*tmp_12*tmp_7*x[1] + tmp_28;
	I tmp_50 = tmp_47 + tmp_48 + tmp_49;
	I tmp_51 = tmp_0*x[2];
	I tmp_52 = 6*tmp_21;
	I tmp_53 = -tmp_19*x[2] + tmp_23 + tmp_52;
	I tmp_54 = 6*tmp_26;
	I tmp_55 = 4*x[2];
	I tmp_56 = tmp_18*tmp_55;
	I tmp_57 = -tmp_12*tmp_56 + tmp_29;
	I tmp_58 = -tmp_25*x[2] + tmp_54 + tmp_57;
	I tmp_59 = -tmp_5*x[0];
	I tmp_60 = tmp_1*x[2] + tmp_59;
	I tmp_61 = tmp_37*x[1];
	I tmp_62 = -tmp_12*tmp_61;
	I tmp_63 = -8*tmp_7*x[1];
	I tmp_64 = 8*tmp_12;
	I tmp_65 = -24*tmp_12*tmp_18*x[2] + tmp_18*tmp_64 + tmp_39;
	I tmp_66 = -8*tmp_12*x[1] + tmp_41 + tmp_65;
	I tmp_67 = 24*tmp_11;
	I tmp_68 = 12*tmp_2 + tmp_67;
	I tmp_69 = 32*tmp_26;
	I tmp_70 = 16*tmp_4;
	I tmp_71 = tmp_18*x[0];
	I tmp_72 = 48*tmp_7;
	I tmp_73 = tmp_12*tmp_18;
	I tmp_74 = tmp_73*x[2];
	I tmp_75 = -32*tmp_74;
	I tmp_76 = -48*tmp_6;
	I tmp_77 = 48*tmp_11 + tmp_75 + tmp_76;
	I tmp_78 = tmp_12*x[1];
	I tmp_79 = 32*tmp_7;
	I tmp_80 = tmp_73*tmp_79;
	I tmp_81 = -tmp_78*tmp_79 + tmp_80;
	I tmp_82 = 48*tmp_21 + tmp_81;
	I tmp_83 = -12*x[0]*x[1]*x[2];
	I tmp_84 = -8*tmp_12*tmp_18*x[2] + tmp_39;
	I tmp_85 = tmp_67 + tmp_83 + tmp_84;
	I tmp_86 = tmp_18*tmp_48 + tmp_62;
	I tmp_87 = tmp_4*tmp_64;
	I tmp_88 = tmp_55*x[1];
	I tmp_89 = tmp_45 + tmp_87 + tmp_88;
	I tmp_90 = -tmp_8*x[1];
	I tmp_91 = -tmp_14*x[1];
	I tmp_92 = -tmp_25*x[1] + tmp_31 + tmp_54 + tmp_91;
	I tmp_93 = 6*tmp_11;
	I tmp_94 = tmp_0*x[1] + tmp_93;
	I tmp_95 = -8*tmp_12*tmp_7*x[1];
	I tmp_96 = tmp_47 + tmp_83 + tmp_95;
	I tmp_97 = tmp_15*tmp_2 + tmp_40;
	I tmp_98 = tmp_25*tmp_4 + tmp_57 + tmp_91;
	I tmp_99 = -tmp_0*tmp_4;
	I tmp_100 = tmp_93 + tmp_99;
	I tmp_101 = -tmp_38*x[0];
	I tmp_102 = -8*tmp_7*x[0];
	I tmp_103 = -8*tmp_18*x[0] + tmp_65;
	I tmp_104 = 16*tmp_3;
	I tmp_105 = 32*tmp_21 + tmp_80;
	I tmp_106 = -tmp_71*tmp_79;
	I tmp_107 = tmp_106 + 48*tmp_26;
	I tmp_108 = tmp_101 + tmp_12*tmp_43;
	I tmp_109 = 8*tmp_21;
	I tmp_110 = tmp_109 + tmp_49 + tmp_55*x[0];
	I tmp_111 = 32*tmp_11;
	I tmp_112 = 16*tmp_12;
	I tmp_113 = -tmp_70*x[0];
	I tmp_114 = tmp_106 + tmp_113 + tmp_12*tmp_70;
	I tmp_115 = tmp_104*tmp_18 + tmp_81;
	I tmp_116 = tmp_2*tmp_27;
	I tmp_117 = tmp_116 + 4*x[0]*x[1];
	I tmp_118 = tmp_2*tmp_37 + tmp_75;
	I tmp_119 = -4*x[0]*x[1]*x[2];
	I tmp_120 = tmp_116 + tmp_119 + tmp_84;
	I tmp_121 = tmp_109 + tmp_95;
	I tmp_122 = tmp_16 - tmp_8*x[0];
	I tmp_123 = -tmp_22*x[0];
	I tmp_124 = tmp_123 - tmp_19*x[0] + tmp_52;
	I tmp_125 = -8*tmp_18*tmp_7*x[0];
	I tmp_126 = tmp_125 + tmp_83;
	I tmp_127 = tmp_123 + tmp_19*tmp_3;
	I tmp_128 = tmp_125 + tmp_87;
	I tmp_129 = tmp_2*tmp_8;
	acc += lagVec[0] * (-tmp_0 - tmp_1 + tmp_17 + 9*tmp_2 + tmp_24 + 9*tmp_3 + tmp_32 + tmp_5 - 27*tmp_6 - 3*x[2] + 1);
	acc += lagVec[1] * (-36*tmp_11 - tmp_13 - tmp_33 - tmp_40 - tmp_42 - tmp_46 - tmp_50 + 12*tmp_7*x[0] + 12*tmp_7*x[1] + 4*x[2]);
	acc += lagVec[2] * (tmp_17 + tmp_51 + tmp_53 + tmp_58 + tmp_60 - x[2]);
	acc += lagVec[3] * (12*tmp_18*x[0] + 12*tmp_18*x[2] - 36*tmp_21 - tmp_30 - tmp_46 - tmp_62 - tmp_63 - tmp_66 - tmp_68 + 4*x[1]);
	acc += lagVec[4] * (-tmp_35 + tmp_38 - tmp_61 + tmp_69 + tmp_70 - tmp_71*tmp_72 + tmp_77 + tmp_82);
	acc += lagVec[5] * (4*tmp_18*x[2] - tmp_63 - tmp_85 - tmp_86 - tmp_89);
	acc += lagVec[6] * (tmp_24 + tmp_60 + tmp_90 + tmp_92 + tmp_94 - x[1]);
	acc += lagVec[7] * (-tmp_33 + 4*tmp_7*x[1] - tmp_89 - tmp_96 - tmp_97);
	acc += lagVec[8] * (tmp_100 + tmp_4 + tmp_53 + tmp_90 + tmp_98);
	acc += lagVec[9] * (-tmp_101 - tmp_102 - tmp_103 + 12*tmp_12*x[1] + 12*tmp_12*x[2] - 4*tmp_12 - 36*tmp_26 - tmp_44 - tmp_50 - tmp_68 + 4*x[0]);
	acc += lagVec[10] * (tmp_104 + tmp_105 + tmp_107 - tmp_12*tmp_34 + tmp_12*tmp_37 - tmp_37*x[0] - tmp_72*tmp_78 + tmp_77);
	acc += lagVec[11] * (-tmp_102 - tmp_108 - tmp_110 + 4*tmp_12*x[2] - tmp_85);
	acc += lagVec[12] * (tmp_107 + tmp_111 + tmp_112*tmp_18 - tmp_112*x[1] + 16*tmp_2 - 16*tmp_71 - 48*tmp_74 + tmp_76 + tmp_82);
	acc += lagVec[13] * (-64*tmp_11 + 64*tmp_12*tmp_18*x[2] + 64*tmp_12*tmp_7*x[1] + 64*tmp_18*tmp_7*x[0] - 64*tmp_21 - 64*tmp_26 - 64*tmp_7*tmp_73 + 64*x[0]*x[1]*x[2]);
	acc += lagVec[14] * (tmp_111 + tmp_114 + tmp_115 + tmp_36);
	acc += lagVec[15] * (-tmp_103 - tmp_108 - tmp_117 + 4*tmp_12*x[1] - tmp_96);
	acc += lagVec[16] * (tmp_105 + tmp_114 + tmp_118 + tmp_62);
	acc += lagVec[17] * (-tmp_101 - tmp_12*tmp_88 - tmp_120 - tmp_121);
	acc += lagVec[18] * (tmp_122 + tmp_124 + tmp_32 + tmp_51 + tmp_59 + tmp_94 - x[0]);
	acc += lagVec[19] * (-tmp_110 - tmp_126 - tmp_42 + 4*tmp_7*x[0] - tmp_97);
	acc += lagVec[20] * (tmp_100 + tmp_122 + tmp_127 + tmp_3 + tmp_58);
	acc += lagVec[21] * (-tmp_117 - tmp_126 + 4*tmp_18*x[0] - tmp_66 - tmp_86);
	acc += lagVec[22] * (tmp_101 + tmp_113 + tmp_115 + tmp_118 + tmp_69);
	acc += lagVec[23] * (-tmp_120 - tmp_128 - tmp_56*x[0] - tmp_62);
	acc += lagVec[24] * (tmp_124 + tmp_129 + tmp_2 + tmp_92 + tmp_99);
	acc += lagVec[25] * (-tmp_119 - tmp_121 - tmp_128 - tmp_13*tmp_2 - tmp_40);
	acc += lagVec[26] * (tmp_127 + tmp_129 - tmp_6 + tmp_98);
	return acc;
}}
#undef R
#undef I