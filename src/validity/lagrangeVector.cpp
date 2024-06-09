#include "lagrangeVector.hpp"

#define R(p, q) Interval(p) / q

namespace element_validity {
template<>
void lagrangeVectorT<1, 1, 1>(const std::vector<fp_t> &cpFP, std::vector<Interval> &out) {
	out.resize(2);
	std::vector<Interval> cp(cpFP.size());
	for (uint i = 0; i < cpFP.size(); ++i) cp[i] = cpFP[i];
	out[0] = -cp[0] + cp[2];
	out[1] = -cp[1] + cp[3];
}

template<>
void lagrangeVectorT<1, 1, 2>(const std::vector<fp_t> &cpFP, std::vector<Interval> &out) {
	out.resize(4);
	std::vector<Interval> cp(cpFP.size());
	for (uint i = 0; i < cpFP.size(); ++i) cp[i] = cpFP[i];
	const Interval tmp_0 = -cp[2];
	const Interval tmp_1 = -cp[3];
	out[0] = -3*cp[0] - cp[4] - tmp_0;
	out[1] = -3*cp[1] - cp[5] - tmp_1;
	out[2] = cp[0] + 3*cp[4] + tmp_0;
	out[3] = cp[1] + 3*cp[5] + tmp_1;
}

template<>
void lagrangeVectorT<1, 1, 3>(const std::vector<fp_t> &cpFP, std::vector<Interval> &out) {
	out.resize(6);
	std::vector<Interval> cp(cpFP.size());
	for (uint i = 0; i < cpFP.size(); ++i) cp[i] = cpFP[i];
	const Interval tmp_0 = -cp[6];
	const Interval tmp_1 = -cp[2] + cp[4];
	const Interval tmp_2 = -cp[7];
	const Interval tmp_3 = -cp[3] + cp[5];
	out[0] = -R(11, 2)*cp[0] - tmp_0 - tmp_1;
	out[1] = -R(11, 2)*cp[1] - tmp_2 - tmp_3;
	out[2] = (R(1, 8))*cp[0] - R(1, 2)*cp[2] + (R(1, 2))*cp[4] + (R(1, 8))*tmp_0;
	out[3] = (R(1, 8))*cp[1] - R(1, 2)*cp[3] + (R(1, 2))*cp[5] + (R(1, 8))*tmp_2;
	out[4] = -cp[0] + (R(11, 2))*cp[6] - tmp_1;
	out[5] = -cp[1] + (R(11, 2))*cp[7] - tmp_3;
}

template<>
void lagrangeVectorT<2, 2, 1>(const std::vector<fp_t> &cpFP, std::vector<Interval> &out) {
	out.resize(3);
	std::vector<Interval> cp(cpFP.size());
	for (uint i = 0; i < cpFP.size(); ++i) cp[i] = cpFP[i];
	const Interval tmp_0 = cp[0]*cp[10] - cp[0]*cp[6] - cp[10]*cp[4] + cp[2]*cp[4] - cp[2]*cp[8] + cp[6]*cp[8];
	const Interval tmp_1 = cp[11]*cp[1] - cp[11]*cp[5] - cp[1]*cp[7] + cp[3]*cp[5] - cp[3]*cp[9] + cp[7]*cp[9];
	out[0] = tmp_0;
	out[1] = (R(1, 4))*cp[0]*cp[11] - R(1, 4)*cp[0]*cp[7] + (R(1, 4))*cp[10]*cp[1] - R(1, 4)*cp[10]*cp[5] - R(1, 4)*cp[11]*cp[4] - R(1, 4)*cp[1]*cp[6] + (R(1, 4))*cp[2]*cp[5] - R(1, 4)*cp[2]*cp[9] + (R(1, 4))*cp[3]*cp[4] - R(1, 4)*cp[3]*cp[8] + (R(1, 4))*cp[6]*cp[9] + (R(1, 4))*cp[7]*cp[8] + (R(1, 4))*tmp_0 + (R(1, 4))*tmp_1;
	out[2] = tmp_1;
}

template<>
void lagrangeVectorT<2, 2, 2>(const std::vector<fp_t> &cpFP, std::vector<Interval> &out) {
	out.resize(18);
	std::vector<Interval> cp(cpFP.size());
	for (uint i = 0; i < cpFP.size(); ++i) cp[i] = cpFP[i];
	const Interval tmp_0 = cp[0]*cp[22];
	const Interval tmp_1 = 3*tmp_0;
	const Interval tmp_2 = cp[0]*cp[6];
	const Interval tmp_3 = cp[12]*cp[2];
	const Interval tmp_4 = cp[0]*cp[14];
	const Interval tmp_5 = cp[20]*cp[2];
	const Interval tmp_6 = cp[2]*cp[4];
	const Interval tmp_7 = cp[2]*cp[8];
	const Interval tmp_8 = cp[0]*cp[10];
	const Interval tmp_9 = -3*tmp_7 + 3*tmp_8;
	const Interval tmp_10 = cp[10]*cp[20];
	const Interval tmp_11 = cp[22]*cp[8];
	const Interval tmp_12 = tmp_10 - tmp_11;
	const Interval tmp_13 = cp[22]*cp[4];
	const Interval tmp_14 = cp[20]*cp[6];
	const Interval tmp_15 = cp[14]*cp[8];
	const Interval tmp_16 = cp[10]*cp[12];
	const Interval tmp_17 = tmp_15 - tmp_16;
	const Interval tmp_18 = tmp_13 - tmp_14 + tmp_17;
	const Interval tmp_19 = cp[12]*cp[6] - cp[14]*cp[4] - tmp_1 + tmp_12 + tmp_18 - 3*tmp_2 - 3*tmp_3 + 3*tmp_4 + 3*tmp_5 + 3*tmp_6 + tmp_9;
	const Interval tmp_20 = cp[10]*cp[21];
	const Interval tmp_21 = cp[11]*cp[20];
	const Interval tmp_22 = cp[11]*cp[21];
	const Interval tmp_23 = cp[15]*cp[5];
	const Interval tmp_24 = cp[22]*cp[9];
	const Interval tmp_25 = cp[23]*cp[8];
	const Interval tmp_26 = cp[23]*cp[9];
	const Interval tmp_27 = cp[0]*cp[23];
	const Interval tmp_28 = 3*tmp_27;
	const Interval tmp_29 = cp[0]*cp[7];
	const Interval tmp_30 = cp[12]*cp[3];
	const Interval tmp_31 = cp[13]*cp[2];
	const Interval tmp_32 = cp[13]*cp[3];
	const Interval tmp_33 = 3*tmp_32;
	const Interval tmp_34 = cp[1]*cp[22];
	const Interval tmp_35 = 3*tmp_34;
	const Interval tmp_36 = cp[1]*cp[6];
	const Interval tmp_37 = cp[1]*cp[7];
	const Interval tmp_38 = 3*tmp_37;
	const Interval tmp_39 = cp[0]*cp[15];
	const Interval tmp_40 = cp[14]*cp[1];
	const Interval tmp_41 = cp[15]*cp[1];
	const Interval tmp_42 = cp[20]*cp[3];
	const Interval tmp_43 = cp[21]*cp[2];
	const Interval tmp_44 = cp[2]*cp[5];
	const Interval tmp_45 = cp[3]*cp[4];
	const Interval tmp_46 = cp[3]*cp[5];
	const Interval tmp_47 = cp[3]*cp[9];
	const Interval tmp_48 = 3*tmp_47;
	const Interval tmp_49 = cp[11]*cp[1];
	const Interval tmp_50 = -tmp_48 + 3*tmp_49;
	const Interval tmp_51 = cp[1]*cp[23];
	const Interval tmp_52 = 3*tmp_51;
	const Interval tmp_53 = cp[21]*cp[3];
	const Interval tmp_54 = -tmp_52 + 3*tmp_53;
	const Interval tmp_55 = cp[2]*cp[9];
	const Interval tmp_56 = cp[3]*cp[8];
	const Interval tmp_57 = cp[0]*cp[11];
	const Interval tmp_58 = cp[10]*cp[1];
	const Interval tmp_59 = -3*tmp_55 - 3*tmp_56 + 3*tmp_57 + 3*tmp_58;
	const Interval tmp_60 = cp[22]*cp[5];
	const Interval tmp_61 = cp[23]*cp[4];
	const Interval tmp_62 = cp[23]*cp[5];
	const Interval tmp_63 = cp[20]*cp[7];
	const Interval tmp_64 = cp[21]*cp[6];
	const Interval tmp_65 = cp[21]*cp[7];
	const Interval tmp_66 = cp[15]*cp[9];
	const Interval tmp_67 = cp[11]*cp[13];
	const Interval tmp_68 = -cp[10]*cp[13] - cp[11]*cp[12] + cp[14]*cp[9] + cp[15]*cp[8] + tmp_66 - tmp_67;
	const Interval tmp_69 = tmp_60 + tmp_61 + tmp_62 - tmp_63 - tmp_64 - tmp_65 + tmp_68;
	const Interval tmp_70 = -tmp_62 + tmp_65;
	const Interval tmp_71 = -tmp_66 + tmp_67;
	const Interval tmp_72 = -tmp_22 + tmp_26;
	const Interval tmp_73 = -3*cp[21]*cp[3] + tmp_52;
	const Interval tmp_74 = (R(1, 2))*tmp_2;
	const Interval tmp_75 = cp[10]*cp[4];
	const Interval tmp_76 = (R(1, 2))*tmp_75;
	const Interval tmp_77 = (R(1, 2))*tmp_3;
	const Interval tmp_78 = cp[18]*cp[8];
	const Interval tmp_79 = (R(1, 2))*tmp_78;
	const Interval tmp_80 = cp[16]*cp[2];
	const Interval tmp_81 = cp[0]*cp[18];
	const Interval tmp_82 = tmp_0 - tmp_5;
	const Interval tmp_83 = -tmp_7 + tmp_8;
	const Interval tmp_84 = (R(1, 2))*tmp_80 - R(1, 2)*tmp_81 + tmp_82 + tmp_83;
	const Interval tmp_85 = tmp_12 + (R(1, 2))*tmp_15 - R(1, 2)*tmp_16;
	const Interval tmp_86 = cp[10]*cp[5];
	const Interval tmp_87 = cp[11]*cp[4];
	const Interval tmp_88 = cp[11]*cp[5];
	const Interval tmp_89 = cp[18]*cp[9];
	const Interval tmp_90 = cp[19]*cp[8];
	const Interval tmp_91 = cp[19]*cp[9];
	const Interval tmp_92 = tmp_80 - tmp_81;
	const Interval tmp_93 = tmp_17 + tmp_92;
	const Interval tmp_94 = cp[17]*cp[3];
	const Interval tmp_95 = cp[19]*cp[1];
	const Interval tmp_96 = tmp_94 - tmp_95;
	const Interval tmp_97 = -cp[0]*cp[19] + cp[16]*cp[3] + cp[17]*cp[2] - cp[18]*cp[1];
	const Interval tmp_98 = tmp_68 + tmp_96 + tmp_97;
	const Interval tmp_99 = 2*tmp_5;
	const Interval tmp_100 = 2*tmp_42;
	const Interval tmp_101 = 2*tmp_43;
	const Interval tmp_102 = 2*tmp_53;
	const Interval tmp_103 = -2*tmp_7;
	const Interval tmp_104 = -2*tmp_55;
	const Interval tmp_105 = -2*tmp_56;
	const Interval tmp_106 = -2*tmp_47;
	const Interval tmp_107 = 2*tmp_8;
	const Interval tmp_108 = 2*tmp_57;
	const Interval tmp_109 = 2*tmp_0;
	const Interval tmp_110 = 2*tmp_27;
	const Interval tmp_111 = 2*tmp_58;
	const Interval tmp_112 = 2*tmp_49;
	const Interval tmp_113 = 2*tmp_34;
	const Interval tmp_114 = 2*tmp_51;
	const Interval tmp_115 = -tmp_100 - tmp_101 - tmp_102 + tmp_103 + tmp_104 + tmp_105 + tmp_106 + tmp_107 + tmp_108 + tmp_109 + tmp_110 + tmp_111 + tmp_112 + tmp_113 + tmp_114 - tmp_99;
	const Interval tmp_116 = 2*tmp_11;
	const Interval tmp_117 = 2*tmp_24;
	const Interval tmp_118 = 2*tmp_25;
	const Interval tmp_119 = 2*tmp_26;
	const Interval tmp_120 = 2*tmp_10;
	const Interval tmp_121 = 2*tmp_20;
	const Interval tmp_122 = 2*tmp_21;
	const Interval tmp_123 = 2*tmp_22;
	const Interval tmp_124 = -tmp_116 - tmp_117 - tmp_118 - tmp_119 + tmp_120 + tmp_121 + tmp_122 + tmp_123;
	const Interval tmp_125 = cp[11]*cp[17];
	const Interval tmp_126 = (R(1, 2))*tmp_95;
	const Interval tmp_127 = cp[7]*cp[9];
	const Interval tmp_128 = (R(1, 2))*tmp_88;
	const Interval tmp_129 = (R(1, 2))*tmp_94;
	const Interval tmp_130 = (R(1, 2))*tmp_91;
	const Interval tmp_131 = tmp_47 - tmp_49;
	const Interval tmp_132 = tmp_131 - R(1, 2)*tmp_66 + (R(1, 2))*tmp_67 + tmp_72;
	const Interval tmp_133 = -R(1, 2)*tmp_32 - R(1, 2)*tmp_37 + (R(1, 2))*tmp_41 + (R(1, 2))*tmp_46;
	const Interval tmp_134 = 3*tmp_10;
	const Interval tmp_135 = cp[10]*cp[16];
	const Interval tmp_136 = cp[6]*cp[8];
	const Interval tmp_137 = -tmp_13 + tmp_14 + tmp_92;
	const Interval tmp_138 = -cp[16]*cp[6] + cp[18]*cp[4] + 3*tmp_11 - tmp_134 + 3*tmp_135 + 3*tmp_136 + tmp_137 - 3*tmp_75 - 3*tmp_78 + tmp_82 + tmp_9;
	const Interval tmp_139 = 3*tmp_20;
	const Interval tmp_140 = 3*tmp_21;
	const Interval tmp_141 = cp[10]*cp[17];
	const Interval tmp_142 = cp[11]*cp[16];
	const Interval tmp_143 = cp[6]*cp[9];
	const Interval tmp_144 = cp[7]*cp[8];
	const Interval tmp_145 = 3*tmp_22;
	const Interval tmp_146 = -tmp_145 + 3*tmp_26;
	const Interval tmp_147 = tmp_51 - tmp_53;
	const Interval tmp_148 = tmp_70 + tmp_96;
	const Interval tmp_149 = -cp[17]*cp[7] + cp[19]*cp[5] + 3*tmp_125 + 3*tmp_127 + tmp_146 + tmp_147 + tmp_148 + tmp_50 - 3*tmp_88 - 3*tmp_91;
	const Interval tmp_150 = -tmp_60 - tmp_61 + tmp_63 + tmp_64 + tmp_97;
	const Interval tmp_151 = (R(1, 2))*tmp_14;
	const Interval tmp_152 = (R(1, 2))*tmp_13;
	const Interval tmp_153 = cp[12]*cp[22];
	const Interval tmp_154 = cp[18]*cp[20];
	const Interval tmp_155 = cp[14]*cp[20];
	const Interval tmp_156 = cp[16]*cp[22];
	const Interval tmp_157 = (R(1, 2))*tmp_153 + (R(1, 2))*tmp_154 - R(1, 2)*tmp_155 - R(1, 2)*tmp_156;
	const Interval tmp_158 = cp[13]*cp[23];
	const Interval tmp_159 = cp[19]*cp[21];
	const Interval tmp_160 = cp[14]*cp[21];
	const Interval tmp_161 = cp[15]*cp[20];
	const Interval tmp_162 = cp[15]*cp[21];
	const Interval tmp_163 = cp[16]*cp[23];
	const Interval tmp_164 = cp[17]*cp[22];
	const Interval tmp_165 = cp[17]*cp[23];
	const Interval tmp_166 = cp[12]*cp[23] + cp[13]*cp[22] + cp[18]*cp[21] + cp[19]*cp[20] + tmp_153 + tmp_154 - tmp_155 - tmp_156 + tmp_158 + tmp_159 - tmp_160 - tmp_161 - tmp_162 - tmp_163 - tmp_164 - tmp_165;
	const Interval tmp_167 = (R(1, 2))*tmp_162;
	const Interval tmp_168 = (R(1, 2))*tmp_165;
	const Interval tmp_169 = -tmp_47 + tmp_49;
	const Interval tmp_170 = tmp_147 - R(1, 2)*tmp_62 + (R(1, 2))*tmp_65;
	const Interval tmp_171 = cp[12]*cp[18] - 3*cp[12]*cp[22] - cp[14]*cp[16] - 3*cp[18]*cp[20] - 3*cp[20]*cp[2] - 3*cp[22]*cp[8] + tmp_1 + tmp_134 + 3*tmp_155 + 3*tmp_156 + tmp_83 + tmp_93;
	const Interval tmp_172 = cp[13]*cp[19];
	const Interval tmp_173 = 3*tmp_162;
	const Interval tmp_174 = 3*tmp_165;
	out[0] = tmp_19;
	out[1] = (R(1, 4))*cp[12]*cp[7] + (R(1, 4))*cp[13]*cp[6] + (R(1, 4))*cp[13]*cp[7] - R(1, 4)*cp[14]*cp[5] - R(1, 4)*cp[15]*cp[4] + (R(1, 4))*tmp_19 + (R(1, 4))*tmp_20 + (R(1, 4))*tmp_21 + (R(1, 4))*tmp_22 - R(1, 4)*tmp_23 - R(1, 4)*tmp_24 - R(1, 4)*tmp_25 - R(1, 4)*tmp_26 - R(1, 4)*tmp_28 - R(3, 4)*tmp_29 - R(3, 4)*tmp_30 - R(3, 4)*tmp_31 - R(1, 4)*tmp_33 - R(1, 4)*tmp_35 - R(3, 4)*tmp_36 - R(1, 4)*tmp_38 + (R(3, 4))*tmp_39 + (R(3, 4))*tmp_40 + (R(3, 4))*tmp_41 + (R(3, 4))*tmp_42 + (R(3, 4))*tmp_43 + (R(3, 4))*tmp_44 + (R(3, 4))*tmp_45 + (R(3, 4))*tmp_46 + (R(1, 4))*tmp_50 + (R(1, 4))*tmp_54 + (R(1, 4))*tmp_59 + (R(1, 4))*tmp_69;
	out[2] = 3*cp[11]*cp[1] + cp[13]*cp[7] + 3*cp[15]*cp[1] + 3*cp[3]*cp[5] - tmp_23 - tmp_33 - tmp_38 - tmp_48 - tmp_70 - tmp_71 - tmp_72 - tmp_73;
	out[3] = (R(1, 2))*cp[0]*cp[14] + (R(1, 2))*cp[10]*cp[16] + (R(1, 2))*cp[2]*cp[4] + (R(1, 2))*cp[6]*cp[8] - tmp_74 - tmp_76 - tmp_77 - tmp_79 - tmp_84 - tmp_85;
	out[4] = (R(1, 8))*cp[0]*cp[14] + (R(1, 8))*cp[0]*cp[15] + (R(1, 8))*cp[10]*cp[16] + (R(1, 8))*cp[10]*cp[17] + (R(1, 8))*cp[11]*cp[16] + (R(1, 8))*cp[11]*cp[17] + (R(1, 8))*cp[14]*cp[1] + (R(1, 8))*cp[15]*cp[1] + (R(1, 8))*cp[2]*cp[4] + (R(1, 8))*cp[2]*cp[5] + (R(1, 8))*cp[3]*cp[4] + (R(1, 8))*cp[3]*cp[5] + (R(1, 8))*cp[6]*cp[8] + (R(1, 8))*cp[6]*cp[9] + (R(1, 8))*cp[7]*cp[8] + (R(1, 8))*cp[7]*cp[9] - R(1, 8)*tmp_115 - R(1, 8)*tmp_124 - R(1, 8)*tmp_2 - R(1, 8)*tmp_29 - R(1, 8)*tmp_3 - R(1, 8)*tmp_30 - R(1, 8)*tmp_31 - R(1, 8)*tmp_32 - R(1, 8)*tmp_36 - R(1, 8)*tmp_37 - R(1, 8)*tmp_75 - R(1, 8)*tmp_78 - R(1, 8)*tmp_86 - R(1, 8)*tmp_87 - R(1, 8)*tmp_88 - R(1, 8)*tmp_89 - R(1, 8)*tmp_90 - R(1, 8)*tmp_91 - R(1, 8)*tmp_93 - R(1, 8)*tmp_98;
	out[5] = (R(1, 2))*tmp_125 + tmp_126 + (R(1, 2))*tmp_127 - tmp_128 - tmp_129 - tmp_130 + tmp_132 + tmp_133 - tmp_51 + tmp_53;
	out[6] = tmp_138;
	out[7] = -R(1, 4)*cp[16]*cp[7] - R(1, 4)*cp[17]*cp[6] + (R(1, 4))*cp[18]*cp[5] + (R(1, 4))*cp[19]*cp[4] + (R(1, 4))*tmp_138 - R(1, 4)*tmp_139 - R(1, 4)*tmp_140 + (R(3, 4))*tmp_141 + (R(3, 4))*tmp_142 + (R(3, 4))*tmp_143 + (R(3, 4))*tmp_144 + (R(1, 4))*tmp_149 + (R(1, 4))*tmp_150 + (R(3, 4))*tmp_24 + (R(3, 4))*tmp_25 + (R(1, 4))*tmp_27 + (R(1, 4))*tmp_34 - R(1, 4)*tmp_42 - R(1, 4)*tmp_43 + (R(1, 4))*tmp_59 - R(3, 4)*tmp_86 - R(3, 4)*tmp_87 - R(3, 4)*tmp_89 - R(3, 4)*tmp_90;
	out[8] = tmp_149;
	out[9] = -tmp_10 + tmp_11 + tmp_151 - tmp_152 + tmp_157 + (R(1, 2))*tmp_4 + (R(1, 2))*tmp_6 - tmp_74 - tmp_77 + tmp_84;
	out[10] = (R(1, 8))*tmp_115 + (R(1, 8))*tmp_116 + (R(1, 8))*tmp_117 + (R(1, 8))*tmp_118 + (R(1, 8))*tmp_119 - R(1, 8)*tmp_120 - R(1, 8)*tmp_121 - R(1, 8)*tmp_122 - R(1, 8)*tmp_123 + (R(1, 8))*tmp_137 + (R(1, 8))*tmp_148 + (R(1, 8))*tmp_150 + (R(1, 8))*tmp_166 - R(1, 8)*tmp_2 - R(1, 8)*tmp_29 - R(1, 8)*tmp_3 - R(1, 8)*tmp_30 - R(1, 8)*tmp_31 - R(1, 8)*tmp_32 - R(1, 8)*tmp_36 - R(1, 8)*tmp_37 + (R(1, 8))*tmp_39 + (R(1, 8))*tmp_4 + (R(1, 8))*tmp_40 + (R(1, 8))*tmp_41 + (R(1, 8))*tmp_44 + (R(1, 8))*tmp_45 + (R(1, 8))*tmp_46 + (R(1, 8))*tmp_6;
	out[11] = -tmp_126 + tmp_129 + tmp_133 + (R(1, 2))*tmp_158 + (R(1, 2))*tmp_159 - tmp_167 - tmp_168 + tmp_169 + tmp_170 + tmp_72;
	out[12] = -tmp_0 + (R(1, 2))*tmp_135 + (R(1, 2))*tmp_136 - tmp_151 + tmp_152 + tmp_157 + tmp_5 - tmp_76 - tmp_79 + tmp_83 + tmp_85;
	out[13] = (R(1, 8))*tmp_100 + (R(1, 8))*tmp_101 + (R(1, 8))*tmp_102 + (R(1, 8))*tmp_103 + (R(1, 8))*tmp_104 + (R(1, 8))*tmp_105 + (R(1, 8))*tmp_106 + (R(1, 8))*tmp_107 + (R(1, 8))*tmp_108 - R(1, 8)*tmp_109 - R(1, 8)*tmp_110 + (R(1, 8))*tmp_111 + (R(1, 8))*tmp_112 - R(1, 8)*tmp_113 - R(1, 8)*tmp_114 + (R(1, 8))*tmp_124 + (R(1, 8))*tmp_125 + (R(1, 8))*tmp_127 + (R(1, 8))*tmp_135 + (R(1, 8))*tmp_136 + (R(1, 8))*tmp_141 + (R(1, 8))*tmp_142 + (R(1, 8))*tmp_143 + (R(1, 8))*tmp_144 + (R(1, 8))*tmp_166 + (R(1, 8))*tmp_18 + (R(1, 8))*tmp_69 - R(1, 8)*tmp_75 - R(1, 8)*tmp_78 - R(1, 8)*tmp_86 - R(1, 8)*tmp_87 - R(1, 8)*tmp_88 - R(1, 8)*tmp_89 - R(1, 8)*tmp_90 - R(1, 8)*tmp_91 + (R(1, 8))*tmp_99;
	out[14] = (R(1, 2))*cp[11]*cp[17] + (R(1, 2))*cp[13]*cp[23] + (R(1, 2))*cp[19]*cp[21] + (R(1, 2))*cp[7]*cp[9] - tmp_128 - tmp_130 - tmp_132 - tmp_167 - tmp_168 - tmp_170;
	out[15] = -tmp_171;
	out[16] = -R(1, 4)*cp[12]*cp[19] + (R(3, 4))*cp[12]*cp[23] - R(1, 4)*cp[13]*cp[18] + (R(3, 4))*cp[13]*cp[22] + (R(3, 4))*cp[13]*cp[23] + (R(1, 4))*cp[14]*cp[17] + (R(1, 4))*cp[15]*cp[16] + (R(1, 4))*cp[15]*cp[17] + (R(3, 4))*cp[18]*cp[21] + (R(3, 4))*cp[19]*cp[20] + (R(3, 4))*cp[19]*cp[21] + (R(3, 4))*cp[20]*cp[3] + (R(3, 4))*cp[21]*cp[2] + (R(3, 4))*cp[22]*cp[9] + (R(3, 4))*cp[23]*cp[8] + (R(3, 4))*cp[23]*cp[9] + (R(1, 4))*cp[2]*cp[9] + (R(1, 4))*cp[3]*cp[8] - R(1, 4)*tmp_139 - R(1, 4)*tmp_140 - R(1, 4)*tmp_145 - R(3, 4)*tmp_160 - R(3, 4)*tmp_161 - R(3, 4)*tmp_163 - R(3, 4)*tmp_164 - R(1, 4)*tmp_169 - R(1, 4)*tmp_171 - R(1, 4)*tmp_172 - R(1, 4)*tmp_173 - R(1, 4)*tmp_174 - R(1, 4)*tmp_28 - R(1, 4)*tmp_35 - R(1, 4)*tmp_57 - R(1, 4)*tmp_58 - R(1, 4)*tmp_73 - R(1, 4)*tmp_98;
	out[17] = cp[15]*cp[17] + tmp_131 + tmp_146 + 3*tmp_158 + 3*tmp_159 - tmp_172 - tmp_173 - tmp_174 + tmp_54 + tmp_71 - tmp_94 + tmp_95;
}

}
