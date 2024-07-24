#include "lagrangeVector.hpp"

#define R(p, q) (Interval(p) / q)
#define I const Interval 

namespace element_validity {
template<>
void lagrangeVectorT<2, 2, 1>(const std::span<const fp_t> cpFP, const std::span<Interval> out) {
	assert(out.size() == 3);
	const uint S = cpFP.size();
	std::vector<Interval> cp(S);
	for (uint i = 0; i < S; ++i) cp[i] = cpFP[i];
	I tmp_0 = cp[0]*cp[10];
	I tmp_1 = cp[2]*cp[4];
	I tmp_2 = cp[6]*cp[8];
	I tmp_3 = cp[11]*cp[1];
	I tmp_4 = cp[3]*cp[5];
	I tmp_5 = cp[7]*cp[9];
	out[0] = cp[0]*cp[6] + cp[10]*cp[4] + cp[2]*cp[8] - tmp_0 - tmp_1 - tmp_2;
	out[1] = -R(1, 4)*cp[0]*cp[11] + (R(1, 4))*cp[0]*cp[6] + (R(1, 4))*cp[0]*cp[7] - R(1, 4)*cp[10]*cp[1] + (R(1, 4))*cp[10]*cp[4] + (R(1, 4))*cp[10]*cp[5] + (R(1, 4))*cp[11]*cp[4] + (R(1, 4))*cp[11]*cp[5] + (R(1, 4))*cp[1]*cp[6] + (R(1, 4))*cp[1]*cp[7] - R(1, 4)*cp[2]*cp[5] + (R(1, 4))*cp[2]*cp[8] + (R(1, 4))*cp[2]*cp[9] - R(1, 4)*cp[3]*cp[4] + (R(1, 4))*cp[3]*cp[8] + (R(1, 4))*cp[3]*cp[9] - R(1, 4)*cp[6]*cp[9] - R(1, 4)*cp[7]*cp[8] - R(1, 4)*tmp_0 - R(1, 4)*tmp_1 - R(1, 4)*tmp_2 - R(1, 4)*tmp_3 - R(1, 4)*tmp_4 - R(1, 4)*tmp_5;
	out[2] = cp[11]*cp[5] + cp[1]*cp[7] + cp[3]*cp[9] - tmp_3 - tmp_4 - tmp_5;
}}
#undef R
#undef I