#include "lagrangeVector.hpp"

#define R(p, q) (Interval(p) / q)

namespace element_validity {
template<>
void lagrangeVectorT<2, 2, 1>(const std::vector<fp_t> &cpFP, std::vector<Interval> &out) {
	out.resize(3);
	const uint S = cpFP.size();
	std::vector<Interval> cp(S);
	for (uint i = 0; i < S; ++i) cp[i] = cpFP[i];
	const Interval tmp_0 = cp[0]*cp[10] - cp[0]*cp[6] - cp[10]*cp[4] + cp[2]*cp[4] - cp[2]*cp[8] + cp[6]*cp[8];
	const Interval tmp_1 = cp[11]*cp[1] - cp[11]*cp[5] - cp[1]*cp[7] + cp[3]*cp[5] - cp[3]*cp[9] + cp[7]*cp[9];
	out[0] = -tmp_0;
	out[1] = -R(1, 4)*cp[0]*cp[11] + (R(1, 4))*cp[0]*cp[7] - R(1, 4)*cp[10]*cp[1] + (R(1, 4))*cp[10]*cp[5] + (R(1, 4))*cp[11]*cp[4] + (R(1, 4))*cp[1]*cp[6] - R(1, 4)*cp[2]*cp[5] + (R(1, 4))*cp[2]*cp[9] - R(1, 4)*cp[3]*cp[4] + (R(1, 4))*cp[3]*cp[8] - R(1, 4)*cp[6]*cp[9] - R(1, 4)*cp[7]*cp[8] - R(1, 4)*tmp_0 - R(1, 4)*tmp_1;
	out[2] = -tmp_1;
}}
