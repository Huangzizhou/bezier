#include "lagrangeVector.hpp"

#define R(p, q) (Interval(p) / q)
#define I const Interval 

namespace element_validity {
template<>
void lagrangeVector<3, 3, 1>(const std::span<const fp_t> cpFP, const std::span<Interval> out) {
	assert(out.size() == 1);
	const uint S = cpFP.size();
	std::vector<Interval> cp(S);
	for (uint i = 0; i < S; ++i) cp[i] = cpFP[i];
	out[0] = cp[0]*(-cp[10]*cp[5] + cp[10]*cp[8] + cp[11]*cp[4] - cp[11]*cp[7] - cp[4]*cp[8] + cp[5]*cp[7]) + cp[10]*(cp[2]*cp[3] - cp[2]*cp[6] - cp[3]*cp[8] + cp[5]*cp[6]) + cp[11]*(-cp[1]*cp[3] + cp[1]*cp[6] + cp[3]*cp[7] - cp[4]*cp[6]) + cp[1]*(cp[3]*cp[8] - cp[5]*cp[6] + cp[5]*cp[9] - cp[8]*cp[9]) + cp[2]*(-cp[3]*cp[7] + cp[4]*cp[6] - cp[4]*cp[9] + cp[7]*cp[9]) + cp[4]*cp[8]*cp[9] - cp[5]*cp[7]*cp[9];
}}
#undef R
#undef I