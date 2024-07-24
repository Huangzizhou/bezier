#include "lagrangeVector.hpp"

#define R(p, q) (Interval(p) / q)
#define I const Interval 

namespace element_validity {
template<>
void lagrangeVector<2, 2, 2>(const std::span<const fp_t> cpFP, const std::span<Interval> out) {
	assert(out.size() == 6);
	const uint S = cpFP.size();
	std::vector<Interval> cp(S);
	for (uint i = 0; i < S; ++i) cp[i] = cpFP[i];
	out[0] = cp[0]*(-12*cp[11] - 3*cp[3] + 3*cp[5] + 12*cp[7]) + cp[10]*(12*cp[1] + 4*cp[3] - 16*cp[7]) + cp[11]*(-4*cp[2] + 16*cp[6]) + cp[1]*(3*cp[2] - 3*cp[4] - 12*cp[6]) + cp[2]*cp[5] - cp[3]*cp[4] + 4*cp[4]*cp[7] - 4*cp[5]*cp[6];
	out[1] = cp[0]*(-2*cp[11] - cp[3] - cp[5] + 2*cp[7] + 2*cp[9]) + cp[10]*(2*cp[1] - 2*cp[5]) + 2*cp[11]*cp[4] + cp[1]*(cp[2] + cp[4] - 2*cp[6] - 2*cp[8]) - cp[2]*cp[5] + cp[3]*cp[4] + cp[4]*(-2*cp[7] - 2*cp[9]) + cp[5]*(2*cp[6] + 2*cp[8]);
	out[2] = cp[0]*(cp[3] + 3*cp[5] - 4*cp[9]) + cp[10]*(-4*cp[3] - 12*cp[5] + 16*cp[9]) + cp[11]*(4*cp[2] + 12*cp[4] - 16*cp[8]) + cp[1]*(-cp[2] - 3*cp[4] + 4*cp[8]) - 3*cp[2]*cp[5] + 3*cp[3]*cp[4] - 12*cp[4]*cp[9] + 12*cp[5]*cp[8];
	out[3] = cp[0]*(-2*cp[11] + cp[3] + cp[5] + 2*cp[7] - 2*cp[9]) + cp[10]*(2*cp[1] - 2*cp[3]) + 2*cp[11]*cp[2] + cp[1]*(-cp[2] - cp[4] - 2*cp[6] + 2*cp[8]) + cp[2]*(-cp[5] - 2*cp[7] + 2*cp[9]) + cp[3]*(cp[4] + 2*cp[6] - 2*cp[8]);
	out[4] = cp[0]*(-cp[3] + cp[5]) + cp[10]*(2*cp[3] - 2*cp[5]) + cp[11]*(-2*cp[2] + 2*cp[4]) + cp[1]*(cp[2] - cp[4]) + cp[2]*(cp[5] - 2*cp[7] + 2*cp[9]) + cp[3]*(-cp[4] + 2*cp[6] - 2*cp[8]) + cp[4]*(2*cp[7] - 2*cp[9]) + cp[5]*(-2*cp[6] + 2*cp[8]);
	out[5] = cp[0]*(-3*cp[3] - cp[5] + 4*cp[9]) + cp[1]*(3*cp[2] + cp[4] - 4*cp[8]) + cp[2]*(-3*cp[5] - 12*cp[7] + 12*cp[9]) + cp[3]*(3*cp[4] + 12*cp[6] - 12*cp[8]) - 4*cp[4]*cp[7] + 4*cp[5]*cp[6] - 16*cp[6]*cp[9] + 16*cp[7]*cp[8];
}}
#undef R
#undef I