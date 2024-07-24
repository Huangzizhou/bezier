#include "lagrangeVector.hpp"

#define R(p, q) (Interval(p) / q)
#define I const Interval 

namespace element_validity {
template<>
void lagrangeVector<2, 1, 1>(const std::span<const fp_t> cpFP, const std::span<Interval> out) {
	assert(out.size() == 4);
	const uint S = cpFP.size();
	std::vector<Interval> cp(S);
	for (uint i = 0; i < S; ++i) cp[i] = cpFP[i];
	out[0] = cp[0]*(cp[3] - cp[7]) + cp[1]*(-cp[2] + cp[6]) + cp[2]*cp[7] - cp[3]*cp[6];
	out[1] = cp[0]*(cp[5] - cp[7]) + cp[1]*(-cp[4] + cp[6]) + cp[4]*cp[7] - cp[5]*cp[6];
	out[2] = cp[0]*(cp[3] - cp[5]) + cp[1]*(-cp[2] + cp[4]) + cp[2]*cp[5] - cp[3]*cp[4];
	out[3] = cp[2]*(cp[5] - cp[7]) + cp[3]*(-cp[4] + cp[6]) + cp[4]*cp[7] - cp[5]*cp[6];
}}
#undef R
#undef I