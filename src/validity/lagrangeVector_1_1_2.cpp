#include "lagrangeVector.hpp"

#define R(p, q) (Interval(p) / q)

namespace element_validity {
template<>
void lagrangeVectorT<1, 1, 2>(const std::span<const fp_t> cpFP, const std::span<Interval> out) {
	assert(out.size() == 4);
	const uint S = cpFP.size();
	std::vector<Interval> cp(S);
	for (uint i = 0; i < S; ++i) cp[i] = cpFP[i];
	const Interval tmp_0 = -4*cp[2];
	const Interval tmp_1 = -4*cp[3];
	out[0] = -3*cp[0] - cp[4] - tmp_0;
	out[1] = -3*cp[1] - cp[5] - tmp_1;
	out[2] = cp[0] + 3*cp[4] + tmp_0;
	out[3] = cp[1] + 3*cp[5] + tmp_1;
}}
#undef R