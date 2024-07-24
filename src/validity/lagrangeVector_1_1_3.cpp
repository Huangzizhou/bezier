#include "lagrangeVector.hpp"

#define R(p, q) (Interval(p) / q)
#define I const Interval 

namespace element_validity {
template<>
void lagrangeVector<1, 1, 3>(const std::span<const fp_t> cpFP, const std::span<Interval> out) {
	assert(out.size() == 3);
	const uint S = cpFP.size();
	std::vector<Interval> cp(S);
	for (uint i = 0; i < S; ++i) cp[i] = cpFP[i];
	out[0] = -R(11, 2)*cp[0] + 9*cp[1] - R(9, 2)*cp[2] + cp[3];
	out[1] = (R(1, 8))*cp[0] - R(27, 8)*cp[1] + (R(27, 8))*cp[2] - R(1, 8)*cp[3];
	out[2] = -cp[0] + (R(9, 2))*cp[1] - 9*cp[2] + (R(11, 2))*cp[3];
}}
#undef R
#undef I