#include "lagrangeVector.hpp"

#define R(p, q) (Interval(p) / q)
#define I const Interval 

namespace element_validity {
template<>
void lagrangeVector<1, 1, 4>(const std::span<const fp_t> cpFP, const std::span<Interval> out) {
	assert(out.size() == 4);
	const uint S = cpFP.size();
	std::vector<Interval> cp(S);
	for (uint i = 0; i < S; ++i) cp[i] = cpFP[i];
	out[0] = -R(25, 3)*cp[0] + 16*cp[1] - 12*cp[2] + (R(16, 3))*cp[3] - cp[4];
	out[1] = -R(7, 81)*cp[0] - R(368, 81)*cp[1] + (R(148, 27))*cp[2] - R(80, 81)*cp[3] + (R(11, 81))*cp[4];
	out[2] = -R(11, 81)*cp[0] + (R(80, 81))*cp[1] - R(148, 27)*cp[2] + (R(368, 81))*cp[3] + (R(7, 81))*cp[4];
	out[3] = cp[0] - R(16, 3)*cp[1] + 12*cp[2] - 16*cp[3] + (R(25, 3))*cp[4];
}}
#undef R
#undef I