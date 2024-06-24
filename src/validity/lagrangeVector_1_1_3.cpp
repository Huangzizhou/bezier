#include "lagrangeVector.hpp"

#define R(p, q) (Interval(p) / q)

namespace element_validity {
template<>
void lagrangeVectorT<1, 1, 3>(const std::vector<fp_t> &cpFP, std::vector<Interval> &out) {
	out.resize(6);
	const uint S = cpFP.size();
	std::vector<Interval> cp(S);
	for (uint i = 0; i < S; ++i) cp[i] = cpFP[i];
	out[0] = -R(11, 2)*cp[0] + 9*cp[2] - R(9, 2)*cp[4] + cp[6];
	out[1] = -R(11, 2)*cp[1] + 9*cp[3] - R(9, 2)*cp[5] + cp[7];
	out[2] = (R(1, 8))*cp[0] - R(27, 8)*cp[2] + (R(27, 8))*cp[4] - R(1, 8)*cp[6];
	out[3] = (R(1, 8))*cp[1] - R(27, 8)*cp[3] + (R(27, 8))*cp[5] - R(1, 8)*cp[7];
	out[4] = -cp[0] + (R(9, 2))*cp[2] - 9*cp[4] + (R(11, 2))*cp[6];
	out[5] = -cp[1] + (R(9, 2))*cp[3] - 9*cp[5] + (R(11, 2))*cp[7];
}}
#undef R