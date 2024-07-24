#include "lagrangeVector.hpp"

#define R(p, q) (Interval(p) / q)
#define I const Interval 

namespace element_validity {
template<>
void lagrangeVector<1, 1, 1>(const std::span<const fp_t> cpFP, const std::span<Interval> out) {
	assert(out.size() == 1);
	const uint S = cpFP.size();
	std::vector<Interval> cp(S);
	for (uint i = 0; i < S; ++i) cp[i] = cpFP[i];
	out[0] = -cp[0] + cp[1];
}}
#undef R
#undef I