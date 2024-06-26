#include "lagrangeVector.hpp"

#define R(p, q) (Interval(p) / q)

namespace element_validity {
#ifdef LAGVEC_GCC_O0
#pragma GCC push_options
#pragma GCC optimize ("-O0")
#endif
template<>
void lagrangeVectorT<1, 1, 1>(const std::vector<fp_t> &cpFP, std::vector<Interval> &out) {
	out.resize(2);
	const uint S = cpFP.size();
	std::vector<Interval> cp(S);
	for (uint i = 0; i < S; ++i) cp[i] = cpFP[i];
	out[0] = -cp[0] + cp[2];
	out[1] = -cp[1] + cp[3];
}}
#ifdef LAGVEC_GCC_O0
#pragma GCC pop_options
#endif
#undef R