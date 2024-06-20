#include "lagrangeVector.hpp"

#define R(p, q) (Interval(p) / q)

namespace element_validity {
template<>
void lagrangeVectorT<1, 1, 4>(const std::vector<fp_t> &cpFP, std::vector<Interval> &out) {
	out.resize(8);
	const uint S = cpFP.size();
	std::vector<Interval> cp(S);
	for (uint i = 0; i < S; ++i) cp[i] = cpFP[i];
	const Interval tmp_0 = 12*cp[4];
	const Interval tmp_1 = 12*cp[5];
	const Interval tmp_2 = -444*cp[4];
	const Interval tmp_3 = -444*cp[5];
	out[0] = -R(25, 3)*cp[0] + 16*cp[2] + (R(16, 3))*cp[6] - cp[8] - tmp_0;
	out[1] = -R(25, 3)*cp[1] + 16*cp[3] + (R(16, 3))*cp[7] - cp[9] - tmp_1;
	out[2] = -R(7, 81)*cp[0] - R(368, 81)*cp[2] - R(80, 81)*cp[6] + (R(11, 81))*cp[8] - R(1, 81)*tmp_2;
	out[3] = -R(7, 81)*cp[1] - R(368, 81)*cp[3] - R(80, 81)*cp[7] + (R(11, 81))*cp[9] - R(1, 81)*tmp_3;
	out[4] = -R(11, 81)*cp[0] + (R(80, 81))*cp[2] + (R(368, 81))*cp[6] + (R(7, 81))*cp[8] + (R(1, 81))*tmp_2;
	out[5] = -R(11, 81)*cp[1] + (R(80, 81))*cp[3] + (R(368, 81))*cp[7] + (R(7, 81))*cp[9] + (R(1, 81))*tmp_3;
	out[6] = cp[0] - R(16, 3)*cp[2] - 16*cp[6] + (R(25, 3))*cp[8] + tmp_0;
	out[7] = cp[1] - R(16, 3)*cp[3] - 16*cp[7] + (R(25, 3))*cp[9] + tmp_1;
}}
