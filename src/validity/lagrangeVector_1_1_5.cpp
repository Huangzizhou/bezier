#include "lagrangeVector.hpp"

#define R(p, q) (Interval(p) / q)

namespace element_validity {
template<>
void lagrangeVectorT<1, 1, 5>(const std::vector<fp_t> &cpFP, std::vector<Interval> &out) {
	out.resize(10);
	const uint S = cpFP.size();
	std::vector<Interval> cp(S);
	for (uint i = 0; i < S; ++i) cp[i] = cpFP[i];
	out[0] = -R(137, 12)*cp[0] + cp[10] + 25*cp[2] - 25*cp[4] + (R(50, 3))*cp[6] - R(25, 4)*cp[8];
	out[1] = cp[11] - R(137, 12)*cp[1] + 25*cp[3] - 25*cp[5] + (R(50, 3))*cp[7] - R(25, 4)*cp[9];
	out[2] = (R(1, 6144))*(-1269*cp[0] - 731*cp[10] - 37575*cp[2] + 51950*cp[4] - 17550*cp[6] + 5175*cp[8]);
	out[3] = (R(1, 6144))*(-731*cp[11] - 1269*cp[1] - 37575*cp[3] + 51950*cp[5] - 17550*cp[7] + 5175*cp[9]);
	out[4] = (R(1, 384))*(-9*cp[0] + 9*cp[10] + 125*cp[2] - 2250*cp[4] + 2250*cp[6] - 125*cp[8]);
	out[5] = (R(1, 384))*(9*cp[11] - 9*cp[1] + 125*cp[3] - 2250*cp[5] + 2250*cp[7] - 125*cp[9]);
	out[6] = (R(1, 6144))*(731*cp[0] + 1269*cp[10] - 5175*cp[2] + 17550*cp[4] - 51950*cp[6] + 37575*cp[8]);
	out[7] = (R(1, 6144))*(1269*cp[11] + 731*cp[1] - 5175*cp[3] + 17550*cp[5] - 51950*cp[7] + 37575*cp[9]);
	out[8] = -cp[0] + (R(137, 12))*cp[10] + (R(25, 4))*cp[2] - R(50, 3)*cp[4] + 25*cp[6] - 25*cp[8];
	out[9] = (R(137, 12))*cp[11] - cp[1] + (R(25, 4))*cp[3] - R(50, 3)*cp[5] + 25*cp[7] - 25*cp[9];
}}
