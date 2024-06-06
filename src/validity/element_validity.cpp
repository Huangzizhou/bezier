#include "element_validity.hpp"

// #define R(p, q) static_cast<Interval>(Rational(p, q))
#define R(p, q) Interval(p) / q

namespace element_validity {
template<>
void lagrangeVector<1, 1, 1, true>(const std::vector<fp_t> &cpFP, std::vector<Interval> &out) {
	out.resize(2);
	std::vector<Interval> cp(cpFP.size());
	for (uint i = 0; i < cpFP.size(); ++i) cp[i] = cpFP[i];
	out[0] = -cp[0] + cp[2];
	out[1] = -cp[1] + cp[3];
}
template<>
void lagrangeVector<1, 1, 2, true>(const std::vector<fp_t> &cpFP, std::vector<Interval> &out) {
	out.resize(4);
	std::vector<Interval> cp(cpFP.size());
	for (uint i = 0; i < cpFP.size(); ++i) cp[i] = cpFP[i];
	const Interval tmp_0 = -cp[2];
	const Interval tmp_1 = -cp[3];
	out[0] = -3*cp[0] - cp[4] - tmp_0;
	out[1] = -3*cp[1] - cp[5] - tmp_1;
	out[2] = cp[0] + 3*cp[4] + tmp_0;
	out[3] = cp[1] + 3*cp[5] + tmp_1;
}
template<>
void lagrangeVector<1, 1, 3, true>(const std::vector<fp_t> &cpFP, std::vector<Interval> &out) {
	out.resize(6);
	std::vector<Interval> cp(cpFP.size());
	for (uint i = 0; i < cpFP.size(); ++i) cp[i] = cpFP[i];
	const Interval tmp_0 = -cp[6];
	const Interval tmp_1 = -cp[2] + cp[4];
	const Interval tmp_2 = -cp[7];
	const Interval tmp_3 = -cp[3] + cp[5];
	out[0] = -R(11, 2)*cp[0] - tmp_0 - tmp_1;
	out[1] = -R(11, 2)*cp[1] - tmp_2 - tmp_3;
	out[2] = (R(1, 8))*cp[0] - R(1, 2)*cp[2] + (R(1, 2))*cp[4] + (R(1, 8))*tmp_0;
	out[3] = (R(1, 8))*cp[1] - R(1, 2)*cp[3] + (R(1, 2))*cp[5] + (R(1, 8))*tmp_2;
	out[4] = -cp[0] + (R(11, 2))*cp[6] - tmp_1;
	out[5] = -cp[1] + (R(11, 2))*cp[7] - tmp_3;
}
template<>
void lagrangeVector<2, 2, 1, true>(const std::vector<fp_t> &cpFP, std::vector<Interval> &out) {
	out.resize(3);
	std::vector<Interval> cp(cpFP.size());
	for (uint i = 0; i < cpFP.size(); ++i) cp[i] = cpFP[i];
	const Interval tmp_0 = cp[0]*cp[10] - cp[0]*cp[6] - cp[10]*cp[4] + cp[2]*cp[4] - cp[2]*cp[8] + cp[6]*cp[8];
	const Interval tmp_1 = cp[11]*cp[1] - cp[11]*cp[5] - cp[1]*cp[7] + cp[3]*cp[5] - cp[3]*cp[9] + cp[7]*cp[9];
	out[0] = tmp_0;
	out[1] = (R(1, 4))*cp[0]*cp[11] - R(1, 4)*cp[0]*cp[7] + (R(1, 4))*cp[10]*cp[1] - R(1, 4)*cp[10]*cp[5] - R(1, 4)*cp[11]*cp[4] - R(1, 4)*cp[1]*cp[6] + (R(1, 4))*cp[2]*cp[5] - R(1, 4)*cp[2]*cp[9] + (R(1, 4))*cp[3]*cp[4] - R(1, 4)*cp[3]*cp[8] + (R(1, 4))*cp[6]*cp[9] + (R(1, 4))*cp[7]*cp[8] + (R(1, 4))*tmp_0 + (R(1, 4))*tmp_1;
	out[2] = tmp_1;
}
}
