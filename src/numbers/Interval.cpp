#include "Interval.hpp"

namespace element_validity {


#ifdef IPRED_ARITHMETIC

void Interval::init() {
	#ifdef IPRED_ARITHMETIC
	static bool initialized = false;
	if (initialized) return;
	initFPU();
	setFPUModeToRoundUP();
	initialized = true;
	#endif
};

#else

const fp_t Interval::POSINF = std::numeric_limits<fp_t>::max();
const fp_t Interval::NEGINF = std::copysign(POSINF, -1);

Interval Interval::pow(uint e) const {
	switch (e) {
	// Base cases
	case 0: return Interval(1);
	case 1: return *this;
	case 2: return sqr();
	default:
		// Exponent is even: divide by 2 and recurse
		if ((e % 2) == 0) return pow(e >> 1).sqr();
		// Exponent is not divisible by 2
		return *this * pow(e - 1);
	}
}

Interval Interval::fromRational(const Rational &rat) {
		const double d = static_cast<double>(rat);
		if (rat < 0) return Interval(roundDn(d), d);
		if (rat > 0) return Interval(d, roundUp(d));
		return Interval(0);
	}
#endif


}
