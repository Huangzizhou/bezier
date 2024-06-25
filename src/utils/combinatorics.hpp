#pragma once
#include "utils/globals.hpp"

namespace element_validity {

	// Binomial coefficient
	constexpr uint binom(uint n, uint k) {
		// edge cases
		if (k>n) return 0;
		if (k==0 || k==n) return 1;
		k = std::min(k, n-k);
		// accumulate results
		uint acc = 1;
		for (uint i=0; i<k; ++i) {
			acc = (acc * (n-i)) / (i+1);
		}
		return acc;
	}

	// Positive integer power
	template<typename bType>
	constexpr bType powi(bType n, uint e) {
		switch(e) {
		// Base cases
		case 0: return 1;
		case 1: return n;
		default:
			// Exponent is even: divide by 2 and recurse
			if ((e % 2) == 0) {
				const bType a = powi(n, e >> 1);
				return a*a;
			}
			// Exponent is not divisible by 2
			return n * powi(n, e - 1);
		}
	}

	// Number of control points
	constexpr uint nControlGeoMap(uint n, uint s, uint p) {
		return binom(p+s, s) * powi(p+1, n - s);
	}

	// Number of control points for the jacobian
	constexpr uint nControlJacobian(uint n, uint s, uint p, bool t) {
		uint L = binom(n*p, s) * powi(n * p, n - s);
		return t ? L * (n+1) : L;
	}
}