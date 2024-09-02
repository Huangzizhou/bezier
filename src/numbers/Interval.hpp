#pragma once

#include <sstream>
#include "Rational.hpp"
#include "utils/globals.hpp"

#ifdef IPRED_ARITHMETIC

#include "numerics.h"
namespace element_validity {

class Interval {
	public:
	static void init();

	private:
	interval_number data;

	inline Interval(const interval_number &o) :
		data(-o.inf(), o.sup()) {}

	Interval fromRational(const Rational &rat) {
		// TODO change this function to not use nextafter
		const double inf = std::numeric_limits<double>::max();
		const double d = static_cast<double>(rat);

		if (rat < 0) return Interval(std::nextafter(d, -inf), d);
		if (rat > 0) return Interval(d, std::nextafter(d, inf));
		return Interval(0);
	}

	public:
	// Constructors
	inline Interval(const double &value_lo, const double &value_hi) :
		data(-value_lo, value_hi) {}
	inline Interval(const double &val) : data(val) {}
	inline Interval() : Interval(0.) {}
	inline Interval(const Rational &rat) : Interval(fromRational(rat)) {}

	// Access
	inline double lower() const { return data.inf(); }
	inline double upper() const { return data.sup(); }

	// Operators
	inline Interval operator-() const { return -data; }
	inline Interval operator+(const Interval &o) const { return data + o.data; }
	inline Interval operator-(const Interval &o) const { return data - o.data; }
	inline Interval operator*(const Interval &o) const { return data * o.data; }

	inline Interval operator+(double o) const { return data + o; }
	inline friend Interval operator+(double o, const Interval m)
		{ return m + o; }
	inline Interval operator-(double o) const { return data - o; }
	inline friend Interval operator-(double o, const Interval m)
		{ return m - o; }
	inline Interval operator*(double o) const { return data * o; }
	inline friend Interval operator*(double o, const Interval m)
		{ return m * o; }
	inline Interval operator/(double o) const { return data / o; }
	inline friend Interval operator/(double o, const Interval m)
		{ return m / o; }

	inline void operator+=(const Interval &o) { data += o.data; }
	inline void operator-=(const Interval &o) { data -= o.data; }
	inline void operator*=(const Interval &o) { data *= o.data; }

	inline void operator+=(double o) { data += o; }
	inline void operator-=(double o) { data -= o; }
	inline void operator*=(double o) { data *= o; }
	inline void operator/=(double o) { data /= o; }

	inline bool operator<(const Interval &x) const 
		{ return data < x.data; }
	inline bool operator<=(const Interval &x) const 
		{ return data <= x.data; }
	inline bool operator>(const Interval &x) const 
		{ return data > x.data; }
	inline bool operator>=(const Interval &x) const 
		{ return data >= x.data; }
	inline bool operator==(const Interval &x) const 
		{ return data == x.data; }
	inline bool operator!=(const Interval &x) const 
		{ return data != x.data; }

	inline bool operator<(double x) const 
		{ return data < x; }
	inline bool operator<=(double x) const 
		{ return data <= x; }
	inline bool operator>(double x) const 
		{ return data > x; }
	inline bool operator>=(double x) const 
		{ return data >= x; }
	inline bool operator==(double x) const 
		{ return data == x; }
	inline bool operator!=(double x) const 
		{ return data != x; }

	inline Interval abs() const { return data.abs(); }
	inline Interval pow(unsigned int e) const { return data.pow(e); }
	// inline int signIsReliable() const { return data.signIsReliable(); }
	inline double width() const { return data.width(); }
	inline double middle() const { return data.getMid(); }
	inline bool contains(double x) const { return (data-x).containsZero(); }
	inline bool contains(const Interval &o) const
		{ return contains(o.lower()) && contains(o.upper()); }

	friend Interval min(const Interval &a, const Interval &b) {
		return min(a.data, b.data);
	}
	friend Interval max(const Interval &a, const Interval &b) {
		return max(a.data, b.data);
	}


	// Conversions
	friend std::ostream& operator<<(std::ostream &ost, const Interval &r) {
		const auto p = std::setprecision(
			std::numeric_limits<double>::max_digits10);
		ost << '[' << p << r.lower() << ',' << p << r.upper() << ']';
        return ost;
    }

	explicit operator std::string() const {
		std::stringstream s;
		s << *this;
		return s.str();
    }

	explicit operator double() const { return middle(); }
};
}
#else

#warning Using the naive interval implementation: computations will be slower!

#include <algorithm>	// for minmax
#include <cmath>
#include <stdexcept>

namespace element_validity {

class Interval {
private:
	fp_t lo, hi;

public:
	static void init() {};
	// Constructors
	Interval(const fp_t &value_lo, const fp_t &value_hi) :
		lo(value_lo), hi(value_hi) {
		if(lo>hi) throw std::invalid_argument(
			"Lower bound cannot be greater than upper bound");
	}
	Interval(const fp_t &value) : lo(value), hi(value)  {}
	Interval() : lo(0.), hi(0.) {}
	Interval(const Interval &o) : lo(o.lo), hi(o.hi) {}
	inline Interval(const Rational &rat) : Interval(fromRational(rat)) {}

	Interval fromRational(const Rational &rat);

	// Access
	inline const fp_t& lower() const { return lo; }
	inline const fp_t& upper() const { return hi; }

	// Opposite
	inline Interval operator-() const { return Interval(-hi, -lo); }

	// Absolute value
	Interval abs() {
		if (hi < 0) return Interval(-hi, -lo);
		if (lo < 0) return Interval(0, std::max(hi, -lo));
		return Interval(lo, hi);
	}

	// Addition
	Interval operator+(const Interval &o) const {
		return Interval(roundDn(lo + o.lo), roundUp(hi + o.hi));
	}
	Interval& operator+=(const Interval& o) {
		lo = roundDn(lo + o.lo);
		hi = roundUp(hi + o.hi);
		return *this;
	}
	Interval operator+(const fp_t &o) const {
		return Interval(roundDn(lo + o), roundUp(hi + o));
	}
	Interval& operator+=(fp_t o) {
		lo = roundDn(lo + o);
		hi = roundUp(hi + o);
		return *this;
	}
	
	// Subtraction
	Interval operator-(const Interval &o) const {
		return Interval(roundDn(lo - o.hi), roundUp(hi - o.lo));
	}
	Interval& operator-=(const Interval& o) {
		lo = roundDn(lo - o.hi);
		hi = roundUp(hi - o.lo);
		return *this;
	}
	Interval operator-(const fp_t &o) const {
		return Interval(roundDn(lo - o), roundUp(hi - o));
	}
	Interval& operator-=(fp_t o) {
		lo = roundDn(lo - o);
		hi = roundUp(hi - o);
		return *this;
	}

	// Multiplication
	Interval operator*(const Interval &o) const {
		const std::pair<fp_t, fp_t> p = std::minmax({
			lo * o.lo, lo * o.hi,
			hi * o.lo, hi * o.hi
		});
		return Interval(roundDn(p.first), roundUp(p.second));
	}
	Interval& operator*=(const Interval& o) {
		const std::pair<fp_t, fp_t> p = std::minmax({
			lo * o.lo, lo * o.hi,
			hi * o.lo, hi * o.hi
		});
		lo = roundDn(p.first);
		hi = roundUp(p.second);
		return *this;
	}
	Interval operator*(const fp_t &o) const {
		const std::pair<fp_t, fp_t> p = std::minmax({lo * o, hi * o});
		return Interval(roundDn(p.first), roundUp(p.second));
	}
	Interval& operator*=(fp_t o) {
		const std::pair<fp_t, fp_t> p = std::minmax({lo * o, hi * o});
		lo = roundDn(p.first);
		hi = roundUp(p.second);
		return *this;
	}

	// Division by a number
	Interval operator/(fp_t o) const {
		if (o == 0) throw std::runtime_error("Division by 0.");
		const std::pair<fp_t, fp_t> p = std::minmax({lo / o, hi / o});
		return Interval(roundDn(p.first), roundUp(p.second));
	}
	Interval& operator/=(fp_t o) {
		if (o == 0) throw std::runtime_error("Division by 0.");
		const std::pair<fp_t, fp_t> p = std::minmax({lo / o, hi / o});
		lo = roundDn(p.first);
		hi = roundUp(p.second);
		return *this;
	}

	// Square
	Interval sqr() const {
		const std::pair<fp_t, fp_t> p = std::minmax(
			{lo * lo, lo * hi, hi * hi});
		return Interval(roundDn(p.first), roundUp(p.second));
	}

	// Generic power
	Interval pow(unsigned int e) const;

	// Whether the interval contains a number or interval
	inline bool contains(fp_t x) const { return (lo <= x && hi >= x); }
	inline bool contains(const Interval &o) const
		{ return contains(o.lower()) && contains(o.upper()); }

	// Comparisons with a floating point number
	inline bool operator<(fp_t x) const { return hi < x; }
	inline bool operator<=(fp_t x) const { return hi <= x; }
	inline bool operator>(fp_t x) const { return lo > x; }
	inline bool operator>=(fp_t x) const { return lo >= x; }
	inline bool operator==(fp_t x) const { return lo == x && hi == x; }
	inline bool operator!=(fp_t x) const { return !contains(x); }

	// Comparisons with another interval
	inline bool operator<(const Interval& o) const { return hi < o.lo; }
	inline bool operator<=(const Interval& o) const { return hi <= o.lo; }
	inline bool operator>(const Interval& o) const { return lo > o.hi; }
	inline bool operator>=(const Interval& o) const { return lo >= o.hi; }
	inline bool operator==(const Interval& o) const
		{ return lo == hi && lo == o.lo && hi == o.hi; }
	inline bool operator!=(const Interval& o) const
		{ return *this < o || *this > o; }

	// Minimum and maximum
	friend inline Interval min(const Interval &a, const Interval &b) {
		return Interval(
			(a.lo < b.lo) ? a.lo : b.lo,
			(a.hi < b.hi) ? a.hi : b.hi
		);
	}
	friend inline Interval max(const Interval &a, const Interval &b) {
		return Interval(
			(a.lo < b.lo) ? b.lo : a.lo,
			(a.hi < b.hi) ? b.hi : a.hi
		);
	}


	// Width or diameter
	inline fp_t width() const { return hi - lo; }

	// Midpoint (approximate)
	inline fp_t middle() const { return lo + (width() / 2); }

	// Conversions
	friend std::ostream& operator<<(std::ostream &ost, const Interval &r) {
		const auto p = std::setprecision(
			std::numeric_limits<double>::max_digits10);
		ost << '[' << p << r.lo << ',' << p << r.hi << ']';
        return ost;
    }

	explicit operator std::string() const {
		std::stringstream s;
		s << *this;
		return s.str();
    }

	explicit operator fp_t() const { return middle(); }

private:
	// Rounding up and down
	static const fp_t POSINF, NEGINF;
	inline static fp_t roundUp(const fp_t &value) {
		return std::nextafter(value, POSINF);
	}
	inline static fp_t roundDn(const fp_t &value) {
		return std::nextafter(value, NEGINF);
	}
};
}

#endif

/*
template<> struct Eigen::NumTraits<Interval> :
	Eigen::GenericNumTraits<Interval>
{
	typedef Interval Real;
	typedef Interval NonInteger;
	typedef Interval Nested;
	static inline Real epsilon() { return 0; }
	static inline Real dummy_precision() { return 0; }
	static inline Real digits10() { return 0; }

	enum {
		IsInteger = 0,
		IsSigned = 1,
		IsComplex = 0,
		RequireInitialization = 1,
		ReadCost = 2,
		AddCost = 4,
		MulCost = 30
	};
};
*/