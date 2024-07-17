#include <catch2/catch_test_macros.hpp>
#include <element_validity.hpp>

using element_validity::Validity;
using element_validity::ValidityChecker;

TEST_CASE("Standard linear tet MTS") {
	ValidityChecker<3, 3, 1> checker;
	const double prec = .1;
    checker.setPrecisionTarget(prec);
    checker.setMaxSubdiv(1);
    const std::vector<double> cp = {
		0.,0.,	0.,0.,	0.,0.,
		1.,1.,	0.,0.,	0.,0.,
		0.,0.,	1.,1.,	0.,0.,
		0.,0.,	0.,0.,	1.,1.,
	};
	std::vector<unsigned> ah;
	const double mts = checker.maxTimeStep(cp, &ah);
	INFO(mts);
	CHECK(mts == 1);
}

TEST_CASE("Standard quadratic tet MTS") {
	ValidityChecker<3, 3, 2> checker;
	const double prec = .1;
    checker.setPrecisionTarget(prec);
    checker.setMaxSubdiv(1);
    const std::vector<double> cp = {
		0.,0.,	0.,0.,	0.,0.,
		1.,1.,	0.,0.,	0.,0.,
		0.,0.,	1.,1.,	0.,0.,
		0.,0.,	0.,0.,	1.,1.,
		.5,.5,	0.,0.,	0.,0.,
		.5,.5,	.5,.5,	0.,0.,
		0.,0.,	.5,.5,	0.,0.,
		0.,0.,	0.,0.,	.5,.5,
		.5,.5,	0.,0.,	.5,.5,
		0.,0.,	.5,.5,	.5,.5,
	};
	std::vector<unsigned> ah;
	const double mts = checker.maxTimeStep(cp, &ah);
	INFO(mts);
	CHECK(mts == 1);
}

TEST_CASE("Invalid quadratic tet MTS") {
	ValidityChecker<3, 3, 2> checker;
	const double prec = .1;
    checker.setPrecisionTarget(prec);
    const std::vector<double> cp = {
		0.,.5,	0.,.5,	0.,.5,
		1.,1.,	0.,0.,	0.,0.,
		0.,0.,	1.,1.,	0.,0.,
		0.,0.,	0.,0.,	1.,1.,
		.5,.5,	0.,0.,	0.,0.,
		.5,.5,	.5,.5,	0.,0.,
		0.,0.,	.5,.5,	0.,0.,
		0.,0.,	0.,0.,	.5,.5,
		.5,.5,	0.,0.,	.5,.5,
		0.,0.,	.5,.5,	.5,.5,
	};
	std::vector<unsigned> ah;
	const double mts = checker.maxTimeStep(cp, &ah);
	INFO(mts);
	CHECK(mts < 1);
}

TEST_CASE("Exact zero at corner quadratic tet MTS") {
	ValidityChecker<3, 3, 2> checker;
	const double prec = .001;
    checker.setPrecisionTarget(prec);
    const std::vector<double> cp = {
		0.,0.,	0.,0.,	0.,0.,
		1.,1.,	0.,0.,	0.,0.,
		0.,0.,	1.,1.,	0.,0.,
		0.,0.,	0.,0.,	1.,1.,
		.5,.5,	0.,0.,	0.,0.,
		.5,.25,	.5,.5,	0.,0.,
		0.,0.,	.5,.5,	0.,0.,
		0.,0.,	0.,0.,	.5,.5,
		.5,.5,	0.,0.,	.5,.5,
		0.,0.,	.5,.5,	.5,.5,
	};
	std::vector<unsigned> ah;
	const double mts = checker.maxTimeStep(cp, &ah);
	INFO(mts);
	CHECK(mts < 1);
	CHECK(mts > 1-prec);
}