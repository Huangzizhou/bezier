#include <catch2/catch_test_macros.hpp>
#include <element_validity.hpp>

using element_validity::Validity;
using element_validity::ValidityChecker;
using element_validity::ContinuousValidityChecker;

// TEST_CASE("Multiple linear tets static validity") {
// 	Eigen::MatrixXd cp(8, 3);
// 	cp <<
// 		.5,.5,.5,
// 		1.,0.,0.,
// 		0.,1.,0.,
// 		0.,0.,1.,
// 		0.,0.,0.,
// 		1.,0.,0.,
// 		0.,1.,0.,
// 		0.,0.,1.;
// 	std::tuple<unsigned, unsigned, unsigned> counter;
// 	element_validity::StaticChecker<3> check(cp, SHAPETET, 1);
// 	const Validity valid = check.isValid(0, nullptr, nullptr, &counter);
// 	CHECK(valid == Validity::Invalid);
// 	CHECK(std::get<0>(counter)==1);
// 	CHECK(std::get<1>(counter)==1);
// 	CHECK(std::get<2>(counter)==0);
// }

TEST_CASE("Standard linear tet validity") {
	ValidityChecker<3, 3, 1> checker;
    checker.setMaxSubdiv(1);
    const std::vector<double> cp = {
		0.,0.,0.,
 		1.,0.,0.,
 		0.,1.,0.,
 		0.,0.,1.,
	};
	std::vector<unsigned> ah;
	const Validity res = checker.isValid(cp, &ah);
	CHECK(res == Validity::valid);
}

TEST_CASE("Invalid linear tet validity") {
	ValidityChecker<3, 3, 1> checker;
    checker.setMaxSubdiv(1);
    const std::vector<double> cp = {
		.5,.5,.5,
 		1.,0.,0.,
 		0.,1.,0.,
 		0.,0.,1.,
	};
	std::vector<unsigned> ah;
	const Validity res = checker.isValid(cp, &ah);
	CHECK(res == Validity::invalid);
}

TEST_CASE("Standard quadratic tet validity") {
	ValidityChecker<3, 3, 2> checker;
    checker.setMaxSubdiv(1);
    const std::vector<double> cp = {
		0.,0.,0.,
		1.,0.,0.,
		0.,1.,0.,
		0.,0.,1.,
		.5,0.,0.,
		.5,.5,0.,
		0.,.5,0.,
		0.,0.,.5,
		.5,0.,.5,
		0.,.5,.5,
	};
	std::vector<unsigned> ah;
	const Validity res = checker.isValid(cp, &ah);
	CHECK(res == Validity::valid);
}


TEST_CASE("Invalid quadratic tet validity") {
	ValidityChecker<3, 3, 2> checker;
    checker.setMaxSubdiv(1);
    const std::vector<double> cp = {
		.5,.5,.5,
		1.,0.,0.,
		0.,1.,0.,
		0.,0.,1.,
		.5,0.,0.,
		.5,.5,0.,
		0.,.5,0.,
		0.,0.,.5,
		.5,0.,.5,
		0.,.5,.5,
	};
	std::vector<unsigned> ah;
	const Validity res = checker.isValid(cp, &ah);
	CHECK(res == Validity::invalid);
}

TEST_CASE("Difficult quadratic tet validity") {
	ValidityChecker<3, 3, 2> checker;
    checker.setMaxSubdiv(1);
    const std::vector<double> cp = {
		-0.42452634383930121231, 	-0.15199946433181563132, 	-0.52560337937709855626,
		-0.44050913176895295508, 	-0.15470197781589292707, 	-0.53435160512323620629,
		-0.42918452277309032671, 	-0.16301150961501897906, 	-0.52689501546960226807,
		-0.43797994733333217576, 	-0.14995527783228912444, 	-0.53289783632873310992,
		-0.43084313006168462045, 	-0.15188149791866326077, 	-0.53189337930183577807,
		-0.43243098252764933242, 	-0.15859376743924027963, 	-0.5309488096272174662,
		-0.42647517075677771059, 	-0.15685809076376799886, 	-0.5264561675879664282,
		-0.4314661156158717814, 	-0.15098007388019624164, 	-0.52468553005828000301,
		-0.43654101656446037127, 	-0.15442916576536502848, 	-0.53114214085428379519,
		-0.43290118494677543026, 	-0.1583252066831398741, 	-0.52536247464088814407,
	};
	std::vector<unsigned> ah;
	const Validity res = checker.isValid(cp, &ah);
	CHECK(res == Validity::uncertain);
}

TEST_CASE("Standard linear tet MTS") {
	ContinuousValidityChecker<3, 3, 1> checker;
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
	ContinuousValidityChecker<3, 3, 2> checker;
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
	ContinuousValidityChecker<3, 3, 2> checker;
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
	ContinuousValidityChecker<3, 3, 2> checker;
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
