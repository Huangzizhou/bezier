#pragma once
#ifdef EIGEN_INTERFACE
#include <Eigen/Core>
#include <vector>
#include <iostream>

namespace element_validity {
std::vector<double> convertEigenMatrix(const Eigen::MatrixXd& mat);
std::vector<double> convertEigenMatrix(
	const Eigen::MatrixXd& mat1,
	const Eigen::MatrixXd& mat2
);
}
#endif