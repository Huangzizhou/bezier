#include "Matrix.hpp"

namespace element_validity {
template class Matrix<Interval>;

#ifndef DATA_DIRECTORY
#warning "DATA_DIRECTORY must be defined"
#endif

/*
Assume that matrices are in a raw text file, each line representing a matrix.
Each line starts with an ID string, followed by a flattened list of tuples
{i,j,numerator,denominator} for the entries.
*/
template<typename T>
void Matrix<T>::read(const std::string &fileID, uint line, uint size) {
	std::string filePath = "DATA_DIRECTORY";
	filePath += fileID;
	std::ifstream in(filePath);
	if (!in)
		throw std::runtime_error("Could not open " + filePath);

	constexpr auto M = std::numeric_limits<std::streamsize>::max();
	for(uint l=1; l<=line; ++l)
		in.ignore(M,'\n');
	std::string data;
	std::getline(in, data);
	std::istringstream datastr(data);

	resize(size);
	uint i, j;
	int num, den;
	while(datastr >> i >> j >> num >> den)
		pushToRow(i, j, static_cast<T>(num) / den);
}

template<typename T>
void Matrix<T>::mult(const std::vector<T> &src, std::vector<T> &dst) const {
	assert(src.size() == s);
	dst.resize(s);
	for(uint i=0; i<s; ++i) {
		T acc = 0;
		for(const auto &p : getRow(i)) acc += p.second * src.at(p.first);
		dst.at(i) = acc;
	}
}
}