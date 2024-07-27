#include "Matrix.hpp"

namespace element_validity {
template class Matrix<Interval>;

/*
Matrices are read from vectors. The first entry is the matrix size,
followed by sequences of {i, j, numerator, denominator} as integers.
*/
template<typename T>
void Matrix<T>::fill(const std::vector<lint> &data) {
	uint s = data.size();
	assert(s>0);
	resize(data.at(0));
	for(uint k=1; k<s; k+=4) {
		assert(data.at(k+0) >= 0 && data.at(k+1) >= 0);
		pushToRow(
			data.at(k+0),
			data.at(k+1),
			static_cast<T>(data.at(k+2)) / data.at(k+3)
		);
	}
}

template<typename T>
void Matrix<T>::mult(const std::span<const T> src, const std::span<T> dst) const {
	assert(src.size() == s);
	assert(dst.size() == s);
	for(uint i=0; i<s; ++i) {
		T acc = 0;
		for(const auto &p : getRow(i)) acc += p.second * src[p.first];
		dst[i] = acc;
	}
}
}