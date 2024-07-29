#pragma once
#include <vector>
#include <string>
#include <span>
#include "Interval.hpp"

/* Simple class for precomputed matrices */

namespace element_validity {

template<typename T>
class Matrix{
	private:
	// Store each row as a list of column indices and values 
	using Row = std::vector<std::pair<uint, T>>;
	std::vector<Row> store;
	uint s;

	private:
	// Write operation that only works in ascending order of column index.
	void pushToRow(uint i, uint j, const T &v) {
		if (store.at(i).size() == 0 || store.at(i).back().first < j)
			store.at(i).emplace_back(j, v);
		else throw std::runtime_error("Attempted out of order insert");
	}

	public:
	Matrix() {}
	Matrix(const std::vector<lint> &data) { fill(data); }
	void fill(const std::vector<lint> &data);
	inline const Row &getRow(uint i) const { return store.at(i); }
	inline void resize(uint size) { s = size; store.resize(s); }

	// Apply matrix to column vector src and directly into dst 
	void mult(const std::span<const T> src, const std::span<T> dst) const;
	void mult(const std::span<const T> src, std::vector<T> &dst) const {
		dst.resize(s);
		mult(src, std::span<T>(dst));
	}

	// Apply matrix to column vector v
	// inline std::vector<T> mult(const std::vector<T> &src) const {
	// 	std::vector<T> res;
	// 	mult(src, res);
	// 	return res;
	// };

	// Stream output
	template<typename U>
	friend std::ostream& operator<<(std::ostream& ost, const Matrix<U>& m) {
		const U zero = 0;
		for(uint i=0; i<m.s; ++i) {
			const typename Matrix<U>::Row &r = m.getRow(i);
			uint j = 0;
			for (const std::pair<uint, U> &p : r) {
				while (j < p.first) { ost << zero << '\t'; ++j; }
				ost << p.second << '\t'; ++j;
			}
			while (j < m.s) { ost << zero << '\t'; ++j; }
			ost << std::endl;
		}
		return ost;
	}
};
}
