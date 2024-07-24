#pragma once
#include "utils/globals.hpp"
#include <vector>

namespace element_validity {
	template<uint n, uint s, uint p>
	void cornerIndices(std::vector<uint> &);

	template<uint n, uint s, uint p>
	void cornerIndicesT(std::vector<uint> &);
}