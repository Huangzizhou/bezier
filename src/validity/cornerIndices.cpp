#include "cornerIndices.hpp"

namespace element_validity {
template<>
void cornerIndicesT<1, 1, 1>(std::vector<uint> &v) { v = {0,1}; }

template<>
void cornerIndicesT<1, 1, 2>(std::vector<uint> &v) { v = {0,2}; }

template<>
void cornerIndicesT<1, 1, 3>(std::vector<uint> &v) { v = {0,3}; }

template<>
void cornerIndicesT<1, 1, 4>(std::vector<uint> &v) { v = {0,4}; }

template<>
void cornerIndicesT<1, 1, 5>(std::vector<uint> &v) { v = {0,5}; }

template<>
void cornerIndicesT<2, 2, 1>(std::vector<uint> &v) { v = {0,1,2}; }

template<>
void cornerIndicesT<2, 2, 2>(std::vector<uint> &v) { v = {0,2,5}; }

template<>
void cornerIndicesT<2, 2, 3>(std::vector<uint> &v) { v = {0,3,9}; }

template<>
void cornerIndicesT<3, 3, 1>(std::vector<uint> &v) { v = {0,1,2,3}; }

template<>
void cornerIndicesT<3, 3, 2>(std::vector<uint> &v) { v = {0,2,5,9}; }


}
