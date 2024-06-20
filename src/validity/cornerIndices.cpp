#include "cornerIndices.hpp"

namespace element_validity {
template<>
void cornerIndicesT<1, 1, 1>(std::vector<uint> &v) { v = {1}; }

template<>
void cornerIndicesT<1, 1, 2>(std::vector<uint> &v) { v = {1,3}; }

template<>
void cornerIndicesT<1, 1, 3>(std::vector<uint> &v) { v = {1,5}; }

template<>
void cornerIndicesT<1, 1, 4>(std::vector<uint> &v) { v = {1,7}; }

template<>
void cornerIndicesT<1, 1, 5>(std::vector<uint> &v) { v = {1,9}; }

template<>
void cornerIndicesT<2, 2, 1>(std::vector<uint> &v) { v = {2}; }

template<>
void cornerIndicesT<2, 2, 2>(std::vector<uint> &v) { v = {2,8,17}; }

template<>
void cornerIndicesT<2, 2, 3>(std::vector<uint> &v) { v = {2,14,44}; }

template<>
void cornerIndicesT<2, 2, 4>(std::vector<uint> &v) { v = {2,20,83}; }

template<>
void cornerIndicesT<3, 3, 1>(std::vector<uint> &v) { v = {3}; }

template<>
void cornerIndicesT<3, 3, 2>(std::vector<uint> &v) { v = {3,15,39,79}; }


}
