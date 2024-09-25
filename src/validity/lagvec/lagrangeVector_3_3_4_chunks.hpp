#pragma once
#include "numbers/Interval.hpp"
#include <array>
namespace element_validity {
void chunk_3_3_4_0(const std::array<Interval, 105> &cp, const span<Interval> &out);
void chunk_3_3_4_1(const std::array<Interval, 105> &cp, const span<Interval> &out);
void chunk_3_3_4_2(const std::array<Interval, 105> &cp, const span<Interval> &out);
void chunk_3_3_4_3(const std::array<Interval, 105> &cp, const span<Interval> &out);
void chunk_3_3_4_4(const std::array<Interval, 105> &cp, const span<Interval> &out);
}
