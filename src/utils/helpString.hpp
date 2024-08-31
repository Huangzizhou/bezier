#pragma once

#include <string>

namespace element_validity {
const std::string helpString = R"HELP(
Usage: bin <inputHDF5> [flags]

List of available flags:
-f Set the first element to be processed (inclusive, defaults to 0)
-g If specified, use global query, otherwise use per-element query
-l Set the maximum level of space subdivisions (defaults to 0=unlimited)
-m Number of threads used for multithreading (defaults to 1=no multithreading)
-n Set the number of elements to be processed (defaults to 0=all elements)
-o Set output path for results in csv format
-p Set the target error % for max valid time step search (defaults to 0.01=1%)
-? Show this help
)HELP";
}