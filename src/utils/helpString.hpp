#pragma once

#include <string>

namespace element_validity {
const std::string helpString = R"HELP(
Usage: bin <inputHDF5> [flags]

List of available flags:
-f Set the first element to be processed (inclusive, defaults to 0=first element)
-g If specified, use global query, otherwise use per-element query
-l Set the maximum level of space subdivisions (defaults to 0=unlimited)
-m Number of threads used for multithreading (defaults to 1=no multithreading)
-n Set the number of elements to be processed (defaults to 0=all elements)
-o Set output path for detailed results in csv format (defaults to cout)
-p Set the target error % for max valid time step search (defaults to 0.01=1%)
-s If specified, run a static pre-check to filter out elements that are not valid at t=0; the argument is the number of maximum subdivisions (defaults to 0=unlimited).
-t If specified, run a static check instead of the dynamic check; the argument is the time in [0,1] to check (defaults to 1=end frame)
-? Show this help
)HELP";
}