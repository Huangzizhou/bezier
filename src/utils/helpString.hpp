#pragma once

#include <string>

namespace element_validity {
const std::string helpString = R"HELP(
Usage: bin resultsFilePath [flags]
Enumerated sequences of paths must be specified using the "#" character as a wildcard.

List of available flags:
-f Set the first element to be processed (inclusive, defaults to 0)
-g If specified, use global query
-i Set the maximum number of space subdivisions (defaults to 0=unlimited)
-l Set the last element to be processed (exclusive, defaults to +inf)
-m Number of threads used for multithreading (defaults to 16)
-o Set output path for results
-p Set the target error % for max valid time step search (defaults to 0.01=1%)
-? Show this help
)HELP";
}