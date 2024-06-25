#pragma once
#include <map>
#include <unordered_map>
#include <string>
#include <iostream>
#include "globals.hpp"
#include "helpString.hpp"

namespace element_validity {
class Settings {
	public:
	std::string filePath = "";
	std::string resultsPath = "";
	bool globalQuery = false;
	uint maxIterations = 0;
	uint numElem = 0;
	uint firstElem = 0;
	uint numThreads = 16;
	fp_t precision = .01;

	// Argument parser
	Settings(int argc, char** argv);
};
}