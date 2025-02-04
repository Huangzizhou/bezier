#include "Settings.hpp"

namespace element_validity {

Settings::Settings(int argc, char** argv) {
	char flag = ' ';
	bool infoHelp = false;
	for (int a = 1; a < argc; ++a) {
		std::string s(argv[a]);	// Get the a-th argument
		if (s.front() == '-') {	// flag arguments
			flag = s.at(1);
			switch (flag) {
			case 'f':
			case 'l':
			case 'm':
			case 'n':
			case 'o':
			case 'p':
				break;
			case 'g': globalQuery = true; break;
			case 's': preCheck = true; break;
			case 't': staticCheck = true; break;
			case '?':
				abort = true;
				std::cout << helpString << std::endl;
				break;
			default: 
				infoHelp = true;
				std::cout << "Ignoring flag \"-" << flag <<
					"\" because I don't recognize it." << std::endl;
				break;
			}
		} else {	// flag additional parameters
			switch (flag) {
			case ' ': filePath = s; break;
			case 'f': firstElem = std::stoi(s); break;
			case 'l': maxIterations = std::stoi(s); break;
			case 'm': numThreads = std::stoi(s); break;
			case 'n': numElem = std::stoi(s); break;
			case 'o': resultsPath = s; break;
			case 'p': precision = std::stod(s); break;
			case 's': preCheckMaxIter = std::stoi(s); break;
			case 't': staticCheckTime = std::stod(s); break;
			default:
				infoHelp = true;
				std::cout << "Ignoring argument \"" << s << "\" because "
					"flag \"-" << flag << "\" doesn't need any more arguments."
					<< std::endl;
				break;
			}
		}
	}
	if (filePath.size() == 0) infoHelp = true;
	if (infoHelp && !abort) {
		std::cout << "No input file provided. Use -? for help." << std::endl;
		abort = true;
	}
}

}