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
			case 'e':
			case 'f':
			case 'i':
			case 'm':
			case 'o':
			case 'p':
				break;
			case 'g': globalQuery = true; break;
			case '?':
				std::cout << helpString << std::endl; break;
			default: 
				infoHelp = true;
				std::cout << "Ignoring flag \"-" << flag <<
					"\" because I don't recognize it." << std::endl;
				break;
			}
		} else {	// flag additional parameters
			switch (flag) {
			case ' ': filePath = s; break;
			case 'e': numElem = std::stoi(s); break;
			case 'f': firstElem = std::stoi(s); break;
			case 'i': maxIterations = std::stoi(s); break;
			case 'm': numThreads = std::stoi(s); break;
			case 'o': resultsPath = s; break;
			case 'p': precision = std::stod(s); break;
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
	if (infoHelp) std::cout << "Use -? for available flags." << std::endl;;
}

}