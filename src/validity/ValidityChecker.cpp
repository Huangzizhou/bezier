#include "ValidityChecker.hpp"

namespace element_validity {

// Priority function
bool Subdomain::operator<(const Subdomain &o) const {
	if (time.lower() != o.time.lower())
		return time.lower() > o.time.lower();
	else return incl.lower() > o.incl.lower();
}


// Debug print
std::ostream& operator<<(std::ostream& ost, const Subdomain &s) {
	ost << "t: " << s.time << std::endl;
	ost << "I: " << s.incl << std::endl;
	ost << "Q: ";
	for (const uint x : s.qSequence) ost << x << ", ";
	ost << std::endl;
	ost << "B: ";
	for (const Interval &b : s.B) ost << '\t' << b << std::endl;

	return ost;
}

// Copy the path I took to get here
void Subdomain::copySequence(std::vector<uint> &dst) const {
	dst.clear();
	dst.reserve(depth());
	std::copy(
		qSequence.cbegin(), qSequence.cend(), std::back_inserter(dst)
	);
}

bool CheckerInfo::success() const {
	switch (status) {
	case Status::completed:
	case Status::reachedTarget:
	case Status::globalCondition:
		return true;
	default:
		return false;
	}
}
std::string CheckerInfo::description() const {
	switch (status) {
	case Status::completed:
		return "Processed all intervals";
	case Status::reachedTarget:
		return "Reached target precision";
	case Status::maxDepth:
		return "User termination condition satisfied";
	case Status::globalCondition:
		return "Global early termination condition satisfied";
	case Status::noSplit:
		return "Cannot split due to machine precision";
	default: return "Something is wrong";
	}
}

}