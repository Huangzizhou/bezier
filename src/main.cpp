#include "globals.hpp"
#include "utils/combinatorics.hpp"
#include "validity/element_validity.hpp"

int main() {
    using element_validity::Validity;
    using element_validity::ValidityChecker;
    // switch () {
    // case Validity::valid: std::cout << "valid"; break;
    // case Validity::invalid: std::cout << "invalid"; break;
    // case Validity::uncertain: std::cout << "uncertain"; break;
    // }
    ValidityChecker<1, 1, 1> v111;
    std::cout << v111.maxTimeStep({-1,-1,1,5}) << std::endl;
    std::cout << std::endl;
    ValidityChecker<1, 1, 2> v112;
    std::cout << v112.maxTimeStep({0,0,2,2,1,1}) << std::endl;
    std::cout << std::endl;
    ValidityChecker<2, 2, 1> v221;
    std::cout << v221.maxTimeStep({-1,-1,-1,-1,0,0,1,2,1,1,0,0}) << std::endl;
    std::cout << std::endl;
    return 0;
}
