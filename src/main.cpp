#include "globals.hpp"
#include "utils/combinatorics.hpp"
#include "validity/element_validity.hpp"

int main() {
    using element_validity::Validity;
    // switch () {
    // case Validity::valid: std::cout << "valid"; break;
    // case Validity::invalid: std::cout << "invalid"; break;
    // case Validity::uncertain: std::cout << "uncertain"; break;
    // }
    element_validity::isValid<1, 1, 1, true>({0,0,1,5});
    std::cout << std::endl;
    element_validity::isValid<1, 1, 2, true>({0,0,2,2,1,2});
    std::cout << std::endl;
    element_validity::isValid<2, 2, 1, true>({0,0,0,0,0,0,1,2,1,1,0,0});
    std::cout << std::endl;
    return 0;
}
