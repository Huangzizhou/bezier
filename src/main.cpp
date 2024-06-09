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
    element_validity::isValidT<1, 1, 1>({-1,-1,1,5});
    std::cout << std::endl;
    element_validity::isValidT<1, 1, 2>({0,0,2,2,1,2});
    std::cout << std::endl;
    element_validity::isValidT<2, 2, 1>({-1,-1,-1,-1,0,0,1,2,1,1,0,0});
    std::cout << std::endl;
    return 0;
}
