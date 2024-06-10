#include "globals.hpp"
#include "numbers/Interval.hpp"
#include "utils/combinatorics.hpp"
#include "validity/element_validity.hpp"

int main() {
    element_validity::Interval().init();
    using element_validity::Validity;
    using element_validity::ValidityChecker;
    // switch () {
    // case Validity::valid: std::cout << "valid"; break;
    // case Validity::invalid: std::cout << "invalid"; break;
    // case Validity::uncertain: std::cout << "uncertain"; break;
    // }
    ValidityChecker<1, 1, 1> v111;
    v111.setPrecisionTarget(.1);
    v111.setMaxSubdiv(20);
    std::cout << v111.maxTimeStep({-1,-1,1,5}) << std::endl;
    std::cout << std::endl;

    ValidityChecker<1, 1, 2> v112;
    v112.setPrecisionTarget(.1);
    v112.setMaxSubdiv(20);
    std::cout << v112.maxTimeStep({0,0,1,1,2,2}) << std::endl;
    std::cout << std::endl;

    ValidityChecker<2, 2, 1> v221;
    v221.setPrecisionTarget(.1);
    v221.setMaxSubdiv(20);
    std::cout << v221.maxTimeStep({-1,-1,-1,-1,0,0,1,2,1,1,0,0}) << std::endl;
    std::cout << std::endl;

    ValidityChecker<2, 2, 2> v222;
    v222.setPrecisionTarget(.1);
    v222.setMaxSubdiv(3);
    std::cout << v222.maxTimeStep({
        0,0, 0,0,
        0,0, 1,1,
        0,0, 2,2,
        1,1, 0,0,
        1,1, 1,1,
        2,2, 0,0
    }) << std::endl;
    std::cout << v222.maxTimeStep({
        0,0, 0,0,
        .43,.43, .25,.25,
        .15,.15, .7,.7,
        .5,.5, .1,.1,
        .58,.58, .46,.46,
        1,1, 0,0 
    }) << std::endl;
    std::cout << v222.maxTimeStep({
        0,0, 0,0,
        .43,.43, .25,.46,
        .15,.85, .7,.7,
        .5,.5, .1,.1,
        .58,.58, .46,.25,
        1,1, 0,0 
    }) << std::endl;
    std::cout << std::endl;
    
    return 0;
}
