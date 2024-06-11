#include "utils/globals.hpp"
#include "numbers/Interval.hpp"
#include "utils/combinatorics.hpp"
#include "validity/element_validity.hpp"
#include "utils/Timer.hpp"

int main() {
    element_validity::Interval().init();
    using element_validity::Validity;
    using element_validity::ValidityChecker;
    // switch () {
    // case Validity::valid: std::cout << "valid"; break;
    // case Validity::invalid: std::cout << "invalid"; break;
    // case Validity::uncertain: std::cout << "uncertain"; break;
    // }

    using element_validity::Timer;
    Timer timer;
    ValidityChecker<2, 2, 2> v222;
    v222.setPrecisionTarget(.01);
    v222.setMaxSubdiv(30);
    timer.start();
    const double t = v222.maxTimeStep({
        0,0, 0,0,
        1,1, 0,0,
        .15,.85, .7,.7,
        .5,.5, -.1,-.1,
        .57,.57, .46,.25,
        .43,.43, .25,.46
    });
    std::cout << "t = " << t << std::endl;
    timer.stop();
    std::cout << timer.read<std::chrono::microseconds>() << "us" << std::endl;
    
    return 0;
}
