#include "globals.hpp"
#include "validity/element_validity.hpp"

int main() {
    using element_validity::Validity;
    switch (element_validity::isValid(2, 2, 1, {0,0,1,0,0,1})) {
    case Validity::valid: std::cout << "valid"; break;
    case Validity::invalid: std::cout << "invalid"; break;
    case Validity::uncertain: std::cout << "uncertain"; break;
    }
    std::cout << std::endl;
    return 0;
}
