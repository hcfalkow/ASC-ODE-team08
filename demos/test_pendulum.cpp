#include <iostream>
#include <array>
#include "pendulum.hpp"

int main() { // set the pendulum up
    Pendulum p(9.81, 1.0); // g & L plugged in
    std::array<double,2> state{0.1, 0.0}; // initial vals uf theta & omega

    auto deriv = p(state); // compute derivatives

    // output derivatives
    std::cout << "dtheta/dt = " << deriv[0] << "\n"; // angular velocity
    std::cout << "domega/dt = " << deriv[1] << "\n"; // angular acceleration

    return 0;
}
