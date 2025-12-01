#ifndef PENDULUM_HPP
#define PENDULUM_HPP

#include <cmath>
#include <array>

// Simple pendulum dynamics
// state[0] = angle theta
// state[1] = angular velocity dtheta/dt
// parameters: g = 9.81, L = length

struct Pendulum {
    double g;
    double L;

    Pendulum(double g_ = 9.81, double L_ = 1.0)
        : g(g_), L(L_) {}

    // Returns time derivative of the state
    std::array<double,2> operator()(const std::array<double,2>& state) const {
        double theta = state[0];
        double omega = state[1];

        double dtheta_dt = omega;
        double domega_dt = -(g/L) * std::sin(theta);

        return {dtheta_dt, domega_dt};
    }
};

#endif // PENDULUM_HPP
