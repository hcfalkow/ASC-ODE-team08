// pendulum_ad_test.cpp
// test templated T_evaluate for pendulum
// evaluate with double and tiny AutoDiff<2>

#include <array>
#include <cmath>
#include <iostream>
#include <iomanip>
#include <vector>
#include <cassert>

using std::size_t;

// small AutoDiff forward-mode for 2 variables:
template <size_t N>
struct AutoDiff {
    double val;
    std::array<double, N> d; // derivatives

    // default initializer:
    AutoDiff() : val(0.0) { d.fill(0.0); }

    // construct from double:
    AutoDiff(double v) : val(v) { d.fill(0.0); }

    // variable constructor: index i --> variable (derivative 1)
    static AutoDiff variable(size_t i, double v = 0.0) {
        AutoDiff a(v);
        if (i < N) a.d[i] = 1.0;
        return a;
    }

    // access derivative:
    const std::array<double, N>& deriv() const { return d; }
};

// arithmetic operators:
template <size_t N>
AutoDiff<N> operator+(AutoDiff<N> a, const AutoDiff<N>& b) {
    AutoDiff<N> r;
    r.val = a.val + b.val;
    for (size_t i=0;i<N;i++) r.d[i] = a.d[i] + b.d[i];
    return r;
}
template <size_t N>
AutoDiff<N> operator-(AutoDiff<N> a, const AutoDiff<N>& b) {
    AutoDiff<N> r;
    r.val = a.val - b.val;
    for (size_t i=0;i<N;i++) r.d[i] = a.d[i] - b.d[i];
    return r;
}
template <size_t N>
AutoDiff<N> operator*(AutoDiff<N> a, const AutoDiff<N>& b) {
    AutoDiff<N> r;
    r.val = a.val * b.val;
    for (size_t i=0;i<N;i++) r.d[i] = a.val * b.d[i] + b.val * a.d[i];
    return r;
}
template <size_t N>
AutoDiff<N> operator/(AutoDiff<N> a, const AutoDiff<N>& b) {
    AutoDiff<N> r;
    r.val = a.val / b.val;
    for (size_t i=0;i<N;i++) r.d[i] = (a.d[i]*b.val - a.val*b.d[i]) / (b.val*b.val);
    return r;
}

// unary functions:
template <size_t N> AutoDiff<N> sin(const AutoDiff<N>& x) {
    AutoDiff<N> r;
    r.val = std::sin(x.val);
    double c = std::cos(x.val);
    for (size_t i=0;i<N;i++) r.d[i] = c * x.d[i];
    return r;
}
template <size_t N> AutoDiff<N> cos(const AutoDiff<N>& x) {
    AutoDiff<N> r;
    r.val = std::cos(x.val);
    double m = -std::sin(x.val);
    for (size_t i=0;i<N;i++) r.d[i] = m * x.d[i];
    return r;
}
template <size_t N> AutoDiff<N> operator-(const AutoDiff<N>& a) {
    AutoDiff<N> r;
    r.val = -a.val;
    for (size_t i=0;i<N;i++) r.d[i] = -a.d[i];
    return r;
}

// small vector container:
template <typename T>
struct TinyVec {
    std::vector<T> data; // --> store elements
    TinyVec(size_t n=0) : data(n) {}
    TinyVec(std::initializer_list<T> init) : data(init) {}
    T& operator()(size_t i) { return data.at(i); }
    const T& operator()(size_t i) const { return data.at(i); }
    size_t size() const { return data.size(); }
};

// PendulumAD class:
class PendulumAD_base {
protected:
    double m_length;
    double m_gravity;
public:
    PendulumAD_base(double length, double gravity=9.81) : m_length(length), m_gravity(gravity) {}
    virtual ~PendulumAD_base() = default;
};

class PendulumAD : public PendulumAD_base {
public:
    PendulumAD(double length, double gravity=9.81) : PendulumAD_base(length, gravity) {}

    size_t dimX() const { return 2; } // input dimension
    size_t dimF() const { return 2; } // output dimension

    // evaluate for raw double:
    void evaluate(const TinyVec<double>& x, TinyVec<double>& f) const {
        T_evaluate<double>(x, f);
    }

    // evaluate derivative using AutoDiff<2>:
    void evaluateDeriv(const TinyVec<double>& x, std::vector<std::vector<double>>& J) const {
        TinyVec< AutoDiff<2> > xa(2), fa(2);
        xa(0) = AutoDiff<2>::variable(0, x(0)); // x0 --> variable
        xa(1) = AutoDiff<2>::variable(1, x(1)); // x1 --> variable
        T_evaluate< AutoDiff<2> >(xa, fa);

        // fill Jacobian:
        J.assign(2, std::vector<double>(2,0.0));
        for (size_t i=0;i<2;i++) for (size_t j=0;j<2;j++) J[i][j] = fa(i).deriv()[j];
    }

    // templated evaluation function
    template <typename T>
    void T_evaluate(const TinyVec<T>& x, TinyVec<T>& f) const {
        // f0 = x1 
        // f1 = - (g / L) * sin(x0)
        f(0) = x(1); // = omega
        f(1) = T(- (m_gravity / m_length)) * sin(x(0)); // = -g/L * sin(theta)
    }
};

// Main:
int main() {
    std::cout << std::fixed << std::setprecision(8);

    PendulumAD p(1.0, 9.81);

    TinyVec<double> x(2);
    x(0) = 0.5; // theta
    x(1) = 0.0; // omega

    TinyVec<double> f(2);
    p.evaluate(x, f); // eval. with double

    std::cout << "Evaluation with double:\n";
    std::cout << "f0 = " << f(0) << "   f1 = " << f(1) << "\n\n";

    // compute Jacobian via AutoDiff:
    std::vector<std::vector<double>> J;
    p.evaluateDeriv(x, J);

    std::cout << "Jacobian from AutoDiff:\n";
    for (size_t i=0;i<2;i++){
        for (size_t j=0;j<2;j++) std::cout << J[i][j] << "   ";
        std::cout << "\n";
    }

    // analytic Jacobian (for comparison):
    double J_analytic[2][2] = {{0.0, 1.0}, { - (9.81/1.0) * std::cos(x(0)), 0.0 }};
    std::cout << "\nAnalytic Jacobian:\n";
    for (size_t i=0;i<2;i++){
        for (size_t j=0;j<2;j++) std::cout << J_analytic[i][j] << "   ";
        std::cout << "\n";
    }

    // check: numeric difference
    std::cout << "\nDifference (AutoDiff - Analytic):\n";
    for (size_t i=0;i<2;i++){
        for (size_t j=0;j<2;j++){
            double diff = J[i][j] - J_analytic[i][j];
            std::cout << diff << "   ";
        }
        std::cout << "\n";
    }

    return 0;
}
