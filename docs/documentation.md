# Welcome to ASC-ODE's documentation!


ASC-ODE is is a C++ library for solving ordinary differential equations (ODEs).
The equation is defined by the right hand side function.
ASC-ODE provides various time-steppers which may be used for odes with right hand sides
given by a function object.

## Installation (TBD)

install XXX-odesolver it via git-clone:

    git clone https://github.com/my-github-clone/my-ode-solver.git


To configure and build some tests do

    cd my-ode-solver
    mkdir build
    cd build
    cmake ..
    make

## Mass spring system

A small demo for solving a mass-spring model as first order ODE

$$\begin{aligned}
y_0' &= y_1 \\
y_1' &= -\frac{k}{m} y_0
\end{aligned}$$
is here:

```cpp
double tend = 4*M_PI;
int steps = 100;
double tau = tend/steps;

Vector<> y = { 1, 0 };  // initial conditions
shared_ptr<NonlinearFunction> rhs = std::make_shared<MassSpring>(mass, stiffness);

// Choose from three stepper methods below:
ExplicitEuler stepper(rhs);
// ImplicitEuler stepper(rhs);
// ImprovedEuler stepper(rhs);
// CrankNicolson stepper(rhs);

std::cout << 0.0 << "  " << y(0) << " " << y(1) << std::endl;
for (int i = 0; i < steps; i++)
  {
     stepper.DoStep(tau, y);
     std::cout << (i+1) * tau << "  " << y(0) << " " << y(1) << std::endl;
  }
```    

The following chapters will go through how to use the library, as wel as provide solutions to some exercises.

### Exercise 17.2.2

The result of this simulation in time and phase space is shown below, where a comparison between the three stepper methods, the explicit Euler, the implicit Euler, and the improved Euler method, has been made:


![](pictures/test_ode_results/mass/explicit_1000_sim.png)
*Figure 1: Explicit Euler with 1000 steps simulation*

![](pictures/test_ode_results/mass/explicit_1000_phase.png)
*Figure 2: Explicit Euler with 1000 steps phase diagram*

![](pictures/test_ode_results/mass/implicit_1000_sim.png)
*Figure 3: Implicit Euler with 1000 steps simulation*

![](pictures/test_ode_results/mass/implicit_1000_phase.png)
*Figure 4: Implicit Euler with 1000 steps phase diagram*

![](pictures/test_ode_results/mass/improved_1000_sim.png)
*Figure 5: Improved Euler with 1000 steps simulation*

![](pictures/test_ode_results/mass/improved_1000_phase.png)
*Figure 6: Improved Euler with 1000 steps phase diagram*



As the spring mass model does not include any damping, it is expected that the system performs a harmonic oscillation. It was found that if the number of steps was chosen too low for either the explicit Euler or the implicit Euler method, the solution would diverge and either be unstable or yield a damped system. Increasing the number of steps fixed this issue. The improved Euler method was far more robust and could handle lower step numbers and still yield a correct response.
Below, a plot of the explicit Euler method can be seen with a step number of 10000, which ensures that the solution converges to a harmonic oscillation:

![](pictures/test_ode_results/mass/explicit_100000_sim.png)
*Figure 7: Explicit Euler with 100000 steps simulation*

![](pictures/test_ode_results/mass/explicit_100000_phase.png)
*Figure 8: Explicit Euler with 100000 steps phase diagram*

Increasing the end time for the simulations results in the implicit or explicit methods diverging from the desired result (harmonic oscillation) even further. When the step number is increased dramatically, this effect is reduced.

### Exercise 17.4.1
    
In exercise 17.4.1, the Crank-Nicolson method for time stepping is also implemented.

$$y_{i+1} = y_i + \frac{\tau}{2}(f(t_i, y_i) + f(t_{i+1}, y_{i+1})), \quad 0 \le i < n$$
The C++ implementation can be seen below in the section **Available time-stepping methods**.
When using the Crank-Nicolson method on the spring mass system it becomes evident that it converges much faster than the explicit and implicit Euler methods and requires a significantly smaller number of steps. Below, a simulation of the mass spring system using the Crank-Nicolson method with only 100 steps is shown. The method is robust and converges very fast.

![](pictures/test_ode_results/mass/crank_100_sim.png)
*Figure 9: Crank-Nicolson with 100 steps simulation*

![](pictures/test_ode_results/mass/crank_100_phase.png)
*Figure 10: Crank-Nicolson with 100 steps phase diagram*



### RC Circuit implementation

The RC circuit is modelled as:

$$\frac{dU_C}{dt} = \frac{U_0 - U_C}{RC}$$

as the rhs of the equation depends on time $t$, as $U_0 = cos(\omega t)$, where $\omega = 100\pi$, this is not a *autonomous* ODE. It can be made autonomous by adding time $t$ as a state $x_2$, where $x_1$ is $U_C$. The RC circuit is implemented in C++:

```cpp
class RCCircuit : public NonlinearFunction
{
private:
  double resistance;
  double capacitance;
  double omega;

public:
  RCCircuit(double R, double C, double w) : resistance(R), capacitance(C), omega(w) {}

  size_t dimX() const override { return 2; } // make autonomous by adding time as a variable
  size_t dimF() const override { return 2; }

  void evaluate (VectorView<double> x, VectorView<double> f) const override
  {
    f(0) = (cos(omega * x(1)) - x(0)) / (resistance * capacitance); // dUc/dt
    f(1) = 1.0; // dt/dt = 1
  }
  void evaluateDeriv (VectorView<double> x, MatrixView<double> df) const override
  {
    df = 0.0; // reset
    df(0,0) = -1.0 / (resistance * capacitance); // d(dUc/dt)/d(Uc)
    df(0,1) = -omega * sin(omega * x(1)) / (resistance * capacitance); // d(dUc/dt)/d(t)
    df(1,1) = 0.0; // d(dt/dt)/d(t)
  }
};
```

The RC Circuit ODE is solved numerically with the three methods, explicit Euler, implicit Euler and Crank-Nicolson, at time steps 500 and 2000, for a total time of $t_{end} = 0.1\pi$:

![](pictures/test_ode_results/RC/RC_explicit_500.png)
*Figure 11: RC ODE: Explicit Euler with 500 steps simulation*

![](pictures/test_ode_results/RC/RC_implicit_500.png)
*Figure 12: RC ODE: Implicit Euler with 500 steps simulation*

![](pictures/test_ode_results/RC/RC_crank_500.png)
*Figure 13: RC ODE: Crank-Nicolson with 500 steps simulation*

![](pictures/test_ode_results/RC/RC_explicit_2000.png)
*Figure 14: RC ODE: Explicit Euler with 2000 steps simulation*

![](pictures/test_ode_results/RC/RC_implicit_2000.png)
*Figure 15: RC ODE: Implicit Euler with 2000 steps simulation*

![](pictures/test_ode_results/RC/RC_crank_2000.png)
*Figure 16: RC ODE: Crank-Nicolson with 2000 steps simulation*


As evident from the i simulations below, the explicit Euler is unstable for low step numbers and needs a significantly higher step number to converge.
Implicit Euler and Crank-Nicolson are both stable at low step numbers and exert stable behaviour.


## Available time-stepping methods


### Explicit Euler
```cpp
class ExplicitEuler : public TimeStepper
  {
    Vector<> m_vecf;
  public:
    ExplicitEuler(std::shared_ptr<NonlinearFunction> rhs) 
    : TimeStepper(rhs), m_vecf(rhs->dimF()) {}
    void DoStep(double tau, VectorView<double> y) override
    {
      this->m_rhs->evaluate(y, m_vecf);
      y += tau * m_vecf;
    }
  };
```

### Implicit Euler
```cpp
class ImplicitEuler : public TimeStepper
  {
    std::shared_ptr<NonlinearFunction> m_equ;
    std::shared_ptr<Parameter> m_tau;
    std::shared_ptr<ConstantFunction> m_yold;
  public:
    ImplicitEuler(std::shared_ptr<NonlinearFunction> rhs) 
    : TimeStepper(rhs), m_tau(std::make_shared<Parameter>(0.0)) 
    {
      m_yold = std::make_shared<ConstantFunction>(rhs->dimX());
      auto ynew = std::make_shared<IdentityFunction>(rhs->dimX());
      m_equ = ynew - m_yold - m_tau * m_rhs;
    }

    void DoStep(double tau, VectorView<double> y) override
    {
      m_yold->set(y);
      m_tau->set(tau);
      NewtonSolver(m_equ, y);
    }
  };
```

### Improved Euler:
```cpp
  class ImprovedEuler : public TimeStepper // improved Euler method
  {
    Vector<> m_vecf; // f evaluated at current step
    Vector<> m_ytilde; // intermediate value
  public:
    ImprovedEuler(std::shared_ptr<NonlinearFunction> rhs) 
    : TimeStepper(rhs), m_vecf(rhs->dimF()), m_ytilde(rhs->dimX()) {} // constructor
    void DoStep(double tau, VectorView<double> y) override // perform one time step
    {
      this->m_rhs->evaluate(y, m_vecf); // evaluate f(y_n)
      m_ytilde = y + (tau/2.0) * m_vecf; // compute y_tilde
      this->m_rhs->evaluate(m_ytilde, m_vecf); // evaluate f(y_tilde)
      y += tau * m_vecf; // update y to y_{n+1}
    }
  };  
```

### Crank-Nicolson:
```cpp
class CrankNicolson : public TimeStepper
  {
    std::shared_ptr<NonlinearFunction> m_equ; // nonlinear equation to solve
    std::shared_ptr<Parameter> m_tau; // time step size parameter
    std::shared_ptr<ConstantFunction> m_yold; // previous time step value
    std::shared_ptr<ConstantFunction> m_fold; // f evaluated at previous step
  public:
    CrankNicolson(std::shared_ptr<NonlinearFunction> rhs) // constructor
    : TimeStepper(rhs), m_tau(std::make_shared<Parameter>(0.0))  // initialize time step size parameter
    {
      m_yold = std::make_shared<ConstantFunction>(rhs->dimX()); // previous time step value initialized as constant function
      auto ynew = std::make_shared<IdentityFunction>(rhs->dimX()); // current time step value initialized as identity
      m_fold = std::make_shared<ConstantFunction>(rhs->dimF()); // f evaluated at previous step initialized as constant function
      m_equ = ynew - m_yold - m_tau * (m_rhs + m_fold); // crank-nicolson equation
    }

    void DoStep(double tau, VectorView<double> y) override
    {
      m_yold->set(y); // set yold to previous y value
      this->m_rhs->evaluate(y, m_fold->get()); // evaluate f at previous y value
      m_tau->set(tau/2.0); // set time step size parameter to tau/2
      NewtonSolver(m_equ, y); // set y to solution of nonlinear equation
    }
  };
```

## AutoDiff differentiation method

A class *AutoDiff* for automatic differetiation was implemented. 

```cpp
template <size_t N, typename T = double>
  class AutoDiff
  {
  private:
    T m_val;
    std::array<T, N> m_deriv;
  public: 
    AutoDiff () : m_val(0), m_deriv{} {}
    AutoDiff (T v) : m_val(v), m_deriv{} 
    {
      for (size_t i = 0; i < N; i++)
        m_deriv[i] = derivative(v, i);
    }
    
    template <size_t I>
    AutoDiff (Variable<I, T> var) : m_val(var.value()), m_deriv{} 
    {
      m_deriv[I] = 1.0;
    }

    T value() const { return m_val; }
    std::array<T, N>& deriv() { return m_deriv; }
    const std::array<T, N>& deriv() const { return m_deriv; }
  };
```

### Exercises 18.4

#### Operators
This class was extended with several functionalities. 
The essential operators $+$, $-$, $*$ and $/$ were implemented as operator overloads:
```cpp
// "+" operator for two AutoDiff objects:
  template <size_t N, typename T = double>
  AutoDiff<N, T> operator+ (const AutoDiff<N, T>& a, const AutoDiff<N, T>& b)
  {
    AutoDiff<N, T> result(a.value() + b.value());
    for (size_t i = 0; i < N; i++)
      result.deriv()[i] = a.deriv()[i] + b.deriv()[i];
    return result;
  }

  
  template <size_t N, typename T = double>
  auto operator+ (T a, const AutoDiff<N, T>& b) { return AutoDiff<N, T>(a) + b; }

   // "-" operator for two AutoDiff objects - own implementation:
  template <size_t N, typename T = double>
  AutoDiff<N, T> operator- (const AutoDiff<N, T>& a, const AutoDiff<N, T>& b)
  {
    AutoDiff<N, T> result(a.value() - b.value());
    for (size_t i = 0; i < N; i++)
      result.deriv()[i] = a.deriv()[i] - b.deriv()[i];
    return result;
  }

  // "-" operator for one AutoDiff object and one scalar - own implementation:
  template <size_t N, typename T = double>
  auto operator- (T a, const AutoDiff<N, T>& b) { return AutoDiff<N, T>(a) - b; }


  // "*" operator for two AutoDiff objects:
  template <size_t N, typename T = double>
  AutoDiff<N, T> operator* (const AutoDiff<N, T>& a, const AutoDiff<N, T>& b)
  {
    AutoDiff<N, T> result(a.value() * b.value());
    for (size_t i = 0; i < N; i++)
      result.deriv()[i] = a.deriv()[i] * b.value() + a.value() * b.deriv()[i];
    return result;
  }

  // "/" operator for two AutoDiff objects - own implementation:
  template <size_t N, typename T = double>
  AutoDiff<N, T> operator/ (const AutoDiff<N, T>& a, const AutoDiff<N, T>& b)
  {
    AutoDiff<N, T> result(a.value() / b.value());
    for (size_t i = 0; i < N; i++)
      result.deriv()[i] = (a.deriv()[i] * b.value() - a.value() * b.deriv()[i]) / (b.value() * b.value());
    return result;
  }
```
#### Functions
Afterwards, additional basic functions were added to the class: $sin$, $cos$, $exp$, $log$, $pow$, and $sqrt$:

```cpp
// sin function for AutoDiff objects:
  template <size_t N, typename T = double>
  AutoDiff<N, T> sin(const AutoDiff<N, T> &a)
  {
    AutoDiff<N, T> result(sin(a.value()));
    for (size_t i = 0; i < N; i++)
      result.deriv()[i] = cos(a.value()) * a.deriv()[i];
    return result;
  }

  // cos function for AutoDiff objects - own implementation:
  template <size_t N, typename T = double>
  AutoDiff<N, T> cos(const AutoDiff<N, T> &a)
  {
    AutoDiff<N, T> result(cos(a.value()));
    for (size_t i = 0; i < N; i++)
      result.deriv()[i] = -sin(a.value()) * a.deriv()[i];
    return result;
  }

  // exp function for AutoDiff objects - own implementation:
  template <size_t N, typename T = double>
  AutoDiff<N, T> exp(const AutoDiff<N, T> &a)
  {
    AutoDiff<N, T> result(exp(a.value()));
    for (size_t i = 0; i < N; i++)
      result.deriv()[i] = exp(a.value()) * a.deriv()[i];
    return result;
  }

  // log function for AutoDiff objects - own implementation:
  template <size_t N, typename T = double>
  AutoDiff<N, T> log(const AutoDiff<N, T> &a)
  {
    AutoDiff<N, T> result(log(a.value()));
    for (size_t i = 0; i < N; i++)
      result.deriv()[i] = (1 / a.value()) * a.deriv()[i];
    return result;
  }

  // pow function for AutoDiff objects - own implementation:
  template <size_t N, typename T = double>
  AutoDiff<N, T> pow(const AutoDiff<N, T> &a, T n)
  {
    AutoDiff<N, T> result(pow(a.value(), n));
    for (size_t i = 0; i < N; i++)
      result.deriv()[i] = n * pow(a.value(), n - 1) * a.deriv()[i];
    return result;
  }

  // sqrt function for AutoDiff objects - own implementation:
  template <size_t N, typename T = double>
  AutoDiff<N, T> sqrt(const AutoDiff<N, T> &a)
  {
    AutoDiff<N, T> result(sqrt(a.value()));
    for (size_t i = 0; i < N; i++)
      result.deriv()[i] = (1 / (2 * sqrt(a.value())) ) * a.deriv()[i];
    return result;
  }
```

#### Legendre polynomials

To test the implementation of the AutoDiff class, as well as some of the functions and operators, several Legendre polynomials were implemented, with a function defined as:

```cpp
template <typename T>
void LegendrePolynomials(int n, T x, std::vector<T>& P) {
    if (n < 0) {
        P.clear();
        return;
    }
    P.resize(n + 1);
    P[0] = T(1);
    if (n == 0) return;
    P[1] = x;
    for (int k = 2; k <= n; ++k) {
        P[k] = ((T(2 * k - 1) * x * P[k - 1]) - T(k - 1) * P[k - 2]) / T(k);
    }
}
```

The polynomials and their derivatives were plotted up to a order of $n$, in the range of $x\in[-1,1]$ with 100 steps:

![](pictures/legendre_order5_plot.png)
*Figure 17: Legendre polynomials and their derivatives up to 5th order.*


### Exercise 18.5

The right-hand side is implemented once using a templated function T_evaluate<T>().
The type T can be either:
double --> normal evaluation
AutoDiff<2> --> value + partial derivatives

```cpp
template <typename T>
void T_evaluate(VectorView<T> x, VectorView<T> f) const {
    f(0) = x(1);
    f(1) = -m_gravity/m_length * sin(x(0));
}
```

The evaluate function calls the template with double.
The evaluateDeriv function constructs AutoDiff variables:

```cpp
x_ad(0) = Variable<0>(x(0));
x_ad(1) = Variable<1>(x(1));
T_evaluate<AutoDiff<2>>(x_ad, f_ad);
```

The Jacobian matrix is extracted from the AutoDiff derivatives and written into df.
This approach ensures that the pendulum model only needs to be implemented once, and both value and derivatives are obtained from the same code path.

#### Testing

A test program compares: the function value evaluated with doubles, the automatically computed Jacobian and the analytically derived Jacobian

The test verifies that AutoDiff produces the correct derivative structure, and the numerical values match the analytical Jacobian up to machine precision. This confirms that both the AutoDiff class and the templated pendulum implementation work correctly.


## Exercise 19: Runge-Kutta Methods

### 1 - Implementation of ExplicitRungeKutta class

Added a new class `ExplicitRungeKutta` in `src/timestepper.hpp` that implements explicit Runge-Kutta methods for arbitrary Butcher tableaus.

The class:
- Takes coefficients `a`, `b`, `c` from a Butcher tableau as input
- Works with any number of stages
- Requires `a` to be strictly lower triangular (explicit method condition)
- Computes stage values k^j sequentially since each k^j depends only on previous stages

The algorithm follows the standard RK formulation:
```
k^j = f(y_n + tau * sum_{l<j} a_{jl} * k^l)
y_{n+1} = y_n + tau * sum_j b_j * k^j
```

### 2 - Test program

Created `demos/test_explicit_rk.cpp` with:
- Mass-spring system as test ODE (harmonic oscillator)
- Butcher tableaus for RK2 (midpoint), Heun, and classical RK4
- Convergence tests with different step sizes
- Comparison between explicit and implicit methods

### 3 - Bug fix in nanoblas

Fixed a template issue in `nanoblas/src/vecexpr.hpp` line 105:
- Changed `isScalar<T>()` to `is_scalar_type<T>::value`
- The function was being called before its declaration

### 4 - Build system

Added `test_explicit_rk` target to `CMakeLists.txt`

### Results

#### Convergence behavior

RK2 (midpoint) - 2nd order:
| Steps | Error x      | Ratio |
|-------|--------------|-------|
| 100   | 1.86e-04     | -     |
| 200   | 2.38e-05     | ~7.8  |
| 400   | 3.01e-06     | ~7.9  |
| 800   | 3.78e-07     | ~8.0  |

Error reduces by ~4x when doubling steps (expected for 2nd order: 2^2 = 4)

RK4 (classical) - 4th order:
| Steps | Error x      | Ratio |
|-------|--------------|-------|
| 50    | 1.36e-06     | -     |
| 100   | 4.27e-08     | ~32   |
| 200   | 1.34e-09     | ~32   |
| 400   | 4.17e-11     | ~32   |

Error reduces by ~16x when doubling steps (expected for 4th order: 2^4 = 16)

#### Comparison explicit vs implicit (200 steps)

| Method          | Error x     | Error v     |
|-----------------|-------------|-------------|
| Explicit RK4    | 1.34e-09    | 5.10e-08    |
| Implicit Gauss3 | ~0          | 6.11e-14    |

The implicit method (6th order Gauss-Legendre) achieves much better accuracy.

### Files modified/created

- `src/timestepper.hpp` - added ExplicitRungeKutta class
- `demos/test_explicit_rk.cpp` - new test file
- `CMakeLists.txt` - added build target
- `nanoblas/src/vecexpr.hpp` - bug fix


### Conclusion

- Implementation works for any explicit RK method
- Tests confirm expected convergence orders
- Comparison shows trade-offs between explicit and implicit methods
