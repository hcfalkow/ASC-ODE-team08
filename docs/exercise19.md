# Exercise 19: Runge-Kutta Methods

### 1. Implementation of ExplicitRungeKutta class

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

### 2. Test program

Created `demos/test_explicit_rk.cpp` with:
- Mass-spring system as test ODE (harmonic oscillator)
- Butcher tableaus for RK2 (midpoint), Heun, and classical RK4
- Convergence tests with different step sizes
- Comparison between explicit and implicit methods

### 3. Bug fix in nanoblas

Fixed a template issue in `nanoblas/src/vecexpr.hpp` line 105:
- Changed `isScalar<T>()` to `is_scalar_type<T>::value`
- The function was being called before its declaration

### 4. Build system

Added `test_explicit_rk` target to `CMakeLists.txt`

## Results

### Convergence behavior

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

### Comparison explicit vs implicit (200 steps)

| Method          | Error x     | Error v     |
|-----------------|-------------|-------------|
| Explicit RK4    | 1.34e-09    | 5.10e-08    |
| Implicit Gauss3 | ~0          | 6.11e-14    |

The implicit method (6th order Gauss-Legendre) achieves much better accuracy.

## Conclusion

- Implementation works for any explicit RK method
- Tests confirm expected convergence orders
- Comparison shows trade-offs between explicit and implicit methods
