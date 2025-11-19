# Welcome to ASC-ODE's documentation!


ASC-ODE is is a C++ library for solving ordinary differential equations (ODEs).
The equation is defined by the right hand side function.
ASC-ODE provides various time-steppers which may be used for odes with right hand sides
given by a function object.

A small demo for solving a mass-spring model as first order ODE
$ \begin{matrix}
y_0^\prime & = & y_1 \\
y_1^\prime & = & -\frac{k}{m} y_0
\end{matrix} $
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

std::cout << 0.0 << "  " << y(0) << " " << y(1) << std::endl;
for (int i = 0; i < steps; i++)
  {
     stepper.DoStep(tau, y);
     std::cout << (i+1) * tau << "  " << y(0) << " " << y(1) << std::endl;
  }
```    


The result of this simulation in phase space is shown below, where a comparison between the three stepper methods, the explicit Euler, the implicit Euler, and the improved Euler method, has been made:

```{image} pictures/massspring_phase.png
:width: 40%
:align: center
```

As the spring mass model does not include any damping, it is expected that the system performs a harmonic oscillation. It was found that if the number of steps was chosen too low for either the explicit Euler or the implicit Euler method, the solution would diverge and either be unstable or yield a damped system. Increasing the number of steps fixed this issue. The improved Euler method was far more robust and could handle lower step numbers and still yield a correct response.



## Installation

install XXX-odesolver it via git-clone:

    git clone https://github.com/my-github-clone/my-ode-solver.git


To configure and build some tests do

    cd my-ode-solver
    mkdir build
    cd build
    cmake ..
    make
    

## Available time-stepping methods are
...




   
