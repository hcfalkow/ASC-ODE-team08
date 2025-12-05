#include <iostream>
#include <fstream>
#include <cmath>

#include <nonlinfunc.hpp>
#include <timestepper.hpp>
#include <implicitRK.hpp>

using namespace ASC_ode;


class MassSpring : public NonlinearFunction
{
  double mass;
  double stiffness;
public:
  MassSpring(double m, double k) : mass(m), stiffness(k) {}
  size_t dimX() const override { return 2; }
  size_t dimF() const override { return 2; }
  
  void evaluate (VectorView<double> x, VectorView<double> f) const override
  {
    f(0) = x(1);
    f(1) = -stiffness/mass*x(0);
  }
  
  void evaluateDeriv (VectorView<double> x, MatrixView<double> df) const override
  {
    df = 0.0;
    df(0,1) = 1;
    df(1,0) = -stiffness/mass;
  }
};


// Butcher tableau for RK2 midpoint: 2nd order, 2 stages
Matrix<double> RK2a { { 0.0, 0.0 }, { 0.5, 0.0 } };
Vector<> RK2b { 0.0, 1.0 };
Vector<> RK2c { 0.0, 0.5 };

// Butcher tableau for classical RK4: 4th order, 4 stages
Matrix<double> RK4a { 
  { 0.0, 0.0, 0.0, 0.0 },
  { 0.5, 0.0, 0.0, 0.0 },
  { 0.0, 0.5, 0.0, 0.0 },
  { 0.0, 0.0, 1.0, 0.0 }
};
Vector<> RK4b { 1.0/6.0, 1.0/3.0, 1.0/3.0, 1.0/6.0 };
Vector<> RK4c { 0.0, 0.5, 0.5, 1.0 };

// Butcher tableau for Heun: 2nd order, 2 stages
Matrix<double> Heuna { { 0.0, 0.0 }, { 1.0, 0.0 } };
Vector<> Heunb { 0.5, 0.5 };
Vector<> Heunc { 0.0, 1.0 };


void run_test(TimeStepper& stepper, const std::string& name, double tend, int steps)
{
  double tau = tend/steps;
  Vector<> y = { 1, 0 };
  
  std::ofstream outfile("output_" + name + ".txt");
  outfile << 0.0 << "  " << y(0) << " " << y(1) << std::endl;

  for (int i = 0; i < steps; i++)
  {
    stepper.DoStep(tau, y);
    outfile << (i+1)*tau << "  " << y(0) << " " << y(1) << std::endl;
  }
  
  // for mass-spring with m=k=1, exact solution is x=cos(t), v=-sin(t)
  double exact_x = cos(tend);
  double exact_v = -sin(tend);
  double err_x = std::abs(y(0) - exact_x);
  double err_v = std::abs(y(1) - exact_v);
  
  std::cout << name << ": steps=" << steps 
            << " err_x=" << err_x 
            << " err_v=" << err_v << std::endl;
}


int main()
{
  double tend = 2*M_PI;
  auto rhs = std::make_shared<MassSpring>(1.0, 1.0);
  
  std::cout << "Testing explicit Runge-Kutta methods\n" << std::endl;
  
  std::cout << "--- RK2 (midpoint) ---" << std::endl;
  for (int steps : {100, 200, 400, 800})
  {
    ExplicitRungeKutta stepper(rhs, RK2a, RK2b, RK2c);
    run_test(stepper, "rk2_" + std::to_string(steps), tend, steps);
  }

  std::cout << "\n--- Heun ---" << std::endl;
  for (int steps : {100, 200, 400, 800})
  {
    ExplicitRungeKutta stepper(rhs, Heuna, Heunb, Heunc);
    run_test(stepper, "heun_" + std::to_string(steps), tend, steps);
  }

  std::cout << "\n--- RK4 (classical) ---" << std::endl;
  for (int steps : {50, 100, 200, 400})
  {
    ExplicitRungeKutta stepper(rhs, RK4a, RK4b, RK4c);
    run_test(stepper, "rk4_" + std::to_string(steps), tend, steps);
  }
  
  std::cout << "\n--- comparison implicit vs explicit ---" << std::endl;
  int steps = 200;
  {
    ExplicitRungeKutta stepper(rhs, RK4a, RK4b, RK4c);
    run_test(stepper, "explicit_rk4", tend, steps);
  }
  {
    // Gauss3c is defined in implicitRK.hpp (ASC_ode namespace)
    auto [Gauss3a,Gauss3b] = ComputeABfromC(Gauss3c);
    ImplicitRungeKutta stepper(rhs, Gauss3a, Gauss3b, Gauss3c);
    run_test(stepper, "implicit_gauss3", tend, steps);
  }

  return 0;
}