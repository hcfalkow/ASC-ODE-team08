#include <iostream>
#include <fstream> 

#include <nonlinfunc.hpp>
#include <timestepper.hpp>

using namespace ASC_ode;


class MassSpring : public NonlinearFunction
{
private:
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

// RC circuit model:
// dUc/dt = (U0 - Uc)/(R*C)
// U0 = cos(omega*t)
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


int main()
{
  double tend = 0.1*M_PI;
  int steps = 2000;
  double tau = tend/steps;

  Vector<> y = { 1, 0 };  // initializer list
  // auto rhs = std::make_shared<MassSpring>(1.0, 1.0);
  auto rhs = std::make_shared<RCCircuit>(100, 1e-6, 100*M_PI);
  // auto rhs = std::make_shared<RCCircuit>(1, 1, 100*M_PI);
  
  // ExplicitEuler stepper(rhs);
  // ImplicitEuler stepper(rhs);
  // ImprovedEuler stepper(rhs);
  CrankNicolson stepper(rhs);

  std::ofstream outfile ("output_test_ode.txt");
  std::cout << 0.0 << "  " << y(0) << " " << y(1) << std::endl;
  outfile << 0.0 << "  " << y(0) << " " << y(1) << std::endl;

  for (int i = 0; i < steps; i++)
  {
     stepper.DoStep(tau, y);

     std::cout << (i+1) * tau << "  " << y(0) << " " << y(1) << std::endl;
     outfile << (i+1) * tau << "  " << y(0) << " " << y(1) << std::endl;
  }
}
