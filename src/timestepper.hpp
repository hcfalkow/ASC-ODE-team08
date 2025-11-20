#ifndef TIMERSTEPPER_HPP
#define TIMERSTEPPER_HPP

#include <functional>
#include <exception>

#include "Newton.hpp"


namespace ASC_ode
{
  
  class TimeStepper
  { 
  protected:
    std::shared_ptr<NonlinearFunction> m_rhs;
  public:
    TimeStepper(std::shared_ptr<NonlinearFunction> rhs) : m_rhs(rhs) {}
    virtual ~TimeStepper() = default;
    virtual void DoStep(double tau, VectorView<double> y) = 0;
  };

  class ExplicitEuler : public TimeStepper
  {
    Vector<> m_vecf; // f evaluated at current step
  public:
    ExplicitEuler(std::shared_ptr<NonlinearFunction> rhs) 
    : TimeStepper(rhs), m_vecf(rhs->dimF()) {}
    void DoStep(double tau, VectorView<double> y) override
    {
      this->m_rhs->evaluate(y, m_vecf);
      y += tau * m_vecf;
    }
  };

  class ImplicitEuler : public TimeStepper
  {
    std::shared_ptr<NonlinearFunction> m_equ; // nonlinear equation to solve
    std::shared_ptr<Parameter> m_tau; // time step size parameter
    std::shared_ptr<ConstantFunction> m_yold; // previous time step value
  public:
    ImplicitEuler(std::shared_ptr<NonlinearFunction> rhs) // constructor
    : TimeStepper(rhs), m_tau(std::make_shared<Parameter>(0.0)) // initialize time step size parameter
    {
      m_yold = std::make_shared<ConstantFunction>(rhs->dimX()); // previous time step value
      auto ynew = std::make_shared<IdentityFunction>(rhs->dimX()); // current time step value
      m_equ = ynew - m_yold - m_tau * m_rhs; // implicit Euler equation
    }

    void DoStep(double tau, VectorView<double> y) override // perform one time step
    {
      m_yold->set(y); // set previous time step value
      m_tau->set(tau); // set time step size
      NewtonSolver(m_equ, y); // solve nonlinear equation for new value
    }
  };

  // Improved Euler:
  // y_tilde = y_n + tau/2 * f(y_n)
  // y_{n+1} = y_n + tau * f(y_tilde)
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

  // Crank-Nicolson:
  // y_{i+1} = y_i + tau/2 * (f(t_i, y_i) + f(t_{i+1}, y_{i+1}))
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


}


#endif
