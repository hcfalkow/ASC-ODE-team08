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

  class ImprovedEuler : public TimeStepper
  {
    Vector<> m_vecf;
    Vector<> m_ytilde;
  public:
    ImprovedEuler(std::shared_ptr<NonlinearFunction> rhs) 
    : TimeStepper(rhs), m_vecf(rhs->dimF()), m_ytilde(rhs->dimX()) {}
    void DoStep(double tau, VectorView<double> y) override
    {
      this->m_rhs->evaluate(y, m_vecf);
      m_ytilde = y + (tau/2.0) * m_vecf;
      this->m_rhs->evaluate(m_ytilde, m_vecf);
      y += tau * m_vecf;
    }
  };

  class CrankNicolson : public TimeStepper
  {
    std::shared_ptr<NonlinearFunction> m_equ;
    std::shared_ptr<Parameter> m_tau;
    std::shared_ptr<ConstantFunction> m_yold;
    std::shared_ptr<ConstantFunction> m_fold;
  public:
    CrankNicolson(std::shared_ptr<NonlinearFunction> rhs)
    : TimeStepper(rhs), m_tau(std::make_shared<Parameter>(0.0))
    {
      m_yold = std::make_shared<ConstantFunction>(rhs->dimX());
      auto ynew = std::make_shared<IdentityFunction>(rhs->dimX());
      m_fold = std::make_shared<ConstantFunction>(rhs->dimF());
      m_equ = ynew - m_yold - m_tau * (m_rhs + m_fold);
    }

    void DoStep(double tau, VectorView<double> y) override
    {
      m_yold->set(y);
      this->m_rhs->evaluate(y, m_fold->get());
      m_tau->set(tau/2.0);
      NewtonSolver(m_equ, y);
    }
  };

  // Explicit RK for arbitrary Butcher tableau (a must be strictly lower triangular)
  class ExplicitRungeKutta : public TimeStepper
  {
    Matrix<> m_a;         // Butcher tableau coefficients
    Vector<> m_b, m_c;
    int m_stages;
    int m_n;              // dimension of ODE system
    Vector<> m_k;         // stores all k^j values concatenated
    Vector<> m_ytmp;
  public:
    ExplicitRungeKutta(std::shared_ptr<NonlinearFunction> rhs,
      const Matrix<> &a, const Vector<> &b, const Vector<> &c) 
    : TimeStepper(rhs), m_a(a), m_b(b), m_c(c),
      m_stages(c.size()), m_n(rhs->dimX()), 
      m_k(m_stages * m_n), m_ytmp(m_n)
    { }

    void DoStep(double tau, VectorView<double> y) override
    {
      m_k = 0.0;
      
      // compute k^j = f(y + tau * sum_{l<j} a_{jl} * k^l)
      for (int j = 0; j < m_stages; j++)
      {
        m_ytmp = y;
        for (int l = 0; l < j; l++)  // lower triangular: sum only over l < j
          m_ytmp += tau * m_a(j,l) * m_k.range(l*m_n, (l+1)*m_n);
        
        m_rhs->evaluate(m_ytmp, m_k.range(j*m_n, (j+1)*m_n));
      }
      
      // y_{n+1} = y_n + tau * sum_j b_j * k^j
      for (int j = 0; j < m_stages; j++)
        y += tau * m_b(j) * m_k.range(j*m_n, (j+1)*m_n);
    }
  };

}

#endif