#ifndef AUTODIFF_HPP
#define AUTODIFF_HPP

#include <cstddef> 
#include <ostream> 
#include <cmath>   
#include <array>  


namespace ASC_ode
{

  template <size_t N, typename T = double>
  class Variable 
  {
    private:
      T m_val;
    public:
      Variable (T v) : m_val(v) {}
      T value() const { return m_val; }
  };

  template <typename T = double>
  auto derivative (T v, size_t /*index*/) { return T(0); } 


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


  template <size_t N, typename T = double>
  auto derivative (AutoDiff<N, T> v, size_t index) 
  {
    return v.deriv()[index];
  }



  template <size_t N, typename T>
  std::ostream & operator<< (std::ostream& os, const AutoDiff<N, T>& ad)
  {
    os << "Value: " << ad.value() << ", Deriv: [";
    for (size_t i = 0; i < N; i++)
    {
      os << ad.deriv()[i];
      if (i < N - 1) os << ", ";
    }
    os << "]";
    return os;
  }

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


  using std::sin;
  using std::cos;

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

} // namespace ASC_ode

#endif
