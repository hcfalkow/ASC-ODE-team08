#include <iostream>
#include <autodiff.hpp>
#include <vector>


using namespace ASC_ode;


// Compute Legendre polynomials up to degree n at point x
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

int main()
{
  double x = 0.5;
  int n = 5;
  AutoDiff<1> adx = Variable<0>(x);

  std::cout << "adx = " << adx << std::endl;

  // Instantiate vector of autodiff objects to hold Legendre polynomials
  std::vector<AutoDiff<1>> adP;
  LegendrePolynomials(n, adx, adP);
  
  for (int i = 0; i <= n; ++i) {
    std::cout << "P_" << i << "(" << x << ") = " << adP[i] << std::endl;
    std::cout << "Derivative P_" << i << "(" << x << ") = " 
              << derivative(adP[i], 0) << std::endl;
  }



  return 0;
}


