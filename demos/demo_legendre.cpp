#include <iostream>
#include <fstream>
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
  // calculate lengendre polynomials of order n and their derivatives at x for x in range [-1,1] with k steps.
  // print the results to a file "legendre_output.txt"
  double x_start = -1.0;
  double x_end = 1.0;
  int n = 5; // order of legendre polynomial
  int k = 100; // number of steps
  double step = (x_end - x_start) / k;
  std::ofstream outfile("legendre_output.txt");
  for (int i = 0; i <= k; ++i) {
    double x = x_start + i * step;
    AutoDiff<1> adx = Variable<0>(x);
    std::vector<AutoDiff<1>> adP;
    LegendrePolynomials(n, adx, adP);
    outfile << x;
    for (int j = 0; j <= n; ++j) {
      outfile << " " << adP[j].value() << " " << derivative(adP[j], 0);
      std::cout << adP[j].value() << "\n";
    }
    outfile << std::endl;
  }
  outfile.close();



  return 0;
}


