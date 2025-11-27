# Python file to plot Legendre polynomials and all their derivatives from files, depending on their order
# This file is used to visualize the results of the C++ code in demo_legendre.cpp
# for n = 3, the file is structured like this:
# x P_0 dP_0 P_1 dP_1 P_2

import numpy as np
import matplotlib.pyplot as plt
def plot_legendre_polynomials(filename, n):
    data = np.loadtxt(filename)
    x = data[:, 0]
    
    plt.figure(figsize=(10, 6))
    
    for i in range(n + 1):
        P = data[:, 1 + 2 * i]
        dP = data[:, 2 + 2 * i]
        
        plt.subplot(2, 1, 1)
        plt.plot(x, P, label=f'P_{i}(x)')
        plt.title('Legendre Polynomials')
        plt.xlabel('x')
        plt.ylabel('P_n(x)')
        plt.legend()
        
        plt.subplot(2, 1, 2)
        plt.plot(x, dP, label=f'dP_{i}/dx', linestyle='--')
        plt.title('Derivatives of Legendre Polynomials')
        plt.xlabel('x')
        plt.ylabel("dP_n/dx")
        plt.legend()
    
    plt.tight_layout()
    plt.show()

# Example usage:
plot_legendre_polynomials('build/legendre_output.txt', n=5) # Change n according to the order used in demo_legendre.cpp when generating the file last time
