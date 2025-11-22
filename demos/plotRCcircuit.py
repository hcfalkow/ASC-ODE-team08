import numpy as np
data = np.loadtxt('build/output_test_ode.txt', usecols=(0, 1, 2))
# print (data)

import matplotlib.pyplot as plt

plt.plot(data[:,0], data[:,1], label='Uc')
# plt.plot(data[:,0], data[:,2], label='t')
plt.xlabel('time')
plt.ylabel('value')
plt.title('RC Circuit Time Evolution')
plt.legend()
plt.grid()
plt.show()


# plt.plot(data[:,1], data[:,2], label='phase plot')
# plt.xlabel('Uc')
# plt.ylabel('t')
# plt.title('RC Circuit Phase Plot')
# plt.legend()
# plt.grid()
# plt.show()

