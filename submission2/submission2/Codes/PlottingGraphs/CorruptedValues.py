import matplotlib.pyplot as plt
import numpy as np

n = np.array([2, 4, 16, 32, 50])

t1 = np.array([6.83, 7.68, 12.124, 17.25, 23.33])

t2 = np.array([6.24, 6.35, 7.902, 8.177, 13.14])


fig = plt.figure()
plt.plot(n,t1)
plt.plot(n,t2)
plt.legend(['PPCA(EM)', 'RPCA'])
plt.suptitle('Comparision of PPCA(EM) and RPCA for Corrupted Values values')
plt.xlabel('% of Corrupted Value')
plt.ylabel('Error in %')
plt.show()
fig.savefig('Corruptedvalues.jpg')
