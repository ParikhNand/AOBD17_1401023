import matplotlib.pyplot as plt
import numpy as np

n = np.array([100, 1000, 1500, 10000,16000])
t1 = np.array([1.56, 4.65, 5.403, 14.04, 17.04,])
t2 = np.array([6.23, 7.5, 8.18, 14.82, 18.87])

n2 = np.array([2, 4, 16, 32, 50])
t3 = np.array([6.83, 7.68, 12.124, 17.25, 23.33])
t4 = np.array([6.24, 6.35, 7.902, 8.177, 13.14])

fig = plt.figure()
plt.subplot(211)
plt.plot(n2,t3)
plt.plot(n2,t4)
plt.legend(['PPCA(EM)', 'RPCA'])
plt.suptitle('Comparision of PPCA(EM) and RPCA for Corrupted Values values')
plt.xlabel('% of Corrupted Value')
plt.ylabel('Error in %')


plt.subplot(222)
plt.plot(n,t1)
plt.plot(n,t2)
plt.legend(['PPCA(EM)', 'RPCA'])
plt.suptitle('Comparision of PPCA(EM) and RPCA for Missing values')
plt.xlabel('Number of missing values')
plt.ylabel('Error in %')
plt.show()
fig.savefig('MissingValues.jpg')

