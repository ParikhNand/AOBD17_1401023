import numpy as np


A = np.load('ppca_model_1000.npy')

iter_1000 = np.load('generator_1000.npy')
iter_1000 = iter_1000[:,:,:,0]
iter_5000 = np.load('generator_5000.npy')
iter_5000 = iter_5000[:,:,:,0]
iter_8000 = np.load('generator_8000.npy')
iter_8000 = iter_8000[:,:,:,0]

A = A[np.random.randint(0,A.shape[0],size=16),:,:]

#As we have 36058 images the probability of any 256 random images for the generator is 256/36058
#From this 256 images the probability of selecting any 16 pc's which resemble the output is 16/256
#Thus we need to multiply the norm (A-iter_x) with 16/36058, for calculating the error in the generated images

norm_1000 = np.linalg.norm(A-iter_1000)*16/36058
norm_5000 = np.linalg.norm(A-iter_5000)*16/36058
norm_8000 = np.linalg.norm(A-iter_8000)*16/36058


print "Norm for 1000 iterations:", norm_1000
print "Norm for 5000 iterations:", norm_5000
print "Norm for 8000 iterations:", norm_8000

