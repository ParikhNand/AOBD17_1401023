import matplotlib.pyplot as plt
import matplotlib.image as mpimg
from numpy.linalg import inv
import numpy as np

import skimage
from skimage import color
from ppca_mv import ppca_mv
from sklearn.metrics import mean_squared_error
from math import sqrt


gray = mpimg.imread("CameraMan.png")
r = np.random.randint(256,size=(5000,2))
r1,c1 = r.shape

for i in range(0,r1):
    gray[r[i,0],r[i,1]] = np.nan

hidden = np.isnan(gray)
missing = np.count_nonzero(hidden)

[a,b,c,d,e]=ppca_mv(gray,20,1)
plt.imshow(e, cmap = plt.cm.Greys_r)
plt.show()


gray[hidden] = 0
rms = sqrt(mean_squared_error(gray, e))

print rms


