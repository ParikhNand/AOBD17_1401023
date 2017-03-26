from keras.models import Sequential
from keras.layers import Convolution2D, MaxPooling2D
from keras.layers import Activation, Dropout, Flatten, Dense
from keras.optimizers import SGD, RMSprop
from keras.utils.np_utils import to_categorical
from keras.models import model_from_json
import numpy as np
from pandas import read_csv
import cPickle as pickle
import json
import time

json_file = open( 'model.json' , 'r' )
loaded_model_json = json_file.read()
json_file.close()
loaded_model = model_from_json(loaded_model_json)
# load weights into new model
loaded_model.load_weights("model.h5")
print("Loaded model from disk")

dataframe = read_csv("fer2013_train.csv", header=None)
dataset = dataframe.values
X  = dataset[0:9,1]
E = []

for i in range(len(X)):
	Y = np.fromstring(X[i], dtype=int, sep = ' ')
	Y = np.reshape(Y,(48, 48))
	E.append(Y)
X_train = np.array(E)
X_train = X_train.reshape(-1, 1, X_train.shape[1], X_train.shape[2])
X_test = X_train.astype('float32')
#print X_train
print X_test.shape
print loaded_model.predict(X_test, batch_size=1, verbose=0)


