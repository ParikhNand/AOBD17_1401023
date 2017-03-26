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

seed = 7
np.random.seed(seed)
dataframe = read_csv("fer2013.csv", header=None)
dataset = dataframe.values
X  = dataset[1:28710,1]
E = []
E_y = []

for i in range(len(X)):
	Y = np.fromstring(X[i], dtype=int, sep = ' ')
	Y = np.reshape(Y,(48, 48))
	E.append(Y)
X_train = np.array(E)
X_train = X_train.reshape(-1, 1, X_train.shape[1], X_train.shape[2])
X_train = X_train.astype('float32')
#print X_train
print X_train.shape

Y = dataset[1:28710,0]

for i in range(len(Y)):
	Y_tr = np.fromstring(Y[i], dtype=int, sep=' ')
	E_y.append(Y_tr)

Y_train = np.array(E_y)
dummy_y = to_categorical(Y_train)


batch_size = 128
nb_epoch = 55


model = Sequential()
model.add(Convolution2D(32, 3, 3, border_mode='same', activation='relu',
                        input_shape=(1, X_train.shape[2], X_train.shape[3])))
model.add(Convolution2D(32, 3, 3, border_mode='same', activation='relu'))
model.add(Convolution2D(32, 3, 3, border_mode='same', activation='relu'))
model.add(MaxPooling2D(pool_size=(2, 2), dim_ordering="th"))

model.add(Convolution2D(64, 3, 3, border_mode='same', activation='relu'))
model.add(Convolution2D(64, 3, 3, border_mode='same', activation='relu'))
model.add(Convolution2D(64, 3, 3, border_mode='same', activation='relu'))
model.add(MaxPooling2D(pool_size=(2, 2),dim_ordering="th"))

model.add(Convolution2D(128, 3, 3, border_mode='same', activation='relu'))
model.add(Convolution2D(128, 3, 3, border_mode='same', activation='relu'))
model.add(Convolution2D(128, 3, 3, border_mode='same', activation='relu'))
model.add(MaxPooling2D(pool_size=(2, 2),dim_ordering="th"))

model.add(Flatten()) 
model.add(Dense(64, activation='relu'))
model.add(Dense(64, activation='relu'))
model.add(Dense(7, activation='sigmoid'))


model.compile(loss='categorical_crossentropy', optimizer='adam', metrics=['accuracy'])
model.fit(X_train, dummy_y, nb_epoch=nb_epoch, batch_size=batch_size, validation_split=0.2, shuffle=True)



model_json = model.to_json()
with open("model.json", "w") as json_file:
	json_file.write(model_json)
# serialize weights to HDF5
model.save_weights("model.h5")
print("Saved model to disk")


