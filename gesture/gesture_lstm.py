
# coding: utf-8

# In[1]:

import os
import numpy as np

import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt


# In[2]:

import keras
from keras.models import Sequential
from keras.layers import LSTM, Dense, Dropout
from keras.utils import np_utils


# In[3]:

import exp_loader


# In[4]:

all_data = exp_loader.load_all_experiments()


# In[6]:

seq_len = 120
num_features = 6 # acc(X, Y, Z), gyro(X, y, Z)
num_samples = len(all_data) # 437


# In[7]:

sign_to_int = {
    's0':0,
    's1':1,
    's2':2,
    's3':3,
    's4':4,
    's5':5,
    's6':6,
    's7':7,
    's8':8,
    's9':9,
    's10':10
}


# ## Prepare data

# In[18]:

data_x = []
data_y = []
for j in range(len(all_data)):
    entry = all_data[j]
    entry_x = [] # should be 120 * 4
    for k in range(120):
        entry_acc =  entry['acc']
        entry_gyr = entry['gyro']
        entry_x.append([entry_acc[k,2], entry_acc[k, 3], entry_acc[k,4],
                        entry_gyr[k,2], entry_gyr[k,3], entry_gyr[k,4]])
    data_x.append(entry_x)
    data_y.append([sign_to_int[entry['sign']]])

# data_x should be 437*120* 4


# In[19]:

X = np.array(data_x)
y = np_utils.to_categorical(np.array(data_y))


# ## Build model

# In[20]:

model = Sequential()
model.add(LSTM(512, input_shape=(120, 6), return_sequences=True))
model.add(Dropout(0.5))
model.add(LSTM(512, return_sequences=False))
model.add(Dropout(0.5))
model.add(Dense(y.shape[1], activation='softmax'))
model.compile(loss='categorical_crossentropy', optimizer='rmsprop', metrics=['accuracy'])


# In[ ]:

model.fit(X, y, nb_epoch=100, verbose=2, validation_split=0.1)


# In[ ]:



