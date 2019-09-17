# -*- coding: utf-8 -*-
"""
Created on Fri Mar  1 19:26:21 2019

@author: ST
"""


from keras.models import Sequential
from keras.layers import Dense, Activation, Dropout, Flatten
from keras.layers import Conv2D, MaxPooling2D
import numpy as np
from keras.utils import to_categorical
from keras.utils import multi_gpu_model
from IPython.display import clear_output
from keras.callbacks import Callback
from matplotlib import pyplot as plt
from keras.callbacks import TensorBoard
from time import time
import pandas as pd
import scipy.io as sio
from sklearn.model_selection import train_test_split
from sklearn.metrics import roc_auc_score
from matplotlib import pyplot
from sklearn.metrics import roc_curve
from keras.utils.vis_utils import plot_model
import os
os.environ["PATH"]+=os.pathsep+'C:\\Program Files (x86)\\graphviz\\release\\bin\\'
from keras.callbacks import EarlyStopping
import h5py
from keras.models import load_model
# load the train, test and cal data from a .mat file
Xy =  h5py.File('D:\\Projects\\P9_EDML\\EDML_trainTestCalData15.mat','r')
arrays ={}
for k,v in Xy.items() :
    arrays[k] = np.array(v)


Xtrain=arrays['trainX']
yTrain=arrays['trainY']
yTrain=yTrain.reshape(yTrain.shape[1],1)
yTrain=to_categorical(yTrain.astype('uint8'),num_classes=2)
Xtrain[np.isnan(Xtrain)]=-2
Xtrain=Xtrain.reshape(Xtrain.shape[0],Xtrain.shape[1],Xtrain.shape[2],1)


XTest=arrays['testX']
yTest=arrays['testY']
yTest=yTest.reshape(yTest.shape[1],1)
yTest=to_categorical(yTest.astype('uint8'),num_classes=2)
XTest[np.isnan(XTest)]=-2
XTest=XTest.reshape(XTest.shape[0],XTest.shape[1],XTest.shape[2],1)

calX=arrays['calX']
calY=arrays['calY']
calY=calY.reshape(calY.shape[1],1)
calY=to_categorical(calY.astype('uint8'),num_classes=2)
calX[np.isnan(calX)]=-2
calX = calX.reshape(calX.shape[0], calX.shape[1], calX.shape[2], 1)

nb_classes=2
batch_size= 128
input_shape=(64,64,1)
nb_neurons=[32,16]
activation='relu'
optimizer='nadam'
count=1
model = Sequential()
    
# Add each layer iteratively.
for i in range(0,2):
    # Need input shape for first layer.
    if i == 0:
        model.add(Conv2D(nb_neurons[i], kernel_size = (3, 3), activation = activation, padding='same', input_shape = input_shape))
    else:
        model.add(Conv2D(nb_neurons[i], kernel_size = (3, 3), activation = activation))
    
    if i < 2: #otherwise we hit zero
        model.add(MaxPooling2D(pool_size=(2, 2)))
    
    model.add(Dropout(0.2))



model.add(Flatten())
model.add(Dense(64, activation = activation))
model.add(Dropout(0.5))
model.add(Dense(nb_classes, activation = 'softmax'))

# uncomment the next line if you want to train the model on gpu
# model=multi_gpu_model(model)
model.compile(loss='categorical_crossentropy',
          optimizer=optimizer,
          metrics=['mse','accuracy'])



print(model.summary())
history=model.fit(Xtrain,yTrain,validation_data=(calX, calY),epochs=500,batch_size=128, callbacks=[EarlyStopping(monitor='val_loss', patience=10, restore_best_weights=True)], verbose=1)
count=count+1
#%% Testing model performance

y_pred_train=model.predict(Xtrain)
fprt, tprt,thresholdst=roc_curve(yTrain[:,1],y_pred_train[:,1])
auct=roc_auc_score(yTrain[:,1], y_pred_train[:,1])

y_pred_cal=model.predict(calX)
fprc, tprc,thresholdsc=roc_curve(calY[:,1],y_pred_cal[:,1])
aucc=roc_auc_score(calY[:,1], y_pred_cal[:,1])

scores=model.evaluate(XTest,yTest)
print(scores[2]*100)

y_pred_keras=model.predict(XTest)
fpr, tpr,thresholds=roc_curve(yTest[:,1],y_pred_keras[:,1])
auc=roc_auc_score(yTest[:,1], y_pred_keras[:,1])
pyplot.plot(fprt,tprt,'k', fpr,tpr,'b',linewidth=1.5)
trstring= "Training auc= {:0.2f} \n N={:d}".format(auct,Xtrain.shape[0])
vstring= "Validation auc= {:0.2f} \n N={:d}".format(auc,XTest.shape[0])
pyplot.legend((trstring,vstring,), loc=3, fontsize=9, frameon=False)
pyplot.xlabel('Probability of false alarm', fontsize=10)
pyplot.ylabel('Probability of detection', fontsize=10)
pyplot.rcParams["figure.figsize"]=(2.5,2.5)
pyplot.xlim(-0.01,1)
pyplot.ylim(0.1)
# saves in your current working directory
pyplot.savefig('MLROC15.tif',dpi=600, bbox_inches='tight')
plot_model(model,to_file='ML_Model15.tif',show_shapes=True,show_layer_names=True)
model.save('MLDeepLearningModel15.h5')

#%% Load independent test data set
Xy =  h5py.File('EDML_Test15.mat','r')
arrays ={}
for k,v in Xy.items() :
    arrays[k] = np.array(v)




Xtest2=arrays['new_testX']
Xtest2[np.isnan(Xtest2)]=-2
Xtest2=Xtest2.reshape(Xtest2.shape[0],Xtest2.shape[1],Xtest2.shape[2],1)
ytest2 = model.predict(Xtest2)


model= load_model('MLDeepLearningModel15.h5')
Xy =  h5py.File('EDML_RestSetTest15.mat','r')
arrays ={}
for k,v in Xy.items() :
    arrays[k] = np.array(v)
Xtest3=arrays['restSet_testX']
Xtest3[np.isnan(Xtest3)]=-2
Xtest3=Xtest3.reshape(Xtest3.shape[0],Xtest3.shape[1],Xtest3.shape[2],1)
ytest3 = model.predict(Xtest3)
