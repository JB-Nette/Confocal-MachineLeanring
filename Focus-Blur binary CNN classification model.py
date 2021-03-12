# -*- coding: utf-8 -*-
# -*- coding: utf-8 -*-
import numpy as np
import os
import cv2
import matplotlib.pyplot as plt
import keras
import time
from sklearn import preprocessing
from sklearn.model_selection import train_test_split
from sklearn.metrics import confusion_matrix

from keras.models import Sequential
from keras.models import load_model
from keras.layers import Dense, Flatten, BatchNormalization
from keras.layers.convolutional import Conv2D, MaxPooling2D

from keras.callbacks import EarlyStopping
from keras.callbacks import ModelCheckpoint

import seaborn as sb
# Mount google drive to access images from Google Drive
# from google.colab import drive
# drive.mount('/content/drive')

# path of different image set and load them in numpy array x and y.
# Focus image is labelled as 1, Blur as 2
# Preprocessing scale is done on all images such that the values are centered to mean with variance 1.

x = []
y = []

path_f = "/content/drive/My Drive/PATH_TO_OUTPUT/20201125-confocal-POR sample/Focus_cy5_resize/"
path_b = "/content/drive/My Drive/PATH_TO_OUTPUT/20201125-confocal-POR sample/Blur_cy5_resize/"
path_f1 = "/content/drive/My Drive/PATH_TO_OUTPUT/20201125-confocal-POR sample/Focus_resize/"
path_b1 = "/content/drive/My Drive/PATH_TO_OUTPUT/20201125-confocal-POR sample/Blur_resize/"
path_f2 = "/content/drive/My Drive/PATH_TO_OUTPUT/20201125-confocal-POR sample/Focus/"
path_b2 = "/content/drive/My Drive/PATH_TO_OUTPUT/20201125-confocal-POR sample/Blur/"
valid_format = ['.tif']

for f in os.listdir(path_f):
    ext = os.path.splitext(f)[1]
    if ext.lower() not in valid_format:
        continue
    img = cv2.imread(os.path.join(path_f, f), -1)
    imgResize = cv2.resize(img, (256, 256))
    img = np.array(imgResize)
    img = preprocessing.scale(img)
    img = img.reshape(256, 256, 1)
    x.append(img)
    y.append(1)

for f in os.listdir(path_b):
    ext = os.path.splitext(f)[1]
    if ext.lower() not in valid_format:
        continue
    img = cv2.imread(os.path.join(path_b, f), -1)
    imgResize = cv2.resize(img, (256, 256))
    img = np.array(imgResize)
    img = preprocessing.scale(img)
    img = img.reshape(256, 256, 1)
    x.append(img)
    y.append(0)

for f in os.listdir(path_f1):
    ext = os.path.splitext(f)[1]
    if ext.lower() not in valid_format:
        continue
    img = cv2.imread(os.path.join(path_f1, f), -1)
    imgResize = cv2.resize(img, (256, 256))
    img = np.array(imgResize)
    img = preprocessing.scale(img)
    img = img.reshape(256, 256, 1)
    x.append(img)
    y.append(1)

for f in os.listdir(path_b1):
    ext = os.path.splitext(f)[1]
    if ext.lower() not in valid_format:
        continue
    img = cv2.imread(os.path.join(path_b1, f), -1)
    imgResize = cv2.resize(img, (256, 256))
    img = np.array(imgResize)
    img = preprocessing.scale(img)
    img = img.reshape(256, 256, 1)
    x.append(img)
    y.append(0)

for f in os.listdir(path_f2):
    ext = os.path.splitext(f)[1]
    if ext.lower() not in valid_format:
        continue
    img = cv2.imread(os.path.join(path_f2, f), -1)
    imgResize = cv2.resize(img, (256, 256))
    img = np.array(imgResize)
    img = preprocessing.scale(img)
    img = img.reshape(256, 256, 1)
    x.append(img)
    y.append(1)

for f in os.listdir(path_b2):
    ext = os.path.splitext(f)[1]
    if ext.lower() not in valid_format:
        continue
    img = cv2.imread(os.path.join(path_b2, f), -1)
    imgResize = cv2.resize(img, (256, 256))
    img = np.array(imgResize)
    img = preprocessing.scale(img)
    img = img.reshape(256, 256, 1)
    x.append(img)
    y.append(0)

x = np.array(x)
y = np.array(y)

# Splitting of train, test and val set randomly
x_train, x_test, y_train, y_test = train_test_split(x, y, test_size=0.2, random_state=42)
x_train, x_val, y_train, y_val = train_test_split(x_train, y_train, test_size=0.1, random_state=42)

# Printing the distribution of Class 0 and 1 for train, test and val set.
(unique, counts) = np.unique(y_train, return_counts=True)
for i in range(len(unique)):
    print(unique[i], ":", counts[i])

(unique, counts) = np.unique(y_test, return_counts=True)
for i in range(len(unique)):
    print(unique[i], ":", counts[i])

(unique, counts) = np.unique(y_val, return_counts=True)
for i in range(len(unique)):
    print(unique[i], ":", counts[i])

# Create CNN model using Sequential()
model = Sequential()
model.add(Conv2D(32, (5, 5), input_shape=(256, 256, 1), activation='relu'))
model.add(MaxPooling2D(pool_size=(3, 3)))
model.add(Conv2D(64, (5, 5), activation='relu'))
model.add(MaxPooling2D(pool_size=(3, 3)))
model.add(Conv2D(64, (3, 3), activation='relu'))
model.add(MaxPooling2D(pool_size=(3, 3)))
model.add(BatchNormalization())
model.add(Flatten())
model.add(Dense(64, activation='relu'))
model.add(Dense(1, activation='sigmoid'))
model.compile(optimizer=keras.optimizers.SGD(lr=0.0001, decay=0.0001), loss='binary_crossentropy', metrics=['accuracy'])
model.summary()

start_time = time.time()
# simple early stopping
es = EarlyStopping(monitor='val_loss', mode='min', verbose=1, patience=30)
mc = ModelCheckpoint('/content/drive/My Drive/PATH_TO_OUTPUT/cdmodelv1', monitor='val_accuracy', mode='max', verbose=1,
                     save_best_only=True)
# fit model
result = model.fit(x_train, y_train, epochs=1000, validation_data=(x_test, y_test), verbose=0, callbacks=[es, mc])
# load the saved model
saved_model = load_model('/content/drive/My Drive/PATH_TO_OUTPUT/cdmodelv1')
# evaluate the model
_, train_acc = saved_model.evaluate(x_train, y_train, verbose=0)
_, test_acc = saved_model.evaluate(x_test, y_test, verbose=0)
print('Train: %.3f, Test: %.3f' % (train_acc, test_acc))
end_time = time.time()
print("--- %s seconds ---" % round((end_time - start_time), 3))


# Plotting validation and training accuracy
_, ax = plt.subplots()
ax.plot(result.history['accuracy'])
ax.plot(result.history['val_accuracy'])
_ = ax.set_xlabel('Epochs')
_ = ax.set_ylabel('Accuracy')
_ = ax.set_title('Deep learning Blur vs Sharp')
_ = ax.legend(["training data", "test data"])
_ = ax.set_ylim([0, 1])
ax.grid()
plt.show()

# Plotting validation and training loss
_, ax = plt.subplots()
ax.plot(result.history['loss'])
ax.plot(result.history['val_loss'])
_ = ax.set_xlabel('Epochs')
_ = ax.set_ylabel('Loss')
_ = ax.set_title('Deep learning Blur vs Sharp')
_ = ax.legend(["training data", "test data"])
_ = ax.set_ylim([0, 1])
ax.grid()
plt.show()

# Evaluate performance on x_val
_, val_acc = saved_model.evaluate(x_val, y_val, verbose=0)
print(val_acc)


# Create the prediction function where prediction of less than 0.5 is classified as blur while more than 0.5 is focus.
def results(mdl, xsample, ysample):
    scores = []
    ypred = mdl.predict(xsample)
    for pred in ypred:
        if pred >= 0.5:
            score = 1
        else:
            score = 0
        scores.append(score)

    scores = np.array(scores)
    prediction = scores == ysample
    truecount = sum(prediction)
    acc = truecount/len(result)
    return acc


# Visualising the confusion matrix
sb.heatmap(confusion_matrix(y_val, results(saved_model, x_val)), annot=True, annot_kws={"size": 16}, fmt='g')

