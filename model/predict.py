import sys
import os
import numpy as np
import pandas as pd
import tensorflow as tf
from glob import glob

class PredictCYP2D6():
  def __init__(self, X=None):
    self.X = X
    self.predictions = None

    self.run()

  def run(self):
    predictions = self.predict(self.X).mean(axis=0)
    self.predictions = self.get_functions(predictions)

  # Cutoff values and logic are from the original repo:
  # https://github.com/gregmcinnes/Hubble2D6
  def get_functions(self, pred):
    cutpoint_1 = 0.4260022
    cutpoint_2 = 0.7360413

    cut1 = np.greater(pred[:, 0], [cutpoint_1])
    cut2 = np.greater(pred[:, 1], [cutpoint_2])

    functions = []
    for i in range(pred.shape[0]):
      if cut1[i] == True and cut2[i] == True:
        functions.append("Normal Function")
      elif cut1[i] == True and cut2[i] == False:
        functions.append("Decreased Function")
      else:
        functions.append("No Function")

    return np.array(functions) 

  def evaluate(self, y):
    labels = self.get_functions(y)

    return np.sum(labels == self.predictions) / len(y)

  def predict(self, X):
    pwd = sys.path[0]
    path = os.path.join(pwd, '../data/ensemble_models/')
    models = glob(path + "*.model.h5")
    
    predictions = []
    for m in models:
      print(f"Predicting using {m}")
      model = tf.keras.models.load_model(m)
      
      pred = model.predict(X)
      predictions.append(pred)

    return np.array(predictions)

if __name__ == "__main__":
  test_X = pd.read_csv('./data/test.samples.csv')
  print(test_X.shape)
  P = PredictCYP2D6()
