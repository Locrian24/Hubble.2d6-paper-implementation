import numpy as np
import pandas as pd
import tensorflow as tf
from glob import glob

class PredictCYP2D6():
  def __init__(self, X=None):
    self.X = X
    self.predictions = self.predict()

  def predict(self):
    dir = './ensemble_models/'
    model_files = glob(dir + "*.model.h5")

    all_predictions = None
    for file in model_files:
      model = tf.keras.models.load_model(file)
      prediction = model.predict(self.X)
      
      if all_predictions = None:
        all_predictions = prediction
      else:
        all_predictions.dstack((prediction, all_predictions))

    return all_predictions

if __name__ == "__main__":
  P = PredictCYP2D6()