import argparse
import sys
import os
import numpy as np

from encode_to_seq import Encode2Seq
from predict import PredictCYP2D6

def main(vcf, verbose=True):
  # External files
  pwd = sys.path[0]
  EMBEDDINGS = os.path.join(pwd, '../data/embeddings.txt')
  ANNOTATIONS = os.path.join(pwd, '../data/gvcf2seq.annotation_embeddings.csv')
  REF = os.path.join(pwd, '../data/ref.seq')

  encoded_data = Encode2Seq(
      vcf=vcf,
      embedding_file=EMBEDDINGS,
      annotation_file=ANNOTATIONS,
      ref_seq=REF
    )

  predict = PredictCYP2D6(encoded_data.X)
  if verbose:
    for i, sample in enumerate(encoded_data.sample_names):
      print(sample, predict.predictions[i])

  return np.c_[encoded_data.sample_names, predict.predictions]

if __name__ == "__main__":
  parser = argparse.ArgumentParser()
  parser.add_argument("-v", "--vcf", default=None)
  options = parser.parse_args()

  if options.vcf is None:
    print("You need to give a vcf file")
    exit(-1)

  main(options.vcf)
