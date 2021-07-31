import sys
import os
import numpy as np
import pandas as pd
# Since simulated diplotypes are based on
# existing star alleles, no need for custom annotation

# take in vcf?
# create embeddings
#   - use *1 as reference

#   - For each sample in the vcf
#   - Look up embeddings for all listed variants - need external file
#   - identify any variants without embeddings
#   - Get the index where that embedding belongs - need external file
#   - Replace the embedding at the location
#   - output seq
#   - pass through model
#   - have a force options that will ignore indels or something.  not needed in v1

# This pre-processing was forked from the original author's published code
# Only main change is the encoding of diplotypes
# https://github.com/gregmcinnes/Hubble2D6/blob/master/bin

class FirstStep2Seq():
  def __init__(self, vcf=None, label_csv=None, annotations=None, embeddings_file=None):
    self.vcf = vcf
    self.label_csv = label_csv
    self.annotations = annotations
    self.embeddings_file = embeddings_file

    self.sample_names = None
    self.X = None
    self.y = None

    self.embeddings = self.create_embedding_matrix()

    seqs = self.build_seqs()
    data = self.format_seqs(seqs)

    self.X = data["X"]
    self.y = data["y"]
    self.sample_names = data["sample_names"]

  def format_seqs(self, seq_data):
    sample_names = []
    seqs = []
    null_vector = 2048

    # Create a full array of encoded seqs for each haplotype per sample
    # Haplotypes per sample are seperated by 32-unit empty space
    for k in seq_data.keys():
      sample_names.append(k)
      seq = seq_data[k]
      null = [null_vector] * 32
      full_seq = seq[0].copy() + null + seq[1].copy()
      seqs.append(full_seq)

    # Turn vector of encodings into one-hot seq w/ annotations
    X_ind = np.array(seqs, dtype=int)
    X = self.indices2embeddings(X_ind)
    y = self.get_labels(sample_names)

    data = {
      "sample_names": sample_names, # Guaranteed to be order of X
      "X": X,
      "y": y
    }

    return data

  def get_labels(self, samples):
    y_df = pd.read_csv(self.label_csv, header=None, index_col=0, usecols=[0, 6], names=['name', 'pt'])
    try:
      y = y_df.loc[samples].values
    except KeyError:
      print("Mismatching labels in X and y files")
      exit(1)

    return y

  def create_embedding_matrix(self):
    headers = pd.read_csv(self.annotations, nrows=0)
    embedding_df = pd.read_csv(self.annotations, usecols=[h for h in headers.columns if h != 'key'], dtype=np.float32)
    return embedding_df.values

  def indices2embeddings(self, data):
    embeddings = np.apply_along_axis(self.embedding_mapper, axis=1, arr=data)
    return embeddings

  def embedding_mapper(self, x):
    return self.embeddings[x]

  # Heavily based on:
  # https://github.com/gregmcinnes/Hubble2D6/blob/master/bin/hubble.py#L77
  def build_seqs(self):
    samples = self.get_samples(self.vcf)
    sample_seqs = self.init_seq_data(samples)
    embeddings = self.precomputed_embeddings()
    ignored_variants, variant_count = 0, 0

    with open(self.vcf) as f:
      for line in f:
        if line.startswith("#"):
          continue

        vcf_row = self.parse_vcf_line(line)
        variant_count += 1

        #! - Alt embedding is array of embeddings
        current_embeddings = self.get_embedding(embeddings, ref=vcf_row['ref_key'], alts=vcf_row['alt_key'])

        if current_embeddings is None:
          # Skip variants with no encodings :(
          ignored_variants += 1
          continue

        # For each row, build embeddings of each sample entry in row
        for i in range(len(samples)):
          s = samples[i]
          gt = vcf_row['diplotypes'][i]
          h1, h2 = [int(i) for i in gt.split('/')]

          h1_position_idx, h1_embedding_idx = current_embeddings[h1]['position_index'], current_embeddings[h1]['embedding_index']
          h2_position_idx, h2_embedding_idx = current_embeddings[h2]['position_index'], current_embeddings[h2]['embedding_index']

          sample_seqs[s][0][h1_position_idx] = h1_embedding_idx
          sample_seqs[s][1][h2_position_idx] = h2_embedding_idx

    print("Ignored %d variants out of %d" % (ignored_variants, variant_count));
    return sample_seqs

  def get_embedding(self, embeddings, ref, alts):
    keys = [ref] + alts
    try:
      if not all (k in embeddings for k in keys):
        raise LookupError(keys)

      return [embeddings[k] for k in keys]
    except LookupError as err:
      return None

  # https://github.com/gregmcinnes/Hubble2D6/blob/master/bin/hubble.py#L188
  def precomputed_embeddings(self):
    wd = sys.path[0]
    file = os.path.join(wd, self.embeddings_file)

    embeddings = {}
    with open(file) as f:
      for line in f:
        fields = line.rstrip().split()
        key = "%s_%s" % (fields[1], fields[2])
        embeddings[key] = {"position": int(fields[1]),
                           "allele": fields[2],
                           "position_index": int(fields[3]),
                           "embedding_index": int(fields[4])}

    return embeddings

  def parse_vcf_line(self, line):
    fields = line.rstrip().split()
    row = {
      'position': int(fields[1]),
      'ref': fields[3],
      'alt': fields[4], 
      'ref_key': '%s_%s' % (fields[1], fields[3]),
      'alt_key': ['%s_%s' % (fields[1], f) for f in fields[4].split(',')],
      'diplotypes': fields[9:]
    }

    return row

  def get_reference_seq(self):
    ref = None
    # Convert embedding indices of reference gene into list
    with open('./data/ref.seq') as f:
      ref = f.readline().rstrip().split()[1:][0].split(',')

    return ref

  def init_seq_data(self, samples):
    sample_seqs = {}
    reference = self.get_reference_seq()

    # There could be a more efficient method than storing two copies of the reference gene for each ht
    # Maybe only store changes?
    for s in samples:
      sample_seqs[s] = [reference, reference]

    return sample_seqs

  def get_samples(self, vcf):
    samples = []

    with open(vcf) as f:
      for line in f:
        if line.startswith("#CHROM"):
          samples = line.rstrip().split()[9:]
          break
      
    return samples

if __name__=='__main__':
  embedding = FirstStep2Seq(vcf='./simulated_cyp2d6_diplotypes/batch_4.vcf', label_csv='./simulated_cyp2d6_diplotypes/batch_4.labels.csv')
  print("Embeddings...")
  print(embedding.X.shape, embedding.y.shape)
