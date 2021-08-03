import pandas as pd
import numpy as np
import argparse

def main(options):
    vcf = options.file
    name = options.name

    # Fetch curated functions per star allele
    print("Fetching curated functions per star allele...", end="")

    # Skip duplicate row in excel
    df = pd.read_excel('https://doi.org/10.1371/journal.pcbi.1008399.s003', usecols=[0, 1])
    print("Done")

    # Get star alleles in file
    samples = []
    with open(vcf) as f:
        for line in f:
            if line.startswith("#CHROM"):
                samples = line.rstrip().split()[9:]
                break

    print("Found %d samples in vcf file" % len(samples))

    #print("Matching samples to curated functions...", end="")
    labels = []
    for sample in samples:
        star = "*" + str(sample.split("_")[1])
        label = df[df["CYP2D6 Star Allele"] == star]["Curated Function"].values[0]
        labels.append([sample, label])

    # Save label csv
    label_df = pd.DataFrame(labels)
    label_df.to_csv(name, index=False, header=None)
    print("Saved labels to %s" % name)

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-n", "--name", help="name of label file")
    parser.add_argument("-f", "--file", help="vcf with CYP2D6 star alleles of interest")
    options = parser.parse_args()
    main(options)
