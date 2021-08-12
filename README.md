# Implementing Hubble.2d6
> Link to project repo: https://github.com/Locrian24/seng474-term-project

This repo contains the report, poster, and source code for an attempted implementation of the Hubble.2D6
tool based on the [original paper](https://www.biorxiv.org/content/10.1101/684357v1). This was a term project for the SENG 474 class at the University of Victoria.

***IMPORTANT***: This is an extremely naive and incomplete implementation of the Hubble.2d6 tool.
This was an undergraduate project and predictions and behaviour do not directly correspond to Hubble.2d6 for the most part and is far from a reliable tool.
Check out the Hubble.2d6 [repo](https://github.com/gregmcinnes/Hubble2D6)

## Acknowledgements

This majority of logic in pre-processing and post-processing of data is taken from the original tool ([here](https://github.com/gregmcinnes/Hubble2D6)).
This also includes supplementary data such as pre-computed embeddings, and specifics in the deep learning networks' architecture.

## How to run

### Creating environment

Using Anaconda:
```python
conda env create --file cannett_474_env.yml
```

Using `pip` with a virtual environment:
```python
python3 -m venv env
source env/bin/activate
pip install -r requirements.txt
```

### Run

Sample file with 3 star alleles:
```python
python3 model/hubble.py -v data/sample.vcf
```

Star alleles from PharmVar:
```python
python3 model/hubble.py -v step3/data/star_samples.vcf
```

## Supplementary material
As well as the source code for the implementation, Google Colab notebooks are included showing the training processes as well as generation of evaluation metrics.

Colab notebooks are found in the `notebooks` directory.

