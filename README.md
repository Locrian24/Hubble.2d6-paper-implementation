# Implementing Hubble.2d6

This repo contains the report, poster, and source code for an attempted implementation of the Hubble.2D6
tool based on the [original paper](https://www.biorxiv.org/content/10.1101/684357v1).

> Link to project repo: https://github.com/Locrian24/seng474-term-project

***IMPORTANT***: This is an extremely naive and incomplete implementation of the Hubble.2d6 tool.
This was an undergraduate project and predictions and behaviour do not directly correspond to Hubble.2d6 for the most part and is far from a reliable tool.
Check out the Hubble.2d6 [repo](https://github.com/gregmcinnes/Hubble2D6)

## Get started

Using Anaconda:
```
conda env create --file cannett_474_env.yml
```

Using `pip` with a virtual environment:
```
python3 -m venv env
source env/bin/activate
pip install -r requirements.txt
```

## Run

Sample file with 3 star alleles:
```
python3 model/hubble.py -v data/sample.vcf
```

Star alleles from PharmVar:
```
python3 model/hubble.py -v step3/data/star_samples.vcf
```
