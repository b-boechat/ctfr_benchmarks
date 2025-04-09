### Benchmarking ctfr implementations

This supplementary repository contains the code to benchmark the CTFR implementations of the article:

> Placeholder

The main repository of the `ctfr` package is available [here](https://github.com/b-boechat/ctfr).

## Instructions

To repeat the experiments, first clone this repository and install it in a virtual environment with:

```shell
git clone git@github.com:b-boechat/ctfr_benchmarks.git
cd ctfr_benchmarks
python -m venv venv
source venv/bin/activate
pip install -e .

```

This is gonna install a clone of the `ctfr` package (version 0.1.0) with the baseline implementations added.

Then, you can run the benchmarks on your local machine running

```shell
python scripts/benchmarks.py
```

You can specify the number of repetitions for each benchmark with a command line argument. For example, for 10 repetitions, you can use:

```shell
python scripts/benchmarks.py 10
```

## Exact dependencies

The exact dependencies' versions used to run the benchmarks are listed in the `requirements.txt` file. You can install the package with them by running:

```shell
pip install -e .[exact]
```

## Plots

To ensure reproducibility, all spectrogram plots in the article can be generated with the `plots.ipynb` Jupyter Notebook located in the `notebooks` folder.

## Implementations code

The code for the implementations of the CTFR algorithms (both baseline and `ctfr`) is available in a [separate documentation](https://ctfr-benchmarks.readthedocs.io/en/latest/).