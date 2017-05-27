This `research` directory contains all the scripts to create the results used in my PhD thesis. Follow the following steps to recreate the benchmark data and/or results.

### Requirements

* installed `tom` toolkit
* Python 3, numpy, scipy, ipython, ipyparallel, jupyter (e.g. from the Anaconda distribution)

### Setup

* Start an ipyparallel cluster from this directory with three engines:
    `ipcluster start --n=3`
* Start a jupyter notebook server from this directory:
    `jupyter notebook`

### 1. Benchmark data

* Enter the `benchmarks` directory.
* Run `python bible_and_ecoli.py` to download and prepare the **real-world benchmarks**.
* The **POMDP benchmarks** are provided ready to use in a json format. To re-create:
  + Install `zmdp` from `https://github.com/trey0/zmdp.git`
  + Edit `ZMDP_LIB_DIR` and `ZMDP_INCLUDE_DIR` in the `Makefile` to point to the correct locations
  + Run `make pomdps` to obtain the original POMDP files, compile `pomdp_to_json` and convert the POMDPs to json format
* Run the `Synthetic OOM Benchmarks` jupyter notebook to compute the **OOM benchmarks**. This also generates the figure `oom_benchmark_selection.pdf` and gives some information about the OOM benchmarks. Computation time is ~1.5 hours on my laptop.

### 2. Comparison of various learning settings




