This `research` directory contains all the scripts to create the results used in my PhD thesis. Follow the following steps to recreate the benchmark data and/or results.

### Requirements

* installed `tom` toolkit (the version tagged as "phd" was used to create the results in the thesis)
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
  + Run `make pomdps` to obtain the original POMDP files, compile `pomdp_to_json` and convert the POMDPs to json
* The **OOM benchmarks** are provided ready to use. To re-create, run the `Synthetic OOM Benchmarks` jupyter notebook. This also generates the figure `oom_benchmark_selection.pdf` and gives some information about the OOM benchmarks. Computation time is ~1.5 hours on my laptop.

### 2. Computation of all research results

* Run the `Compute` jupyter notebook. Computation time is approx 2.5 days on my laptop.

### 3. Generate plots

* Run the `Plots` jupyter notebook.
