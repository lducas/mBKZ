# Module-BKZ
Artifact for the paper _Predicting Module-Lattice Reduction_ by Léo Ducas, Lynn Engelberts, and Paola de Perthuis. 

This paper is accepted for publication in ASIACRYPT 2025, with the full version available on the IACR Cryptology ePrint Archive: https://eprint.iacr.org/2025/1904. 
Slides presenting this work may be found in the `slides/` directory. 

This artifact provides an implementation of the module-BKZ algorithm and the prediction framework from the paper. 
It allows users to: 
- run the module-BKZ algorithm on random q-ary module lattices over cyclotomic number fields of low degree;
- reproduce the paper's concrete predictions for slopes and blocksizes, including the individual terms contributing to the slopes;
- reproduce the paper's asymptotic predictions for the blocksizes.

The repository also includes the CSV data files used to generate the figures in the paper. 

## Integrated artifact

The paper contains clickable hyperlinks marked [.py], pointing to specific lines in our artifact corresponding to the formulae referenced at that point in the paper and used to generate the figures.

## Installation
The following packages should be installed system-wide: 
`autoconf, automake, libtool, virtualenv, libgmp-dev, libmpfr-dev, libqd-dev`. 

The file `requirements.txt` lists the Python dependencies. 

The following command will install libraries `g6k`, `fpylll`, and `fplll` in a virtual environment:

```
git clone --recurse-submodules https://github.com/lducas/mBKZ.git
cd g6k
PYTHON=python3 ./bootstrap.sh
cd ..
```
For SSH access, replace the first line with: 
```
git clone --recurse-submodules git@github.com:lducas/mBKZ.git
``` 
Note that the execution of `./bootstrap.sh` may take a long time. It creates a virtual environment `g6k-env` inside the `g6k` directory. 

To run the Python scripts, you will need to activate the virtual environment `g6k-env` using:
```
cd g6k
source ./activate
cd ..
```

## Running experiments, examples

Please note that the examples provided here may take hours to run. 

To generate a csv file with the profile of mBKZ over Q(ω_3), Q-dimension 120, after 5 tours of BKZ-12, run:
```
python exp_profile.py c=3 m=120 beta=12 tours=5
```

To measure the slope of mBKZ over Q(ω_3) for increasing blocksizes, Q-dimension 80, run:
```
python exp_slope.py c=3 m=80
```

To measure and predict the skewness over Q(ω_5) and index for increasing rank, run:
```
python exp_skewness_index.py c=5 samples=50
```

To measure and predict the Gaussian Heuristic over Q(ω_5) and index for increasing rank, run:
```
python exp_module_gh.py c=5 samples=50
```

## Reproducing predictions and experiments from the paper

The predicted gains on the blocksize (Fig. 2) are generated with the script `pred_gain.py`.

The extensive experiments from the paper are launched via the `bash` scripts 
`allexp_slope.sh` (Fig. 1), `allexp_profile.sh` (Fig. 4), `allexp_module_gh.sh` (Fig. 5), and `allexp_skewness_index.sh` (Fig. 7 & 8), which will store data in a `data` subfolder. For the profile and slope data, one datafile per experiment is generated, and can be averaged using a script `average.py` in the corresponding `data` subfolder. The slope averaging script also appends predictions to the experimental average for comparison.

These bash scripts will launch up to 35 parallel processes, be sure to use them on a computational server that can handle such load! 

## Core implementation files

The file `cyclotomics.py` contains utility functions to deal with cyclotomic fields, for both predictions and experiments.

The file `modlatred.py` contains our implementation of module-BKZ.

The file `predictions.py` contains the prediction formulae from our paper.

The file `asympt_gain.py` contains our formula for the asymptotic relationship between mBKZ and BKZ blocksizes. 
