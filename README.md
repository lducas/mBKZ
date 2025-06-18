# Module-BKZ
Artifact for the paper _Predicting Module-Lattice Reduction_
by Léo Ducas, Lynn Engelberts and Paola de Perthuis. 

Pre-pre-print included on this repository (`mBKZ.pdf`). Pre-print soon available on IACR eprint archive. Slides presenting this work may be found in the slides directory. 

## Integrated artifacts

The paper contains clickable hyperlinks marked [.py], which link directly to specific files and lines of our artifacts corresponding to the current formula, feature or technicality discussed at that point of the paper.

## Installation
The following packages should be installed systemwide
`autoconf, automake, libtool, virtualenv, libgmp-dev, libmpfr-dev, libqd-dev`.

The following will install libraries `g6k`, `fpylll` and `fplll` in a virtual environment:

```
git submodule add "https://github.com/fplll/g6k"
cd g6k
PYTHON=python3 ./bootstrap.sh
cd ..
```
To run the python scripts, you will need to activate the virtual environment using:
```
cd g6k
source ./activate
cd ..
```

## Running experiments, examples

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

The file `predictions.py` contains the prediction formulas from our paper.
