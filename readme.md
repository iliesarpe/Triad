## This folder contains the code to reproduce and run the results for our VLDB paper (Efficient and Adaptive Estimation of Local Triadic Coefficients)


## Required libraries:
1. https://github.com/BlueBrain/HighFive
2. https://github.com/nlohmann/json 
3. [Gurobi](https://www.gurobi.com/) with license

Note that for 1. and 2. we provide header files in the current folder (/statLibs).

_We will soon release a container to build easily the entire code._

__STEP 1a__ (_Build the executable_)

- Properly set the path to your Gurobi installation in the file "Makefile"
- `GUROBI_INCL=""` (Should be `.../gurobiVersion/linux64/include`)
- `GUROBI_LIB=""` (Should be `.../gurobiVersion/linux64/lib`)

__STEP 1b__ (_Compilation_)

Run `make` (that should build `.clustcoeff`), now you are all set-up!

# Reproducing the results

__STEP 2__ 

Download the datasets (e.g., use `wget`) and make sure each dataset is in in the (src dst) txt format (see folder `/data`)

__STEP 3__ 

Set up folder structure create the following folders (e.g., use `mkdir`):

`["data/proc", "data/raw/", "results/qval/", "results/ouralg/fixedss/", "data/raw/labelled/", "data/raw/txt/", "results/adapt/"]`

__STEP 4__
Place the raw data (from STEP 2) in the folder `data/raw/`

__STEP 5__

`python driver.py tohdf5`
`python driver.py exactTriCounts`
`python driver.py modePartition`

Now you are set up and you can replicate all the experiments of the paper!

__STEP 6__

Use the python file `driver.py` to replicate all the experiments: modes are: `modeTestQ`, `modeFixedSampleSizeOurs`, `modeCompAdapt` (see main paper)

__DBLP CASE STUDY__

1. run the script `downDblp.py` in the folder `/scripts`
2. run `python driver.py dblpTemporal`
