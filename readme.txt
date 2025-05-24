This folder contains the code to reproduce and run the results for our VLDB submission.
Required libraries:
1. https://github.com/BlueBrain/HighFive
2. https://github.com/nlohmann/json 
3. Gurobi with license

Note that for 1. and 2. we provide header files in the current folder (/statLibs).
In the final version of the code we will also provide a container to easily build the executable and reproduce the results.

STEP 1: Build the executable
Properly set the path to your Gurobi installation in the file "Makefile"
GUROBI_INCL="" (Should be .../,gurobiVersion/linux64/include)
GUROBI_LIB="" (Should be .../,gurobiVersion/linux64/lib)

Run "make"
This should build the executable ".clustcoeff"

STEP 2: Download the datasets and make sure they are in the (src dst) txt format (see folder /data)
STEP 3: Set up folder structure
create the following folders:
"data/proc"
"data/raw/"
"results/qval/"
"results/ouralg/fixedss/"
"data/raw/labelled/"
"data/raw/txt/"
"results/adapt/"

STEP 4: Place the raw data (from STEP 2) in the folder "data/raw/"
STEP 5: Run 
"python driver.py tohdf5"
"python driver.py exactTriCounts"
"python driver.py modePartition"
Now you are set up and you can replicate all the experiments of the paper
STEP 6:
Use the python file "driver.py" to replicate all the experiments: modes are: modeTestQ, modeFixedSampleSizeOurs, modeCompAdapt

DBLP CASE STUDY:
STEP 1: run the script "downDblp.py" in the folder "/scripts"
STEP 2: run "python driver.py dblpTemporal"
