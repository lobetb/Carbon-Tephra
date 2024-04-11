[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.10958530.svg)](https://doi.org/10.5281/zenodo.10958530)

# What is **Carbon-Tephra** ?
## Purpose
**Carbon-Tephra** is a Python script to compute the amount of carbon trapped under [tephra](https://en.wikipedia.org/wiki/Tephra) deposited after volcano eruptions. It takes two input files : a MAT file describing isoprobabilities of tephra deposit (from [Biass, S., Bonadonna, C., Connor, L. et al. TephraProb: a Matlab package for probabilistic hazard assessments of tephra fallout. J Appl. Volcanol. 5, 10 (2016).](https://doi.org/10.1186/s13617-016-0050-5) ). and an xls file from the Global Volcanism Program (https://volcano.si.edu/).
Before use, the xls file needs to be opened in Excel, the first line removed and the file resaved (something makes the file unreadable directly by Python).

## How it works
**Carbon-Tephra** uses the data from the GVP to compute the volcanic eruptions of VEI 4, 5 and 6 during the holocène. VEI 4 eruptions are too small to leave geological traces, and as such we only know about those that occurred since the beginning of systematic logging of volcanic eruptions. The script aims at computing the amount of carbon accumulated under tephra during the holocène, way before that systematic logging started.

The script includes three distincts modes to infer dates and locations of eruptions during the holocène :

1. The "stochastic" mode calculates the probability of each VEI for each volcano and makes a random draw at each timestep to determine if an eruption occurs for each volcano-vei duo.
2. The "mixed" mode uses geological data for VEIs 5 and 6 and the stochastic approach for VEI 4 eruptions only.
3. The "sequential" mode changes the workflow : at each timestep, computes if there's an eruption, then its VEI, then the volcano according to their respective probabilities.

After all the eruptions are attributed, the coordinate points at which the tephra cover is sufficient to trap accumulated carbon are extracted from the MAT file. For each point, the amount of trapped C is computed between the eruption dates. An option also allows to compute the C accumulated on the surface since the last eruption.

While first versions of the script used an UTM coordinates system, this version use lat/lon (conversions are done using the [UTM package](https://github.com/Turbo87/utm). This allows for the script to be used on any country, or the whole world.

Batch execution is possible and code is included in this project for that purpose.

## Outputs
The script outputs a number of files :
1. in a "Runs" subfolder, a file "runX.csv" containing lat/lon pairs and the amount of stored C and surface C (if the option is enabled)
2. in a "logC" subfolder, a file "logCrunX.csv" containing the total amount of C stored as a function of time of simulation.
3. in a "logErupt" subfolder, a file containing a line per eruption, with the volcano, the VEI and the year of eruption.

