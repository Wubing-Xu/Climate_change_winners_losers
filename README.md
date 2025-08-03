# Climate_change_winners_losers

## Data
This metacommunity time series analyzed were selected from four databases: BioTIME, RivFishTIME, InsectChange, and Metacommunity Resurveys. The BioTIME data can be accessed on Zenodo (https://doi.org/10.5281/zenodo.2602708); the RivFishTIME data can be accessed through the iDiv Biodiversity Portal (https://doi.org/10.25829/idiv.1873-10-4000); the InsectChange data can be accessed on KNB (https://doi.org/10.5063/F11V5C9V) or through the data paper (http://onlinelibrary.wiley.com/doi/10.1002/ecy.3354/suppinfo); the 'Metacommunity Resurveys' data was compiled using the R code available here (https://github.com/chase-lab/metacommunity_surveys/tree/german_resurvey) and can be accessed through the iDiv Biodiversity Portal (https://doi.org/10.25829/idiv.3503-jevu6s).

In this repository, the selected metacommunity time series can be accessed:
```
  data/Combined_assemblages.RDATA
```

The first RDATA file contains filtered dataset based on grid-based approach, and the second contains only sites in same geographic coordinates sampled across years (see manuscript for more details).

## R Analysis Files
These scripts were used to prepare the data for analysis, calculate species' thermal limits and positions, calculate temporal changes in occupancy and abundance, fit models and produce figures.

To run first few R scripts (named with prefix1-5), raw assemblage data are required to download from links described above. The output is the selected metacommunity time series and estimated species' thermal limits (the file `Species_thermal_limits.RDATA` in the folder data).

Please note that some of the code in this repository was written to run on a HPC cluster.
