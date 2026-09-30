

# Prior Smoothing for Multivariate Disease Mapping Models

This repository contains the R code to fit in NIMBLE the code to replicate and reproduce the simulation study and real-data illustration of the paper entitled *"Prior Smoothing for Multivariate Disease Mapping Models"* [(Retegui et al., 2026)](https://doi.org/10.1016/j.spasta.2026.101055).

## Table of contents

1.  [Data](#data)
2.  [R code](#r-code)
    1.  Within Simulation Study
    2.  Across Simulation Study
    3.  RealData Illustration
3.  [Acknowledgements](#acknowledgements)
4.  [References](#references)

## Data 

This folder contains the datasets and cartography files used in the simulation studies and data illustrations presented in the work.

-   `Data_Spain_47areas.Rdata`, `Data_Spain_100areas.Rdata` and `Data_Spain_300areas.Rdata`: These datasets include counts for lung cancer and the corresponding population at risk, focusing on aggregate data for females in Spain from 2019 to 2021 divided into 47, 100 and 300 areas, respectively.

-   `Carto_Spain_47areas.Rdata`, `Carto_Spain_100areas.Rdata` and `Carto_Spain_300areas.Rdata`: These cartography files include the cartography of peninsular Spain divided into 47, 100 and 300 areas, respectively.

-   `Data_SimulationStudy_S1_47areas.Rdata`, `Data_SimulationStudy_S1_100areas.Rdata`, and `Data_SimulationStudy_S1_300areas.Rdata`: These datasets include the simulated data for Scenario 1 when Spain is divided into 47, 100 and 300 areas, respectively. Similar file naming structure, replacing `S1` with `S2`, `S3` or `S4` can be found for Scenarios 2, 3 and 4, respectively.

## R code 

This folder contains the R code to replicate and reproduce the within prior and across priors simulation studies, as well as the data illustration described in the paper. The code is organized into three subfolders, each corresponding to a specific part of the study:

1.  Within Prior Simulation Study (available [here](https://github.com/spatialstatisticsupna/Multivariate_Smoothing/tree/main/R/Within_SimulationStudy)). 

    -   This folder includes code to fit the three spatial priors discussed in Section 3.2 of the paper: iCAR, LCAR, and L$_j$CAR.

    -   Model-fitting scripts follow the format: `Code_Within_iCAR.R`, `Code_Within_LCAR.R` and `Code_Within_LjCAR.R`.

    -   Note: the LCAR model implementation is based on Beltrán-Sánchez et al (2024).

    -   Results are generated using `Code_ResultsWithin.R`.

2.  Across Priors Simulation Study (available [here](https://github.com/spatialstatisticsupna/Multivariate_Smoothing/tree/main/R/Across_SimulationStudy)). 

    -   This folder includes both model-fitting scripts and result generation code.

    -   Model-fitting scripts follow the format: `Code_SimulationStudy_iCAR.R`, `Code_SimulationStudy_LCAR.R` and `Code_SimulationStudy_LjCAR.R`.

    -   Results are generated using `Code_EmpiricalMetrics.R`.

3.  Data Illustration (available [here](https://github.com/spatialstatisticsupna/Multivariate_Smoothing/tree/main/R/RealData_Illustration)).

    -   This folder contains the code used to fit the models and produce the results for the real data applications.

    -   Separate scripts are provided for the pairwise analysis and joint analysis.

## Acknowledgements

The work was supported by Project PID2020-113125RB-I00/MCIN/AEI/10.13039/501100011033, PID2024-155382OB-I00 funded by MICIU/AEI/10.13039/501100011033 and FEDER, UE and BIOSTATNET - PROYECTOS REDES DE INVESTIGACIÓN 2024 - RED2024-153680-T/MICIU/AEI/.

![plot](https://github.com/spatialstatisticsupna/Multivariate_Smoothing/blob/main/miciu-aei.jpg)

## References

[Retegui, G.,Ugarte, M.D., Etxeberria, J. and Gelfand, A.E. (2026). Prior smoothing for multivariate disease mapping models. *Spatial Statistics*, in press, https://doi.org/10.1016/j.spasta.2026.101055](https://doi.org/10.1016/j.spasta.2026.101055)
