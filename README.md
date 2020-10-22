# Non-Linear Mixed Effect Models
## Comparision of Pseudolikelihood and Quadrature Methods in SAS

Carrie Brown - January 2020

Analyses for each model can be found within the respective folders:
 - Logistic: `./logistic/`
 - Michaelis-Menten: `./michaelis-menten/`

These analyses are written to be ran in R and SAS on an HPC Cluster with the SLURM scheduler. Post analysis requires the R packages `plyr`, `tidyverse`, and `latex2exp`.

To begin an analysis, run the `start` executable within the desired model's directory.

For example, generating 1000 simulations for the Logistic NLME model can be done with the command:

`cd ./logistic`

`./start <folder_name>`

where `<folder_name>` is replaced with the desired name for the output directory

### References:

 - Pinheiro, J.C. and Bates, D.M. (1995). "Approximations to the Log-Likelihood Function in the Nonlinear Mixed-Effects Model." *Journal of Computational and Graphical Statistics* 4:12-35
 - Stroup, W. W., Milliken, G. A., Claassen, E. A., & Wolfinger, R. D. (2018). *SAS for mixed models: introduction and basic applications*. Cary, NC: SAS.
 - Stroup, Walter W. (2013) *Generalized Linear Mixed Models: Modern Concepts, Methods and Applications*. Boca Faton, FL: CRC Press
 - Wolfringer, R.D. and O'Connell, M.A. (1993). "Generalized Linear Mixed Models: A Pseudo-likelihood Approach." *Journal of Statistical Computation and Simulation* 48:233-243


The code for this project was generated using SAS software, Version 9.4 of the SAS System for Linux. Copyright © 2016 SAS Institute Inc. SAS and all other SAS Institute Inc. product or service names are registered trademarks or trademarks of SAS Institute Inc., Cary, NC, USA. 
