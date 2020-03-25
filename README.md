# Non-Linear Mixed Models
## Comparision of Pseudolikelihood and Quadrature Methods in SAS

Carrie Brown - January 2020

Analysis for each model can be found within the respective folders:
 - Logistic NLMM: `./logistic/`
 - Michaelis-Menton: `./michaelis-mention/`

These analyses are written to be ran on an HPC Cluster with the SLURM scheduler.

To begin an analysis, run the `start` executable within the desired model's directory.

For example, generating 1000 simulations for the Logistic NLMM can be done with the command:
`./logistic/start <folder_name>`
where `<folder_name>` is replaced with the desired name for the output directory.
