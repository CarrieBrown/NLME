# Non-Linear Mixed Models
## Comparision of Pseudolikelihood and Approximation Methods in SAS

Carrie Brown - January 2020

Analysis for each model can be found within the respective folders:
 - Logistic NLMM: `./logistic/`

![equation](https://latex.codecogs.com/svg.download?%5Clarge%20%5Cmu_%7Bij%7D%3D%5Calpha%20+%20a_i%20+%20%5Cfrac%7B%5Cbeta%20+%20b_i%7D%7B%5Cexp%7B-%28%5Cgamma%20+%20c_i%20+%20%28%5Cdelta%20+%20d_i%29X_%7Bij%7D%29%7D%7D)

 - Michaelis-Menton: `./michaelis-mention/`

These analyses are written to be ran on an HPC Cluster with the SLURM scheduler.

To begin an analysis, run the `start` executable within the desired model's directory.

For example, generating 1000 simulations for the Logistic NLMM can be done with the command:
`./logistic/start <folder_name>`
where `<folder_name>` is replaced with the desired name for the output directory.


LaTeX equation generator used: https://www.codecogs.com/latex/eqneditor.php
