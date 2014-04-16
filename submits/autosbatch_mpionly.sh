#!/bin/bash

# Submit jobs. Assume 16 cores per proc. Make sure line numbers are the same.
for (( i = 1; i<= 16 ; i++ ))
do
    sed -i "/--nodes=/c #SBATCH --nodes=$i                      # number of nodes requested (n)" submits/slurm_noomp.sh
    sed -i "/--ntasks=/c #SBATCH --ntasks=$((16*$i))                    # required number of CPUs (n)" submits/slurm_noomp.sh
    sbatch submits/slurm_noomp.sh
done

