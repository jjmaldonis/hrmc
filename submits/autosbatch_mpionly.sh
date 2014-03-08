#!/bin/bash

LINE1=$(sed -n "/--nodes=/=" slurm.sh)
LINE2=$(sed -n "/--ntasks=/=" slurm.sh)

# Submit jobs. Assume 16 cores per proc. Make sure line numbers are the same.
for (( i = 1; i<= 16 ; i++ ))
do
    sed -i "/--nodes=/c #SBATCH --nodes=$i                      # number of nodes requested (n)" slurm.sh
    sed -i "/--ntasks=/c #SBATCH --ntasks=$i                    # required number of CPUs (n)" slurm.sh
    sbatch slurm.sh
done

