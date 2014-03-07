#!/bin/bash
# For Odie
# Submit jobs. Assume 8 cores per proc. Make sure line numbers are the same.
# Change 14 to 16 when I fix the memory problem. DONE
#for (( i = 1; i<= 16 ; i++ ))
for (( i = 1; i<= 32 ; i++ ))
do
    sed -i "/#$ -pe orte/c #$ -pe orte $(($i*4))" submit.sh
    qsub submit.sh
done

