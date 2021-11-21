#!/bin/sh

for J in 20 50 100 200 400 600 800 1000 1200 1400 1600 1800 2000
do
    printf "Submitting job for J = $J\n"
    $(qsub -F "0.5 1 $J 10000 5000 100 --binX1" /lustre/haven/user/crock2/AdM/pimc/Fermionic-PIMC-toy/submission_julia.pbs);
done
