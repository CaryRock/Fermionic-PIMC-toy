#! /bin/bash
#SBATCH --partition=workq
#SBATCH --nodes=1
#SBATCH --ntasks=17
#SBATCH --mem=17GB
#SBATCH --job-name=Boson_u_N_comparison
#SBATCH --output=Boson_u_N_comparison.log
#SBATCH --error=Boson_u_N_comparison_error.log

#cd /home/crock2/Projects/Fermionic-PIMC-toy/JULIA/base/ 
pimc.e -T 1.0 -t 0.002 -L 20.0 -X harmonic -I free -m 48.48 -N 2 -u -0.5 -M 256 --canonical --relaxmu --window=1 --gaussian_window_width=0.5 --relax --no_save_state --bin_size=100 -E 100000 -S 1000 --estimator=virial --estimator="linear density rho"
