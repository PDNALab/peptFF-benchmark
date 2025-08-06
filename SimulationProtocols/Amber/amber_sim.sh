#!/bin/bash
#SBATCH --job-name=minim
#SBATCH --output=min.out
#SBATCH --error=min.err
#SBATCH --mail-type=ALL
#SBATCH --time=6:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --ntasks-per-node=1
#SBATCH --distribution=cyclic:cyclic
#SBATCH --mem-per-cpu=500
#SBATCH --qos=alberto.perezant

cd $SLURM_SUBMIT_DIR

module purge

module load  cuda/11.1.0 nvhpc/20.11 openmpi/4.0.5 amber/20

## MINIMIZATION STAGE
srun sander -O -i min1.in -o min1.out -p peptide_sol.prmtop -c peptide_sol.rst7 -r min1.rst -inf min1.mdinfo -ref peptide_sol.rst7
srun sander -O -i min2.in -o min2.out -p peptide_sol.prmtop -c min1.rst -r min2.rst -inf min2.mdinfo -ref min1.rst
srun sander -O -i min3.in -o min3.out -p peptide_sol.prmtop -c min2.rst -r min3.rst -inf min3.mdinfo -ref min2.rst
srun sander -O -i min4.in -o min4.out -p peptide_sol.prmtop -c min3.rst -r min4.rst -inf min4.mdinfo -ref min3.rst
srun sander -O -i min5.in -o min5.out -p peptide_sol.prmtop -c min4.rst -r min5.rst -inf min5.mdinfo -ref min4.rst
srun sander -O -i min6.in -o min6.out -p peptide_sol.prmtop -c min5.rst -r min6.rst -inf min6.mdinfo -ref min5.rst

## EQUILIBRATION STAGE
pmemd.cuda -O -i mdt.in -o heat.out -p peptide_sol.prmtop -c ../minim/min6.rst -r heat.rst -x heat.nc -inf heat.mdinfo -ref ../minim/min6.rst
pmemd.cuda -O -i npt1.in -o npt1.out -p peptide_sol.prmtop -c heat.rst -r npt1.rst -x npt1.nc -inf npt1.mdinfo -ref heat.rst
pmemd.cuda -O -i npt2.in -o npt2.out -p peptide_sol.prmtop -c npt1.rst -r npt2.rst -x npt2.nc -inf npt2.mdinfo -ref npt1.rst
pmemd.cuda -O -i md.in -o md1.out -p peptide_hmass.prmtop -c npt2.rst -r md1.rst -x md1.nc -inf md1.mdinfo -ref npt2.rst #200ns simulation

## PRODUCTION STAGE
# Set the initial step number
initial_step=2

# Set the number of steps
num_steps=49

# Set the input file names
input_file="md.in"
output_prefix="md"

# Perform the MD simulations
for ((step=initial_step; step<=(initial_step+num_steps-1); step++))
do
    # Define the output file names for the current step
    output_file="${output_prefix}${step}.out"
    restart_file="${output_prefix}${step}.rst"
    trajectory_file="${output_prefix}${step}.nc"
    info_file="${output_prefix}${step}.mdinfo"
    
    # Define the input file from the previous step
    prev_step=$((step-1))
    prev_restart_file="${output_prefix}${prev_step}.rst"
    prev_output_file="${output_prefix}${prev_step}.out"
    
    # Run the MD simulation
    pmemd.cuda -O -i $input_file -o $output_file -p peptide_hmass.prmtop -c $prev_restart_file -r $restart_file -x $trajectory_file -inf $info_file -ref $prev_restart_file
    
    # Check if the previous output file exists and remove it
    if [ -f $prev_output_file ]; then
        rm $prev_output_file
    fi
done



