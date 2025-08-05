#!/bin/bash
#SBATCH --job-name=a99SBdisp
#SBATCH --output=P_md.out
#SBATCH --error=P_md.err
#SBATCH --time=7-00:00:00
#SBATCH --ntasks=1
#SBATCH --mem-per-cpu=2000
#SBATCH --cpus-per-gpu=8
#SBATCH --gpus-per-task=1
#SBATCH --partition=gpu
#SBATCH --account=alberto.perezant
#SBATCH --qos=alberto.perezant
pwd; hostname; date

module purge
module load cuda/10.0.130  intel/2018.1.163  openmpi/4.0.0  gromacs/2019.2

export OMP_NUM_THREADS=$SLURM_CPUS_PER_GPU

pwd

############################################################ Subroutines: START  ############################################################
minimization()
{
	##MINIMISATION USING STEEPEST GRADIENT##
	gmx grompp -f min1.mdp -c peptide_solv_ions.gro -p topol.top -o min1.tpr
	srun --mpi=pmi2 --accel-bind=g --ntasks=$SLURM_NTASKS gmx mdrun -v -deffnm min1 -nb gpu

	##MINIMISATION USING CONJUGATE GRADIENT##
	gmx grompp -f min2.mdp -c min1.gro -p topol.top -o min2.tpr
	srun --mpi=pmi2 --accel-bind=g --ntasks=$SLURM_NTASKS gmx mdrun -v -deffnm min2 -nb gpu
}

equilibration()
{
	##EQUILIBRATION (NVT for 10ns)
	gmx grompp -f nvt.mdp -c min2.gro -r min2.gro -p topol.top -o nvt.tpr
	srun --mpi=pmi2 --accel-bind=g --ntasks=$SLURM_NTASKS gmx mdrun -deffnm nvt -nb gpu

	##EQUILIBRATION (NPT for 20ns)
	gmx grompp -f npt.mdp -c nvt.gro -r nvt.gro -t nvt.cpt -p topol.top -o npt.tpr
	srun --mpi=pmi2 --accel-bind=g --ntasks=$SLURM_NTASKS gmx mdrun -deffnm npt -nb gpu
}

production()
{	
	cd production/

	##PRODUCTION RUN (each step 200ns)
	gmx grompp -f md.mdp -c ../npt.gro -r ../npt.gro -t ../npt.cpt -p topol_hmass.top -o md_1.tpr -maxwarn 1
	srun --mpi=pmi2 --accel-bind=g --ntasks=$SLURM_NTASKS gmx mdrun -deffnm md_1 -nb gpu

        for i in {2..50}
        do
                j=`expr $i - 1`
                gmx grompp -f md.mdp -c md_$j.gro -t md_$j.cpt -p topol_hmass.top -o md_$i.tpr
                srun --mpi=pmi2 --accel-bind=g --ntasks=$SLURM_NTASKS gmx mdrun -deffnm md_$i -nb gpu
                sleep 120
    
        done
}

############################################################ Subroutines: END  ############################################################

############################################################ MAIN SCRIPT BODY: START #######################################################

###################    Minimization   #################
minimization
#exit
###################    Equilibration  #################
equilibration
#exit
###################    Production     #################
production

############################################################ MAIN SCRIPT BODY: END ########################################################

exit
