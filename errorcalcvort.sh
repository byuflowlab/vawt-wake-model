#!/bin/bash

#Submit this script with: sbatch thefilename


#SBATCH --time=24:00:00		# walltime
#SBATCH --ntasks=64		# number of processor cores (i.e. tasks)
#SBATCH --mem-per-cpu=2G	# memory per CPU
#SBATCH --nodes=4		# number of nodes
#SBATCH -J "errorcheck"   # job name
#SBATCH -C 'ib'


# Compatibility variables for PBS. Delete if not needed.
# export PBS_NODEFILE=`/fslapps/fslutils/generate_pbs_nodefile`
# export PBS_JOBID=SLURM_JOB_ID
# export PBS_O_WORKDIR="SLURM_SUBMIT_DIR"
# export PBS_QUEUE=batch

#cd ~/fslhome/ebtingey/compute
printf "Job debugging info\n___\n"
printf "Node List"
cat $PBS_NODEFILE
printf "___\n"
printf "SIM_FILE variable: $1\n"
printf "___\n"
printf "Current Directory: $PWD\n"
printf "___\n"
printf "Starting STAR-CCM+ run at "
date
printf "\n___\n"

module purge
module load python/2/7
module load compiler_gnu/4.9.2

python errorcalc_full.py

printf "___\n"
printf "Ending starccm+run at "
date
printf "\___\n"
printf "Removing temp files\n"
rm -vf ./*~
printf "___\n"
exit 0
