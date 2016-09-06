#!/bin/bash

#SBATCH --time=24:00:00         # walltime
#SBATCH --ntasks=1              # number of processor cores (i.e. tasks)
#SBATCH --mem-per-cpu=24G        # memory per CPU
#SBATCH --nodes=1               # number of nodes
#SBATCH -J "cvwakesheet"   # job name


# Compatibility variables for PBS. Delete if not needed.
export PBS_NODEFILE=`/fslapps/fslutils/generate_pbs_nodefile`
export PBS_JOBID=SLURM_JOB_ID
export PBS_O_WORKDIR="SLURM_SUBMIT_DIR"
export PBS_QUEUE=batch

printf "Job debugging info\n___\n"
printf "Node List"
cat $PBS_NODEFILE
printf "___\n"
printf "SIM_FILE variable: $1\n"
printf "___\n"
printf "Current Directory: $PWD\n"
printf "___\n"
printf "Starting optimization at "
date
printf "\n___\n"

module purge
module load python/2/7

args=("$@")
rdata=${args[0]}
cvon=${args[1]}
fitopt=${args[2]}
spr1=${args[3]}
pow1=${args[4]}
pow2=${args[5]}
pow3=${args[6]}
spr2=${args[7]}
skw=${args[8]}
scl1=${args[9]}
scl2=${args[10]}
scl3=${args[11]}

python QME_cvx_cv_curvefit.py $rdata $cvon $fitopt $spr1 $pow1 $pow2 $pow3 $spr2 $skw $scl1 $scl2 $scl3

printf "___\n"
printf "Ending starccm+run at "
date
printf "\___\n"
printf "Removing temp files\n"
rm -vf ./*~
printf "___\n"
exit 0
