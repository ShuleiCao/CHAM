#!/bin/bash
#SBATCH -J match_redmapper_with_gold
#SBATCH -o "slurm_outputs/match_redmapper_with_gold_%j.out"
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --mail-type=END,FAIL,REQUEUE
#SBATCH --mail-user=shuleic@smu.edu

##SBATCH --cpus-per-task=128
##SBATCH -p "standard-s"
##SBATCH -t 00-23:59:59
##SBATCH --mem=500GB

##SBATCH --cpus-per-task=128
##SBATCH --mem=2000000MB
#SBATCH --cpus-per-task=64
#SBATCH --mem=1000000MB
#SBATCH -p "highmem"
#SBATCH -t 04-23:59:59

# module load slurm/22.05.8-zglugih
# module load gsl/2.7.1-oyt4rtp
# module purge
module load gcc/11.2.0
module load openmpi/4.1.6-vfi4iwj
source /users/shuleic/miniconda3/etc/profile.d/conda.sh
conda activate /lustre/work/client/users/shuleic/.conda/envs/redmapper_env

python /users/shuleic/scripts/pythonscripts/All_other_useful_tools/match_redmapper_Cardinal/match_redmapper_with_gold_process.py --mismatched
# python /users/shuleic/scripts/pythonscripts/All_other_useful_tools/match_redmapper_Cardinal/match_redmapper_with_gold_process_noshare.py
