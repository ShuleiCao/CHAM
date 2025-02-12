#!/bin/bash
#SBATCH -J match_cluster_with_halo
#SBATCH -o "slurm_outputs/match_cluster_with_halo_%j.out"
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --mail-type=END,FAIL,REQUEUE
#SBATCH --mail-user=YourEmail
#SBATCH --cpus-per-task=64
#SBATCH --mem=500GB
#SBATCH -p YourPartition
#SBATCH -t 04-23:59:59 # Time limit

# load modules and activate conda environment
module load gcc/11.2.0
module load openmpi/4.1.6-vfi4iwj
source /users/shuleic/miniconda3/etc/profile.d/conda.sh
conda activate /lustre/work/client/users/shuleic/.conda/envs/redmapper_env

python /path/to/your/script/CHM.py --keys px,py,pz,haloid,m200,r200,mem_match_id,coadd_object_id \
--member_path /path/to/your_cluster_member/member_data_lgt20.fits --lambda_cut_suffix _lgt20 \
--halo_path /path/to/your_halo/halo_data.fits --output_loc /path/to/your_output/ \
--temp_dir /path/to/your_temp_folder/

