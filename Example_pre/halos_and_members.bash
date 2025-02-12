#!/bin/bash
#SBATCH -J halos_and_members
#SBATCH -o "slurm_outputs/halos_and_members_%j.out"
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
source /path/to/your/miniconda3/etc/profile.d/conda.sh
conda activate /path/to/your/.conda/envs/redmapper_env

python /path/to/your/script/halos_and_members.py --keys px,py,pz,haloid,m200,r200,mem_match_id,coadd_object_id,id,ra,dec,z \
--member_path /path/to/your_cluster_member/original_member_data_lgt20.fit --lambda_cut_suffix _lgt20 \
--mock_path /path/to/your_mock/Cardinal-3_v2.0_Y6a_gold.h5 --output_loc /path/to/your_output/ \
--redshift_path /path/to/your_mock/Cardinal-3_v2.0_Y6a_bpz.h5 --temp_dir /path/to/your_temp_folder/
