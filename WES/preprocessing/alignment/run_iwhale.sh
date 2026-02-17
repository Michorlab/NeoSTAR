#!/bin/bash
#SBATCH --job-name=reruns
#SBATCH -n 16
#SBATCH --mem=30G       # total memory need
#SBATCH --time=7-00:00:00
#SBATCH --mail-type=END,FAIL # email notification when job ends/fails
#SBATCH --mail-user=jacobg@ds.dfci.harvard.edu # email to notify
##SBATCH --error=/michorlab/jacobg/dryrun.error      # error file
##SBATCH --output=/michorlab/jacobg/dryrun.out      # output file

# NOTE TO OTHERS: VARIANT CALLING FROM THIS PIPELINE WAS NOT USED FOR ANALYSYS

module load singularity

# Change to appropriate directory !
dir=/liulab/jacobg/Ellisen/iwhale_rundirs/reruns
cd $dir
mkdir ${dir}/Variants # needed because we cant make these in singularity because we dont mount /working itself
mkdir ${dir}/VCF

declare -a BIND_ARGS=()
MOUNT_DIR="/working"

# Find all symlinks and construct bind arguments
while IFS= read -r -d '' symlink; do
    relative_path=${symlink#$dir/}
    internal_path="$MOUNT_DIR/$relative_path"
    target=$(readlink -f "$symlink")
    BIND_ARGS+=("--bind $target:$internal_path")
done < <(find "$dir" -type l -print0)

echo ${BIND_ARGS[@]}
echo singularity exec --bind ${dir}/configuration.py:/working/configuration.py --bind ${dir}/tumor_control_samples.txt:/working/tumor_control_samples.txt ${BIND_ARGS[@]} --bind /michorlab/jacobg/Michalina/iwhale/annotations:/annotations --bind /michorlab/jacobg/Michalina/iwhale/iwhale_sandbox/tmp:/tmp  --bind /liulab/jacobg/Ellisen/iwhale/iWhale/iwhale.py:/iwhale.py /michorlab/jacobg/Michalina/iwhale/iwhale_custom.sif iwhale -s /tmp


singularity exec --pwd /working --bind ${dir}/configuration.py:/working/configuration.py --bind ${dir}/tumor_control_samples.txt:/working/tumor_control_samples.txt ${BIND_ARGS[@]} --bind /michorlab/jacobg/Michalina/iwhale/annotations:/annotations --bind /michorlab/jacobg/Michalina/iwhale/iwhale_sandbox/tmp:/tmp --bind /liulab/jacobg/Ellisen/iwhale/iWhale/iwhale.py:/iwhale.py --bind ${dir}/Variants:/working/Variants --bind ${dir}/VCF:/working/VCF /michorlab/jacobg/Michalina/iwhale/iwhale_custom.sif /iwhale.py -s /tmp



