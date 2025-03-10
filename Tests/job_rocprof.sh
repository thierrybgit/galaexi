#!/bin/bash -l
#SBATCH --job-name=galaexi
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --time=01:30:00
#SBATCH --partition=dev-g
#SBATCH --hint=nomultithread
#SBATCH --account=project_462000031
#SBATCH --gpus-per-node=1

module load LUMI PrgEnv-cray partition/G buildtools cray-hdf5-parallel rocm
module li

export LD_LIBRARY_PATH=$CRAY_LD_LIBRARY_PATH:$LD_LIBRARY_PATH
GALAEXI_ROOT=../galaexi
export LD_LIBRARY_PATH=$GALAEXI_ROOT/build/lib:$LD_LIBRARY_PATH
EXE=$GALAEXI_ROOT/build/bin/galaexi

rm -rf TGV_test_* userblock*
srun rocprof --stats --basenames on -o stats-only.csv $EXE parameter_flexi.ini
rm -rf TGV_test_* userblock*
srun rocprof --stats --hsa-trace --basenames on -o stat+hsatrace.csv $EXE parameter_flexi.ini
rm -rf TGV_test_* userblock*
srun rocprof --stats --hip-trace --basenames on -o stat+hiptrace.csv $EXE parameter_flexi.ini
rm -rf TGV_test_* userblock*
srun rocprof --stats --sys-trace --basenames on -o stat+systrace.csv $EXE parameter_flexi.ini
