#!/bin/csh



# Project:
#SBATCH --account=nn9469k

# Wall clock limit:
##SBATCH --time=48:00:00

# Max memory usage:
##SBATCH --mem-per-cpu=64G --partition=hugemem
#SBATCH --mem-per-cpu=8G

#SBATCH --ntasks=1
#SBATCH --ntasks-per-node=1
#SBATCH --nodes=1

##SBATCH --qos=lowpri

set model = $argv[1]
set analysis = $argv[2]
set imp = $argv[3]
set chain = $argv[4]

module purge
module load R/3.5.0.gnu
module load gcc/6.3.0

Rscript run_zinbx+raw.R $model $analysis $imp $chain > a{$analysis}_i{$imp}_c{$chain}_m{$model}x.Rout

