#!/bin/csh



# Project:
#SBATCH --account=nn9469k

# Wall clock limit:
#SBATCH --time=48:00:00

# Max memory usage:
##SBATCH --mem-per-cpu=64G --partition=hugemem
#SBATCH --mem-per-cpu=4G

#SBATCH --ntasks=1
#SBATCH --ntasks-per-node=1
#SBATCH --nodes=1

##SBATCH --qos=lowpri

set model = $argv[1]
set analysis = $argv[2]
set imp = $argv[3]
set chain = $argv[4]

module load R

Rscript run_zinb.R $model $analysis $imp $chain > a{$analysis}_i{$imp}_c{$chain}_m{$model}s.Rout

