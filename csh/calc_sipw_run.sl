#!/bin/csh

# Project:
#SBATCH --account=nn9469k

# Wall clock limit:
#SBATCH --time=18:00:00

# Max memory usage:
##SBATCH --mem-per-cpu=64G --partition=hugemem
#SBATCH --mem-per-cpu=4G

#SBATCH --ntasks=4
##SBATCH --ntasks-per-node=1
#SBATCH --nodes=1

##SBATCH --qos=lowpri

set i = $argv[1]

Rscript clalc_IPW_pADHD_2SB_v9_sl.R $i > ipw_i{$i}.Rout

