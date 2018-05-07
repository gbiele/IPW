#!/bin/csh

# Project:
#SBATCH --account=nn9469k

# Wall clock limit:
#SBATCH --time=18:00:00

# Max memory usage:
##SBATCH --mem-per-cpu=64G --partition=hugemem
#SBATCH --mem-per-cpu=8G

#SBATCH --ntasks=1
#SBATCH --ntasks-per-node=1
#SBATCH --nodes=1

##SBATCH --qos=lowpri

set stan_model = 2
set analysis = $argv[2]
set imp = $argv[3]
set chain = $argv[4]

set base_dir = /cluster/home/guidopb/programming/R/representativeness/IPW/

set datf = $base_dir/data/data_a{$analysis}_i{$imp}_pADHD_v9.rdump

cd $base_dir/stan


if ($stan_model == 1) then
  set program = zinb_deltaARCC
  set of = $base_dir/samples/samples_m{$stan_model}_a{$analysis}_i{$imp}_c{$chain}_pADHD_v9.csv
else
  set program = zinb_deltaARCClog
  set of = $base_dir/samples/pADHD/samples_m{$stan_model}_a{$analysis}_i{$imp}_c{$chain}_log_v9.csv
endif

./$program sample num_samples=500 num_warmup=250 save_warmup=1 algorithm=hmc init=0 id=$chain data file=$datf output file=$of refresh=10


