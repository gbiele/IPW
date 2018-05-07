#!/bin/csh 

set analysis = $argv[1]
set start_chain = $argv[2]
set end_chain = $argv[3]

if ($#argv == 4) then
  set imputations = $argv[4]
else
  set imputations = 1
endif

set cn = $start_chain
set imp = 1

while ($imp <= $imputations)
  while ($cn <= $end_chain) 
    set outfile = /cluster/home/guidopb/programming/R/representativeness/IPW/csh/slurm-%j-a{$analysis}-i{$imp}-c{$cn}.out
    sbatch -J a{$analysis}-i{$imp}-c{$cn} -o $outfile run_zinbR.sl 1 $analysis $imp $cn    
    echo $outfile 
    @ cn++
    sleep 3
  end
  set cn = $start_chain
  @ imp++
end


