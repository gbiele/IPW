#!/bin/csh 

set imps = $argv[1]
set imp = 6

while ($imp <= $imps)
  echo $imp
  sbatch -J calc_sipw_i{$imp} calc_sipw_run.sl $imp   
  sleep 3
  @ imp++
end


