
set files = `ls -lt samples/samples_m1_a3_* | awk '$0 !~ /+/ {print $9}'`
set files = ()

foreach f ($files) 

 set finished = `cat $f | grep Elapsed`

 if ($#finished > 0) then
  set ff = `echo $f | sed 's/samples\/samples_m1_//' | sed ' s/_init0_v9.csv//' | sed 's/_/-/g'`
  set jobid =  `queue | grep $ff | awk '{print $1}'`
  if ($#jobid > 0) then
    set rmf = `echo $f | sed 's/v9.csv/v9+.csv/'`
    echo $jobid $ff $f $rmf done
    #scancel $jobid
    #rm $rmf
  endif
  
  if ($#finished < 1) then
    set ff = `echo $f | sed 's/samples\/samples_m1_//' | sed ' s/_init0_v9+.csv//' | sed 's/_/-/g'`
    set jobid =  `queue | grep $ff | awk '{print $1}'`
    echo $jobid $ff $f
  endif
  
 endif

end 



set files = `ls slurm*`

foreach f ($files) 
set finished = `cat $f | grep Elapsed`
  if ($#finished > 0) then
    echo $f $#finished
    rm $f
  endif

end 
