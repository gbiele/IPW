#!/bin/sh

# Project:
##SBATCH --account=uio
#SBATCH --account=NN9266K


# Wall clock limit:
#SBATCH --time=48:00:00

# Max memory usage:
#SBATCH --mem-per-cpu=32G

#SBATCH --ntasks=1
#SBATCH --ntasks-per-node=1
#SBATCH --nodes=1

module load plink2


export PATH="/usit/abel/u1/guidopb/software/ldsc/:$PATH"
export PATH="/usit/abel/u1/guidopb/software/ldsc/ldsc:$PATH"



vars=($(cut var_N.txt -d " " -f1 ) )

LDroot=eur_w_ld_chr/

l=${#vars[@]}
nvars=`expr $l - 1`


#for t1 in `seq 0 $nvars`;
#do
t1=11
	for t2 in `seq 0 $nvars`;
        do 
        	if [ "$t1" -gt "$t2" ]; then
			#echo $t1 $t2
			v1=${vars[$t1]}
			v2=${vars[$t2]}
        		echo $v1 $v2
			ldsc.py --rg $v1.sumstats.gz,$v2.sumstats.gz --ref-ld-chr $LDroot --w-ld-chr $LDroot --out rg/$v1"_"$v2
			#fn=rg/$v1"_"$v2.log
			#if [ ! -f $fn ]; then
			#	echo $v1 $v2 not done
			#	ldsc.py --rg $v1.sumstats.gz,$v2.sumstats.gz --ref-ld-chr $LDroot --w-ld-chr $LDroot --out rg/$v1"_"$v2
			#fi
        	fi
        done    	
#done    


#vars=($(cut var_N.txt -d " " -f1 ) )


#for i in "${vars[@]}"
#do
#	ldsc.py --h2 $i.sumstats.gz --ref-ld-chr $LDroot --w-ld-chr $LDroot --out h2/$i.h2	
#done


