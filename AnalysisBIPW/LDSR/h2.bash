vars=($(cut var_N.txt -d " " -f1 ) )

LDdir=/usit/abel/u1/guidopb/programming/R/representativeness/IPW/GWASresults/LDscores_broadinsititute/1000G_Phase3_baselineLD_ldscores/

for i in "${vars[@]}"
do
	ldsc.py --h2 $i.1000G.EUR.sumstats.gz --ref-ld-chr $LDdir --w-ld-chr $LDdir --out h2/$i.1000G.EUR.h2	
done


