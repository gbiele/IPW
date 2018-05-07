vars=($(cut var_N.txt -d " " -f1 ) )

for i in "${vars[@]}"
do
	N=`awk -v v="$i" '{if ($1==v) print $2}' var_N.txt`
	echo $i $N
	munge_sumstats.py --sumstats $i.1000G.EUR.txt --N $N --out $i.1000G.EUR --merge-alleles sniplist.1000G.EUR.QC.txt
done


