filenames='c2670 c3540 c5315 c6288 c7552 s13207 s15850 s35932'

for file in $filenames
do
	if [[ $file == c* ]];
	then
		type='comb'
	else
		type='seq'
	fi
	./single_bench_run.sh $file $type > Results/$file.log &
done

wait
