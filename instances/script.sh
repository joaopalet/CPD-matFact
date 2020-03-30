#|/bin/bash

for f in ./*.in; do
	echo "============================"
	out=${f::-3}
	for i in 1 2 4 8; do
		eval "$(./aux.sh $i )"
		echo "NUM_THREADS =" $OMP_NUM_THREADS
		outfile="./outputs/${out}_time_$OMP_NUM_THREADS.out"
		output=$( { time ../openmp/matFact-omp $f; } 2>&1)
		echo $output &>> $outfile
	done
done