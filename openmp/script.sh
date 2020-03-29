#|/bin/bash


tests="$1"
for f in ${tests}/*.in; do
	echo "============================"
	out=${f::-3}
	for i in 1 2 4 8; do
		eval "$(./aux.sh $i )"
		echo "NUM_THREADS =" $OMP_NUM_THREADS
		outfile="${out}_time_$OMP_NUM_THREADS.out"
		exec 3>&1 4>&2
		output=$( { time ./matFact-omp $f 1>&3 2>&4; } 2>&1)
		exec 3>&- 4>&-
		echo $output > $outfile
	done
done