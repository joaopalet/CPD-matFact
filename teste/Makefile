# Makefile

main:
	gcc -o matFact-omp matrix.c matFact-omp.c -fopenmp -g

clean:
	rm -f *.o matFact-omp

run:
	# ./matFact-omp ../instances/inst0.in
	./matFact-omp ../instances/inst500-500-20-2-100.in
	# ./matFact-omp ../instances/inst30-40-10-2-10.in
