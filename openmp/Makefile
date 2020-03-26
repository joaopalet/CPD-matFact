# Makefile

main:
	gcc -o matFact-omp matrix.c matFact-omp.c -fopenmp

clean:
	rm -f *.o matFact-omp

run:
	./matFact-omp ../instances/inst500-500-20-2-100.in
