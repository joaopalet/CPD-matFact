# Makefile

main:
	mpicc -o matFact-mpi matrix.c matFact-mpi.c -g

clean:
	rm -f *.o matFact-mpi

run:
	mpirun -n 4 ./matFact-mpi ../instances/inst30-40-10-2-10.in
