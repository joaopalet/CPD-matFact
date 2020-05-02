# Makefile

main:
	mpicc -o matFact-mpi matrix.c matFact-mpi.c -g

clean:
	rm -f *.o matFact-mpi

run:
	mpirun -n 3 ./matFact-mpi ../instances/inst0.in
