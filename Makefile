# Makefile

main:
	gcc -o matFact matFact.c

clean:
	rm -f *.o matFact

run:
	./matFact < instances/inst0.in
