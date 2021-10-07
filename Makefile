# Makefile for MPI Game of Life
#
# Jim Teresco, CS 338, Williams College
# CS 341, Mount Holyoke College
# CS 400/335, Siena College
#
CFILES=mpilife.c
OFILES=$(CFILES:.c=.o)
CC=mpicc

mpilife:	$(OFILES)
	$(CC) -o mpilife $(OFILES)

clean::
	/bin/rm -f mpilife $(OFILES)
