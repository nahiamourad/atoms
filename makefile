nanoOBJ = program.o  main.o test.o perturbation.o bisection.o  accidental-pd.o accidental.o accidentalplus.o golden.o occupation-C.o occupationnbplus.o rayleighqoutient.o oda.o occupationnb.o shapefunction.o legendre.o exchange-matrix.o integration.o gauleg.o assembling-matrices.o used-functions.o Hartree-potential.o poisson-equation.o energy.o eigenfunctions.o      #ssdrvfinal.o #dsdrv33.o
LIBS = -llapack -lblas -larpack

OPT = -O  -Werror -Wall -pedantic -std=f95 

#FF77 = fort77
#CC = gcc

FF77 = gfortran
CC = gcc

#FF77 = pgf77
#CC = gcc
all: nano
nano	: $(nanoOBJ)
	$(FF77) $(OPT)  -o nano $(nanoOBJ) $(LIBS)
main.o: main.f95
	$(FF77)  $(OPT) -c $^
test.o: test.f95
	$(FF77)  $(OPT) -c $^	
eigenfunctions.o: eigenfunctions.f95
	$(FF77)  $(OPT) -c $^
perturbation.o: perturbation.f95
	$(FF77)  $(OPT) -c $^
accidental.o: accidental.f95
	$(FF77)  $(OPT) -c $^
accidental-pd.o: accidental-pd.f95
	$(FF77)  $(OPT) -c $^	
accidentalplus.o: accidentalplus.f95
	$(FF77)  $(OPT) -c $^	
program.o: program.f95
	$(FF77)  $(OPT) -c $^
bisection.o: bisection.f95
	$(FF77)  $(OPT) -c $^ 	
golden.o: golden.f90
	$(FF77)  $(OPT) -c $^	
occupation-C.o: occupation-C.f95
	$(FF77)  $(OPT) -c $^	
occupationnb.o: occupationnb.f95
	$(FF77)  $(OPT) -c $^
occupationnbplus.o: occupationnbplus.f95
	$(FF77)  $(OPT) -c $^	
rayleighqoutient.o: rayleighqoutient.f95
	$(FF77)  $(OPT) -c $^ 	
oda.o: oda.f95
	$(FF77)  $(OPT) -c $^	
shapefunction.o: shapefunction.f95
	$(FF77)   $(OPT) -c $^
legendre.o: legendre.f95
	$(FF77) $(OPT) -c $^
exchange-matrix.o: exchange-matrix.f95
	$(FF77) $(OPT)  -c $^
integration.o: integration.f95
	$(FF77) $(OPT)  -c $^
gauleg.o: gauleg.f95
	$(FF77)  $(OPT) -c $^
assembling-matrices.o: assembling-matrices.f95
	$(FF77) $(OPT) -c $^	
used-functions.o: used-functions.f95
	$(FF77) $(OPT) -c $^
Hartree-potential.o: Hartree-potential.f95
	$(FF77) $(OPT) -c $^
poisson-equation.o: poisson-equation.f95
	$(FF77) $(OPT) -c $^	
energy.o: energy.f95
	$(FF77) $(OPT) -c $^	

#dsdrv33.o:dsdrv33.f
	$(FF77) $(OPT) -c $^
#ssdrvfinal.o:ssdrvfinal.f
	$(FF77) $(OPT) -c $^

.PHONY: clean
clean: 
	rm *.o 	*.txt
