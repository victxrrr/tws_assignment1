FC=gfortran
EXT=f90

FFLAGS_g95      = -g -pedantic -Wall -fbounds-check -ftrace=full
FFLAGS_gfortran = -g -pedantic -Wall -Wimplicit-interface -Wunderflow -fbounds-check -fimplicit-none
FFLAGS_ifort    = -g -debug full -implicitnone -check -warn -free
FFLAGS_nagfor   = -g -gline -u -info -colour


# Select the right flags for the current compiler
FFLAGS=$(FFLAGS_$(FC))

all: exp

exp.o: ./exp.$(EXT)
	$(FC) $(FFLAGS) -o $@ -c $^

exp: exp.o
	$(FC) $(FFLAGS) -o $@ $^

solvers.o: ./solvers.$(EXT)
	$(FC) $(FFLAGS) -o $@ -c $^

siqrd.o: ./siqrd.$(EXT)
	$(FC) $(FFLAGS) -o $@ -c $^

siqrd: solvers.o siqrd.o solver_gfortran.o
	$(FC) $(FFLAGS) -o $@ $^ 

# to avoid errors when using solver_gfortran.o
%.o: %.mod

plot: plot.tex
	pdflatex plot.tex >/dev/null

clean:
	@ rm -f exp.o siqrd.o solvers.o solvers.mod exp siqrd *.aux *.log plot.pdf