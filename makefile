#If compiled with fujitsu f90
#include fujitsu.def

#If compiled with pgf90
#include pgf90.def

##If compiled with ifort
#include ifort.def

#If compiled with gfortran
include gfortran.def

# If compiled with gfortran
#FC=gfortran
#FFLAGS= -O2 -Wtabs  -ffixed-form # -fbounds-check
#FFLAGS1= -O2 -Wtabs  -ffixed-form # -fbounds-check

LBITS := $(shell getconf LONG_BIT)
ifeq ($(LBITS),64)
# Compile with 64-bits lapack
#LFLAGS= -L./lapack64 -llapack -lrefblas
else
# Compile with 32-bits lapack
#LFLAGS= -L./lapack -llapack -lrefblas
endif



# BINDIR= directory to install executables
BINDIR=$(PWD)/bin
TOPDIR=$(PWD)
SRC=$(TOPDIR)

.SUFFIXES: .c .o

OBJ=modules.o coul90.o scatcc.o  thox.o basis.o hdiag.o ho.o lst-amos.o utils.o \
    pauli.o sort.o nag2.o whittaker.o  \
     continuum.o belam.o ham.o  transition_targdef.o\
        transition.o   clebsg.o  projpot.o fragpot.o \
        rmatel.o solvecc.o lapack.o  rmatrix.o ccbins.o xsections.o \
	bincc2.o ceigen.o readwf.o
#coulfg4.f 
all: $(OBJ)
	$(FC) -o thox $(OBJ) $(LFLAGS) $(PARALLEL)  
#	$(FC) -g -o thox $(OBJ) $(LFLAGS) -L./lapack64 -llapack -lrefblas
#	$(FC) -g -o thox $(OBJ) $(LFLAGS) -L./lapack -llapack -lrefblas
modules.o:modules.f90 
	$(FC) $(FFLAGS) -c modules.f90
thox.o:thox.f90 modules.f90
	$(FC) $(FFLAGS) -c thox.f90
ho.o:ho.f90 modules.f90
	$(FC) $(FFLAGS) -c ho.f90
lst-amos.o:lst-amos.f90 modules.f90
	$(FC) $(FFLAGS) -c lst-amos.f90
basis.o:basis.f90 modules.f90
	$(FC) $(FFLAGS) -c basis.f90
ham.o:ham.f90 modules.f90
	$(FC) $(FFLAGS) -c ham.f90
hdiag.o:hdiag.f90 modules.f90
	$(FC) $(FFLAGS) -c hdiag.f90
utils.o:utils.f90 modules.f90
	$(FC) $(FFLAGS) -c utils.f90
pauli.o:pauli.f90 modules.f90
	$(FC) $(FFLAGS) -c pauli.f90
sort.o:sort.f90 modules.f90
	$(FC) $(FFLAGS) -c sort.f90
nag2.o:nag2.f90 modules.f90
	$(FC) $(FFLAGS) -c nag2.f90
whittaker.o:whittaker.f90 modules.f90
	$(FC) $(FFLAGS) -c whittaker.f90
coul90.o:coul90.f modules.f90
	$(FC) $(FFLAGS1) -c coul90.f
ccbins.o:ccbins.f90 modules.f90
	$(FC) $(FFLAGS1) -c ccbins.f90
continuum.o:continuum.f90 modules.f90
	$(FC) $(FFLAGS1) -c continuum.f90
scatcc.o:scatcc.f90 modules.f90
	$(FC) $(FFLAGS1) -c scatcc.f90 
coulfg4.o:coulfg4.f modules.f90
	$(FC) $(FFLAGS1) -c coulfg4.f
belam.o:belam.f90 modules.f90
	$(FC) $(FFLAGS1) -c belam.f90
transition.o:transition.f90 modules.f90
	$(FC) $(FFLAGS1)  -c transition.f90
transition_targdef.o:transition_targdef.f90 modules.f90
	$(FC) $(FFLAGS1)  -c transition_targdef.f90
fragpot.o:fragpot.f90 modules.f90
	$(FC) $(FFLAGS1) -c fragpot.f90
projpot.o:projpot.f90 modules.f90
	$(FC) $(FFLAGS1) -c projpot.f90
rmatel.o:rmatel.f90 modules.f90
	$(FC) $(FFLAGS1) -c rmatel.f90
read_frad.o:read_frad.f90 modules.f90
	$(FC) $(FFLAGS1) -c read_frad.f90
clebsg.o:clebsg.f90 modules.f90
	$(FC) $(FFLAGS1) -c clebsg.f90
solvecc.o:solvecc.f90 modules.f90
	$(FC) $(FFLAGS) -c solvecc.f90
xsections.o:xsections.f90
	$(FC) $(FFLAGS) $(PARALLEL) -c  xsections.f90
lapack.o:lapack.f 
	$(FC) $(FFLAGS) -c lapack.f
dsyev-lapack.o:dsyev-lapack.f 
	$(FC) $(FFLAGS) -c dsyev-lapack.f
rmatrix.o:rmatrix.f
	$(FC) $(FFLAGS) -c rmatrix.f
bincc2.o:bincc2.f
	$(FC) $(FFLAGS) -c bincc2.f
ceigen.o:ceigen.F
	$(FC) $(FFLAGS2) -c ceigen.F
readwf.o:readwf.f90
	$(FC) $(FFLAGS) -c readwf.f90

install: all
#%if BINDIR
	-mkdir $(BINDIR)
	 @echo 
	 @echo Installing executables on $(BINDIR) 
	 @echo -------------------------------------------------
	-cp $(SRC)/thox $(BINDIR)/
	@echo ---------------------------------------------------

clean:
	rm -f *.o core
	rm -f *.mod
# DO NOT DELETE



