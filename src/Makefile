####################################################################
#
#	Makefile to compile, print or check the 2DVM programs that
#	compute bending vibrational spectra for linear/bent
#	triatomic molecules using the U(3) dynamical algebra. 
#       2DVM = 2 Dimensional limit of the Vibron Model.
#
#       by Currix TM.
#
#
BINPATH = ../bin
#
VERSION = "3.0"
#
####################################################################
####################################################################
#
# ALGEBRAIC MODEL ENERGIES AND EIGENVECTORS (MODEL HAMILTONIAN)
#
AVECSRC = nrtype.f90 defparam_2dvm.f90 u3_2dvm_mod.f90 avalavec_modelH_2dvm.f90 
AVECOBJ = $(AVECSRC:.f90=.o)
#
####################################################################
####################################################################
FC = gfortran
#
FOPT	= -c 
FOPT   += -Wall
MODLIB   = -I/usr/local/lib/lapack95/lapack95_modules -I/usr/local/include/lapack95_modules
FOPT   += $(MODLIB)
#
LAPACK95 = -L/usr/local/lib/lapack95 -llapack95 -L/usr/local/lib -llapack95
## Cluster
LAPACK77 = -llapack
BLAS     = -lblas
#
LIBS     = $(LAPACK95)  $(LAPACK77) $(BLAS)
#
######################################################################
######################################################################
#
all: avecmodh
#
.PHONY : all
#
# Standard compilation rules
#
%.o: %.f90
	$(FC) $(FOPT) $<
#
%.o: %.F90
	$(FC) $(FOPT) $<
#
%.o: %.f
	$(FC) $(FOPT) $<
#################################################################
# HAMILTONIAN EIGENVALUES AND EIGENVECTORS (MODEL HAMILTONIAN)
avecmodh: $(AVECOBJ)
	$(FC) -o $(BINPATH)/avec_U3_mh $(AVECOBJ) $(LIBS)
#################################################################
clean:
	@echo "Cleaning object files"
	@rm -f $(AVECOBJ)




