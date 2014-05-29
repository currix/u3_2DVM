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
AVECMODSRC = nrtype.f90 defparam_2dvm.f90 u3_2dvm_mod.f90 avalavec_modelH_2dvm.f90 
AVECMODOBJ = $(AVECMODSRC:.f90=.o)
#
#
# ALGEBRAIC MODEL INVERSE PARTICIPATION RATIO (MODEL HAMILTONIAN)
#
IPRMODSRC = nrtype.f90 defparam_2dvm.f90 u3_2dvm_mod.f90 ipr_modelH_2dvm.f90 
IPRMODOBJ = $(IPRMODSRC:.f90=.o)
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
all: avecmodh iprmodh
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
avecmodh: $(AVECMODOBJ)
	$(FC) -o $(BINPATH)/avec_U3_mh $(AVECMODOBJ) $(LIBS)
#################################################################
# HAMILTONIAN EIGENVALUES AND EIGENVECTORS (MODEL HAMILTONIAN)
iprmodh: $(IPRMODOBJ)
	$(FC) -o $(BINPATH)/ipr_U3_mh $(IPRMODOBJ) $(LIBS)
#################################################################
clean:
	@echo "Cleaning object files"
	@rm -f $(AVECMODOBJ) $(IPRMODOBJ)




