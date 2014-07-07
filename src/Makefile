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
# Check how to reduce the number of invocations of loggamma.f90 to the Husimi case.
#
# ALGEBRAIC MODEL ENERGIES AND EIGENVECTORS (MODEL HAMILTONIAN) Chain I
#
AVECMODSRC = nrtype.f90 defparam_2dvm.f90 u3_2dvm_mod.f90 loggamma.f90 avalavec_modelH_2dvm.f90 
AVECMODOBJ = $(AVECMODSRC:.f90=.o)
#
#
# ALGEBRAIC MODEL ENERGIES AND EIGENVECTORS (MODEL HAMILTONIAN) Chain II
#
AVECMODIISRC = nrtype.f90 defparam_2dvm.f90 u3_2dvm_mod.f90 loggamma.f90 avalavec_modelH_2dvm_II.f90 
AVECMODIIOBJ = $(AVECMODIISRC:.f90=.o)
#
#
# ALGEBRAIC MODEL INVERSE PARTICIPATION RATIO (Husimi distribution MODEL HAMILTONIAN)
#
IPRHMODSRC = nrtype.f90 defparam_2dvm.f90 u3_2dvm_mod.f90 loggamma.f90 ipr_Husimi_modelH_2dvm.f90 
IPRHMODOBJ = $(IPRHMODSRC:.f90=.o)
#
#
# ALGEBRAIC MODEL INVERSE PARTICIPATION RATIO (MODEL HAMILTONIAN) Chain I
#
IPRMODSRC = nrtype.f90 defparam_2dvm.f90 u3_2dvm_mod.f90 loggamma.f90 ipr_modelH_2dvm.f90 
IPRMODOBJ = $(IPRMODSRC:.f90=.o)
#
#
# ALGEBRAIC MODEL INVERSE PARTICIPATION RATIO (MODEL HAMILTONIAN) Chain II
#
IPRMODIISRC = nrtype.f90 defparam_2dvm.f90 u3_2dvm_mod.f90 loggamma.f90 ipr_modelH_2dvm_II.f90 
IPRMODIIOBJ = $(IPRMODIISRC:.f90=.o)
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
all: avecmodh avecIImodh iprHmodh iprmodh iprIImodh
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
# HAMILTONIAN EIGENVALUES AND EIGENVECTORS (MODEL HAMILTONIAN, Chain I)
avecmodh: $(AVECMODOBJ)
	$(FC) -o $(BINPATH)/avec_U3_mh $(AVECMODOBJ) $(LIBS)
#################################################################
# HAMILTONIAN EIGENVALUES AND EIGENVECTORS (MODEL HAMILTONIAN, Chain II)
avecIImodh: $(AVECMODIIOBJ)
	$(FC) -o $(BINPATH)/avec_U3_mh_SO3 $(AVECMODIIOBJ) $(LIBS)
#################################################################
# HAMILTONIAN EIGENVALUES AND EIGENVECTORS (MODEL HAMILTONIAN)
iprHmodh: $(IPRHMODOBJ)
	$(FC) -o $(BINPATH)/ipr_Husimi_U3_mh $(IPRHMODOBJ) $(LIBS)
#################################################################
# HAMILTONIAN EIGENVALUES AND EIGENVECTORS (MODEL HAMILTONIAN, Chain I)
iprmodh: $(IPRMODOBJ)
	$(FC) -o $(BINPATH)/ipr_U3_mh $(IPRMODOBJ) $(LIBS)
#################################################################
# HAMILTONIAN EIGENVALUES AND EIGENVECTORS (MODEL HAMILTONIAN, Chain II)
iprIImodh: $(IPRMODIIOBJ)
	$(FC) -o $(BINPATH)/ipr_U3_mh_SO3 $(IPRMODIIOBJ) $(LIBS)
#################################################################
clean:
	@echo "Cleaning object files"
	@rm -f *.o *.mod




