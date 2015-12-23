####################################################################
#
#	Makefile to compile, print or check 2DVM programs 
#       2DVM stands for 2 Dimensional limit of the Vibron Model.
#
#       by Currix TM.
#
#
BINPATH = ../bin
#
FC = gfortran
#
##OPT = -O3 -ftree-vectorize -funroll-loops
OPT=-O3
#
ARCH = -march=native
#
#
# Module directories
MODLIB   = -I/usr/local/lib/lapack95/lapack95_modules -I/usr/local/include/lapack95_modules
#
OPTFLAGS = $(OPT)
OPTFLAGS += $(ARCH)
OPTFLAGS += $(MODLIB)
#
FOPT	= -c 
##FOPT   += -Wall
FOPT   += $(OPTFLAGS)
#
# Libraries
LAPACK95 = -L/usr/local/lib/lapack95 -llapack95 -L/usr/local/lib -llapack95
## Debian Lapack and Blas Libraries
LAPACK77 = -llapack
BLAS     = -lblas
#
LIBS     = $(LAPACK95)  $(LAPACK77) $(BLAS)
####################################################################
####################################################################
#
# To Do : Check how to reduce the number of invocations of loggamma.f90 to the Husimi case.
#
# ALGEBRAIC MODEL ENERGIES AND EIGENVECTORS (MODEL HAMILTONIAN) Chain I
#
AVECMODSRC = nrtype.f90 defparam_2dvm.f90 u3_2dvm_mod.f90 loggamma.f90 avalavec_modelH_2dvm.f90 
#
# ALGEBRAIC MODEL ENERGIES AND EIGENVECTORS (MODEL HAMILTONIAN) Chain II
#
AVECMODIISRC = nrtype.f90 defparam_2dvm.f90 u3_2dvm_mod.f90 loggamma.f90 avalavec_modelH_2dvm_II.f90 
#
#
# ALGEBRAIC MODEL INVERSE PARTICIPATION RATIO (Husimi distribution MODEL HAMILTONIAN)
#
IPRHMODSRC = nrtype.f90 defparam_2dvm.f90 u3_2dvm_mod.f90 loggamma.f90 ipr_Husimi_modelH_2dvm.f90 
#
#
# ALGEBRAIC MODEL INVERSE PARTICIPATION RATIO (MODEL HAMILTONIAN) Chain I
#
IPRMODSRC = nrtype.f90 defparam_2dvm.f90 u3_2dvm_mod.f90 loggamma.f90 ipr_modelH_2dvm.f90 
#
#
# ALGEBRAIC MODEL INVERSE PARTICIPATION RATIO (MODEL HAMILTONIAN) Chain II
#
IPRMODIISRC = nrtype.f90 defparam_2dvm.f90 u3_2dvm_mod.f90 loggamma.f90 ipr_modelH_2dvm_II.f90 
#
####################################################################
####################################################################
#
######################################################################
######################################################################
#
all: avec_U3_mh avec_U3_mh_SO3 ipr_Husimi_U3_mh ipr_U3_mh  ipr_U3_mh_SO3
#
.PHONY : all
#
#################################################################
# HAMILTONIAN EIGENVALUES AND EIGENVECTORS (MODEL HAMILTONIAN, Chain I)
avec_U3_mh:  $(AVECMODSRC) Makefile .PHONY 
	$(FC) -o $(BINPATH)/$@ $(AVECMODSRC) $(OPTFLAGS) $(LIBS) 
#################################################################
# HAMILTONIAN EIGENVALUES AND EIGENVECTORS (MODEL HAMILTONIAN, Chain II)
avec_U3_mh_SO3: $(AVECMODIISRC) Makefile  .PHONY 
	$(FC) -o $(BINPATH)/$@ $(AVECMODIISRC) $(OPTFLAGS) $(LIBS)
#################################################################
# HAMILTONIAN EIGENVALUES AND EIGENVECTORS (MODEL HAMILTONIAN)
ipr_Husimi_U3_mh: $(IPRHMODSRC) Makefile  .PHONY 
	$(FC) -o $(BINPATH)/$@ $(IPRHMODSRC) $(OPTFLAGS) $(LIBS)
#################################################################
# HAMILTONIAN EIGENVALUES AND EIGENVECTORS (MODEL HAMILTONIAN, Chain I)
ipr_U3_mh: $(IPRMODSRC) Makefile  .PHONY 
	$(FC) -o $(BINPATH)/$@ $(IPRMODSRC) $(OPTFLAGS) $(LIBS)
#################################################################
# HAMILTONIAN EIGENVALUES AND EIGENVECTORS (MODEL HAMILTONIAN, Chain II)
ipr_U3_mh_SO3: $(IPRMODIISRC) Makefile  .PHONY 
	$(FC) -o $(BINPATH)/$@ $(IPRMODIISRC) $(OPTFLAGS) $(LIBS)
#################################################################
clean:
	@echo "Cleaning object files"
	@rm -f *.o *.mod




