####################################################################
#
#	Makefile to compile, print or check 2DVM programs
#
#       gfortran compiler
#
#       2DVM stands for 2 Dimensional limit of the Vibron Model.
#
#       by Currix TM.
#
#
BINPATH = ../bin
####################################################################
####################################################################
#
# ALGEBRAIC MODEL ENERGIES AND EIGENVECTORS (MODEL HAMILTONIAN) Chain (i) U(3) -> U(2)
#
AVEC_modH_u2_SRC = nrtype.f90 defparam_2dvm.f90 u3_2dvm_mod.f90 loggamma.f90 avalavec_modelH_2dvm_u2.f90 
AVEC_modH_u2_OBJ = $(AVEC_modH_u2_SRC:.f90=.o)
#AVEC_modH_u2_OPENMP_OBJ = $(AVEC_modH_u2_SRC:.f90=.oOMP)
#
# ALGEBRAIC MODEL ENERGIES AND EIGENVECTORS (MODEL HAMILTONIAN) Chain (ii) U(3) -> SO(3)
#
AVEC_modH_so3_SRC = nrtype.f90 defparam_2dvm.f90 u3_2dvm_mod.f90 loggamma.f90 avalavec_modelH_2dvm_so3.f90 
AVEC_modH_so3_OBJ = $(AVEC_modH_so3_SRC:.f90=.o)
#AVEC_modH_so3_OPENMP_OBJ = $(AVEC_modH_so3_SRC:.f90=.oOMP)
#
#
# ALGEBRAIC MODEL INVERSE PARTICIPATION RATIO (Husimi distribution MODEL HAMILTONIAN)
#
IPR_modH_Hus_u2_SRC = nrtype.f90 defparam_2dvm.f90 u3_2dvm_mod.f90 loggamma.f90 ipr_Husimi_modelH_2dvm_u2.f90 
IPR_modH_Hus_u2_OBJ = $(IPR_modH_Hus_u2_SRC:.f90=.o)
#IPR_modH_Hus_u2_OPENMP_OBJ = $(IPR_modH_Hus_u2_SRC:.f90=.oOMP)
#
#
# ALGEBRAIC MODEL INVERSE PARTICIPATION RATIO (MODEL HAMILTONIAN) Chain  (i) U(3) -> U(2)
#
IPR_modH_u2_SRC = nrtype.f90 defparam_2dvm.f90 u3_2dvm_mod.f90 loggamma.f90 ipr_modelH_2dvm_u2.f90 
IPR_modH_u2_OBJ = $(IPR_modH_u2_SRC:.f90=.o)
#IPR_modH_u2_OPENMP_OBJ = $(IPR_modH_u2_SRC:.f90=.oOMP)
#
#
# ALGEBRAIC MODEL INVERSE PARTICIPATION RATIO (MODEL HAMILTONIAN) Chain  (ii) U(3) -> SO(3)
#
IPR_modH_so3_SRC = nrtype.f90 defparam_2dvm.f90 u3_2dvm_mod.f90 loggamma.f90 ipr_modelH_2dvm_so3.f90 
IPR_modH_so3_OBJ = $(IPR_modH_so3_SRC:.f90=.o)
#IPR_modH_so3_OPENMP_OBJ = $(IPR_modH_so3_SRC:.f90=.oOMP)
#
# ALGEBRAIC MODEL ENERGIES AND EIGENVECTORS (GENERAL HAMILTONIAN) Chain  (i) U(3) -> U(2)
#
AVEC_u2_SRC = nrtype.f90 defparam_2dvm.f90 u3_2dvm_mod.f90 loggamma.f90 avalavec_2dvm_u2.f90 
AVEC_u2_OBJ = $(AVEC_u2_SRC:.f90=.o)
#AVEC_u2_OPENMP_OBJ = $(AVEC_u2_SRC:.f90=.oOMP)
#
# ALGEBRAIC MODEL ENERGIES AND EIGENVECTORS (GENERAL HAMILTONIAN) Chain  (ii) U(3) -> SO(3)
#
AVEC_so3_SRC = nrtype.f90 defparam_2dvm.f90 u3_2dvm_mod.f90 loggamma.f90 avalavec_2dvm_so3.f90 
AVEC_so3_OBJ = $(AVEC_so3_SRC:.f90=.o)
#AVEC_so3_OPENMP_OBJ = $(AVEC_so3_SRC:.f90=.oOMP)
#
# ALGEBRAIC MODEL PARTICIPATION RATIO (GENERAL HAMILTONIAN) Chain  (i) U(3) -> U(2)
#
IPR_u2_SRC = nrtype.f90 defparam_2dvm.f90 u3_2dvm_mod.f90 loggamma.f90 ipr_2dvm_u2.f90 
IPR_u2_OBJ = $(IPR_u2_SRC:.f90=.o)
#IPR_u2_OPENMP_OBJ = $(IPR_u2_SRC:.f90=.oOMP)
#
#
# ALGEBRAIC MODEL PARTICIPATION RATIO (GENERAL HAMILTONIAN) Chain  (ii) U(3) -> SO(3)
#
IPR_so3_SRC = nrtype.f90 defparam_2dvm.f90 u3_2dvm_mod.f90 loggamma.f90 ipr_2dvm_so3.f90 
IPR_so3_OBJ = $(IPR_so3_SRC:.f90=.o)
#IPR_so3_OPENMP_OBJ = $(IPR_so3_SRC:.f90=.oOMP)
#
# ALGEBRAIC MODEL PARTICIPATION RATIO (GENERAL 4-Body HAMILTONIAN) Chain (i) U(3) -> U(2)
#
IPR_4b_u2_SRC = nrtype.f90 defparam_2dvm.f90 u3_2dvm_mod.f90 loggamma.f90 ipr_2dvm_4b_u2.f90 
IPR_4b_u2_OBJ = $(IPR_4b_u2_SRC:.f90=.o)
#IPR_4b_u2_OPENMP_OBJ = $(IPR_4b_u2_SRC:.f90=.oOMP)
#
#
#
# ALGEBRAIC MODEL PARTICIPATION RATIO (GENERAL 4-Body HAMILTONIAN) Chain (i) U(3) -> SO(3)
#
IPR_4b_so3_SRC = nrtype.f90 defparam_2dvm.f90 u3_2dvm_mod.f90 loggamma.f90 ipr_2dvm_4b_so3.f90 
IPR_4b_so3_OBJ = $(IPR_4b_so3_SRC:.f90=.o)
#IPR_4b_so3_OPENMP_OBJ = $(IPR_4b_so3_SRC:.f90=.oOMP)
#
#	
# ALGEBRAIC MODEL Chi^2 TO EXPERIMENT (GENERAL 4BODY HAMILTONIAN) Chains  (i) U(3) -> U(2) and (ii) U(3) -> SO(3)
#
CHI2_SRC = nrtype.f90 loggamma.f90 defparam_2dvm.f90 u3_2dvm_mod.f90 fit_2dvm.f90 chi_sqre_prg.f90 
CHI2_OBJ = $(CHI2_SRC:.f90=.o)
#CHI2_OPENMP_OBJ = $(CHI2_SRC:.f90=.oOMP)
#
#
#Minimization with Minuit-77
Min_minuit_SRC = nrtype.f90 loggamma.f90 defparam_2dvm.f90 u3_2dvm_mod.f90 fit_2dvm.f90 minuit-cern.f FCN_func.f90 minimization.f90 
Min_minuit_OBJ = nrtype.o loggamma.o defparam_2dvm.o u3_2dvm_mod.o fit_2dvm.o minuit-cern.o FCN_func.o minimization.o

####################################################################
####################################################################
#
.SUFFIXES: 
.SUFFIXES: .o .f90 .oOMP .f
#
FC = gfortran
FC_Path = /usr/bin/gfortran
#
FOPT=-O3 -Wall
#
PREP = -cpp
#
ARCH = -march=native
#
#
# Module directories
MODLIB   = -I/usr/local/lib/lapack95/lapack95_modules -I/usr/local/include/lapack95_modules
#
OPTFLAGS = $(FOPT)
OPTFLAGS += $(ARCH)
OPTFLAGS += $(MODLIB)
#
# Libraries
LAPACK95 = -L/usr/local/lib/lapack95 -llapack95 -L/usr/local/lib -llapack95
## Debian Lapack and Blas Libraries
LAPACK77 = -llapack
BLAS     = -lblas
#
LIBS     = $(LAPACK95)  $(LAPACK77) $(BLAS)
#
######################################################################
######################################################################
#
all: avec_U3_U2 avec_U3_U2_mh avec_U3_SO3  avec_U3_SO3_mh ipr_Husimi_U3_U2_mh ipr_U3_U2_mh ipr_U3_SO3_mh ipr_U3_U2 ipr_U3_SO3 ipr_4b_U3_U2 ipr_4b_U3_SO3 chi2_U3 Min_minuit
#
.PHONY : all
#################################################################
.f90.o:
	$(info )
	$(info Compiling single thread object file:)
	$(FC_Path) -c $(PREP) $(OPTFLAGS)  "$<"
#################################################################
.f.o:
	$(info )
	$(info Compiling single thread fortran 77 object file:)
	$(FC_Path) -c $(PREP) $(OPTFLAGS) "$<"
#################################################################
.f90.oOMP:
	$(info )
	$(info Compiling multi thread object file:)
	$(FC_Path) -c -fopenmp $(PREP) $(OPTFLAGS) -o "$@" "$<"
#################################################################
#
#################################################################
# HAMILTONIAN EIGENVALUES AND EIGENVECTORS (2 BODY HAMILTONIAN, Chain (i)
avec_U3_U2:  $(AVEC_u2_OBJ) Makefile 
	$(info )
	$(info Linking single thread avec_U3_U2 executable:)
	$(FC_Path) -o $(BINPATH)/$@_$(FC) $(AVEC_u2_OBJ) $(OPTFLAGS) $(LIBS) 
# avec_U3_U2_openmp:  $(AVEC_u2_OPENMP_OBJ) Makefile .PHONY 
# 	$(info )
# 	$(info Linking multithread avec_U3_U2 executable:)
# 	$(FC) -fopenmp -o $(BINPATH)/$@_$(FC) $(AVEC_u2_OPENMP_OBJ) $(OPTFLAGS) $(LIBS) 
#################################################################
# HAMILTONIAN EIGENVALUES AND EIGENVECTORS (2 BODY HAMILTONIAN, Chain II)
avec_U3_SO3: $(AVEC_so3_OBJ) Makefile
	$(info )
	$(info Linking single thread avec_U3_SO3 executable:)
	$(FC_Path) -o $(BINPATH)/$@_$(FC) $(AVEC_so3_OBJ) $(OPTFLAGS) $(LIBS)
# avec_U3_SO3_openmp: $(AVEC_so3_OPENMP_OBJ) Makefile .PHONY 
# 	$(info )
# 	$(info Linking multi thread avec_U3_SO3 executable:)
# 	$(FC) -fopenmp -o $(BINPATH)/$@_$(FC) $(AVEC_so3_OPENMP_OBJ) $(OPTFLAGS) $(LIBS)
#################################################################
# HAMILTONIAN EIGENVALUES AND EIGENVECTORS (MODEL HAMILTONIAN, Chain (i)
avec_U3_U2_mh:  $(AVEC_modH_u2_OBJ) Makefile
	$(info )
	$(info Linking single thread avec_U3_U2_mh executable:)
	$(FC_Path) -o $(BINPATH)/$@_$(FC) $(AVEC_modH_u2_OBJ) $(OPTFLAGS) $(LIBS) 
# avec_U3_U2_mh_openmp:  $(AVEC_modH_u2_OPENMP_OBJ) Makefile .PHONY 
# 	$(info )
# 	$(info Linking multithread avec_U3_U2_mh executable:)
# 	$(FC) -fopenmp -o $(BINPATH)/$@_$(FC) $(AVEC_modH_u2_OPENMP_OBJ) $(OPTFLAGS) $(LIBS) 
#################################################################
# HAMILTONIAN EIGENVALUES AND EIGENVECTORS (MODEL HAMILTONIAN, Chain (i)
avec_U3_SO3_mh: $(AVEC_modH_so3_OBJ) Makefile
	$(info )
	$(info Linking single thread avec_U3_SO3_mh executable:)
	$(FC_Path) -o $(BINPATH)/$@_$(FC) $(AVEC_modH_so3_OBJ) $(OPTFLAGS) $(LIBS)
# avec_U3_SO3_mh_openmp: $(AVEC_modH_so3_OPENMP_OBJ) Makefile .PHONY 
# 	$(info )
# 	$(info Linking multi thread avec_U3_SO3_mh executable:)
# 	$(FC) -fopenmp -o $(BINPATH)/$@_$(FC) $(AVEC_modH_so3_OPENMP_OBJ) $(OPTFLAGS) $(LIBS)
#################################################################
# MODEL HAMILTONIAN HUSIMI IBR Chain (i)
ipr_Husimi_U3_U2_mh: $(IPR_modH_Hus_u2_OBJ) Makefile 
	$(info )
	$(info Linking single thread $@ executable:)
	$(FC_Path) -o $(BINPATH)/$@_$(FC) $(IPR_modH_Hus_u2_OBJ) $(OPTFLAGS) $(LIBS)
# ipr_Husimi_U3_U2_mh_openmp: $(IPR_modH_Hus_u2_OPENMP_OBJ) Makefile  .PHONY
# 	$(info )
# 	$(info Linking multi thread $@ executable:)
# 	$(FC) -fopenmp -o $(BINPATH)/$@_$(FC) $(IPR_modH_Hus_u2_OPENMP_OBJ) $(OPTFLAGS) $(LIBS)
#################################################################
# HAMILTONIAN EIGENVALUES AND EIGENVECTORS. MODEL HAMILTONIAN, Chain (i)
ipr_U3_U2_mh: $(IPR_modH_u2_OBJ) Makefile
	$(info )
	$(info Linking single thread $@ executable:)
	$(FC_Path) -o $(BINPATH)/$@_$(FC) $(IPR_modH_u2_OBJ) $(OPTFLAGS) $(LIBS)
# ipr_U3_U2_mh_openmp: $(IPR_modH_u2_OPENMP_OBJ) Makefile  .PHONY 
# 	$(info )
# 	$(info Linking multi thread $@ executable:)
# 	$(FC) -fopenmp -o $(BINPATH)/$@_$(FC) $(IPR_modH_u2_OPENMP_OBJ) $(OPTFLAGS) $(LIBS)
#################################################################
# HAMILTONIAN EIGENVALUES, EIGENVECTORS, AND PR. MODEL HAMILTONIAN, Chain (ii)
ipr_U3_SO3_mh: $(IPR_modH_so3_OBJ) Makefile
	$(info )
	$(info Linking single thread $@ executable:)
	$(FC_Path) -o $(BINPATH)/$@_$(FC) $(IPR_modH_so3_OBJ) $(OPTFLAGS) $(LIBS)
# ipr_U3_SO3_mh_openmp: $(IPR_modH_so3_OPENMP_OBJ) Makefile  .PHONY 
# 	$(info )
# 	$(info Linking multi thread $@ executable:)
# 	$(FC) -fopenmp -o $(BINPATH)/$@_$(FC) $(IPR_modH_so3_OPENMP_OBJ) $(OPTFLAGS) $(LIBS)
#################################################################
# HAMILTONIAN EIGENVALUES, EIGENVECTORS, AND PR (GENERAL 2 BODY HAMILTONIAN, Chain I)
ipr_U3_U2: $(IPR_u2_OBJ) Makefile 
	$(info )
	$(info Linking single thread $@ executable:)
	$(FC_Path) -o $(BINPATH)/$@_$(FC) $(IPR_u2_OBJ) $(OPTFLAGS) $(LIBS)
# ipr_U3_U2_openmp: $(IPR_u2_OPENMP_OBJ) Makefile  .PHONY 
# 	$(info )
# 	$(info Linking multi thread $@ executable:)
# 	$(FC) -fopenmp -o $(BINPATH)/$@_$(FC) $(IPR_u2_OPENMP_OBJ) $(OPTFLAGS) $(LIBS)
#################################################################
# HAMILTONIAN EIGENVALUES, EIGENVECTORS, AND PR (GENERAL 2 BODY HAMILTONIAN, Chain II)
ipr_U3_SO3: $(IPR_so3_OBJ) Makefile
	$(info )
	$(info Linking single thread $@ executable:)
	$(FC_Path) -o $(BINPATH)/$@_$(FC) $(IPR_so3_OBJ) $(OPTFLAGS) $(LIBS)
# ipr_U3_SO3_openmp: $(IPR_so3_OPENMP_OBJ) Makefile  .PHONY 
# 	$(info )
# 	$(info Linking multi thread $@ executable:)
# 	$(FC) -fopenmp -o $(BINPATH)/$@_$(FC) $(IPR_so3_OPENMP_OBJ) $(OPTFLAGS) $(LIBS)
#################################################################
# HAMILTONIAN EIGENVALUES, EIGENVECTORS, AND PR (GENERAL 4 BODY HAMILTONIAN, Chain I)
ipr_4b_U3_U2: $(IPR_4b_u2_OBJ) Makefile
	$(info )
	$(info Linking single thread $@ executable:)
	$(FC_Path) -o $(BINPATH)/$@_$(FC) $(IPR_4b_u2_OBJ) $(OPTFLAGS) $(LIBS)
# ipr_4b_U3_U2_openmp: $(IPR_4b_u2_OPENMP_OBJ) Makefile  .PHONY 
# 	$(info )
# 	$(info Linking multi thread $@ executable:)
# 	$(FC) -fopenmp -o $(BINPATH)/$@_$(FC) $(IPR_4b_u2_OPENMP_OBJ) $(OPTFLAGS) $(LIBS)
#################################################################
# HAMILTONIAN EIGENVALUES, EIGENVECTORS, AND PR (GENERAL 4 BODY HAMILTONIAN, Chain II)
ipr_4b_U3_SO3: $(IPR_4b_so3_OBJ) Makefile
	$(info )
	$(info Linking single thread $@ executable:)
	$(FC_Path) -o $(BINPATH)/$@_$(FC) $(IPR_4b_so3_OBJ) $(OPTFLAGS) $(LIBS)
# ipr_4b_U3_SO3_openmp: $(IPR_4b_so3_OPENMP_OBJ) Makefile  .PHONY 
# 	$(info )
# 	$(info Linking multi thread $@ executable:)
# 	$(FC) -fopenmp -o $(BINPATH)/$@_$(FC) $(IPR_4b_so3_OPENMP_OBJ) $(OPTFLAGS) $(LIBS)
#################################################################
#################################################################
# SQUARE CHI VALUE (GENERAL 4 BODY HAMILTONIAN, Chain I)
chi2_U3: $(CHI2_OBJ) Makefile
	$(info )
	$(info Linking single thread $@ executable:)
	$(FC_Path) -Wall -o $(BINPATH)/$@_$(FC) $(CHI2_OBJ) $(OPTFLAGS) $(LIBS)
# chi2_U3_openmp: $(CHI2_OPENMP_OBJ) Makefile  .PHONY 
# 	$(info )
# 	$(info Linking multi thread $@ executable:)
# 	$(FC) -fopenmp -o $(BINPATH)/$@_$(FC) $(CHI2_OPENMP_OBJ) $(OPTFLAGS) $(LIBS)
#################################################################
#Minuit minimization:
Min_minuit: $(Min_minuit_OBJ) Makefile
		$(info )
		$(info Linking single thread $@ executable:)
		$(FC_Path) -o $(BINPATH)/$@_$(FC) $(Min_minuit_OBJ) $(OPTFLAGS) $(LIBS)
#Minuit minimization:
Min_minuit_static: $(Min_minuit_OBJ) Makefile
		$(info )
		$(info Linking single thread $@ executable:)
		$(FC_Path) -static -o $(BINPATH)/$@_$(FC) $(Min_minuit_OBJ) $(OPTFLAGS) $(LIBS)
#################################################################

clean:
	@echo "Cleaning object files"
	@rm -f *.o *.mod




