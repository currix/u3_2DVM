
PARALLEL=y

BINPATH	= ../bin

##################################
#**********GNU compiler**********#
#************gfortran************#
##################################

# FC_path	= 
# FC	= gfortran
# OPTS	= -O3 -Wall
# PREP	= -cpp
# ARCHFLAGS	= -march=native
# #Lapack95 modules:
# MODLIB   = -I/home/jamil/usr/lib/lapack95_modules/

# #Lapack95, Lapack77 and Blas libs
# LAPACK95 = /home/jamil/usr/lib/lapack95.a
# LAPACK77 = /home/jamil/usr/lib/lapack.a
# BLAS     = /home/jamil/usr/lib/blas.a
# TMGLIB	 = /home/jamil/usr/lib/tmglib.a

##################################
#**********Intel compiler********#
#**************ifort*************#
##################################

FC_path	= 
FC	= ifort
OPTS	= -O3 
PREP	= -cpp
ARCHFLAGS	= 

#MKL modules
MODLIB	= -I${MKLROOT}/include/intel64/lp64/

#Lapack95:

LAPACK95 = $(MKLROOT)/lib/intel64/libmkl_lapack95_lp64.a -L${MKLROOT}/lib/intel64 -lmkl_intel_lp64 -mkl
