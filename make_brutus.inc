#################################################
#  PCMALIB make include file 
#################################################

# Name of the program
LIB_CMA := libpcma


######### SHELL OPTIONS ################

SHELL   := /bin/sh

######## OUTPUT Folder ################
# will create directories BUILDDIR/bin, BUILDDIR/objects and BUILDDIR/include here
BUILDDIR:= .

####### INPUT Folder ##################
# expects all the source files in that folder and its subfolders
SRCDIR := .


######### COMPILER AND COMPILER OPTIONS #########

# If you want to keep the pre-compiled sources! (comment if not)

# C compiler - needed if BBOB is encluded
CC = icc
# Fortran Compiler
FC = ifort
# -O3 is optimization level and -heap-arrays puts local variables on the heap instead of the stack
CFLAGS = -O3 -heap-arrays 0
# One for dependency generation
DEPCPP  := cpp
# and one for F90 code generation
CPP     := cpp


##### comment or set to empty string if not running debug
DEBUG  := -g




######## LINKER FLAGS #############################

LINKER = ifort
MPI_LINKER = /cluster/apps/openmpi/1.3.3/x86_64/intel_11.0.081/bin/mpif90


################### MATLAB ######################
# Declare if Matlab is present on your machine (0=no,1=yes)
# If yes, the program will produce a .mat output file.
# WATCH OUT NOT TO HAVE A SPACE AFTER 1/0

HAS_MAT = 1
# If Matlab is present, please specify the appropriate include paths

MATINC = -I/cluster/apps/matlab/7.6/x86_64/extern/include
MATLIB = -L/cluster/apps/matlab/7.6/x86_64/bin/glnxa64


################### LAPACK ######################
# You might need to tell the linker where to find your LAPACK library if not
# already included in your PATH variable

LAPACKLIB  = -lmkl_intel_lp64 -lmkl_sequential -lmkl_core -lmkl_lapack

################## PRECISION ####################
# Set the precision that should be used for computation
# (Single Precision = 1, Double Precision = 2)
# Please do a make clean whenever you change the precision
# THE SINGLE SETTING MOST LIKLY WONT WORK ANY MORE!
# this is due to 3rd party code included in pcmalib
PREC = 2


##################### MPI #######################
# If MPI is installed, set MPI = 1, otherwise MPI = 0
# WATCH OUT NOT TO HAVE A SPACE AFTER 1/0
# Define where to find the MPI include files needed for compilation
# And tell the linker where to find the MPI-Libraries if they are not
# already included in your PATH variable

HAS_MPI = 0
MPIINC = -I/cluster/apps/openmpi/1.3.3/x86_64/intel_11.0.081/include/
MPILIB = -L/cluster/apps/openmpi/1.3.3/x86_64/intel_11.0.081/lib/

##################### BBOB #######################
# include the BBOB benchmarks in the code? (yes = 1)
# WATCH OUT NOT TO HAVE A SPACE AFTER 1/0
# you need a working C compiler to include this
BBOB = 0
