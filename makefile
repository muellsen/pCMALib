#------------------------------------------------------------------------
# pCMALib: a parallel fortran 90 library for the evolution strategy with
#          covariance matrix adaptation
# Christian L. Mueller, Benedikt Baumgartner, Georg Ofenbeck
# MOSAIC group, ETH Zurich, Switzerland
#-------------------------------------------------------------------------

include make_brutus.inc

# Default value for the dependency pre-processor
# = same as code-generating pre-processor
DEPCPP ?= $(CPP)
#$(CPP) = fpp

# Sources dir
CMA_DIR := $(SRCDIR)/libcma
# Where to put .o and .d files
OBJ_DIR := $(BUILDDIR)/objects
# Target
TARGET  := $(OBJ_DIR)/$(LIB_CMA)
# Installation directories
INSTALL_DIR := $(BUILDDIR)/bin
MODULES_DIR := $(BUILDDIR)/include

CEC2005_DIR := $(SRCDIR)/libtestfcns
MPI_DIR := $(SRCDIR)/libcma/mpi
MATLAB_DIR := $(SRCDIR)/libcma/matlab
RANDGEN_DIR := $(SRCDIR)/librng
BBOB_DIR := $(SRCDIR)/BBOB
LJ_DIR := $(SRCDIR)/energy_landscapes/LJ
BFGS_DIR := $(SRCDIR)/bfgs
WATER_DIR := $(SRCDIR)/energy_landscapes/TIPnP
USER_DIR := $(SRCDIR)/user
# List the source files
SOURCES :=  $(notdir $(wildcard $(CEC2005_DIR)/*.f90))
SOURCES +=  $(notdir $(wildcard $(CMA_DIR)/*.f90))
SOURCES +=  $(notdir $(wildcard $(RANDGEN_DIR)/*.f90))
SOURCES +=  $(notdir $(wildcard $(LJ_DIR)/*.f90))
SOURCES +=  $(notdir $(wildcard $(BFGS_DIR)/*.f90))
SOURCES +=  $(notdir $(wildcard $(WATER_DIR)/*.f90))
SOURCES +=  $(notdir $(wildcard $(USER_DIR)/*.f90))
SOURCES_F :=  $(notdir $(wildcard $(RANDGEN_DIR)/*.f))
SOURCES_F +=  $(notdir $(wildcard $(BFGS_DIR)/*.f))
SOURCES_F +=  $(notdir $(wildcard $(LJ_DIR)/*.f))
SOURCES_F +=  $(notdir $(wildcard $(WATER_DIR)/*.f))

#Where to search for the .f90 files
.PATH.f90 = $(CEC2005_DIR);$(CMA_DIR);$(RANDGEN_DIR);$(MATLAB_DIR);$(LJ_DIR);$(BFGS_DIR);$(WATER_DIR);$(MPI_DIR);$(USER_DIR)


# Now list all the paths to check
VPATH := $(CMA_DIR):
VPATH += $(MODULES_DIR):
VPATH += $(CEC2005_DIR):
VPATH += $(RANDGEN_DIR):
VPATH += $(LJ_DIR):
VPATH += $(BFGS_DIR):
VPATH += $(WATER_DIR):
VPATH += $(USER_DIR):


###########################################################
# LIBS to be linked
###########################################################

LIBS = $(LAPACKLIB)

ifeq "$(HAS_MAT)" "1"
LIBS +=  $(MATLIB) -lmat -lmx -leng
DEFINE += -D__HAVE_MATLAB__
VPATH += $(MATLAB_DIR):
INCLS += $(MATINC)
SOURCES +=  $(notdir $(wildcard $(MATLAB_DIR)/*.f90))
endif


ifeq "$(HAS_MPI)" "1"
DEFINE += -D__HAVE_MPI__
INCLS += $(MPIINC)
LIBS += $(MPILIB) #-lmpi -lmpifarg
VPATH += $(MPI_DIR):
LINKER = $(MPI_LINKER)
SOURCES +=  $(notdir $(wildcard $(MPI_DIR)/*.f90))
endif

ifeq "$(BBOB)" "1"
SOURCES_C :=  $(notdir $(wildcard $(BBOB_DIR)/*.c))
VPATH += $(BBOB_DIR):
DEFINE += -D__BBOB
endif



# Useful derived macros
INCLS += $(patsubst %,-I%,$(subst :, ,$(VPATH)))
OBJECTS := $(SOURCES:%.f90=$(OBJ_DIR)/%.o)
OBJECTS += $(SOURCES_F:%.f=$(OBJ_DIR)/%.o)
OBJECTS += $(SOURCES_C:%.c=$(OBJ_DIR)/%.o)
#OBJECTS += prng_r.o  randomGenerator.o
MODULES := $(SOURCES:%.f90=$(MODULES_DIR)/%.mod)
MODSRCS := $(SOURCES:%.f90=$(MODULES_DIR)/__%.f90)
DEPENDENCIES := $(SOURCES:%.f90=$(OBJ_DIR)/%.d)

# If you want to keep the pre-compiled sources!
ifdef KEEPCPP
RMCPP := 'touch'
else
RMCPP := 'rm'
endif


$(warning Checking for directories...)
$(shell mkdir $(OBJ_DIR))
$(shell mkdir $(MODULES_DIR))


###########################################################
# Set Precision and additional compiler flags
###########################################################
ifeq ($(PREC), 1)
        DEFINE += -D__SP $()
else
        DEFINE += -D__DP $()
endif


###########################################################
# Set some variables concerning compilation with Matlab support
###########################################################
# check if Matlab is present and set libs and make targets appropriatly


.DEFAULT: ;
# Default target is everything
all: $(TARGET)


# final linking
$(TARGET): $(OBJECTS)
#	echo $(SOURCES)
#	echo "$(FC) -o $@ $(OBJECTS)"
	$(LINKER)  -o $@ $(OBJECTS) $(LIBS)
#	echo $(OBJECTS)
#	$(AR) $(ARFLAGS) $@ $(OBJECTS)



# How to make dependency files
# Some explanations: 
# 1) we use the preprocessor to find the includes
# 2) we add the dependency file as a target
# 3) find all the INCLUDE statements that are not inside a comment
# (actually we just check for lines starting with !)
# We exclude mpif.h from this search. 
# If make does not find mpif.h (user did not set the MPIINCLDIR correctly), 
# this could lead to an infinite loop for the dependcy files.
# 4) find all the USE statements that are not inside a comment
# 5) append a bogus line
$(OBJ_DIR)/%.d: %.f90
	@echo '[ quietly making' $@ ']'
	@$(DEPCPP) $(INCLS) -M $< | \
	sed -e 's#$*.o#$(OBJ_DIR)/$*.o $(OBJ_DIR)/$*.d#' \
	    -e 's#$$#\\#' -e's#\\\\#\\#' > $@
	@$(DEPCPP) -P $(DEFINE) $(INCLS) $< > $(OBJ_DIR)/__$*.f90
	@grep "INCLUDE " $(OBJ_DIR)/__$*.f90 | \
        sed -e 's#^[ \t]*##;s#[ \t]*##' \
	    -e  '/^!/d;/mpif.h/d' \
		-e  '/^!/d;/fftw3.f90/d' \
            -e 's#INCLUDE ##;s#$$#\\#' \
            -e 's#"##g' -e "s#'##g" -e 's# ##g' | \
        tr -d '\011' >> $@
	@echo '# end of source dependencies for .o and .d files' >> $@
	@echo ''$(OBJ_DIR)/$*.o ':\' >> $@ 
	@grep "USE " $(OBJ_DIR)/__$*.f90 | \
	sed -e 's#^[ \t]*##;s#[ \t]*##' \
	    -e '/^!/d' \
	    -e 's#,.*##' | \
	sed -e 's#USE #$(OBJ_DIR)/#' \
	    -e 's# ##g;s#$$#.o\\#' | \
	tr -d '\011' >> $@
	@echo '# end of module dependencies for .o file' >> $@
#	@rm $(OBJ_DIR)/__$*.f90




# How to make .o and .mod files out of f90
# Note how we move the .mod file: we use * and not $* because
# compilers lower-case the file name, thus creating a mismatch
# in certain cases, e.g. ppm_module_neighlist_MkNeighIdx.f
$(OBJ_DIR)/%.o: %.f90
	@echo "creating pp file __"
	$(CPP) -P -C $(DEFINE) $(INCLS) $< > $(OBJ_DIR)/__$*.f90
	@echo "compiling $(FC) $(INCLS) $(FFLAGS) -c -o $@ $(OBJ_DIR)/__$*.f90"
	$(FC) $(INCLS) $(FFLAGS) -c $(DEBUG) -o $@ $(OBJ_DIR)/__$*.f90
	@echo "moving"
	-@mv $(SRCDIR)/*.mod $(MODULES_DIR)
	@echo "deleting"
	@$(RMCPP) $(OBJ_DIR)/__$*.f90
#	@rm $(OBJ_DIR)/__$*.f90


# How to make .o and .mod files out of f90
$(OBJ_DIR)/%.o: %.f
	@echo "creating pp file __"
	$(CPP) -traditional-cpp -P -C $(DEFINE) $(INCLS) $< > $(OBJ_DIR)/__$*.f
	@echo "compiling $(FC) $(INCLS) $(FFLAGS) -c -extend-source -o $@ $(OBJ_DIR)/__$*.f"
	$(FC) $(INCLS) $(FFLAGS) $(DEBUG) -c -extend-source -ffixed-line-length-132   -o $@ $(OBJ_DIR)/__$*.f
	@echo "moving"
	-@mv $(SRCDIR)/*.mod $(MODULES_DIR)
	@echo "deleting"
	@$(RMCPP) $(OBJ_DIR)/__$*.f
#	@rm $(OBJ_DIR)/__$*.f

# How to make .o and .mod files out of f90
$(OBJ_DIR)/%.o: %.c
	@echo "creating objects out of c files_"
	$(CC) -c $< -o $(OBJ_DIR)/$*.o 



# at the install part we copy the input files and the programm itself to 
install:
	$(shell mkdir $(INSTALL_DIR))
	$(shell mv $(TARGET) $(INSTALL_DIR))
	$(shell cp -r $(CEC2005_DIR)/supportData $(INSTALL_DIR))


# Get rid of the .d's too. Notice this may produce weird behavior during
# the next invocation, since the makefile itself needs the .d's.
reallyclean:
	rm -f $(OBJECTS)
	rm -f $(MODULES)
	rm -f $(TARGET)
	rm -f $(DEPENDENCIES)

# Make new
new: reallyclean all install

# Do not run this makefile until all dependency files are up-to-date
-include $(DEPENDENCIES)
