########################################################################
# A Matter Flow Method for Staggered Lagrangian Hydrodynamics 
# on Triangular Grids (MATTERFLOW)
# 
# This is the Makefile for MATTERFLOW test! 
#
########################################################################

#==============================================================================#
# User compilers                                                               #
# MATTERFLOW has been tested with many different compilers.                    #
#==============================================================================#
CC  = icc 
CPP = icpc 
FC  = ifort -nofor_main
AR  = ar ruc

########################################################################      
# Compiling options                                                             
########################################################################        
BOPT=-O2


########################################################################
# Root directory for MATTERFLOW package
########################################################################
INCLUDE = -I ./include -I ./src 

#==============================================================================#
# User preprocessing definitions                                               #
#==============================================================================#
CDEFS=


COPTS=$(BOPT)
CINCLUDES=$(INCLUDE)
CFLAGS=$(CDEFS) $(COPTS) $(CINCLUDES)

FOPTS=$(BOPT)
FDEFS=$(CDEFS)
FINCLUDES=$(CINCLUDES)
FFLAGS=$(FDEFS) $(FOPTS) $(FINCLUDES)

########################################################################
# Link options
########################################################################
LINKOPTS=$(BOPT)

LIBS=$(TESTLIB) 

CLFLAGS=-lstdc++ $(LINKOPTS) $(LIBS)
FLFLAGS=-lm $(LINKOPTS) $(LIBS)

CSRCDIR = ./src
FSRCDIR = ./src
TESTLIB = ./lib/libmatterflow.a

########################################################################
# Rules for compiling the source files
########################################################################
.SUFFIXES: .c .cc .cpp .for .f .f77 .f90 .f95
#
FSRC := $(foreach dir,$(FSRCDIR),$(wildcard $(FSRCDIR)/*.for))
FSRC += $(foreach dir,$(FSRCDIR),$(wildcard $(FSRCDIR)/*.f))
FSRC += $(foreach dir,$(FSRCDIR),$(wildcard $(FSRCDIR)/*.f77))
FSRC += $(foreach dir,$(FSRCDIR),$(wildcard $(FSRCDIR)/*.f90))
FSRC += $(foreach dir,$(FSRCDIR),$(wildcard $(FSRCDIR)/*.f95))
CSRC := $(foreach dir,$(CSRCDIR),$(wildcard $(CSRCDIR)/*.c))
CSRC += $(foreach dir,$(EXTRDIR),$(wildcard $(EXTRDIR)/*.c))
#
OBJSF := $(patsubst %.for,%.o,$(FSRC))
OBJSF += $(patsubst %.f,%.o,$(FSRC))
OBJSF += $(patsubst %.f77,%.o,$(FSRC))
OBJSF += $(patsubst %.f90,%.o,$(FSRC))
OBJSF += $(patsubst %.f95,%.o,$(FSRC))
OBJSC := $(patsubst %.c,%.o,$(CSRC))
#
.for.o:
	@$(FC) -c $< -o $@ $(FFLAGS)
	@echo 'Building F object $@'
	@$(AR) $(TESTLIB) $@
#
.f.o:
	@$(FC) -c $< -o $@ $(FFLAGS)
	@echo 'Building F object $@'
	@$(AR) $(TESTLIB) $@
#
.f90.o:
	@$(FC) -c $< -o $@ $(FFLAGS)
	@echo 'Building F90 object $@'
	@$(AR) $(TESTLIB) $@
#
.f95.o:
	@$(FC) -c $< -o $@ $(FFLAGS)
	@echo 'Building F95 object $@'
	@$(AR) $(TESTLIB) $@
#
.c.o:
	@$(CC) -c $< -o $@ $(CFLAGS)
	@echo 'Building C object $@'
	@$(AR) $(TESTLIB) $@
#
.cpp.o:
	@$(CPP) -c $< -o $@ $(CFLAGS)
	@echo 'Building CPP object $@'
	@$(AR) $(TESTLIB) $@
#
########################################################################
# List of all programs to be compiled
########################################################################

# Everything
ALLPROG=$(TESTLIB)

########################################################################
# Link
########################################################################

all: $(ALLPROG) test

Default:
	test

headers: 
	cat $(CSRCDIR)/*.c \
	| awk -v name="matterflow_functs.h" -f ./util/mkheaders.awk > ./include/matterflow_functs.h

$(TESTLIB): $(OBJSC) $(OBJSF)
	@ranlib $(TESTLIB)
	@echo 'Generating library $@'

lib: $(OBJSC) $(OBJSF)
	ranlib $(TESTLIB)
	@echo 'Generating library $@'

########################################################################
# Some test problems
########################################################################

test:
	@$(CC) $(CFLAGS) -c main/Main.c -o main/Main.o
	@$(FC) $(LOPT) main/Main.o $(FLFLAGS) -o test.ex
	@echo 'Building executable $@'


########################################################################
# Clean up
########################################################################

.PHONY : clean distclean help

clean:
	@rm -f $(CSRCDIR)/*.o
	@rm -f $(FSRCDIR)/*.o
	@rm -f main/*.o

distclean:
	@make clean
	@rm -f lib/*.a
	@rm -f *~ *.ex *.out
	@rm -f $(CSRCDIR)/*~
	@rm -f $(FSRCDIR)/*~
	@rm -f main/*~

help:
	@echo "======================================================"
	@echo " A Matter Flow Method for SGH on Triangular Grids     "
	@echo "======================================================"
	@echo " "
	@echo " make            : build all exe files "
	@echo " make headers    : build the header file automatically"
	@echo " make clean      : clean all obj files "
	@echo " make distclean  : clean all obj, exe, bak, out files "
	@echo " make help       : show this screen "
	@echo " "
