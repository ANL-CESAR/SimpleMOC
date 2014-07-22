#===============================================================================
# User Options
#===============================================================================

COMPILER    = gnu
MPI         = no
OPENMP      = yes
OPTIMIZE    = no
DEBUG       = yes
PROFILE     = no
PAPI        = no
BENCHMARK   = no

#===============================================================================
# Program name & source code list
#===============================================================================

program = SimpleMOC

source = \
main.c \
comms.c \
source.c \
test.c \
init.c \
io.c \
tracks.c \
utils.c \
solver.c

obj = $(source:.c=.o)

#===============================================================================
# Sets Flags
#===============================================================================

# Regular gcc Compiler
ifeq ($(COMPILER),gnu)
  CC = gcc
endif

# intel Compiler
ifeq ($(COMPILER),intel)
  CC = icc
endif

# BG/Q gcc Cross-Compiler
ifeq ($(MACHINE),bluegene)
  CC = mpicc
endif

# Standard Flags
CFLAGS := -std=gnu99

# Linker Flags
LDFLAGS = -lm

# Debug Flags
ifeq ($(DEBUG),yes)
  CFLAGS += -g
endif

# Profiling Flags
ifeq ($(PROFILE),yes)
  CFLAGS += -pg
endif

# Optimization Flags
ifeq ($(OPTIMIZE),yes)
  CFLAGS += -O3
endif

# PAPI source (you may need to provide -I and -L pointing
# to PAPI depending on your installation
ifeq ($(PAPI),yes)
  source += papi.c
  CFLAGS += -DPAPI
  #CFLAGS += -I/soft/apps/packages/papi/papi-5.1.1/include
  #LDFLAGS += -L/soft/apps/packages/papi/papi-5.1.1/lib -lpapi
  LDFLAGS += -lpapi
endif

# MPI
ifeq ($(MPI),yes)
  CC = mpicc
  CFLAGS += -DMPI
endif

# OpenMP
ifeq ($(OPENMP),yes)
	CFLAGS += -fopenmp
endif


#===============================================================================
# Targets to Build
#===============================================================================

$(program): $(obj) SimpleMOC_header.h
	$(CC) $(CFLAGS) $(obj) -o $@ $(LDFLAGS)

%.o: %.c
	$(CC) $(CFLAGS) -c $< -o $@

clean:
	rm -rf $(program) $(obj)

edit:
	vim -p $(source) SimpleMOC_header.h

run:
	./$(program)
