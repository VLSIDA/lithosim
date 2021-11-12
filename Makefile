# Copyright (c) 2003 - 2007 ACM/SIGDA
#
# Written by Florian Krohm (fkrohm@us.ibm.com)
#
# This Makefile will compile all C (*.c) and C++ (*.C) sources in 
# this directory and link them into an executable specified by variable TARGET.
# In order not to have to model the exact dependencies on header files this
# Makefile assumes that all C/C++ files depend on all header files.
# Bison and flex input files, if any, are handled as expected.
#
# make debug   - compile and link to produce debuggable executable
# make opt     - compile and link to produce optimized executable
# make clean   - remove all generated files
# make test    - run all testcases
# make submit  - copy relevant files to solution directory
#
# You may change the value of SUBMIT_FILES to fit your needs
# You *must not* modify TARGET or SUBMIT_DIR.
#

TARGET = lithosim

SUBMIT_DIR   = ../solution
SUBMIT_FILES = $(wildcard *.[Cchyl]) Makefile lib.defs

#
# Tools used
#
CC    = gcc -std=c99 
CXX   = g++ 
FLEX  = flex
BISON = bison
COMPARE_RESULT = diff-float

#
# Extra libraries
#
LIBS = -lm

#
# Define LD_RUN_PATH such that ld.so will find the shared objcets
#export OA_HOME=/mada/software/oa2.2.6
#export LD_RUN_PATH=$(OA_HOME)/lib
#export LD_LIBRARY_PATH = $(OA_HOME)/lib/linux_rhel30_64/opt

#
# Bison, Flex input if any 
# (Only one grammar/scanner spec at this point)
#
Y_FILE := $(wildcard *.y)
Y_SRCS := $(subst .y,.tab.c,$(Y_FILE))
Y_HDRS := $(subst .y,.tab.h,$(Y_FILE))
L_FILE := $(wildcard *.l)
L_SRCS := $(subst .l,.yy.c,$(L_FILE))

#
# Assemble sources, objects, and headers
#
C_SRCS  := $(wildcard *.c)
C_SRCS  := $(subst $(L_SRCS),,$(C_SRCS))
C_SRCS  := $(subst $(Y_SRCS),,$(C_SRCS))
C_SRCS  += $(L_SRCS)
C_DOBJS := $(C_SRCS:.c=.o)
C_SRCS  += $(Y_SRCS)
C_OBJS  := $(C_SRCS:.c=.o)
L_OBJS  := $(L_SRCS:.c=.o)

CXX_SRCS  := $(wildcard *.C)
CXX_OBJS  := $(CXX_SRCS:.C=.o)

HDRS   := $(wildcard *.h)
HDRS   := $(subst $(Y_HDRS),,$(HDRS))
HDRS   += $(Y_HDRS)

# Choose suitable commandline flags 
#
ifeq "$(MAKECMDGOALS)" "opt"
#CFLAGS   = -O3 -march=athlon64
CFLAGS   = -O3 -mssse3
CXXFLAGS = -O3 -mssse3
else
CFLAGS   = -g -W -Wall -pedantic 
CXXFLAGS = -g 
#-W -Wall -pedantic 
#-pg
#-W -Wall -pedantic
endif

CFLAGS   +=  -I. 
CXXFLAGS +=  -I.

.PHONY:	clean test debug opt submit

.SECONDARY: $(Y_SRCS) $(L_SRCS)

debug opt: $(TARGET)

$(TARGET):  $(C_OBJS) $(CXX_OBJS) 
	$(CXX) $(CXXFLAGS) $(LIBS) -o $(TARGET) $(CXX_OBJS) 

%.yy.c:%.l
	$(FLEX) -l -o$(notdir $*).yy.c $<

%.tab.c %.tab.h:%.y
	$(BISON) -d --output=$*.tab.c $<

%.o:%.c
	$(CC) -c $(CFLAGS) $<

%.o:%.C
	$(CXX) -c $(CXXFLAGS) $< 

$(C_DOBJS) $(CXX_OBJS): $(HDRS)

#
# Flex generates code that causes harmless warnings; suppress those
#
$(L_OBJS): $(L_SRCS) $(HDRS)
	$(CC) -c $(CFLAGS) -Wno-unused-function -Wno-unused-label \
                 -Wno-implicit-function-declaration $<

bar: $(TARGET)
	./$(TARGET) 8 tests/2_thick_bar7.pbm results/bar-aerial.pnm results/bar-contours.pnm

baropc: $(TARGET)
	./$(TARGET) 8 tests/2_thick_bar7.pbm results/bar-opc-aerial.pnm results/bar-opc-contours.pnm results/bar-opc.pnm

via: $(TARGET)
	./$(TARGET) 8 tests/via.pbm results/via-aerial.pnm results/via-contours.pnm

viaopc: $(TARGET)
	./$(TARGET) 8 tests/via.pbm results/via-opc-aerial.pnm results/via-opc-contours.pnm results/via-opc.pnm

tiny: $(TARGET)
#	./$(TARGET) 16 tests/tiny.pbm tiny-aerial-180.pnm tiny-contours-180.pnm
#	./$(TARGET) 12 tests/tiny.pbm tiny-aerial-130.pnm tiny-contours-130.pnm
	./$(TARGET) 8 tests/tiny.pbm results/tiny-aerial-90.pnm results/tiny-contours-90.pnm
#	./$(TARGET) 6 tests/tiny.pbm tiny-aerial-65.pnm tiny-contours-65.pnm
#	./$(TARGET) 4 tests/tiny.pbm tiny-aerial-45.pnm tiny-contours-45.pnm

tinyopc: $(TARGET)
	./$(TARGET) 8 tests/tiny.pbm results/tiny-opc-aerial-90.pnm results/tiny-opc-contours-90.pnm results/tiny-opc-90.pnm

small: $(TARGET)
	./$(TARGET) 8 tests/small.pbm results/small-aerial.pnm results/small-contours.pnm

smallopc: $(TARGET)
	./$(TARGET) 8 tests/small.pbm results/small-opc-aerial.pnm results/small-opc-contours.pnm results/bar-opc.pnm

#medium: $(TARGET)
#	./$(TARGET) tests/medium.pbm results/medium-aerial.pnm results/medium-contours.pnm


clean : 
	rm -f *.o $(TARGET) $(Y_HDRS) $(Y_SRCS) $(L_SRCS) results/*

submit:
	$(HOME)/bin/submit $(SUBMIT_DIR) $(SUBMIT_FILES)
