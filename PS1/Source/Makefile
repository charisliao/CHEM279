# Problem Set 1 Makefile for source codes 
# Date: 9/10/2023 
# Creator: Charis Liao 
# University of California, Berkeley 
# CHEM 279 


# GNU C++ Compiler
CPP         = g++      # C++ Compuler

# compiler flags -g debug, -c creates a library object
CPPFLAGS    = -g -c -std=c++11

# The directory where the include files needed to create the library objects are
INC_DIR = ../Include

# The directory where all object files are going to moved afer creation
LIB_DIR     = ../Libdir

# The directory path for armadillo library
ARMADILLO_DIR = /Users/charisliao/Desktop/armadillo-12.6.2

# The next line contains the list of object files created by this Makefile.
ENERGY_OBJS = lennard_jones.o molecule.o steepest_descent.o

ENERGY_LIB = lennard_jones.a

lennard_jones.a: $(ENERGY_OBJS)
# create a library file from the object files
	ar ru $(ENERGY_LIB) $(ENERGY_OBJS)
# check the library
	ranlib $(ENERGY_LIB)
	mv $(ENERGY_LIB) $(LIB_DIR)/.
	rm $(ENERGY_OBJS)


lennard_jones.o: $(INC_DIR)/lennard_jones.h
	$(CPP) $(CPPFLAGS) lennard_jones.cpp -I$(INC_DIR) -I$(ARMADILLO_DIR)/include

molecule.o: $(INC_DIR)/molecule.h
	$(CPP) $(CPPFLAGS) molecule.cpp -I$(INC_DIR) -I$(ARMADILLO_DIR)/include

steepest_descent.o: $(INC_DIR)/steepest_descent.h
	$(CPP) $(CPPFLAGS) steepest_descent.cpp -I$(INC_DIR) -I$(ARMADILLO_DIR)/include


all: $(ENERGY_LIB)

cleanall: 
	cd $(LIB_DIR); rm $(ENERGY_LIB)