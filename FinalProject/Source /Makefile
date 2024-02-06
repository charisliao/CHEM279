# Problem Set 1 Makefile for source codes
# Date: 9/10/2023
# Creator: Charis Liao
# University of California, Berkeley
# CHEM 279

# GNU C++ Compiler
CPP = g++  # C++ Compiler

# Update CPPFLAGS to include -std=c++11
CPPFLAGS = -g -c -std=c++11 -pg

# Update LDFLAGS to include -fopenmp
LDFLAGS = -pg -fopenmp


# The directory where the include files needed to create the library objects are
INC_DIR = ../Include

# The directory where all object files are going to be moved after creation
LIB_DIR = ../Libdir

# The directory path for armadillo library
ARMADILLO_DIR = /Users/charisliao/Desktop/armadillo-12.6.2

# The next line contains the list of object files created by this Makefile.
INTEGRATION_OBJS = AO.o utils.o cndo.o
INTEGRATION_LIB = AO.a

AO.a: $(INTEGRATION_OBJS)
	# create a library file from the object files
	ar ru $(INTEGRATION_LIB) $(INTEGRATION_OBJS)
	# check the library
	ranlib $(INTEGRATION_LIB)
	mv $(INTEGRATION_LIB) $(LIB_DIR)/.
	rm $(INTEGRATION_OBJS)

AO.o: $(INC_DIR)/AO.h $(INC_DIR)/utils.h
	$(CPP) $(CPPFLAGS) AO.cpp -I$(INC_DIR) -I$(ARMADILLO_DIR)/include

utils.o: $(INC_DIR)/utils.h
	$(CPP) $(CPPFLAGS) utils.cpp -I$(INC_DIR) -I$(ARMADILLO_DIR)/include

cndo.o: $(INC_DIR)/CNDO.h $(INC_DIR)/AO.h $(INC_DIR)/utils.h
	clang -Xclang -fopenmp -L/opt/homebrew/opt/libomp/lib -I/opt/homebrew/opt/libomp/include -lomp -o cndo.o -c cndo.cpp -std=c++11 -I$(INC_DIR) -I$(ARMADILLO_DIR)/include

hamiltonian.o: $(INC_DIR)/hamiltonian.h $(INC_DIR)/AO.h $(INC_DIR)/utils.h
	$(CPP) $(CPPFLAGS) hamiltonian.cpp -I$(INC_DIR) -I$(ARMADILLO_DIR)/include

all: $(INTEGRATION_LIB)

cleanall:
	cd $(LIB_DIR); rm $(INTEGRATION_LIB)
