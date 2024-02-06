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
INTEGRATION_OBJS = gaussian.o numerical-integration.o shell.o
INTEGRATION_LIB = numerical-integration.a

numerical-integration.a: $(INTEGRATION_OBJS)
# create a library file from the object files
	ar ru $(INTEGRATION_LIB) $(INTEGRATION_OBJS)
# check the library
	ranlib $(INTEGRATION_LIB)
	mv $(INTEGRATION_LIB) $(LIB_DIR)/. 
	rm $(INTEGRATION_OBJS)


gaussian.o: $(INC_DIR)/gaussian.h
	$(CPP) $(CPPFLAGS) gaussian.cpp -I$(INC_DIR) -I$(ARMADILLO_DIR)/include

numerical-integration.o: $(INC_DIR)/numerical-integration.h $(INC_DIR)/gaussian.h
	$(CPP) $(CPPFLAGS) numerical-integration.cpp -I$(INC_DIR) -I$(ARMADILLO_DIR)/include

shell.o: $(INC_DIR)/shell.h
	$(CPP) $(CPPFLAGS) shell.cpp -I$(INC_DIR) -I$(ARMADILLO_DIR)/include


all: $(INTEGRATION_LIB)

cleanall: 
	cd $(LIB_DIR); rm $(INTEGRATION_LIB)