# Problem Set 1 Makefile for tests 
# Creator: Charis Liao 
# Date: 9/10/2023 

# This makefi;e creates the text executable files that test the 
# functionalities in molecular and lennerd_jones files. 

# GNU C++ Compiler
CPP         = g++      # C++ Compuler

# compiler flags -g debug, -o creates an executable file
CPPFLAGS    = -g -std=c++11 -o

# The directory where the include files needed to create the library objects are
INC_DIR = ../Include

# The directory where all object files are going to be moved afer creation
LIB_DIR     = ../Libdir

# The directory where all executable files are going to be moved afer creation
BIN_DIR     = ../Bindir

# The directory path for armadillo library
ARMADILLO_DIR = /Users/charisliao/Desktop/armadillo-12.6.2

# The next line contains the list of executable files that are created by the Makefile
EXECUTABLES = test_2.1 test_2.2

#The object file for insertion sort. The file is saved in Libdir.

test_2.1: $(INC_DIR)/gaussian.h $(INC_DIR)/numerical-integration.h
	$(CPP) $(CPPFLAGS) test_2.1 test_2.1.cpp \
	-I$(INC_DIR) -I$(ARMADILLO_DIR)/include \
	$(LIB_DIR)/numerical-integration.a
	mv test_2.1 $(BIN_DIR)

test_2.2: $(INC_DIR)/shell.h
	$(CPP) $(CPPFLAGS) test_2.2 test_2.2.cpp \
	-I$(INC_DIR) -I$(ARMADILLO_DIR)/include \
	$(LIB_DIR)/numerical-integration.a
	mv test_2.2 $(BIN_DIR)



all: $(EXECUTABLES)

cleanall:
	cd $(BIN_DIR); rm $(EXECUTABLES)