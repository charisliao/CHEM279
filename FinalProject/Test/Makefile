# Problem Set 1 Makefile for tests
# Creator: Charis Liao
# Date: 9/10/2023

# This makefile creates the text executable files that test the
# functionalities in molecular and lennerd_jones files.

# GNU C++ Compiler
CPP = g++      # C++ Compiler

# Compiler flags -g for debug, -o for creating an executable file
CPPFLAGS = -g -std=c++11 -pg -o  # Added -pg for profiling

# The directory where the include files needed to create the library objects are
INC_DIR = ../Include

# The directory where all object files are going to be moved after creation
LIB_DIR = ../Libdir

# The directory where all executable files are going to be moved after creation
BIN_DIR = ../Bindir

# The directory path for Armadillo library
ARMADILLO_DIR = /Users/charisliao/Desktop/armadillo-12.6.2

# Additional linker flags
LDFLAGS = -L/opt/homebrew/opt/veclibfort/lib -lveclibfort

# The next line contains the list of executable files that are created by the Makefile
EXECUTABLES = test_h test_h2 test_c2h4 test_hf+ test_n2 test_ot test_modified_runtime generate_molecules


# The object file for insertion sort. The file is saved in Libdir.

test_h: $(INC_DIR)/AO.h $(INC_DIR)/utils.h $(INC_DIR)/CNDO.h
	$(CPP) $(CPPFLAGS) test_h test_h.cpp \
	-I$(INC_DIR) -I$(ARMADILLO_DIR)/include \
	$(LIB_DIR)/AO.a $(LDFLAGS)
	mv test_h $(BIN_DIR)

test_h2: $(INC_DIR)/AO.h $(INC_DIR)/utils.h $(INC_DIR)/CNDO.h
	$(CPP) $(CPPFLAGS) test_h2 test_h2.cpp \
	-I$(INC_DIR) -I$(ARMADILLO_DIR)/include \
	$(LIB_DIR)/AO.a $(LDFLAGS)
	mv test_h2 $(BIN_DIR)

test_hf+: $(INC_DIR)/AO.h $(INC_DIR)/utils.h $(INC_DIR)/CNDO.h
	$(CPP) $(CPPFLAGS) test_hf+ test_hf+.cpp \
	-I$(INC_DIR) -I$(ARMADILLO_DIR)/include \
	$(LIB_DIR)/AO.a $(LDFLAGS)
	mv test_hf+ $(BIN_DIR)

test_c2h4: $(INC_DIR)/AO.h $(INC_DIR)/utils.h $(INC_DIR)/CNDO.h
	$(CPP) $(CPPFLAGS) test_c2h4 test_c2h4.cpp \
	-I$(INC_DIR) -I$(ARMADILLO_DIR)/include \
	$(LIB_DIR)/AO.a $(LDFLAGS)
	mv test_c2h4 $(BIN_DIR)

test_n2: $(INC_DIR)/AO.h $(INC_DIR)/utils.h $(INC_DIR)/CNDO.h
	$(CPP) $(CPPFLAGS) test_n2 test_n2.cpp \
	-I$(INC_DIR) -I$(ARMADILLO_DIR)/include \
	$(LIB_DIR)/AO.a $(LDFLAGS)
	mv test_n2 $(BIN_DIR)


test_ot: $(INC_DIR)/AO.h $(INC_DIR)/utils.h $(INC_DIR)/CNDO.h
	$(CPP) $(CPPFLAGS) test_ot test_ot.cpp \
	-I$(INC_DIR) -I$(ARMADILLO_DIR)/include \
	$(LIB_DIR)/AO.a $(LDFLAGS)
	mv test_ot $(BIN_DIR)

test_modified_runtime: $(INC_DIR)/AO.h $(INC_DIR)/utils.h $(INC_DIR)/CNDO.h
	$(CPP) $(CPPFLAGS) test_modified_runtime test_modified_runtime.cpp \
	-I$(INC_DIR) -I$(ARMADILLO_DIR)/include \
	$(LIB_DIR)/AO.a $(LDFLAGS)
	mv test_modified_runtime $(BIN_DIR)

generate_molecules: $(INC_DIR)/AO.h $(INC_DIR)/utils.h $(INC_DIR)/CNDO.h
	$(CPP) $(CPPFLAGS) generate_molecules generate_molecules.cpp \
	-I$(INC_DIR) -I$(ARMADILLO_DIR)/include \
	$(LIB_DIR)/AO.a $(LDFLAGS)
	mv generate_molecules $(BIN_DIR)

all: $(EXECUTABLES)

cleanall:
	cd $(BIN_DIR); rm $(EXECUTABLES)
