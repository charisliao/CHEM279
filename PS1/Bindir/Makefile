# Problem Set 1 Makefile for executables generated from tests  
# Date: 9/10/2023 
# Creator: Charis Liao 
# University of California, Berkeley 
# CHEM 279 

# GNU C++ Compiler
CPP         = g++      # C++ Compuler

# compiler flags -g debug, -c creates a library object
CPPFLAGS    = -g -c    

# The next line contains the list of executables run by this Makefile.
EXECUTABLES = run_test_1 run_test_2 run_test_2.1 run_test_3

EXECUTABLES_NAME = test_1 test_2 test_2.1 test_3

UTILS_DIR = ../Utils/Plots/Data

# UTILS_DIR = ../Utils/Plots/Data

# The next line will run the battery_model executable
run_test_1:
	./test_1 

run_test_2:
	./test_2 

run_test_2.1:
	./test_2.1
	mv truncation_error.txt truncation_error_2.txt $(UTILS_DIR)/.

run_test_3:
	./test_3
	
# The next line will run all the executables
all: $(EXECUTABLES)

# The next line will remove all the executables 
cleanall:
	rm -f problem_2.txt 
	rm -f $(EXECUTABLES_NAME)
	




