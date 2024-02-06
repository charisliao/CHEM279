/**
 * @file molecule.cpp
 * @author your name (you@domain.com)
 * @brief Create a class called molecule 
 * @version 0.1
 * @date 2023-09-09
 * 
 * @copyright Copyright (c) 2023
 * 
 */

#include "molecule.h"
#include <iostream>
#include <fstream>


using namespace std;

// Constructor
molecule::molecule(const char *filename){

    ifstream infile(filename);
    if (!infile.is_open()) {
        cerr << "Error: Unable to open file " << filename << endl;
    // Add additional error handling if needed
    } else {
    // Proceed with reading the file

        infile >> natoms;
        atomic_numbers.resize(natoms); 
        coord.resize(natoms, 3);  // this is an armadiilo matrix

        for (int i = 0; i < natoms; i++){
            infile >> atomic_numbers[i] >> coord(i, 0) >> coord(i, 1) >> coord(i, 2);
        }
    }

    // check if the file has any atoms other than gold or other unacceptable inputs
    for (int atom = 0; atom < natoms; atom++){
        if (atomic_numbers[atom] != 79){
            cerr << "Error: Atomic number " << atomic_numbers[atom] << " is not gold." << endl;
        }
    }

     


    infile.close();
}