/**
 * @file test1.cpp
 * @author Charis Liao 
 * @brief Test for problem 1 - The Lennard-Jones energy of a set of atoms
 * @version 0.1
 * @date 2023-09-10
 * 
 * @copyright Copyright (c) 2023
 * 
 */


#include "molecule.h"
#include "lennard_jones.h"
#include "utils.h"
#include <iostream>
#include <armadillo>
#include <vector>
#include <string>
#include <cassert>
#include <iomanip>
#include <fstream>
#include <cmath>


using namespace std;
using namespace arma;

int main() {
    // Read in the file into our attribute
    molecule mol("test_case_1.txt");
    molecule mol2("test_case_2.txt");
    LennardJones LJ_energy(mol);
    LennardJones LJ_energy_2(mol2);

    // Since we're dealing with all gold here, the value for
    // epsilon and sigma are the same for all atoms.

    double epsilon_i = 5.29; // kcal/mol
    double sigma_i = 2.951; // Angstrom
    double epsilon_j = 5.29; // kcal/mol
    double sigma_j = 2.951; // Angstrom

    // Calculate the total Lennard-Jones potential of the molecule.
    double total_LJ_1 = LJ_energy.calculate_total_LJ(mol, epsilon_i, epsilon_j, sigma_i, sigma_j);
    double total_LJ_2 = LJ_energy_2.calculate_total_LJ(mol2, epsilon_i, epsilon_j, sigma_i, sigma_j);
    // Print the total Lennard-Jones potential of the molecule.
    cout << "The total Lennard-Jones potential of the molecule for test case 1 is " << total_LJ_1 << " kcal/mol." << endl;
    cout << "The total Lennard-Jones potential of the molecule for test case 2 is " << total_LJ_2 << " kcal/mol." << endl;
    return 0;   
};