/**
 * @file test1.cpp
 * @author Charis Liao 
 * @brief Test for problem 2 - Finite Difference Versus Analytical Forces 
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

    molecule mol("test_case_1.txt");
    molecule mol2("test_case_2.txt");
    
    LennardJones LJ_forces_1(mol);
    LennardJones LJ_forces_2(mol2);
    

    // Since we're dealing with all gold here, the value for
    // epsilon and sigma are the same for all atoms.

    double epsilon_i = 5.29; // kcal/mol
    double sigma_i = 2.951; // Angstrom
    double epsilon_j = 5.29; // kcal/mol
    double sigma_j = 2.951; // Angstrom
    
    
    // Calculate the analytical forces of the molecule.
    mat LJ_forces1 = LJ_forces_1.calculate_LJ_forces(mol, epsilon_i, epsilon_j, sigma_i, sigma_j);
    mat LJ_forces2 = LJ_forces_2.calculate_LJ_forces(mol2, epsilon_i, epsilon_j, sigma_i, sigma_j);
    
    // Calulate the norm of the analytical forces
    double norm1 = LJ_forces_1.calculate_force_norm(LJ_forces1);
    double norm2 = LJ_forces_2.calculate_force_norm(LJ_forces2);

    // Print analytical forces calculated and force norm
    cout << "The analytical forces calculated for test case 1 are: " << endl;
    LJ_forces1.print();
    cout << "The force norm for test case 1 is: " << norm1 << endl;



    cout << "The analytical forces calculated for test case 2 are: " << endl;
    LJ_forces2.print();
    cout << "The force norm for test case 2 is: " << norm2 << endl;



    // Approximate the analytical forces using forward difference approximation 
    double h = 0.01; 
    mat LJ_forces_forward1 = LJ_forces_1.calculate_LJforces_forward(mol, h, epsilon_i, epsilon_j, sigma_i, sigma_j);
    mat LJ_forces_forward2 = LJ_forces_2.calculate_LJforces_forward(mol2, h, epsilon_i, epsilon_j, sigma_i, sigma_j);

    // Print forward difference approximation forces calculated
    cout << "The forward difference approximation forces calculated for test case 1 are: " << endl;
    LJ_forces_forward1.print();
    cout << "The forward difference approximation forces calculated for test case 2 are: " << endl;
    LJ_forces_forward2.print();

    // Approximate the analytical forces using central difference approximation
    mat LJ_forces_central1 = LJ_forces_1.calculate_LJforces_central(mol, h, epsilon_i, epsilon_j, sigma_i, sigma_j);
    mat LJ_forces_central2 = LJ_forces_2.calculate_LJforces_central(mol2, h, epsilon_i, epsilon_j, sigma_i, sigma_j);

    // Print central difference approximation forces calculated
    cout << "The central difference approximation forces calculated for test case 1 are: " << endl;
    LJ_forces_central1.print();
    cout << "The central difference approximation forces calculated for test case 2 are: " << endl;
    LJ_forces_central2.print();

    // write all the print statements to a file 
    ofstream problem2;
    problem2.open("problem_2.txt");
   
    problem2 << "The analytical forces calculated for test case 1 are: " << endl;
    LJ_forces1.print(problem2);
    problem2 << "The force norm for test case 1 is: " << norm1 << endl;
    problem2 << "The analytical forces calculated for test case 2 are: " << endl;
    LJ_forces2.print(problem2);
    problem2 << "The force norm for test case 2 is: " << norm2 << endl;
    problem2 << "The forward difference approximation forces calculated for test case 1 are: " << endl;
    LJ_forces_forward1.print(problem2);
    problem2 << "The forward difference approximation forces calculated for test case 2 are: " << endl;
    LJ_forces_forward2.print(problem2);
    problem2 << "The central difference approximation forces calculated for test case 1 are: " << endl;
    LJ_forces_central1.print(problem2);
    problem2 << "The central difference approximation forces calculated for test case 2 are: " << endl;
    LJ_forces_central2.print(problem2);

    problem2.close();


    return 0;
}