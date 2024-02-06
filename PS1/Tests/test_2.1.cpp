/**
 * @file test1.cpp
 * @author Charis Liao 
 * @brief Test for problem 2.1 - Finite Difference Versus Analytical Forces 
 * We're going to investivate the truncation error of the forward and 
 * central difference expressions as a function of h, trying h = 0.1, 0.01, 0.001, 0.0001.
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

        vector<double> step_sizes = {0.1, 0.01, 0.001, 0.0001};
        vector<double> truncation_error_forward;
        vector<double> truncation_error_central;
        vector<double> truncation_error_forward_2;
        vector<double> truncation_error_central_2;

        // Truncation error for test case 1
        for (double step : step_sizes) {
            // Calculate forward and central difference forces 
            mat LJ_forces_forward = LJ_forces_1.calculate_LJforces_forward(mol, step, epsilon_i, epsilon_j, sigma_i, sigma_j);
            mat LJ_forces_central = LJ_forces_1.calculate_LJforces_central(mol, step, epsilon_i, epsilon_j, sigma_i, sigma_j);

            // Calculate the truncation error 
            mat TruncationError_central = abs(LJ_forces_central - LJ_forces1);
            mat TruncationError_forward = abs(LJ_forces_forward - LJ_forces1);

            // Calculate the norm of the truncation error 
            double norm_truncation_error_central = LJ_forces_1.calculate_force_norm(TruncationError_central);
            double norm_truncation_error_forward = LJ_forces_1.calculate_force_norm(TruncationError_forward);

            // Store the norm of the truncation error in a vector
            truncation_error_central.push_back(norm_truncation_error_central);
            truncation_error_forward.push_back(norm_truncation_error_forward);
        }

        // Write the truncation error to a file
        ofstream truncation_error;
        truncation_error.open("truncation_error.txt");
        truncation_error <<"0.1, 0.01, 0.001, 0.0001" << endl;
        truncation_error << "Truncation error for forward difference approximation for test case 1: " << endl;
        for (double error : truncation_error_forward) {
            truncation_error << error << ", ";
        }
        truncation_error << endl;
        truncation_error << "Truncation error for central difference approximation for test case 1: " << endl;
        for (double error : truncation_error_central) {
            truncation_error << error << ", ";
        }
        truncation_error << endl;
        truncation_error.close();


        // Truncation error for test case 2
        for (double step : step_sizes) {
            // Calculate forward and central difference forces 
            mat LJ_forces_forward = LJ_forces_2.calculate_LJforces_forward(mol2, step, epsilon_i, epsilon_j, sigma_i, sigma_j);
            mat LJ_forces_central = LJ_forces_2.calculate_LJforces_central(mol2, step, epsilon_i, epsilon_j, sigma_i, sigma_j);

            // Calculate the truncation error 
            mat TruncationError_central = abs(LJ_forces_central - LJ_forces2);
            mat TruncationError_forward = abs(LJ_forces_forward - LJ_forces2);

            // Calculate the norm of the truncation error 
            double norm_truncation_error_central = LJ_forces_2.calculate_force_norm(TruncationError_central);
            double norm_truncation_error_forward = LJ_forces_2.calculate_force_norm(TruncationError_forward);

            // Store the norm of the truncation error in a vector
            truncation_error_central_2.push_back(norm_truncation_error_central);
            truncation_error_forward_2.push_back(norm_truncation_error_forward);
        }

        // Write the truncation error 2 to a file
        ofstream truncation_error_2;
        truncation_error_2.open("truncation_error_2.txt");
        truncation_error_2 <<"0.1, 0.01, 0.001, 0.0001" << endl;
        truncation_error_2 << "Truncation error for forward difference approximation for test case 2: " << endl;
        for (double error : truncation_error_forward_2) {
            truncation_error_2 << error << ", ";
        }
        truncation_error_2 << endl;
        truncation_error_2 << "Truncation error for central difference approximation for test case 2: " << endl;
        for (double error : truncation_error_central_2) {
            truncation_error_2 << error << ", ";
        }
        truncation_error_2 << endl;
        truncation_error_2.close();

}
