/**
 * @file steepest_descend.cpp
 * @author Charis Liao 
 * @brief Implement the steepest descent algorithm with line search to find the minimum of the Lennard-Jones potential
 * @version 0.1
 * @date 2023-09-11
 * 
 * @copyright Copyright (c) 2023
 * 
 */

#include "molecule.h"
#include "lennard_jones.h"
#include "steepest_descent.h"
#include <iostream>
#include <armadillo>
#include <vector>

using namespace std;
using namespace arma;

// Constructor
SteepestDescent::SteepestDescent(const molecule& mol) : mol_(mol) {
    cout << "Steepest descent initialized." << endl;
}

void SteepestDescent::steepest_descent_optimizer(molecule mol, double stepsize, double tol, double epsilon_i, double epsilon_j, double sigma_i, double sigma_j) {
    int count = 0; // count the number of iterations
    LennardJones LJ_forces(mol);
    mat gradient = -1 * (LJ_forces.calculate_LJ_forces(mol, epsilon_i, epsilon_j, sigma_i, sigma_j));
    mat new_geom = mol.coord; 

    while (LJ_forces.calculate_force_norm(gradient) > tol && count < 100000)  {
        molecule old_mol = mol;
        new_geom = mol.coord + stepsize * gradient;
        mol.coord = new_geom;
        gradient = -1 * (LJ_forces.calculate_LJ_forces(mol, epsilon_i, epsilon_j, sigma_i, sigma_j));

        cout << "The potential energy at iteration " << count << " is " << LJ_forces.calculate_total_LJ(mol, epsilon_i, epsilon_j, sigma_i, sigma_j) << " kcal/mol." << endl;
        if (LJ_forces.calculate_total_LJ(mol, epsilon_i, epsilon_j, sigma_i, sigma_j) > LJ_forces.calculate_total_LJ(old_mol, epsilon_i, epsilon_j, sigma_i, sigma_j)) {
            stepsize = 1.2 * stepsize;
        } else {
            stepsize = 0.5 * stepsize;
        }
        count++;
    }
    cout << "The minimum of the Lennard-Jones potential is " << LJ_forces.calculate_total_LJ(mol, epsilon_i, epsilon_j, sigma_i, sigma_j) << " kcal/mol." << endl;
    cout << "It takes " << count << " iterations to reach the minimum." << endl;

}