/**
 * @file lennard-jones.cpp
 * @author Charis Liao
 * @brief This file implements the Lennard-Jones potential of a set of atoms.
 * @version 0.1
 * @date 2023-09-08
 * 
 * @copyright Copyright (c) 2023
 * 
 */


#include <iostream>
#include <cmath>

#include "lennard_jones.h"
#include "molecule.h"

using namespace std;

// Constructor
LennardJones::LennardJones(const molecule& mol) : mol_(mol) {
    cout << "Lennard-Jones potential initialized." << endl;
}


// Calculate the epsilon_ij between atom i and atom j for the Lennard-Jones potential.
double LennardJones::calculate_epsilon_ij(double epsilon_i, double epsilon_j) {
    return sqrt(epsilon_i * epsilon_j);
}


// Calculate the sigma_ij between atom i and atom j for the Lennard-Jones potential.
double LennardJones::calculate_sigma_ij(double sigma_i, double sigma_j) {
    return sqrt(sigma_i * sigma_j);
}

/**
 * @brief Calculate the Lennard-Jones potential between atom i and atom j.
 * 
 * @param epsilon_ij: geometric mean of epsilon_i and epsilon_j
 * @param sigma_ij: geometric mean of sigma_i and sigma_j
 * @param r_ij: distance between atom i and atom j
 * @return double 
 */
double LennardJones::calculate_LJ(double epsilon_ij, double sigma_ij, double r_ij) {
    return epsilon_ij * (pow(sigma_ij / r_ij, 12) - 2 * pow(sigma_ij / r_ij, 6));
}

// Calculate the distance between atom i and atom j.
double LennardJones::calculate_distance(double x_i, double y_i, double z_i, double x_j, double y_j, double z_j) {
    return sqrt(pow(x_i - x_j, 2) + pow(y_i - y_j, 2) + pow(z_i - z_j, 2));
}

// Calculate the total Lennard-Jones potential of a molecule.
double LennardJones::calculate_total_LJ(molecule mol, double epsilon_i, double epsilon_j, double sigma_i, double sigma_j) {
    double total_LJ = 0.0;
    double epsilon_ij, sigma_ij, r_ij;
    int natoms = mol.natoms;
    arma::mat geom = mol.coord;
    vector<int> atomic_numbers = mol.atomic_numbers;

    for (int i = 0; i < natoms; i++) {
        for (int j = i + 1; j < natoms; j++) {
            epsilon_ij = calculate_epsilon_ij(epsilon_i, epsilon_j);
            sigma_ij = calculate_sigma_ij(sigma_i, sigma_j);
            r_ij = calculate_distance(geom(i, 0), geom(i, 1), geom(i, 2), geom(j, 0), geom(j, 1), geom(j, 2));
            total_LJ += calculate_LJ(epsilon_ij, sigma_ij, r_ij);
        }
    }

    return total_LJ;
}

// Calculate the derivative of the Lennard-Jones potential with respect to position (x or y or z).
double LennardJones::calculate_LJ_force(double epsilon_ij, double sigma_ij, double r_ij, double x_i, double x_j) {
    double partial_derivative = (x_i - x_j) / r_ij;
    return -1 * (epsilon_ij * (-12 * pow(sigma_ij, 12) * pow(1 / r_ij, 13) + 12 * pow(sigma_ij, 6) * pow(1 / r_ij, 7)) * partial_derivative);
}

// Calculate the analytical forces of a molecule 
arma::mat LennardJones::calculate_LJ_forces(molecule mol, double epsilon_i, double epsilon_j, double sigma_i, double sigma_j) {
    double LJ_force = 0.0; 
    double epsilon_ij, sigma_ij, r_ij, force_norm;
    int natoms = mol.natoms;
    arma::mat geom = mol.coord;
    arma::mat forces = arma::zeros(natoms, 3);
    

    for (int i = 0; i < natoms; i++) {
        for (int j = i + 1; j < natoms; j++) {
            
            epsilon_ij = calculate_epsilon_ij(epsilon_i, epsilon_i);
            sigma_ij = calculate_sigma_ij(sigma_i, sigma_i);
            r_ij = calculate_distance(mol.coord(i, 0), mol.coord(i, 1), mol.coord(i, 2), mol.coord(j, 0), mol.coord(j, 1), mol.coord(j, 2));

            // Calculate the LJ force for each atom
            // arma::vec3 rij_vector = mol.coord.row(i) - mol.coord.row(j);
            arma::vec3 rij_vector = {mol.coord(i, 0) - mol.coord(j, 0), mol.coord(i, 1) - mol.coord(j, 1), mol.coord(i, 2) - mol.coord(j, 2)};

            for (int dim = 0; dim < 3; dim++) {
                LJ_force = calculate_LJ_force(epsilon_ij, sigma_ij, r_ij, mol.coord(i, dim), mol.coord(j, dim));
                forces(i, dim) += LJ_force;
                forces(j, dim) -= LJ_force;
            }
          


            // update forces on atom i and atom j
            // forces.row(i) += LJ_force * (rij_vector / r_ij).t();
            // forces.row(j) -= LJ_force * (rij_vector / r_ij).t();

            // forces.row(i) += (LJ_force / r_ij * rij_vector).t();
            // forces.row(j) -= (LJ_force / r_ij * rij_vector).t();
        }
    }
    
    


    return forces;
}

double LennardJones::calculate_force_norm(arma::mat forces) {
    double force_norm = 0.0;
    for (int i = 0; i < forces.n_rows; i++) {
        for (int dim = 0; dim < 3; dim++) {
            force_norm += pow(forces(i, dim), 2);
        }
    }
    return sqrt(force_norm);
}


// Approximate the analytical forces of a molecule using forward difference approximation. 
arma::mat LennardJones::calculate_LJforces_forward(molecule mol, double h, double epsilon_i, double epsilon_j, double sigma_i, double sigma_j) {
    arma::mat forces = arma::zeros(mol.natoms, 3); // initialize the forces matrix
    arma::mat geom = mol.coord; 
    int natoms = mol.natoms;

    for (int i = 0; i < natoms; i++) {
        for (int dim = 0; dim < 3; dim++) {
            double original_coord = mol.coord(i, dim);
            mol.coord(i, dim) += h;

            // Calculate the energy after adding h 
            double energy_plus_h = calculate_total_LJ(mol, epsilon_i, epsilon_j, sigma_i, sigma_j);

            // restore the original coordinate
            mol.coord(i, dim) = original_coord;

            double original_energy = calculate_total_LJ(mol, epsilon_i, epsilon_j, sigma_i, sigma_j);

            double force = - (energy_plus_h - original_energy) / h;

            // update the force matrix 
            forces(i, dim) = force;
        }
    }
    return forces;
}

// Approximate the analytical forces of a molecule using central difference approximation.
arma::mat LennardJones::calculate_LJforces_central(molecule mol, double h, double epsilon_i, double epsilon_j, double sigma_i, double sigma_j) {
    arma::mat forces = arma::zeros(mol.natoms, 3); // initialize the forces matrix
    arma::mat geom = mol.coord; 
    int natoms = mol.natoms;

    for (int i = 0; i < natoms; i++) {
        for (int dim = 0; dim < 3; dim++) {
            double original_coord = mol.coord(i, dim);
            mol.coord(i, dim) += h;

            // Calculate the energy after adding h 
            double energy_plus_h = calculate_total_LJ(mol, epsilon_i, epsilon_j, sigma_i, sigma_j);

            // restore the original coordinate
            mol.coord(i, dim) = original_coord;

            // Calculate the energy after subtracting h 
            mol.coord(i, dim) -= h;
            double energy_minus_h = calculate_total_LJ(mol, epsilon_i, epsilon_j, sigma_i, sigma_j);

            // calculate the force with respect to the position 
            double force = - (energy_plus_h - energy_minus_h) / (2 * h);

            // update the force matrix 
            forces(i, dim) = force;
        }
    }
    return forces;
}


