/**
 * @file lennard-jones.h
 * @author Charis Liao
 * @brief The header file for the Lennard Jones class 
 * @version 0.1
 * @date 2023-09-10
 * 
 * @copyright Copyright (c) 2023
 * 
 */

#pragma once // Ensures that this file is only included once during compilation
#include "molecule.h" // molecule class


class LennardJones {
    private: 
        // constants 
        double epsilon_i, epsilon_j, sigma_i, sigma_j;

        // functions
        double calculate_epsilon_ij(double epsilon_i, double epsilon_j);
        double calculate_sigma_ij(double sigma_i, double sigma_j);
        double calculate_LJ(double epsilon_ij, double sigma_ij, double r_ij);
        double calculate_distance(double x_i, double y_i, double z_i, double x_j, double y_j, double z_j);
        double calculate_LJ_force(double epsilon_ij, double sigma_ij, double r_ij, double x_i, double x_j);

        // molecule object
        molecule mol_; 
    public: 
        LennardJones(const molecule& mol);
        double calculate_total_LJ(molecule mol, double epsilon_i, double epsilon_j, double sigma_i, double sigma_j);
        arma::mat calculate_LJ_forces(molecule mol, double epsilon_i, double epsilon_j, double sigma_i, double sigma_j);
        arma::mat calculate_LJforces_forward(molecule mol, double h, double epsilon_i, double epsilon_j, double sigma_i, double sigma_j);
        arma::mat calculate_LJforces_central(molecule mol, double h, double epsilon_i, double epsilon_j, double sigma_i, double sigma_j);
        double calculate_force_norm(arma::mat forces);
};