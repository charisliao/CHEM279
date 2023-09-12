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
#include "lennard_jones.h"

class SteepestDescent {
    private: 
        // molecule object
        molecule mol_; 
    public: 
        SteepestDescent(const molecule& mol);
        void steepest_descent_optimizer(molecule mol, double stepsize, double tol, double epsilon_i, double epsilon_j, double sigma_i, double sigma_j);
};
