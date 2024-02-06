#include <iostream>
#include <armadillo>
#include <vector>
#include <string>

#pragma once // Ensures that this file is only included once during compilation

// Create a class called molecule 
class molecule{
    public: 

    int natoms; 
    arma::mat coord;
    std::vector<int> atomic_numbers;

    // Constructor
    molecule(const char *filename);
};