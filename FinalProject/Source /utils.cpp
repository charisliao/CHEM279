
#include "utils.h"
#include <armadillo> 
#include <cmath>
#include <chrono>
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>



double factorial(double n) {
    if (n == 0 || n == 1) {
        return 1;
    }
    else {
        return n * factorial(n - 1);
    }
}

// double factorial function
double double_factorial(double n) {
    if (n <= 1) {
        return 1;
    }
    else {
        return n * double_factorial(n - 2);
    }
}

// binomial coefficient function
double binomial(int m, int n) {
    if (m < n) {
        return 0.0;
    }
    return factorial(m) / (factorial(n) * factorial(m - n));
}

// Center of the Product 
double center_of_product(double center_vec1, double alpha1, double center_vec2, double alpha2) {
    return (center_vec1 * alpha1 + center_vec2 * alpha2) / (alpha1 + alpha2);
}

// exponential prefector 
double prefactor(double center_vec1, double alpha1, double center_vec2, double alpha2) {
    double r = center_vec1 - center_vec2;
    return exp(-(alpha1 * alpha2 * r * r) / (alpha1 + alpha2)) * pow(M_PI / (alpha1 + alpha2), 0.5);
}

double I2e_pG(arma::rowvec &Ra, arma::rowvec &Rb, double sigmaA, double sigmaB) {
    double U = pow(M_PI * sigmaA, 1.5) * pow(M_PI * sigmaB, 1.5); //eq 3.8 and 3.11
    double V2 = 1.0 / (sigmaA+sigmaB); //eq. 3.9 

    double Rd = arma::norm(Ra - Rb, 2);
    
    if (Rd == 0.0) {
        //use eq 3.15
        return U * sqrt(2*V2) * sqrt(2/M_PI);
    }
    //eq. 3.14 need sqrt T
    double srT = sqrt(V2) * Rd; 

    // eq 3.14
    double result = U / Rd * erf(srT); 
    return result; 

}


// Function to measure the execution time of a given function
double measureExecutionTime(const std::function<void()>& func, const std::string& functionName) {
    

    auto start = std::chrono::high_resolution_clock::now();
    func();
    auto end = std::chrono::high_resolution_clock::now();

    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();
    // Print details to the terminal
    if (!functionName.empty()) {
        std::cout << "   " << functionName << ": " << duration << " milliseconds." << std::endl;
    }

    return static_cast<double>(duration);
}

void generateMoleculeFile(int num_C, int num_H, const std::string& filename) {
    std::ofstream outputFile(filename);
    if (!outputFile.is_open()) {
        std::cerr << "Error: Unable to open file " << filename << " for writing." << std::endl;
        return;
    }

    std::vector<std::string> molecule;

    // Add C atoms
    for (int i = 0; i < num_C; ++i) {
        molecule.push_back("C    " + std::to_string(1.5 * i) + "    0.000    0.0");
    }

    // Add H atoms
    for (int i = 0; i < num_H; ++i) {
        molecule.push_back("H    " + std::to_string(1.5 * num_C + 1.0 * i) + "    0.000    0.0");
    }

    // Combine all lines into a single string and write to file
    outputFile << num_C + num_H << " 0\n";
    for (const auto& line : molecule) {
        outputFile << line << "\n";
    }

    std::cout << "Molecule generated and saved to " << filename << std::endl;
}