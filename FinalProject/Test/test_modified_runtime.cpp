#include <iostream>
#include <fstream>
#include <vector>
#include <armadillo>
#include <functional>
#include <chrono>
#include "AO.h"
#include "utils.h"
#include "CNDO.h"

using namespace std;

int main() {
    // File names
    vector<string> fileNames = {"H.txt", "C2H4.txt", "C50H50.txt", "C250H250.txt", "C500H500.txt"};

    // Open the output file
    std::ofstream outputFile("output_optimized.txt");
    if (!outputFile.is_open()) {
        std::cerr << "Error: Unable to open output file." << std::endl;
        return 1;
    }

    for (const auto& fileName : fileNames) {
        // Read in file
        AO ao(fileName);

        // Compute gamma and hcore matrices
        CNDO cndo;
        int natoms = ao.get_natoms();
        vector<string> atom_types = ao.get_atom_types();
        vector<BasisFunction> basis_set = ao.basis_set;
        std::map<int, BasisFunction> basis_map; 
        for (int i = 0; i < basis_set.size(); i++) {
            if (basis_set[i].AO_type.find("s") != std::string::npos) {
                basis_map[basis_set[i].atom_index] = basis_set[i];
            }
        }
        arma::mat gamma = cndo.computeGammaMatrix(natoms, ao.basis_set);
        arma::mat S = overlap_matrix(basis_set);
        arma::mat Hcore = cndo.computeCoreHamiltonianMatrixOptimized(S, basis_map, atom_types, basis_set);
        arma::mat p = arma::zeros(basis_set.size(), basis_set.size());
        arma::vec Ptotal = arma::zeros(ao.get_natoms());

        // Measure and store execution time for S
        double durationS = measureExecutionTime([&]() {
            vector<BasisFunction> basis_set = ao.basis_set;
            arma::mat S = overlap_matrix(basis_set);
        }, "Compute Overlap Matrix (S)");

        // Measure and store execution time for H
        double durationH = measureExecutionTime([&]() {
            CNDO cndo;
            arma::mat Hcore = cndo.computeCoreHamiltonianMatrixOptimized(S, basis_map, atom_types, basis_set);
        }, "Compute Core Hamiltonian Matrix (H)");

        double durationGamma = measureExecutionTime([&]() {
            CNDO cndo;
            int natoms = ao.get_natoms();
            arma::mat gamma = cndo.computeGammaMatrix(natoms, ao.basis_set);
        }, "Compute Gamma Matrix");

        // Measure and store execution time for F
        double durationF = measureExecutionTime([&]() {
            CNDO cndo;
            arma::mat Fock = cndo.computeFockMatrixOptimized(basis_map, atom_types, basis_set, S, Hcore, p, Ptotal);
        }, "Compute Fock Matrix (F)");

        // Measure and store execution time for P
        double durationP = measureExecutionTime([&]() {
            CNDO cndo;
            cndo.updateDensityMatrixOptimized(basis_map, ao, S, Hcore, gamma, atom_types, basis_set);
        }, "Update Density Matrix (P)");

        outputFile << fileName << ":" << std::endl;

        outputFile << durationS << " " << durationH << " " << durationGamma << " " << durationF << " " << durationP << std::endl;
    }

    outputFile.close();

    return 0;
}
