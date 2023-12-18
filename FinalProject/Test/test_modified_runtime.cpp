#include <iostream>
#include <fstream>
#include <vector>
#include <armadillo>
#include <functional>
#include <chrono>
#include <filesystem>  // Include filesystem header for directory traversal
#include <algorithm>   // Include algorithm header for sorting
#include <omp.h>
#include "AO.h"
#include "utils.h"
#include "CNDO.h"

using namespace std;
namespace fs = std::filesystem;

// Function to extract the numeric part of the filename
int extractNumber(const std::string& fileName) {
    size_t found = fileName.find_first_of("0123456789");
    return std::stoi(fileName.substr(found));
}

int main() {
    // Directory containing molecule files
    const std::string moleculesDir = "molecules";

    // Open the output file
    std::ofstream outputFile("optimized_output.txt");
    if (!outputFile.is_open()) {
        std::cerr << "Error: Unable to open output file." << std::endl;
        return 1;
    }

    // Vector to store sorted file entries
    std::vector<std::pair<std::string, int>> sortedEntries;

    for (const auto& entry : fs::directory_iterator(moleculesDir)) {
        const std::string fileName = entry.path().filename().string();

        // Skip non-text files
        if (fileName.find(".txt") == std::string::npos)
            continue;

        // Extract the numeric part of the filename
        int fileNumber = extractNumber(fileName);

        // Add to the vector for sorting
        sortedEntries.emplace_back(fileName, fileNumber);
    }

    // Sort the entries based on the numeric part
    std::sort(sortedEntries.begin(), sortedEntries.end(), [](const auto& a, const auto& b) {
        return a.second < b.second;
    });

    // Process files in sorted order
    for (const auto& entry : sortedEntries) {
        const std::string fileName = entry.first;
        const std::string filePath = moleculesDir + "/" + fileName;

        // Read in file
        AO ao(filePath);

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
            // Set the number of threads for the following parallel region
            #pragma omp parallel
            {
                int numThreads = omp_get_num_threads();
                // std::cout << "Number of threads for updateDensityMatrixOptimized: " << numThreads << std::endl;

                // Inside the parallel region
                #pragma omp single
                {
                    // Set the desired number of threads
                    omp_set_num_threads(numThreads);

                    cndo.updateDensityMatrixOptimized(basis_map, ao, S, Hcore, gamma, atom_types, basis_set);
                }
            }
        }, "Update Density Matrix (P)");

        std::cout << fileName << " processed." << std::endl;
        outputFile << ao.get_natoms() << " " << durationS << " " << durationH << " " << durationGamma << " " << durationF << " " << durationP << std::endl;
    }

    outputFile.close();

    return 0;
}
