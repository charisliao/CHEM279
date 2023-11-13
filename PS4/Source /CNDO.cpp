/**
 * @file CNDO.cpp
 * @author Charis Liao (charisliao@berkeley.edu)
 * @brief 
 * @version 0.1
 * @date 2023-11-05
 * 
 * @copyright Copyright (c) 2023
 * 
 */


#include <iostream>
#include <armadillo>
#include <cmath>
#include <string>
#include <vector>
#include <map>
#include "AO.h"
#include "util.h"
#include "CNDO.h"




// create CNDO constructor 
CNDO::CNDO() {
        // Define CNDO parameters for different elements
        CNDO_parameters H_para{{{"s", 7.176}}, 9.};
        CNDO_parameters C_para{{{"s", 14.051}, {"p", 5.572}}, 21.};
        CNDO_parameters N_para{{{"s", 19.316}, {"p", 7.275}}, 25.};
        CNDO_parameters O_para{{{"s", 25.390}, {"p", 9.111}}, 31.};
        CNDO_parameters F_para{{{"s", 32.272}, {"p", 11.080}}, 39.};

        // Populate CNDO_para_map with the parameters
        CNDO_para_map = {
            {"H", H_para},
            {"C", C_para},
            {"N", N_para},
            {"O", O_para},
            {"F", F_para}
        };
    }


double CNDO::getBeta(const string& element) {
    auto it = CNDO_para_map.find(element);
    if (it != CNDO_para_map.end()) {
        return it->second.beta;
    } else {
        cout << "This is not a valid element" << endl;
        return 0.0;
    }
}

double CNDO::getIA(const std::string& element, const std::string& orbitalType) const {
    auto it = CNDO_para_map.find(element);
    if (it != CNDO_para_map.end()) {
        // Treat "px", "py", and "pz" as "p"
        std::string orbitalKey = (orbitalType == "px" || orbitalType == "py" || orbitalType == "pz") ? "p" : orbitalType;

        auto orbitalIt = it->second.IA.find(orbitalKey);
        if (orbitalIt != it->second.IA.end()) {
            return orbitalIt->second;
        } else {
            std::cerr << "Invalid orbital type: " << orbitalType << std::endl;
            return 0.0;
        }
    } else {
        std::cerr << "Element not found: " << element << std::endl;
        return 0.0;
    }
}

double CNDO::gamma(BasisFunction& basisFunction1, BasisFunction& basisFunction2) {
    arma::vec la = basisFunction1.lmn, lb = basisFunction2.lmn;
    if(!(arma::accu(la) == 0 && arma::accu(lb) == 0)) {
        std::cout << "Only s orbitals allowed" << std::endl;
    }
    //information from basis set
    arma::vec da=basisFunction1.contracted_coefficients % basisFunction1.normalization_constants, db = basisFunction2.contracted_coefficients % basisFunction2.normalization_constants;
    // cout << "Done initializing da and db" <<endl;
    arma::vec alphaa = basisFunction1.exponents, alphab = basisFunction2.exponents;
    // cout << "Done initializing alphaa and alphab"<< endl;
    arma::rowvec Ra = basisFunction1.center, Rb = basisFunction2.center;
    // cout << "Done initializing Ra and Rb" << endl;
    int len = basisFunction1.exponents.size();
    // cout << "Initialization complete" << endl;

    double sum = 0;
    for (int k1=0; k1 < len; k1++) {
        for (int k2=0; k2 < len; k2++) {
            double sigmaA = 1.0/(alphaa(k1) + alphaa(k2)); //eq 3.10

            for (int j1 =0; j1 < len; j1++) {
                for (int j2 = 0; j2 < len; j2++) {
                    double sigmaB = 1.0/(alphab(j1) + alphab(j2)); //eq 3.10
                    double I2e = I2e_pG(Ra, Rb, sigmaA, sigmaB); //eq 3.14

                    sum += da(k1) * da(k2) * db(j1) * db(j2) * I2e; //eq 3.13
                    // cout << "sum: " << sum << endl;
                }
            }
        }

    }
    // cout << "sum: " << sum << endl;
    // cout << "hartree_to_ev: " << hartree_to_ev << endl;
    // cout << "converted sum: " << sum * hartree_to_ev << endl;
    return sum * hartree_to_ev; 
}

arma::mat CNDO::computeGammaMatrix(int natoms, vector<BasisFunction>& basis_set) {
    arma::mat gamma_mat = arma::zeros(natoms, natoms);
    cout << "gamma matrix initialized" << endl;
    int valid_i = 0; // Counter for valid 'i' indices
    int valid_j = 0; // Counter for valid 'j' indices

    // Loop through the basis set and only consider S orbitals 
    for (int i = 0; i < basis_set.size(); i++) {
        if (basis_set[i].AO_type.find("s") != std::string::npos) { // Ensure it's an s orbital
            valid_j = 0; // Reset valid_j for each valid i
            for (int j = 0; j < basis_set.size(); j++) {
                if (basis_set[j].AO_type.find("s") != std::string::npos) { // Ensure it's an s orbital
                    // cout << "inside gamma matrix call gamma" << endl;
                    // cout << "i" << i << ","  << " j: " << j << " ," << "gamma: " << gamma(basis_set[i], basis_set[j]) << endl;
                    // Insert gamma into the matrix at the valid indices
                    gamma_mat(valid_i, valid_j) = gamma(basis_set[i], basis_set[j]);
                    valid_j++; // Increment valid_j for each valid j
                }
            }
            valid_i++; // Increment valid_i for each valid i
        }
    }
    return gamma_mat;
}


double CNDO::diagCoreH(double gammaAA, BasisFunction& basisFunctionA, double ZB_gamma_AB) {
    double diag_element = 0.0; 
    // cout << "element" << extractElement(basisFunctionA.AO_type) << endl;
    // cout << "orbital" << extractOrbital(basisFunctionA.AO_type) << endl;
    double IA = getIA(extractElement(basisFunctionA.AO_type), extractOrbital(basisFunctionA.AO_type));
    // cout << "Calculate diag element" << endl;
    // cout << "IA: " << IA << endl; 
    // cout << "i" << i << ","  << " j: " << j << " ," << "gammaAA: " << gammaAA << endl;
    // cout << "ZB_gamma_AB: " << ZB_gamma_AB << endl;
    diag_element = -IA - (basisFunctionA.valence_e - 0.5) * gammaAA - ZB_gamma_AB;
    
    return diag_element;

}

std::string CNDO::extractElement(const std::string& AO_type) {
    return AO_type.substr(0, AO_type.find("_"));
}

std::string CNDO::extractOrbital(const std::string& AO_type) {
    return AO_type.substr(AO_type.find("_") + 1);
}

double CNDO::computeOffDiagMat(int i, int j, double Beta_A, double Beta_B, arma::mat overlapMatrix) {
    double off_diag_element = 0.0;
    // cout << " Inside Compute Off Diag Mat" << endl;
    // cout << "i: " << i << " j: " << j << endl;
    // cout << "Beta_A: " << Beta_A << " Beta_B: " << Beta_B << endl;
    // cout << "overlapMatrix(i, j): " << overlapMatrix(i, j) << endl;
    off_diag_element = -1 * 0.5 * (Beta_A + Beta_B) * overlapMatrix(i, j) ;
    return off_diag_element;
}

arma::mat CNDO::computeCoreHamiltonianMatrix(vector<string> atom_types, vector<BasisFunction>& basis_set) {
    arma::mat H = arma::zeros(basis_set.size(), basis_set.size());
    arma::mat gamma_mat = computeGammaMatrix(atom_types.size(), basis_set);
    arma::mat overlap_mat = overlap_matrix(basis_set);
    int atom_index = 0;

     // Print the vector elements
    // std::cout << "Vector contents: ";
    // for (const auto& element : atom_types) {
    //     std::cout << element << " ";
    // }
    // std::cout << std::endl;
    
    //create a map to store atom index and the s orbitle basis functions 
    std::map<int, BasisFunction> basis_map; 
    for (int i = 0; i < basis_set.size(); i++) {
        if (basis_set[i].AO_type.find("s") != std::string::npos) {
            atom_index = basis_set[i].atom_index;
            basis_map[atom_index] = basis_set[i];
        }
    }
    
    for (int i = 0; i < basis_set.size(); i++) {
        for (int j = 0; j < basis_set.size(); j++) {
            
            if (i == j) {
                int atom_indexA = basis_set[i].atom_index;
                double gammaAA = 0.0;
                BasisFunction basisFunctionA = basis_map[atom_indexA];
                gammaAA = gamma(basisFunctionA, basisFunctionA);

                cout << "gammaAA: " << gammaAA << endl;
                double ZB_gamma_AB = 0.0;

                for (int index = 0; index < atom_types.size(); index++) {
                    if (atom_indexA == index) {
                        continue;
                    } else {
                        BasisFunction basisFunctionB = basis_map[index];
                        double ZB = basisFunctionB.valence_e;
                        double gammaAB = gamma(basisFunctionA, basisFunctionB);
                        ZB_gamma_AB += ZB * gammaAB;
                    }
                }
                cout << "ZB_gamma_AB: " << ZB_gamma_AB << endl;
                H(i, j) = diagCoreH(gammaAA, basis_set[i], ZB_gamma_AB);
                cout << "Current atom index: " << atom_index <<  " Current atom type: " << atom_types[atom_index] << endl;
                
            } else {
                
                double Beta_A = getBeta(extractElement(basis_set[i].AO_type));
                double Beta_B = getBeta(extractElement(basis_set[j].AO_type));
                H(i, j) = computeOffDiagMat(i, j, Beta_A, Beta_B, overlap_mat);
            }
        }
    }
    return H;
}



