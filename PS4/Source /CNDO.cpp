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
#include <fstream>
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
    // if(!(arma::accu(la) == 0 && arma::accu(lb) == 0)) {
    //     std::cout << "Only s orbitals allowed" << std::endl;
    // }
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
    // cout << "gamma matrix initialized" << endl; 
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
                

                // cout << "gammaAA: " << gammaAA << endl;
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
                // cout << "ZB_gamma_AB: " << ZB_gamma_AB << endl;
                H(i, j) = diagCoreH(gammaAA, basis_set[i], ZB_gamma_AB);
                // cout << "Current atom index: " << atom_index <<  " Current atom type: " << atom_types[atom_index] << endl;
                
            } else {
                
                double Beta_A = getBeta(extractElement(basis_set[i].AO_type));
                double Beta_B = getBeta(extractElement(basis_set[j].AO_type));
                H(i, j) = computeOffDiagMat(i, j, Beta_A, Beta_B, overlap_mat);
            }
        }
    }
    return H;
}
 
double CNDO::offDiagonalFockElement(double BetaA, double BetaB, double overlap, double p, double gammaAB) {
    return 0.5 * (-BetaA - BetaB) * overlap - p * gammaAB;
}

double CNDO::diagonalFockElement(double IA, double ZA, double gammaAA, double pTotalA, double pMuMu, double pZBgammaAB) {
    return -IA + ((pTotalA - ZA) - (pMuMu - 0.5)) * gammaAA + pZBgammaAB;
}

arma::mat CNDO::computeFockMatrix(vector<string> atom_types, vector<BasisFunction>& basis_set, arma::mat& overlap_mat, arma::mat& Hcore_mat, arma::mat& p, arma::vec& Ptotal) {
    std::map<int, BasisFunction> basis_map; 
    int atom_index;
    for (int i = 0; i < basis_set.size(); i++) {
        if (basis_set[i].AO_type.find("s") != std::string::npos) {
            atom_index = basis_set[i].atom_index;
            basis_map[atom_index] = basis_set[i];
        }
    }
    arma::mat F = arma::zeros(basis_set.size(), basis_set.size());
    
    for (int i = 0; i < basis_set.size(); i++) {
        for (int j = 0; j < basis_set.size(); j++) {
            if (i == j) {
                double IA = getIA(extractElement(basis_set[i].AO_type), extractOrbital(basis_set[i].AO_type));
                int atom_indexA = basis_set[i].atom_index;
                double gammaAA = 0.0;
                BasisFunction basisFunctionA = basis_map[atom_indexA];
                gammaAA = gamma(basisFunctionA, basisFunctionA);
                double ZA = basisFunctionA.valence_e;
                double pTotalA = Ptotal(atom_indexA);

                // cout << "gammaAA: " << gammaAA << endl;
                double pZB_gamma_AB = 0.0;
                double pMuMu = p(i, j);

                for (int index = 0; index < atom_types.size(); index++) {
                    if (atom_indexA == index) {
                        continue;
                    } else {
                        BasisFunction basisFunctionB = basis_map[index];
                        double ZB = basisFunctionB.valence_e;
                        double gammaAB = gamma(basisFunctionA, basisFunctionB);
                        double pTotalB = Ptotal(index);
                        pZB_gamma_AB += (pTotalB - ZB) * gammaAB;
                    }
                }
                // cout << "pZB_gamma_AB: " << pZB_gamma_AB << endl;
                F(i, j) = diagonalFockElement(IA, ZA, gammaAA, pTotalA, pMuMu, pZB_gamma_AB);
                
                
            } else {
                
                double Beta_A = getBeta(extractElement(basis_set[i].AO_type));
                double Beta_B = getBeta(extractElement(basis_set[j].AO_type));
                double overlap = overlap_mat(i, j);
                double density = p(i, j);
                double gammaAB = gamma(basis_set[i], basis_set[j]);
                F(i, j) = offDiagonalFockElement(Beta_A, Beta_B, overlap, density, gammaAB);
            }
        }
    }
    return F;
}

/**
 * @brief 
 * 
 * @param AO_object 
 * Step1: Guess p_alpha = p_beta = 0
 * Step2: Build f_alpha and f_beta
 * Step3: Solve the eigenvalue problems to obtain new MO coefficients, C_alpha and C_beta and 
 * corresponding eigenvalues epsilon_alpha and epsilon_beta.
 * Step4: Copy the old density matrices to a new variable, P_alpha_old and P_beta_old
 * Step5: Build the new density matrices, P_alpha and P_beta by occupying the p lowest energy alpha
 * and q lowest energy beta MOs, respectively.
 * Step6: If the maximum magnitude of the change in the alpha na dbeta density matrices is less than 
 * a tolerance (10e-6), the we have converged. Otherwise, go back to step 2.
 */
arma::vec CNDO::updatePtotal(vector<string> atom_types, vector<BasisFunction> basis_set, arma::mat& densityMat) {
    int num_atoms = atom_types.size();
    arma::vec newPtotal = arma::zeros(num_atoms);
    for (int i = 0; i < basis_set.size(); i++) {
        for (int j = 0; j < basis_set.size(); j++) {
            if (i == j && basis_set[i].atom_index == basis_set[j].atom_index) {
                newPtotal(basis_set[i].atom_index) += densityMat(i, j);
            }
        }
    }
    return newPtotal;

}

double CNDO::calculateNuclearRepulsionEnergy(vector<string> atom_types, AO AO_object) {
    double nuclear_repulsion_energy = 0.0;
    arma::mat coord = AO_object.get_coord();
    // create valence electron map 
    map<string, int> valence_e_map;
    valence_e_map["H"] = 1.0;
    valence_e_map["C"] = 4.0;
    valence_e_map["N"] = 5.0;
    valence_e_map["O"] = 6.0;
    valence_e_map["F"] = 7.0;

    for (int i = 0; i < atom_types.size(); i++) {
        for (int j = 0; j < i; j++) {
            double ZA = valence_e_map[atom_types[i]];
            double ZB = valence_e_map[atom_types[j]];
            double RAB = arma::norm(coord.row(i) - coord.row(j));
            nuclear_repulsion_energy += ZA * ZB / RAB;
        }
    }
    return nuclear_repulsion_energy * hartree_to_ev;
}

double CNDO::calculateElectronEnergy(arma::mat& densityMat, arma::mat& Hcore_mat, arma::mat& Fock_mat) {
    double electron_energy = 0.0;
    for (int i = 0; i < densityMat.n_rows; i++) {
        for (int j = 0; j < densityMat.n_cols; j++) {
            electron_energy += densityMat(i, j) * (Hcore_mat(i, j) + Fock_mat(i, j));
        }
    }
    return 0.5 * electron_energy;
}

void CNDO::updateDensityMatrix(AO AO_object, std::string output_file_name) {
    std::ofstream outputFile(output_file_name);
    if (!outputFile.is_open()) {
        std::cerr << "Error: Unable to open output file " << output_file_name << std::endl;
        return;
    }
    // Redirect cout to the output file
    std::streambuf* coutBuffer = std::cout.rdbuf();
    std::cout.rdbuf(outputFile.rdbuf());

    // std::cout << "State Updating Density Matrix" << std::endl;
    arma::mat densityAlpha = arma::zeros(AO_object.basis_set.size(), AO_object.basis_set.size());
    arma::mat densityBeta = arma::zeros(AO_object.basis_set.size(), AO_object.basis_set.size());
    arma::mat coefficientAlpha;
    arma::mat coefficientBeta;
    arma::vec pTotal = arma::zeros(AO_object.get_natoms());


    arma::mat totalDensityMat; 
    //initialize eigan values for alpha and beta
    arma::vec epsilonAlpha;
    arma::vec epsilonBeta;
    arma::mat FockAlpha;
    arma::mat FockBeta;


    int p = AO_object.get_p();
    int q = AO_object.get_q();

    vector<string> atom_types = AO_object.get_atom_types();
    vector<BasisFunction> basis_set = AO_object.basis_set;
    arma::mat overlap_mat = overlap_matrix(basis_set);
    arma::mat Hcore_mat = computeCoreHamiltonianMatrix(atom_types, basis_set);

    bool converged = false;

    int iteration = 0;
    while (converged != true) {
        std::cout << "Iteration: " << iteration << std::endl;

        std::cout << "Fock Matrix for Alpha: " << std::endl;
        FockAlpha = computeFockMatrix(atom_types, basis_set, overlap_mat, Hcore_mat, densityAlpha, pTotal);
        FockAlpha.print();

        std::cout << "Fock Matrix for Beta: " << std::endl;
        FockBeta = computeFockMatrix(atom_types, basis_set, overlap_mat, Hcore_mat, densityBeta, pTotal);
        FockBeta.print();

        // solve eigan value problem for alpha and beta
        arma::eig_sym(epsilonAlpha, coefficientAlpha, FockAlpha);
        arma::eig_sym(epsilonBeta, coefficientBeta, FockBeta);
        std::cout << "after solving eigen quation: " << iteration << std::endl;

        std::cout << "MO Coefficient Matrix for Alpha: " << std::endl;
        coefficientAlpha.print();

        std::cout << "MO Coefficient Matrix for Beta: " << std::endl;
        coefficientBeta.print();

        std::cout << "p: " << p << " q: " << q << std::endl;

        // copy old density matrix to new density matrix
        arma::mat densityAlpha_old = densityAlpha;
        arma::mat densityBeta_old = densityBeta;

        // build new density matrix
        densityAlpha = arma::zeros(densityAlpha.n_rows, densityAlpha.n_cols);
        densityBeta = arma::zeros(densityBeta.n_rows, densityBeta.n_cols);

        for (int i = 0; i < densityAlpha.n_rows; i++) {
            for (int j = 0; j < densityAlpha.n_cols; j++) {
                // densityAlpha(i, j) = 0.0
                // densityBeta(i, j) = 0.0;
                for (int k = 0; k < p; k++) {
                    densityAlpha(i, j) += coefficientAlpha(i, k) * coefficientAlpha(j, k);
                }
                for (int k = 0; k < q; k++) {
                    densityBeta(i, j) += coefficientBeta(i, k) * coefficientBeta(j, k);
                }
                
            }
        }

        std::cout << "Density Matrix for Alpha: " << std::endl;
        densityAlpha.print();

        std::cout << "Density Matrix for Beta: " << std::endl;
        densityBeta.print();

        totalDensityMat = densityAlpha + densityBeta;
        std::cout << "Total Density Matrix: " << std::endl;
        totalDensityMat.print();

        pTotal = updatePtotal(atom_types, basis_set, totalDensityMat);
        std::cout << "Total Density Vector" << std::endl;
        pTotal.print();

        // / 5) Check for convergence (tolerance = 1e-6)
        if (abs(densityAlpha - densityAlpha_old).max() < 1e-6 && \
            abs(densityBeta - densityBeta_old).max() < 1e-6) {
                std::cout << "Converged!" << std::endl;
                converged = true;
                std::cout << "Ea" << std::endl;
                epsilonAlpha.print();
                std::cout << "Eb" << std::endl;
                epsilonBeta.print();
                std::cout << "Ca" << std::endl;
                coefficientAlpha.print();
                std::cout << "Cb" << std::endl;
                coefficientBeta.print();

                double nuclear_repulsion_energy = calculateNuclearRepulsionEnergy(atom_types, AO_object);
                std::cout << "Nuclear Repulsion Energy: " << nuclear_repulsion_energy << "eV" << std::endl;
                double electron_energy = calculateElectronEnergy(densityAlpha, Hcore_mat, FockAlpha);
                electron_energy += calculateElectronEnergy(densityBeta, Hcore_mat, FockBeta);
                std::cout << "Electron Energy: " << electron_energy << "eV" << std::endl;
                std::cout << "Total Energy: " << electron_energy + nuclear_repulsion_energy << "eV" << std::endl;

        }

        

        iteration++;
    }
    // Redirect cout back to the console
    std::cout.rdbuf(coutBuffer);

    // Close the output file
    outputFile.close();
}






