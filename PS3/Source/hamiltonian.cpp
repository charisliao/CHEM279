/**
 * @file hamiltonianEnergy.cpp
 * @author Charis Liao (charisliao@berkeley.edu)
 * @brief 
 * @version 0.1
 * @date 2023-10-16
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
#include "hamiltonian.h"



arma::mat createhamiltonianEnergy(vector<BasisFunction>& basis_set) {
    
    
    arma::mat S = overlap_matrix(basis_set);
    arma::mat H = arma::zeros(basis_set.size(), basis_set.size());

    // Create a that stores the orbital and it's corresponding energy 
    map<string, double> hamiltonianEnergy;
    hamiltonianEnergy.insert({"H_s", -13.6});
    hamiltonianEnergy.insert({"C_s", -21.4});
    hamiltonianEnergy.insert({"C_px", -11.4});
    hamiltonianEnergy.insert({"C_py", -11.4});
    hamiltonianEnergy.insert({"C_pz", -11.4});
    
    // symmetric orthogoalization of S 
    for (int i = 0; i < basis_set.size(); i++) {
        for (int j = 0; j < basis_set.size(); j++) {
            
            string AO_type1 = basis_set[i].AO_type;
            string AO_type2 = basis_set[j].AO_type;
            double H_i = hamiltonianEnergy[AO_type1];
            double H_j = hamiltonianEnergy[AO_type2];

            if (i == j) {
                H(i, j) = H_i;
            } else {
                H(i, j) = K / 2 * (H_i + H_j) * S(i, j);
            }
        }
    }
    return H;    
}


arma::mat createX(arma::mat overlap_matrix) {
    
   
    // get the eigenvector and matrix 
    arma::vec s;
    arma::mat U;
    arma::eig_sym(s, U, overlap_matrix);

    // s.print("S Eigenvalues: ");
    // U.print("S Eigenvectors: ");
// get inverse sqrt of eigenvalues and put them into a diagonal matrix 
    for (int i = 0; i <s.size(); i++) {
        s(i) = 1.0 / std::sqrt(s(i));
    }

    arma::mat s_diagonal_matrix = arma::diagmat(s);
    arma::mat X = U * s_diagonal_matrix * U.t();

    
    return X;
}

arma::mat fancyH(arma::mat X, arma::mat H) {
    arma::mat fancyH = arma::zeros(H.n_rows, H.n_cols);

    fancyH = X.t() * H * X;
    return fancyH;
}

arma::mat MO_coefficients(arma::mat X, arma::mat Fancy_H) {
    arma::vec eiganval_energy; 
    arma::mat V; 
    arma::eig_sym(eiganval_energy, V, Fancy_H);

    arma::mat MO_Coefficients = X * V;
    return MO_Coefficients;
}

double calculateEnergy(arma::mat X, arma::mat Fancy_H, int num_electrons) {
    arma::vec eiganval_energy; 
    arma::mat V; 
    arma::eig_sym(eiganval_energy, V, Fancy_H);

    arma::mat MO_Coefficients = X * V;

    double energy = 0.0; 
    for (int i = 0; i < num_electrons; i++) {
        energy += eiganval_energy(i) * 2;
    }

    return energy;
}

double calculateHamiltonianEnergy(AO AO_object, arma::mat Overlap_matrix) {
    vector<BasisFunction> basis_set = AO_object.basis_set;
    
    arma::mat H = createhamiltonianEnergy(basis_set);
    arma::mat X = createX(Overlap_matrix);
    arma::mat Fancy_H = fancyH(X, H);
    // cout << "Number of electrons: " << AO_object.num_electrons << endl;
    double energy = calculateEnergy(X, Fancy_H, AO_object.num_electrons);
    return energy;
}




