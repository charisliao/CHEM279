
#include <iostream>
#include <armadillo>
#include "AO.h"
#include "utils.h"
#include "CNDO.h"


using namespace std;


int main() {

    // read in file named "C2H2.txt" 
    // create AO object
    // print out the number of basis functions, number of electrons, and number of atoms
    // print out the basis set
    
    AO H2_ao("H2.txt");

    vector<BasisFunction> basis_set = H2_ao.basis_set;
    CNDO CNDO_H2;
    
    
// arma::mat createhamiltonianEnergy(vector<BasisFunction>& basis_set);
// arma::mat fancyH(arma::mat X, arma::mat H);
// arma::mat createX(arma::mat overlap_matrix);
// arma::mat MO_coefficients(arma::mat X, arma::mat Fancy_H);
// double calculateEnergy(arma::mat X, arma::mat Fancy_H, int num_electrons);
// double calculateHamiltonianMatrix(AO AO_object, arma::mat Overlap_matrix);
    // compute the gamma matrix

    std::map<int, BasisFunction> basis_map; 
    for (int i = 0; i < basis_set.size(); i++) {
        if (basis_set[i].AO_type.find("s") != std::string::npos) {
            basis_map[basis_set[i].atom_index] = basis_set[i];
        }
    }

    cout << "Gamma Matrix for H2: " << endl;
    int natoms = H2_ao.get_natoms();
    arma::mat gamma = CNDO_H2.computeGammaMatrix(natoms, basis_set);
    gamma.print();
    
    cout << "Overlap Matrix for H2: " << endl;
    arma::mat S = overlap_matrix(basis_set);
    S.print();

    cout << "Compute Core Hamiltonian Matrix for H: " << endl;
    vector<string> atom_types = H2_ao.get_atom_types();
    arma::mat Hcore = CNDO_H2.computeCoreHamiltonianMatrixOptimized(S, basis_map, atom_types, basis_set);
    Hcore.print();

    cout << "Compute Fock Matrix for H: " << endl;
    arma::mat density_matrix = arma::zeros(basis_set.size(), basis_set.size());
    arma::vec Ptotal = arma::zeros(natoms);
    arma::mat Fock = CNDO_H2.computeFockMatrixOptimized(basis_map, atom_types, basis_set, S, Hcore, density_matrix, Ptotal);
    Fock.print();

    CNDO_H2.updateDensityMatrixOptimized(basis_map, H2_ao, S, Hcore, gamma, atom_types, basis_set);

    return 0;
}