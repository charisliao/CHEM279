
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
    cout << "H: " << endl;
    AO H_ao("H.txt");
    cout << "Done initializing H" << endl;
    vector<BasisFunction> basis_set = H_ao.basis_set;

    std::map<int, BasisFunction> basis_map; 
    for (int i = 0; i < basis_set.size(); i++) {
        if (basis_set[i].AO_type.find("s") != std::string::npos) {
            basis_map[basis_set[i].atom_index] = basis_set[i];
        }
    }

    cout << "Overlap Matrix for H: " << endl;
    
    

    arma::mat S = overlap_matrix(basis_set);
    S.print();

    // create a CNDO instance 

    CNDO CNDO_H;
    // compute the gamma matrix 
    cout << "Compute Gamma Matrix for H: " << endl;
    int natoms = H_ao.get_natoms();
    arma::mat gamma = CNDO_H.computeGammaMatrix(natoms, basis_set);
    gamma.print();

    cout << "Compute Core Hamiltonian Matrix for H: " << endl;
    vector<string> atom_types = H_ao.get_atom_types();
    arma::mat Hcore = CNDO_H.computeCoreHamiltonianMatrixOptimized(S, basis_map, atom_types, basis_set);
    Hcore.print();

    cout << "Compute Fock Matrix for H: " << endl;
    arma::mat density_matrix = arma::zeros(basis_set.size(), basis_set.size());
    arma::vec Ptotal = arma::zeros(natoms);
    arma::mat Fock = CNDO_H.computeFockMatrixOptimized(basis_map, atom_types, basis_set, S, Hcore, density_matrix, Ptotal);
    Fock.print();
    
    CNDO_H.updateDensityMatrixOptimized(basis_map, H_ao, S, Hcore, gamma, atom_types, basis_set);



    return 0;
}