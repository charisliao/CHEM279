
#include <iostream>
#include <armadillo>
#include "AO.h"
#include "utils.h"
#include "CNDO.h"



using namespace std;


int main() {
    
    AO HF_ao("HF+.txt");

    cout << "Overlap Matrix for HF+: " << endl;
    vector<BasisFunction> basis_set = HF_ao.basis_set;

    std::map<int, BasisFunction> basis_map; 
    for (int i = 0; i < basis_set.size(); i++) {
        if (basis_set[i].AO_type.find("s") != std::string::npos) {
            basis_map[basis_set[i].atom_index] = basis_set[i];
        }
    }

    arma::mat S = overlap_matrix(basis_set);
    S.print();
    CNDO CNDO_HF;
    // compute the gamma matrix
    cout << "Gamma Matrix for HF+: " << endl;
    int natoms = HF_ao.get_natoms();
    arma::mat gamma = CNDO_HF.computeGammaMatrix(natoms, basis_set);
    gamma.print();

    cout << "Compute Core Hamiltonian Matrix for H: " << endl;
    vector<string> atom_types = HF_ao.get_atom_types();
    arma::mat Hcore = CNDO_HF.computeCoreHamiltonianMatrixOptimized(S, basis_map, atom_types, basis_set);
    Hcore.print();

    cout << "Compute Fock Matrix for H: " << endl;
    arma::mat density_matrix = arma::zeros(basis_set.size(), basis_set.size());
    arma::vec Ptotal = arma::zeros(natoms);
    arma::mat Fock = CNDO_HF.computeFockMatrixOptimized(basis_map, atom_types, basis_set, S, Hcore, density_matrix, Ptotal);
    Fock.print();

    CNDO_HF.updateDensityMatrixOptimized(basis_map, HF_ao, S, Hcore, gamma, atom_types, basis_set);



    return 0;
}