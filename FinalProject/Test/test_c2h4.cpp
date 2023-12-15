
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
    
    AO C2H4_ao("C2H4.txt");

    cout << "Overlap Matrix for H2: " << endl;
    vector<BasisFunction> basis_set = C2H4_ao.basis_set;
    
    std::map<int, BasisFunction> basis_map; 
    for (int i = 0; i < basis_set.size(); i++) {
        if (basis_set[i].AO_type.find("s") != std::string::npos) {
            basis_map[basis_set[i].atom_index] = basis_set[i];
        }
    }
 
    arma::mat S = overlap_matrix(basis_set);
    S.print();

    // create a CNDO instance
    CNDO CNDO_C2H4;
    // compute the gamma matrix
    cout << "Compute Gamma Matrix for C2H4: " << endl;
    int natoms = C2H4_ao.get_natoms();
    cout << "natoms: " << natoms << endl;
    arma::mat gamma = CNDO_C2H4.computeGammaMatrix(natoms, basis_set);
    gamma.print();

    cout << "Compute Core Hamiltonian Matrix for H: " << endl;
    vector<string> atom_types = C2H4_ao.get_atom_types();
    arma::mat Hcore = CNDO_C2H4.computeCoreHamiltonianMatrixOptimized(S, basis_map, atom_types, basis_set);
    Hcore.print();

    cout << "Compute Fock Matrix for H: " << endl;
    arma::mat density_matrix = arma::zeros(basis_set.size(), basis_set.size());
    arma::vec Ptotal = arma::zeros(natoms);
    arma::mat Fock = CNDO_C2H4.computeFockMatrixOptimized(basis_map, atom_types, basis_set, S, Hcore, density_matrix, Ptotal);
    Fock.print();

    CNDO_C2H4.updateDensityMatrixOptimized(basis_map, C2H4_ao, S, Hcore, gamma, atom_types, basis_set);


    return 0;
}