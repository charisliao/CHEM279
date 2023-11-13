
#include <iostream>
#include <armadillo>
#include "AO.h"
#include "utils.h"
#include "CNDO.h"



using namespace std;


int main() {

    // read in file named "C2HF+.txt" 
    // create AO object
    // print out the number of basis functions, number of electrons, and number of atoms
    // print out the basis set
    
    AO HF_ao("HF+.txt");

    cout << "Overlap Matrix for HF+: " << endl;
    vector<BasisFunction> basis_set = HF_ao.basis_set;
    
// arma::mat createhamiltonianEnergy(vector<BasisFunction>& basis_set);
// arma::mat fancyH(arma::mat X, arma::mat H);
// arma::mat createX(arma::mat overlap_matrix);
// arma::mat MO_coefficients(arma::mat X, arma::mat Fancy_H);
// double calculateEnergy(arma::mat X, arma::mat Fancy_H, int num_electrons);
// double calculateHamiltonianMatrix(AO AO_object, arma::mat Overlap_matrix);
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
    arma::mat Hcore = CNDO_HF.computeCoreHamiltonianMatrix(atom_types, basis_set);
    Hcore.print();

    



    return 0;
}