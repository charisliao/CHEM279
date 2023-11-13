
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

    cout << "Overlap Matrix for H2: " << endl;
    vector<BasisFunction> basis_set = H2_ao.basis_set;
    
// arma::mat createhamiltonianEnergy(vector<BasisFunction>& basis_set);
// arma::mat fancyH(arma::mat X, arma::mat H);
// arma::mat createX(arma::mat overlap_matrix);
// arma::mat MO_coefficients(arma::mat X, arma::mat Fancy_H);
// double calculateEnergy(arma::mat X, arma::mat Fancy_H, int num_electrons);
// double calculateHamiltonianMatrix(AO AO_object, arma::mat Overlap_matrix);
    arma::mat S = overlap_matrix(basis_set);
    S.print();
    CNDO CNDO_H2;
    // compute the gamma matrix
    cout << "Gamma Matrix for H2: " << endl;
    int natoms = H2_ao.get_natoms();
    arma::mat gamma = CNDO_H2.computeGammaMatrix(natoms, basis_set);
    gamma.print();
    
    cout << "Compute Core Hamiltonian Matrix for H: " << endl;
    vector<string> atom_types = H2_ao.get_atom_types();
    arma::mat Hcore = CNDO_H2.computeCoreHamiltonianMatrix(atom_types, basis_set);
    Hcore.print();


    return 0;
}