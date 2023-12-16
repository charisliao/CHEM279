
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

    cout << "Overlap Matrix for H: " << endl;
    vector<BasisFunction> basis_set = H_ao.basis_set;
    
// arma::mat createhamiltonianEnergy(vector<BasisFunction>& basis_set);
// arma::mat fancyH(arma::mat X, arma::mat H);
// arma::mat createX(arma::mat overlap_matrix);
// arma::mat MO_coefficients(arma::mat X, arma::mat Fancy_H);
// double calculateEnergy(arma::mat X, arma::mat Fancy_H, int num_electrons);
// double calculateHamiltonianMatrix(AO AO_object, arma::mat Overlap_matrix);
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
    arma::mat Hcore = CNDO_H.computeCoreHamiltonianMatrix(atom_types, basis_set);
    Hcore.print();
    
    CNDO_H.updateDensityMatrix(H_ao, "myH.txt");



    return 0;
}