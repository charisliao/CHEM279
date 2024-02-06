
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
    
    AO N2_ao("N2.txt");

    N2_ao.set_p(4);
    N2_ao.set_q(1);

    cout << "Overlap Matrix for N2: " << endl;
    vector<BasisFunction> basis_set = N2_ao.basis_set;
    
// arma::mat createhamiltonianEnergy(vector<BasisFunction>& basis_set);
// arma::mat fancyH(arma::mat X, arma::mat H);
// arma::mat createX(arma::mat overlap_matrix);
// arma::mat MO_coefficients(arma::mat X, arma::mat Fancy_H);
// double calculateEnergy(arma::mat X, arma::mat Fancy_H, int num_electrons);
// double calculateHamiltonianMatrix(AO AO_object, arma::mat Overlap_matrix);
    arma::mat S = overlap_matrix(basis_set);
    S.print();

    // create a CNDO instance
    CNDO CNDO_N2;
    // compute the gamma matrix
    cout << "Compute Gamma Matrix for N2: " << endl;
    int natoms = N2_ao.get_natoms();
    
    cout << "natoms: " << natoms << endl;
    arma::mat gamma = CNDO_N2.computeGammaMatrix(natoms, basis_set);
    gamma.print();

    cout << "Compute Core Hamiltonian Matrix for N2: " << endl;
    vector<string> atom_types = N2_ao.get_atom_types();
    arma::mat Hcore = CNDO_N2.computeCoreHamiltonianMatrix(atom_types, basis_set);
    Hcore.print();

    CNDO_N2.updateDensityMatrix(N2_ao, "myN2.txt");

    return 0;
}