/**
 * @file hamiltonian.h
 * @author Charis Liao (charisliao@berkeley.edu)
 * @brief 
 * @version 0.1
 * @date 2023-10-17
 * 
 * @copyright Copyright (c) 2023
 * 
 */

const double K = 1.75; 

arma::mat createhamiltonianEnergy(vector<BasisFunction>& basis_set);
arma::mat fancyH(arma::mat X, arma::mat H);
arma::mat createX(arma::mat overlap_matrix);
arma::mat MO_coefficients(arma::mat X, arma::mat Fancy_H);
double calculateEnergy(arma::mat X, arma::mat Fancy_H, int num_electrons);
double calculateHamiltonianEnergy(AO AO_object, arma::mat Overlap_matrix);


