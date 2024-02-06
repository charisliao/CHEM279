#include <iostream> 
#include <armadillo>

int main() {

    // caculate S and H for H2 
    arma::mat S;
    S << 1.0 << 0.6559 << arma::endr 
    << 0.6559 << 1.0 << arma::endr;

    arma::mat H;
    H << -13.6 << -15.705 << arma::endr
    << -15.705 << -13.6 << arma::endr;

    // symmetric orthogoalization 
    // use armadillo! 

    // diagonalize overlap matrix -> gives us eigenvectors 
    // (matrix U) and the eigenvalues (vector s) 
    arma::vec s;
    arma::mat U;

    // eigen decomposition of a symmetric matrix 
    arma::eig_sym(s, U, S);

    s.print("S Eigenvalues: ");
    U.print("S Eigenvectors: ");

    // get inverse sqrt of eigenvalues and put them into a diagonal matrix 
    for (int i = 0; i <s.size(); i++) {
        s(i) = 1.0 / std::sqrt(s(i));
    }


    arma::mat s_diag_mat = arma::diagmat(s);
    arma::mat X = U * s_diag_mat * U.t();
    X.print("X matrix: ");

    arma::mat test = arma::trans(X) * S * X;
    test.print("Unit matrix?  ");

    // transform H to H_prime = transpose(X) * H * X
    // get eigenvalues and eigenvectors of H_prime form eig_sym(eigenvalues, eigenvectors, H_prime)
    // eigenvalues -> orbital energies, add up the energies of the occupied orbitals (lowest energy)
    // transform back the eigenvectors to get the MOs: C = X*C_prime (C_prime is eigenvectors from eig_sym)
    return 0;
}