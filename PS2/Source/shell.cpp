
#include <armadillo>
#include <cmath>
#include "shell.h"

// Constructor
Shell::Shell(arma::vec center_vec, double alpha, int angular_momentum) {
    this->center_vec = center_vec;
    this->alpha = alpha;
    this->angular_momentum = angular_momentum;
}

arma::vec Shell::get_center_vec() {
    return center_vec;
}

double Shell::get_alpha() {
    return alpha;
}

arma::mat Shell::get_angular_momentum() {
    if (angular_momentum == 0) {
        return arma::zeros(1, 3);
    } else if(angular_momentum == 1) {
        arma::mat identity_matrix = arma::eye(3, 3);
        return identity_matrix;
    } else {
        return arma::mat();
    }
}

// Constructor
ShellPair::ShellPair(Shell s1, Shell s2): s1(s1), s2(s2) {
    arma::vec s1_center_vec = s1.get_center_vec();
    arma::vec s2_center_vec = s2.get_center_vec();
    double s1_alpha = s1.get_alpha();
    double s2_alpha = s2.get_alpha();
    arma::mat s1_angular_momentum = s1.get_angular_momentum();
    arma::mat s2_angular_momentum = s2.get_angular_momentum();

    this-> s1_angular_momentum = s1_angular_momentum;
    this-> s2_angular_momentum = s2_angular_momentum;
    this-> s1_center_vec = s1_center_vec;
    this-> s2_center_vec = s2_center_vec;
    this-> s1_alpha = s1_alpha;
    this-> s2_alpha = s2_alpha;
}

double ShellPair::factorial(double n) {
    if (n == 0 || n == 1) {
        return 1;
    }
    else {
        return n * factorial(n - 1);
    }
}

// double factorial function
double ShellPair::double_factorial(double n) {
    if (n == 0 || n == 1) {
        return 1;
    }
    else {
        return n * double_factorial(n - 2);
    }
}

// binomial coefficient function
double ShellPair::binomial(int m, int n) {
    return factorial(m) / (factorial(n) * factorial(m - n));
}

// Center of the Product 
arma::vec ShellPair::center_of_product() {
    return (s1_center_vec * s1_alpha + s2_center_vec * s2_alpha) / (s1_alpha + s2_alpha);
}

// exponential prefector 
double ShellPair::prefector(int dimension) {
    arma::vec r = s1_center_vec - s2_center_vec;
    return exp(-(s1_alpha * s2_alpha * r(dimension) * r(dimension)) / (s1_alpha + s2_alpha)) * pow(M_PI / (s1_alpha + s2_alpha), 0.5);
}

double ShellPair::overlap_integral_1D(int coord, int angular_momentum1, int angular_momentum2) {

    double prefactor_1D = prefector(coord);

    // check how many rows the angular momentum matrix has
    int num_rows1 = s1_angular_momentum.n_rows;
    int num_rows2 = s2_angular_momentum.n_rows;

    arma::mat double_summation_matrix;

    if (num_rows1 * num_rows2 == 1) { // both S-orbitals
        double_summation_matrix = arma::zeros(1, 1);
    } else if (num_rows1 * num_rows2 == 3) { // one S-orbital and one P-orbital
        double_summation_matrix = arma::zeros(1, 3);
    } else if (num_rows1 * num_rows2 == 9) { // both P-orbitals
        double_summation_matrix = arma::zeros(3, 3);
    } else {
        double_summation_matrix = arma::zeros(1, 1);
    }

    arma::vec center_product = center_of_product();

    for (int i = 0; i <= angular_momentum1; i++) {
        for (int j = 0; j <= angular_momentum2; j++) {
            if ((i+j) % 2 == 0) {
                result += binomial(angular_momentum1, i) * binomial(angular_momentum2, j) \
                        * (double_factorial(i+j-1) * pow(center_product(coord) - s1_center_vec(coord), angular_momentum1 - i) \
                        * pow(center_product(coord) - s2_center_vec(coord), angular_momentum2 - j)) / pow(2 * (s1_alpha + s2_alpha), (i+j)/2);
            }   
        }
    }
    return result;
    
    }

    return double_summation_matrix;

}




    










    // for (int row1 = 0; row1 < num_rows1; row1++) {
    //     for (int row2 = 0; row2 < num_rows2; row2++) {
    //         arma::rowvec s1_angular_momentum_row = s1_angular_momentum.row(row1);
    //         arma::rowvec s2_angular_momentum_row = s2_angular_momentum.row(row2);

    //         // Loop through dimensions (x, y, z)
    //         for (int dimension = 0; dimension < 3; dimension++) {
    //             double result = 0.0;

    //             for (int i = 0; i <= s1_angular_momentum_row(dimension); i++) {
    //                 for (int j = 0; j <= s2_angular_momentum_row(dimension); j++) {
    //                     if ((i + j) % 2 == 0) {
    //                         result += binomial(s1_angular_momentum_row(dimension), i) * binomial(s2_angular_momentum_row(dimension), j) *
    //                                   (double_factorial(i + j - 1) *
    //                                    pow(center_product(dimension) - s1_center_vec(dimension), s1_angular_momentum_row(dimension) - i) *
    //                                    pow(center_product(dimension) - s2_center_vec(dimension), s2_angular_momentum_row(dimension) - j)) /
    //                                   pow(2 * (s1_alpha + s2_alpha), (i + j) / 2);
    //                     }
    //                 }
    //             }

    //             // Store the result in the appropriate row and dimension of the result matrix
    //             double_summation_matrix(row1 * num_rows2 + row2, dimension) = result;
    //             // multiply by the prefector
    //             double_summation_matrix(row1 * num_rows2 + row2, dimension) *= prefector(dimension);
    //         }

    //     }

