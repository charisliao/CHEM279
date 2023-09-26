
#include <armadillo>
#include <cmath>

#pragma once



class Shell {
    private:
        arma::vec center_vec;
        double alpha;
        int angular_momentum;
    public:
        Shell(arma::vec center_vec, double alpha, int angular_momentum);

        arma::vec get_center_vec();
        double get_alpha();
        arma::mat get_angular_momentum();
};

class ShellPair {

    private:
        Shell s1;
        Shell s2;
        arma::vec s1_center_vec;
        arma::vec s2_center_vec;
        double s1_alpha;
        double s2_alpha;
        arma::mat s1_angular_momentum;
        arma::mat s2_angular_momentum;

    public:

        ShellPair(Shell s1, Shell s2);
        double factorial(double n);
        double double_factorial(double n);
        double binomial(int m, int n);
        arma::vec center_of_product();
        double prefector(int dimension);
        double overlap_integral_1D(int dimension, int angular_momentum1, int angular_momentum2);
        arma::mat overlap_integral_3D();
};

