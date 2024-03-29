/**
 * @file gaussian.h
 * @author Charis Liao (charisliao@berkeley.edu)
 * @brief The header file for the Gaussian class 
 * @version 0.1
 * @date 2023-09-25
 * 
 * @copyright Copyright (c) 2023
 * 
 */

#pragma once
#include <armadillo>

class Gaussian {
    public:
        double center; 
        double alpha; 
        double angular_momentum; 
        arma::vec center_vec;
        arma::mat angular_momentum_mat;

        Gaussian(double center, double alpha, double angular_momentum);

        double gaussian(double point);


};

class OverlapIntegral {
    public:
        OverlapIntegral(Gaussian g1, Gaussian g2);
        double operator()(double point);
    
    private: 
        Gaussian g1; 
        Gaussian g2;
};