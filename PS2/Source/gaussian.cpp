
// Create a gaussian class that allows users to define center, alpha, and angular momentum 

#include "gaussian.h"
#include <cmath>



// Constructor
Gaussian::Gaussian(double center, double alpha, double angular_momentum) {
    this->center = center; 
    this->alpha = alpha; 
    this->angular_momentum = angular_momentum; 
}

// Gaussian function 
double Gaussian::gaussian(double point) {
    double r = point - center;        // distance between two points
    double result = exp(-alpha * r * r);
    return pow(r, angular_momentum) * result;
}

// overlap integral constructor
OverlapIntegral::OverlapIntegral(Gaussian g1, Gaussian g2): g1(g1), g2(g2) {
}

// overlap integral operator
double OverlapIntegral::operator()(double point) {
    return g1.gaussian(point) * g2.gaussian(point);
}





