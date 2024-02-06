
#include "utils.h"
#include <armadillo> 
#include <cmath>



double factorial(double n) {
    if (n == 0 || n == 1) {
        return 1;
    }
    else {
        return n * factorial(n - 1);
    }
}

// double factorial function
double double_factorial(double n) {
    if (n <= 1) {
        return 1;
    }
    else {
        return n * double_factorial(n - 2);
    }
}

// binomial coefficient function
double binomial(int m, int n) {
    if (m < n) {
        return 0.0;
    }
    return factorial(m) / (factorial(n) * factorial(m - n));
}

// Center of the Product 
double center_of_product(double center_vec1, double alpha1, double center_vec2, double alpha2) {
    return (center_vec1 * alpha1 + center_vec2 * alpha2) / (alpha1 + alpha2);
}

// exponential prefector 
double prefactor(double center_vec1, double alpha1, double center_vec2, double alpha2) {
    double r = center_vec1 - center_vec2;
    return exp(-(alpha1 * alpha2 * r * r) / (alpha1 + alpha2)) * pow(M_PI / (alpha1 + alpha2), 0.5);
}

