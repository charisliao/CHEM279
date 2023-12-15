#include <armadillo> 
#include <cmath>
#include <chrono>
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>

double factorial(double n); 
double double_factorial(double n);
double binomial(int m, int n);
double center_of_product(double center_vec1, double alpha1, double center_vec2, double alpha2);
double prefactor(double center_vec1, double alpha1, double center_vec2, double alpha2);
double I2e_pG(arma::rowvec &Ra, arma::rowvec &Rb, double sigmaA, double sigmaB);
double measureExecutionTime(const std::function<void()>& func, const std::string& functionName);
void generateMoleculeFile(int num_C, int num_H, const std::string& filename);