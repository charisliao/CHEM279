/**
 * @file AO.cpp
 * @author Charis Liao (charisliao@berkeley.edu)
 * @brief 
 * @version 0.1
 * @date 2023-10-15
 * 
 * @copyright Copyright (c) 2023
 * 
 */


#include <armadillo>
#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <iomanip>
#include <cmath>
#include <algorithm>
#include <numeric>
#include "AO.h"
#include "utils.h"

using namespace std;  



//constructor
AO::AO(const char *filename) {

    ifstream infile(filename);
    vector<BasisFunction> basis_set;
    if (!infile.is_open()) {
        cerr << "Error: Unable to open file " << filename << endl;
    // Add additional error handling if needed
    } else {
    // Proceed with reading the file

        infile >> natoms;
        // skip next number
        infile.ignore(256, '\n');
        atomic_numbers.resize(natoms); 
        coord.resize(natoms, 3);  // this is an armadiilo matrix

        for (int i = 0; i < natoms; i++){
            infile >> atomic_numbers[i] >> coord(i, 0) >> coord(i, 1) >> coord(i, 2);
            // cout << "Atomic number: " << atomic_numbers[i] << endl;
            if (atomic_numbers[i] == 1) {
                num_hydrogen += 1;
                BasisFunction hydrogen_basis = createHydrogenBasis(coord.row(i));
                // cout << "Hydrogen basis: " << hydrogen_basis.center << endl;
                basis_set.push_back(hydrogen_basis);
            } else if (atomic_numbers[i] == 6) {
                num_carbon += 1;
                vector<BasisFunction> carbon_basis = createCarbonBasis(coord.row(i));
                basis_set.insert(basis_set.end(), carbon_basis.begin(), carbon_basis.end());
            }
        }
        num_basis = count_basis(num_carbon, num_hydrogen);
        num_electrons = count_electrons(num_carbon, num_hydrogen);
    }

    updateBasisSet(basis_set);
    this->basis_set = basis_set;
    
    

    infile.close();
}

AO::~AO() {
    // destructor
}   

int AO::count_basis(int num_carbon, int num_hydrogen) {
    int num_basis = 0;
    num_basis = 4 * num_carbon + num_hydrogen;
    return num_basis;
}

int AO::count_electrons(int num_carbon, int num_hydrogen) {
    
    // throw error if num_electrons is not an integer
    if ((2 * num_carbon + num_hydrogen) % 2 != 0) {
        cerr << "Error: Number of electrons is not an integer" << endl;
        throw std::runtime_error("Number of electrons is not an integer");
    } 
    int num_electrons = 2 * num_carbon + (num_hydrogen / 2);
    
    return num_electrons;
}



BasisFunction AO::createHydrogenBasis(arma::rowvec center) {
    return {"H_s", {3.42525091, 0.62391373, 0.16885540}, {0.15432897, 0.53532814, 0.44463454}, center, {0, 0, 0}, {0.0, 0.0, 0.0}};
}

vector<BasisFunction> AO::createCarbonBasis(arma::rowvec center) {
    vector<BasisFunction> carbonBasis; 
    carbonBasis.push_back({"C_s", {2.94124940, 0.68348310, 0.22228990}, {-0.09996723, 0.39951283, 0.70011547}, center, {0, 0, 0}, {0.0, 0.0, 0.0}});
    carbonBasis.push_back({"C_px", {2.94124940, 0.68348310, 0.22228990}, {0.15591627, 0.60768372, 0.39195739}, center, {1, 0, 0}, {0.0, 0.0, 0.0}});
    carbonBasis.push_back({"C_py", {2.94124940, 0.68348310, 0.22228990}, {0.15591627, 0.60768372, 0.39195739}, center, {0, 1, 0}, {0.0, 0.0, 0.0}});
    carbonBasis.push_back({"C_px", {2.94124940, 0.68348310, 0.22228990}, {0.15591627, 0.60768372, 0.39195739}, center, {0, 0, 1}, {0.0, 0.0, 0.0}});
    return carbonBasis;
}



void AO::updateNormalizationConstants(BasisFunction& basisFunction) {

    double omega1 = overlap_integral_3D(basisFunction.center, basisFunction.center, basisFunction.exponents(0), basisFunction.exponents(0), basisFunction.lmn, basisFunction.lmn);
    double omega2 = overlap_integral_3D(basisFunction.center, basisFunction.center, basisFunction.exponents(1), basisFunction.exponents(1), basisFunction.lmn, basisFunction.lmn);
    double omega3 = overlap_integral_3D(basisFunction.center, basisFunction.center, basisFunction.exponents(2), basisFunction.exponents(2), basisFunction.lmn, basisFunction.lmn);
    
    // cout << "omega1: " << omega1 << endl;
    // cout << "omega2: " << omega2 << endl;
    // cout << "omega3: " << omega3 << endl;

    // cout << "omega1_"



    basisFunction.normalization_constants(0) = pow(omega1, -0.5);
    basisFunction.normalization_constants(1) = pow(omega2, -0.5);
    basisFunction.normalization_constants(2) = pow(omega3, -0.5);

    // cout << "normalization_constants: " << basisFunction.normalization_constants << endl;

}

void AO::updateBasisSet(vector<BasisFunction>& basis_set) {
    for (int i = 0; i < basis_set.size(); i++) {
        updateNormalizationConstants(basis_set[i]);
    }
}



double overlap_integral_1D(double center1, double center2, double alpha1, double alpha2, double l1, double l2) {
   
    double prefactor_1D = prefactor(center1, alpha1, center2, alpha2);
    double center_product = center_of_product(center1, alpha1, center2, alpha2);
    double double_sum = 0.0;
    
    for (int i = 0; i <= l1; i++) {
        for (int j = 0; j <= l2; j++) {
            if ((i+j) % 2 == 0) {
                double_sum += binomial(l1, i) * binomial(l2, j) * ((double_factorial(i+j-1) * pow(center_product - center1, l1 - i) * pow(center_product - center2, l2 - j)) / pow(2 * (alpha1 + alpha2), double(i+j)/2));
            }   
        }
    }
    
    return double_sum * prefactor_1D;
}


double overlap_integral_3D(arma::rowvec center1, arma::rowvec center2, double alpha1, double alpha2, arma::vec lmn1, arma::vec lmn2) {
    
    double overlap_orbital = overlap_integral_1D(center1(0), center2(0), alpha1, alpha2, lmn1(0), lmn2(0)) \
                            * overlap_integral_1D(center1(1), center2(1), alpha1, alpha2, lmn1(1), lmn2(1)) \
                            * overlap_integral_1D(center1(2), center2(2), alpha1, alpha2, lmn1(2), lmn2(2));
    
    // cout << "used in normalization: " << overlap_orbital << endl;
    return overlap_orbital;

}
double unnormOverlapIntegral(BasisFunction& basisFunction1, BasisFunction& basisFunction2, int i, int j) {
    
    double overlap_orbital = overlap_integral_1D(basisFunction1.center(0), basisFunction2.center(0), basisFunction1.exponents(i), basisFunction2.exponents(j), basisFunction1.lmn(0), basisFunction2.lmn(0)) \
                            * overlap_integral_1D(basisFunction1.center(1), basisFunction2.center(1), basisFunction1.exponents(i), basisFunction2.exponents(j), basisFunction1.lmn(1), basisFunction2.lmn(1)) \
                            * overlap_integral_1D(basisFunction1.center(2), basisFunction2.center(2), basisFunction1.exponents(i), basisFunction2.exponents(j), basisFunction1.lmn(2), basisFunction2.lmn(2));
    // cout << "Unnormalized overlap_orbital: " << overlap_orbital << endl;
    return overlap_orbital;
}

double overlapIntegral(BasisFunction& basisFunction1, BasisFunction& basisFunction2) {
    
    double overlap_orbital = 0.0;
    // double unnormalized_overlap = unnormOverlapIntegral(basisFunction1, basisFunction2);
    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
            // cout << "normalization_constants(i): " << basisFunction1.normalization_constants(i) << endl;
            // cout << "normalization_constants(j): " << basisFunction2.normalization_constants(j) << endl;
            // cout << "contracted_coefficients(i): " << basisFunction1.contracted_coefficients(i) << endl;
            // cout << "contracted_coefficients(j): " << basisFunction2.contracted_coefficients(j) << endl;
            // cout << "unnormalized_overlap: " << unnormalized_overlap << endl;

            overlap_orbital += basisFunction1.normalization_constants(i) * basisFunction2.normalization_constants(j) * 
                             basisFunction1.contracted_coefficients(i) * basisFunction2.contracted_coefficients(j) * unnormOverlapIntegral(basisFunction1, basisFunction2, i, j);
        }
    }


    // cout << "Normalized overlap_orbital: " << overlap_orbital << endl;
    
    return overlap_orbital;
}

arma::mat overlap_matrix(vector<BasisFunction>& basis_set) {
    arma::mat S = arma::zeros(basis_set.size(), basis_set.size());
    for (int i = 0; i < basis_set.size(); i++) {
        for (int j = 0; j < basis_set.size(); j++) {
            // cout << "inside overlap matrix call overlap integral" << endl;
            // cout << "basis_set[i]: " << basis_set[i].normalization_constants<< endl;
            // cout << "basis_set[j]: " << basis_set[j].normalization_constants << endl;
            S(i, j) = overlapIntegral(basis_set[i], basis_set[j]);
            // cout << "inside overlap matrix call overlap integral. DONE!" << endl;
        }
    }
    return S;
}




    