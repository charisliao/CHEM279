/**
 * @file golden-section.cpp
 * @author Charis Liao 
 * @brief 
 * @version 0.1
 * @date 2023-09-21
 * 
 * @copyright Copyright (c) 2023
 * 
 */

// assume we have a funciton: LJ energy, we have a function for analytical gradient. 
// assume we have the molecule coordinates and constants 

// we want to optimize our geometry with golden section search 

void goldenSectionSearch(arma::mat cood, double a, double b, double tol) {
    const double goldenRatio = (sqrt(5.0) - 1.0) / 2.0; 

    while ((b - a) > tol) { //eg. 10e-5 Angstrom

        double x1 = b - goldenRatio*(b-a);
        double x2 = a - goldenRatio*(b-a);

        // decide on two atoms that we are looking at 
        arma::Row candidate1 = geom(0); // as in the first row of our geometry matrix
        arma::Row candidate2 = geom(1);

        // compute the gradient on both of the atoms 
        arma::Row grad1, grad2;

        computeGradient(candidate1, grad1);
        computeGradient(candidate2, grad2);

        // update the coordinates 
        for (int i=0; i < 3; i++) {
            candidate1(i) -= grad1(i) * (x1-x2);
            candidate2(i) -= grad2(i) * (x2-x1);
        }

        // after this loop, we've updated x, y, z coordinates
        // of atom 1 and atom 2 within golden ratio interval. 

        double energy1 = LJ(candidate1);
        double energy2 = LJ(candidate2);

        if (energy1 < energy2) {
            b = x2;
            geom(0) = candidate1;
        } else {
            a=x1;
            geom(0) = candidate2;
        }

    }
}