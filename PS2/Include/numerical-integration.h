/**
 * @file 1D-numerical-integration.h
 * @author Charis Liao (charisliao@berkeley.edu)
 * @brief This file conatins the function declaration for 1D numerical
 *  integration using the trapezoidal rule.
 * @version 0.1
 * @date 2023-09-21
 * 
 * @copyright Copyright (c) 2023
 * 
 */


#include <armadillo>
#include "gaussian.h"
#pragma once

// Forward declaration of overlapIntegral
class OverlapIntegral;

class numericalIntegration {
public:
    numericalIntegration(double lower_bound, double upper_bound, int num_steps, OverlapIntegral& overlap);
    double trapezoidalRule();

private:
    double lower_bound;
    double upper_bound;
    int num_steps;
    OverlapIntegral& overlap; 
};
