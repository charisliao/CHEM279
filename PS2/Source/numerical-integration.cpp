/**
 * @file 1D-numerical-integration.cpp
 * @author Charis Liao (charisliao@berkeley.edu)
 * @brief 
 * @version 0.1
 * @date 2023-09-21
 * 
 * @copyright Copyright (c) 2023
 * 
 */

#include <iostream>
#include <armadillo>
#include <vector>
#include <string>
#include <cmath>
#include <functional>
#include "numerical-integration.h"
#include "gaussian.h"

// Constructor
numericalIntegration::numericalIntegration(double lower_bound, double upper_bound, int num_steps, OverlapIntegral& overlap)
    : lower_bound(lower_bound), upper_bound(upper_bound), num_steps(num_steps), overlap(overlap) {
}

double numericalIntegration::trapezoidalRule() {
    double h = (upper_bound - lower_bound) / num_steps; 
    double sum = 0.0; 
    for (int i=1; i < num_steps; i++) {
        sum += overlap(lower_bound + i*h);
    }
    sum += 0.5 * (overlap(lower_bound) + overlap(upper_bound));
    return sum * h; 
}










