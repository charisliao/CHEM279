/**
 * @file main.cpp
 * @author Charis Liao 
 * @brief 
 * @version 0.1
 * @date 2023-09-23
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


void NumericalIntegrationTest(int num_steps) {

    Gaussian s1(0.0, 1.0, 0.0);   // center at 0, alpha = 1, angular momentum = 0
    Gaussian s2(1.0, 1.0, 0.0);   // center at 1, alpha = 1, angular momentum = 0
    Gaussian p1(0.0, 1.0, 1.0);   // center at 0, alpha = 1, angular momentum = 1
    Gaussian p2(1.0, 1.0, 1.0);   // center at 1, alpha = 1, angular momentum = 1

    OverlapIntegral overlap_s1s1(s1, s1);
    OverlapIntegral overlap_s1p1(s1, p1);
    OverlapIntegral overlap_s1s2(s1, s2);
    OverlapIntegral overlap_s1p2(s1, p2);

    numericalIntegration integral_s1s1(-5.0, 5.0, num_steps, overlap_s1s1);
    numericalIntegration integral_s1p1(-5.0, 5.0, num_steps, overlap_s1p1);
    numericalIntegration integral_s1s2(-5.0, 5.0, num_steps, overlap_s1s2);
    numericalIntegration integral_s1p2(-5.0, 5.0, num_steps, overlap_s1p2);


    std::cout << "1D overalp integral test with " << num_steps << " steps" << std::endl;
    std::cout << "--------------------------------------------------" << std::endl;
    std::cout << "Overlap integral between s1 and s1: " << integral_s1s1.trapezoidalRule() << std::endl;
    std::cout << "Overlap integral between s1 and p1: " << integral_s1p1.trapezoidalRule() << std::endl;
    std::cout << "Overlap integral between s1 and s2: " << integral_s1s2.trapezoidalRule() << std::endl;
    std::cout << "Overlap integral between s1 and p2: " << integral_s1p2.trapezoidalRule() << std::endl;
    std::cout << std::endl;

}

int main() {

    std::cout << "Start of integration test with different steps" << std::endl;

    NumericalIntegrationTest(5);
    NumericalIntegrationTest(10);
    NumericalIntegrationTest(1000);
    NumericalIntegrationTest(10000);

    return 0;

}