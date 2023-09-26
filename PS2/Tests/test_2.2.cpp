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
#include "shell.h"



int main() {

    Shell s1(arma::vec ({0,0,0}), 1.0, 0);   // center at 0, alpha = 1, angular momentum = 0
    Shell s2(arma::vec ({1,0,0}), 1.0, 0);   // center at 1, alpha = 1, angular momentum = 0
    Shell p1(arma::vec ({0,0,0}), 1.0, 1);   // center at 0, alpha = 1, angular momentum = 1
    Shell p2(arma::vec ({1,0,0}), 1.0, 1);   // center at 1, alpha = 1, angular momentum = 1

    ShellPair s1s1(s1, s1);
    ShellPair s1p1(s1, p1);
    ShellPair s1s2(s1, s2);
    ShellPair s1p2(s1, p2);



    std::cout << "Analytical integral test" << std::endl;
    std::cout << "--------------------------------------------------" << std::endl;
    std::cout << "Overlap integral between s1 and s1: " << s1s1.overlap_integral() << std::endl;
    std::cout << "Overlap integral between s1 and p1: " << s1p1.overlap_integral() << std::endl;
    std::cout << "Overlap integral between s1 and s2: " << s1s2.overlap_integral() << std::endl;
    std::cout << "Overlap integral between s1 and p2: " << s1p2.overlap_integral() << std::endl;
    std::cout << std::endl;
  
    return 0;

}