#include "AO.h"
#include <stdexcept>
#include <fstream>
#include <stdio.h>
#include <math.h>
#include <sstream>
#include <cassert>
#include <string>
#include "util.h"
#include <vector>

using namespace std;

//AO functions
void AO::printinfo() {
    printf("This AO infor: %s, R( %1.2f, %1.2f, %1.2f), with angular momentum: %lld %lld %lld\n", label.c_str().
    R0(0), R0(1), R0(2), lmn(0), lmn(1), lmn(2));

    d_coe.print("d_coe");
    alpha.print("alpha");
}

AO::AO(arma::vec &R0_input, arma::vec &alpha_input, arma::vec &d_input, arma::uvec &lmn_input, string label_input):
R0(R0_input), alpha(alpha_input), d_coe(d_input), lmn(lmn_input), label(label_input) {
    assert(R0.n_elem == 3);
    assert(lmn.n_elem == 3);
    len = alpha.n_elem;
    assert(d_coe.n_elem == len);
    for (size_t k = 0; k < len; k++) {
        double Overlap_Self = Overlap_3d(R0, R0, alpha(k), alpha(k), lmn, lmn);
        d_coe(k) /= sqrt(Overlap_Self);
    }
}



double CNDO::getEnergy() {
    arma::mat Ptotal = Pa + Pb;
    Ee = arma::dot(Pa, Ga) / 2. + arma::dot(Pb, Gb) / 2.;
    Ee += arma::dot(Ptotal, Hcore);
    Etotal = Ee + Ec;
    cout << "Nuclear Repulsion Energy is " << Ec << " eV." <<endl;
    cout << "Electronic Energy is " << Ee << " eV." << endl;
    return Etotal;
}


/**
 * @brief Calculate the core Hamiltonian matrix.
 * 
 * This function calculates the core Hamiltonian matrix using the semi-empirical
 * parameters and the overlap matrix.
 */
mat CNDO::calcHCoreMat() {
    mat hCoreMat = arma::zeros(molecule.nBasisFunctions, molecule.nBasisFunctions);
    int mu = 0;      // Index for all atomA AOs
    // Loop over all atoms in molecule (atomA)
    for (int A = 0; A < molecule.nAtoms; A++) {
        string chemSymA = molecule.atomicSymbols[A];
        double ZA = molecule.atomValences(A);
        double gammaAA = gammaMatrix(A, A);
        int numAOs_A = calcNumAOs(chemSymA);

        // Loop over all AOs in atomA
        for (int i = 0; i < numAOs_A; i++) {
            // Calculate the diagonal elements of the matrix
            hCoreMat(mu, mu) = -diagCNDOPara[chemSymA][molecule.basisFunctionsList[mu].AO_type] \
                               - (ZA - 0.5) * gammaAA;

            int nu = 0;   // Index for all atomB AOs
            // Loop over all atoms in molecule (atomB)
            for (int B = 0; B < molecule.nAtoms; B++) {
                string chemSymB = molecule.atomicSymbols[B];
                double ZB = molecule.atomValences(B);
                double gammaAB = gammaMatrix(A, B);
                int numAOs_B = calcNumAOs(chemSymB);

                // Subtract from the diagonal elements of the matrix when A != B
                if (A != B) {
                    hCoreMat(mu, mu) -= ZB * gammaAB;
                }

                // Loop over all AOs in atomB
                for (int j = 0; j < numAOs_B; j++) {
                    if (mu != nu) {
                        // Calculate the off-diagonal elements of the matrix
                        hCoreMat(mu, nu) = -(offDiagCNDOPara[chemSymA] + offDiagCNDOPara[chemSymB]) \
                                            / 2.0 * overlapMatrix(mu, nu);
                    }
                    nu++;   // Increment with each AO in atomB
                }
            }
            mu++;  // Increment with each AO in atomA
        }
    }
    return hCoreMat;
}