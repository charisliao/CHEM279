#include "CNDO.h"
#include <stdio.h>
#include <math.h>
#include <cassert>
#include <unordered_map>
#include <armadillo>
#include <iostream>
#include <vector>
#include <iomanip>

using namespace std;

double hartree_to_ev  = 27.211396641308;
struct CNDO_para {
    std::map<std::string, double> IA;
    double beta;
} C_para = {{{"s", 14.051}, {"p", 5.572}}, 21.}, H_para{{{"s", 7.176}}, 9.};
CNDO_para N_para = {{{"s", 19.316}, {"p", 7.275}}, 25.};
CNDO_para O_para = {{{"s", 25.390}, {"p", 9.111}}, 31.};
CNDO_para F_para = {{{"s", 32.272}, {"p", 11.080}}, 39.};

std::map<std::string, CNDO_para> CNDO_para_map = {{"H", H_para}, {"C", C_para}, {"N", N_para}, {"O", O_para}, {"F", F_para}};


double Eval_2eI_sAO(AO &ao1, AO &ao2) {
    arma::vec la = ao1.get_lmn(); lb = ao2.get_lmn();
    if(!(arma::accu(la) == 0 && arma::accu(lb) == 0)) {
        std::cout << "Only s orbitals allowed" << std::endl;
    }
    //information from basis set
    arma::vec da=ao1.get_d_coe(), db = ao2.get_d_coe();
    arma::vec alphaa = ao1.get_alpha(), alphab = ao2.get_alpha();
    arma::vec Ra = ao1.get_R0(), Rb = ao2.get_R0();

    int len = ao1.get_len();

    double sum = 0;
    for (int k1=0; k1 < len; k1++) {
        for (int k2=0; k2 < len; k2++) {
            double sigmaA = 1.0/(alphaa(k1) + alphaa(k2)); //eq 3.10

            for (int j1 =0; j1 < len; j1++) {
                for (int j2 = 0; j2 < len; j2++) {
                    double sigmaB = 1.0/(alphab(j1) + alphab(j2)); //eq 3.10
                    double I2e = I2e_pG(Ra, Rb, sigmaA, sigmaB); //eq 3.14
                    sum += da(k1) * da(k2) * db(j1) * db(j2) * I2e; //eq 3.13
                }
            }
        }

    }
    return sum; 
}

// she put this one in utils.cpp (to evaluate 3.14)
double I2e_pG(arma::vec &Ra, arma::vec &Rb, double sigmaA, double sigmaB) {
    double U = pow(M_PI * sigmaA, 1.5) * pow(M_Pi * sigmaB, 1.5); //eq 3.8 and 3.11
    double V2 = 1.0 / (sigmaA+sigmaB); //eq. 3.9 

    double Rd = arma::norm(Ra - Rb, 2);
    
    if (Rd == 0.0) {
        //use eq 3.15
        return U * sqrt(2*V2) * sqrt(2/M_PI);
    }
    //eq. 3.14 need sqrt T
    double srT = sqrt(V2) * Rd; 

    // eq 3.14
    double result = U / Rd * erf(srT); 
    return result; 

}

//H1 C1 C2 H2
//gamma H1H1 "Gamma AA", gamma H1C1, gamma H1C2, gamma H1H2


arma::mat CNDO::computeCoreHamiltonianMatrix(vector<string> atom_types, vector<BasisFunction>& basis_set) {
    arma::mat H = arma::zeros(basis_set.size(), basis_set.size());
    arma::mat gamma_mat = computeGammaMatrix(atom_types.size(), basis_set);
    arma::mat overlap_mat = overlap_matrix(basis_set);
    int atom_index = 0;

    for (int i = 0; i < basis_set.size(); i++) {

        for (int j = 0; j < basis_set.size(); j++) {
            
            if (i == j) {
                double gammaAA = 0.0;
                BasisFunction basisFunctionA = basis_set[i];

                // Loop through the basis functions associated with the current atom
                for (int basis = 0; basis < basis_set.size(); basis++) {
                    if (basis_set[basis].atom_index == atom_index) {
                        if (basis_set[basis].AO_type.find("s") == std::string::npos) {
                            // find the basis function that has the same element but in S orbital 
                            for (int basisS = 0; basisS < basis_set.size(); basisS++) {
                                if (basis_set[basisS].atom_index == atom_index &&
                                    basis_set[basisS].AO_type.find(extractElement(basis_set[basis].AO_type)) != std::string::npos) {
                                    basisFunctionA = basis_set[basisS];
                                }
                            }
                        }
                        gammaAA += gamma(basisFunctionA, basisFunctionA);
                    }
                }

                double ZB_gamma_AB = 0.0;
                num_atoms = atom_types.size()-1;
                for (int k = 0; k < basis_set.size(); k++) {
                    if (i == k) {
                        continue;
                    } else if (basis_set[k].atom_index == atom_index) {
                        
                        BasisFunction basisFunctionB = basis_set[k];
                        if (basis_set[k].AO_type.find("s") == std::string::npos) {
                            // find the basis function that has the same element but in S orbital
                            for (int basisS = 0; basisS < basis_set.size(); basisS++) {
                                if (basis_set[basisS].atom_index == atom_index &&
                                    basis_set[basisS].AO_type.find(extractElement(basis_set[k].AO_type)) != std::string::npos) {
                                    basisFunctionB = basis_set[basisS];
                                }
                            }
                            double ZB = basisFunctionB.valence_e;
                            double gammaAB = gamma(basisFunctionA, basisFunctionB);
                            ZB_gamma_AB += ZB * gammaAB;
                        }
                    }
                }

                H(i, j) = diagCoreH(gammaAA, basisFunctionA, ZB_gamma_AB);
                cout << "Current atom index: " << atom_index <<  " Current atom type: " << atom_types[atom_index] << endl;

            } else {
                double Beta_A = getBeta(extractElement(basis_set[i].AO_type));
                double Beta_B = getBeta(extractElement(basis_set[j].AO_type));
                H(i, j) = computeOffDiagMat(i, j, Beta_A, Beta_B, overlap_mat);
            }
        }
    }
    return H;
}
