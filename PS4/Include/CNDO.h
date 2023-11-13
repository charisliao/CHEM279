
#include <cmath>
#include <iostream>
#include <armadillo>
#include <string>
#include <vector>
#include <map>
#include "AO.h"
#include "utils.h"


struct CNDO_parameters {
    std::map<std::string, double> IA;
    double beta;
};


class CNDO {
    private:
        std::map<std::string, CNDO_parameters> CNDO_para_map;
        

    public: 
        
        CNDO();
        double getBeta(const string& element);

        double getIA(const std::string& element, const std::string& orbitalType) const;
        double gamma(BasisFunction& basisFunction1, BasisFunction& basisFunction2);
        arma::mat computeGammaMatrix(int natoms, vector<BasisFunction>& basis_set);
        double hartree_to_ev  = 27.211396641308;
        double diagCoreH(double gammaAA, BasisFunction& basisFunctionA, double ZB_gamma_AB);
        double computeOffDiagMat(int i, int j, double Beta_A, double Beta_B, arma::mat overlapMatrix);
        arma::mat computeCoreHamiltonianMatrix(vector<string> atom_types, vector<BasisFunction>& basis_set);
        std::string extractElement(const std::string& AO_type);

        std::string extractOrbital(const std::string& AO_type);

};