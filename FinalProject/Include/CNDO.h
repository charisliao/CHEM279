
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
        double offDiagonalFockElement(double BetaA, double BetaB, double overlap, double p, double gammaAB);
        double diagonalFockElement(double IA, double ZA, double gammaAA, double pTotalA, double pMuMu, double pZBgammaAB);
        arma::mat computeFockMatrix(vector<string> atom_types, vector<BasisFunction>& basis_set, arma::mat& overlap_mat, arma::mat& Hcore_mat, arma::mat& p, arma::vec& Ptotal);
        arma::vec updatePtotal(vector<string> atom_types, vector<BasisFunction> basis_set, arma::mat& densityMat);
        double calculateNuclearRepulsionEnergy(vector<string> atom_types, AO AO_object);
        double calculateElectronEnergy(arma::mat& densityMat, arma::mat& Hcore_mat, arma::mat& Fock_mat);
        void updateDensityMatrix(AO AO_object, std::string output_file_name);

        // optimized functions for final project 
        arma::mat computeCoreHamiltonianMatrixOptimized(arma::mat& overlap_mat, std::map<int, BasisFunction> basis_map, vector<string> atom_types, vector<BasisFunction>& basis_set);
        arma::mat computeFockMatrixOptimized(std::map<int, BasisFunction> basis_map, vector<string> atom_types, vector<BasisFunction>& basis_set, arma::mat& overlap_mat, arma::mat& Hcore_mat, arma::mat& p, arma::vec& Ptotal);
        arma::mat updateDensityMatrixOptimized(std::map<int, BasisFunction> basis_map, AO AO_object, arma::mat overlap_mat, arma::mat Hcore_mat, arma::mat gamma_mat, vector<string> atom_types, vector<BasisFunction> basis_set);
};  