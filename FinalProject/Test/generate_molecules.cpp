
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <iomanip>
#include <filesystem>
#include "utils.h"


void generateMoleculeFiles(const std::string& folderPath) {
    // Create the folder if it doesn't exist
    std::filesystem::create_directories(folderPath);

    for (int numAtoms = 1; numAtoms <= 335; ++numAtoms) {
        std::ostringstream filenameStream;
        filenameStream << folderPath << "/" << numAtoms << ".txt";
        std::string filename = filenameStream.str();

        generateMoleculeFile(numAtoms, filename);
    }
}

int main(){
    generateMoleculeFiles("molecules");
    return 0;
}