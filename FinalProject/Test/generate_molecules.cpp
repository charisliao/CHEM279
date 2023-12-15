
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include "utils.h"


int main(){
    generateMoleculeFile(50, 50, "C50H50.txt");
    generateMoleculeFile(250, 250, "C250H250.txt");
    generateMoleculeFile(500, 500, "C500H500.txt");
    return 0;
}