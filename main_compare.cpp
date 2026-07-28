#include <iostream>
#include <string>
#include "Plotter.h"

int main(int argc, char** argv) {
    if (argc < 4) {
        std::cout << "Usage: ./draw_delta <data_F2.txt> <mixed_F2.txt> <output_basename>" << std::endl;
        return 1;
    }

    std::string dataFile = argv[1];
    std::string mixedFile = argv[2];
    std::string outBase = argv[3];

    std::cout << "Calculating Delta F2..." << std::endl;
    std::cout << "Data:  " << dataFile << std::endl;
    std::cout << "Mixed: " << mixedFile << std::endl;

    Plotter::DrawDeltaF2(dataFile, mixedFile, outBase);

    return 0;
}