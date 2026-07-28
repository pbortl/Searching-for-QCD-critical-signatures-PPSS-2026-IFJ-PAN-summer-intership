#include <iostream>
#include <string>
#include "Analysis_mixed.h"

int main(int argc, char** argv) {
    std::string fileListPath = "/home/p/Desktop/pathtofiles.txt";
    
    if (argc > 1) {
        fileListPath = argv[1];
    }

    std::cout << "Starting MIXED EVENTS analysis using list: " << fileListPath << std::endl;

    Analysis_mixed myAnalysis(fileListPath);
    myAnalysis.Run();

    return 0;
}