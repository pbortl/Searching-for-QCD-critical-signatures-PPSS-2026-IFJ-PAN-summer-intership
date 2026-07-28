#include <iostream>
#include <string>
#include "Analysis.h"

int main(int argc, char** argv) {
    std::string fileListPath = "/home/p/Desktop/pathtofiles.txt";
    
    if (argc > 1) { 
        fileListPath = argv[1];
    }

    std::cout << "Starting analysis using list: " << fileListPath << std::endl;

    Analysis myAnalysis(fileListPath); // Create an instance of the Analysis class with the provided file list path
    myAnalysis.Run();

    return 0;
}