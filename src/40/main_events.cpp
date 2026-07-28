#include <iostream>
#include <string>
#include "Analysis_Events.h"

int main(int argc, char** argv) {
    std::string fileListPath = "/home/p/Desktop/pathtofiles.txt";
    
    if (argc > 1) { 
        fileListPath = argv[1];
    }

    std::cout << "Starting FAST EVENT-ONLY analysis using list: " << fileListPath << std::endl;

    Analysis_Events myAnalysis(fileListPath);
    myAnalysis.Run();

    return 0;
}