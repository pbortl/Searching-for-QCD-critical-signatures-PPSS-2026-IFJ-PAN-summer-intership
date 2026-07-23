#pragma once

#include <string>

class Analysis {
private:
    std::string fileList;
    void LoadLibraries();

public:
    
    Analysis(const std::string& fileListPath);

    
    void Run();
};
