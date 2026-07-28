#pragma once

#include <string>

class Analysis_mixed {
public:
    Analysis_mixed(const std::string& fileListPath);
    void Run();

private:
    std::string fileList;
    void LoadLibraries();
};