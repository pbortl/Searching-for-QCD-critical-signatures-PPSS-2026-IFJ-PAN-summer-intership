#ifndef ANALYSIS_EVENTS_H
#define ANALYSIS_EVENTS_H

#include <string>

class Analysis_Events {
public:
    Analysis_Events(const std::string& fileListPath);
    void Run();

private:
    void LoadLibraries();
    std::string fileList;
};

#endif // ANALYSIS_EVENTS_H