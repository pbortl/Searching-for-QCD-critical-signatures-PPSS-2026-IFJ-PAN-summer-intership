#pragma once
#include <string>
#include "HistogramManager.h"

class Plotter {
public:
    
    static void DrawAndSaveAll(HistogramManager& hists, const std::string& baseFileName, const std::vector<double>& cent_limits, const std::vector<int>& cent_events);
};
