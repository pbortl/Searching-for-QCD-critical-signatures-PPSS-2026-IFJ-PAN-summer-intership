#ifndef EVENT_PLOTTER_H
#define EVENT_PLOTTER_H

#include "HistogramManager.h"
#include <string>
#include <vector>

class EventPlotter {
public:
    static void DrawAndSaveEventCuts(HistogramManager& hists, 
                                     const std::string& baseFileName, 
                                     const std::vector<double>& cent_limits, 
                                     const std::vector<int>& cent_events);
};

#endif // EVENT_PLOTTER_H