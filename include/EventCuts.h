#pragma once
#include "EventXeLaMag.h"
#include "Config.h"

class EventCuts {
public:
    // Return true if the event passes 
    static bool PassVertexZ(double vZ);
    static bool PassTracksRatio(int nTracks, int nTracksInFit);
    static bool PassPSDEnergy(double energyPSD, double energyPSDPeripheral);
};
