#include "EventCuts.h"

bool EventCuts::PassVertexZ(double vZ) {// Check if the vertex Z position is within the specified range defined in Config.h
    if (vZ < Config::z_min || vZ > Config::z_max) return false; // Check if the vertex Z position is within the specified range
    return true;
}

bool EventCuts::PassPSDEnergy(double energyPSD, double energyPSDPeripheral) { // Check if the PSD energy is above the specified peripheral energy cut defined in Config.h
    (void)energyPSD;
    if (energyPSDPeripheral < Config::PSD_per_cut) return false;
    return true;
}

bool EventCuts::PassTracksRatio(int nTracks, int nTracksInFit) { // Check if the ratio of fitted tracks to total tracks
    double x = static_cast<double>(nTracks);
    double y = static_cast<double>(nTracksInFit);

    double D1 = (Config::P2_x - Config::P1_x) * (y - Config::P1_y) -
    (Config::P2_y - Config::P1_y) * (x - Config::P1_x); // Calculate the determinant to determine the position of the point (x, y) relative to the line segment P1-P2

    double D2 = (Config::P4_x - Config::P3_x) * (y - Config::P3_y) -
    (Config::P4_y - Config::P3_y) * (x - Config::P3_x);

    if (D1 < 0.0 && x < Config::P2_x) return false;
    if (D2 < 0.0 && x >= Config::P3_x) return false;

    return true;
} 