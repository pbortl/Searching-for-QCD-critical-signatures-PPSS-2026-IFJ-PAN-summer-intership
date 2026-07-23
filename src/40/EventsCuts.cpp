#include "EventCuts.h"

bool EventCuts::PassVertexZ(double vZ) {
    if (vZ < Config::z_min || vZ > Config::z_max) return false;
    return true;
}

bool EventCuts::PassPSDEnergy(double energyPSD, double energyPSDPeripheral) {
    (void)energyPSD;
    if (energyPSDPeripheral < Config::PSD_per_cut) return false;
    return true;
}

bool EventCuts::PassTracksRatio(int nTracks, int nTracksInFit) {
    double x = static_cast<double>(nTracks);
    double y = static_cast<double>(nTracksInFit);

    double D1 = (Config::P2_x - Config::P1_x) * (y - Config::P1_y) -
    (Config::P2_y - Config::P1_y) * (x - Config::P1_x);

    double D2 = (Config::P4_x - Config::P3_x) * (y - Config::P3_y) -
    (Config::P4_y - Config::P3_y) * (x - Config::P3_x);

    if (D1 < 0.0 && x < Config::P2_x) return false;
    if (D2 < 0.0 && x >= Config::P3_x) return false;

    return true;
} 