#pragma once

namespace Config {
    constexpr int BeamMomentum = 40;
    constexpr double PSD_per_cut = 300.0;

    constexpr double z_peak = -602.9;
    constexpr double z_min = z_peak - 2.0;
    constexpr double z_max = z_peak + 2.0;

    constexpr unsigned int VTPC_sum_min = 15;
    constexpr unsigned int ClustersAll_min = 30;
    constexpr unsigned int dEdx_min = 30;
    constexpr double Ratio_min = 0.5;
    constexpr double Ratio_max = 1.1;
    
    constexpr double ImpactParam_bx_max = 4.0; 
    constexpr double ImpactParam_by_max = 2.0; 

    constexpr double P1_x = 0.0;    constexpr double P1_y = 0.0;
    constexpr double P2_x = 350.0;  constexpr double P2_y = 52.5;
    constexpr double P3_x = 350.0;  constexpr double P3_y = 52.5;
    constexpr double P4_x = 1350.0; constexpr double P4_y = 350.0;
}