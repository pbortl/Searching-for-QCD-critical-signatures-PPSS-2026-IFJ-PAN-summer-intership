#pragma once
#include "TH1F.h"
#include "TH2F.h"
#include "EventXeLaMag.h"

class HistogramManager {
public:
    HistogramManager();
    ~HistogramManager();

    TH1F* h_Vz_all;
    TH1F* h_Vz_cut;
    TH2F* h2_Tracks_vs_PSD_all;
    TH2F* h2_Tracks_vs_PSD_cut;
    TH2F* h2_PSD_Peripheral_vs_Selected;
    TH2F* h2_PSD_Peripheral_vs_Selected_cut;
    TH1F* h_PSDperipheral_all;
    TH1F* h_PSDperipheral_cut;
    TH2F* h2_TracksRatio_all;
    TH2F* h2_TracksRatio_cut;
    
    TH2F* h2_TracksInFit_vs_PSD_all;
    TH2F* h2_TracksInFit_vs_PSD_after_PSD;
    TH2F* h2_TracksInFit_vs_PSD_after_Vz;
    TH2F* h2_TracksInFit_vs_PSD_cut;
    TH2F* h_dedx_ptot_pos;
    TH2F* h_dedx_ptot_neg;
    TH1F* h_PSD_T2;
    TH2F* h2_TracksInFit_vs_PSD_Central;

    TH1F* h_clusters_dEdx_all;
    TH1F* h_clusters_dEdx_cut;
    TH1F* h_clusters_VTPC_sum_all;
    TH1F* h_clusters_VTPC_sum_cut;
    TH1F* h_clusters_PotAll_all;
    TH1F* h_clusters_PotAll_cut;
    TH1F* h_clusters_All_all;
    TH1F* h_clusters_All_cut;
    TH1F* h_clusters_Ratio_all;
    TH1F* h_clusters_Ratio_cut;
    TH2F* h2_bx_by_all;
    TH2F* h2_bx_by_cut;
    
    TH2F* h_dedx_ptot_protons;
    TH2F* h_dedx_ptot_neg_protons;
    TH2F* h2_px_py_pos;
    TH2F* h2_px_py_neg;
    TH1F* h_Y_CM_tracks;
    
    TH1F* hist_events;
    TH1F* hist_tracks;
};