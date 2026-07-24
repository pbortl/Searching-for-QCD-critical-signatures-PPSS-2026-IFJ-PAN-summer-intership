#include "HistogramManager.h"
#include "Config.h"

HistogramManager::HistogramManager() {
    h_Vz_all = new TH1F("h_Vz_all", "VertexZ (All); Vertex Z [cm]; Entries", 500, -613, -593);
    h_Vz_cut = new TH1F("h_Vz_cut", "VertexZ (Cut); Vertex Z [cm]; Entries", 500, -613, -593);

    h2_Tracks_vs_PSD_all = new TH2F("h2_Tracks_vs_PSD_all", "Tracks vs PSD (All); PSD Energy [GeV]; Total Tracks", 350, 0, 7000, 150, 0, 1500);
    h2_Tracks_vs_PSD_cut = new TH2F("h2_Tracks_vs_PSD_cut", "Tracks vs PSD (Cut); PSD Energy [GeV]; Total Tracks", 350, 0, 7000, 150, 0, 1500);

    h2_PSD_Peripheral_vs_Selected = new TH2F("h2_PSD_Peripheral_vs_Selected", "PSD energy peripheral vs selected (All); PSD Selected Energy [GeV]; PSD Peripheral Energy [GeV]", 350, 0, 7000, 100, 0, 2000);
    h2_PSD_Peripheral_vs_Selected_cut = new TH2F("h2_PSD_Peripheral_vs_Selected_cut", "PSD energy after cuts (Cut); PSD Selected Energy [GeV]; PSD Peripheral Energy [GeV]", 350, 0, 7000, 100, 0, 2000);

    h_PSDperipheral_all = new TH1F("h_PSDperipheral_all", "PSD peripheral energy (All); Peripheral PSD Energy [GeV]; Entries", 500, 0, 1500);
    h_PSDperipheral_cut = new TH1F("h_PSDperipheral_cut", "PSD peripheral energy (Cut); Peripheral PSD Energy [GeV]; Entries", 400, Config::PSD_per_cut, 1500);

    h2_TracksRatio_all = new TH2F("h2_TracksRatio_all", "All tracks vs tracks in fit (All); nTracksAll; nTracksInFit", 150, 0, 1500, 350, 0, 350);
    h2_TracksRatio_cut = new TH2F("h2_TracksRatio_cut", "All tracks vs tracks in fit (Cut); nTracksAll; nTracksInFit", 150, 0, 1500, 350, 0, 350);

    h2_TracksInFit_vs_PSD_all = new TH2F("h2_TracksInFit_vs_PSD_all", "Selected energy vs tracks in fit (All); PSD Energy (selected) [GeV]; nTracks in fit", 350, 0, 7000, 350, 0, 350);
    h2_TracksInFit_vs_PSD_after_PSD = new TH2F("h2_TracksInFit_vs_PSD_after_PSD", "Selected energy vs tracks in fit (After PSD Cut); PSD Energy (selected) [GeV]; nTracks in fit", 350, 0, 7000, 350, 0, 350);
    h2_TracksInFit_vs_PSD_after_Vz = new TH2F("h2_TracksInFit_vs_PSD_after_Vz", "Selected energy vs tracks in fit (After Vz Cut); PSD Energy (selected) [GeV]; nTracks in fit", 350, 0, 7000, 350, 0, 350);
    h2_TracksInFit_vs_PSD_cut = new TH2F("h2_TracksInFit_vs_PSD_cut", "Selected energy vs tracks in fit (Cut); PSD Energy (selected) [GeV]; nTracks in fit", 350, 0, 7000, 350, 0, 350);
    
    h2_TracksInFit_vs_PSD_Central = new TH2F("h2_TracksInFit_vs_PSD_Central", "Selected energy vs tracks in fit (0-20% Central); PSD Energy (selected) [GeV]; nTracks in fit", 350, 0, 7000, 350, 0, 350);

    h_PSD_T2 = new TH1F("h_PSD_T2", "PSD energy T2; Energy [GeV]; Events", 500, 0, 5000);

    h_clusters_dEdx_all = new TH1F("h_clusters_dEdx_all", "Clusters dE/dx (All); #clusters dE/dx; # tracks", 250, 0, 250);
    h_clusters_dEdx_cut = new TH1F("h_clusters_dEdx_cut", "Clusters dE/dx (Cut); #clusters dE/dx; # tracks", 250, 0, 250);

    h_clusters_VTPC_sum_all = new TH1F("h_clusters_VTPC_sum_all", "Clusters VTPC1+VTPC2 (All); VTPC sum; # tracks", 190, 0, 190);
    h_clusters_VTPC_sum_cut = new TH1F("h_clusters_VTPC_sum_cut", "Clusters VTPC1+VTPC2 (Cut); VTPC sum; # tracks", 190, 0, 190);

    h_clusters_PotAll_all = new TH1F("h_clusters_PotAll_all", "Clusters Potential All (All); Pot. All; # tracks", 300, 0, 300);
    h_clusters_PotAll_cut = new TH1F("h_clusters_PotAll_cut", "Clusters Potential All (Cut); Pot. All; # tracks", 300, 0, 300);

    h_clusters_All_all = new TH1F("h_clusters_All_all", "Clusters All (All); All clusters; # tracks", 300, 0, 300);
    h_clusters_All_cut = new TH1F("h_clusters_All_cut", "Clusters All (Cut); All clusters; # tracks", 300, 0, 300);

    h_clusters_Ratio_all = new TH1F("h_clusters_Ratio_all", "Ratio Points/Pot.points (All); Ratio; # tracks", 200, 0, 2.0);
    h_clusters_Ratio_cut = new TH1F("h_clusters_Ratio_cut", "Ratio Points/Pot.points (Cut); Ratio; # tracks", 200, 0, 2.0);

    h2_bx_by_all = new TH2F("h2_bx_by_all", "b_{x} vs b_{y} (All); b_{x} [cm]; b_{y} [cm]", 200, -10.0, 10.0, 200, -10.0, 10.0);
    h2_bx_by_cut = new TH2F("h2_bx_by_cut", "b_{x} vs b_{y} (Cut); b_{x} [cm]; b_{y} [cm]", 200, -10.0, 10.0, 200, -10.0, 10.0);
    
    h_dedx_ptot_pos = new TH2F("h_dedx_ptot_pos", "dE/dx vs log(p_tot) positive tracks; Log_{10}(p_{tot} / [GeV/c]); dE/dx (MIP)", 300, -0.6, 2.2, 300, 0.4, 2.0);
    h_dedx_ptot_neg = new TH2F("h_dedx_ptot_neg", "dE/dx vs log(p_tot) negative tracks; Log_{10}(p_{tot} / [GeV/c]); dE/dx (MIP)", 300, -0.6, 2.2, 300, 0.4, 2.0);
    
    h_dedx_ptot_protons = new TH2F("h_dedx_ptot_protons", "dE/dx vs log(p_tot) selected protons; Log_{10}(p_{tot} / [GeV/c]); dE/dx (MIP)", 300, -0.6, 2.2, 300, 0.4, 2.0);
    h_dedx_ptot_neg_protons = new TH2F("h_dedx_ptot_neg_protons", "dE/dx vs log(p_tot) selected anti-protons; Log_{10}(p_{tot} / [GeV/c]); dE/dx (MIP)", 300, -0.6, 2.2, 300, 0.4, 2.0);
    
    h2_px_py_pos = new TH2F("h2_px_py_pos", "Positive charge p_{y} vs p_{x}; p_{x} [GeV/c]; p_{y} [GeV/c]", 300, -1.5, 1.5, 300, -1.5, 1.5);
    h2_px_py_neg = new TH2F("h2_px_py_neg", "Negative charge p_{y} vs p_{x}; p_{x} [GeV/c]; p_{y} [GeV/c]", 300, -1.5, 1.5, 300, -1.5, 1.5);
    h_Y_CM_tracks = new TH1F("h_Y_CM_tracks", "Selected Poxitive and Negative tracks y^{CM}_{track}; counts", 200, -3.0, 3.0);

    TH2F* h2_array[] = {h2_Tracks_vs_PSD_all, h2_Tracks_vs_PSD_cut, 
                        h2_PSD_Peripheral_vs_Selected, h2_PSD_Peripheral_vs_Selected_cut, 
                        h2_TracksRatio_all, h2_TracksRatio_cut, 
                        h2_TracksInFit_vs_PSD_all, h2_TracksInFit_vs_PSD_after_PSD,
                        h2_TracksInFit_vs_PSD_after_Vz, h2_TracksInFit_vs_PSD_cut, 
                        h2_TracksInFit_vs_PSD_Central,
                        h2_bx_by_all, h2_bx_by_cut};
    
    for (auto* h2 : h2_array) {
        h2->SetOption("COLZ");
        h2->SetContour(255);
    }

    hist_events = new TH1F("hist_events", "Event Cuts Statistics; Cut Step; Number of Events", 5, 0.5, 5.5);
    hist_events->GetXaxis()->SetBinLabel(1, "before cuts");
    hist_events->GetXaxis()->SetBinLabel(2, "Vertex Z");
    hist_events->GetXaxis()->SetBinLabel(3, "PSD Energy");
    hist_events->GetXaxis()->SetBinLabel(4, "Tracks Ratio");
    hist_events->GetXaxis()->SetBinLabel(5, "0-20% Centrality");
    hist_events->SetStats(0);

    hist_tracks = new TH1F("hist_tracks", "Track Cuts Statistics; Cut Step; Number of Tracks", 7, 0.5, 7.5);
    hist_tracks->GetXaxis()->SetBinLabel(1, "before track cuts");
    hist_tracks->GetXaxis()->SetBinLabel(2, "VTPC1+2 clusters > 15");
    hist_tracks->GetXaxis()->SetBinLabel(3, "All clusters > 30");
    hist_tracks->GetXaxis()->SetBinLabel(4, "dE/dx clusters > 30");
    hist_tracks->GetXaxis()->SetBinLabel(5, "0.5 < points/pot. points < 1.1");
    hist_tracks->GetXaxis()->SetBinLabel(6, "Impact Parameter");
    hist_tracks->GetXaxis()->SetBinLabel(7, "Protons Selected");
    hist_tracks->SetStats(0);
}

HistogramManager::~HistogramManager() {
    delete h_Vz_all; delete h_Vz_cut;
    delete h2_Tracks_vs_PSD_all; delete h2_Tracks_vs_PSD_cut;
    delete h2_PSD_Peripheral_vs_Selected; delete h2_PSD_Peripheral_vs_Selected_cut;
    delete h_PSDperipheral_all; delete h_PSDperipheral_cut;
    delete h2_TracksRatio_all; delete h2_TracksRatio_cut;
    delete h2_TracksInFit_vs_PSD_all; 
    delete h2_TracksInFit_vs_PSD_after_PSD;
    delete h2_TracksInFit_vs_PSD_after_Vz;
    delete h2_TracksInFit_vs_PSD_cut;
    delete h2_TracksInFit_vs_PSD_Central;
    delete h_PSD_T2;

    delete h_clusters_dEdx_all; delete h_clusters_dEdx_cut;
    delete h_clusters_VTPC_sum_all; delete h_clusters_VTPC_sum_cut;
    delete h_clusters_PotAll_all; delete h_clusters_PotAll_cut;
    delete h_clusters_All_all; delete h_clusters_All_cut;
    delete h_clusters_Ratio_all; delete h_clusters_Ratio_cut;
    delete h2_bx_by_all; delete h2_bx_by_cut;
    delete h_dedx_ptot_pos;
    delete h_dedx_ptot_neg;
    delete h_dedx_ptot_protons; 
    delete h_dedx_ptot_neg_protons;
    delete h2_px_py_pos;
    delete h2_px_py_neg;
    delete h_Y_CM_tracks;
    delete hist_events;
    delete hist_tracks;
}