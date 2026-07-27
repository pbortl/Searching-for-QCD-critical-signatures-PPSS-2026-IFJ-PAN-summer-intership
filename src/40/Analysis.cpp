#include "Analysis.h"
#include "Config.h"
#include "HistogramManager.h"
#include "EventCuts.h"
#include "Plotter.h"
#include "BetheBloch/src/include/BetheBlochWrapper.h"

#include <iostream>
#include <fstream>
#include <string>
#include <chrono>
#include <ctime>
#include <iomanip>
#include <sstream>
#include <cmath>
#include <vector>
#include <utility>

#include "TROOT.h"
#include "TChain.h"
#include "TSystem.h"
#include "TError.h"
#include "TFile.h"

Analysis::Analysis(const std::string& fileListPath) : fileList(fileListPath) {}

void Analysis::LoadLibraries() {
    static bool initialized = false;
    if (!initialized) {
        std::cout << "Initializing libraries..." << std::endl;
        const char* libs[] = {
            "/home/p/Documents/PPSS_tools/BetheBloch/BetheBlochWrapper.so",
            "/home/p/Documents/PPSS_tools/Tools/Event.so",
            "/home/p/Documents/PPSS_tools/Tools/CutsMap.noDict.so",
            "/home/p/Documents/PPSS_tools/Tools/EventMixer.noDict.so",
            "/home/p/Documents/PPSS_tools/Tools/EventXeLa.so",
            "/home/p/Documents/PPSS_tools/Tools/EventXeLaMag.so"
        };
        for (const auto& lib : libs) {
            if (gSystem->Load(lib, "", kTRUE) < 0) {
                std::cerr << "ERROR: Failed to load: " << lib << std::endl;
            }
        }
        initialized = true;
    }
}

void Analysis::Run() {
    gROOT->SetBatch(kTRUE);
    gErrorIgnoreLevel = kError;
    LoadLibraries();

    std::string outputDir;
    std::cout << "Enter the path to the output directory (e.g., /home/user/results): ";
    std::getline(std::cin, outputDir);

    if (!outputDir.empty() && outputDir.back() == '/') {
        outputDir.pop_back();
    }

    std::cout << "--- Applying cuts for Xe+La " << Config::BeamMomentum << "A GeV/c ---" << std::endl;

    auto now = std::chrono::system_clock::now();
    std::time_t now_c = std::chrono::system_clock::to_time_t(now);
    std::tm* now_tm = std::localtime(&now_c);
    std::stringstream ss_time;

    ss_time << std::put_time(now_tm, "%H_%M_%S_%d_%m_%Y");

    std::string eventName = "EventXeLaMag_" + std::to_string(Config::BeamMomentum) + "A";
    std::string baseFileName = outputDir + "/" + ss_time.str() + "_" + eventName;
    
    std::string rapidityTxtPath = baseFileName + "_rapidity_details.txt";
    std::ofstream rapFile(rapidityTxtPath);
    int rap_print_counter = 0; 

    TChain chain("event_tree");
    std::ifstream file(fileList);
    if (!file.is_open()) {
        std::cerr << "ERROR: Cannot open the file with paths list!" << std::endl;
        return;
    }

    std::string filePath;
    while (std::getline(file, filePath)) {
        if (!filePath.empty()) chain.Add(filePath.c_str());
    }
    file.close();

    chain.SetBranchStatus("*", 1);
    EventClass::EventXeLaMag* pevent = nullptr;
    chain.SetBranchAddress("EventXeLaMag", &pevent);

    Long64_t nentries = chain.GetEntries();
    std::cout << "Analyzing " << nentries << " events..." << std::endl;

    HistogramManager hists;
    Dedx::BetheBlochWrapper bbWrapper;

    unsigned int nEvents = 0;
    unsigned int nEvents_PSD = 0;
    unsigned int nEvents_VertexZ = 0;
    unsigned int nEvents_tracksratio = 0;

    unsigned int nTracks = 0;
    unsigned int nTracks_VTPC12 = 0;
    unsigned int nTracks_ClustersAll = 0;
    unsigned int nTracks_dEdx = 0;
    unsigned int nTracks_pointRatio = 0;
    unsigned int nTracks_ImpactParameter = 0;
    unsigned int nTracks_protons = 0;

    unsigned int nEvents_before_centrality = 0;
    unsigned int nTracksFit_before_centrality = 0;
    unsigned int nEvents_after_centrality = 0;
    unsigned int nTracksFit_after_centrality = 0;

    double mp = 0.93827;
    double pL = Config::BeamMomentum; //40
    double E_L = std::sqrt(mp*mp + pL*pL);
    double E_CM = 0.5 * std::sqrt(2.0 * mp * (mp + E_L));
    double Y_beam_shift = std::acosh(E_CM / mp);
    //rapility in CM frame: Y_CM = Y_LAB - Y_beam_shift first part

    // F2(M) I
    const int M_MAX = 150;
    const double L_window = 3.0; // px, py (from -1.5 to 1.5)
    const double R_0 = L_window / std::sqrt(M_PI);
    const double m_f = 0.27;

    std::vector<double> sum_N2(M_MAX, 0.0);
    std::vector<double> sum_N2_sq(M_MAX, 0.0);
    double sum_mul = 0.0;
    double sum_mul_sq = 0.0;
    long long f2_valid_events = 0;
   

    if (rapFile.is_open()) {
        rapFile << "=BEAM KINEMATICS =\n";
        rapFile << "Rest mass (mp): " << mp << " GeV/c^2\n";
        rapFile << "Beam momentum (pL): " << pL << " A GeV/c\n";
        rapFile << "Beam energy LAB (E_L): " << E_L << " GeV\n";
        rapFile << "CM energy (E_CM): " << E_CM << " GeV\n";
        rapFile << "CM rapidity shift (Y_beam_shift): " << Y_beam_shift << "\n\n";
        rapFile << "= TRACK RAPIDITY (Y_LAB - Y_beam_shift = Y_CM) =\n";
        rapFile << "Printing values for the first 1000 tracks:\n\n";
        rapFile << std::setw(15) << "Y_LAB" << std::setw(15) << "Y_beam_shift" << std::setw(15) << "Y_CM" << "\n";
    }

    std::cout << "dynamic centrality limit" << std::endl;
    for (Long64_t ievents = 0; ievents < nentries; ievents++) {
        chain.GetEntry(ievents);

        if (ievents > 0 && ievents % 100000 == 0) {
            double percent = (static_cast<double>(ievents) / nentries) * 100.0;
            std::cout << "Phase 1 Processed: " << ievents << " / " << nentries << " (" << static_cast<int>(percent) << "%)" << std::endl;
        }

        if (!pevent) continue;

        if (pevent->run_number == 34948 || pevent->run_number == 34976 ||
            pevent->run_number == 35093 || pevent->run_number == 35142) {
            continue;
        }
        hists.h_PSD_T2->Fill(pevent->energyPSDSelectedModules);
    }

    //counting centrality limits and events for each centrality class
    std::vector<double> cent_limits(4, 0.0);
    std::vector<int> cent_events(4, 0);
    double fracs[4] = {0.05, 0.10, 0.15, 0.20};

    TH1* cumul = hists.h_PSD_T2->GetCumulative();
    TH1* cumul_raw = hists.h_PSD_T2->GetCumulative();

    if (cumul->GetMaximum() > 0) {
        cumul->Scale(0.423 / cumul->GetMaximum());
    }

    for (int i = 0; i < 4; ++i) {
        int bin = cumul->FindFirstBinAbove(fracs[i]);
        if (bin > 0) {
            cent_limits[i] = cumul->GetXaxis()->GetBinCenter(bin);
            cent_events[i] = cumul_raw->GetBinContent(bin);
        }
    }
    double dynamic_limit_20 = cent_limits[3];
    delete cumul;
    delete cumul_raw;

    std::cout << ">> Dynamic 20% centrality threshold calculated as: " << dynamic_limit_20 << " GeV\n" << std::endl;

    std::cout << "appling cuts" << std::endl;
    for (Long64_t ievents = 0; ievents < nentries; ievents++) {
        chain.GetEntry(ievents);

        if (ievents > 0 && ievents % 100000 == 0) {
            double percent = (static_cast<double>(ievents) / nentries) * 100.0;
            std::cout << "Processed: " << ievents << " / " << nentries << " (" << static_cast<int>(percent) << "%)" << std::endl;
        }

        if (!pevent) continue;

        if (pevent->run_number == 34948 || pevent->run_number == 34976 ||
            pevent->run_number == 35093 || pevent->run_number == 35142) {
            continue;
        }

        nEvents++;

        hists.h_Vz_all->Fill(pevent->VertexZ);
        hists.h2_TracksInFit_vs_PSD_all->Fill(pevent->energyPSDSelectedModules, pevent->nTracksFit);
        hists.h2_Tracks_vs_PSD_all->Fill(pevent->energyPSDSelectedModules, pevent->nTracksAll);

        if (EventCuts::PassVertexZ(pevent->VertexZ)) {
            nEvents_VertexZ++;
            hists.h_Vz_cut->Fill(pevent->VertexZ);
            hists.h2_TracksInFit_vs_PSD_after_Vz->Fill(pevent->energyPSDSelectedModules, pevent->nTracksFit);
            hists.h2_PSD_Peripheral_vs_Selected->Fill(pevent->energyPSDSelectedModules, pevent->energyPSDPeripheralModules);
            hists.h_PSDperipheral_all->Fill(pevent->energyPSDPeripheralModules);

            if (pevent->energyPSDPeripheralModules > Config::PSD_per_cut) {
                nEvents_PSD++;
                hists.h_PSDperipheral_cut->Fill(pevent->energyPSDPeripheralModules);
                hists.h2_PSD_Peripheral_vs_Selected_cut->Fill(pevent->energyPSDSelectedModules, pevent->energyPSDPeripheralModules);
                hists.h2_TracksInFit_vs_PSD_after_PSD->Fill(pevent->energyPSDSelectedModules, pevent->nTracksFit);
                hists.h2_TracksRatio_all->Fill(pevent->nTracksAll, pevent->nTracksFit);

                if (EventCuts::PassTracksRatio(pevent->nTracksAll, pevent->nTracksFit)) {
                    nEvents_tracksratio++;
                    hists.h2_TracksRatio_cut->Fill(pevent->nTracksAll, pevent->nTracksFit);
                    hists.h2_Tracks_vs_PSD_cut->Fill(pevent->energyPSDSelectedModules, pevent->nTracksAll);
                    hists.h2_TracksInFit_vs_PSD_cut->Fill(pevent->energyPSDSelectedModules, pevent->nTracksFit);

                    nEvents_before_centrality++;
                    nTracksFit_before_centrality += pevent->nTracksFit;

                    if (pevent->energyPSDSelectedModules < dynamic_limit_20) {

                        hists.h2_TracksInFit_vs_PSD_Central->Fill(pevent->energyPSDSelectedModules, pevent->nTracksFit);

                        nEvents_after_centrality++;
                        nTracksFit_after_centrality += pevent->nTracksFit;

                        
                        std::vector<std::pair<double, double>> event_protons;// Store accepted positive protons for F2 calculation
                        // Loop over tracks in the event

                        for (const auto& tracks : pevent->tracks) {
                            nTracks++;

                            unsigned int cl_dEdx = tracks.clustersdEdx;
                            unsigned int cl_VTPC_sum = tracks.clustersVTPC1 + tracks.clustersVTPC2;
                            unsigned int cl_Pot = tracks.clustersPotentialAll;
                            unsigned int cl_All = tracks.clustersAll;
                            double ratio = (cl_Pot > 0) ? static_cast<double>(cl_All) / static_cast<double>(cl_Pot) : 0.0;

                            hists.h_clusters_VTPC_sum_all->Fill(cl_VTPC_sum);
                            if (cl_VTPC_sum <= Config::VTPC_sum_min) continue;
                            nTracks_VTPC12++;
                            hists.h_clusters_VTPC_sum_cut->Fill(cl_VTPC_sum);

                            hists.h_clusters_All_all->Fill(cl_All);
                            if (cl_All <= Config::ClustersAll_min) continue;
                            nTracks_ClustersAll++;
                            hists.h_clusters_All_cut->Fill(cl_All);

                            hists.h_clusters_dEdx_all->Fill(cl_dEdx);
                            if (cl_dEdx <= Config::dEdx_min) continue;
                            nTracks_dEdx++;
                            hists.h_clusters_dEdx_cut->Fill(cl_dEdx);

                            hists.h_clusters_Ratio_all->Fill(ratio);
                            hists.h_clusters_PotAll_all->Fill(cl_Pot);
                            if (ratio < Config::Ratio_min || ratio > Config::Ratio_max) continue;
                            nTracks_pointRatio++;
                            hists.h_clusters_Ratio_cut->Fill(ratio);
                            hists.h_clusters_PotAll_cut->Fill(cl_Pot);

                            hists.h2_bx_by_all->Fill(tracks.bx, tracks.by);
                            if (std::abs(tracks.bx) > Config::ImpactParam_bx_max || std::abs(tracks.by) > Config::ImpactParam_by_max) continue;
                            nTracks_ImpactParameter++;
                            hists.h2_bx_by_cut->Fill(tracks.bx, tracks.by);

                            //vector momentum components all
                            double ptot = std::sqrt(tracks.px * tracks.px + tracks.py * tracks.py + tracks.pz * tracks.pz);
                            double log_ptot = std::log10(ptot);

                            double bb_proton_dedx = bbWrapper(3, ptot);
                            double bb_kaon_dedx = bbWrapper(2, ptot);
                            double delta_kp = bb_kaon_dedx - bb_proton_dedx;

                            //counting of rapidity in centre
                            double E_track = std::sqrt(mp*mp + ptot*ptot);
                            double Y_LAB_track = 0.5 * std::log((E_track + tracks.pz) / (E_track - tracks.pz)); //Frame of reference
                            double Y_CM_track = Y_LAB_track - Y_beam_shift;

                            if (rapFile.is_open() && rap_print_counter < 1000) {
                                rapFile << std::setw(15) << Y_LAB_track 
                                        << std::setw(15) << Y_beam_shift 
                                        << std::setw(15) << Y_CM_track << "\n";
                                rap_print_counter++;
                            }

                            //cuttings
                            bool pass_log_ptot = (log_ptot >= 0.55 && log_ptot <= 2.0);
                            bool pass_px = (tracks.px > -1.5 && tracks.px < 1.5);
                            bool pass_py = (tracks.py > -1.5 && tracks.py < 1.5);
                            bool pass_rapidity = (std::abs(Y_CM_track) <= 0.75);

                            //bloch 
                            double upper_limit = bb_proton_dedx + 0.15 * delta_kp;

                            if (tracks.dEdx > 0) {
                                hists.h_dedx_ptot_pos->Fill(log_ptot, tracks.dEdx);
                                bool pass_dedx = (tracks.dEdx <= upper_limit);

                                if (pass_log_ptot && pass_px && pass_py && pass_dedx && pass_rapidity) {
                                    hists.h_dedx_ptot_protons->Fill(log_ptot, tracks.dEdx);
                                    nTracks_protons++;
                                    //moemntum components for positive tracks
                                    hists.h2_px_py_pos->Fill(tracks.px, tracks.py);
                                    hists.h_Y_CM_tracks->Fill(Y_CM_track);

                                    // Store accepted positive protons for F2 calculation
                                    event_protons.push_back({tracks.px, tracks.py});
                                }

                            } else if (tracks.dEdx < 0) {
                                double abs_dedx = std::abs(tracks.dEdx);
                                hists.h_dedx_ptot_neg->Fill(log_ptot, abs_dedx);
                                bool pass_dedx = (abs_dedx <= upper_limit);

                                if (pass_log_ptot && pass_px && pass_py && pass_dedx && pass_rapidity) {
                                    hists.h_dedx_ptot_neg_protons->Fill(log_ptot, abs_dedx);
                                    //moemntum components for negative tracks
                                    hists.h2_px_py_neg->Fill(tracks.px, tracks.py);
                                    hists.h_Y_CM_tracks->Fill(Y_CM_track);
                                }
                            }
                        }

                        // F2(M) II 
                        double current_mul = event_protons.size();
                        sum_mul += current_mul;
                        sum_mul_sq += current_mul * current_mul;
                        f2_valid_events++;

                        std::vector<double> event_N2(M_MAX, 0.0);

                        for (size_t i = 0; i < event_protons.size(); ++i) {
                            for (size_t j = i + 1; j < event_protons.size(); ++j) {
                                double dx = event_protons[i].first - event_protons[j].first;
                                double dy = event_protons[i].second - event_protons[j].second;
                                double dist = std::sqrt(dx*dx + dy*dy);

                                if (dist > 0.0) {
                                    int idx = std::floor(R_0 / dist + m_f);
                                    if (idx > M_MAX) idx = M_MAX;
                                    
                                    for (int m = 0; m < idx; ++m) {
                                        event_N2[m] += 1.0;
                                    }
                                }
                            }
                        }

                        for (int m = 0; m < M_MAX; ++m) {
                            sum_N2[m] += event_N2[m];
                            sum_N2_sq[m] += event_N2[m] * event_N2[m];
                        }
                        
                    }
                }
            }
        }
    }

    hists.hist_events->SetBinContent(1, nEvents);
    hists.hist_events->SetBinContent(2, nEvents_VertexZ);
    hists.hist_events->SetBinContent(3, nEvents_PSD);
    hists.hist_events->SetBinContent(4, nEvents_tracksratio);
    hists.hist_events->SetBinContent(5, nEvents_after_centrality);

    hists.hist_tracks->SetBinContent(1, nTracks);
    hists.hist_tracks->SetBinContent(2, nTracks_VTPC12);
    hists.hist_tracks->SetBinContent(3, nTracks_ClustersAll);
    hists.hist_tracks->SetBinContent(4, nTracks_dEdx);
    hists.hist_tracks->SetBinContent(5, nTracks_pointRatio);
    hists.hist_tracks->SetBinContent(6, nTracks_ImpactParameter);
    hists.hist_tracks->SetBinContent(7, nTracks_protons);
    hists.hist_tracks->SetBinContent(8, nTracksFit_before_centrality);

    if (rapFile.is_open()) {
        rapFile.close();
    }

    
    std::string f2TxtPath = baseFileName + "_F2_results.txt";
    std::ofstream f2File(f2TxtPath);
    // F2(M) III 
    if (f2File.is_open() && f2_valid_events > 0) {
        double avg_mul = sum_mul / f2_valid_events;
        double avg_mul_sq = sum_mul_sq / f2_valid_events;// Calculate variance and error for multiplicity
        
        double var_mul = avg_mul_sq - (avg_mul * avg_mul);// Variance of multiplicity
        if (var_mul < 0) var_mul = 0.0; // Ensure non-negative variance
        double err_mul = std::sqrt(var_mul) / std::sqrt(f2_valid_events);

        f2File << "M\tF2_M\terr_F2_M\tavg_N2\terr_N2\tavg_mul\terr_mul\n";// Header for the output file

        for (int m = 0; m < M_MAX; ++m) {
            int M_val = m + 1; 
            
            double avg_N2 = sum_N2[m] / f2_valid_events;
            double avg_N2_sq = sum_N2_sq[m] / f2_valid_events;
            
            double var_N2 = avg_N2_sq - (avg_N2 * avg_N2);
            if (var_N2 < 0) var_N2 = 0.0;
            double err_N2 = std::sqrt(var_N2) / std::sqrt(f2_valid_events);

            double F2_M = 0.0;
            double err_F2_M = 0.0;

            if (avg_mul > 0 && avg_N2 > 0) {
                F2_M = (2.0 * M_val * M_val * avg_N2) / (avg_mul * avg_mul);
                
                double rel_err_N2 = err_N2 / avg_N2;
                double rel_err_mul = err_mul / avg_mul;
                
                err_F2_M = F2_M * std::sqrt(rel_err_N2 * rel_err_N2 + 4.0 * rel_err_mul * rel_err_mul);
            }

            f2File << M_val << "\t" 
                   << F2_M << "\t" 
                   << err_F2_M << "\t" 
                   << avg_N2 << "\t" 
                   << err_N2 << "\t" 
                   << avg_mul << "\t" 
                   << err_mul << "\n";
        }
        f2File.close();
        std::cout << "=> F2(M) calculation finished! Results saved to: " << f2TxtPath << std::endl;
    }
    

    std::cout << "Analysis finished successfully!" << std::endl;
    std::cout << "=> Rapidity parameters saved to file: " << rapidityTxtPath << std::endl;

    std::cout << "\n=== Centrality Stats ===" << std::endl;
    std::cout << "Events surviving 0-20% cut: " << nEvents_after_centrality << std::endl;
    if (nEvents_before_centrality > 0)
        std::cout << "Avg tracks/event (before cut): " << (double)nTracksFit_before_centrality / nEvents_before_centrality << std::endl;
    if (nEvents_after_centrality > 0)
        std::cout << "Avg tracks/event (after cut):  " << (double)nTracksFit_after_centrality / nEvents_after_centrality << std::endl;
    std::cout << "========================\n" << std::endl;

    Plotter::DrawAndSaveAll(hists, baseFileName, cent_limits, cent_events);
}