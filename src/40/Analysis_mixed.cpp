#include "Analysis_mixed.h"
#include "Config.h"
#include "HistogramManager.h"
#include "EventCuts.h"
#include "Plotter.h" // Required for PDF/ROOT saving at the end
#include "BetheBloch/src/include/BetheBlochWrapper.h"

// Mixer libraries
#include "/home/p/Documents/PPSS_tools/Tools/src/include/EventMixerT.h"
#include "TRandom3.h"

#include <iostream>
#include <fstream>
#include <string>
#include <chrono>
#include <ctime>
#include <iomanip>
#include <sstream>
#include <cmath>
#include <vector>

#include "TROOT.h"
#include "TChain.h"
#include "TSystem.h"
#include "TError.h"
#include "TFile.h"

Analysis_mixed::Analysis_mixed(const std::string& fileListPath) : fileList(fileListPath) {} 

void Analysis_mixed::LoadLibraries() { 
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

void Analysis_mixed::Run() {
    gROOT->SetBatch(kTRUE);
    gErrorIgnoreLevel = kError;
    LoadLibraries();

    std::string outputDir;
    std::cout << "Enter the path to the output directory (e.g., /home/user/results): ";
    std::getline(std::cin, outputDir);
    if (!outputDir.empty() && outputDir.back() == '/') outputDir.pop_back();

    auto now = std::chrono::system_clock::now();
    std::time_t now_c = std::chrono::system_clock::to_time_t(now);
    std::tm* now_tm = std::localtime(&now_c);
    std::stringstream ss_time;
    ss_time << std::put_time(now_tm, "%H_%M_%S_%d_%m_%Y");

    std::string baseFileName = outputDir + "/" + ss_time.str() + "_EventXeLaMag_MIXED_" + std::to_string(Config::BeamMomentum) + "A";

    // Open a text file to export mixed proton px and py for inspection if needed
    std::string pxpyTxtPath = baseFileName + "_protons_pxpy.txt";
    std::ofstream pxpyFile(pxpyTxtPath);
    if (pxpyFile.is_open()) {
        pxpyFile << "px\tpy\n";
    }

    TChain chain("event_tree");
    std::ifstream file(fileList);
    std::string filePath;
    while (std::getline(file, filePath)) {
        if (!filePath.empty()) chain.Add(filePath.c_str());
    }
    file.close();

    chain.SetBranchStatus("*", 1);
    EventClass::EventXeLaMag* pevent = nullptr;
    chain.SetBranchAddress("EventXeLaMag", &pevent);

    Long64_t nentries = chain.GetEntries();
    
    HistogramManager hists;
    double mp = 0.93827;
    double pL = Config::BeamMomentum;
    double E_L = std::sqrt(mp*mp + pL*pL);
    double E_CM = 0.5 * std::sqrt(2.0 * mp * (mp + E_L));
    
    /* 
     * --- KINEMATICS: LAB to CM Lorentz Transformation ---
     * 1. Mandelstam s: (E_L + m_N)^2 - p_L^2 = 4(E_CM)^2
     * 2. CM Energy:    \sqrt{S_{NN}} = 2(E_CM) = \sqrt{2m_N(m_N + E_L)}
     * 3. Beam Rapidity: y_{beam}^{CM} = acosh(E_CM / m_N)
     * -> Track CM rapidity is then: y_CM = y_LAB - y_beam_CM
     */
    double Y_beam_shift = std::acosh(E_CM / mp);

    // Calculate centrality limit on CLEAN events (fills h_PSD_T2 for mixed)
    std::cout << "Step 1: Calculating dynamic centrality limit on CLEAN events..." << std::endl;
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

        // Apply event cuts BEFORE calculating centrality limit
        if (!EventCuts::PassVertexZ(pevent->VertexZ)) continue;
        if (pevent->energyPSDPeripheralModules <= Config::PSD_per_cut) continue;
        if (!EventCuts::PassTracksRatio(pevent->nTracksAll, pevent->nTracksFit)) continue;

        hists.h_PSD_T2->Fill(pevent->energyPSDSelectedModules);
    }

    std::vector<double> cent_limits(4, 0.0);
    std::vector<int> cent_events(4, 0); // Added for plot consistency
    double fracs[4] = {0.05, 0.10,0.15,0.20};
    
    TH1* cumul = hists.h_PSD_T2->GetCumulative();
    TH1* cumul_raw = hists.h_PSD_T2->GetCumulative();

    if (cumul->GetMaximum() > 0) cumul->Scale(0.423 / cumul->GetMaximum());
    for (int i = 0; i < 4; ++i) {
        int bin = cumul->FindFirstBinAbove(fracs[i]);
        if (bin > 0) {
            cent_limits[i] = cumul->GetXaxis()->GetBinCenter(bin);
            cent_events[i] = cumul_raw->GetBinContent(bin);
        }
    }
    double dynamic_limit_20 = cent_limits[1];
    delete cumul;
    delete cumul_raw;
    std::cout << ">> Dynamic 20% centrality threshold calculated as: " << dynamic_limit_20 << " GeV\n" << std::endl;

    // Mixer settings create
    std::cout << "Step 2: Starting Event Mixer..." << std::endl;
    std::string mixerCutsPath = outputDir + "/mixer_cuts.txt";
    std::ofstream mcf(mixerCutsPath);
    mcf << "bufferSize 20000\nmixTimes 10\n";
    mcf.close();

    TRandom3 rng(12345);

    // Event cuts lambda
    auto event_cuts = [dynamic_limit_20](const EventClass::EventXeLaMag& ev) {
        if (ev.run_number == 34948 || ev.run_number == 34976 ||
            ev.run_number == 35093 || ev.run_number == 35142) return false;
        if (!EventCuts::PassVertexZ(ev.VertexZ)) return false;
        if (ev.energyPSDPeripheralModules <= Config::PSD_per_cut) return false;
        if (!EventCuts::PassTracksRatio(ev.nTracksAll, ev.nTracksFit)) return false;
        if (ev.energyPSDSelectedModules >= dynamic_limit_20) return false;
        return true;
    };

    // Track cuts lambda (protons)
    Dedx::BetheBlochWrapper bbWrap;
    auto track_cuts = [&bbWrap, Y_beam_shift, mp](const EventClass::Track& t) mutable {
        unsigned int cl_VTPC_sum = t.clustersVTPC1 + t.clustersVTPC2;
        unsigned int cl_Pot = t.clustersPotentialAll;
        unsigned int cl_All = t.clustersAll;
        double ratio = (cl_Pot > 0) ? static_cast<double>(cl_All) / static_cast<double>(cl_Pot) : 0.0;

        if (cl_VTPC_sum <= Config::VTPC_sum_min) return false;
        if (cl_All <= Config::ClustersAll_min) return false;
        if (t.clustersdEdx <= Config::dEdx_min) return false;
        if (ratio < Config::Ratio_min || ratio > Config::Ratio_max) return false;
        if (std::abs(t.bx) > Config::ImpactParam_bx_max || std::abs(t.by) > Config::ImpactParam_by_max) return false;

        double ptot = std::sqrt(t.px * t.px + t.py * t.py + t.pz * t.pz);
        double log_ptot = std::log10(ptot);
        if (log_ptot < 0.55 || log_ptot > 2.0) return false;
        if (t.px <= -1.5 || t.px >= 1.5) return false;
        if (t.py <= -1.5 || t.py >= 1.5) return false;

        double E_track = std::sqrt(mp*mp + ptot*ptot);
        double Y_LAB_track = 0.5 * std::log((E_track + t.pz) / (E_track - t.pz));
        double Y_CM_track = Y_LAB_track - Y_beam_shift;
        if (std::abs(Y_CM_track) > 0.75) return false;

        if (t.dEdx <= 0) return false; // Only positive tracks
        
        /*
         * --- PARTICLE ID: Proton Selection via Bethe-Bloch (BB) ---
         * Particles lose energy (dE/dx) in the TPC. We select protons by placing
         * a cut slightly above the theoretical proton BB curve:
         * Cut = BB_proton + 0.15 * (BB_kaon - BB_proton)
         * This empirical 0.15 margin safely separates protons from kaons.
         */
        double bb_proton_dedx = bbWrap(3, ptot);
        double bb_kaon_dedx = bbWrap(2, ptot);
        double delta_kp = bb_kaon_dedx - bb_proton_dedx;
        if (t.dEdx > bb_proton_dedx + 0.15 * delta_kp) return false;

        return true; // Passed all cuts
    };

    EventMixerT<EventClass::EventXeLaMag> mixer(chain, rng, mixerCutsPath, event_cuts, track_cuts);
    auto eventStream = mixer.randomStream();

    /*
     * --- F2(M) CALCULATION: Correlation Integral (Moving Circles) ---
     * Instead of a rigid M x M grid, we use distance-based pair counting:
     * 1. Calculate distance between track pairs in p_x, p_y space.
     * 2. N_2(M) = number of pairs within radius R_M.
     * 3. F_2(M) relates to the correlation integral C(R):
     *    F_2(M) = (2 * M^2 * <N_2(M)>) / <mul>^2
     */
    // F2(M) INITIALIZATION
    const int M_MAX = 150;
    const double R_0 = 3.0 / std::sqrt(M_PI);
    const double m_f = 0.27;

    std::vector<double> sum_N2(M_MAX, 0.0);
    std::vector<double> sum_N2_sq(M_MAX, 0.0);
    double sum_mul = 0.0;
    double sum_mul_sq = 0.0;
    long long f2_valid_events = 0;

    int mixed_event_counter = 0;

    // Mixed events loop
    for (const auto& mixed_event : eventStream) {
        mixed_event_counter++;
        if (mixed_event_counter % 50000 == 0) {
            std::cout << "Mixed events processed: " << mixed_event_counter << std::endl;
        }

        double current_mul = mixed_event.tracks.size();
        sum_mul += current_mul;
        sum_mul_sq += current_mul * current_mul;
        f2_valid_events++;

        std::vector<double> event_N2(M_MAX, 0.0);

        for (size_t i = 0; i < mixed_event.tracks.size(); ++i) {
            const auto& t1 = mixed_event.tracks[i];
            
            // Fill histograms for mixed events
            double ptot = std::sqrt(t1.px * t1.px + t1.py * t1.py + t1.pz * t1.pz);
            double log_ptot = std::log10(ptot);
            double E_track = std::sqrt(mp*mp + ptot*ptot);
            double Y_LAB_track = 0.5 * std::log((E_track + t1.pz) / (E_track - t1.pz));
            double Y_CM_track = Y_LAB_track - Y_beam_shift;

            // Mixed events contain only final selected protons
            // Replaced the old histogram pointers with the new ones ending in "_rap"
            hists.h_dedx_ptot_protons_rap->Fill(log_ptot, t1.dEdx);
            hists.h2_px_py_protons_rap->Fill(t1.px, t1.py);
            hists.h_Y_CM_protons_rap->Fill(Y_CM_track);
            
            unsigned int cl_VTPC_sum = t1.clustersVTPC1 + t1.clustersVTPC2;
            unsigned int cl_Pot = t1.clustersPotentialAll;
            unsigned int cl_All = t1.clustersAll;
            double ratio = (cl_Pot > 0) ? static_cast<double>(cl_All) / static_cast<double>(cl_Pot) : 0.0;

            hists.h_clusters_VTPC_sum_cut->Fill(cl_VTPC_sum);
            hists.h_clusters_All_cut->Fill(cl_All);
            hists.h_clusters_dEdx_cut->Fill(t1.clustersdEdx);
            hists.h_clusters_Ratio_cut->Fill(ratio);
            hists.h_clusters_PotAll_cut->Fill(cl_Pot);
            hists.h2_bx_by_cut->Fill(t1.bx, t1.by);

            // Optional export of mixed proton px and py to text file for inspection
            if (pxpyFile.is_open()) {
                pxpyFile << t1.px << "\t" << t1.py << "\n";
            }

            for (size_t j = i + 1; j < mixed_event.tracks.size(); ++j) {
                const auto& t2 = mixed_event.tracks[j];
                double dx = t1.px - t2.px;
                double dy = t1.py - t2.py;
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

    if (pxpyFile.is_open()) {
        pxpyFile.close();
    }

    // Save mixed results to TXT
    std::string f2TxtPath = baseFileName + "_F2_results.txt"; // Name matched for Plotter consistency
    std::ofstream f2File(f2TxtPath);

    /*
     * --- ERROR PROPAGATION: Statistical Uncertainty for F2(M) ---
     * 1. Variance: Var(x) = <x^2> - <x>^2  (for both multiplicity and pairs)
     * 2. Standard Error: err_x = sqrt(Var(x) / N_events)
     * 3. F2 Error: Propagated from relative errors of pairs (N_2) and multiplicity (mul):
     *    dF2 / F2 = sqrt( (dN_2 / N_2)^2 + 4 * (dmul / mul)^2 )
     */
    if (f2File.is_open() && f2_valid_events > 0) {
        double avg_mul = sum_mul / f2_valid_events;
        double avg_mul_sq = sum_mul_sq / f2_valid_events;
        double var_mul = std::max(0.0, avg_mul_sq - (avg_mul * avg_mul));
        double err_mul = std::sqrt(var_mul) / std::sqrt(f2_valid_events);

        f2File << "M\tF2_M\terr_F2_M\tavg_N2\terr_N2\tavg_mul\terr_mul\n";

        for (int m = 0; m < M_MAX; ++m) {
            int M_val = m + 1; 
            double avg_N2 = sum_N2[m] / f2_valid_events;
            double avg_N2_sq = sum_N2_sq[m] / f2_valid_events;
            double var_N2 = std::max(0.0, avg_N2_sq - (avg_N2 * avg_N2));
            double err_N2 = std::sqrt(var_N2) / std::sqrt(f2_valid_events);

            double F2_M = 0.0, err_F2_M = 0.0;

            if (avg_mul > 0 && avg_N2 > 0) {
                F2_M = (2.0 * M_val * M_val * avg_N2) / (avg_mul * avg_mul);
                double rel_err_N2 = err_N2 / avg_N2;
                double rel_err_mul = err_mul / avg_mul;
                err_F2_M = F2_M * std::sqrt(rel_err_N2 * rel_err_N2 + 4.0 * rel_err_mul * rel_err_mul);
            }

            f2File << M_val << "\t" << F2_M << "\t" << err_F2_M << "\t" 
                   << avg_N2 << "\t" << err_N2 << "\t" << avg_mul << "\t" << err_mul << "\n";
        }
        f2File.close();
        std::cout << "=> F2(M) MIXED calculation finished! Results saved to: " << f2TxtPath << std::endl;
        std::cout << "=> Mixed proton px and py coordinates exported to file: " << pxpyTxtPath << std::endl;
    }

    // Generate histograms for mixed events
    std::cout << "Generating PDF plots for Mixed Events..." << std::endl;
    Plotter::DrawAndSaveAll(hists, baseFileName, cent_limits, cent_events);
}