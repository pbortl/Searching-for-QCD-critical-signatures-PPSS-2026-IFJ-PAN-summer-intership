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
    // Ensure the output directory exists
    std::cout << "--- Applying cuts for Xe+La " << Config::BeamMomentum << "A GeV/c ---" << std::endl;

    auto now = std::chrono::system_clock::now();
    std::time_t now_c = std::chrono::system_clock::to_time_t(now);
    std::tm* now_tm = std::localtime(&now_c);
    std::stringstream ss_time;

    ss_time << std::put_time(now_tm, "%H_%M_%S_%d_%m_%Y");
    // Construct the base file name for output files
    std::string eventName = "EventXeLaMag_" + std::to_string(Config::BeamMomentum) + "A";
    std::string baseFileName = outputDir + "/" + ss_time.str() + "_" + eventName;
    // Open a text file to log rapidity details for the first 1000 tracks
    std::string rapidityTxtPath = baseFileName + "_rapidity_details.txt";
    std::ofstream rapFile(rapidityTxtPath);
    int rap_print_counter = 0; 
    
    // Open a text file to export proton px and py for inspection if needed
    std::string pxpyTxtPath = baseFileName + "_protons_pxpy.txt";
    std::ofstream pxpyFile(pxpyTxtPath);
    if (pxpyFile.is_open()) {
        pxpyFile << "px\tpy\n";
    }

    // Ensure the rapidity file is open before writing
    TChain chain("event_tree");
    std::ifstream file(fileList);
    if (!file.is_open()) {
        std::cerr << "ERROR: Cannot open the file with paths list!" << std::endl;
        return;
    }
    // Read each line from the file and add it to the TChain
    std::string filePath; // Variable to hold each file path
    while (std::getline(file, filePath)) {
        if (!filePath.empty()) chain.Add(filePath.c_str());
    }
    file.close();
    // Set all branches to active for reading
    chain.SetBranchStatus("*", 1); // Enable all branches for reading
    EventClass::EventXeLaMag* pevent = nullptr; // Pointer to hold the event data
    chain.SetBranchAddress("EventXeLaMag", &pevent);
    // Get the total number of entries in the TChain
    Long64_t nentries = chain.GetEntries(); // Total number of events in the chain
    std::cout << "Analyzing " << nentries << " events..." << std::endl; 
    // Initialize histogram manager and Bethe-Bloch wrapper
    HistogramManager hists;
    Dedx::BetheBlochWrapper bbWrapper; // Initialize Bethe-Bloch wrapper for dE/dx calculations
    // Initialize counters for events and tracks
    unsigned int nEvents = 0; // Total number of events processed
    unsigned int nEvents_PSD = 0; // Number of events passing the PSD cut
    unsigned int nEvents_VertexZ = 0; // Number of events passing the Vertex Z cut
    unsigned int nEvents_tracksratio = 0; // Number of events passing the tracks ratio cut

    unsigned int nTracks = 0; // Total number of tracks processed   
    unsigned int nTracks_VTPC12 = 0; // Number of tracks passing the VTPC1+VTPC2 clusters cut
    unsigned int nTracks_ClustersAll = 0; // Number of tracks passing the total clusters cut
    unsigned int nTracks_dEdx = 0; // Number of tracks passing the dE/dx clusters cut
    unsigned int nTracks_pointRatio = 0; // Number of tracks passing the potential clusters ratio cut
    unsigned int nTracks_ImpactParameter = 0; // Number of tracks passing the impact parameter cut
    unsigned int nTracks_protons = 0; // Number of tracks identified as protons

    unsigned int nEvents_before_centrality = 0; // Number of events before applying the centrality cut
    unsigned int nTracksFit_before_centrality = 0; // Number of fitted tracks before applying the centrality cut
    unsigned int nEvents_after_centrality = 0; // Number of events after applying the centrality cut
    unsigned int nTracksFit_after_centrality = 0; // Number of fitted tracks after applying the centrality cut

    double mp = 0.93827;
    double pL = Config::BeamMomentum; // Beam momentum in A GeV/c
    double E_L = std::sqrt(mp*mp + pL*pL); // Beam energy in the LAB frame (E_L = sqrt(m^2 + p^2))
    double E_CM = 0.5 * std::sqrt(2.0 * mp * (mp + E_L));

    /* 
     * KINEMATICS STEP-BY-STEP: Lorentz Transformation from LAB to CM Frame
     * 
     * 1. We use the Mandelstam invariant 's', which is conserved across reference frames:
     *    s = (p_beam + p_target)^2 = (E_L + m_N)^2 - (p_L + 0)^2 = 4 * (E_CM)^2
     *
     * 2. Expanding the LAB side, recalling the energy-momentum relation (E_L^2 - p_L^2 = m_N^2):
     *    E_L^2 - p_L^2 + m_N^2 + 2 * m_N * E_L = 4 * (E_CM)^2
     *    m_N^2 + m_N^2 + 2 * m_N * E_L = 4 * (E_CM)^2
     *    2 * m_N * (m_N + E_L) = 4 * (E_CM)^2
     *
     * 3. The total CM energy per nucleon pair (\sqrt{S_{NN}}) is therefore:
     *    \sqrt{S_{NN}} = 2 * E_CM = \sqrt{2 * m_N * (m_N + E_L)}
     *
     * 4. The Lorentz boost factor \gamma (from LAB to CM) is defined as E_CM / m_N.
     *
     * 5. Rapidity (y) is related to \gamma by \gamma = \cosh(y). Thus, the beam rapidity 
     *    in the CM frame is:
     *    y_{beam}^{CM} = \text{acosh}(E_CM / m_N)
     *
     * 6. Because rapidity is simply additive under Lorentz boosts along the z-axis, 
     *    we can translate any particle's LAB rapidity to CM rapidity by simply subtracting:
     *    y_{CM} = y_{LAB} - y_{beam}^{CM}
     */
    double Y_beam_shift = std::acosh(E_CM / mp); // Beam rapidity shift in the CM frame

    // F2(M) I
    const int M_MAX = 150; // Maximum number of bins for F2(M) calculation
    const double L_window = 3.0; // The characteristic length scale for pair counting in p_x, p_y space
    const double R_0 = L_window / std::sqrt(M_PI); // Radius corresponding to the area of a circle with area L_window^2
    const double m_f = 0.27; // Fractional offset to adjust the binning for pair counting
    

    std::vector<double> sum_N2(M_MAX, 0.0); // Accumulates the total number of pairs N2(M) across all events
    std::vector<double> sum_N2_sq(M_MAX, 0.0); // Accumulates the square of N2(M) for variance calculation
    double sum_mul = 0.0; // Accumulates the total multiplicity (number of tracks) across all events
    double sum_mul_sq = 0.0;// Accumulates the square of multiplicity for variance calculation
    long long f2_valid_events = 0; // Counts the number of events that pass all cuts and are valid for F2(M) calculation
   
    // Save rapidity details to TXT
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
    
    /* Save rapidity details to TXT */
    std::cout << "Step 1: Calculating dynamic centrality limit on CLEAN events..." << std::endl; // Calculate dynamic centrality limit based on PSD energy distribution
    for (Long64_t ievents = 0; ievents < nentries; ievents++) { 
        chain.GetEntry(ievents); // Load the event data for the current entry in the TChain

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

        hists.h_PSD_T2->Fill(pevent->energyPSDSelectedModules); // Fill the histogram for PSD energy distribution to determine centrality thresholds
    }

    std::vector<double> cent_limits(4, 0.0); // Vector to hold the calculated centrality limits for 5%, 10%, 15 % and 20 % centrality
    std::vector<int> cent_events(4, 0); // Vector to hold the number of events corresponding to each centrality limit for plotting consistency
    double fracs[4] = {0.05, 0.10, 0.15, 0.20}; // Fractions for 5%, 10%, 15% and 20% centrality thresholds

    TH1* cumul = hists.h_PSD_T2->GetCumulative(); // Get the cumulative distribution of the PSD energy histogram
    TH1* cumul_raw = hists.h_PSD_T2->GetCumulative(); // Get a raw copy of the cumulative distribution for event counting

    if (cumul->GetMaximum() > 0) {
        cumul->Scale(0.423 / cumul->GetMaximum()); // Scale the cumulative histogram to a maximum of 0.423 for normalization
    }

    // Calculate centrality limits based on the cumulative distribution of PSD energy
    for (int i = 0; i < 4; ++i) {
        int bin = cumul->FindFirstBinAbove(fracs[i]); // Find the first bin in the cumulative histogram that exceeds the specified fraction (5%, 10%, 15 and 20 %)
        if (bin > 0) {
            cent_limits[i] = cumul->GetXaxis()->GetBinCenter(bin); // Store the corresponding PSD energy value for the centrality limit
            cent_events[i] = cumul_raw->GetBinContent(bin); // Store the number of events corresponding to that centrality limit for plotting consistency
        }
    }
    double dynamic_limit_20 = cent_limits[1]; 
    delete cumul;
    delete cumul_raw;

    std::cout << ">> Dynamic 20% centrality threshold calculated as: " << dynamic_limit_20 << " GeV\n" << std::endl;

    std::cout << "applying cuts" << std::endl;
    for (Long64_t ievents = 0; ievents < nentries; ievents++) {
        chain.GetEntry(ievents);
        //copy of the run number filter from above to avoid processing unwanted runs
        if (ievents > 0 && ievents % 100000 == 0) {
            double percent = (static_cast<double>(ievents) / nentries) * 100.0;
            std::cout << "Processed: " << ievents << " / " << nentries << " (" << static_cast<int>(percent) << "%)" << std::endl;
        }

        if (!pevent) continue;
        //copy of the run number filter from above to avoid processing unwanted runs
        if (pevent->run_number == 34948 || pevent->run_number == 34976 ||
            pevent->run_number == 35093 || pevent->run_number == 35142) {
            continue;
        }

        nEvents++; // Increment the total number of events processed

        /* Fill histograms for all events */
        hists.h_Vz_all->Fill(pevent->VertexZ); 
        hists.h2_TracksInFit_vs_PSD_all->Fill(pevent->energyPSDSelectedModules, pevent->nTracksFit);
        hists.h2_Tracks_vs_PSD_all->Fill(pevent->energyPSDSelectedModules, pevent->nTracksAll);
        
        // LEVEL 1: EVENT CUTS
        if (EventCuts::PassVertexZ(pevent->VertexZ)) { // Check if the event passes the Vertex Z cut    
            nEvents_VertexZ++; // Increment the count of events passing the Vertex Z cut
            hists.h_Vz_cut->Fill(pevent->VertexZ); 
            hists.h2_TracksInFit_vs_PSD_after_Vz->Fill(pevent->energyPSDSelectedModules, pevent->nTracksFit);
            hists.h2_PSD_Peripheral_vs_Selected->Fill(pevent->energyPSDSelectedModules, pevent->energyPSDPeripheralModules);
            hists.h_PSDperipheral_all->Fill(pevent->energyPSDPeripheralModules);
            // Fill histograms for events passing the Vertex Z cut
            
            if (pevent->energyPSDPeripheralModules > Config::PSD_per_cut) { // Check if the event passes the PSD peripheral energy cut
                nEvents_PSD++;
                hists.h_PSDperipheral_cut->Fill(pevent->energyPSDPeripheralModules);
                hists.h2_PSD_Peripheral_vs_Selected_cut->Fill(pevent->energyPSDSelectedModules, pevent->energyPSDPeripheralModules);
                hists.h2_TracksInFit_vs_PSD_after_PSD->Fill(pevent->energyPSDSelectedModules, pevent->nTracksFit);
                hists.h2_TracksRatio_all->Fill(pevent->nTracksAll, pevent->nTracksFit);
                // Fill histograms for events passing the PSD peripheral energy cut
                
                if (EventCuts::PassTracksRatio(pevent->nTracksAll, pevent->nTracksFit)) {
                    nEvents_tracksratio++;
                    hists.h2_TracksRatio_cut->Fill(pevent->nTracksAll, pevent->nTracksFit);
                    hists.h2_Tracks_vs_PSD_cut->Fill(pevent->energyPSDSelectedModules, pevent->nTracksAll);
                    hists.h2_TracksInFit_vs_PSD_cut->Fill(pevent->energyPSDSelectedModules, pevent->nTracksFit);

                    nEvents_before_centrality++; // Increment the count of events before applying the centrality cut
                    nTracksFit_before_centrality += pevent->nTracksFit;

                    if (pevent->energyPSDSelectedModules < dynamic_limit_20) { // Check if the event passes the dynamic centrality cut (120% most central events)

                        hists.h2_TracksInFit_vs_PSD_Central->Fill(pevent->energyPSDSelectedModules, pevent->nTracksFit);// Fill histogram for events passing the dynamic centrality cut

                        nEvents_after_centrality++;
                        nTracksFit_after_centrality += pevent->nTracksFit;

                        std::vector<std::pair<double, double>> event_protons; //vector instead txt file

                        for (const auto& tracks : pevent->tracks) { // Loop over all tracks in the event
                            nTracks++;

                            // LEVEL 2: TRACK QUALITY CUTS
                            unsigned int cl_dEdx = tracks.clustersdEdx; // Number of clusters used for dE/dx measurement
                            unsigned int cl_VTPC_sum = tracks.clustersVTPC1 + tracks.clustersVTPC2; // Sum of clusters in VTPC1 and VTPC2
                            unsigned int cl_Pot = tracks.clustersPotentialAll; // Number of potential clusters in all TPCs
                            unsigned int cl_All = tracks.clustersAll; // Total number of clusters in all TPCs
                            double ratio = (cl_Pot > 0) ? static_cast<double>(cl_All) / static_cast<double>(cl_Pot) : 0.0; // Calculate the ratio of total clusters to potential clusters, ensuring no division by zero

                            hists.h_clusters_VTPC_sum_all->Fill(cl_VTPC_sum); 
                            if (cl_VTPC_sum <= Config::VTPC_sum_min) continue; // Check if the sum of clusters in VTPC1 and VTPC2 is above the minimum threshold
                            nTracks_VTPC12++;
                            hists.h_clusters_VTPC_sum_cut->Fill(cl_VTPC_sum); // Fill histogram for tracks passing the VTPC1+VTPC2 clusters cut

                            hists.h_clusters_All_all->Fill(cl_All);
                            if (cl_All <= Config::ClustersAll_min) continue; // Check if the total number of clusters is above the minimum threshold
                            nTracks_ClustersAll++;
                            hists.h_clusters_All_cut->Fill(cl_All); // Fill histogram for tracks passing the total clusters cut

                            hists.h_clusters_dEdx_all->Fill(cl_dEdx);
                            if (cl_dEdx <= Config::dEdx_min) continue; // Check if the number of clusters used for dE/dx measurement is above the minimum threshold
                            nTracks_dEdx++;
                            hists.h_clusters_dEdx_cut->Fill(cl_dEdx); // Fill histogram for tracks passing the dE/dx clusters cut

                            hists.h_clusters_Ratio_all->Fill(ratio);
                            hists.h_clusters_PotAll_all->Fill(cl_Pot);
                            if (ratio < Config::Ratio_min || ratio > Config::Ratio_max) continue; // Check if the ratio of total clusters to potential clusters is within the specified range
                            nTracks_pointRatio++;
                            hists.h_clusters_Ratio_cut->Fill(ratio); 
                            hists.h_clusters_PotAll_cut->Fill(cl_Pot);

                            hists.h2_bx_by_all->Fill(tracks.bx, tracks.by);
                            if (std::abs(tracks.bx) > Config::ImpactParam_bx_max || std::abs(tracks.by) > Config::ImpactParam_by_max) continue; // Check if the impact parameters (bx, by) are within the specified maximum limits
                            nTracks_ImpactParameter++;
                            hists.h2_bx_by_cut->Fill(tracks.bx, tracks.by);

                            double ptot = std::sqrt(tracks.px * tracks.px + tracks.py * tracks.py + tracks.pz * tracks.pz);
                            double log_ptot = std::log10(ptot); // Calculate the logarithm of the total momentum for dE/dx vs log10(p_tot) analysis

                            /*
                             * : Proton Selection via Bethe-Bloch (BB) Curves
                             * 
                             * 1. Charged particles ionize the TPC gas, losing energy (dE/dx) which is 
                             *    measured alongside their total momentum (p_tot).
                             *
                             * 2. When plotting dE/dx vs \log_{10}(p_{tot}), particles form distinct bands 
                             *    (pions, kaons, protons) that follow theoretical Bethe-Bloch functions.
                             *
                             * 3. At higher momenta, the proton and kaon/pion bands begin to overlap. To cleanly 
                             *    select protons and reject kaons, we define a dynamic upper threshold.
                             *
                             * 4. This threshold is positioned as a fraction of the distance between the theoretical 
                             *    proton curve (BB_p) and kaon curve (BB_K):
                             *    Upper Limit = BB_p + 0.15 * (BB_K - BB_p)
                             *
                             * 5. The coefficient 0.15 is chosen empirically. If a particle has a positive dE/dx 
                             *    that falls strictly below this upper limit, it is accepted as a proton candidate.
                             */
                            double bb_proton_dedx = bbWrapper(3, ptot); // Calculate the theoretical dE/dx for protons using the Bethe-Bloch formula
                            double bb_kaon_dedx = bbWrapper(2, ptot); // Calculate the theoretical dE/dx for kaons using the Bethe-Bloch formula
                            double delta_kp = bb_kaon_dedx - bb_proton_dedx; // Calculate the difference between the kaon and proton dE/dx values

                            double E_track = std::sqrt(mp*mp + ptot*ptot);
                            double Y_LAB_track = 0.5 * std::log((E_track + tracks.pz) / (E_track - tracks.pz)); // Calculate the rapidity of the track in the LAB frame
                            double Y_CM_track = Y_LAB_track - Y_beam_shift; // Transform the track's rapidity to the CM frame by subtracting the beam rapidity shift

                            if (rapFile.is_open() && rap_print_counter < 1000) { // Log the rapidity details for the first 1000 tracks to the text file
                                rapFile << std::setw(15) << Y_LAB_track //std::setw(15)  - Set the width of the output field to 15 characters for alignment
                                        << std::setw(15) << Y_beam_shift 
                                        << std::setw(15) << Y_CM_track << "\n";
                                rap_print_counter++; // Increment the counter for printed rapidity details
                            }

                            // Fill histograms after quality cuts
                            if (tracks.dEdx > 0) {
                                hists.h_dedx_ptot_pos_qual->Fill(log_ptot, tracks.dEdx);
                                hists.h2_px_py_pos_qual->Fill(tracks.px, tracks.py);
                            } else if (tracks.dEdx < 0) {
                                hists.h_dedx_ptot_neg_qual->Fill(log_ptot, std::abs(tracks.dEdx));
                                hists.h2_px_py_neg_qual->Fill(tracks.px, tracks.py);
                            }

                            // LEVEL 3: TRACK MOMENTUM CUTS
                            bool pass_log_ptot = (log_ptot >= 0.55 && log_ptot <= 2.0); // Check if the logarithm of the total momentum is within the specified range for proton selection
                            bool pass_px = (tracks.px > -1.5 && tracks.px < 1.5); // Check if the track's x-component of momentum is within the specified range
                            bool pass_py = (tracks.py > -1.5 && tracks.py < 1.5); // Check if the track's y-component of momentum is within the specified range

                            if (!pass_log_ptot || !pass_px || !pass_py) continue;

                            // Fill histograms after momentum cuts
                            if (tracks.dEdx > 0) {
                                hists.h_dedx_ptot_pos_mom->Fill(log_ptot, tracks.dEdx);
                                hists.h2_px_py_pos_mom->Fill(tracks.px, tracks.py);
                            } else if (tracks.dEdx < 0) {
                                hists.h_dedx_ptot_neg_mom->Fill(log_ptot, std::abs(tracks.dEdx));
                                hists.h2_px_py_neg_mom->Fill(tracks.px, tracks.py);
                            }

                            // LEVEL 4: PROTON IDENTIFICATION CUTS PID
                            if (tracks.dEdx <= 0) continue; // Only consider tracks with positive dE/dx values for proton selection

                            double upper_limit = bb_proton_dedx + 0.15 * delta_kp;
                            if (tracks.dEdx > upper_limit) continue; // Drop track if above proton band

                            nTracks_protons++; // Increment the proton track count

                            // Fill histograms after PID
                            hists.h_dedx_ptot_protons_id->Fill(log_ptot, tracks.dEdx);
                            hists.h2_px_py_protons_id->Fill(tracks.px, tracks.py);
                            hists.h_Y_CM_protons_id->Fill(Y_CM_track);

                            // LEVEL 5: RAPIDITY CUT  _rap (FINAL)
                            bool pass_rapidity = (std::abs(Y_CM_track) <= 0.75); // Check if the track's rapidity in the CM frame is within the specified range
                            if (!pass_rapidity) continue;

                            // Fill histograms after rapidity cut
                            hists.h_dedx_ptot_protons_rap->Fill(log_ptot, tracks.dEdx); // Fill the histogram for dE/dx vs log10(p_tot) for selected proton candidates
                            hists.h2_px_py_protons_rap->Fill(tracks.px, tracks.py);
                            hists.h_Y_CM_protons_rap->Fill(Y_CM_track);

                            // Optional export of final proton px and py to text file for inspection
                            if (pxpyFile.is_open()) {
                                pxpyFile << tracks.px << "\t" << tracks.py << "\n";
                            }

                            event_protons.push_back({tracks.px, tracks.py}); // I don't use txt file as primary container because vector form is way faster
                        }

                        // F2(M) II
                        /*
                         *  Correlation Integral ("Moving Circles" Method)
                         * 
                         * Intermittency analysis searches for scaling behavior in the 2nd factorial moment.
                         *
                         * 1. Standard approach: Divide the (p_x, p_y) space into a fixed M \times M grid. 
                         *    This suffers from binning artifacts (grid-boundary effects).
                         *
                         * 2. Our approach (Correlation Integral C(R)): Instead of rigid grids, we measure 
                         *    the density of pairs as a function of their continuous distance R.
                         *
                         * 3. For every unique pair of protons in the event, we compute the 2D Euclidean 
                         *    distance in transverse momentum space: 
                         *    dist = \sqrt{\Delta p_x^2 + \Delta p_y^2}
                         *
                         * 4. We map this 'dist' directly to an equivalent grid subdivision M index:
                         *    idx = \lfloor R_0 / dist + m_f \rfloor
                         *    (where R_0 scales with the phase space window, and m_f is an offset factor).
                         *
                         * 5. We increment N_2(M), which counts how many pairs fall within a radius R_M.
                         *
                         * 6. F_2(M) is later analytically recovered from this pair count via the relation:
                         *    F_2(M) = \frac{2 \cdot M^2}{\langle mul \rangle^2} \langle N_2(M) \rangle
                         */
                        double current_mul = event_protons.size(); // Current event multiplicity (number of selected protons)
                        sum_mul += current_mul;
                        sum_mul_sq += current_mul * current_mul; // Accumulate the square of the multiplicity for variance calculation
                        f2_valid_events++;

                        std::vector<double> event_N2(M_MAX, 0.0); // Initialize the pair count N2(M) for the current event

                        for (size_t i = 0; i < event_protons.size(); ++i) { // Loop over all unique pairs of protons in the event to calculate their distances and update N2(M)
                            for (size_t j = i + 1; j < event_protons.size(); ++j) {
                                double dx = event_protons[i].first - event_protons[j].first; // Difference in p_x between the two protons
                                double dy = event_protons[i].second - event_protons[j].second; // Difference in p_y between the two protons
                                double dist = std::sqrt(dx*dx + dy*dy);

                                if (dist > 0.0) {
                                    int idx = std::floor(R_0 / dist + m_f); // Map the distance to an equivalent M index for pair counting
                                    if (idx > M_MAX) idx = M_MAX; // Ensure the index does not exceed the maximum M value
                                    
                                    for (int m = 0; m < idx; ++m) { // Increment the pair count for all M values up to the calculated index
                                        event_N2[m] += 1.0; // Each pair contributes to all smaller M bins, reflecting the cumulative nature of the correlation integral
                                    }
                                }
                            }
                        }

                        for (int m = 0; m < M_MAX; ++m) { // Accumulate the total pair counts and their squares for variance calculation across all events
                            sum_N2[m] += event_N2[m]; // Accumulate the total number of pairs N2(M) across all events
                            sum_N2_sq[m] += event_N2[m] * event_N2[m]; // Accumulate the square of N2(M) for variance calculation
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
    
    if (pxpyFile.is_open()) {
        pxpyFile.close();
    }

    std::string f2TxtPath = baseFileName + "_F2_results.txt";
    std::ofstream f2File(f2TxtPath);
    
    // F2(M) III 
    /*
     * ERROR PROPAGATION: Statistical Uncertainty for F2(M)
     * 
     * F_2(M) is calculated as a ratio of two averaged random variables: the average 
     * number of pairs \langle N_2 \rangle and the square of the average multiplicity \langle mul \rangle^2.
     * We must properly propagate the statistical errors of both variables.
     *
     * 1. Calculate the variance for both event multiplicity (mul) and pair counts (N_2):
     *    Var(X) = \langle X^2 \rangle - \langle X \rangle^2
     *
     * 2. The standard error of the mean (for N valid events) is:
     *    err_X = \sqrt{ Var(X) / N_{events} }
     *
     * 3. According to the standard error propagation formula for a ratio F = A / B^2 
     *    (treating A and B as independent variables for simplicity):
     *    \frac{err_F}{F} = \sqrt{ \left(\frac{err_A}{A}\right)^2 + \left(2 \cdot \frac{err_B}{B}\right)^2 }
     *
     *    where:
     *    A = \langle N_2 \rangle
     *    B = \langle mul \rangle
     *
     * 4. This ensures that the final error bars on the F_2(M) plot accurately reflect 
     *    the statistical fluctuations of the recorded event sample.
     */
    if (f2File.is_open() && f2_valid_events > 0) { // Ensure the output file is open and there are valid events for F2(M) calculation
        double avg_mul = sum_mul / f2_valid_events;
        double avg_mul_sq = sum_mul_sq / f2_valid_events;
        
        double var_mul = avg_mul_sq - (avg_mul * avg_mul);
        if (var_mul < 0) var_mul = 0.0; 
        double err_mul = std::sqrt(var_mul) / std::sqrt(f2_valid_events);

        f2File << "M\tF2_M\terr_F2_M\tavg_N2\terr_N2\tavg_mul\terr_mul\n";

        for (int m = 0; m < M_MAX; ++m) { // Loop over all M values to calculate F2(M) and its statistical error
            int M_val = m + 1; 
            
            double avg_N2 = sum_N2[m] / f2_valid_events; // Calculate the average number of pairs N2(M) for the current M value
            double avg_N2_sq = sum_N2_sq[m] / f2_valid_events; // Calculate the average of the square of N2(M) for variance calculation
            
            double var_N2 = avg_N2_sq - (avg_N2 * avg_N2);
            if (var_N2 < 0) var_N2 = 0.0;
            double err_N2 = std::sqrt(var_N2) / std::sqrt(f2_valid_events);

            double F2_M = 0.0; // Initialize F2(M) for the current M value
            double err_F2_M = 0.0;

            if (avg_mul > 0 && avg_N2 > 0) {
                F2_M = (2.0 * M_val * M_val * avg_N2) / (avg_mul * avg_mul);
                
                double rel_err_N2 = err_N2 / avg_N2;
                double rel_err_mul = err_mul / avg_mul;
                
                err_F2_M = F2_M * std::sqrt(rel_err_N2 * rel_err_N2 + 4.0 * rel_err_mul * rel_err_mul); // Propagate the errors for F2(M) using the standard error propagation formula for a ratio
            }
            // Write the results for the current M value to the output file, including F2(M), its error, average N2, its error, average multiplicity, and its error
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
    std::cout << "=> Proton px and py coordinates exported to file: " << pxpyTxtPath << std::endl;

    std::cout << "\n=== Centrality Stats ===" << std::endl;
    std::cout << "Events surviving 0-20% cut: " << nEvents_after_centrality << std::endl;
    if (nEvents_before_centrality > 0)
        std::cout << "Avg tracks/event (before cut): " << (double)nTracksFit_before_centrality / nEvents_before_centrality << std::endl;
    if (nEvents_after_centrality > 0)
        std::cout << "Avg tracks/event (after cut):  " << (double)nTracksFit_after_centrality / nEvents_after_centrality << std::endl;
    std::cout << "========================\n" << std::endl;

    Plotter::DrawAndSaveAll(hists, baseFileName, cent_limits, cent_events); // Call the Plotter to draw and save all histograms and plots based on the analysis results
}