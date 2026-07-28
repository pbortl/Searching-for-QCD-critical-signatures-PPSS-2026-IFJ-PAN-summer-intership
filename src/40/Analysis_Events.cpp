#include "Analysis_Events.h"
#include "Config.h"
#include "HistogramManager.h"
#include "EventCuts.h"
#include "EventPlotter.h"

#include <iostream>
#include <fstream>
#include <string>
#include <chrono>
#include <ctime>
#include <iomanip>
#include <sstream>
#include <vector>

#include "TROOT.h"
#include "TChain.h"
#include "TSystem.h"
#include "TError.h"
#include "TFile.h"

Analysis_Events::Analysis_Events(const std::string& fileListPath) : fileList(fileListPath) {}

void Analysis_Events::LoadLibraries() {
    static bool initialized = false;
    if (!initialized) {
        std::cout << "Initializing libraries..." << std::endl;
        const char* libs[] = {
            "/home/p/Documents/PPSS_tools/Tools/Event.so",
            "/home/p/Documents/PPSS_tools/Tools/CutsMap.noDict.so",
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

void Analysis_Events::Run() {
    gROOT->SetBatch(kTRUE);
    gErrorIgnoreLevel = kError; 
    LoadLibraries();

    std::string outputDir;
    std::cout << "Enter the path to the output directory (e.g., /home/user/results): "; 
    std::getline(std::cin, outputDir);

    if (!outputDir.empty() && outputDir.back() == '/') { 
        outputDir.pop_back();
    }
    std::cout << "--- Applying ONLY EVENT CUTS for Xe+La " << Config::BeamMomentum << "A GeV/c ---" << std::endl;

    auto now = std::chrono::system_clock::now();
    std::time_t now_c = std::chrono::system_clock::to_time_t(now);
    std::tm* now_tm = std::localtime(&now_c);
    std::stringstream ss_time;
    ss_time << std::put_time(now_tm, "%H_%M_%S_%d_%m_%Y");
    
    std::string baseFileName = outputDir + "/" + ss_time.str() + "_EventXeLaMag_EventCutsOnly_" + std::to_string(Config::BeamMomentum) + "A";

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
    std::cout << "Analyzing " << nentries << " events (Fast Mode)..." << std::endl; 

    HistogramManager hists;

    unsigned int nEvents = 0;
    unsigned int nEvents_VertexZ = 0;
    unsigned int nEvents_PSD = 0;
    unsigned int nEvents_tracksratio = 0;
    unsigned int nEvents_after_centrality = 0;

    // STEP 1: CALCULATE CENTRALITY ON CLEAN EVENTS
    std::cout << "Step 1: Calculating dynamic centrality limit on CLEAN events..." << std::endl;
    for (Long64_t ievents = 0; ievents < nentries; ievents++) { 
        chain.GetEntry(ievents);

        if (!pevent || pevent->run_number == 34948 || pevent->run_number == 34976 ||
            pevent->run_number == 35093 || pevent->run_number == 35142) continue;

        if (!EventCuts::PassVertexZ(pevent->VertexZ)) continue;
        if (pevent->energyPSDPeripheralModules <= Config::PSD_per_cut) continue;
        if (!EventCuts::PassTracksRatio(pevent->nTracksAll, pevent->nTracksFit)) continue;
        //first cuttings
        hists.h_PSD_T2->Fill(pevent->energyPSDSelectedModules);
    }

    std::vector<double> cent_limits(4, 0.0);
    std::vector<int> cent_events(4, 0);
    double fracs[4] = {0.05, 0.10, 0.15, 0.20};

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
    double dynamic_limit_20 = cent_limits[3]; 
    delete cumul; delete cumul_raw;

    std::cout << ">> Dynamic 20% centrality threshold: " << dynamic_limit_20 << " GeV\n" << std::endl;

    // STEP 2: FILL EVENT HISTOGRAMS
    std::cout << "Step 2: Filling event histograms..." << std::endl;
    for (Long64_t ievents = 0; ievents < nentries; ievents++) {
        chain.GetEntry(ievents);

        if (!pevent || pevent->run_number == 34948 || pevent->run_number == 34976 ||
            pevent->run_number == 35093 || pevent->run_number == 35142) continue;

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

                    if (pevent->energyPSDSelectedModules < dynamic_limit_20) {
                        nEvents_after_centrality++;
                        hists.h2_TracksInFit_vs_PSD_Central->Fill(pevent->energyPSDSelectedModules, pevent->nTracksFit);
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

    std::cout << "Event analysis finished! Calling EventPlotter..." << std::endl;

    EventPlotter::DrawAndSaveEventCuts(hists, baseFileName, cent_limits, cent_events);
}