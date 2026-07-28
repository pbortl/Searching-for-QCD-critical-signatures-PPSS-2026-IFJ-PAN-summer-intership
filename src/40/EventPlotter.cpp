#include "EventPlotter.h"
#include "Config.h"

#include <iostream>
#include <fstream>
#include <cmath>
#include <functional>
#include <vector>
#include <string>
#include <iomanip>

#include "TCanvas.h"
#include "TF1.h"
#include "TLine.h"
#include "TFile.h"
#include "TStyle.h"
#include "TPaveText.h"
#include "TAxis.h"

namespace {
    void DrawTracksRatioLines() {
        TLine lineT1(Config::P1_x, Config::P1_y, Config::P2_x, Config::P2_y);
        lineT1.SetLineColor(kRed); lineT1.SetLineWidth(2); lineT1.DrawClone("SAME");
        
        TLine lineT2(Config::P2_x, Config::P1_y, Config::P2_x, Config::P2_y);
        lineT2.SetLineColor(kRed); lineT2.SetLineWidth(2); lineT2.DrawClone("SAME");
        
        TLine lineT3(Config::P1_x, Config::P1_y, Config::P2_x, Config::P1_y);
        lineT3.SetLineColor(kRed); lineT3.SetLineWidth(2); lineT3.DrawClone("SAME");

        double slope2 = (Config::P4_y - Config::P3_y) / (Config::P4_x - Config::P3_x);
        double end_x2 = Config::P3_x + (350.0 - Config::P3_y) / slope2;
        double end_y2 = 350.0;

        if (end_x2 > 1500.0) {
            end_x2 = 1500.0;
            end_y2 = Config::P3_y + slope2 * (1500.0 - Config::P3_x);
        }
        
        TLine lineVector(Config::P3_x, Config::P3_y, end_x2, end_y2);
        lineVector.SetLineColor(kRed); lineVector.SetLineWidth(2); lineVector.DrawClone("SAME");
    }

    void DrawCentralityLines1D(double yMax, const std::vector<double>& cent_limits, const std::vector<int>& cent_events, int total_evts) {
        for(size_t i = 0; i < cent_limits.size(); i++) {
            if(cent_limits[i] > 0) {
                TLine line(cent_limits[i], 0, cent_limits[i], yMax);
                line.SetLineColor(kRed);
                line.SetLineStyle(2); 
                line.SetLineWidth(2);
                line.DrawClone("SAME");
            }
        }
        
        TPaveText pt(0.65, 0.35, 0.85, 0.65, "NDC");
        pt.SetFillColor(kWhite);
        pt.SetBorderSize(1);
        pt.SetTextAlign(12);
        pt.AddText(Form("total evets: %d", total_evts));
        if (cent_events.size() >= 1) pt.AddText(Form("0-5%%: %d", cent_events[0]));
        if (cent_events.size() >= 2) pt.AddText(Form("5-10%%: %d", cent_events[1] - cent_events[0]));
        pt.DrawClone("SAME");
    }

    void DrawCentralityLines2D(double yMax, const std::vector<double>& cent_limits, const std::vector<int>& cent_events, int total_evts) {
        for(size_t i = 0; i < cent_limits.size(); i++) {
            if(cent_limits[i] > 0) {
                TLine line(cent_limits[i], 0, cent_limits[i], yMax);
                line.SetLineColor(kRed);
                line.SetLineStyle(2); 
                line.SetLineWidth(2);
                line.DrawClone("SAME");
            }
        }

        TPaveText pt(0.65, 0.65, 0.88, 0.88, "NDC");
        pt.SetFillColor(kWhite);
        pt.SetBorderSize(1);
        pt.SetTextAlign(12);
        pt.AddText(Form("total evets: %d", total_evts));
        if (cent_events.size() >= 1) pt.AddText(Form("0-5%%: %d", cent_events[0]));
        if (cent_events.size() >= 2) pt.AddText(Form("5-10%%: %d", cent_events[1] - cent_events[0]));
        pt.DrawClone("SAME");
    }
}

void EventPlotter::DrawAndSaveEventCuts(HistogramManager& hists, 
                                        const std::string& baseFileName, 
                                        const std::vector<double>& cent_limits, 
                                        const std::vector<int>& cent_events) {
    std::cout << "Generating Event Cuts multi-page PDF and saving to ROOT file..." << std::endl;

    gStyle->SetOptStat(1111);

    std::string txtFilePath = baseFileName + "_vertex_fit.txt";
    std::string pdfFilePath = baseFileName + "_event_cuts.pdf";
    std::string rootOutputFilePath = baseFileName + "_event_cuts.root";

    TFile *fOut = new TFile(rootOutputFilePath.c_str(), "RECREATE");

    hists.hist_events->SetMinimum(0);

    std::string suffix = " for Xe+La " + std::to_string(Config::BeamMomentum) + "A GeV/c";

    // Set histogram titles matching Plotter style
    hists.h_Vz_all->SetTitle(("distrib vz before the cut" + suffix).c_str());
    hists.h_Vz_cut->SetTitle(("distrib vz cutted" + suffix).c_str());
    hists.h2_PSD_Peripheral_vs_Selected->SetTitle(("PSD energy peripheral modules before cut" + suffix).c_str());
    hists.h2_PSD_Peripheral_vs_Selected_cut->SetTitle(("peripheral energy after the cut" + suffix).c_str());
    hists.h2_TracksRatio_all->SetTitle(("ntracks in fit vs ntracks all before cut" + suffix).c_str());
    hists.h2_TracksRatio_cut->SetTitle(("ntracks in fit vs ntracks after linear cut" + suffix).c_str());
    
    hists.h2_TracksInFit_vs_PSD_all->SetTitle(("PSD energy before any cuts" + suffix).c_str());
    hists.h2_TracksInFit_vs_PSD_after_Vz->SetTitle(("PSD energy after vertexz cut" + suffix).c_str());
    hists.h2_TracksInFit_vs_PSD_after_PSD->SetTitle(("PSD energy after vertexz cut and periph energy cut" + suffix).c_str());
    
    hists.h_PSD_T2->SetTitle(("events vs E_{PSD} (cut)" + suffix).c_str());
    hists.h2_TracksInFit_vs_PSD_cut->SetTitle(("track (fit) vs E_{PSD} (cut)" + suffix).c_str());
    hists.h2_TracksInFit_vs_PSD_Central->SetTitle(("fitted tracks vs PSD energy AFTER 0-10% cut" + suffix).c_str());
    
    hists.hist_events->SetTitle(("nb of events after each event cut" + suffix).c_str());

    // Gaussian fit for Vertex Z
    TF1 gaussFit("gaussFit", "gaus", -606, -600);
    gaussFit.SetParameters(100000, Config::z_peak, 1.0);
    hists.h_Vz_all->Fit(&gaussFit, "R0");

    std::ofstream txtOut(txtFilePath);
    if (txtOut.is_open()) {
        txtOut << "--- GAUSSIAN FIT RESULTS (VertexZ) ---\n";
        txtOut << "Amplitude: " << gaussFit.GetParameter(0) << "\n";
        txtOut << "Mean:      " << gaussFit.GetParameter(1) << " cm\n";
        txtOut << "Sigma:     " << gaussFit.GetParameter(2) << " cm\n";
        txtOut.close();
    }

    TCanvas *c = new TCanvas("c", "Event Cuts Plots", 1000, 700);
    c->Print((pdfFilePath + "[").c_str());

    auto printPage = [&](TH1* h, const char* opt, bool logY, bool logZ, std::function<void()> drawExtra) {
        c->Clear();
        gPad->SetLogy(logY ? 1 : 0);
        gPad->SetLogz(logZ ? 1 : 0);
        gPad->SetGrid(0, 0);
        
        gPad->SetBottomMargin(0.12);
        gPad->SetLeftMargin(0.12);
        gPad->SetRightMargin(logZ ? 0.15 : 0.05);

        if (std::string(opt).find("TEXT") != std::string::npos) {
            gPad->SetGridy(1);
            gPad->SetBottomMargin(0.25);
        }

        if (std::string(opt) == "HIST") {
            h->SetLineColor(kBlue+2);
            h->SetLineWidth(2);
        }
        
        h->Draw(opt);
        c->Update(); 
        
        if (drawExtra) drawExtra();
        
        c->Print(pdfFilePath.c_str());

        if (fOut) {
            fOut->cd();
            c->Write((std::string(h->GetName()) + "_canvas").c_str());
        }
    };

    // ==========================================
    // STEP 1: CUT 1 - VERTEX Z
    // ==========================================
    printPage(hists.h_Vz_all, "HIST", true, false, nullptr);
    printPage(hists.h_Vz_cut, "HIST", true, false, nullptr);
    printPage(hists.h2_TracksInFit_vs_PSD_all, "COLZ", false, true, nullptr);
    printPage(hists.h2_TracksInFit_vs_PSD_after_Vz, "COLZ", false, true, nullptr);

    // ==========================================
    // STEP 2: CUT 2 - PSD PERIPHERAL ENERGY
    // ==========================================
    printPage(hists.h_PSDperipheral_all, "HIST", true, false, nullptr);
    printPage(hists.h_PSDperipheral_cut, "HIST", true, false, nullptr);
    printPage(hists.h2_PSD_Peripheral_vs_Selected, "COLZ", false, true, nullptr);
    printPage(hists.h2_PSD_Peripheral_vs_Selected_cut, "COLZ", false, true, nullptr);
    printPage(hists.h2_TracksInFit_vs_PSD_after_PSD, "COLZ", false, true, nullptr);

    // ==========================================
    // STEP 3: CUT 3 - TRACKS RATIO
    // ==========================================
    printPage(hists.h2_TracksRatio_all, "COLZ", false, true, [](){ DrawTracksRatioLines(); });
    printPage(hists.h2_TracksRatio_cut, "COLZ", false, true, [](){ DrawTracksRatioLines(); });

    // ==========================================
    // STEP 4: CENTRALITY SELECTION
    // ==========================================
    int total_evts = hists.h_PSD_T2->GetEntries();
    double yMax1D = hists.h_PSD_T2->GetMaximum();
    hists.h_PSD_T2->SetFillColor(kAzure+1);
    printPage(hists.h_PSD_T2, "HIST", false, false, [&](){ DrawCentralityLines1D(yMax1D, cent_limits, cent_events, total_evts); });
    
    double yMax2D = hists.h2_TracksInFit_vs_PSD_cut->GetYaxis()->GetXmax();
    printPage(hists.h2_TracksInFit_vs_PSD_cut, "COLZ", false, true, [&](){ DrawCentralityLines2D(yMax2D, cent_limits, cent_events, total_evts); });
    printPage(hists.h2_TracksInFit_vs_PSD_Central, "COLZ", false, true, nullptr);

    // ==========================================
    // STEP 5: EVENT STATISTICS SUMMARY
    // ==========================================
    hists.hist_events->SetFillColor(kBlue-7);
    hists.hist_events->GetXaxis()->SetLabelSize(0.045);
    printPage(hists.hist_events, "HIST TEXT0", false, false, nullptr);

    c->Print((pdfFilePath + "]").c_str());

    delete c;

    fOut->cd();
    hists.hist_events->Write(); 
    hists.h_Vz_all->Write(); 
    hists.h_Vz_cut->Write();
    hists.h2_Tracks_vs_PSD_all->Write(); 
    hists.h2_Tracks_vs_PSD_cut->Write();
    hists.h2_PSD_Peripheral_vs_Selected->Write(); 
    hists.h2_PSD_Peripheral_vs_Selected_cut->Write();
    hists.h_PSDperipheral_all->Write(); 
    hists.h_PSDperipheral_cut->Write();
    hists.h2_TracksRatio_all->Write(); 
    hists.h2_TracksRatio_cut->Write();
    hists.h2_TracksInFit_vs_PSD_all->Write(); 
    hists.h2_TracksInFit_vs_PSD_after_Vz->Write(); 
    hists.h2_TracksInFit_vs_PSD_after_PSD->Write();
    hists.h2_TracksInFit_vs_PSD_cut->Write();
    hists.h2_TracksInFit_vs_PSD_Central->Write();
    hists.h_PSD_T2->Write();
    fOut->Close();

    std::cout << "=> Successfully saved event cut plots into: " << pdfFilePath << " and " << rootOutputFilePath << std::endl;
}