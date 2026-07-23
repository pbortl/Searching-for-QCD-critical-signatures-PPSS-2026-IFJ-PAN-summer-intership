#include "Plotter.h"
#include "Config.h"
#include "BetheBloch/src/include/BetheBlochWrapper.h"

#include <iostream>
#include <fstream>
#include <cmath>
#include <functional>
#include <vector>
#include <string>

#include "TCanvas.h"
#include "TF1.h"
#include "TLine.h"
#include "TLegend.h"
#include "TFile.h"
#include "TGaxis.h"
#include "TStyle.h"
#include "TPaveText.h"
#include "TText.h"

void Plotter::DrawAndSaveAll(HistogramManager& hists, const std::string& baseFileName, const std::vector<double>& cent_limits, const std::vector<int>& cent_events) {
    std::cout << "Generating multi-page PDF and saving to ROOT files..." << std::endl;

    gStyle->SetOptStat(1111);

    std::string txtFilePath = baseFileName + ".txt";
    std::string pdfFilePath = baseFileName + "_histograms.pdf";
    std::string rootOutputFilePath = baseFileName + "_histograms.root";
    std::string rootProtonsFilePath = baseFileName + "_protons.root";

    TFile *fOut = new TFile(rootOutputFilePath.c_str(), "RECREATE");

    hists.hist_events->SetMinimum(0);
    hists.hist_tracks->SetMinimum(0);

    std::string suffix = " for Xe+La " + std::to_string(Config::BeamMomentum) + "A GeV/c";

    hists.h_Vz_all->SetTitle(("distrib vz before the cut" + suffix).c_str());
    hists.h_Vz_cut->SetTitle(("distrib vz cutted" + suffix).c_str());
    hists.h2_PSD_Peripheral_vs_Selected->SetTitle(("PSD energy peripheral modules before cut" + suffix).c_str());
    hists.h2_PSD_Peripheral_vs_Selected_cut->SetTitle(("peripheral energy after the cut" + suffix).c_str());
    hists.h2_TracksRatio_all->SetTitle(("ntracks in fit vs ntracks all before cut" + suffix).c_str());
    hists.h2_TracksRatio_cut->SetTitle(("ntracks in fit vs ntracks after linear cut" + suffix).c_str());
    hists.h_clusters_VTPC_sum_all->SetTitle(("clusters VTPC 1 and 2 before the cut" + suffix).c_str());
    hists.h_clusters_VTPC_sum_cut->SetTitle(("clusters VTPC 1 and 2 cutted" + suffix).c_str());
    hists.h_clusters_All_all->SetTitle(("clusters all before cut" + suffix).c_str());
    hists.h_clusters_All_cut->SetTitle(("clusters after cut" + suffix).c_str());
    hists.h_clusters_dEdx_all->SetTitle(("clusters dE/dx before cut" + suffix).c_str());
    hists.h_clusters_dEdx_cut->SetTitle(("clusters dE/dx after cut" + suffix).c_str());
    hists.h_clusters_Ratio_all->SetTitle(("distrib ratio points/potential points before" + suffix).c_str());
    hists.h_clusters_Ratio_cut->SetTitle(("distrib ratio points/potential points after cut" + suffix).c_str());
    
    hists.h2_TracksInFit_vs_PSD_all->SetTitle(("PSD energy before any cuts" + suffix).c_str());
    hists.h2_TracksInFit_vs_PSD_after_Vz->SetTitle(("PSD energy after vertexz cut" + suffix).c_str());
    hists.h2_TracksInFit_vs_PSD_after_PSD->SetTitle(("PSD energy after vertexz cut and periph energy cut" + suffix).c_str());
    
    hists.h_PSD_T2->SetTitle(("events vs E_{PSD} (cut)" + suffix).c_str());
    hists.h2_TracksInFit_vs_PSD_cut->SetTitle(("track (fit) vs E_{PSD} (cut)" + suffix).c_str());
    hists.h2_TracksInFit_vs_PSD_Central->SetTitle(("fitted tracks vs PSD energy AFTER 0-20% cut" + suffix).c_str());
    
    hists.h2_bx_by_all->SetTitle(("bxby before" + suffix).c_str());
    hists.h2_bx_by_cut->SetTitle(("bxby cutted" + suffix).c_str());
    
    hists.hist_events->SetTitle(("nb of events after each event cut" + suffix).c_str());
    hists.hist_tracks->SetTitle(("(total) nb of tracks after each track cut" + suffix).c_str());

    TF1 *gaussFit = new TF1("gaussFit", "gaus", -606, -600);
    gaussFit->SetParameters(100000, Config::z_peak, 1.0);
    hists.h_Vz_all->Fit(gaussFit, "R0");

    std::ofstream txtOut(txtFilePath);
    if (txtOut.is_open()) {
        txtOut << "--- GAUSSIAN FIT RESULTS (VertexZ) ---\n";
        txtOut << "Amplitude: " << gaussFit->GetParameter(0) << "\n";
        txtOut << "Mean:      " << gaussFit->GetParameter(1) << " cm\n";
        txtOut << "Sigma:     " << gaussFit->GetParameter(2) << " cm\n";
        txtOut.close();
    }

    auto drawTracksRatioLines = [&]() {
        TLine *lineT1 = new TLine(Config::P1_x, Config::P1_y, Config::P2_x, Config::P2_y);
        lineT1->SetLineColor(kRed); lineT1->SetLineWidth(2); lineT1->Draw("SAME");
        TLine *lineT2 = new TLine(Config::P2_x, Config::P1_y, Config::P2_x, Config::P2_y);
        lineT2->SetLineColor(kRed); lineT2->SetLineWidth(2); lineT2->Draw("SAME");
        TLine *lineT3 = new TLine(Config::P1_x, Config::P1_y, Config::P2_x, Config::P1_y);
        lineT3->SetLineColor(kRed); lineT3->SetLineWidth(2); lineT3->Draw("SAME");

        double slope2 = (Config::P4_y - Config::P3_y) / (Config::P4_x - Config::P3_x);
        double end_x2 = Config::P3_x + (350.0 - Config::P3_y) / slope2;
        double end_y2 = 350.0;

        if (end_x2 > 1500.0) {
            end_x2 = 1500.0;
            end_y2 = Config::P3_y + slope2 * (1500.0 - Config::P3_x);
        }
        TLine *lineVector = new TLine(Config::P3_x, Config::P3_y, end_x2, end_y2);
        lineVector->SetLineColor(kRed); lineVector->SetLineWidth(2); lineVector->Draw("SAME");
    };

    auto drawCentralityLines1D = [&]() {
        double yMax = hists.h_PSD_T2->GetMaximum();
        for(int i=0; i<4; i++) {
            if(cent_limits[i] > 0) {
                TLine *line = new TLine(cent_limits[i], 0, cent_limits[i], yMax);
                line->SetLineColor(kRed);
                line->SetLineStyle(2); 
                line->SetLineWidth(2);
                line->Draw("SAME");
            }
        }
        
        TPaveText *pt = new TPaveText(0.65, 0.35, 0.85, 0.65, "NDC");
        pt->SetFillColor(kWhite);
        pt->SetBorderSize(1);
        pt->SetTextAlign(12);
        
        int total_evts = hists.h_PSD_T2->GetEntries();
        pt->AddText(Form("total evets: %d", total_evts));
        
        int n0_5 = cent_events[0];
        int n5_10 = cent_events[1] - cent_events[0];
        int n10_15 = cent_events[2] - cent_events[1];
        int n15_20 = cent_events[3] - cent_events[2];
        
        pt->AddText(Form("0-5%%: %d", n0_5));
        pt->AddText(Form("5-10%%: %d", n5_10));
        pt->AddText(Form("10-15%%: %d", n10_15));
        pt->AddText(Form("15-20%%: %d", n15_20));
        pt->Draw("SAME");
    };

    auto drawCentralityLines2D = [&]() {
        double yMax = hists.h2_TracksInFit_vs_PSD_cut->GetYaxis()->GetXmax();
        for(int i=0; i<4; i++) {
            if(cent_limits[i] > 0) {
                TLine *line = new TLine(cent_limits[i], 0, cent_limits[i], yMax);
                line->SetLineColor(kRed);
                line->SetLineStyle(2); 
                line->SetLineWidth(2);
                line->Draw("SAME");
            }
        }

        TPaveText *pt = new TPaveText(0.65, 0.65, 0.88, 0.88, "NDC");
        pt->SetFillColor(kWhite);
        pt->SetBorderSize(1);
        pt->SetTextAlign(12);
        
        int total_evts = hists.h_PSD_T2->GetEntries();
        pt->AddText(Form("total evets: %d", total_evts));
        
        int n0_5 = cent_events[0];
        int n5_10 = cent_events[1] - cent_events[0];
        int n10_15 = cent_events[2] - cent_events[1];
        int n15_20 = cent_events[3] - cent_events[2];
        
        pt->AddText(Form("0-5%%: %d", n0_5));
        pt->AddText(Form("5-10%%: %d", n5_10));
        pt->AddText(Form("10-15%%: %d", n10_15));
        pt->AddText(Form("15-20%%: %d", n15_20));
        pt->Draw("SAME");
    };

    TCanvas *c = new TCanvas("c", "Plots", 1000, 700);
    c->Print((pdfFilePath + "[").c_str());

    auto printPage = [&](TH1* h, const char* opt, bool logY, bool logZ, std::function<void()> drawExtra, TFile* targetFile) {
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

        if (targetFile) {
            targetFile->cd();
            c->Write((std::string(h->GetName()) + "_canvas").c_str());
        }
    };

    printPage(hists.h_Vz_all, "HIST", false, false, nullptr, fOut);
    printPage(hists.h_Vz_cut, "HIST", false, false, nullptr, fOut);

    printPage(hists.h_PSDperipheral_all, "HIST", true, false, nullptr, fOut);
    printPage(hists.h_PSDperipheral_cut, "HIST", true, false, nullptr, fOut);

    printPage(hists.h2_PSD_Peripheral_vs_Selected, "COLZ", false, true, nullptr, fOut);
    printPage(hists.h2_PSD_Peripheral_vs_Selected_cut, "COLZ", false, true, nullptr, fOut);

    printPage(hists.h2_TracksRatio_all, "COLZ", false, true, drawTracksRatioLines, fOut);
    printPage(hists.h2_TracksRatio_cut, "COLZ", false, true, drawTracksRatioLines, fOut);

    printPage(hists.h_clusters_VTPC_sum_all, "HIST", true, false, nullptr, fOut);
    printPage(hists.h_clusters_VTPC_sum_cut, "HIST", true, false, nullptr, fOut);
    
    printPage(hists.h_clusters_All_all, "HIST", true, false, nullptr, fOut);
    printPage(hists.h_clusters_All_cut, "HIST", true, false, nullptr, fOut);
    
    printPage(hists.h_clusters_dEdx_all, "HIST", true, false, nullptr, fOut);
    printPage(hists.h_clusters_dEdx_cut, "HIST", true, false, nullptr, fOut);
    
    printPage(hists.h_clusters_Ratio_all, "HIST", true, false, nullptr, fOut);
    printPage(hists.h_clusters_Ratio_cut, "HIST", true, false, nullptr, fOut);

    printPage(hists.h2_TracksInFit_vs_PSD_all, "COLZ", false, true, nullptr, fOut);
    printPage(hists.h2_TracksInFit_vs_PSD_after_Vz, "COLZ", false, true, nullptr, fOut);
    printPage(hists.h2_TracksInFit_vs_PSD_after_PSD, "COLZ", false, true, nullptr, fOut);
    
    hists.h_PSD_T2->SetFillColor(kAzure+1);
    printPage(hists.h_PSD_T2, "HIST", false, false, drawCentralityLines1D, fOut);
    
    printPage(hists.h2_TracksInFit_vs_PSD_cut, "COLZ", false, true, drawCentralityLines2D, fOut);
    printPage(hists.h2_TracksInFit_vs_PSD_Central, "COLZ", false, true, nullptr, fOut);

    printPage(hists.h2_bx_by_all, "COLZ", false, true, nullptr, fOut);
    printPage(hists.h2_bx_by_cut, "COLZ", false, true, nullptr, fOut);

    hists.hist_events->SetFillColor(kBlue-7);
    hists.hist_events->GetXaxis()->SetLabelSize(0.045);
    printPage(hists.hist_events, "HIST TEXT0", false, false, nullptr, fOut);

    hists.hist_tracks->SetFillColor(kGreen-7);
    hists.hist_tracks->GetXaxis()->SetLabelSize(0.045);
    printPage(hists.hist_tracks, "HIST TEXT0", false, false, nullptr, fOut);

    Dedx::BetheBlochWrapper bbWrapper;
    double xMin = -0.5;
    double xMax = 2.2;

    TF1 *bbElectron = new TF1("bbElectron", [&](double *x, double *){ return bbWrapper(0, std::pow(10, x[0])); }, xMin, xMax, 0);
    bbElectron->SetLineColor(kMagenta);
    bbElectron->SetLineStyle(7);
    bbElectron->SetLineWidth(2);

    TF1 *bbPion = new TF1("bbPion", [&](double *x, double *){ return bbWrapper(1, std::pow(10, x[0])); }, xMin, xMax, 0);
    bbPion->SetLineColor(kBlue);
    bbPion->SetLineStyle(7);
    bbPion->SetLineWidth(2);

    TF1 *bbKaon = new TF1("bbKaon", [&](double *x, double *){ return bbWrapper(2, std::pow(10, x[0])); }, xMin, xMax, 0);
    bbKaon->SetLineColor(kGreen);
    bbKaon->SetLineStyle(7);
    bbKaon->SetLineWidth(2);

    TF1 *bbProton = new TF1("bbProton", [&](double *x, double *){ return bbWrapper(3, std::pow(10, x[0])); }, xMin, xMax, 0);
    bbProton->SetLineColor(kRed);
    bbProton->SetLineStyle(7);
    bbProton->SetLineWidth(2);

    TF1 *bbDeuteron = new TF1("bbDeuteron", [&](double *x, double *){ return bbWrapper(4, std::pow(10, x[0])); }, xMin, xMax, 0);
    bbDeuteron->SetLineColor(kOrange+1);
    bbDeuteron->SetLineStyle(7);
    bbDeuteron->SetLineWidth(2);

    TLegend *bbLegend = new TLegend(0.65, 0.15, 0.85, 0.40);
    bbLegend->AddEntry(bbElectron, "Electron", "l");
    bbLegend->AddEntry(bbPion, "Pion", "l");
    bbLegend->AddEntry(bbKaon, "Kaon", "l");
    bbLegend->AddEntry(bbProton, "Proton", "l");
    bbLegend->AddEntry(bbDeuteron, "Deuteron", "l");
    bbLegend->SetBorderSize(1);
    bbLegend->SetFillColor(kWhite);

    auto drawBBCurves = [&]() {
        bbElectron->Draw("SAME");
        bbPion->Draw("SAME");
        bbKaon->Draw("SAME");
        bbProton->Draw("SAME");
        bbDeuteron->Draw("SAME");
        bbLegend->Draw("SAME");
    };

    printPage(hists.h_dedx_ptot_pos, "COLZ", false, true, drawBBCurves, fOut);
    printPage(hists.h_dedx_ptot_neg, "COLZ", false, true, drawBBCurves, fOut);
    
    TFile *fOutProtons = new TFile(rootProtonsFilePath.c_str(), "RECREATE");

    hists.h_dedx_ptot_protons->SetTitle(("dE/dx vs log(p_{tot}) selected protons positive tracks " + suffix + "; Log_{10}(p_{tot} / [GeV/c]); dE/dx (MIP)").c_str());
    hists.h_dedx_ptot_neg_protons->SetTitle(("dE/dx vs log(p_{tot}) selected protons negative tracks " + suffix + "; Log_{10}(p_{tot} / [GeV/c]); dE/dx (MIP)").c_str());

    printPage(hists.h_dedx_ptot_protons, "COLZ", false, true, drawBBCurves, fOutProtons);
    printPage(hists.h_dedx_ptot_neg_protons, "COLZ", false, true, drawBBCurves, fOutProtons);

    c->Print((pdfFilePath + "]").c_str());
    
    delete bbElectron; delete bbPion; delete bbKaon; delete bbProton; delete bbDeuteron;
    delete bbLegend;
    delete c;

    fOut->cd();
    hists.hist_events->Write(); 
    hists.hist_tracks->Write();
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
    hists.h_clusters_dEdx_all->Write(); 
    hists.h_clusters_dEdx_cut->Write();
    hists.h_clusters_VTPC_sum_all->Write(); 
    hists.h_clusters_VTPC_sum_cut->Write();
    hists.h_clusters_PotAll_all->Write(); 
    hists.h_clusters_PotAll_cut->Write();
    hists.h_clusters_All_all->Write(); 
    hists.h_clusters_All_cut->Write();
    hists.h_clusters_Ratio_all->Write(); 
    hists.h_clusters_Ratio_cut->Write();
    hists.h2_bx_by_all->Write(); 
    hists.h2_bx_by_cut->Write();
    hists.h_dedx_ptot_pos->Write();
    hists.h_dedx_ptot_neg->Write();
    hists.h_dedx_ptot_protons->Write();
    fOut->Close();

    fOutProtons->cd();
    hists.h_dedx_ptot_protons->Write();
    hists.h_dedx_ptot_neg_protons->Write();
    fOutProtons->Close();

    std::cout << "Successfully saved all plots into: " << pdfFilePath << std::endl;
}