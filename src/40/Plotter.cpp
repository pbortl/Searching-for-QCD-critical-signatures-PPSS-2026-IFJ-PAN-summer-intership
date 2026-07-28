#include "Plotter.h"
#include "Config.h"
#include "BetheBloch/src/include/BetheBlochWrapper.h"

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
#include "TLegend.h"
#include "TFile.h"
#include "TGaxis.h"
#include "TStyle.h"
#include "TPaveText.h"
#include "TText.h"
#include "TGraphErrors.h"
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
        for(size_t i = 0; i < cent_limits.size(); i++) { // ZMIANA: dynamiczny rozmiar
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
        pt.AddText(Form("0-5%%: %d", cent_events[0]));
        pt.AddText(Form("5-10%%: %d", cent_events[1] - cent_events[0]));
        pt.AddText(Form("10-15%%: %d", cent_events[2] - cent_events[1]));
        pt.AddText(Form("15-20%%: %d", cent_events[3] - cent_events[2]));
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
        pt.AddText(Form("0-5%: %d", cent_events[0]));
        pt.AddText(Form("5-10%: %d", cent_events[1] - cent_events[0]));
        pt.AddText(Form("10-15%: %d", cent_events[2] - cent_events[1]));
        pt.AddText(Form("15-20%: %d", cent_events[3] - cent_events[2]));
        pt.DrawClone("SAME");
    }
}

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

    // VZ LOG SCALE REQUESTED BY NIKOS:
    printPage(hists.h_Vz_all, "HIST", true, false, nullptr, fOut);
    printPage(hists.h_Vz_cut, "HIST", true, false, nullptr, fOut);

    printPage(hists.h_PSDperipheral_all, "HIST", true, false, nullptr, fOut);
    printPage(hists.h_PSDperipheral_cut, "HIST", true, false, nullptr, fOut);

    printPage(hists.h2_PSD_Peripheral_vs_Selected, "COLZ", false, true, nullptr, fOut);
    printPage(hists.h2_PSD_Peripheral_vs_Selected_cut, "COLZ", false, true, nullptr, fOut);

    printPage(hists.h2_TracksRatio_all, "COLZ", false, true, [](){ DrawTracksRatioLines(); }, fOut);
    printPage(hists.h2_TracksRatio_cut, "COLZ", false, true, [](){ DrawTracksRatioLines(); }, fOut);

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
    
    int total_evts = hists.h_PSD_T2->GetEntries();
    double yMax1D = hists.h_PSD_T2->GetMaximum();
    hists.h_PSD_T2->SetFillColor(kAzure+1);
    printPage(hists.h_PSD_T2, "HIST", false, false, [&](){ DrawCentralityLines1D(yMax1D, cent_limits, cent_events, total_evts); }, fOut);
    
    double yMax2D = hists.h2_TracksInFit_vs_PSD_cut->GetYaxis()->GetXmax();
    printPage(hists.h2_TracksInFit_vs_PSD_cut, "COLZ", false, true, [&](){ DrawCentralityLines2D(yMax2D, cent_limits, cent_events, total_evts); }, fOut);
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

  //track cut
    
    // 1. After Quality Cuts
    printPage(hists.h_dedx_ptot_pos_qual, "COLZ", false, true, drawBBCurves, fOut);
    printPage(hists.h_dedx_ptot_neg_qual, "COLZ", false, true, drawBBCurves, fOut);
    printPage(hists.h2_px_py_pos_qual, "COLZ", false, false, nullptr, fOut);
    printPage(hists.h2_px_py_neg_qual, "COLZ", false, false, nullptr, fOut);

    // 2. After Momentum Cuts
    printPage(hists.h_dedx_ptot_pos_mom, "COLZ", false, true, drawBBCurves, fOut);
    printPage(hists.h_dedx_ptot_neg_mom, "COLZ", false, true, drawBBCurves, fOut);
    printPage(hists.h2_px_py_pos_mom, "COLZ", false, false, nullptr, fOut);
    printPage(hists.h2_px_py_neg_mom, "COLZ", false, false, nullptr, fOut);

    // 3. After Proton ID
    printPage(hists.h_dedx_ptot_protons_id, "COLZ", false, true, drawBBCurves, fOut);
    printPage(hists.h2_px_py_protons_id, "COLZ", false, false, nullptr, fOut);
    hists.h_Y_CM_protons_id->SetFillColor(kBlue-7);
    printPage(hists.h_Y_CM_protons_id, "HIST", false, false, nullptr, fOut);

    // 4. After Rapidity (Final output for F2(M))
    TFile *fOutProtons = new TFile(rootProtonsFilePath.c_str(), "RECREATE");
    printPage(hists.h_dedx_ptot_protons_rap, "COLZ", false, true, drawBBCurves, fOutProtons);
    printPage(hists.h2_px_py_protons_rap, "COLZ", false, false, nullptr, fOutProtons);
    hists.h_Y_CM_protons_rap->SetFillColor(kBlue-7);
    printPage(hists.h_Y_CM_protons_rap, "HIST", false, false, nullptr, fOutProtons);

   

    std::string f2TxtPath = baseFileName + "_F2_results.txt";
    std::ifstream f2FileIn(f2TxtPath);
    TGraphErrors* gr_F2 = nullptr;
    //creating F2 file
    if (f2FileIn.is_open()) {
        std::string header;
        std::getline(f2FileIn, header);

        gr_F2 = new TGraphErrors();
        gr_F2->SetName("g_F2_M");
        gr_F2->SetTitle(("Raw F_{2}(M)" + suffix + "; M^{2} (Number of cells); F_{2}(M)").c_str());

        int M;
        double F2_M, err_F2_M, avg_N2, err_N2, avg_mul, err_mul;
        int pointIdx = 0;

        while (f2FileIn >> M >> F2_M >> err_F2_M >> avg_N2 >> err_N2 >> avg_mul >> err_mul) {
            double M2 = M * M;
            gr_F2->SetPoint(pointIdx, M2, F2_M);
            gr_F2->SetPointError(pointIdx, 0.0, err_F2_M);
            pointIdx++;
        }
        f2FileIn.close();

        if (pointIdx > 0) {
            c->Clear();
            gPad->SetLogy(0);
            gPad->SetLogz(0);
            gPad->SetGrid(1, 1);
            gPad->SetBottomMargin(0.12);
            gPad->SetLeftMargin(0.12);
            
            c->SetLogx(1);

            gr_F2->SetMarkerStyle(20);
            gr_F2->SetMarkerSize(0.6); 
            gr_F2->SetMarkerColor(kBlue+2);
            gr_F2->SetLineColor(kBlue+2);
            gr_F2->SetLineWidth(1); 
            
            gr_F2->GetHistogram()->GetXaxis()->SetLimits(1.0, 25000.0);

            TAxis* xAxis = gr_F2->GetXaxis();
            xAxis->SetMoreLogLabels(kFALSE);
            xAxis->SetNoExponent(kTRUE);
            xAxis->SetMaxDigits(4);
            xAxis->SetLabelSize(0.035);
            
            gStyle->SetEndErrorSize(3);
            
            gr_F2->Draw("AP"); 
            c->Update();
            c->Print(pdfFilePath.c_str());
            
            fOut->cd();
            gr_F2->Write("F2_M_Graph");
        }
    }

    c->Print((pdfFilePath + "]").c_str());
    
    delete bbElectron; delete bbPion; delete bbKaon; delete bbProton; delete bbDeuteron;
    delete bbLegend;
    delete c;
    if(gr_F2) delete gr_F2;

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
    fOut->Close();

    fOutProtons->Close();

    std::cout << "Successfully saved all plots into: " << pdfFilePath << std::endl;
}

void Plotter::DrawDeltaF2(const std::string& dataFile, const std::string& mixedFile, const std::string& outputBaseName) {
    gStyle->SetOptTitle(1);
    gStyle->SetPadGridX(kTRUE);
    gStyle->SetPadGridY(kTRUE);
    gStyle->SetEndErrorSize(3);

    std::ifstream inData(dataFile);
    std::ifstream inMixed(mixedFile);

    if (!inData.is_open() || !inMixed.is_open()) {
        std::cerr << "Error: Cannot open one of the TXT files!" << std::endl;
        return;
    }

    std::string header;
    std::getline(inData, header);
    std::getline(inMixed, header);

    TGraphErrors* grData = new TGraphErrors();
    grData->SetName("g_Data_F2");
    
    TGraphErrors* grMixed = new TGraphErrors();
    grMixed->SetName("g_Mixed_F2");

    TGraphErrors* grDelta = new TGraphErrors();
    grDelta->SetName("g_Delta_F2");

    int M_d, M_m;
    double F2_d, errF2_d, avgN2_d, errN2_d, avgMul_d, errMul_d;
    double F2_m, errF2_m, avgN2_m, errN2_m, avgMul_m, errMul_m;
    int pt = 0;

    std::cout << "\nReading data for comparison and delta plot..." << std::endl;

    while ((inData >> M_d >> F2_d >> errF2_d >> avgN2_d >> errN2_d >> avgMul_d >> errMul_d) &&
           (inMixed >> M_m >> F2_m >> errF2_m >> avgN2_m >> errN2_m >> avgMul_m >> errMul_m)) {
        
        if (M_d != M_m) {
            std::cerr << "Error: M step mismatch between files!" << std::endl;
            break;
        }

        double M2 = M_d * M_d;

        grData->SetPoint(pt, M2, F2_d);
        grData->SetPointError(pt, 0.0, errF2_d);

        grMixed->SetPoint(pt, M2, F2_m);
        grMixed->SetPointError(pt, 0.0, errF2_m);

        double delta_val = F2_d - F2_m;
        double delta_err = std::sqrt(errF2_d * errF2_d + errF2_m * errF2_m);
        grDelta->SetPoint(pt, M2, delta_val);
        grDelta->SetPointError(pt, 0.0, delta_err);

        pt++;
    }

    inData.close();
    inMixed.close();

    TCanvas* c1 = new TCanvas("c1", "F2 Comparison and Delta", 1000, 700);
    std::string pdfOut = outputBaseName + "_F2_Comparison.pdf";
    
    c1->Print((pdfOut + "[").c_str());

    c1->Clear();
    gPad->SetBottomMargin(0.12);
    gPad->SetLeftMargin(0.12);
    gPad->SetLogx(1);
    gPad->SetGridx(1);
    gPad->SetGridy(1);

    grData->SetTitle("F_{2}(M) Comparison; M^{2} (Number of cells); F_{2}(M)");
    grData->SetMarkerStyle(20);
    grData->SetMarkerSize(0.6); 
    grData->SetMarkerColor(kBlue+2);
    grData->SetLineColor(kBlue+2);
    grData->SetLineWidth(1); 
    grData->GetHistogram()->GetXaxis()->SetLimits(1.0, 25000.0);

    TAxis* xAxisTop = grData->GetXaxis();
    xAxisTop->SetMoreLogLabels(kFALSE);
    xAxisTop->SetNoExponent(kTRUE);
    xAxisTop->SetMaxDigits(4);
    xAxisTop->SetLabelSize(0.04); 
    xAxisTop->SetTitleSize(0.045);
    
    TAxis* yAxisTop = grData->GetYaxis();
    yAxisTop->SetLabelSize(0.04);
    yAxisTop->SetTitleSize(0.045);

    grMixed->SetMarkerStyle(20);
    grMixed->SetMarkerSize(0.6); 
    grMixed->SetMarkerColor(kRed+1);
    grMixed->SetLineColor(kRed+1);
    grMixed->SetLineWidth(1); 

    grData->Draw("AP");
    grMixed->Draw("P SAME");

    TLegend* leg = new TLegend(0.65, 0.15, 0.85, 0.30);
    leg->AddEntry(grData, "Original Events", "lep");
    leg->AddEntry(grMixed, "Mixed Events", "lep");
    leg->SetBorderSize(1);
    leg->SetFillColor(kWhite);
    leg->Draw("SAME");

    c1->Print(pdfOut.c_str()); 

    c1->Clear();
    gPad->SetBottomMargin(0.12);
    gPad->SetLeftMargin(0.12);
    gPad->SetLogx(1);
    gPad->SetGridx(1);
    gPad->SetGridy(1);

    grDelta->SetTitle("Difference #Delta F_{2}(M) (Original - Mixed); M^{2} (Number of cells); #Delta F_{2}(M)");
    grDelta->SetMarkerStyle(20);
    grDelta->SetMarkerSize(0.6);
    grDelta->SetMarkerColor(kBlack);
    grDelta->SetLineColor(kBlack);
    grDelta->SetLineWidth(1);
    grDelta->GetHistogram()->GetXaxis()->SetLimits(1.0, 25000.0);

    TAxis* xAxisBot = grDelta->GetXaxis();
    xAxisBot->SetMoreLogLabels(kFALSE);
    xAxisBot->SetNoExponent(kTRUE);
    xAxisBot->SetMaxDigits(4);
    xAxisBot->SetLabelSize(0.04);
    xAxisBot->SetTitleSize(0.045);

    TAxis* yAxisBot = grDelta->GetYaxis();
    yAxisBot->SetLabelSize(0.04);
    yAxisBot->SetTitleSize(0.045);
    yAxisBot->SetTitleOffset(1.2);

    grDelta->Draw("AP");

    c1->Print(pdfOut.c_str()); 

    c1->Print((pdfOut + "]").c_str());

    std::string rootOut = outputBaseName + "_F2_Comparison.root";
    TFile* fOut = new TFile(rootOut.c_str(), "RECREATE");
    grData->Write("Data_F2_Graph");
    grMixed->Write("Mixed_F2_Graph");
    grDelta->Write("Delta_F2_Graph");
    
    TCanvas* cData = new TCanvas("cData", "Comparison");
    cData->cd();
    gPad->SetLogx(1);
    gPad->SetGridx(1);
    gPad->SetGridy(1);
    grData->Draw("AP");
    grMixed->Draw("P SAME");
    leg->Draw("SAME");
    cData->Write("c_F2_Comparison_Canvas");
    
    TCanvas* cDelta = new TCanvas("cDelta", "Delta");
    cDelta->cd();
    gPad->SetLogx(1);
    gPad->SetGridx(1);
    gPad->SetGridy(1);
    grDelta->Draw("AP");
    cDelta->Write("c_F2_Delta_Canvas");
    
    fOut->Close();

    std::cout << "=> Successfully saved F2 comparison and delta plots to " << pdfOut << " (2 pages) and " << rootOut << std::endl;
    
    delete leg;
    delete cData;
    delete cDelta;
    delete c1;
    delete grData;
    delete grMixed;
    delete grDelta;
}