#include <iostream>
#include <fstream>
#include <string>
#include <cmath>

// Standard ROOT headers
#include "TCanvas.h"
#include "TGraphErrors.h"
#include "TAxis.h"
#include "TLegend.h"
#include "TStyle.h"
#include "TFile.h"
#include "TH1.h"

// Zintegrowana funkcja rysująca (nie potrzebuje Plotter.h ani EventXeLaMag.h!)
void DrawDeltaF2(const std::string& dataFile, const std::string& mixedFile, const std::string& outputBaseName) {
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

    // --- STRONA 1: Porównanie ---
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

    // --- STRONA 2: Różnica (Delta) ---
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

    // Zapis do ROOT
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

int main(int argc, char** argv) {
    if (argc < 4) {
        std::cout << "Usage: ./draw_delta <data_F2.txt> <mixed_F2.txt> <output_basename>" << std::endl;
        return 1;
    }

    std::string dataFile = argv[1];
    std::string mixedFile = argv[2];
    std::string outBase = argv[3];

    std::cout << "Calculating Delta F2..." << std::endl;
    std::cout << "Data:  " << dataFile << std::endl;
    std::cout << "Mixed: " << mixedFile << std::endl;

    DrawDeltaF2(dataFile, mixedFile, outBase);

    return 0;
}