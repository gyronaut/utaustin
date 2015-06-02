#include "TLegend.h"
#include <string>
#include <sstream>

void plot_K_phi(string inputName){
    gROOT->Reset();
    gROOT->SetStyle("Plain");
    gStyle->SetPalette(1);
    gStyle->SetCanvasColor(kWhite);
    gStyle->SetCanvasBorderMode(0);
    gStyle->SetPadBorderMode(0);
    gStyle->SetTitleBorderSize(0);
    gStyle->SetOptStat(0);
    gStyle->SetErrorX(0);

    TFile *histoFile = new TFile(inputName.c_str());
    TH1F *hadronPt = (TH1F *)histoFile->Get("HadronPt");
    hadronPt->SetTitle("p_{T}^{h} Distribution");
    TH1F *phiPt = (TH1F *)histoFile->Get("PhiPt");
    phiPt->SetTitle("p_{T}^{#phi} Distribution");
    THnSparseF *DphiHPhi = (THnSparseF *)histoFile->Get("DphiHPhi"); 
    THnSparseF *DphiHK0 = (THnSparseF *)histoFile->Get("DphiHK0");
    THnSparseF *DphiPiPhi = (THnSparseF *)histoFile->Get("DphiPiPhi");
    THnSparseF *DphiKPhi = (THnSparseF *)histoFile->Get("DphiKPhi");
    THnSparseF *DphiPiK0 = (THnSparseF *)histoFile->Get("DphiPiK0");

    TH1F *hadronPtPEREVENT = (TH1F*) hadronPt->Clone("HadronpTperEvent");
    TH1F *phiPtPEREVENT = (TH1F*) phiPt->Clone("PhipTperEvent");
    //set-up "per event" histograms
    for(int i = 0; i < hadronPtPEREVENT->GetNbinsX(); i++){
	    hadronPtPEREVENT->SetBinContent(i, (float)hadronPtPEREVENT->GetBinContent(i)/hadronPtPEREVENT->GetEntries());
	    phiPtPEREVENT->SetBinContent(i, (float)phiPtPEREVENT->GetBinContent(i)/phiPtPEREVENT->GetEntries());
    }

    fprintf(stdout, "About to project sparse! \n");

    TH3D* HPhiDphi = DphiHPhi->Projection(0, 1, 2);
    TH3D* HK0Dphi = DphiHK0->Projection(0, 1, 2);
    TH3D* PiPhiDphi = DphiPiPhi->Projection(0, 1, 2);
    TH3D* KPhiDphi = DphiKPhi->Projection(0,1,2);
    TH3D* PiKDphi = DphiPiK0->Projection(0,1,2);

    fprintf(stdout, "done projecting sparse! \n");

    TH1D *fHPhiDphi[3];
    fHPhiDphi[0] = HPhiDphi->ProjectionZ("ptH4_1ptPhi2", 80, 1000, 20, 40);
    fHPhiDphi[0]->SetTitle("p_{T}^{h} > 4,   1 < p_{T}^{#phi} < 2 GeV/c");
    fHPhiDphi[1] = HPhiDphi->ProjectionZ("ptH4_2ptPhi4", 80, 1000, 41, 80);
    fHPhiDphi[1]->SetTitle("p_{T}^{h} > 4,   2 < p_{T}^{#phi} < 4 GeV/c");
    fHPhiDphi[2] = HPhiDphi->ProjectionZ("ptH4_4ptPhi", 80, 1000, 81, 1000);
    fHPhiDphi[2]->SetTitle("p_{T}^{h} > 4,   p_{T}^{#phi} > 4 GeV/c");

    TH1D *fHK0Dphi[3];
    fHK0Dphi[0] = HK0Dphi->ProjectionZ("ptH4_1ptK02", 80, 1000, 20, 40);
    fHK0Dphi[0]->SetTitle("p_{T}^{h} > 4,   1 < p_{T}^{K} < 2 GeV/c");
    fHK0Dphi[1] = HK0Dphi->ProjectionZ("ptH4_2ptK04", 80, 1000, 41, 80);
    fHK0Dphi[1]->SetTitle("p_{T}^{h} > 4,   2 < p_{T}^{K} < 4 GeV/c");
    fHK0Dphi[2] = HK0Dphi->ProjectionZ("ptH4_4ptK0", 80, 1000, 81, 1000);
    fHK0Dphi[2]->SetTitle("p_{T}^{h} > 4,   p_{T}^{K} > 4 GeV/c");

    TH1D *fPiPhiDphi[3];
    fPiPhiDphi[0] = PiPhiDphi->ProjectionZ("ptPi4_1ptPhi2", 80, 1000, 20, 40);
    fPiPhiDphi[0]->SetTitle("p_{T}^{h} > 4,   1 < p_{T}^{#phi} < 2 GeV/c");
    fPiPhiDphi[1] = PiPhiDphi->ProjectionZ("ptPi4_2ptPhi4", 80, 1000, 41, 80);
    fPiPhiDphi[1]->SetTitle("p_{T}^{h} > 4,   2 < p_{T}^{#phi} < 4 GeV/c");
    fPiPhiDphi[2] = PiPhiDphi->ProjectionZ("ptPi4_4ptPhi", 80, 1000, 81, 1000);
    fPiPhiDphi[2]->SetTitle("p_{T}^{#pi} > 4,   p_{T}^{#phi} > 4 GeV/c");

    TH1D *fKPhiDphi[3];
    fKPhiDphi[0] = KPhiDphi->ProjectionZ("ptK4_1ptPhi2", 80, 1000, 20, 40);
    fKPhiDphi[0]->SetTitle("p_{T}^{K} > 4,   1 < p_{T}^{#phi} < 2 GeV/c");
    fKPhiDphi[1] = KPhiDphi->ProjectionZ("ptK4_2ptPhi4", 80, 1000, 41, 80);
    fKPhiDphi[1]->SetTitle("p_{T}^{K} > 4,   2 < p_{T}^{#phi} < 4 GeV/c");
    fKPhiDphi[2] = KPhiDphi->ProjectionZ("ptK4_4ptPhi", 80, 1000, 81, 1000);
    fKPhiDphi[2]->SetTitle("p_{T}^{K} > 4,   p_{T}^{#phi} > 4 GeV/c");

    TH1D *fPiKDphi[3];
    fPiKDphi[0] = PiKDphi->ProjectionZ("ptPi4_1ptK2", 80, 1000, 20, 40);
    fPiKDphi[0]->SetTitle("p_{T}^{#pi} > 4,   1 < p_{T}^{K} < 2 GeV/c");
    fPiKDphi[1] = PiKDphi->ProjectionZ("ptPi4_2ptK4", 80, 1000, 41, 80);
    fPiKDphi[1]->SetTitle("p_{T}^{#pi} > 4,   2 < p_{T}^{K} < 4 GeV/c");
    fPiKDphi[2] = PiKDphi->ProjectionZ("ptPi4_4ptK", 80, 1000, 81, 1000);
    fPiKDphi[2]->SetTitle("p_{T}^{#pi} > 4,   p_{T}^{K} > 4 GeV/c");

 
    TH1D* HPhiDphiInclusive = HPhiDphi->ProjectionZ("dphiptPhiInc",80,1000,40,1000);
    HPhiDphiInclusive->SetTitle("#Delta#phi for hadron-#phi(1020), p_{T}^{h} > 4 GeV/c, p_{T}^{#phi} > 2 GeV/c");

    TH1D* HK0DphiInclusive = HK0Dphi->ProjectionZ("dphiptKInc,80, 1000, 40, 1000");
    HK0DphiInclusive->SetTitle("#Delta#phi for hadron-K^{0}, p_{T}^{h} > 4 GeV/c, p_{T}^{K} > 2 GeV/c");

    Double_t xbins[19] = {-1.57079, -1.22173, -0.872662, -0.523596, -0.17453, 0.174536, \
                        0.523602, 0.87267, 1.22173, 1.5708, 1.91987, 2.26893, 2.618, 2.96706, \
                        3.31613, 3.66519, 4.01426, 4.36333, 4.71239}; 
    TH1D* ScaledHPhiDphiInclusive = HPhiDphiInclusive->Rebin(18, "ScaledDphiHPhiInclusive", xbins);
    TH1D* ScaledHK0DphiInclusive = HK0DphiInclusive->Rebin(18, "ScaledHK0DphiInclusive", xbins);

    //Setting up histogram for Christina's data:
    Double_t data[18] = {287.344, 263.394, 340.795, 444.768, 1162.75, 526.076, 285.806, \
                   230.348, 300.421, 271.289, 283.492, 354.686, 423.02, 497.954, 424.278, \
                   399.197, 312.045, 339.193};
    Double_t error[18] = {19.7231, 19.1311, 21.6102, 28.5307, 49.8197, 29.0345, 20.445, 19, \
                  20.025, 19.8997, 21.0713, 23.4094, 27.0924, 29.1204, 26.5707, 24.3516, \
                  21.166, 21.5639};
    TH1D* HPhiDphiData = ScaledHPhiDphiInclusive->Clone();
    for(int i = 0; i < 18; i++){
        HPhiDphiData->SetBinContent(i+1,data[i]);
        HPhiDphiData->SetBinError(i+1, error[i]);
    }
    //Checking data
    TCanvas *cData = new TCanvas("dphidata", "#Delta#phi data", 50, 50, 600, 400);
    cData->cd();
    HPhiDphiData->SetLineColor(1);
    HPhiDphiData->Draw("HIST");

    //Checking rebinned data
    TCanvas *cRebin = new TCanvas("dphirebinned", "#Delta#phi rebinned", 50, 50, 600, 400);
    cRebin->cd();
    ScaledHPhiDphiInclusive->SetLineColor(1);
    ScaledHPhiDphiInclusive->Draw("HIST");

    //Checking rebinned K0
    TCanvas *cRebinK0 = new TCanvas("dphirebinnedK0", "#Delta#phi K0 rebinned", 50, 50, 600, 400);
    cRebinK0->cd();
    ScaledHK0DphiInclusive->SetLineColor(1);
    ScaledHK0DphiInclusive->Draw("HIST");

/*
    //Setting up the scaling factor between the simulation/data:
    Double_t scale_factor = 0;
    for(int i = 1; i < 19; i++){
        HPhiDphiData->SetBinContent(i, HPhiDphiData->GetBinContent(i)-287.674);
        scale_factor += ((HPhiDphiData->GetBinContent(i))/(ScaledHPhiDphiInclusive->GetBinContent(i)+1));
    }
    scale_factor = scale_factor/18.0;
    ScaledHPhiDphiInclusive->Scale(scale_factor);

    //Graph of scaled data
    TCanvas *cScaled = new TCanvas("dphiscaled", "#Delta#phi scaled", 50, 50, 600, 400);
    cScaled->cd();
    ScaledHPhiDphiInclusive->SetLineColor(1);
    ScaledHPhiDphiInclusive->Draw("SAME E1");
    HPhiDphiData->SetLineColor(2);
    HPhiDphiData->Draw("SAME E1");
*/
    // Graph of dphi all pT
    TCanvas *cInc = new TCanvas("dphiInclusive", "#Delta#phi inclusive", 50, 50, 600, 400);
    cInc->cd();
    HPhiDphiInclusive->SetLineColor(1);
    HPhiDphiInclusive->SetFillColor(16);
    HPhiDphiInclusive->DrawCopy();
    HPhiDphiInclusive->Sumw2();
    HPhiDphiInclusive->SetMarkerStyle(21);
    HPhiDphiInclusive->SetMarkerSize(2);
    HPhiDphiInclusive->GetXaxis()->SetTitle("#Delta#phi");
    HPhiDphiInclusive->Draw("SAME E1");

    TCanvas *cInc = new TCanvas("dphiKInclusive", "#Delta#phi inclusive", 65, 65, 600, 400);
    cInc->cd();
    HK0DphiInclusive->SetLineColor(1);
    HK0DphiInclusive->SetFillColor(16);
    HK0DphiInclusive->DrawCopy();
    HK0DphiInclusive->Sumw2();
    HK0DphiInclusive->SetMarkerStyle(21);
    HK0DphiInclusive->SetMarkerSize(2);
    HK0DphiInclusive->GetXaxis()->SetTitle("#Delta#phi");
    HK0DphiInclusive->Draw("SAME E1");
    
    TF1* P = new TF1("P", "pol0", 1, 2);
    P->SetParameter(0, HPhiDphiInclusive->GetBinContent(1.5));
    HPhiDphiInclusive->Fit(P, "N", "", 1, 2);
    TF1* fit = new TF1("fit1", "pol0(0)+gaus(1)+gaus(4)", -1.5, 4.7);
    fit->SetParameter(0, P->GetParameter(0));
//    fit->SetParameter(1, 2000);
    fit->SetParameter(2, 0);
    fit->SetParameter(3, 1);
    fit->SetParameter(5, 3.1416);
    fit->SetParameter(6, 1.5);

    fit->SetLineColor(2);
    fit->SetLineWidth(3);
    HPhiDphiInclusive->Fit(fit, "", "SAME C");
 
    // Graph of Phi(1020) pT Distribution
    TCanvas *c2 = new TCanvas("PhiPtCanvas", "Pt distribution for Phi mesons", 60, 60, 600, 400);
    c2->cd();
    c2->SetLogy();
    phiPt->GetXaxis()->SetRangeUser(0, 15);
    phiPt->GetXaxis()->SetTitle("p_{T}^{#phi}");
    phiPt->SetLineColor(1);
    phiPt->SetFillColor(16);
    phiPt->Draw();
 
    Double_t numphiPt1_2 = phiPt->Integral(20, 40);
    Double_t numphiPt2_4 = phiPt->Integral(41, 80);
    Double_t numphiPtgt4 = phiPt->Integral(81, 1000);   
    
    printf("Number of phi1_2 = %d\nNumber of phi2_4 = %d\nNumber of phigt4 = %d\n", numphiPt1_2, numphiPt2_4, numphiPtgt4);
    
    TH1F* phiPt1_2 = (TH1F*) phiPt->Clone();
    phiPt1_2->GetXaxis()->SetRangeUser(1, 2);
    phiPt1_2->SetLineColor(1);
    phiPt1_2->SetFillColor(2);
    phiPt1_2->Draw("SAME");

    TH1F* phiPt2_4 = (TH1F*) phiPt->Clone();
    phiPt2_4->GetXaxis()->SetRangeUser(2, 4);
    phiPt2_4->SetLineColor(1);
    phiPt2_4->SetFillColor(3);
    phiPt2_4->Draw("SAME");

    TH1F* phiPtgt4 = (TH1F*) phiPt->Clone();
    phiPtgt4->GetXaxis()->SetRangeUser(4, 100);
    phiPtgt4->SetLineColor(1);
    phiPtgt4->SetFillColor(4);
    phiPtgt4->Draw("SAME");

    TLatex label;
    label.SetTextSize(0.04);

    Int_t maxHeight = phiPt->GetBinContent(phiPt->GetMaximumBin());

    label.DrawLatex(8, maxHeight*0.9, "Counts of #phi(1020) mesons:");

    string label1_2 =  "#color[2]{1 < p_{T}^{#phi} < 2 GeV/c} : ";
    stringstream ss;
    ss << numphiPt1_2;
    string number1_2 = ss.str();
    string str_full_label1_2 = label1_2 + number1_2;
    char* full_label1_2 = str_full_label1_2.c_str();
    label.DrawLatex(8.5, maxHeight*.7, full_label1_2);

    ss.str(std::string());
    string label2_4 =  "#color[3]{2 < p_{T}^{#phi} < 4 GeV/c} : ";
    ss << numphiPt2_4;
    string number2_4 = ss.str();
    string str_full_label2_4 = label2_4 + number2_4;
    char* full_label2_4 = str_full_label2_4.c_str();
    label.DrawLatex(8.5, maxHeight*.7*.7, full_label2_4);

    ss.str(std::string());
    string labelgt4 = "#color[4]{p_{T}^{#phi} > 4 GeV/c} : ";
    ss << numphiPtgt4;
    string numbergt4 = ss.str();
    string str_full_labelgt4 = labelgt4 + numbergt4;
    char* full_labelgt4 = str_full_labelgt4.c_str();
    label.DrawLatex(8.5, maxHeight*.7*.7*.7, full_labelgt4);

    // Graph of Phi(1020) pT Distribution for hadron pT > 4 GeV/c
    TCanvas *c3 = new TCanvas("PhiPtgt5Canvas", "Pt distribution for Phi mesons", 60, 60, 600, 400);
    c3->cd();
    TH1D* phiPt_Hgt5 = HPhiDphi->ProjectionY("phiPtgt5", 80, 1000, 0, 256);
    phiPt_Hgt5->SetTitle("p_{T}^{#phi} Distribution for events with a hadron with p_{T}^{h} > 5 GeV/c");
    phiPt_Hgt5->GetXaxis()->SetRangeUser(0, 15);
    phiPt_Hgt5->GetXaxis()->SetTitle("p_{T}^{#phi}");
    phiPt_Hgt5->SetLineColor(1);
    phiPt_Hgt5->SetFillColor(16);
    phiPt_Hgt5->Draw();
 
    Double_t numphiPt_Hgt5_1_2 = phiPt_Hgt5->Integral(20, 40);
    Double_t numphiPt_Hgt5_2_4 = phiPt_Hgt5->Integral(41, 80);
    Double_t numphiPt_Hgt5_gt4 = phiPt_Hgt5->Integral(81, 1000);   
   
    printf("Number of phi1_2 = %d\nNumber of phi2_4 = %d\nNumber of phigt4 = %d\n", numphiPt_Hgt5_1_2, numphiPt_Hgt5_2_4, numphiPt_Hgt5_gt4);

    // Graph of high-pt hadron (> 4 Gev/c) phi(1020) correlation
    TCanvas *c6 = new TCanvas("H-Phi correlation for pTh>4, ptPhi ranges", "H-Phi correlation for pTh>4, ptPhi ranges", 50, 50, 600, 400);
    c6->Divide(2,2);
    for(int i=-1; i<3; i++){
        c6->cd(i+2);
        if(i+2==1){
            TLatex heading;
            heading.SetTextAlign(22);
            heading.SetTextSize(0.1);
            heading.DrawLatex(0.5, 0.5, "#splitline{#Delta#phi correlation between}{#splitline{high p_{T} hadrons (> 4 GeV/c)}{and #phi(1020) mesons}}");
        }else{
            fHPhiDphi[i]->GetXaxis()->SetTitle("#Delta#phi");
            fHPhiDphi[i]->SetLineWidth(0);
            fHPhiDphi[i]->SetLineColor(0);
            fHPhiDphi[i]->SetFillColor(16);
            fHPhiDphi[i]->DrawCopy("HIST");
            fHPhiDphi[i]->Sumw2();
            fHPhiDphi[i]->SetMarkerStyle(21);
            fHPhiDphi[i]->SetMarkerSize(1);
            fHPhiDphi[i]->SetMarkerColor(1);
            fHPhiDphi[i]->SetLineColor(1);
            fHPhiDphi[i]->SetLineWidth(1);
            fHPhiDphi[i]->Draw("SAME E1");
        }
    }

    // Graph of high-pt hadron (> 4 Gev/c) K0 correlation
    TCanvas *c9 = new TCanvas("H-K0 correlation for pTh>4, ptK0 ranges", "H-K0 correlation for pTh>4, ptK0 ranges", 50, 50, 600, 400);
    c9->Divide(2,2);
    for(int i=-1; i<3; i++){
        c9->cd(i+2);
        if(i+2==1){
            TLatex heading;
            heading.SetTextAlign(22);
            heading.SetTextSize(0.1);
            heading.DrawLatex(0.5, 0.5, "#splitline{#Delta#phi correlation between}{#splitline{high p_{T} hadrons (> 4 GeV/c)}{and K^{0} mesons}}");
        }else{
            fHK0Dphi[i]->GetXaxis()->SetTitle("#Delta#phi");
            fHK0Dphi[i]->SetLineWidth(0);
            fHK0Dphi[i]->SetLineColor(0);
            fHK0Dphi[i]->SetFillColor(16);
            fHK0Dphi[i]->DrawCopy("HIST");
            fHK0Dphi[i]->Sumw2();
            fHK0Dphi[i]->SetMarkerStyle(21);
            fHK0Dphi[i]->SetMarkerSize(1);
            fHK0Dphi[i]->SetMarkerColor(1);
            fHK0Dphi[i]->SetLineColor(1);
            fHK0Dphi[i]->SetLineWidth(1);
            fHK0Dphi[i]->Draw("SAME E1");
        }
    }

    // Graph of high-pt Pion (> 4 Gev/c) Phi correlation
    TCanvas *c12 = new TCanvas("Pi-Phi correlation for pTh>4, ptPhi ranges", "Pi-Phi correlation for pTh>4, ptPhi ranges", 50, 50, 600, 400);
    c12->Divide(2,2);
    for(int i=-1; i<3; i++){
        c12->cd(i+2);
        if(i+2==1){
            TLatex heading;
            heading.SetTextAlign(22);
            heading.SetTextSize(0.1);
            heading.DrawLatex(0.5, 0.5, "#splitline{#Delta#phi correlation between}{#splitline{high p_{T} Pions (> 4 GeV/c)}{and #phi mesons}}");
        }else{
            fPiPhiDphi[i]->GetXaxis()->SetTitle("#Delta#phi");
            fPiPhiDphi[i]->SetLineWidth(0);
            fPiPhiDphi[i]->SetLineColor(0);
            fPiPhiDphi[i]->SetFillColor(16);
            fPiPhiDphi[i]->DrawCopy("HIST");
            fPiPhiDphi[i]->Sumw2();
            fPiPhiDphi[i]->SetMarkerStyle(21);
            fPiPhiDphi[i]->SetMarkerSize(1);
            fPiPhiDphi[i]->SetMarkerColor(1);
            fPiPhiDphi[i]->SetLineColor(1);
            fPiPhiDphi[i]->SetLineWidth(1);
            fPiPhiDphi[i]->Draw("SAME E1");
        }
    }

    // Graph of high-pt hadron (> 4 Gev/c) K0 correlation
    TCanvas *c13 = new TCanvas("K0-Phi correlation for pTh>4, ptPhi ranges", "K-Phi correlation for pTh>4, ptPhi ranges", 50, 50, 600, 400);
    c13->Divide(2,2);
    for(int i=-1; i<3; i++){
        c13->cd(i+2);
        if(i+2==1){
            TLatex heading;
            heading.SetTextAlign(22);
            heading.SetTextSize(0.1);
            heading.DrawLatex(0.5, 0.5, "#splitline{#Delta#phi correlation between}{#splitline{high p_{T} Kaons (> 4 GeV/c)}{and #phi mesons}}");
        }else{
            fKPhiDphi[i]->GetXaxis()->SetTitle("#Delta#phi");
            fKPhiDphi[i]->SetLineWidth(0);
            fKPhiDphi[i]->SetLineColor(0);
            fKPhiDphi[i]->SetFillColor(16);
            fKPhiDphi[i]->DrawCopy("HIST");
            fKPhiDphi[i]->Sumw2();
            fKPhiDphi[i]->SetMarkerStyle(21);
            fKPhiDphi[i]->SetMarkerSize(1);
            fKPhiDphi[i]->SetMarkerColor(1);
            fKPhiDphi[i]->SetLineColor(1);
            fKPhiDphi[i]->SetLineWidth(1);
            fKPhiDphi[i]->Draw("SAME E1");
        }
    }
 
    // Graph of high-pt hadron (> 4 Gev/c) K0 correlation
    TCanvas *c14 = new TCanvas("Pi-K0 correlation for pTPi>4, ptK0 ranges", "Pi-K0 correlation for pTPi>4, ptK0 ranges", 50, 50, 600, 400);
    c14->Divide(2,2);
    for(int i=-1; i<3; i++){
        c14->cd(i+2);
        if(i+2==1){
            TLatex heading;
            heading.SetTextAlign(22);
            heading.SetTextSize(0.1);
            heading.DrawLatex(0.5, 0.5, "#splitline{#Delta#phi correlation between}{#splitline{high p_{T} Pions (> 4 GeV/c)}{and K^{0} mesons}}");
        }else{
            fPiKDphi[i]->GetXaxis()->SetTitle("#Delta#phi");
            fPiKDphi[i]->SetLineWidth(0);
            fPiKDphi[i]->SetLineColor(0);
            fPiKDphi[i]->SetFillColor(16);
            fPiKDphi[i]->DrawCopy("HIST");
            fPiKDphi[i]->Sumw2();
            fPiKDphi[i]->SetMarkerStyle(21);
            fPiKDphi[i]->SetMarkerSize(1);
            fPiKDphi[i]->SetMarkerColor(1);
            fPiKDphi[i]->SetLineColor(1);
            fPiKDphi[i]->SetLineWidth(1);
            fPiKDphi[i]->Draw("SAME E1");
        }
    }
}
