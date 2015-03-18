#include "TLegend.h"
#include <string>
#include <sstream>

void plot_phi_histo(string inputName){
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
    TH2F *DphiPhiH = (TH2F *)histoFile->Get("DphiPhiH");
    TH3F *HPhiDphi = (TH3F *)histoFile->Get("DphiHPhi"); 

    TH1F *hadronPtPEREVENT = (TH1F*) hadronPt->Clone("HadronpTperEvent");
    TH1F *phiPtPEREVENT = (TH1F*) phiPt->Clone("PhipTperEvent");
    //set-up "per event" histograms
    for(int i = 0; i < hadronPtPEREVENT->GetNbinsX(); i++){
	hadronPtPEREVENT->SetBinContent(i, (float)hadronPtPEREVENT->GetBinContent(i)/hadronPtPEREVENT->GetEntries());
	phiPtPEREVENT->SetBinContent(i, (float)phiPtPEREVENT->GetBinContent(i)/phiPtPEREVENT->GetEntries());
    }

    TH1D *fHPhiDphi[9];
    fHPhiDphi[0] = HPhiDphi->ProjectionZ("ptH1ptPhiInc", 20, 1000, 0, 1000);
    fHPhiDphi[0]->SetTitle("p_{T}^{h} > 1 GeV/c");
    fHPhiDphi[1] = HPhiDphi->ProjectionZ("ptH5ptPhiInc", 100, 1000, 0, 1000);
    fHPhiDphi[1]->SetTitle("p_{T}^{h} > 5 GeV/c");
    fHPhiDphi[2] = HPhiDphi->ProjectionZ("ptH7ptPhiInc", 140, 1000, 0, 1000);
    fHPhiDphi[2]->SetTitle("p_{T}^{h} > 7 GeV/c");
    fHPhiDphi[3] = HPhiDphi->ProjectionZ("ptH5_1ptPhi2", 100, 1000, 20, 40);
    fHPhiDphi[3]->SetTitle("p_{T}^{h} > 5,   1 < p_{T}^{#phi} < 2 GeV/c");
    fHPhiDphi[4] = HPhiDphi->ProjectionZ("ptH5_2ptPhi4", 100, 1000, 41, 80);
    fHPhiDphi[4]->SetTitle("p_{T}^{h} > 5,   2 < p_{T}^{#phi} < 4 GeV/c");
    fHPhiDphi[5] = HPhiDphi->ProjectionZ("ptH5_4ptPhi", 100, 1000, 81, 1000);
    fHPhiDphi[5]->SetTitle("p_{T}^{h} > 5,   p_{T}^{#phi} > 4 GeV/c");
    fHPhiDphi[6] = HPhiDphi->ProjectionZ("ptH7_1ptPhi2", 140, 1000, 20, 40);
    fHPhiDphi[6]->SetTitle("p_{T}^{h} > 7,   1 < p_{T}^{#phi} < 2 GeV/c");
    fHPhiDphi[7] = HPhiDphi->ProjectionZ("ptH7_2ptPhi4", 140, 1000, 41, 80);
    fHPhiDphi[7]->SetTitle("p_{T}^{h} > 7,   2 < p_{T}^{#phi} < 4 GeV/c");
    fHPhiDphi[8] = HPhiDphi->ProjectionZ("ptH7_4ptPhi", 140, 1000, 81, 1000);
    fHPhiDphi[8]->SetTitle("p_{T}^{h} > 7,   p_{T}^{#phi} > 4 GeV/c");

    TH1D* HPhiDphiInclusive = HPhiDphi->ProjectionZ("dphiptPhiInc",80,1000,40,1000);
    HPhiDphiInclusive->SetTitle("#Delta#phi for hadron-#phi(1020), p_{T}^{h} > 4 GeV/c,  p_{T}^{#phi} > 2 GeV/c");

    // Graph of dphi all pT
    TCanvas *cInc = new TCanvas("dphiInclusive", "#Delta#phi inclusive", 50, 50, 1000, 1000);
    cInc->cd();
    HPhiDphiInclusive->SetLineColor(1);
    HPhiDphiInclusive->SetFillColor(16);
    HPhiDphiInclusive->DrawCopy();
    HPhiDphiInclusive->Sumw2();
    HPhiDphiInclusive->SetMarkerStyle(21);
    HPhiDphiInclusive->SetMarkerSize(2);
    HPhiDphiInclusive->GetXaxis()->SetTitle("#Delta#phi");
    HPhiDphiInclusive->Draw("SAME E1");
    
    TF1* P = new TF1("P", "pol0", 1, 2);
    P->SetParameter(0, HPhiDphiInclusive->GetBinContent(1.5));
    HPhiDphiInclusive->Fit(P, "N", "", 1, 2);
    TF1* fit = new TF1("fit1", "pol0(0)+gaus(1)+gaus(4)", -1.5, 4.7);
    fit->SetParameter(0, P->GetParameter(0));
//    fit->SetParameter(1, 2000);
    fit->FixParameter(2, 0);
    fit->SetParameter(3, 1);
    fit->FixParameter(5, 3.1416);
    fit->SetParameter(6, 1.5);

    fit->SetLineColor(2);
    fit->SetLineWidth(10);
    HPhiDphiInclusive->Fit(fit, "", "SAME C");

    // Graph of Hadron pT Distribution
    TCanvas *c1 = new TCanvas("HadronPtCanvas", "Pt distribution for hadrons", 50, 50, 1000, 1000);
    c1->cd();
    c1->SetLogy();
    hadronPt->GetXaxis()->SetRangeUser(0, 15);
    hadronPt->SetLineColor(1);
    hadronPt->SetFillColor(16);
    hadronPt->Draw();

    // Graph of Phi(1020) pT Distribution
    TCanvas *c2 = new TCanvas("PhiPtCanvas", "Pt distribution for Phi mesons", 60, 60, 1000, 1000);
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

    // Graph of Phi(1020) pT Distribution for hadron pT > 5 GeV/c
    TCanvas *c3 = new TCanvas("PhiPtgt5Canvas", "Pt distribution for Phi mesons", 60, 60, 1000, 1000);
    c3->cd();
    TH1D* phiPt_Hgt5 = DphiHPhi->ProjectionY("phiPtgt5", 100, 1000, 0, 64);
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
    
    TH1F* phiPt_Hgt5_1_2 = (TH1F*) phiPt_Hgt5->Clone();
    phiPt_Hgt5_1_2->GetXaxis()->SetRangeUser(1, 2);
    phiPt_Hgt5_1_2->SetLineColor(1);
    phiPt_Hgt5_1_2->SetFillColor(2);
    phiPt_Hgt5_1_2->Draw("SAME");

    TH1F* phiPt_Hgt5_2_4 = (TH1F*) phiPt_Hgt5->Clone();
    phiPt_Hgt5_2_4->GetXaxis()->SetRangeUser(2, 4);
    phiPt_Hgt5_2_4->SetLineColor(1);
    phiPt_Hgt5_2_4->SetFillColor(3);
    phiPt_Hgt5_2_4->Draw("SAME");

    TH1F* phiPt_Hgt5_gt4 = (TH1F*) phiPt_Hgt5->Clone();
    phiPt_Hgt5_gt4->GetXaxis()->SetRangeUser(4, 15);
    phiPt_Hgt5_gt4->SetLineColor(1);
    phiPt_Hgt5_gt4->SetFillColor(4);
    phiPt_Hgt5_gt4->Draw("SAME");

    TLatex label_Hgt5;
    label_Hgt5.SetTextSize(0.03);

    Int_t maxHeight_Hgt5 = phiPt_Hgt5->GetBinContent(phiPt_Hgt5->GetMaximumBin());
    label_Hgt5.DrawLatex(10, maxHeight_Hgt5*.9, "Counts of #phi(1020) mesons:");

    string labelH_Hgt5_1_2 =  "#color[2]{1 < p_{T}^{#phi} < 2 GeV/c} : ";
    stringstream ss;
    ss << numphiPt_Hgt5_1_2;
    string number_Hgt5_1_2 = ss.str();
    string str_full_labelH_Hgt5_1_2 = labelH_Hgt5_1_2 + number_Hgt5_1_2;
    char* full_labelH_Hgt5_1_2 = str_full_labelH_Hgt5_1_2.c_str();
    label_Hgt5.DrawLatex(10, maxHeight_Hgt5*.82, full_labelH_Hgt5_1_2);

    ss.str(std::string());
    string labelH_Hgt5_2_4 =  "#color[3]{2 < p_{T}^{#phi} < 4 GeV/c} : ";
    ss << numphiPt_Hgt5_2_4;
    string number_Hgt5_2_4 = ss.str();
    string str_full_labelH_Hgt5_2_4 = labelH_Hgt5_2_4 + number_Hgt5_2_4;
    char* full_labelH_Hgt5_2_4 = str_full_labelH_Hgt5_2_4.c_str();
    label_Hgt5.DrawLatex(10, maxHeight_Hgt5*.74, full_labelH_Hgt5_2_4);

    ss.str(std::string());
    string labelH_Hgt5_gt4 = "#color[4]{p_{T}^{#phi} > 4 GeV/c} : ";
    ss << numphiPt_Hgt5_gt4;
    string number_Hgt5_gt4 = ss.str();
    string str_full_labelH_Hgt5_gt4 = labelH_Hgt5_gt4 + number_Hgt5_gt4;
    char* full_labelH_Hgt5_gt4 = str_full_labelH_Hgt5_gt4.c_str();
    label_Hgt5.DrawLatex(10, maxHeight_Hgt5*.66, full_labelH_Hgt5_gt4);

    // Graph of Phi(1020) pT Distribution for hadron pT > 7 GeV/c
    TCanvas *c4 = new TCanvas("PhiPtgt7Canvas", "Pt distribution for Phi mesons", 60, 60, 1000, 1000);
    c4->cd();
    TH1D* phiPt_Hgt7 = DphiHPhi->ProjectionY("phiPtgt7", 140, 1000, 0, 64);
    phiPt_Hgt7->SetTitle("p_{T}^{#phi} Distribution for events with a hadron with p_{T}^{h} > 7 GeV/c");
    phiPt_Hgt7->GetXaxis()->SetRangeUser(0, 15);
        phiPt_Hgt7->GetXaxis()->SetTitle("p_{T}^{#phi}");
    phiPt_Hgt7->SetLineColor(1);
    phiPt_Hgt7->SetFillColor(16);
    phiPt_Hgt7->Draw();
 
    Double_t numphiPt_Hgt7_1_2 = phiPt_Hgt7->Integral(20, 40);
    Double_t numphiPt_Hgt7_2_4 = phiPt_Hgt7->Integral(41, 80);
    Double_t numphiPt_Hgt7_gt4 = phiPt_Hgt7->Integral(81, 1000);   
    
    printf("Number of phi1_2 = %d\nNumber of phi2_4 = %d\nNumber of phigt4 = %d\n", numphiPt_Hgt7_1_2, numphiPt_Hgt7_2_4, numphiPt_Hgt7_gt4);
    
    TH1F* phiPt_Hgt7_1_2 = (TH1F*) phiPt_Hgt7->Clone();
    phiPt_Hgt7_1_2->GetXaxis()->SetRangeUser(1, 2);
    phiPt_Hgt7_1_2->SetLineColor(1);
    phiPt_Hgt7_1_2->SetFillColor(2);
    phiPt_Hgt7_1_2->Draw("SAME");

    TH1F* phiPt_Hgt7_2_4 = (TH1F*) phiPt_Hgt7->Clone();
    phiPt_Hgt7_2_4->GetXaxis()->SetRangeUser(2, 4);
    phiPt_Hgt7_2_4->SetLineColor(1);
    phiPt_Hgt7_2_4->SetFillColor(3);
    phiPt_Hgt7_2_4->Draw("SAME");

    TH1F* phiPt_Hgt7_gt4 = (TH1F*) phiPt_Hgt7->Clone();
    phiPt_Hgt7_gt4->GetXaxis()->SetRangeUser(4, 15);
    phiPt_Hgt7_gt4->SetLineColor(1);
    phiPt_Hgt7_gt4->SetFillColor(4);
    phiPt_Hgt7_gt4->Draw("SAME");

    TLatex label_Hgt7;
    label_Hgt7.SetTextSize(0.03);

    Int_t maxHeight_Hgt7 = phiPt_Hgt7->GetBinContent(phiPt_Hgt7->GetMaximumBin());
    label_Hgt7.DrawLatex(10, .9*maxHeight_Hgt7, "Counts of #phi(1020) mesons:");

    string labelH_Hgt7_1_2 =  "#color[2]{1 < p_{T}^{#phi} < 2 GeV/c} : ";
    stringstream ss;
    ss << numphiPt_Hgt7_1_2;
    string number_Hgt7_1_2 = ss.str();
    string str_full_labelH_Hgt7_1_2 = labelH_Hgt7_1_2 + number_Hgt7_1_2;
    char* full_labelH_Hgt7_1_2 = str_full_labelH_Hgt7_1_2.c_str();
    label_Hgt7.DrawLatex(10, .82*maxHeight_Hgt7, full_labelH_Hgt7_1_2);

    ss.str(std::string());
    string labelH_Hgt7_2_4 =  "#color[3]{2 < p_{T}^{#phi} < 4 GeV/c} : ";
    ss << numphiPt_Hgt7_2_4;
    string number_Hgt7_2_4 = ss.str();
    string str_full_labelH_Hgt7_2_4 = labelH_Hgt7_2_4 + number_Hgt7_2_4;
    char* full_labelH_Hgt7_2_4 = str_full_labelH_Hgt7_2_4.c_str();
    label_Hgt7.DrawLatex(10, .74*maxHeight_Hgt7, full_labelH_Hgt7_2_4);

    ss.str(std::string());
    string labelH_Hgt7_gt4 = "#color[4]{p_{T}^{#phi} > 4 GeV/c} : ";
    ss << numphiPt_Hgt7_gt4;
    string number_Hgt7_gt4 = ss.str();
    string str_full_labelH_Hgt7_gt4 = labelH_Hgt7_gt4 + number_Hgt7_gt4;
    char* full_labelH_Hgt7_gt4 = str_full_labelH_Hgt7_gt4.c_str();
    label_Hgt7.DrawLatex(10, .66*maxHeight_Hgt7, full_labelH_Hgt7_gt4);


    // Graph of H-Phi correlation for pT-Phi inclusive    
    TCanvas *c5 = new TCanvas("H-Phi correlation for pTh range, ptPhi inclusive", "H-Phi correlation for pTh range, ptPhi inclusive", 50, 50, 1000, 1000);
    c5->Divide(2,2);
    for(int i=-1; i<3; i++){
        c5->cd(i+2);
        if(i+2 == 1){
            TLatex heading;
            heading.SetTextAlign(22);
            heading.SetTextSize(0.1);
            heading.DrawLatex(0.5, 0.5, "#splitline{#Delta#phi correlation between}{#splitline{hadrons and #phi(1020) mesons,}{p_{T}^{#phi} inclusive}}");            
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

    // Graph of high-pt hadron (> 5 Gev/c) phi(1020) correlation
    TCanvas *c6 = new TCanvas("H-Phi correlation for pTh>5, ptPhi ranges", "H-Phi correlation for pTh>5, ptPhi ranges", 50, 50, 1000, 1000);
    c6->Divide(2,2);
    for(int i=2; i<6; i++){
        c6->cd(i-1);
        if(i-1==1){
            TLatex heading;
            heading.SetTextAlign(22);
            heading.SetTextSize(0.1);
            heading.DrawLatex(0.5, 0.5, "#splitline{#Delta#phi correlation between}{#splitline{high p_{T} hadrons (> 5 GeV/c)}{and #phi(1020) mesons}}");
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

    // Graph of high pT Hadron (> 7 GeV/c) and phi-meson correlation
    TCanvas *c7 = new TCanvas("H-Phi correlation for pTh>7, ptPhi ranges", "H-Phi correlation for pTh>7, ptPhi ranges", 50, 50, 1000, 1000);
    c7->Divide(2,2);
    for(int i=5; i<9; i++){
        c7->cd(i-4);
        if(i-4 == 1){
            TLatex heading;
            heading.SetTextAlign(22);
            heading.SetTextSize(0.1);
            heading.DrawLatex(0.5, 0.5, "#splitline{#Delta#phi correlation between}{#splitline{high p_{T} hadrons (> 7 GeV/c)}{and #phi(1020) mesons}}");
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
}
