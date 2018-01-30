#include <stdio.h>

KKInvMass(){
    gStyle->SetOptStat(0);
    
    TFile* file = new TFile("phiCorrelations_mult_20_50.root");
    //file->cd("PhiReconstruction");
    //THnSparseF* kkUSDist = (THnSparseF*)InvMass->FindObject("fkkUSDist");
    //THnSparseF* kkLSDist = (THnSparseF*)InvMass->FindObject("fkkLSDist");
    TList* list = (TList*) file->Get("phiCorr_mult_20_50");
    THnSparseF* kkUSDist= (THnSparseF*)list->FindObject("fkkUSDist");
    THnSparseF* kkLSDist= (THnSparseF*)list->FindObject("fkkLSDist");

    kkUSDist->GetAxis(0)->SetRangeUser(2.0, 4.0);
    kkLSDist->GetAxis(0)->SetRangeUser(2.0, 4.0);
    
    TH1D* USInvMass = kkUSDist->Projection(1);
    TH1D* LSInvMass = kkLSDist->Projection(1);
    
    float sidebandUS = (float)USInvMass->Integral(USInvMass->GetXaxis()->FindBin(1.04), USInvMass->GetXaxis()->FindBin(1.06));
    float sidebandLS = (float)LSInvMass->Integral(LSInvMass->GetXaxis()->FindBin(1.04), LSInvMass->GetXaxis()->FindBin(1.06));
    float scale = sidebandUS/sidebandLS;
    
    LSInvMass->Scale(scale);

    LSInvMass->SetLineColor(kRed);
    LSInvMass->SetLineWidth(2);
    USInvMass->SetLineWidth(3);
    
    USInvMass->SetTitle("");
    USInvMass->GetXaxis()->SetTitle("m_{KK} (GeV/c^{2})");
    USInvMass->GetXaxis()->SetTitleSize(0.05);
    USInvMass->SetLineColor(kBlack);
   
    TH1D* corrected = (TH1D*)USInvMass->Clone("corrected");
    corrected->Add(LSInvMass, -1.0);
    corrected->SetLineWidth(2);
    TF1* fit = new TF1("fit",  "[0]*TMath::Voigt(x - [1], [2], [3], 4) + pol2(4)",0.99, 1.1);
    fit->SetParameter(1, 1.020);
    fit->SetParameter(2, 0.0002);
    fit->SetParameter(0, 600);
    fit->FixParameter(3, 0.00426);
    fit->SetParLimits(1, 1.010, 1.030);
    fit->SetLineColor(kBlue);
    fit->SetLineStyle(7);
    fit->SetLineWidth(6);
    corrected->Fit(fit, "R");

    TLine* sbLine1 = new TLine(1.04, 0, 1.04, 44000);
    TLine* sbLine2 = new TLine(1.06, 0, 1.06, 44000);

    sbLine1->SetLineStyle(4);
    sbLine2->SetLineStyle(4);
    sbLine1->SetLineColor(kBlack);
    sbLine2->SetLineColor(kBlack);
    sbLine1->SetLineWidth(2);
    sbLine2->SetLineWidth(2);

    TLine* peak1 = new TLine(1.01, 0, 1.01, 44000);
    TLine* peak2 = new TLine(1.03, 0, 1.03, 44000);

    peak1->SetLineStyle(9);
    peak2->SetLineStyle(9);
    peak1->SetLineColor(kBlack);
    peak2->SetLineColor(kBlack);
    peak1->SetLineWidth(2);
    peak2->SetLineWidth(2);

    TLegend *lineLeg = new TLegend(0.4564, 0.5550, 0.8792, 0.6632);
    lineLeg->AddEntry(sbLine1, "Sideband Region", "l");
    //lineLeg->AddEntry(peak1, "#phi(1020) Peak Region", "l");

    TLegend *leg = new TLegend(0.4581, 0.3927, 0.8809, 0.5637);
    leg->AddEntry(USInvMass, "US Kaon Pairs");
    leg->AddEntry(LSInvMass, "#splitline{Est. BG}{(using LS Kaon Pairs)}");

    TPaveText *text = new TPaveText(0.5537, 0.7102, 0.8741, 0.8726, "NDC");
    text->AddText("ALICE");
    text->AddText("Work In Progress");
    text->AddText("p-Pb #sqrt{s_{NN}} = 5 TeV");
    text->SetFillColor(kWhite);
    text->SetBorderSize(0);

    TPaveText *pTText = new TPaveText(0.4564, 0.5550, 0.8792, 0.6632, "NDC");
    pTText->AddText("2.0 < p_{T}^{KK} < 4.0 GeV/c");
    pTText->SetFillColor(kWhite);
    pTText->SetBorderSize(0);
   
    TCanvas *c = new TCanvas("c", "c", 50, 50, 600, 600);
    c->cd();
    USInvMass->Draw();
    LSInvMass->Draw("SAME");
    sbLine1->Draw("SAME");
    sbLine2->Draw("SAME");
    //peak1->Draw("SAME");
    //peak2->Draw("SAME");
    //lineLeg->Draw();
    leg->Draw();
    text->Draw();
    pTText->Draw();

    TCanvas *c2 = new TCanvas("c2", "c2", 50, 50, 600, 600);
    c2->cd();
    corrected->Draw();
    text->Draw();
    pTText->Draw();

}
