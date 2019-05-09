#include <stdio.h>

void KKInvMass(){
    gStyle->SetOptStat(0);
    gStyle->SetOptFit(0);

    TFile* file = new TFile("~/phiStudies/results_onlineEff/Combined/TOTAL_hphi_0_20_50_80.root");
    TList* list = (TList*) file->Get("phiCorr_mult_0_20");
    THnSparseF* kkUSDist= (THnSparseF*)list->FindObject("fkkUSDist");
    THnSparseF* kkLSDist= (THnSparseF*)list->FindObject("fkkLSDist");

    kkUSDist->GetAxis(0)->SetRangeUser(2.0, 4.0);
    kkLSDist->GetAxis(0)->SetRangeUser(2.0, 4.0);
    
    TH1D* USInvMass = kkUSDist->Projection(1);
    TH1D* LSInvMass = kkLSDist->Projection(1);
    TH1D* origLSInvMass = (TH1D*)LSInvMass->Clone("origLSInvMass");
    origLSInvMass->SetLineColor(kRed);
    origLSInvMass->SetLineWidth(2);
    
    Double_t sidebandUS = (Double_t)(USInvMass->Integral(USInvMass->GetXaxis()->FindBin(1.04), USInvMass->GetXaxis()->FindBin(1.06)) + USInvMass->Integral(USInvMass->GetXaxis()->FindBin(0.995), USInvMass->GetXaxis()->FindBin(1.005)));
    Double_t sidebandLS = (Double_t)(LSInvMass->Integral(LSInvMass->GetXaxis()->FindBin(1.04), LSInvMass->GetXaxis()->FindBin(1.06)) + LSInvMass->Integral(LSInvMass->GetXaxis()->FindBin(0.995), LSInvMass->GetXaxis()->FindBin(1.005)));
    Double_t scale = sidebandUS/sidebandLS;
    
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
    corrected->SetMarkerSize(2);
    corrected->SetMarkerStyle(34);
    TF1* fit = new TF1("fit",  "[0]*TMath::Voigt(x - [1], [2], [3], 4) + pol2(4)",0.99, 1.07);
    fit->SetParameter(1, 1.020);
    fit->SetParameter(2, 0.0002);
    fit->SetParameter(0, 600);
    fit->FixParameter(3, 0.00426);
    fit->SetParLimits(1, 1.010, 1.030);
    fit->SetLineColor(kBlue);
    fit->SetLineStyle(7);
    fit->SetLineWidth(6);
    corrected->Fit(fit, "R");

    TF1* voigtFit = new TF1("voigtFit", "[0]*TMath::Voigt(x - [1], [2], [3], 4)", 0.99, 1.07);
    voigtFit->SetParameters(fit->GetParameter(0), fit->GetParameter(1), fit->GetParameter(2), fit->GetParameter(3));

    TF1* bgFit = new TF1("bgFit", "pol2(0)", 0.99, 1.07);
    bgFit->SetParameters(fit->GetParameter(4), fit->GetParameter(5), fit->GetParameter(6));

    Float_t fullInt = voigtFit->Integral(0.99, 1.07);
    Float_t onlyHistInt = corrected->Integral(corrected->GetXaxis()->FindBin(0.9901), corrected->GetXaxis()->FindBin(1.06999), "width");
    Float_t onlyBGInt = bgFit->Integral(0.99, 1.07);
    Float_t fullHistInt = corrected->Integral(corrected->GetXaxis()->FindBin(0.9901), corrected->GetXaxis()->FindBin(1.06999), "width") - bgFit->Integral(0.99, 1.07);
    Float_t widePct = voigtFit->Integral(1.01, 1.03);
    widePct = widePct/fullInt;
    Float_t peakInt = voigtFit->Integral(1.014, 1.026);
    Float_t bgpeakInt = bgFit->Integral(1.014, 1.026);
    Float_t bg2peak = bgpeakInt/peakInt;
    
    printf("=================\n\nfit integral: %e,    hist integral: %e,   only hist: %e,    only BG: %e,    ratio to residual: %4.2f%% \n=====================\n", fullInt, fullHistInt, onlyHistInt, onlyBGInt, bg2peak*100);

    Float_t diffPct[11];
    Float_t diffInt[11];
    Float_t sigBG[11];
    Float_t purity[11];
    TH1D* effHist = new TH1D("effHist", "Efficiency & Purity of Different Mass Ranges;Invariant Mass Range (MeV/c^{2});Percent", 11, 0., 11.);
    TH1D* purityHist = new TH1D("purityHist", "Purity;Invariant Mass Range;#frac{Signal}{Signal + BG}", 11, 0., 11.);

    for(int i = 0; i < 11; i++){
        diffPct[i] = Float_t(voigtFit->Integral(1.008 + 0.001*i, 1.032 - 0.001*i))/fullInt;
        diffInt[i] = Float_t(corrected->Integral(corrected->GetXaxis()->FindBin(1.0081 + 0.001*i), corrected->GetXaxis()->FindBin(1.032 - 0.0001 - 0.001*i), "width") - bgFit->Integral(1.008 + 0.001*i, 1.032 - 0.001*i))/fullHistInt;
        sigBG[i] = Float_t(corrected->Integral(corrected->GetXaxis()->FindBin(1.0081 + 0.001*i), corrected->GetXaxis()->FindBin(1.032 - 0.0001 - 0.001*i), "width") - bgFit->Integral(1.008 + 0.001*i, 1.032 - 0.001*i))/Float_t(LSInvMass->Integral(LSInvMass->GetXaxis()->FindBin(1.0081+0.001*i), LSInvMass->GetXaxis()->FindBin(1.0319 - 0.001*i),"width") + bgFit->Integral(1.008 + 0.001*i, 1.032 - 0.001*i));
        purity[i] = Float_t(corrected->Integral(corrected->GetXaxis()->FindBin(1.0081 + 0.001*i), corrected->GetXaxis()->FindBin(1.032 - 0.0001 - 0.001*i), "width") - bgFit->Integral(1.008 + 0.001*i, 1.032 - 0.001*i))/Float_t(USInvMass->Integral(USInvMass->GetXaxis()->FindBin(1.0081+0.001*i), USInvMass->GetXaxis()->FindBin(1.0319 - 0.001*i),"width") + bgFit->Integral(1.008 + 0.001*i, 1.032 - 0.001*i));
        effHist->SetBinContent(i+1, diffInt[i]*100.0);
        effHist->GetXaxis()->SetBinLabel(i+1, Form("%i", (int)(((1.032 - 0.001*i) - (1.008 + 0.001*i))*1000)));
        purityHist->SetBinContent(i+1, sigBG[i]);
    }


    TCanvas *ceff = new TCanvas("ceff", "ceff", 50, 50, 900, 600);
    ceff->cd();
    effHist->GetYaxis()->SetRangeUser(0.0, 100.0);
    effHist->SetLineWidth(2);
    effHist->Draw();
    ceff->Update();

    Float_t rightmax = 1.1*purityHist->GetMaximum();
    Float_t axisscale = gPad->GetUymax()/rightmax;
    purityHist->SetLineColor(kRed);
    purityHist->SetLineWidth(2);
    purityHist->Scale(axisscale);
    purityHist->Draw("HIST SAME");

    TGaxis *axis = new TGaxis(gPad->GetUxmax(), gPad->GetUymin(), gPad->GetUxmax(), gPad->GetUymax(),0,rightmax,510,"+L");
    axis->SetLineColor(kRed);
    axis->SetTextColor(kRed);
    axis->SetLabelColor(kRed);
    axis->SetLabelSize(0.035);
    axis->SetLabelFont(42);
    axis->SetTitle("#frac{Signal}{BG}");
    axis->Draw();
    ceff->Update();

    for(int j = 0; j < 11; j++){
        printf(" fit range %4.3f to %4.3f: %4.2f%%\n", 1.008+0.001*j, 1.032 -0.001*j, diffPct[j]*100.0);
        printf("hist range %4.3f to %4.3f: %4.2f%%, Signal/(Signal+ BG): %4.2f Signal/BG: %4.2f\n\n", 1.008+0.001*j, 1.032 -0.001*j, diffInt[j]*100.0, purity[j], sigBG[j]);
    }



    TLine* sbLine1 = new TLine(1.04, 0, 1.04, 44000);
    TLine* sbLine2 = new TLine(1.06, 0, 1.06, 44000);

    sbLine1->SetLineStyle(4);
    sbLine2->SetLineStyle(4);
    sbLine1->SetLineColor(kBlack);

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
    //text->AddText("Work In Progress");
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


    TLegend *corrleg = new TLegend(0.46, 0.39, 0.88, 0.57);
    corrleg->AddEntry(corrected, "Corrected US Inv. Mass", "p");
    corrleg->AddEntry(fit, "Inv. Mass Fit", "l");
    corrleg->AddEntry(bgFit, "Residual BG Fit", "l");

    TCanvas *c2 = new TCanvas("c2", "c2", 50, 50, 600, 600);
    c2->cd();
    corrected->Draw();
    bgFit->Draw("SAME");
    text->Draw();
    pTText->Draw();
    corrleg->Draw();

    TH1D* wideInvMass = (TH1D*)USInvMass->Clone("wideInvMass");
    wideInvMass->SetFillColor(36);
    TH1D* narrowInvMass=(TH1D*)USInvMass->Clone("narrowInvMass");
    narrowInvMass->SetFillColor(30);
    TH1D* LSBinvmass = (TH1D*)USInvMass->Clone("LSBinvmass");
    LSBinvmass->SetFillColor(kGray+2);
    TH1D* RSBinvmass = (TH1D*)USInvMass->Clone("RSBinvmass");
    RSBinvmass->SetFillColor(kGray+2);

    for(int i = 1; i<= wideInvMass->GetXaxis()->GetNbins(); i++){
        if(TMath::Abs(wideInvMass->GetXaxis()->GetBinCenter(i) - 1.020) > 0.010){
            wideInvMass->SetBinContent(i, 0.);
        }
        if(TMath::Abs(narrowInvMass->GetXaxis()->GetBinCenter(i) - 1.020) > 0.006){
            narrowInvMass->SetBinContent(i, 0.);
        }
        if(narrowInvMass->GetBinCenter(i) < 0.995 || narrowInvMass->GetBinCenter(i) > 1.005){
            LSBinvmass->SetBinContent(i, 0.);
        }
        if(narrowInvMass->GetBinCenter(i) < 1.040 || narrowInvMass->GetBinCenter(i) > 1.060){
            RSBinvmass->SetBinContent(i, 0.);
        }
    }

    TCanvas* c3 = new TCanvas("c3", "c3", 50, 50, 600, 600);
    c3->cd();
    USInvMass->Draw();
    wideInvMass->Draw("SAME");
    USInvMass->Draw("SAME");
    LSInvMass->Draw("SAME");
    leg->Draw();
    
    TCanvas* c4 = new TCanvas("c4", "c4", 50, 50, 600, 600);
    c4->cd();
    USInvMass->Draw();
    narrowInvMass->Draw("SAME");
    LSBinvmass->Draw("SAME");
    RSBinvmass->Draw("SAME");
    USInvMass->Draw("SAME");
    LSInvMass->Draw("SAME");
    leg->Draw();

    TLegend *leg2 = new TLegend(0.4581, 0.3927, 0.8809, 0.5637);
    leg2->AddEntry(USInvMass, "US Kaon Pairs", "l");
    
    TCanvas* c5 = new TCanvas("c5", "c5", 50, 50, 600, 600);
    c5->cd();
    USInvMass->Draw();
    narrowInvMass->Draw("SAME");
    LSBinvmass->Draw("SAME");
    RSBinvmass->Draw("SAME");
    USInvMass->Draw("SAME");
    leg2->Draw();

    TLegend *leg3 = new TLegend(0.4581, 0.3927, 0.8809, 0.5637);
    leg3->AddEntry(USInvMass, "US Kaon Pairs", "l");
    leg3->AddEntry(origLSInvMass, "LS Kaon Pairs", "l");

    TCanvas* c6 = new TCanvas("c6", "c6", 50, 50, 600, 600);
    c6->cd();
    USInvMass->Draw();
    narrowInvMass->Draw("SAME");
    LSBinvmass->Draw("SAME");
    RSBinvmass->Draw("SAME");
    USInvMass->Draw("SAME");
    origLSInvMass->Draw("SAME");
    leg3->Draw();

}
