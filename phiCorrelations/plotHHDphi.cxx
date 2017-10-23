plotHHDPhi(string inputfile){
    gStyle->SetOptStat(0);
    gStyle->SetOptFit(0);

    TFile* filehh = new TFile(inputfile.c_str());    
    string output = inputfile.substr(0,22);
    output+= "dphi.pdf";
    TH2D* hh2D_0_20 = (TH2D*)hh2D->Clone("hhdphi");
    TH1D* hhdphi_0_20 = (TH1D*)hh2D_0_20->ProjectionY("hhdphi_50_100", hh2D_0_20->GetXaxis()->FindBin(-1.2), hh2D_0_20->GetXaxis()->FindBin(1.2));
    hhdphi_0_20->Rebin();
    hhdphi_0_20->SetLineWidth(4);
    hhdphi_0_20->SetLineColor(kBlue+2);
    hhdphi_0_20->SetMarkerColor(kBlue+2);
    hhdphi_0_20->SetMarkerStyle(21);
    hhdphi_0_20->SetMarkerSize(2);
    hhdphi_0_20->GetXaxis()->SetTitle("#Delta#varphi");
    hhdphi_0_20->SetTitle("");
    hhdphi_0_20->Scale(1.0/(hhdphi_0_20->Integral()));
    hhdphi_0_20->GetYaxis()->SetTitle("Normalized Arb. Units");
    hhdphi_0_20->GetYaxis()->SetTitleOffset(1.60);
    hhdphi_0_20->GetXaxis()->SetTitleSize(0.05);
    hhdphi_0_20->GetXaxis()->SetTitleOffset(0.90);

    TFile* file0_20 = new TFile("~/phiStudies/LHC16q_FAST_V0A_moremixing/LS_trig_4_8_assoc_2_4_mixcorr_phiCorrelations_mult_0_20.root");    
    TH2D* hPhi2D_0_20 = (TH2D*)RLSsubhPhi2Dpeak->Clone("hPhidphi");
    TH1D* hPhidphi_0_20 = (TH1D*)hPhi2D_0_20->ProjectionY("hPhidphi_50_100", hPhi2D_0_20->GetXaxis()->FindBin(-1.2), hPhi2D_0_20->GetXaxis()->FindBin(1.2));
    //hPhidphi_0_20->Rebin();
    hPhidphi_0_20->SetLineWidth(4);
    hPhidphi_0_20->SetLineColor(kRed+2);
    hPhidphi_0_20->SetMarkerColor(kRed+2);
    hPhidphi_0_20->SetMarkerStyle(22);
    hPhidphi_0_20->SetMarkerSize(2);
    hPhidphi_0_20->GetXaxis()->SetTitle("#Delta#varphi");
    hPhidphi_0_20->SetTitle("");
    hPhidphi_0_20->Scale(1.0/(hPhidphi_0_20->Integral()));
    hPhidphi_0_20->GetYaxis()->SetTitle("Normalized Arb. Units");
    hPhidphi_0_20->GetYaxis()->SetTitleOffset(1.20);
    hPhidphi_0_20->GetXaxis()->SetTitleSize(0.05);
    hPhidphi_0_20->GetXaxis()->SetTitleOffset(0.90);

   /* TFile* file20_50 = new TFile("LS_eta20_corrected_phiCorrelations_mult_20_50.root");
    TH2D* hPhi2D_20_50 = (TH2D*)RLSsubhPhi2Dpeak->Clone("hphi2D_20_50");
    TH1D* hPhidphi_20_50 = (TH1D*)hPhi2D_20_50->ProjectionY("hPhidphi_20_50", hPhi2D_20_50->GetXaxis()->FindBin(-1.0), hPhi2D_20_50->GetXaxis()->FindBin(1.0));
    hPhidphi_20_50->Rebin();
    hPhidphi_20_50->SetLineWidth(2);
    hPhidphi_20_50->GetXaxis()->SetTitle("#Delta#varphi");
    hPhidphi_20_50->SetTitle("");

    TFile* file50_100 = new TFile("LS_eta20_corrected_phiCorrelations_mult_50_100.root");
    TH2D* hPhi2D_50_100 = (TH2D*)RLSsubhPhi2Dpeak->Clone("hphi2D_50_100");
    TH1D* hPhidphi_50_100 = (TH1D*)hPhi2D_50_100->ProjectionY("hPhidphi_50_100", hPhi2D_50_100->GetXaxis()->FindBin(-1.0), hPhi2D_50_100->GetXaxis()->FindBin(1.0));
    hPhidphi_50_100->SetLineWidth(2);
    hPhidphi_50_100->GetXaxis()->SetTitle("#Delta#varphi");
    hPhidphi_50_100->SetTitle("");
*/
    TF1 *corrFit = new TF1("corrFit", "gaus(0) + gaus(3) + pol0(6)", -1.4, 4.6);
    corrFit->FixParameter(6, 0.25*(hhdphi_0_20->GetBinContent(8)+hhdphi_0_20->GetBinContent(9)+hhdphi_0_20->GetBinContent(16)+hhdphi_0_20->GetBinContent(1)));
    //corrFit->SetParLimits(6, hhdphi_0_20->GetBinContent(8)*0.9, hhdphi_0_20->GetBinContent(8)*1.1);
    corrFit->SetParameter(0, hhdphi_0_20->GetBinContent(hhdphi_0_20->GetXaxis()->FindBin(0)) - corrFit->GetParameter(6));
    corrFit->SetParLimits(0, corrFit->GetParameter(0)*0.5, corrFit->GetParameter(0)*1.5);
    corrFit->SetParameter(1, 0.0);
    corrFit->SetParLimits(1, -0.5, 0.5);
    corrFit->SetParameter(2, 0.5);
    corrFit->SetParLimits(2, 0.2, 0.9);
    corrFit->SetParameter(3, hhdphi_0_20->GetBinContent(hhdphi_0_20->GetXaxis()->FindBin(3.14)) - corrFit->GetParameter(6));
    corrFit->SetParLimits(3, corrFit->GetParameter(3)*0.5, corrFit->GetParameter(3)*1.5);
    corrFit->SetParameter(4, 3.14);
    corrFit->SetParLimits(4, 3.0, 3.25);
    corrFit->SetParameter(5, 0.5);
    corrFit->SetParLimits(5, 0.2, 0.9);

    corrFit->SetParNames("peak1", "mean1", "sigma1", "peak2", "mean2", "sigma2", "BG");

    corrFit->SetLineColor(kBlue);
    corrFit->SetLineWidth(6);
    corrFit->SetLineStyle(7);

    TF1 *corrFit2 = new TF1("corrFit2", "gaus(0) + gaus(3) + pol0(6)", -1.4, 4.6);
    corrFit2->FixParameter(6, 0.25*(hPhidphi_0_20->GetBinContent(8)+hPhidphi_0_20->GetBinContent(9)+hPhidphi_0_20->GetBinContent(16)+hPhidphi_0_20->GetBinContent(1)));
    //corrFit->SetParLimits(6, hPhidphi_0_20->GetBinContent(8)*0.9, hPhidphi_0_20->GetBinContent(8)*1.1);
    corrFit2->SetParameter(0, hPhidphi_0_20->GetBinContent(hPhidphi_0_20->GetXaxis()->FindBin(0)) - corrFit->GetParameter(6));
    corrFit2->SetParLimits(0, corrFit->GetParameter(0)*0.5, corrFit->GetParameter(0)*1.5);
    corrFit2->SetParameter(1, 0.0);
    corrFit2->SetParLimits(1, -0.5, 0.5);
    corrFit2->SetParameter(2, 0.5);
    corrFit2->SetParLimits(2, 0.2, 0.9);
    corrFit2->SetParameter(3, hPhidphi_0_20->GetBinContent(hPhidphi_0_20->GetXaxis()->FindBin(3.14)) - corrFit->GetParameter(6));
    corrFit2->SetParLimits(3, corrFit->GetParameter(3)*0.5, corrFit->GetParameter(3)*1.5);
    corrFit2->SetParameter(4, 3.14);
    corrFit2->SetParLimits(4, 3.0, 3.25);
    corrFit2->SetParameter(5, 0.5);
    corrFit2->SetParLimits(5, 0.2, 0.9);

    corrFit2->SetParNames("peak1", "mean1", "sigma1", "peak2", "mean2", "sigma2", "BG");

    corrFit2->SetLineColor(kRed);
    corrFit2->SetLineWidth(6);
    corrFit2->SetLineStyle(7);

    /*
    TF1* corrFit_20_50 = (TF1*)corrFit->Clone("corrFit_20_50");
    corrFit_20_50->FixParameter(6, 0.25*(hPhidphi_20_50->GetBinContent(8)+hPhidphi_20_50->GetBinContent(9)+hPhidphi_20_50->GetBinContent(16)+hPhidphi_20_50->GetBinContent(1)));
    corrFit_20_50->SetLineColor(kBlue-5);
    //corrFit_20_50->SetParLimits(6, hPhidphi_20_50->GetBinContent(8)*0.9, hPhidphi_20_50->GetBinContent(8)*1.1);
    corrFit_20_50->SetParameter(0, hPhidphi_20_50->GetMaximum());
    TF1* corrFIt_50_100 = (TF1*)corrFit->Clone("corrFit_50_100");
    corrFit_50_100->FixParameter(6, 0.25*(hPhidphi_50_100->GetBinContent(8)+hPhidphi_50_100->GetBinContent(9)+hPhidphi_50_100->GetBinContent(16)+hPhidphi_50_100->GetBinContent(1)));
    corrFit_50_100->SetLineColor(kViolet+6);
    //corrFit_50_100->SetParLimits(6, hPhidphi_50_100->GetBinContent(8)*0.9, hPhidphi_50_100->GetBinContent(8)*1.1);
    corrFit_50_100->SetParameter(0, hPhidphi_50_100->GetMaximum());
*/
    hhdphi_0_20->Fit("corrFit", "R");
    hPhidphi_0_20->Fit("corrFit2", "R");
//    hPhidphi_20_50->Fit("corrFit_20_50", "R");
//    hPhidphi_50_100->Fit("corrFit_50_100", "R");

    TF1 *gaus1 = new TF1("gaus1", "gaus(0)", -1.4, 4.6);
    TF1 *gaus2 = new TF1("gaus2", "gaus(0)", -1.4, 4.6);
    TF1 *bg = new TF1("bg", "pol0(0)", -1.4, 4.6);
    gaus1->SetParameter(0, corrFit->GetParameter(0));
    gaus1->SetParameter(1, corrFit->GetParameter(1));
    gaus1->SetParameter(2, corrFit->GetParameter(2));
    gaus2->SetParameter(0, corrFit->GetParameter(3));
    gaus2->SetParameter(1, corrFit->GetParameter(4));
    gaus2->SetParameter(2, corrFit->GetParameter(5));
    bg->SetParameter(0, corrFit->GetParameter(6));    

    //float signal1 = gaus1->Integral(gaus1->GetParameter(1) - (gaus1->GetParameter(2)*2), gaus1->GetParameter(1) - (gaus1->GetParameter(2)*2));
    float signal1 = gaus1->Integral(-1.4, 4.6);
    //float signal2 = gaus2->Integral(gaus2->GetParameter(1) - (gaus2->GetParameter(2)*2), gaus2->GetParameter(1) - (gaus2->GetParameter(2)*2));
    float signal2 = gaus2->Integral(-1.4, 4.6);
    //float bg1 = bg->Integral(gaus1->GetParameter(1) - (gaus1->GetParameter(2)*2), gaus1->GetParameter(1) - (gaus1->GetParameter(2)*2));
    //float bg2 = bg->Integral(gaus2->GetParameter(1) - (gaus2->GetParameter(2)*2), gaus2->GetParameter(1) - (gaus2->GetParameter(2)*2));
    float bg1 = bg->Integral(-1.4, 4.6);  

    float signalOverBG = (signal1+signal2)/(bg1);

    printf("\n================\nsignal over BG: %e\n===================\n", signalOverBG);

    TPaveText *text = new TPaveText(0.4815, 0.7056, 0.8658, 0.8551, "NDC");
    text->AddText("ALICE Work in Progress");
    text->AddText("p-Pb #sqrt{#it{s}_{NN}} = 5 TeV");
    text->AddText("0%-20% Multiplicity");
    text->SetTextSizePixels(20);
    text->SetFillColor(kWhite);

    TPaveText *text2 = new TPaveText(0.6, 0.9, 0.85, 0.85, "NDC");
    text2->AddText("trigger: 2.0 < #it{p}_{T}^{h} < 4.0 GeV/c");
    text2->AddText("assoc: 1.0 < #it{p}_{T}^{a} < 2.0 GeV/c");
    text2->SetTextSizePixels(18);
    text2->SetFillColor(kWhite);

    TLegend  *legend = new TLegend(0.3791, 0.1518, 0.8772, 0.2688);
    legend->SetMargin(0.15);
    legend->AddEntry(corrFit2, "Hadron-#phi(1020) Correlation", "l");
    legend->AddEntry(corrFit, "Hadron-hadron Correlations", "l");
    TCanvas *c0_20 = new TCanvas("c0_20", "c0_20", 50, 50, 550, 600);
    c0_20->cd();
    c0_20->SetMargin(0.12, 0.05, 0.1, 0.05);
    hhdphi_0_20->GetYaxis()->SetRangeUser(0.03, 0.11);
    hhdphi_0_20->Draw("E0 X0");
    hPhidphi_0_20->Draw("E0 X0 SAME");
    legend->Draw();
    text->Draw("SAME");
    text2->Draw("SAME");
    //c0_20->Print(output.c_str(),"pdf");

  /*  TCanvas *c20_50 = new TCanvas("c20_50", "c20_50", 50, 50, 600, 600);
    c20_50->cd();
    hPhidphi_20_50->Draw();

    TCanvas *c50_100 = new TCanvas("c50_100", "c50_100", 50, 50, 600, 600);
    c50_100->cd();
    hPhidphi_50_100->Draw();
    */
}
