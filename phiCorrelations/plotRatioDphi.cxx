plotRatioDphi(){
    gStyle->SetOptStat(0);
    gStyle->SetOptFit(0);

    TFile *hhFile = new TFile("~/phiStudies/LHC16q_FAST_hh_0_20/trig_4_8_assoc_2_4_hh_phiCorrelations_mult_0_20.root");
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
    hhdphi_0_20->GetYaxis()->SetTitleOffset(1.60);
    hhdphi_0_20->GetXaxis()->SetTitleSize(0.05);
    hhdphi_0_20->GetXaxis()->SetTitleOffset(0.90);


    TFile* phiFile = new TFile("~/phiStudies/LHC16q_FAST_0_20_widecuts/LS_trig_4_8_assoc_2_4_mixcorr_phiCorrelations_mult_0_20.root");    
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
    hPhidphi_0_20->GetYaxis()->SetTitleOffset(1.20);
    hPhidphi_0_20->GetXaxis()->SetTitleSize(0.05);
    hPhidphi_0_20->GetXaxis()->SetTitleOffset(0.90);

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
    corrFit2->SetParLimits(5, 0.1, 0.9);

    corrFit2->SetParNames("peak1", "mean1", "sigma1", "peak2", "mean2", "sigma2", "BG");

    corrFit2->SetLineColor(kRed);
    corrFit2->SetLineWidth(6);
    corrFit2->SetLineStyle(7);

    hhdphi_0_20->Fit("corrFit", "R0");
    hPhidphi_0_20->Fit("corrFit2", "R0");

    TF1 *hhBG = new TF1("hhBG", "pol0(0)", -1.4, 4.6);
    hhBG->SetParLimits(0, 0.00001, 10000000.0);
    hhBG->SetParameter(0, -1.0*corrFit->GetParameter(6));

    TF1 *hphiBG = new TF1("hphiBG", "pol0(0)", -1.4, 4.6);
    hphiBG->SetParLimits(0, 0.00001, 10000000.0);
    hphiBG->SetParameter(0, -1.0*corrFit2->GetParameter(6));


    /*hhdphi_0_20->Add(hhBG);
    hPhidphi_0_20->Add(hphiBG);
    //re-fit the BG subtracted distributions for plotting
    hhdphi_0_20->Fit("corrFit", "R");
    hPhidphi_0_20->Fit("corrFit2", "R");
    */

    TH1D *ratio = hPhidphi_0_20->Clone("ratio");
    ratio->Divide(hhdphi_0_20);
    

    TLegend  *legend = new TLegend(0.3791, 0.1518, 0.8772, 0.2688);
    legend->SetMargin(0.15);
    legend->AddEntry(corrFit2, "Hadron-#phi(1020) Correlation", "l");
    legend->AddEntry(corrFit, "Hadron-hadron Correlations", "l");
    
    TCanvas *c0_20 = new TCanvas("c0_20", "c0_20", 50, 50, 550, 600);
    c0_20->cd();
    c0_20->SetMargin(0.12, 0.05, 0.1, 0.05);
    //hhdphi_0_20->GetYaxis()->SetRangeUser(0.03, 0.11);
    hhdphi_0_20->Draw("E0 X0");
    corrFit->Draw();
    corrFit2->Draw("SAME");
    hPhidphi_0_20->Draw("E0 X0 SAME");
    legend->Draw();

    TCanvas *cratio = new TCanvas("cratio", "cratio", 50, 50, 550, 600);
    cratio->cd();
    cratio->SetMargin(0.12, 0.05, 0.1, 0.05);
    ratio->Draw("E0 X0");


}

