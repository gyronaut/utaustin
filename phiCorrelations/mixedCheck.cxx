void mixedCheck(){

    TFile* corrFile = new TFile("eta20_corrected_phiCorrelations_mult_0_20.root");
    
    TH2D* corrhUSPeak = (TH2D*)hPhi2Dpeak;
    TH1D* corrhUSPeak_dphi = corrhUSPeak->ProjectionY("corrhUSPeak_dphi", corrhUSPeak->GetXaxis()->FindBin(-1.2), corrhUSPeak->GetXaxis()->FindBin(1.2));
    corrhUSPeak_dphi->SetLineWidth(2);
    corrhUSPeak_dphi->SetStats(kFALSE);
    corrhUSPeak_dphi->SetTitle("");
    corrhUSPeak_dphi->GetXaxis()->SetTitle("#Delta#varphi");
    TH2D* uncorrhUSPeak = (TH2D*)uncorrhPhi2Dpeak;
    uncorrhUSPeak->Sumw2();
    TH1D* uncorrhUSPeak_dphi = uncorrhUSPeak->ProjectionY("uncorrhUSPeak_dphi", uncorrhUSPeak->GetXaxis()->FindBin(-1.2), uncorrhUSPeak->GetXaxis()->FindBin(1.2));
    uncorrhUSPeak_dphi->SetLineColor(kRed);
    uncorrhUSPeak_dphi->SetLineWidth(2);

    TH2D* corrhLSPeak = (TH2D*)hKK2Dpeak;
    TH1D* corrhLSPeak_dphi = corrhLSPeak->ProjectionY("corrhLSPeak_dphi", corrhLSPeak->GetXaxis()->FindBin(-1.2), corrhLSPeak->GetXaxis()->FindBin(1.2));
    corrhLSPeak_dphi->SetLineWidth(2);
    corrhLSPeak_dphi->SetStats(kFALSE);
    corrhLSPeak_dphi->SetTitle("");
    corrhLSPeak_dphi->GetXaxis()->SetTitle("");
    TH2D* uncorrhLSPeak = (TH2D*)uncorrhKK2Dpeak;
    uncorrhLSPeak->Sumw2();
    TH1D* uncorrhLSPeak_dphi = uncorrhLSPeak->ProjectionY("uncorrhLSPeak_dphi", uncorrhLSPeak->GetXaxis()->FindBin(-1.2), uncorrhLSPeak->GetXaxis()->FindBin(1.2));
    uncorrhLSPeak_dphi->SetLineColor(kRed);
    uncorrhLSPeak_dphi->SetLineWidth(2);

    TH1D* USPeakRatio = corrhUSPeak_dphi->Clone("USPeakRatio");
    uncorrhUSPeak_dphi->Rebin();
    USPeakRatio->Divide(uncorrhUSPeak_dphi);
    USPeakRatio->Scale((uncorrhUSPeak_dphi->Integral())/(corrhUSPeak_dphi->Integral()));
    USPeakRatio->SetLineColor(kViolet);
    USPeakRatio->SetStats(kFALSE);
    USPeakRatio->SetTitle("");
    USPeakRatio->GetYaxis()->SetTitle("#frac{Corrected}{Uncorrected}");
    USPeakRatio->GetYaxis()->SetTitleOffset(1.3);

    TH1D* LSPeakRatio = corrhLSPeak_dphi->Clone("LSPeakRatio");
    uncorrhLSPeak_dphi->Rebin();
    LSPeakRatio->Divide(uncorrhLSPeak_dphi);
    LSPeakRatio->Scale((uncorrhLSPeak_dphi->Integral())/(corrhLSPeak_dphi->Integral()));
    LSPeakRatio->SetLineColor(kBlue);
    LSPeakRatio->SetStats(kFALSE);

    
    TLegend* uspeaklegend = new TLegend(0.4513, 0.6108, 0.8775, 0.7504);
    uspeaklegend->SetMargin(0.15);
    uspeaklegend->AddEntry(corrhUSPeak_dphi, "Acceptance Corrected", "le");
    uspeaklegend->AddEntry(uncorrhUSPeak_dphi, "Uncorrected", "le");

    TPaveText* uspeaktext = new TPaveText(0.1359, 0.7854, 0.5822, 0.8743, "NDC");
    uspeaktext->AddText("h-(K^{+}K^{-}) #Delta#varphi Correlation");
    uspeaktext->AddText("1.010 < m_{KK} < 1.030 GeV/c^{2}");
    uspeaktext->SetBorderSize(0);
    uspeaktext->SetFillColor(kWhite);

    TCanvas* cUSpeak = new TCanvas("cUSpeak", "cUSpeak", 50, 50, 600, 600);
    cUSpeak->cd();
    corrhUSPeak_dphi->Draw("H");
    uncorrhUSPeak_dphi->Draw("H SAME");
    uspeaklegend->Draw();
    uspeaktext->Draw();

    TLegend* lspeaklegend = new TLegend(0.4513, 0.6108, 0.8775, 0.7504);
    lspeaklegend->SetMargin(0.15);
    lspeaklegend->AddEntry(corrhLSPeak_dphi, "Acceptance Corrected", "le");
    lspeaklegend->AddEntry(uncorrhLSPeak_dphi, "Uncorrected", "le");

    TPaveText* lspeaktext = new TPaveText(0.1359, 0.7854, 0.5822, 0.8743, "NDC");
    lspeaktext->AddText("h-(K^{#pm}K^{#pm}) #Delta#varphi Correlation");
    lspeaktext->AddText("1.010 < m_{KK} < 1.030 GeV/c^{2}");
    lspeaktext->SetBorderSize(0);
    lspeaktext->SetFillColor(kWhite);

    TCanvas* cLSpeak = new TCanvas("cLSpeak", "cLSpeak", 50, 50, 600, 600);
    cLSpeak->cd();
    corrhLSPeak_dphi->Draw("H");
    uncorrhLSPeak_dphi->Draw("H SAME");
    lspeaklegend->Draw();
    lspeaktext->Draw();

    TLegend* ratiopeaklegend = new TLegend(0.4513, 0.7108, 0.8775, 0.8504);
    ratiopeaklegend->AddEntry(USPeakRatio, "Unlike-sign", "le");
    ratiopeaklegend->AddEntry(LSPeakRatio, "Like-sign", "le"); 

    TCanvas* cratiopeak = new TCanvas("cratiopeak", "cratiopeak", 50, 50, 600, 600);
    cratiopeak->cd();
    USPeakRatio->Draw("H");
    LSPeakRatio->Draw("H SAME");
    ratiopeaklegend->Draw();

    //Histograms for Corr/Uncorr for just center Delta-eta region (middle 2 bins);
    TH1D* USPeakRatioCenter = corrhUSPeak->ProjectionY("USPeakRatioCenter", corrhUSPeak->GetXaxis()->FindBin(-0.01), corrhUSPeak->GetXaxis()->FindBin(0.01));
    TH1D* uncorr= uncorrhUSPeak->ProjectionY("uncorr", uncorrhUSPeak->GetXaxis()->FindBin(-0.01), uncorrhUSPeak->GetXaxis()->FindBin(0.01));
    uncorr->Rebin();
    Float_t scale = uncorr->Integral()/USPeakRatioCenter->Integral();
    USPeakRatioCenter->Divide(uncorr);
    USPeakRatioCenter->Scale(scale);
    USPeakRatioCenter->SetLineWidth(2);
    USPeakRatioCenter->SetLineColor(kViolet);
    USPeakRatioCenter->SetTitle("");
    USPeakRatioCenter->SetStats(kFALSE);
    TH1D* LSPeakRatioCenter = corrhLSPeak->ProjectionY("LSPeakRatioCenter", corrhLSPeak->GetXaxis()->FindBin(-0.01), corrhLSPeak->GetXaxis()->FindBin(0.01));
    uncorr= uncorrhLSPeak->ProjectionY("uncorr", uncorrhLSPeak->GetXaxis()->FindBin(-0.01), uncorrhLSPeak->GetXaxis()->FindBin(0.01));
    uncorr->Rebin();
    Float_t scale = uncorr->Integral()/LSPeakRatioCenter->Integral();
    LSPeakRatioCenter->Divide(uncorr);
    LSPeakRatioCenter->Scale(scale);
    LSPeakRatioCenter->SetLineWidth(2);
    LSPeakRatioCenter->SetLineColor(kBlue);

    TCanvas* cratiopeakcenter = new TCanvas("cratiopeakcenter", "cratiopeakcenter", 50, 50, 600, 600);
    cratiopeakcenter->cd();
    USPeakRatioCenter->Draw("H");
    LSPeakRatioCenter->Draw("H SAME");
    ratiopeaklegend->Draw();


    TH2D* corrhUSRside = (TH2D*)hPhi2DRside;
    TH1D* corrhUSRside_dphi = corrhUSRside->ProjectionY("corrhUSRside_dphi", corrhUSRside->GetXaxis()->FindBin(-1.2), corrhUSRside->GetXaxis()->FindBin(1.2));
    corrhUSRside_dphi->SetLineWidth(2);
    corrhUSRside_dphi->SetStats(kFALSE);
    corrhUSRside_dphi->SetTitle("");
    corrhUSRside_dphi->GetXaxis()->SetTitle("#Delta#varphi");
    corrhUSRside_dphi->SetLineColor(kBlue-1);
    TH2D* uncorrhUSRside = (TH2D*)uncorrhPhi2DRside;
    uncorrhUSRside->Sumw2();
    TH1D* uncorrhUSRside_dphi = uncorrhUSRside->ProjectionY("uncorrhUSRside_dphi", uncorrhUSRside->GetXaxis()->FindBin(-1.2), uncorrhUSRside->GetXaxis()->FindBin(1.2));
    uncorrhUSRside_dphi->SetLineColor(kRed-2);
    uncorrhUSRside_dphi->SetLineWidth(2);

    TH2D* corrhLSRside = (TH2D*)hKK2DRside;
    TH1D* corrhLSRside_dphi = corrhLSRside->ProjectionY("corrhLSRside_dphi", corrhLSRside->GetXaxis()->FindBin(-1.2), corrhLSRside->GetXaxis()->FindBin(1.2));
    corrhLSRside_dphi->SetLineWidth(2);
    corrhLSRside_dphi->SetStats(kFALSE);
    corrhLSRside_dphi->SetTitle("");
    corrhLSRside_dphi->GetXaxis()->SetTitle("");
    corrhLSRside_dphi->SetLineColor(kBlue-1);
    TH2D* uncorrhLSRside = (TH2D*)uncorrhKK2DRside;
    uncorrhLSRside->Sumw2();
    TH1D* uncorrhLSRside_dphi = uncorrhLSRside->ProjectionY("uncorrhLSRside_dphi", uncorrhLSRside->GetXaxis()->FindBin(-1.2), uncorrhLSRside->GetXaxis()->FindBin(1.2));
    uncorrhLSRside_dphi->SetLineColor(kRed-2);
    uncorrhLSRside_dphi->SetLineWidth(2);

    TH1D* USRsideRatio = corrhUSRside_dphi->Clone("USRsideRatio");
    uncorrhUSRside_dphi->Rebin();
    USRsideRatio->Divide(uncorrhUSRside_dphi);
    USRsideRatio->Scale((uncorrhUSRside_dphi->Integral())/(corrhUSRside_dphi->Integral()));
    USRsideRatio->SetLineColor(kViolet-1);
    USRsideRatio->SetStats(kFALSE);
    USRsideRatio->SetTitle("");
    USRsideRatio->GetYaxis()->SetTitle("#frac{Corrected}{Uncorrected}");
    USRsideRatio->GetYaxis()->SetTitleOffset(1.3);

    TH1D* LSRsideRatio = corrhLSRside_dphi->Clone("LSRsideRatio");
    uncorrhLSRside_dphi->Rebin();
    LSRsideRatio->Divide(uncorrhLSRside_dphi);
    LSRsideRatio->Scale((uncorrhLSRside_dphi->Integral())/(corrhLSRside_dphi->Integral()));
    LSRsideRatio->SetLineColor(kBlue-2);
    LSRsideRatio->SetStats(kFALSE);

    
    TLegend* usrsidelegend = new TLegend(0.4513, 0.6108, 0.8775, 0.7504);
    usrsidelegend->SetMargin(0.15);
    usrsidelegend->AddEntry(corrhUSRside_dphi, "Acceptance Corrected", "le");
    usrsidelegend->AddEntry(uncorrhUSRside_dphi, "Uncorrected", "le");

    TPaveText* usrsidetext = new TPaveText(0.1359, 0.7854, 0.5822, 0.8743, "NDC");
    usrsidetext->AddText("h-(K^{+}K^{-}) #Delta#varphi Correlation");
    usrsidetext->AddText("1.040 < m_{KK} < 1.060 GeV/c^{2}");
    usrsidetext->SetBorderSize(0);
    usrsidetext->SetFillColor(kWhite);

    TCanvas* cUSrside = new TCanvas("cUSrside", "cUSrside", 50, 50, 600, 600);
    cUSrside->cd();
    corrhUSRside_dphi->Draw("H");
    uncorrhUSRside_dphi->Draw("H SAME");
    usrsidelegend->Draw();
    usrsidetext->Draw();

    TLegend* lsrsidelegend = new TLegend(0.4513, 0.6108, 0.8775, 0.7504);
    lsrsidelegend->SetMargin(0.15);
    lsrsidelegend->AddEntry(corrhLSRside_dphi, "Acceptance Corrected", "le");
    lsrsidelegend->AddEntry(uncorrhLSRside_dphi, "Uncorrected", "le");

    TPaveText* lsrsidetext = new TPaveText(0.1359, 0.7854, 0.5822, 0.8743, "NDC");
    lsrsidetext->AddText("h-(K^{#pm}K^{#pm}) #Delta#varphi Correlation");
    lsrsidetext->AddText("1.040 < m_{KK} < 1.060 GeV/c^{2}");
    lsrsidetext->SetBorderSize(0);
    lsrsidetext->SetFillColor(kWhite);

    TCanvas* cLSrside = new TCanvas("cLSrside", "cLSrside", 50, 50, 600, 600);
    cLSrside->cd();
    corrhLSRside_dphi->Draw("H");
    uncorrhLSRside_dphi->Draw("H SAME");
    lsrsidelegend->Draw();
    lsrsidetext->Draw();

    TLegend* ratiorsidelegend = new TLegend(0.4513, 0.7108, 0.8775, 0.8504);
    ratiorsidelegend->AddEntry(USRsideRatio, "Unlike-sign", "le");
    ratiorsidelegend->AddEntry(LSRsideRatio, "Like-sign", "le"); 

    TCanvas* cratiorside = new TCanvas("cratiorside", "cratiorside", 50, 50, 600, 600);
    cratiorside->cd();
    USRsideRatio->Draw("H");
    LSRsideRatio->Draw("H SAME");
    ratiorsidelegend->Draw();
 
}
