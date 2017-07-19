void plotUncorrected(){
    TFile* eta20File = new TFile("eta20_corrected_phiCorrelations_mult_50_100.root");

    TH2D* eta20peak = uncorrhPhi2Dpeak->Clone("eta20peak");
    TH2D* eta20RSB = uncorrhPhi2DRside->Clone("eta20RSB");
    TH2D* eta20LSB = uncorrhPhi2DLside->Clone("eta20LSB");

    eta20peak->GetXaxis()->SetTitle("#Delta#eta");
    eta20peak->GetXaxis()->SetTitleSize(0.05);
    eta20peak->GetXaxis()->SetTitleOffset(1.3);
    eta20peak->GetYaxis()->SetTitle("#Delta#varphi");
    eta20peak->GetYaxis()->SetTitleSize(0.05);
    eta20peak->GetYaxis()->SetTitleOffset(1.3);
    eta20peak->SetTitle("");
    eta20peak->SetStats(kFALSE);

    eta20RSB->GetXaxis()->SetTitle("#Delta#eta");
    eta20RSB->GetXaxis()->SetTitleSize(0.05);
    eta20RSB->GetXaxis()->SetTitleOffset(1.3);
    eta20RSB->GetYaxis()->SetTitle("#Delta#varphi");
    eta20RSB->GetYaxis()->SetTitleSize(0.05);
    eta20RSB->GetYaxis()->SetTitleOffset(1.3);
    eta20RSB->SetTitle("");
    eta20RSB->SetStats(kFALSE);

    eta20LSB->GetXaxis()->SetTitle("#Delta#eta");
    eta20LSB->GetXaxis()->SetTitleSize(0.05);
    eta20LSB->GetXaxis()->SetTitleOffset(1.3);
    eta20LSB->GetYaxis()->SetTitle("#Delta#varphi");
    eta20LSB->GetYaxis()->SetTitleSize(0.05);
    eta20LSB->GetYaxis()->SetTitleOffset(1.3);
    eta20LSB->SetTitle("");
    eta20LSB->SetStats(kFALSE);

    TH1D* eta20peakEta = eta20peak->ProjectionX("eta20peakEta");
    eta20peakEta->GetXaxis()->SetTitleOffset(1.0);
    eta20peakEta->SetStats(kFALSE);
    TH1D* eta20peakPhi = eta20peak->ProjectionY("eta20peakPhi", eta20peak->GetXaxis()->FindBin(-1.5), eta20peak->GetXaxis()->FindBin(1.5));
    eta20peakPhi->SetLineColor(kBlue);
    eta20peakPhi->SetStats(kFALSE);
    TH1D* eta20peakPhiNarrow = eta20peak->ProjectionY("eta20peakPhiNarrow", eta20peak->GetXaxis()->FindBin(-1.2), eta20peak->GetXaxis()->FindBin(1.2));
    eta20peakPhiNarrow->SetLineColor(kViolet);
    TH1D* eta20peakPhiNarrowest = eta20peak->ProjectionY("eta20peakPhiNarrowest", eta20peak->GetXaxis()->FindBin(-1.0), eta20peak->GetXaxis()->FindBin(1.0));
    eta20peakPhiNarrowest->SetLineColor(kRed);

    TH1D* eta20RSBEta = eta20RSB->ProjectionX("eta20RSBEta");
    eta20RSBEta->GetXaxis()->SetTitleOffset(1.0);
    eta20RSBEta->SetStats(kFALSE);
    TH1D* eta20RSBPhi = eta20RSB->ProjectionY("eta20RSBPhi", eta20RSB->GetXaxis()->FindBin(-1.5), eta20RSB->GetXaxis()->FindBin(1.5));
    eta20RSBPhi->SetLineColor(kBlue);
    eta20RSBPhi->SetStats(kFALSE);
    TH1D* eta20RSBPhiNarrow = eta20RSB->ProjectionY("eta20RSBPhiNarrow", eta20RSB->GetXaxis()->FindBin(-1.2), eta20RSB->GetXaxis()->FindBin(1.2));
    eta20RSBPhiNarrow->SetLineColor(kViolet);
    TH1D* eta20RSBPhiNarrowest = eta20RSB->ProjectionY("eta20RSBPhiNarrowest", eta20RSB->GetXaxis()->FindBin(-1.0), eta20RSB->GetXaxis()->FindBin(1.0));
    eta20RSBPhiNarrowest->SetLineColor(kRed);

    TH1D* eta20LSBEta = eta20LSB->ProjectionX("eta20LSBEta");
    TH1D* eta20LSBPhi = eta20LSB->ProjectionY("eta20LSBPhi", eta20LSB->GetXaxis()->FindBin(-1.5), eta20LSB->GetXaxis()->FindBin(1.5));
    eta20LSBPhi->SetLineColor(kBlue);
    eta20LSBPhi->SetStats(kFALSE);
    TH1D* eta20LSBPhiNarrow = eta20LSB->ProjectionY("eta20LSBPhiNarrow", eta20LSB->GetXaxis()->FindBin(-1.2), eta20LSB->GetXaxis()->FindBin(1.2));
    eta20LSBPhiNarrow->SetLineColor(kViolet);
    TH1D* eta20LSBPhiNarrowest = eta20LSB->ProjectionY("eta20LSBPhiNarrowest", eta20LSB->GetXaxis()->FindBin(-1.0), eta20LSB->GetXaxis()->FindBin(1.0));
    eta20LSBPhiNarrowest->SetLineColor(kRed);


    //reset eta range to narrow view for 2D plotting and rebin
    eta20peak->Rebin2D(2,2);
    eta20RSB->Rebin2D(2,2);
    eta20LSB->Rebin2D(2,2);
    eta20peak->GetXaxis()->SetRangeUser(-1.2, 1.2);
    eta20RSB->GetXaxis()->SetRangeUser(-1.2, 1.2);
    eta20LSB->GetXaxis()->SetRangeUser(-1.2, 1.2);


    TH2D* LSeta20peak = uncorrhKK2Dpeak->Clone("LSeta20peak");
    TH2D* Lseta20RSB = uncorrhKK2DRside->Clone("LSeta20RSB");
    TH2D* LSeta20LSB = uncorrhKK2DLside->Clone("LSeta20LSB");

    LSeta20peak->GetXaxis()->SetTitle("#Delta#eta");
    LSeta20peak->GetXaxis()->SetTitleSize(0.05);
    LSeta20peak->GetXaxis()->SetTitleOffset(1.3);
    LSeta20peak->GetYaxis()->SetTitle("#Delta#varphi");
    LSeta20peak->GetYaxis()->SetTitleSize(0.05);
    LSeta20peak->GetYaxis()->SetTitleOffset(1.3);
    LSeta20peak->SetTitle("");
    LSeta20peak->SetStats(kFALSE);

    LSeta20RSB->GetXaxis()->SetTitle("#Delta#eta");
    LSeta20RSB->GetXaxis()->SetTitleSize(0.05);
    LSeta20RSB->GetXaxis()->SetTitleOffset(1.3);
    LSeta20RSB->GetYaxis()->SetTitle("#Delta#varphi");
    LSeta20RSB->GetYaxis()->SetTitleSize(0.05);
    LSeta20RSB->GetYaxis()->SetTitleOffset(1.3);
    LSeta20RSB->SetTitle("");
    LSeta20RSB->SetStats(kFALSE);

    LSeta20LSB->GetXaxis()->SetTitle("#Delta#eta");
    LSeta20LSB->GetXaxis()->SetTitleSize(0.05);
    LSeta20LSB->GetXaxis()->SetTitleOffset(1.3);
    LSeta20LSB->GetYaxis()->SetTitle("#Delta#varphi");
    LSeta20LSB->GetYaxis()->SetTitleSize(0.05);
    LSeta20LSB->GetYaxis()->SetTitleOffset(1.3);
    LSeta20LSB->SetTitle("");
    LSeta20LSB->SetStats(kFALSE);

    LSeta20peak->Rebin2D(2,2);
    LSeta20RSB->Rebin2D(2,2);
    LSeta20LSB->Rebin2D(2,2);
    LSeta20peak->GetXaxis()->SetRangeUser(-1.2, 1.2);
    LSeta20RSB->GetXaxis()->SetRangeUser(-1.2, 1.2);
    LSeta20LSB->GetXaxis()->SetRangeUser(-1.2, 1.2);




    TCanvas* ceta20peak = new TCanvas("ceta20peak", "ceta20peak", 50, 50, 800, 800);
    ceta20peak->Divide(2,2);
    ceta20peak->cd(1)->SetTheta(50);
    ceta20peak->cd(1)->SetPhi(50);
    eta20peak->Draw("SURF1");
    ceta20peak->cd(2);
    eta20peakEta->Draw("H");
    ceta20peak->cd(3);
    eta20peakPhi->Draw("H");
    eta20peakPhiNarrow->Draw("H SAME");
    eta20peakPhiNarrowest->Draw("H SAME");

    TCanvas* ceta20RSB = new TCanvas("ceta20RSB", "ceta20RSB", 60, 60, 800, 800);
    ceta20RSB->Divide(2,2);
    ceta20RSB->cd(1)->SetTheta(50);
    ceta20RSB->cd(1)->SetPhi(50);
    eta20RSB->Draw("SURF1");
    ceta20RSB->cd(2);
    eta20RSBEta->Draw("H");
    ceta20RSB->cd(3);
    eta20RSBPhi->Draw("H");
    eta20RSBPhiNarrow->Draw("H SAME");
    eta20RSBPhiNarrowest->Draw("H SAME");  

    TCanvas* ceta20LSB = new TCanvas("ceta20LSB", "ceta20LSB", 70, 70, 800, 800);
    ceta20LSB->Divide(2,2);
    ceta20LSB->cd(1)->SetTheta(50);
    ceta20LSB->cd(1)->SetPhi(50);
    eta20LSB->Draw("SURF1");
    ceta20LSB->cd(2);
    eta20LSBEta->Draw("H");
    ceta20LSB->cd(3);
    eta20LSBPhi->Draw("H");
    eta20LSBPhiNarrow->Draw("H SAME");
    eta20LSBPhiNarrowest->Draw("H SAME");

   
    TCanvas* cuncorrUSpeak = new TCanvas("cuncorrUSpeak", "cuncorrUSpeak", 50, 50, 800, 800);
    cuncorrUSpeak->cd()->SetTheta(50);
    cuncorrUSpeak->cd()->SetPhi(50);
    eta20peak->Draw("SURF1");

    TCanvas* cuncorrUSRside = new TCanvas("cuncorrUSRside", "cuncorrUSRside", 50, 50, 800, 800);
    cuncorrUSRside->cd()->SetTheta(50);
    cuncorrUSRside->cd()->SetPhi(50);
    eta20RSB->Draw("SURF1");

    TCanvas* cuncorrUSLside = new TCanvas("cuncorrUSLside", "cuncorrUSLside", 50, 50, 800, 800);
    cuncorrUSLside->cd()->SetTheta(50);
    cuncorrUSLside->cd()->SetPhi(50);
    eta20LSB->Draw("SURF1");


    TCanvas* cuncorrLSpeak = new TCanvas("cuncorrLSpeak", "cuncorrLSpeak", 50, 50, 800, 800);
    cuncorrLSpeak->cd()->SetTheta(50);
    cuncorrLSpeak->cd()->SetPhi(50);
    LSeta20peak->Draw("SURF1");

    TCanvas* cuncorrLSRside = new TCanvas("cuncorrLSRside", "cuncorrLSRside", 50, 50, 800, 800);
    cuncorrLSRside->cd()->SetTheta(50);
    cuncorrLSRside->cd()->SetPhi(50);
    LSeta20RSB->Draw("SURF1");

    TCanvas* cuncorrLSLside = new TCanvas("cuncorrLSLside", "cuncorrLSLside", 50, 50, 800, 800);
    cuncorrLSLside->cd()->SetTheta(50);
    cuncorrLSLside->cd()->SetPhi(50);
    LSeta20LSB->Draw("SURF1");


    TH2D* mixedUSpeak = uncorrhPhiMixed2Dpeak->Clone("MixedUSpeak");
    mixedUSpeak->SetStats(kFALSE);
    mixedUSpeak->SetTitle("");
    Float_t scale = 0.5*(mixedUSpeak->GetBinContent(mixedUSpeak->GetXaxis()->FindBin(-0.01), mixedUSpeak->GetYaxis()->FindBin(0.0)) + mixedUSpeak->GetBinContent(mixedUSpeak->GetXaxis()->FindBin(0.01), mixedUSpeak->GetYaxis()->FindBin(0.0)));
    mixedUSpeak->Scale(1.0/scale);
    TCanvas* cMixedUSpeak = new TCanvas("cMixedUSpeak", "cMixedUSpeak", 50, 50, 800, 800);
    cMixedUSpeak->cd()->SetTheta(50);
    cMixedUSpeak->cd()->SetPhi(50);
    mixedUSpeak->Draw("SURF1");

    TH2D* mixedLSpeak = uncorrhKKMixed2Dpeak->Clone("MixedLSpeak");
    mixedLSpeak->SetStats(kFALSE);
    mixedLSpeak->SetTitle("");
    TCanvas* cMixedLSpeak = new TCanvas("cMixedLSpeak", "cMixedLSpeak", 50, 50, 800, 800);
    cMixedLSpeak->cd()->SetTheta(50);
    cMixedLSpeak->cd()->SetPhi(50);
    mixedLSpeak->Draw("SURF1");

    TCanvas* cmixedratioRSB = new TCanvas("cmixedratioRSB", "cmixedratioRSB", 70, 70, 800, 800);
    cmixedratioRSB->Divide(2,2);
    cmixedratioRSB->cd(1)->SetTheta(50);
    cmixedratioRSB->cd(1)->SetPhi(50);
    mixedratioRSB->Draw("SURF1");
    cmixedratioRSB->cd(2);
    mixedratioRSBdeta->Draw("H");
    cmixedratioRSB->cd(3);
    mixedratioRSBdphi->Draw("H");

    TCanvas* cmixedratioPeak = new TCanvas("cmixedratioPeak", "cmixedratioPeak", 70, 70, 800, 800);
    cmixedratioPeak->Divide(2,2);
    cmixedratioPeak->cd(1)->SetTheta(50);
    cmixedratioPeak->cd(1)->SetPhi(50);
    mixedratioPeak->Draw("SURF1");
    cmixedratioPeak->cd(2);
    mixedratioPeakdeta->Draw("H");
    cmixedratioPeak->cd(3);
    mixedratioPeakdphi->Draw("H");

    
    TH2D* sameratioRSB = uncorrhPhi2DRside->Clone("sameratioRSB");
    sameratioRSB->SetTitle("");
    sameratioRSB->SetStats(kFALSE);
    Float_t scale = uncorrhKK2DRside->Integral()/sameratioRSB->Integral();
    sameratioRSB->Divide(uncorrhKK2DRside);
    sameratioRSB->Scale(scale);
    TH1D* sameratioRSBdeta = sameratioRSB->ProjectionX("sameratioRSBdeta");
    sameratioRSBdeta->SetStats(kFALSE);
    sameratioRSBdeta->SetTitle("");
    TH1D* sameratioRSBdphi = sameratioRSB->ProjectionY("sameratioRSBdphi", sameratioRSB->GetXaxis()->FindBin(-1.2), sameratioRSB->GetXaxis()->FindBin(1.2));
    sameratioRSBdphi->SetStats(kFALSE);
    sameratioRSBdphi->SetTitle("");

    TCanvas* csameratioRSB = new TCanvas("csameratioRSB", "csameratioRSB", 70, 70, 800, 800);
    csameratioRSB->Divide(2,2);
    csameratioRSB->cd(1)->SetTheta(50);
    csameratioRSB->cd(1)->SetPhi(50);
    sameratioRSB->Draw("SURF1");
    csameratioRSB->cd(2);
    sameratioRSBdeta->Draw("H");
    csameratioRSB->cd(3);
    sameratioRSBdphi->Draw("H");

    TH2D* sameratiopeak = uncorrhPhi2Dpeak->Clone("sameratiopeak");
    sameratiopeak->SetTitle("");
    sameratiopeak->SetStats(kFALSE);
    Float_t scale = uncorrhKK2Dpeak->Integral()/sameratiopeak->Integral();
    sameratiopeak->Divide(uncorrhKK2Dpeak);
    sameratiopeak->Scale(scale);
    TH1D* sameratiopeakdeta = sameratiopeak->ProjectionX("sameratiopeakdeta");
    sameratiopeakdeta->SetStats(kFALSE);
    sameratiopeakdeta->SetTitle("");
    TH1D* sameratiopeakdphi = sameratiopeak->ProjectionY("sameratiopeakdphi", sameratiopeak->GetXaxis()->FindBin(-1.2), sameratiopeak->GetXaxis()->FindBin(1.2));
    sameratiopeakdphi->SetStats(kFALSE);
    sameratiopeakdphi->SetTitle("");

    TCanvas* csameratiopeak = new TCanvas("csameratiopeak", "csameratiopeak", 70, 70, 800, 800);
    csameratiopeak->Divide(2,2);
    csameratiopeak->cd(1)->SetTheta(50);
    csameratiopeak->cd(1)->SetPhi(50);
    sameratiopeak->Draw("SURF1");
    csameratiopeak->cd(2);
    sameratiopeakdeta->Draw("H");
    csameratiopeak->cd(3);
    sameratiopeakdphi->Draw("H");

}
