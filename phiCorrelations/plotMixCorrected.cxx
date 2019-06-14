void plotMixCorrected(string inputfile){
    TFile* eta20File = new TFile(inputfile.c_str());

    TH2D* hPhi2Dpeak = (TH2D*)eta20File->Get("hPhi2Dpeak");
    TH2D* hPhi2DRside = (TH2D*)eta20File->Get("hPhi2DRside");
    TH2D* hPhi2DLside = (TH2D*)eta20File->Get("hPhi2DLside");

    TH2D* eta20peak = (TH2D*)hPhi2Dpeak->Clone("eta20peak");
    TH2D* eta20RSB = (TH2D*)hPhi2DRside->Clone("eta20RSB");
    TH2D* eta20LSB = (TH2D*)hPhi2DLside->Clone("eta20LSB");

    Float_t epsilon = 0.001;

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
    TH1D* eta20peakPhi = eta20peak->ProjectionY("eta20peakPhi", eta20peak->GetXaxis()->FindBin(-1.5 + epsilon), eta20peak->GetXaxis()->FindBin(1.5 - epsilon));
    eta20peakPhi->SetLineColor(kBlue);
    eta20peakPhi->SetStats(kFALSE);
    TH1D* eta20peakPhiNarrow = eta20peak->ProjectionY("eta20peakPhiNarrow", eta20peak->GetXaxis()->FindBin(-1.2 + epsilon), eta20peak->GetXaxis()->FindBin(1.2 - epsilon));
    eta20peakPhiNarrow->SetLineColor(kViolet);
    TH1D* eta20peakPhiNarrowest = eta20peak->ProjectionY("eta20peakPhiNarrowest", eta20peak->GetXaxis()->FindBin(-1.0 + epsilon), eta20peak->GetXaxis()->FindBin(1.0 - epsilon));
    eta20peakPhiNarrowest->SetLineColor(kRed);

    TH1D* eta20RSBEta = eta20RSB->ProjectionX("eta20RSBEta");
    eta20RSBEta->GetXaxis()->SetTitleOffset(1.0);
    eta20RSBEta->SetStats(kFALSE);
    TH1D* eta20RSBPhi = eta20RSB->ProjectionY("eta20RSBPhi", eta20RSB->GetXaxis()->FindBin(-1.5 + epsilon), eta20RSB->GetXaxis()->FindBin(1.5 - epsilon));
    eta20RSBPhi->SetLineColor(kBlue);
    eta20RSBPhi->SetStats(kFALSE);
    TH1D* eta20RSBPhiNarrow = eta20RSB->ProjectionY("eta20RSBPhiNarrow", eta20RSB->GetXaxis()->FindBin(-1.2 + epsilon), eta20RSB->GetXaxis()->FindBin(1.2 - epsilon));
    eta20RSBPhiNarrow->SetLineColor(kViolet);
    TH1D* eta20RSBPhiNarrowest = eta20RSB->ProjectionY("eta20RSBPhiNarrowest", eta20RSB->GetXaxis()->FindBin(-1.0 + epsilon), eta20RSB->GetXaxis()->FindBin(1.0 - epsilon));
    eta20RSBPhiNarrowest->SetLineColor(kRed);

    TH1D* eta20LSBEta = eta20LSB->ProjectionX("eta20LSBEta");
    TH1D* eta20LSBPhi = eta20LSB->ProjectionY("eta20LSBPhi", eta20LSB->GetXaxis()->FindBin(-1.5 + epsilon), eta20LSB->GetXaxis()->FindBin(1.5 - epsilon));
    eta20LSBPhi->SetLineColor(kBlue);
    eta20LSBPhi->SetStats(kFALSE);
    TH1D* eta20LSBPhiNarrow = eta20LSB->ProjectionY("eta20LSBPhiNarrow", eta20LSB->GetXaxis()->FindBin(-1.2 + epsilon), eta20LSB->GetXaxis()->FindBin(1.2 - epsilon));
    eta20LSBPhiNarrow->SetLineColor(kViolet);
    TH1D* eta20LSBPhiNarrowest = eta20LSB->ProjectionY("eta20LSBPhiNarrowest", eta20LSB->GetXaxis()->FindBin(-1.0 + epsilon), eta20LSB->GetXaxis()->FindBin(1.0 - epsilon));
    eta20LSBPhiNarrowest->SetLineColor(kRed);


    //reset eta range to narrow view for 2D plotting and rebin
//    eta20peak->Rebin2D(2,2);
//    eta20RSB->Rebin2D(2,2);
//    eta20LSB->Rebin2D(2,2);
    eta20peak->GetXaxis()->SetRangeUser(-1.2 + epsilon, 1.2 - epsilon);
    eta20RSB->GetXaxis()->SetRangeUser(-1.2 + epsilon, 1.2 - epsilon);
    eta20LSB->GetXaxis()->SetRangeUser(-1.2 + epsilon, 1.2 - epsilon);


    TH2D* hKK2Dpeak = (TH2D*)eta20File->Get("hKK2Dpeak");
    TH2D* hKK2DRside = (TH2D*)eta20File->Get("hKK2DRside");
    TH2D* hKK2DLside = (TH2D*)eta20File->Get("hKK2DLside");

    TH2D* LSeta20peak = (TH2D*)hKK2Dpeak->Clone("LSeta20peak");
    TH2D* LSeta20RSB = (TH2D*)hKK2DRside->Clone("LSeta20RSB");
    TH2D* LSeta20LSB = (TH2D*)hKK2DLside->Clone("LSeta20LSB");

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

//    LSeta20peak->Rebin2D(2,2);
//    LSeta20RSB->Rebin2D(2,2);
//    LSeta20LSB->Rebin2D(2,2);
    LSeta20peak->GetXaxis()->SetRangeUser(-1.2 + epsilon, 1.2);
    LSeta20RSB->GetXaxis()->SetRangeUser(-1.2 + epsilon, 1.2);
    LSeta20LSB->GetXaxis()->SetRangeUser(-1.2 + epsilon, 1.2);




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

    TPaveText *data = new TPaveText(0.593, 0.805, 0.992, 0.990, "NDC");
    data->AddText("ALICE Preliminary");
    data->AddText("p-Pb #sqrt{s_{NN}} = 5.02 TeV");
    data->AddText("0-20% Multiplicity Class (V0A)");
    data->SetFillStyle(0);
    data->SetBorderSize(0);
    data->SetTextFont(42);

    TPaveText *othertextpeak = new TPaveText(0.009, 0.805, 0.333, 0.990, "NDC");
    othertextpeak->AddText("Mass Peak Region (Signal + BG)");
    othertextpeak->AddText("4.0 < #it{p}^{h}_{T,trig} < 8.0 GeV/#it{c}");
    othertextpeak->AddText("2.0 < #it{p}^{#phi}_{T,assoc} < 4.0 GeV/#it{c}");
    othertextpeak->SetFillStyle(0);
    othertextpeak->SetBorderSize(0);
    othertextpeak->SetTextFont(42);
    
    TPaveText *othertextLSB = new TPaveText(0.009, 0.805, 0.333, 0.990, "NDC");
    othertextLSB->AddText("Left Sideband");
    othertextLSB->AddText("4.0 < #it{p}^{h}_{T,trig} < 8.0 GeV/#it{c}");
    othertextLSB->AddText("2.0 < #it{p}^{#phi}_{T,assoc} < 4.0 GeV/#it{c}");
    othertextLSB->SetFillStyle(0);
    othertextLSB->SetBorderSize(0);
    othertextLSB->SetTextFont(42);

    TPaveText *othertextRSB = new TPaveText(0.009, 0.805, 0.333, 0.990, "NDC");
    othertextRSB->AddText("Right Sideband");
    othertextRSB->AddText("4.0 < #it{p}^{h}_{T,trig} < 8.0 GeV/#it{c}");
    othertextRSB->AddText("2.0 < #it{p}^{#phi}_{T,assoc} < 4.0 GeV/#it{c}");
    othertextRSB->SetFillStyle(0);
    othertextRSB->SetBorderSize(0);
    othertextRSB->SetTextFont(42);

    TPaveText *othertext = new TPaveText(0.009, 0.805, 0.333, 0.990, "NDC");
    othertext->AddText("4.0 < #it{p}^{h}_{T,trig} < 8.0 GeV/#it{c}");
    othertext->AddText("2.0 < #it{p}^{#phi}_{T,assoc} < 4.0 GeV/#it{c}");
    othertext->SetFillStyle(0);
    othertext->SetBorderSize(0);
    othertext->SetTextFont(42);

   
    TCanvas* ccorrUSpeak = new TCanvas("ccorrUSpeak", "ccorrUSpeak", 50, 50, 800, 800);
    ccorrUSpeak->cd()->SetTheta(16.2);
    ccorrUSpeak->cd()->SetPhi(-65);
    ccorrUSpeak->cd()->SetMargin(0.15, 0.05, 0.12, 0.18);
    eta20peak->GetZaxis()->SetTitle("1/#it{N}_{trig} d^{2}#it{N}_{assoc}/d#it{#Delta#phi}d#it{#Delta#eta}");
    eta20peak->GetZaxis()->SetMaxDigits(2);
    eta20peak->GetZaxis()->SetTitleOffset(1.2);
    eta20peak->GetZaxis()->SetTitleSize(0.04);
    eta20peak->Scale(1.0/(eta20peak->GetXaxis()->GetBinWidth(1)*eta20peak->GetYaxis()->GetBinWidth(1)*0.49*0.82));
    eta20peak->GetXaxis()->CenterTitle();
    eta20peak->GetYaxis()->CenterTitle();
    eta20peak->Draw("SURF1");
    data->Draw();
    othertextpeak->Draw();

    TCanvas* ccorrUSRside = new TCanvas("ccorrUSRside", "ccorrUSRside", 50, 50, 800, 800);
    ccorrUSRside->cd()->SetTheta(16.2);
    ccorrUSRside->cd()->SetPhi(-65);
    ccorrUSRside->cd()->SetMargin(0.15, 0.05, 0.12, 0.18);
    eta20RSB->GetZaxis()->SetTitle("1/#it{N}_{trig} d^{2}#it{N}_{assoc}/d#it{#Delta#phi}d#it{#Delta#eta}");
    eta20RSB->GetZaxis()->SetMaxDigits(2);
    eta20RSB->GetZaxis()->SetTitleOffset(1.2);
    eta20RSB->GetZaxis()->SetTitleSize(0.04);
    eta20RSB->Scale(1.0/(eta20RSB->GetXaxis()->GetBinWidth(1)*eta20RSB->GetYaxis()->GetBinWidth(1)*0.49*0.82));
    eta20RSB->GetXaxis()->CenterTitle();
    eta20RSB->GetYaxis()->CenterTitle();
    eta20RSB->Draw("SURF1");
    data->Draw();
    othertextRSB->Draw();

    TCanvas* ccorrUSLside = new TCanvas("ccorrUSLside", "ccorrUSLside", 50, 50, 800, 800);
    ccorrUSLside->cd()->SetTheta(16.2);
    ccorrUSLside->cd()->SetPhi(-65);
    ccorrUSLside->cd()->SetMargin(0.15, 0.05, 0.12, 0.18);
    eta20LSB->GetZaxis()->SetTitle("1/#it{N}_{trig} d^{2}#it{N}_{assoc}/d#it{#Delta#phi}d#it{#Delta#eta}");
    eta20LSB->GetZaxis()->SetMaxDigits(2);
    eta20LSB->GetZaxis()->SetTitleOffset(1.2);
    eta20LSB->GetZaxis()->SetTitleSize(0.04);
    eta20LSB->Scale(1.0/(eta20LSB->GetXaxis()->GetBinWidth(1)*eta20LSB->GetYaxis()->GetBinWidth(1)*0.49*0.82));
    eta20LSB->GetXaxis()->CenterTitle();
    eta20LSB->GetYaxis()->CenterTitle();
    eta20LSB->Draw("SURF1");
    data->Draw();
    othertextLSB->Draw();

    TCanvas* ccorrLSpeak = new TCanvas("ccorrLSpeak", "ccorrLSpeak", 50, 50, 800, 800);
    ccorrLSpeak->cd()->SetTheta(50);
    ccorrLSpeak->cd()->SetPhi(50);
    LSeta20peak->Draw("SURF1");

    TCanvas* ccorrLSRside = new TCanvas("ccorrLSRside", "ccorrLSRside", 50, 50, 800, 800);
    ccorrLSRside->cd()->SetTheta(50);
    ccorrLSRside->cd()->SetPhi(50);
    LSeta20RSB->Draw("SURF1");

    TCanvas* ccorrLSLside = new TCanvas("ccorrLSLside", "ccorrLSLside", 50, 50, 800, 800);
    ccorrLSLside->cd()->SetTheta(50);
    ccorrLSLside->cd()->SetPhi(50);
    LSeta20LSB->Draw("SURF1");

    //plotting sideband distributions
    TH1D* dphiLSB = (TH1D*)eta20LSBPhiNarrow->Clone("dphiLSB");
    dphiLSB->Scale(1.0/(dphiLSB->Integral()*dphiLSB->GetXaxis()->GetBinWidth(1)));
    TH1D* dphiRSB = (TH1D*)eta20RSBPhiNarrow->Clone("dphiRSB");
    dphiRSB->Scale(1.0/(dphiRSB->Integral()*dphiRSB->GetXaxis()->GetBinWidth(1)));

    TLegend *sblegend = new TLegend(0.4, 0.6, 0.7, 0.9);
    sblegend->AddEntry(dphiLSB, "Left Sideband Region", "le");
    sblegend->AddEntry(dphiRSB, "Right Sideband Region", "le");
    sblegend->SetLineWidth(0);

    TCanvas* cSBcompare = new TCanvas("cSBcompare", "cSBcompare", 50, 50, 600, 600);
    cSBcompare->cd()->SetMargin(0.12, 0.04, 0.1, 0.04);;
    dphiLSB->SetLineWidth(2);
    dphiLSB->SetLineColor(kViolet+1);
    dphiLSB->GetXaxis()->SetTitle("#Delta #varphi");
    dphiLSB->GetYaxis()->SetTitle("1/#it{N}_{entries} d#it{N}_{assoc}/d#it{#Delta#varphi}");
    dphiLSB->GetYaxis()->SetMaxDigits(3);
    dphiLSB->SetStats(kFALSE);
    dphiLSB->GetXaxis()->SetTitleOffset(0.8);
    dphiLSB->Draw("H");
    dphiRSB->SetLineWidth(2);
    dphiRSB->SetLineColor(kRed+1);
    dphiRSB->Draw("H SAME");
    sblegend->Draw();
    data->Draw();
    othertext->Draw();

    //plot mix corrected LSB region (delta phi)
    TPaveText* LSBText = new TPaveText(0.51, 0.67, 0.88, 0.87, "NDC");
    LSBText->SetBorderSize(0);
    LSBText->SetFillColor(kWhite);
    LSBText->AddText("Left Sideband");

    TCanvas *cLSB = new TCanvas("cLSB", "cLSB", 50, 50, 600, 600);
    cLSB->cd();
    eta20LSBPhiNarrow->GetYaxis()->SetMaxDigits(3);
    eta20LSBPhiNarrow->GetYaxis()->SetTitle("Per-Trigger Pair Yield");
    eta20LSBPhiNarrow->SetStats(kFALSE);
    eta20LSBPhiNarrow->SetMarkerSize(2);
    eta20LSBPhiNarrow->SetMarkerStyle(33);
    eta20LSBPhiNarrow->GetXaxis()->SetTitleOffset(0.8);
    eta20LSBPhiNarrow->Draw("H");
    LSBText->Draw();

    //plot mix corrected RSB region (delta phi)
    TPaveText* RSBText = new TPaveText(0.51, 0.67, 0.88, 0.87, "NDC");
    RSBText->SetBorderSize(0);
    RSBText->SetFillColor(kWhite);
    RSBText->AddText("Right Sideband");

    TCanvas *cRSB = new TCanvas("cRSB", "cRSB", 50, 50, 600, 600);
    cRSB->cd();
    eta20RSBPhiNarrow->GetYaxis()->SetMaxDigits(3);
    eta20RSBPhiNarrow->GetYaxis()->SetTitle("Per-Trigger Pair Yield");
    eta20RSBPhiNarrow->SetStats(kFALSE);
    eta20RSBPhiNarrow->SetMarkerSize(2);
    eta20RSBPhiNarrow->SetMarkerStyle(33);
    eta20RSBPhiNarrow->GetXaxis()->SetTitleOffset(0.8);
    eta20RSBPhiNarrow->Draw("H");
    RSBText->Draw();

    //plot mix corrected mass peak region (delta phi)
    TPaveText* peakText = new TPaveText(0.51, 0.67, 0.88, 0.87, "NDC");
    peakText->SetBorderSize(0);
    peakText->SetFillColor(kWhite);
    peakText->AddText("Mass Peak Region");

    TCanvas *cpeak = new TCanvas("cpeak", "cpeak", 50, 50, 600, 600);
    cpeak->cd();
    eta20peakPhiNarrow->GetYaxis()->SetMaxDigits(3);
    eta20peakPhiNarrow->GetYaxis()->SetTitle("Per-Trigger Pair Yield");
    eta20peakPhiNarrow->SetStats(kFALSE);
    eta20peakPhiNarrow->SetMarkerSize(2);
    eta20peakPhiNarrow->SetMarkerStyle(33);
    eta20peakPhiNarrow->GetXaxis()->SetTitleOffset(0.8);
    eta20peakPhiNarrow->Draw("H");
    peakText->Draw();


}
