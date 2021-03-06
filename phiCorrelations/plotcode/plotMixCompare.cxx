void plotMixCompare(){

    TFile* simplecorrFile = new TFile("~/phiStudies/LHC16q_FAST_V0A_20170925/trig_4_8_assoc_2_4_simplecorr_phiCorrelations_mult_0_20.root");
    TH2D* eta20peak = hPhi2Dpeak->Clone("eta20peak");
    TH2D* eta20RSB = hPhi2DRside->Clone("eta20RSB");
    TH2D* eta20LSB = hPhi2DLside->Clone("eta20LSB");
    TH2D* LSeta20peak = hKK2Dpeak->Clone("LSeta20peak");
    TH2D* LSeta20RSB = hKK2DRside->Clone("LSeta20RSB");
    TH2D* LSeta20LSB = hKK2DLside->Clone("LSeta20LSB");

    //(Simple - Mix)/Mix
    TFile* mixcorrFile = new TFile("~/phiStudies/LHC16q_FAST_V0A_20170925/trig_4_8_assoc_2_4_mixcorr_phiCorrelations_mult_0_20.root");
    eta20peak->Add(hPhi2Dpeak, -1.0);
    eta20RSB->Add(hPhi2DRside, -1.0);
    eta20LSB->Add(hPhi2DLside, -1.0);
    LSeta20peak->Add(hKK2Dpeak, -1.0);
    LSeta20RSB->Add(hKK2DRside, -1.0);
    LSeta20LSB->Add(hKK2DLside, -1.0);

    eta20peak->Divide(hPhi2Dpeak);
    eta20RSB->Divide(hPhi2DRside);
    eta20LSB->Divide(hPhi2DLside);
    LSeta20peak->Divide(hKK2Dpeak);
    LSeta20RSB->Divide(hKK2DRside);
    LSeta20LSB->Divide(hKK2DLside);

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

    TH1D* eta20peakEta = eta20peak->ProjectionX("eta20peakEta", 1, eta20peak->GetYaxis()->GetNbins());
    //eta20peakEta->Scale(1.0/16.0);
    eta20peakEta->GetXaxis()->SetTitleOffset(1.0);
    eta20peakEta->SetStats(kFALSE);
    TH1D* eta20peakPhi = eta20peak->ProjectionY("eta20peakPhi", eta20peak->GetXaxis()->FindBin(-1.5), eta20peak->GetXaxis()->FindBin(1.5));
    //eta20peakPhi->Scale(1.0/(float(eta20peak->GetXaxis()->FindBin(1.5) - eta20peak->GetXaxis()->FindBin(-1.5) + 1)));
    eta20peakPhi->SetLineColor(kBlue);
    eta20peakPhi->SetStats(kFALSE);
    TH1D* eta20peakPhiNarrow = eta20peak->ProjectionY("eta20peakPhiNarrow", eta20peak->GetXaxis()->FindBin(-1.2), eta20peak->GetXaxis()->FindBin(1.2));
    //eta20peakPhiNarrow->Scale(1.0/(float(eta20peak->GetXaxis()->FindBin(1.2) - eta20peak->GetXaxis()->FindBin(-1.2) + 1)));
    eta20peakPhiNarrow->SetLineColor(kViolet);
    TH1D* eta20peakPhiNarrowest = eta20peak->ProjectionY("eta20peakPhiNarrowest", eta20peak->GetXaxis()->FindBin(-1.0), eta20peak->GetXaxis()->FindBin(1.0));
    //eta20peakPhiNarrowest->Scale(1.0/(float(eta20peak->GetXaxis()->FindBin(1.1) - eta20peak->GetXaxis()->FindBin(-1.1) + 1)));
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


    eta20peak->GetXaxis()->SetRangeUser(-1.2, 1.2);
    eta20RSB->GetXaxis()->SetRangeUser(-1.2, 1.2);
    eta20LSB->GetXaxis()->SetRangeUser(-1.2, 1.2);

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
    //eta20peakPhi->Draw("H");
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

   
    TCanvas* ccorrUSpeak = new TCanvas("ccorrUSpeak", "ccorrUSpeak", 50, 50, 800, 800);
    ccorrUSpeak->cd()->SetTheta(50);
    ccorrUSpeak->cd()->SetPhi(50);
    eta20peak->Draw("SURF1");

    TCanvas* ccorrUSRside = new TCanvas("ccorrUSRside", "ccorrUSRside", 50, 50, 800, 800);
    ccorrUSRside->cd()->SetTheta(50);
    ccorrUSRside->cd()->SetPhi(50);
    eta20RSB->Draw("SURF1");

    TCanvas* ccorrUSLside = new TCanvas("ccorrUSLside", "ccorrUSLside", 50, 50, 800, 800);
    ccorrUSLside->cd()->SetTheta(50);
    ccorrUSLside->cd()->SetPhi(50);
    eta20LSB->Draw("SURF1");

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


}
