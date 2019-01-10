void plotMixCorrected(string inputfile){
    TFile* eta20File = new TFile(inputfile.c_str());

    TH2D* hPhi2Dpeak = (TH2D*)eta20File->Get("hPhi2Dpeak");
    TH2D* hPhi2DRside = (TH2D*)eta20File->Get("hPhi2DRside");
    TH2D* hPhi2DLside = (TH2D*)eta20File->Get("hPhi2DLside");

    TH2D* eta20peak = (TH2D*)hPhi2Dpeak->Clone("eta20peak");
    TH2D* eta20RSB = (TH2D*)hPhi2DRside->Clone("eta20RSB");
    TH2D* eta20LSB = (TH2D*)hPhi2DLside->Clone("eta20LSB");

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
//    eta20peak->Rebin2D(2,2);
//    eta20RSB->Rebin2D(2,2);
//    eta20LSB->Rebin2D(2,2);
    eta20peak->GetXaxis()->SetRangeUser(-1.2, 1.2);
    eta20RSB->GetXaxis()->SetRangeUser(-1.2, 1.2);
    eta20LSB->GetXaxis()->SetRangeUser(-1.2, 1.2);


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
