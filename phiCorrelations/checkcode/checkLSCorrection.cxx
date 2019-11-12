void checkLSCorrection(string inputFile){
    TFile* input = new TFile(inputFile.c_str());

    TH2D* hPhi2DLside = (TH2D*)input->Get("hPhi2DLside");
    hPhi2DLside->GetXaxis()->SetRangeUser(-1.2, 1.2);
    TH2D* hPhi2DRside = (TH2D*)input->Get("hPhi2DRside");
    hPhi2DRside->GetXaxis()->SetRangeUser(-1.2, 1.2);
    TH2D* hKK2Dpeak = (TH2D*)input->Get("hKK2Dpeak");
    hKK2Dpeak->GetXaxis()->SetRangeUser(-1.2, 1.2);
    TH2D* hKK2DLside = (TH2D*)input->Get("hKK2DLside");
    hKK2DLside->GetXaxis()->SetRangeUser(-1.2, 1.2);
    TH2D* hKK2DRside = (TH2D*)input->Get("hKK2DRside");
    hKK2DRside->GetXaxis()->SetRangeUser(-1.2, 1.2);

    //Normalize everything in the smaller delta-eta range to make comparisons easier
    hKK2Dpeak->Scale(1.0/hKK2Dpeak->Integral(hKK2Dpeak->GetXaxis()->FindBin(-1.2), hKK2Dpeak->GetXaxis()->FindBin(1.2), 1, hKK2Dpeak->GetYaxis()->GetNbins()));
    hKK2DLside->Scale(1.0/hKK2DLside->Integral(hKK2DLside->GetXaxis()->FindBin(-1.2), hKK2DLside->GetXaxis()->FindBin(1.2), 1, hKK2DLside->GetYaxis()->GetNbins()));
    hKK2DRside->Scale(1.0/hKK2DRside->Integral(hKK2DRside->GetXaxis()->FindBin(-1.2), hKK2DRside->GetXaxis()->FindBin(1.2), 1, hKK2DRside->GetYaxis()->GetNbins()));
    hPhi2DLside->Scale(1.0/hPhi2DLside->Integral(hPhi2DLside->GetXaxis()->FindBin(-1.2), hPhi2DLside->GetXaxis()->FindBin(1.2), 1, hPhi2DLside->GetYaxis()->GetNbins()));
    hPhi2DRside->Scale(1.0/hPhi2DRside->Integral(hPhi2DRside->GetXaxis()->FindBin(-1.2), hPhi2DRside->GetXaxis()->FindBin(1.2), 1, hPhi2DRside->GetYaxis()->GetNbins()));

    //compute average of left and right side to compare to LS peak region
    TH2D* avgLSsideband = (TH2D*)hKK2DLside->Clone("avgLSsideband");
    avgLSsideband->Add(hKK2DRside);
    avgLSsideband->Scale(0.5);

    TH2D* avgUSsideband = (TH2D*)hPhi2DLside->Clone("avgUSsideband");
    avgUSsideband->Add(hPhi2DRside);
    avgUSsideband->Scale(0.5);

    //compute residuals for the different scenarios (US and LS sidebands)
    TH2D* resUSLeftRight = (TH2D*)hPhi2DRside->Clone("resUSLeftRight");
    resUSLeftRight->Add(hPhi2DLside, -1.0);
    resUSLeftRight->Divide(hPhi2DRside);

    TH2D* resLSLeftRight = (TH2D*)hKK2DRside->Clone("resLSLeftRight");
    resLSLeftRight->Add(hKK2DLside, -1.0);
    resLSLeftRight->Divide(hKK2DRside);

    TH2D* resLSPeakAvg = (TH2D*)hKK2Dpeak->Clone("resLSPeakAvg");
    resLSPeakAvg->Add(avgLSsideband, -1.0);
    resLSPeakAvg->Divide(hKK2Dpeak);

    TH2D* resUSPeakAvg = (TH2D*)hKK2Dpeak->Clone("resUSPeakAvg");
    resUSPeakAvg->Add(avgUSsideband, -1.0);
    resUSPeakAvg->Divide(hKK2Dpeak);

    //set up the 1D projections for each of these 2D plots
    TH1D* hKK2DLsidePhi = (TH1D*)hKK2DLside->ProjectionY("hKK2DLsidePhi", hKK2DLside->GetXaxis()->FindBin(-1.2), hKK2DLside->GetXaxis()->FindBin(1.2));
    TH1D* hKK2DLsideEta = (TH1D*)hKK2DLside->ProjectionX("hKK2DLsideEta", 1, hKK2DLside->GetYaxis()->GetNbins());

    TH1D* hKK2DRsidePhi = (TH1D*)hKK2DRside->ProjectionY("hKK2DRsidePhi", hKK2DRside->GetXaxis()->FindBin(-1.2), hKK2DRside->GetXaxis()->FindBin(1.2));
    TH1D* hKK2DRsideEta = (TH1D*)hKK2DRside->ProjectionX("hKK2DRsideEta", 1, hKK2DRside->GetYaxis()->GetNbins());

    TH1D* hPhi2DLsidePhi = (TH1D*)hPhi2DLside->ProjectionY("hPhi2DLsidePhi", hPhi2DLside->GetXaxis()->FindBin(-1.2), hPhi2DLside->GetXaxis()->FindBin(1.2));
    TH1D* hPhi2DLsideEta = (TH1D*)hPhi2DLside->ProjectionX("hPhi2DLsideEta", 1, hPhi2DLside->GetYaxis()->GetNbins());

    TH1D* hPhi2DRsidePhi = (TH1D*)hPhi2DRside->ProjectionY("hPhi2DRsidePhi", hPhi2DRside->GetXaxis()->FindBin(-1.2), hPhi2DRside->GetXaxis()->FindBin(1.2));
    TH1D* hPhi2DRsideEta = (TH1D*)hPhi2DRside->ProjectionX("hPhi2DRsideEta", 1, hPhi2DRside->GetYaxis()->GetNbins());

    TH1D* hKK2DpeakPhi = (TH1D*)hKK2Dpeak->ProjectionY("hKK2DpeakPhi", hKK2Dpeak->GetXaxis()->FindBin(-1.2), hKK2Dpeak->GetXaxis()->FindBin(1.2));
    TH1D* hKK2DpeakEta = (TH1D*)hKK2Dpeak->ProjectionX("hKK2DpeakEta", 1, hKK2Dpeak->GetYaxis()->GetNbins());

    TH1D* avgLSsidebandPhi = (TH1D*)avgLSsideband->ProjectionY("avgLSsidebandPhi", avgLSsideband->GetXaxis()->FindBin(-1.2), avgLSsideband->GetXaxis()->FindBin(1.2));
    TH1D* avgLSsidebandEta = (TH1D*)avgLSsideband->ProjectionX("avgLSsidebandEta", 1, avgLSsideband->GetYaxis()->GetNbins());

    TH1D* avgUSsidebandPhi = (TH1D*)avgUSsideband->ProjectionY("avgUSsidebandPhi", avgUSsideband->GetXaxis()->FindBin(-1.2), avgUSsideband->GetXaxis()->FindBin(1.2));
    TH1D* avgUSsidebandEta = (TH1D*)avgUSsideband->ProjectionX("avgUSsidebandEta", 1, avgUSsideband->GetYaxis()->GetNbins());


/*    
    TH1D* resLSLeftRightPhi = (TH1D*)resLSLeftRight->ProjectionY("resLSLeftRightPhi", resLSLeftRight->GetXaxis()->FindBin(-1.2), resLSLeftRight->GetXaxis()->FindBin(1.2));
    TH1D* resLSLeftRightEta = (TH1D*)resLSLeftRight->ProjectionX("resLSLeftRightEta", 1, resLSLeftRight->GetYaxis()->GetNbins());

    TH1D* resLSPeakAvgPhi = (TH1D*)resLSPeakAvg->ProjectionY("resLSPeakAvgPhi", resLSPeakAvg->GetXaxis()->FindBin(-1.2), resLSPeakAvg->GetXaxis()->FindBin(1.2));
    TH1D* resLSPeakAvgEta = (TH1D*)resLSPeakAvg->ProjectionX("resLSPeakAvgEta", 1, resLSPeakAvg->GetYaxis()->GetNbins());
*/
    TH1D* resLSLeftRightPhi = (TH1D*)hKK2DRsidePhi->Clone("resLSLeftRightPhi");
    resLSLeftRightPhi->Add(hKK2DLsidePhi, -1.0);
    resLSLeftRightPhi->Divide(hKK2DRsidePhi);

    TH1D* resLSLeftRightEta = (TH1D*)hKK2DRsideEta->Clone("resLSLeftRightEta");
    resLSLeftRightEta->Add(hKK2DLsideEta, -1.0);
    resLSLeftRightEta->Divide(hKK2DRsideEta);

    TH1D* resLSPeakAvgPhi = (TH1D*)hKK2DpeakPhi->Clone("resLSPeakAvgPhi");
    resLSPeakAvgPhi->Add(avgLSsidebandPhi, -1.0);
    resLSPeakAvgPhi->Divide(hKK2DpeakPhi);

    TH1D* resLSPeakAvgEta = (TH1D*)hKK2DpeakEta->Clone("resLSPeakAvgEta");
    resLSPeakAvgEta->Add(avgLSsidebandEta, -1.0);
    resLSPeakAvgEta->Divide(hKK2DpeakEta);

    TH1D* resUSLeftRightPhi = (TH1D*)hPhi2DRsidePhi->Clone("resUSLeftRightPhi");
    resUSLeftRightPhi->Add(hPhi2DLsidePhi, -1.0);
    resUSLeftRightPhi->Divide(hPhi2DRsidePhi);

    TH1D* resUSLeftRightEta = (TH1D*)hPhi2DRsideEta->Clone("resUSLeftRightEta");
    resUSLeftRightEta->Add(hPhi2DLsideEta, -1.0);
    resUSLeftRightEta->Divide(hPhi2DRsideEta);

    TH1D* resUSPeakAvgPhi = (TH1D*)hKK2DpeakPhi->Clone("resUSPeakAvgPhi");
    resUSPeakAvgPhi->Add(avgUSsidebandPhi, -1.0);
    resUSPeakAvgPhi->Divide(hKK2DpeakPhi);

    TH1D* resUSPeakAvgEta = (TH1D*)hKK2DpeakEta->Clone("resUSPeakAvgEta");
    resUSPeakAvgEta->Add(avgUSsidebandEta, -1.0);
    resUSPeakAvgEta->Divide(hKK2DpeakEta);


    //Draw
    
    TCanvas* cLS = new TCanvas("cLS", "cLS", 50, 50, 1200, 1200);
    cLS->Divide(3,3);
    cLS->cd(1);
    hKK2DLside->SetTitle("h-KK LS correlation for Left Sideband");
    hKK2DLside->SetStats(0);
    hKK2DLside->Draw("SURF1");
    gPad->SetTheta(22.5);
    gPad->SetPhi(-76);
    gPad->Update();
    cLS->cd(2);
    resLSLeftRight->SetTitle("Residual for LS sidebands: (Right - Left)/Right");
    resLSLeftRight->SetStats(0);
    resLSLeftRight->Draw("SURF1");;
    gPad->SetTheta(22.5);
    gPad->SetPhi(-76);
    gPad->Update();
    cLS->cd(3);
    hKK2DRside->SetTitle("h-KK LS correlation for Right Sideband");
    hKK2DRside->SetStats(0);
    hKK2DRside->Draw("SURF1");;
    gPad->SetTheta(22.5);
    gPad->SetPhi(-76);
    gPad->Update();
    cLS->cd(4);
    hKK2DLsideEta->Draw("HE");
    cLS->cd(5);
    resLSLeftRightEta->Draw("HE");
    cLS->cd(6);
    hKK2DRsideEta->Draw("HE");
    cLS->cd(7);
    hKK2DLsidePhi->Draw("HE");
    cLS->cd(8);
    resLSLeftRightPhi->Draw("HE");
    cLS->cd(9);
    hKK2DRsidePhi->Draw("HE");

    TCanvas* cUS = new TCanvas("cUS", "cUS", 50, 50, 1200, 1200);
    cUS->Divide(3,3);
    cUS->cd(1);
    hPhi2DLside->SetTitle("h-Phi US correlation for Left Sideband");
    hPhi2DLside->SetStats(0);
    hPhi2DLside->Draw("SURF1");
    gPad->SetTheta(22.5);
    gPad->SetPhi(-76);
    gPad->Update();
    cUS->cd(2);
    resUSLeftRight->SetTitle("Residual for US sidebands: (Right - Left)/Right");
    resUSLeftRight->SetStats(0);
    resUSLeftRight->Draw("SURF1");;
    gPad->SetTheta(22.5);
    gPad->SetPhi(-76);
    gPad->Update();
    cUS->cd(3);
    hPhi2DRside->SetTitle("h-Phi US correlation for Right Sideband");
    hPhi2DRside->SetStats(0);
    hPhi2DRside->Draw("SURF1");;
    gPad->SetTheta(22.5);
    gPad->SetPhi(-76);
    gPad->Update();
    cUS->cd(4);
    hPhi2DLsideEta->Draw("HE");
    cUS->cd(5);
    resUSLeftRightEta->Draw("HE");
    cUS->cd(6);
    hPhi2DRsideEta->Draw("HE");
    cUS->cd(7);
    hPhi2DLsidePhi->Draw("HE");
    cUS->cd(8);
    resUSLeftRightPhi->Draw("HE");
    cUS->cd(9);
    hPhi2DRsidePhi->Draw("HE");

    TCanvas* cAvg = new TCanvas("cAvg", "cAvg", 55, 55, 1200, 1200);
    cAvg->Divide(3,3);
    cAvg->cd(1);
    hKK2Dpeak->SetTitle("h-KK LS correlation for mass peak region");
    hKK2Dpeak->SetStats(0);
    hKK2Dpeak->Draw("SURF1");
    gPad->SetTheta(22.5);
    gPad->SetPhi(-76);
    gPad->Update();
    cAvg->cd(2);
    resLSPeakAvg->SetTitle("Residual For LS avg: (Peak region - Left&Right)/Peak");
    resLSPeakAvg->SetStats(0);
    resLSPeakAvg->Draw("SURF1");;
    gPad->SetTheta(22.5);
    gPad->SetPhi(-76);
    gPad->Update();
    cAvg->cd(3);
    avgLSsideband->SetTitle("h-KK LS correlation Average between Left & Right sidebands");
    avgLSsideband->SetStats(0);
    avgLSsideband->Draw("SURF1");;
    gPad->SetTheta(22.5);
    gPad->SetPhi(-76);
    gPad->Update();
    cAvg->cd(4);
    hKK2DpeakEta->Draw("HE");
    cAvg->cd(5);
    resLSPeakAvgEta->Draw("HE");
    cAvg->cd(6);
    avgLSsidebandEta->Draw("HE");
    cAvg->cd(7);
    hKK2DpeakPhi->Draw("HE");
    cAvg->cd(8);
    resLSPeakAvgPhi->Draw("HE");
    cAvg->cd(9);
    avgLSsidebandPhi->Draw("HE");

    TCanvas* cUSAvg = new TCanvas("cUSAvg", "cUSAvg", 55, 55, 1200, 1200);
    cUSAvg->Divide(3,3);
    cUSAvg->cd(1);
    hKK2Dpeak->SetTitle("h-KK LS correlation for mass peak region");
    hKK2Dpeak->SetStats(0);
    hKK2Dpeak->Draw("SURF1");
    gPad->SetTheta(22.5);
    gPad->SetPhi(-76);
    gPad->Update();
    cUSAvg->cd(2);
    resUSPeakAvg->SetTitle("Residual For US avg: (LS Peak region - US Left&Right)/Peak");
    resUSPeakAvg->SetStats(0);
    resUSPeakAvg->Draw("SURF1");;
    gPad->SetTheta(22.5);
    gPad->SetPhi(-76);
    gPad->Update();
    cUSAvg->cd(3);
    avgUSsideband->SetTitle("h-Phi US correlation Average between Left & Right sidebands");
    avgUSsideband->SetStats(0);
    avgUSsideband->Draw("SURF1");;
    gPad->SetTheta(22.5);
    gPad->SetPhi(-76);
    gPad->Update();
    cUSAvg->cd(4);
    hKK2DpeakEta->Draw("HE");
    cUSAvg->cd(5);
    resUSPeakAvgEta->Draw("HE");
    cUSAvg->cd(6);
    avgUSsidebandEta->Draw("HE");
    cUSAvg->cd(7);
    hKK2DpeakPhi->Draw("HE");
    cUSAvg->cd(8);
    resUSPeakAvgPhi->Draw("HE");
    cUSAvg->cd(9);
    avgUSsidebandPhi->Draw("HE");

    TCanvas* cLSnoRes = new TCanvas("cLSnoRes", "cLSnoRes", 55, 55, 1200, 500);
    cLSnoRes->Divide(3,1);
    cLSnoRes->cd(1);
    hKK2DLside->SetTitle("h-KK LS correlation for Left Sideband");
    hKK2DLside->SetStats(0);
    hKK2DLside->Draw("SURF1");
    gPad->SetTheta(22.5);
    gPad->SetPhi(-76);
    gPad->Update();
    cLSnoRes->cd(2);
    hKK2Dpeak->SetTitle("h-KK LS correlation for mass peak region");
    hKK2Dpeak->SetStats(0);
    hKK2Dpeak->Draw("SURF1");
    gPad->SetTheta(22.5);
    gPad->SetPhi(-76);
    gPad->Update();
    cLSnoRes->cd(3);
    hKK2DRside->SetTitle("h-KK LS correlation for Right Sideband");
    hKK2DRside->SetStats(0);
    hKK2DRside->Draw("SURF1");;
    gPad->SetTheta(22.5);
    gPad->SetPhi(-76);
    gPad->Update();

}
