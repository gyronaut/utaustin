void makeUSCorrections(string inputFile){
    TFile* input = new TFile(inputFile.c_str());
    TH2D* hPhi2Dpeak = (TH2D*)input->Get("hPhi2Dpeak");
    TH2D* hPhi2DLside = (TH2D*)input->Get("hPhi2DLside");
    TH2D* hPhi2DRside = (TH2D*)input->Get("hPhi2DRside");
    TH2D* hKK2Dpeak = (TH2D*)input->Get("hKK2Dpeak");
    TH2D* hKK2DLside = (TH2D*)input->Get("hKK2DLside");
    TH2D* hKK2DRside = (TH2D*)input->Get("hKK2DRside");

    TH2D* trigDistSameUS = (TH2D*)input->Get("fTrigSameUSDist");
    TH2D* trigDistSameLS = (TH2D*)input->Get("fTrigSameLSDist");

    hPhi2Dpeak->SetName("uncorrectedhPhi2Dpeak");

    Float_t totalTrigUS = trigDistSameUS->Integral(trigDistSameUS->GetXaxis()->FindBin(4.0), trigDistSameUS->GetXaxis()->FindBin(10.0));
    Float_t totalTrigLS = trigDistSameLS->Integral(trigDistSameLS->GetXaxis()->FindBin(4.0), trigDistSameLS->GetXaxis()->FindBin(10.0));
    hPhi2Dpeak->Scale(1.0/totalTrigUS);
    hPhi2DLside->Scale(1.0/totalTrigUS);
    hPhi2DRside->Scale(1.0/totalTrigUS);
    hKK2Dpeak->Scale(1.0/totalTrigLS);
    hKK2DLside->Scale(1.0/totalTrigLS);
    hKK2DRside->Scale(1.0/totalTrigLS);

    //Using US sideband regions to estimate the BG under the peak region
    //Now using the (unweighted) average of the Left and Right side sideband (correct? needs additional checks?)
    TH2D* hPhiBGPeakRegionL = (TH2D*)hPhi2DLside->Clone("hPhiBGPeakRegionL");
    hPhiBGPeakRegionL->Scale(1.0/(hPhi2DLside->Integral(hPhi2DLside->GetXaxis()->FindBin(-1.2), hPhi2DLside->GetXaxis()->FindBin(1.2), 1, hPhi2DLside->GetYaxis()->GetNbins())));
    hPhiBGPeakRegionL->GetXaxis()->SetRangeUser(-1.2, 1.2);
    TH1D* hPhiBGPeakRegionL_deta = (TH1D*)hPhiBGPeakRegionL->ProjectionX("hPhiBGPeakRegionL_deta", 1, hPhiBGPeakRegionL->GetYaxis()->GetNbins());
    TH1D* hPhiBGPeakRegionL_dphi = (TH1D*)hPhiBGPeakRegionL->ProjectionY("hPhiBGPeakRegionL_dphi", hPhiBGPeakRegionL->GetXaxis()->FindBin(-1.2), hPhiBGPeakRegionL->GetXaxis()->FindBin(1.2));

    TH2D* hPhiBGPeakRegionR = (TH2D*)hPhi2DRside->Clone("hPhiBGPeakRegionR");
    hPhiBGPeakRegionR->Scale(1.0/(hPhi2DRside->Integral(hPhi2DRside->GetXaxis()->FindBin(-1.2), hPhi2DRside->GetXaxis()->FindBin(1.2), 1, hPhi2DRside->GetYaxis()->GetNbins())));
    hPhiBGPeakRegionR->GetXaxis()->SetRangeUser(-1.2, 1.2);
    TH1D* hPhiBGPeakRegionR_deta = (TH1D*)hPhiBGPeakRegionR->ProjectionX("hPhiBGPeakRegionR_deta", 1, hPhiBGPeakRegionR->GetYaxis()->GetNbins());
    TH1D* hPhiBGPeakRegionR_dphi = (TH1D*)hPhiBGPeakRegionR->ProjectionY("hPhiBGPeakRegionR_dphi", hPhiBGPeakRegionR->GetXaxis()->FindBin(-1.2), hPhiBGPeakRegionR->GetXaxis()->FindBin(1.2));

    TH2D* hPhiBGPeakRegion = (TH2D*)hPhiBGPeakRegionL->Clone("hPhiBGPeakregion");
    hPhiBGPeakRegion->Add(hPhiBGPeakRegionR);
    hPhiBGPeakRegion->Scale(0.5);
    hPhiBGPeakRegion->GetXaxis()->SetRangeUser(-1.2, 1.2);
    TH1D* hPhiBGPeakRegion_deta = (TH1D*)hPhiBGPeakRegion->ProjectionX("hPhiBGPeakRegion_deta", 1, hPhiBGPeakRegion->GetYaxis()->GetNbins());
    TH1D* hPhiBGPeakRegion_dphi = (TH1D*)hPhiBGPeakRegion->ProjectionY("hPhiBGPeakRegion_dphi", hPhiBGPeakRegion->GetXaxis()->FindBin(-1.2), hPhiBGPeakRegion->GetXaxis()->FindBin(1.2));


    //US residual checks between SB average and the Left and Right separately
    TH2D* resLeftVsAvg = (TH2D*)hPhiBGPeakRegionL->Clone("resLeftVsAbg");
    resLeftVsAvg->Add(hPhiBGPeakRegion, -1.0);
    resLeftVsAvg->Divide(hPhiBGPeakRegionL);
    resLeftVsAvg->GetXaxis()->SetRangeUser(-1.2, 1.2);
    TH1D* resLeftVsAvg_deta = (TH1D*)hPhiBGPeakRegionL_deta->Clone("resLeftVsAvg_deta");
    resLeftVsAvg_deta->Add(hPhiBGPeakRegion_deta, -1.0);
    resLeftVsAvg_deta->Divide(hPhiBGPeakRegionL_deta);
    resLeftVsAvg_deta->GetXaxis()->SetRangeUser(-1.2, 1.2);
    TH1D* resLeftVsAvg_dphi = (TH1D*)hPhiBGPeakRegionL_dphi->Clone("resLeftVsAvg_dphi");
    resLeftVsAvg_dphi->Add(hPhiBGPeakRegion_dphi, -1.0);
    resLeftVsAvg_dphi->Divide(hPhiBGPeakRegionL_dphi);

    TH2D* resRightVsAvg = (TH2D*)hPhiBGPeakRegionR->Clone("resRightVsAbg");
    resRightVsAvg->Add(hPhiBGPeakRegion, -1.0);
    resRightVsAvg->Divide(hPhiBGPeakRegionR);
    resRightVsAvg->GetXaxis()->SetRangeUser(-1.2, 1.2);
    TH1D* resRightVsAvg_deta = (TH1D*)hPhiBGPeakRegionR_deta->Clone("resRightVsAvg_deta");
    resRightVsAvg_deta->Add(hPhiBGPeakRegion_deta, -1.0);
    resRightVsAvg_deta->Divide(hPhiBGPeakRegionR_deta);
    resRightVsAvg_deta->GetXaxis()->SetRangeUser(-1.2, 1.2);
    TH1D* resRightVsAvg_dphi = (TH1D*)hPhiBGPeakRegionR_dphi->Clone("resRightVsAvg_dphi");
    resRightVsAvg_dphi->Add(hPhiBGPeakRegion_dphi, -1.0);
    resRightVsAvg_dphi->Divide(hPhiBGPeakRegionR_dphi);



    Float_t leftscale = hPhi2DLside->Integral(hPhi2DLside->GetXaxis()->FindBin(-1.2), hPhi2DLside->GetXaxis()->FindBin(1.2), 1, hPhi2DLside->GetYaxis()->GetNbins())/hKK2DLside->Integral(hPhi2DLside->GetXaxis()->FindBin(-1.2), hPhi2DLside->GetXaxis()->FindBin(1.2), 1, hPhi2DLside->GetYaxis()->GetNbins());
    
    TH2D* LLSsubhPhi2DLside = (TH2D*)hPhi2DLside->Clone("LLSsubhPhi2DLside");
    TH2D* LLSsubhPhi2Dpeak = (TH2D*)hPhi2Dpeak->Clone("LLSsubhPhi2Dpeak");
    LLSsubhPhi2DLside->Add(hKK2DLside, -1.0*leftscale);
    //LLSsubhPhi2DLside->Divide(hKK2DLside);
    //LLSsubhPhi2DLside->Scale(1.0/leftscale);
    LLSsubhPhi2Dpeak->Add(hKK2Dpeak, -1.0*leftscale);

    TH1D* LLSsubhPhi2DLside_deta = LLSsubhPhi2DLside->ProjectionX("LLSsubhPhi2DLside_deta", 1, LLSsubhPhi2DLside->GetYaxis()->GetNbins());
    TH1D* LLSsubhPhi2DLside_dphi = LLSsubhPhi2DLside->ProjectionY("LLSsubhPhi2DLside_dphi", LLSsubhPhi2DLside->GetXaxis()->FindBin(-1.2), LLSsubhPhi2DLside->GetXaxis()->FindBin(1.2));

    //Float_t rightscale = hPhi2DRside->Integral(1, hPhi2DRside->GetXaxis()->GetNbins(), 1, hPhi2DRside->GetYaxis()->GetNbins())/hKK2DRside->Integral(1, hPhi2DRside->GetXaxis()->GetNbins(), 1, hPhi2DRside->GetYaxis()->GetNbins());
    Float_t rightscale = hPhi2DRside->Integral(hPhi2DRside->GetXaxis()->FindBin(-1.2), hPhi2DRside->GetXaxis()->FindBin(1.2), 1, hPhi2DRside->GetYaxis()->GetNbins())/hKK2DRside->Integral(hPhi2DRside->GetXaxis()->FindBin(-1.2), hPhi2DRside->GetXaxis()->FindBin(1.2), 1, hPhi2DRside->GetYaxis()->GetNbins());
    TH2D* RLSsubhPhi2DRside = (TH2D*)hPhi2DRside->Clone("RLSsubhPhi2DRside");
    TH2D* RLSsubhPhi2Dpeak = (TH2D*)hPhi2Dpeak->Clone("RLSsubhPhi2Dpeak");
    RLSsubhPhi2DRside->Add(hKK2DRside, -1.0*rightscale);
    //RLSsubhPhi2DRside->Divide(hKK2DRside);
    //RLSsubhPhi2DRside->Scale(1.0/rightscale);
    RLSsubhPhi2Dpeak->Add(hKK2Dpeak, -1.0*rightscale);
    RLSsubhPhi2Dpeak->GetXaxis()->SetRangeUser(-1.2, 1.2);

    TH1D* RLSsubhPhi2DRside_deta = RLSsubhPhi2DRside->ProjectionX("RLSsubhPhi2DRside_deta", 1, RLSsubhPhi2DRside->GetYaxis()->GetNbins());
    TH1D* RLSsubhPhi2DRside_dphi = RLSsubhPhi2DRside->ProjectionY("RLSsubhPhi2DRside_dphi", RLSsubhPhi2DRside->GetXaxis()->FindBin(-1.2), RLSsubhPhi2DRside->GetXaxis()->FindBin(1.2));

    TH1D* RLSsubhPhi2Dpeak_deta = RLSsubhPhi2Dpeak->ProjectionX("RLSsubhPhi2Dpeak_deta", 1, RLSsubhPhi2Dpeak->GetYaxis()->GetNbins());
    TH1D* RLSsubhPhi2Dpeak_dphi = RLSsubhPhi2Dpeak->ProjectionY("RLSsubhPhi2Dpeak_dphi", RLSsubhPhi2Dpeak->GetXaxis()->FindBin(-1.2), RLSsubhPhi2Dpeak->GetXaxis()->FindBin(1.2));


    TH1D* scales = new TH1D("scales", "scales", 2, -1, 1);
    scales->SetBinContent(1, leftscale);
    scales->SetBinContent(2, rightscale);

    TH2D* rebinRLSsubhPhi2Dpeak = (TH2D*)RLSsubhPhi2Dpeak->Clone("rebinRLSsubhPhi2Dpeak");
    rebinRLSsubhPhi2Dpeak->Rebin2D(2, 2);

    //Using US estimate for BG to subtract off the from the peak region:

    Float_t scaleUS = (rightscale)*hKK2Dpeak->Integral(hKK2Dpeak->GetXaxis()->FindBin(-1.2), hKK2Dpeak->GetXaxis()->FindBin(1.2), 1, hKK2Dpeak->GetYaxis()->GetNbins());
    Float_t scaletest = (rightscale)*hPhi2Dpeak->Integral(hPhi2Dpeak->GetXaxis()->FindBin(-1.2), hPhi2Dpeak->GetXaxis()->FindBin(1.2), 1, hPhi2Dpeak->GetYaxis()->GetNbins());


    printf("\n\nscaleUS = %e\n\ntestscale = %e \n\n", scaleUS, scaletest);

    TH2D* AvgUSsubhPhi2Dpeak = (TH2D*)hPhi2Dpeak->Clone("AvgUSsubhPhi2Dpeak");
    AvgUSsubhPhi2Dpeak->Add(hPhiBGPeakRegion, -1.0*scaleUS);
    AvgUSsubhPhi2Dpeak->GetXaxis()->SetRangeUser(-1.2, 1.2);
    TH1D* AvgUSsubhPhi2Dpeak_deta = (TH1D*)AvgUSsubhPhi2Dpeak->ProjectionX("AvgUSsubhPhi2Dpeak_deta", 1, AvgUSsubhPhi2Dpeak->GetYaxis()->GetNbins());
    TH1D* AvgUSsubhPhi2Dpeak_dphi = (TH1D*)AvgUSsubhPhi2Dpeak->ProjectionY("AvgUSsubhPhi2Dpeak_dphi", AvgUSsubhPhi2Dpeak->GetXaxis()->FindBin(-1.2), AvgUSsubhPhi2Dpeak->GetXaxis()->FindBin(1.2));

    TH2D* resUSvsLS = (TH2D*)AvgUSsubhPhi2Dpeak->Clone("resUSvsLS");
    resUSvsLS->Add(RLSsubhPhi2Dpeak, -1.0);
    resUSvsLS->Divide(AvgUSsubhPhi2Dpeak);
    resUSvsLS->GetXaxis()->SetRangeUser(-1.2, 1.2);
    TH1D* resUSvsLS_deta = (TH1D*)AvgUSsubhPhi2Dpeak_deta->Clone("resUSvsLS_deta");
    resUSvsLS_deta->Add(RLSsubhPhi2Dpeak_deta, -1.0);
    resUSvsLS_deta->Divide(AvgUSsubhPhi2Dpeak_deta);
    resUSvsLS_deta->GetXaxis()->SetRangeUser(-1.2, 1.2);
    TH1D* resUSvsLS_dphi = (TH1D*)AvgUSsubhPhi2Dpeak_dphi->Clone("resUSvsLS_dphi");
    resUSvsLS_dphi->Add(RLSsubhPhi2Dpeak_dphi, -1.0);
    resUSvsLS_dphi->Divide(AvgUSsubhPhi2Dpeak_dphi);



    TFile* output = new TFile(Form("US_%s", inputFile.c_str()), "RECREATE");
    LLSsubhPhi2DLside->Write();
    LLSsubhPhi2DLside_deta->Write();
    LLSsubhPhi2DLside_dphi->Write();
    LLSsubhPhi2Dpeak->Write();
    RLSsubhPhi2DRside->Write();
    RLSsubhPhi2DRside_deta->Write();
    RLSsubhPhi2DRside_dphi->Write();
    RLSsubhPhi2Dpeak->Write();
    RLSsubhPhi2Dpeak_deta->Write();
    RLSsubhPhi2Dpeak_dphi->Write();
    rebinRLSsubhPhi2Dpeak->Write();
    scales->Write();
    hPhi2Dpeak->Write();
    hPhiBGPeakRegionL->Write();
    hPhiBGPeakRegionL_deta->Write();
    hPhiBGPeakRegionL_dphi->Write();
    hPhiBGPeakRegionR->Write();
    hPhiBGPeakRegionR_deta->Write();
    hPhiBGPeakRegionR_dphi->Write();
    hPhiBGPeakRegion->Write();
    hPhiBGPeakRegion_deta->Write();
    hPhiBGPeakRegion_dphi->Write();
    resLeftVsAvg->Write();
    resLeftVsAvg_deta->Write();
    resLeftVsAvg_dphi->Write();
    resRightVsAvg->Write();
    resRightVsAvg_deta->Write();
    resRightVsAvg_dphi->Write();
    AvgUSsubhPhi2Dpeak->Write();
    AvgUSsubhPhi2Dpeak_deta->Write();
    AvgUSsubhPhi2Dpeak_dphi->Write();
    resUSvsLS->Write();
    resUSvsLS_deta->Write();
    resUSvsLS_dphi->Write();
}
