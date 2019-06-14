void makeUSCorrections(string inputFile, float peakLow = 1.014, float peakHigh = 1.026){
    TFile* input = new TFile(inputFile.c_str());
    TH2D* hPhi2Dpeak = (TH2D*)input->Get("hPhi2Dpeak");
    TH2D* hPhi2DLside = (TH2D*)input->Get("hPhi2DLside");
    TH2D* hPhi2DRside = (TH2D*)input->Get("hPhi2DRside");
    TH2D* hKK2Dpeak = (TH2D*)input->Get("hKK2Dpeak");
    TH2D* hKK2DLside = (TH2D*)input->Get("hKK2DLside");
    TH2D* hKK2DRside = (TH2D*)input->Get("hKK2DRside");

    TH2D* trigDistSameUS = (TH2D*)input->Get("fTrigSameUSDist");
    TH2D* trigDistSameLS = (TH2D*)input->Get("fTrigSameLSDist");

    TH1D* USinvmass = (TH1D*)input->Get("USinvmassTotal");
    TH1D* LSinvmass = (TH1D*)input->Get("LSinvmassTotal");

    hPhi2Dpeak->SetName("uncorrectedhPhi2Dpeak");

    Double_t weight = 1.0;
    TString filename(inputFile.c_str());
    if(filename.Contains("PIDcut28")){
           weight = 0.995/(0.99488*0.99488);
    }else if(filename.Contains("PIDcut26")){
           weight = 0.995/(0.99098*0.99098);
    }else if(filename.Contains("PIDcut24")){
            weight = 0.995/(0.98366*0.98366);
    }else if(filename.Contains("PIDcut22")){
            weight = 0.995/(0.97220*0.97220);
    }else if(filename.Contains("PIDcut20")){
            weight = 0.995/(0.95445*0.95445);
    }
/*
    if(filename.Contains("TOFcut275")){
           weight = 0.995/(0.994*0.994);
    }else if(filename.Contains("TOFcut25")){
           weight = 0.995/(0.9756*0.9756);
    }else if(filename.Contains("TOFcut225")){
            weight = 0.995/(0.9566*0.9566);
    }
*/    
    printf("weight: %f\n", weight);

/*    Float_t totalTrigUS = trigDistSameUS->Integral(trigDistSameUS->GetXaxis()->FindBin(4.0), trigDistSameUS->GetXaxis()->FindBin(10.0));
    Float_t totalTrigLS = trigDistSameLS->Integral(trigDistSameLS->GetXaxis()->FindBin(4.0), trigDistSameLS->GetXaxis()->FindBin(10.0));
    hPhi2Dpeak->Scale(1.0/totalTrigUS);
    hPhi2DLside->Scale(1.0/totalTrigUS);
    hPhi2DRside->Scale(1.0/totalTrigUS);
    hKK2Dpeak->Scale(1.0/totalTrigLS);
    hKK2DLside->Scale(1.0/totalTrigLS);
    hKK2DRside->Scale(1.0/totalTrigLS);
*/
    //Using US sideband regions to estimate the BG under the peak region
    //Now using the (unweighted) average of the Left and Right side sideband (correct? needs additional checks?)
    
    Float_t epsilon = 0.001;

    TH2D* hPhiBGPeakRegionL = (TH2D*)hPhi2DLside->Clone("hPhiBGPeakRegionL");
    hPhiBGPeakRegionL->Scale(1.0/(hPhi2DLside->Integral(hPhi2DLside->GetXaxis()->FindBin(-1.2 + epsilon), hPhi2DLside->GetXaxis()->FindBin(1.2 - epsilon), 1, hPhi2DLside->GetYaxis()->GetNbins())));
    hPhiBGPeakRegionL->GetXaxis()->SetRangeUser(-1.2 + epsilon, 1.2 - epsilon);
    TH1D* hPhiBGPeakRegionL_deta = (TH1D*)hPhiBGPeakRegionL->ProjectionX("hPhiBGPeakRegionL_deta", 1, hPhiBGPeakRegionL->GetYaxis()->GetNbins());
    TH1D* hPhiBGPeakRegionL_dphi = (TH1D*)hPhiBGPeakRegionL->ProjectionY("hPhiBGPeakRegionL_dphi", hPhiBGPeakRegionL->GetXaxis()->FindBin(-1.2 + epsilon), hPhiBGPeakRegionL->GetXaxis()->FindBin(1.2 - epsilon));

    TH2D* hPhiBGPeakRegionR = (TH2D*)hPhi2DRside->Clone("hPhiBGPeakRegionR");
    hPhiBGPeakRegionR->Scale(1.0/(hPhi2DRside->Integral(hPhi2DRside->GetXaxis()->FindBin(-1.2 + epsilon), hPhi2DRside->GetXaxis()->FindBin(1.2 - epsilon), 1, hPhi2DRside->GetYaxis()->GetNbins())));
    hPhiBGPeakRegionR->GetXaxis()->SetRangeUser(-1.2 + epsilon, 1.2 - epsilon);
    TH1D* hPhiBGPeakRegionR_deta = (TH1D*)hPhiBGPeakRegionR->ProjectionX("hPhiBGPeakRegionR_deta", 1, hPhiBGPeakRegionR->GetYaxis()->GetNbins());
    TH1D* hPhiBGPeakRegionR_dphi = (TH1D*)hPhiBGPeakRegionR->ProjectionY("hPhiBGPeakRegionR_dphi", hPhiBGPeakRegionR->GetXaxis()->FindBin(-1.2 + epsilon), hPhiBGPeakRegionR->GetXaxis()->FindBin(1.2 - epsilon));

    TH2D* hPhiBGPeakRegion = (TH2D*)hPhiBGPeakRegionL->Clone("hPhiBGPeakregion");
    hPhiBGPeakRegion->Add(hPhiBGPeakRegionR);
    hPhiBGPeakRegion->Scale(0.5);
    hPhiBGPeakRegion->GetXaxis()->SetRangeUser(-1.2 + epsilon, 1.2 - epsilon);
    TH1D* hPhiBGPeakRegion_deta = (TH1D*)hPhiBGPeakRegion->ProjectionX("hPhiBGPeakRegion_deta", 1, hPhiBGPeakRegion->GetYaxis()->GetNbins());
    TH1D* hPhiBGPeakRegion_dphi = (TH1D*)hPhiBGPeakRegion->ProjectionY("hPhiBGPeakRegion_dphi", hPhiBGPeakRegion->GetXaxis()->FindBin(-1.2 + epsilon), hPhiBGPeakRegion->GetXaxis()->FindBin(1.2 - epsilon));


    //US residual checks between SB average and the Left and Right separately
    TH2D* resLeftVsAvg = (TH2D*)hPhiBGPeakRegionL->Clone("resLeftVsAbg");
    resLeftVsAvg->Add(hPhiBGPeakRegion, -1.0);
    resLeftVsAvg->Divide(hPhiBGPeakRegionL);
    resLeftVsAvg->GetXaxis()->SetRangeUser(-1.2 + epsilon, 1.2 - epsilon);
    TH1D* resLeftVsAvg_deta = (TH1D*)hPhiBGPeakRegionL_deta->Clone("resLeftVsAvg_deta");
    resLeftVsAvg_deta->Add(hPhiBGPeakRegion_deta, -1.0);
    resLeftVsAvg_deta->Divide(hPhiBGPeakRegionL_deta);
    resLeftVsAvg_deta->GetXaxis()->SetRangeUser(-1.2 + epsilon, 1.2 - epsilon);
    TH1D* resLeftVsAvg_dphi = (TH1D*)hPhiBGPeakRegionL_dphi->Clone("resLeftVsAvg_dphi");
    resLeftVsAvg_dphi->Add(hPhiBGPeakRegion_dphi, -1.0);
    resLeftVsAvg_dphi->Divide(hPhiBGPeakRegionL_dphi);

    TH2D* resRightVsAvg = (TH2D*)hPhiBGPeakRegionR->Clone("resRightVsAbg");
    resRightVsAvg->Add(hPhiBGPeakRegion, -1.0);
    resRightVsAvg->Divide(hPhiBGPeakRegionR);
    resRightVsAvg->GetXaxis()->SetRangeUser(-1.2 + epsilon, 1.2 - epsilon);
    TH1D* resRightVsAvg_deta = (TH1D*)hPhiBGPeakRegionR_deta->Clone("resRightVsAvg_deta");
    resRightVsAvg_deta->Add(hPhiBGPeakRegion_deta, -1.0);
    resRightVsAvg_deta->Divide(hPhiBGPeakRegionR_deta);
    resRightVsAvg_deta->GetXaxis()->SetRangeUser(-1.2 + epsilon, 1.2 - epsilon);
    TH1D* resRightVsAvg_dphi = (TH1D*)hPhiBGPeakRegionR_dphi->Clone("resRightVsAvg_dphi");
    resRightVsAvg_dphi->Add(hPhiBGPeakRegion_dphi, -1.0);
    resRightVsAvg_dphi->Divide(hPhiBGPeakRegionR_dphi);



    Double_t intbuff1 = 0.0;
    Double_t intbuff2 = 0.0;
    Double_t errbuff1 = 0.0;
    Double_t errbuff2 = 0.0;
    Double_t leftscale = hPhi2DLside->Integral(hPhi2DLside->GetXaxis()->FindBin(-1.2 + epsilon), hPhi2DLside->GetXaxis()->FindBin(1.2 - epsilon), 1, hPhi2DLside->GetYaxis()->GetNbins())/hKK2DLside->Integral(hPhi2DLside->GetXaxis()->FindBin(-1.2 + epsilon), hPhi2DLside->GetXaxis()->FindBin(1.2 - epsilon), 1, hPhi2DLside->GetYaxis()->GetNbins());
    intbuff1 = hPhi2DLside->IntegralAndError(hPhi2DLside->GetXaxis()->FindBin(-1.2 + epsilon), hPhi2DLside->GetXaxis()->FindBin(1.2 - epsilon), 1, hPhi2DLside->GetYaxis()->GetNbins(), errbuff1);
    intbuff2 = hKK2DLside->IntegralAndError(hKK2DLside->GetXaxis()->FindBin(-1.2 + epsilon), hKK2DLside->GetXaxis()->FindBin(1.2 - epsilon), 1, hKK2DLside->GetYaxis()->GetNbins(), errbuff2);
    Double_t lefterror = (intbuff1/intbuff2)*TMath::Sqrt(TMath::Power(errbuff1/intbuff1, 2) + TMath::Power(errbuff2/intbuff2, 2));


    TH2D* LLSsubhPhi2DLside = (TH2D*)hPhi2DLside->Clone("LLSsubhPhi2DLside");
    TH2D* LLSsubhPhi2Dpeak = (TH2D*)hPhi2Dpeak->Clone("LLSsubhPhi2Dpeak");
    LLSsubhPhi2DLside->Add(hKK2DLside, -1.0*leftscale);
    //LLSsubhPhi2DLside->Divide(hKK2DLside);
    //LLSsubhPhi2DLside->Scale(1.0/leftscale);
    LLSsubhPhi2Dpeak->Add(hKK2Dpeak, -1.0*leftscale);

    TH1D* LLSsubhPhi2DLside_deta = LLSsubhPhi2DLside->ProjectionX("LLSsubhPhi2DLside_deta", 1, LLSsubhPhi2DLside->GetYaxis()->GetNbins());
    TH1D* LLSsubhPhi2DLside_dphi = LLSsubhPhi2DLside->ProjectionY("LLSsubhPhi2DLside_dphi", LLSsubhPhi2DLside->GetXaxis()->FindBin(-1.2 + epsilon), LLSsubhPhi2DLside->GetXaxis()->FindBin(1.2 - epsilon));

    Double_t rightscale = hPhi2DRside->Integral(hPhi2DRside->GetXaxis()->FindBin(-1.2 + epsilon), hPhi2DRside->GetXaxis()->FindBin(1.2 - epsilon), 1, hPhi2DRside->GetYaxis()->GetNbins())/hKK2DRside->Integral(hPhi2DRside->GetXaxis()->FindBin(-1.2 + epsilon), hPhi2DRside->GetXaxis()->FindBin(1.2 - epsilon), 1, hPhi2DRside->GetYaxis()->GetNbins());
    intbuff1 = hPhi2DRside->IntegralAndError(hPhi2DRside->GetXaxis()->FindBin(-1.2 + epsilon), hPhi2DRside->GetXaxis()->FindBin(1.2 - epsilon), 1, hPhi2DRside->GetYaxis()->GetNbins(), errbuff1);
    intbuff2 = hKK2DRside->IntegralAndError(hKK2DRside->GetXaxis()->FindBin(-1.2 + epsilon), hKK2DRside->GetXaxis()->FindBin(1.2 - epsilon), 1, hKK2DRside->GetYaxis()->GetNbins(), errbuff2);
    Double_t righterror = (intbuff1/intbuff2)*TMath::Sqrt(TMath::Power(errbuff1/intbuff1, 2) + TMath::Power(errbuff2/intbuff2, 2));


    TH2D* RLSsubhPhi2DRside = (TH2D*)hPhi2DRside->Clone("RLSsubhPhi2DRside");
    TH2D* RLSsubhPhi2Dpeak = (TH2D*)hPhi2Dpeak->Clone("RLSsubhPhi2Dpeak");
    RLSsubhPhi2DRside->Add(hKK2DRside, -1.0*rightscale);
    //RLSsubhPhi2DRside->Divide(hKK2DRside);
    //RLSsubhPhi2DRside->Scale(1.0/rightscale);
    RLSsubhPhi2Dpeak->Add(hKK2Dpeak, -1.0*rightscale);
    RLSsubhPhi2Dpeak->GetXaxis()->SetRangeUser(-1.2 + epsilon, 1.2 - epsilon);

    TH1D* RLSsubhPhi2DRside_deta = RLSsubhPhi2DRside->ProjectionX("RLSsubhPhi2DRside_deta", 1, RLSsubhPhi2DRside->GetYaxis()->GetNbins());
    TH1D* RLSsubhPhi2DRside_dphi = RLSsubhPhi2DRside->ProjectionY("RLSsubhPhi2DRside_dphi", RLSsubhPhi2DRside->GetXaxis()->FindBin(-1.2 + epsilon), RLSsubhPhi2DRside->GetXaxis()->FindBin(1.2 - epsilon));

    TH1D* RLSsubhPhi2Dpeak_deta = RLSsubhPhi2Dpeak->ProjectionX("RLSsubhPhi2Dpeak_deta", 1, RLSsubhPhi2Dpeak->GetYaxis()->GetNbins());
    TH1D* RLSsubhPhi2Dpeak_dphi = RLSsubhPhi2Dpeak->ProjectionY("RLSsubhPhi2Dpeak_dphi", RLSsubhPhi2Dpeak->GetXaxis()->FindBin(-1.2 + epsilon), RLSsubhPhi2Dpeak->GetXaxis()->FindBin(1.2 - epsilon));


    TH2D* rebinRLSsubhPhi2Dpeak = (TH2D*)RLSsubhPhi2Dpeak->Clone("rebinRLSsubhPhi2Dpeak");
    rebinRLSsubhPhi2Dpeak->Rebin2D(2, 2);

    //Using US estimate for BG to subtract off the from the peak region:

    Double_t avgscale = 0.5*(leftscale+rightscale);
    Double_t avgerror = TMath::Sqrt(TMath::Power(lefterror, 2) + TMath::Power(righterror, 2));

    Double_t peakerror = 0.0;
    Double_t peakscale = hKK2Dpeak->IntegralAndError(hKK2Dpeak->GetXaxis()->FindBin(-1.2 + epsilon), hKK2Dpeak->GetXaxis()->FindBin(1.2 - epsilon), 1, hKK2Dpeak->GetYaxis()->GetNbins(), peakerror);
    Double_t scaleUS = (rightscale)*peakscale;
    Double_t fullscaleright = rightscale*peakscale;
    Double_t fullscaleleft = leftscale*peakscale;
    Double_t fullscaleavg = avgscale*peakscale;
    Double_t fullerrorright = scaleUS*TMath::Sqrt(TMath::Power(righterror/rightscale, 2) + TMath::Power(peakerror/peakscale, 2));
    Double_t fullerrorleft = (leftscale)*peakscale*TMath::Sqrt(TMath::Power(lefterror/leftscale, 2) + TMath::Power(peakerror/peakscale, 2));
    Double_t fullerroravg = 0.5*(avgscale)*peakscale*TMath::Sqrt(TMath::Power(avgerror/avgscale, 2) + TMath::Power(peakerror/peakscale, 2));

    TH1D* scales = new TH1D("scales", "scales", 6, 0, 6);
    scales->SetBinContent(1, leftscale);
    scales->SetBinError(1, lefterror);
    scales->SetBinContent(2, rightscale);
    scales->SetBinError(2, righterror);
    scales->SetBinContent(3, avgscale);
    scales->SetBinError(3, avgerror);
    scales->SetBinContent(4, fullscaleright);
    scales->SetBinError(4, fullerrorright);
    scales->SetBinContent(5, fullscaleleft);
    scales->SetBinError(5, fullerrorleft);
    scales->SetBinContent(6, fullscaleavg);
    scales->SetBinError(6, fullerroravg);


    Double_t scaletest = (rightscale)*hPhi2Dpeak->Integral(hPhi2Dpeak->GetXaxis()->FindBin(-1.2 + epsilon), hPhi2Dpeak->GetXaxis()->FindBin(1.2 - epsilon), 1, hPhi2Dpeak->GetYaxis()->GetNbins());

    printf("\n\nscaleUS = %e\n\ntestscale = %e \n\n", scaleUS, scaletest);

    //do inv. mass reconstruction efficiency
    TH1D* corrMass = (TH1D*)USinvmass->Clone("correctedMass");
    corrMass->Add(LSinvmass, -1.0*rightscale);

    TF1 *massfitR = new TF1("massfitR", "[0]*TMath::Voigt(x-[1], [2], [3], 4) + pol2(4)", 0.99, 1.07);
    massfitR->SetParameter(0, 1000);
    massfitR->SetParameter(1, 1.020);
    massfitR->SetParameter(2, 0.001);
    massfitR->FixParameter(3, 0.00426);
    massfitR->SetParLimits(1, 1.015, 1.025);
    massfitR->SetParLimits(2, 0.0001, 0.005);

    corrMass->Fit(massfitR);

    TF1 *massBGR = new TF1("massBGR", "pol2(0)", 0.99, 1.07);
    massBGR->SetParameters(massfitR->GetParameter(4), massfitR->GetParameter(5), massfitR->GetParameter(6));
    Double_t fullfitIntR = massfitR->Integral(1.0, 1.05) - massBGR->Integral(1.0, 1.05);
    Double_t peakfitIntR = massfitR->Integral(peakLow, peakHigh) - massBGR->Integral(peakLow, peakHigh);
    Double_t fitratioR = peakfitIntR/fullfitIntR;

    Double_t peakregion = corrMass->Integral(corrMass->GetXaxis()->FindBin(peakLow + 0.0001), corrMass->GetXaxis()->FindBin(peakHigh - 0.0001));
    Double_t totalregion = corrMass->Integral(corrMass->GetXaxis()->FindBin(1.0001), corrMass->GetYaxis()->FindBin(1.04999));
    Double_t fraction = peakregion/totalregion;

    printf("peak: %f, total: %f, fraction: %f\n", peakregion, totalregion, fraction); 
    
    TH1D* corrMassL = (TH1D*)USinvmass->Clone("correctedMassL");
    corrMassL->Add(LSinvmass, -1.0*leftscale);

    TF1 *massfitL = new TF1("massfitL", "[0]*TMath::Voigt(x-[1], [2], [3], 5) + pol2(4)", 0.99, 1.07);
    massfitL->SetParameter(0, 1000);
    massfitL->SetParameter(1, 1.020);
    massfitL->SetParameter(2, 0.0002);
    massfitL->FixParameter(3, 0.00426);
    massfitL->SetParLimits(1, 1.015, 1.025);

    corrMassL->Fit(massfitL);

    TF1 *massBGL = new TF1("massBGL", "pol2(0)", 0.99, 1.07);
    massBGL->SetParameters(massfitL->GetParameter(4), massfitL->GetParameter(5), massfitL->GetParameter(6));
    Double_t fullfitIntL = massfitL->Integral(1.0, 1.05) - massBGL->Integral(1.0, 1.05);
    Double_t peakfitIntL = massfitL->Integral(peakLow, peakHigh) - massBGL->Integral(peakLow, peakHigh);
    Double_t fitratioL = peakfitIntL/fullfitIntL;


    TH1D* corrMassAvg = (TH1D*)USinvmass->Clone("correctedMassAvg");
    corrMassAvg->Add(LSinvmass, -1.0*(avgscale+avgerror));

    TF1 *massfitAvg = new TF1("massfitAvg", "[0]*TMath::Voigt(x-[1], [2], [3], 5) + pol2(4)", 0.99, 1.07);
    massfitAvg->SetParameter(0, 1000);
    massfitAvg->SetParameter(1, 1.020);
    massfitAvg->SetParameter(2, 0.0002);
    massfitAvg->FixParameter(3, 0.00426);
    massfitAvg->SetParLimits(1, 1.015, 1.025);

    corrMassAvg->Fit(massfitAvg);

    TF1 *massBGAvg = new TF1("massBGAvg", "pol2(0)", 0.99, 1.07);
    massBGAvg->SetParameters(massfitAvg->GetParameter(4), massfitAvg->GetParameter(5), massfitAvg->GetParameter(6));
    Double_t fullfitIntAvg = massfitAvg->Integral(1.0, 1.05) - massBGAvg->Integral(1.0, 1.05);
    Double_t peakfitIntAvg = massfitAvg->Integral(peakLow, peakHigh) - massBGAvg->Integral(peakLow, peakHigh);
    Double_t fitratioAvg = peakfitIntAvg/fullfitIntAvg;


    printf("fitpeakR: %f, fittotalR: %f, fitfractionR: %f\n", peakfitIntR, fullfitIntR, fitratioR);
    printf("fitpeakL: %f, fittotalL: %f, fitfractionL: %f\n", peakfitIntL, fullfitIntL, fitratioL);
    printf("fitpeakAvg: %f, fittotalAvg: %f, fitfractionAvg: %f\n", peakfitIntAvg, fullfitIntAvg, fitratioAvg);

    fitratioR = fitratioR*(0.49); //account for branching ratio fraction for KK

    //avg of right and left US sideband tests
    TH2D* AvgUSsubhPhi2Dpeak = (TH2D*)hPhi2Dpeak->Clone("AvgUSsubhPhi2Dpeak");
    AvgUSsubhPhi2Dpeak->Add(hPhiBGPeakRegion, -1.0*scaleUS);
    AvgUSsubhPhi2Dpeak->GetXaxis()->SetRangeUser(-1.2 + epsilon, 1.2 - epsilon);
    AvgUSsubhPhi2Dpeak->Scale(weight*1.0/fitratioR);
    TH1D* AvgUSsubhPhi2Dpeak_deta = (TH1D*)AvgUSsubhPhi2Dpeak->ProjectionX("AvgUSsubhPhi2Dpeak_deta", 1, AvgUSsubhPhi2Dpeak->GetYaxis()->GetNbins());
    TH1D* AvgUSsubhPhi2Dpeak_dphi = (TH1D*)AvgUSsubhPhi2Dpeak->ProjectionY("AvgUSsubhPhi2Dpeak_dphi", AvgUSsubhPhi2Dpeak->GetXaxis()->FindBin(-1.2 + epsilon), AvgUSsubhPhi2Dpeak->GetXaxis()->FindBin(1.2 - epsilon));

    //using right scale +- errors
    TH2D* AvgUSsubhPhi2Dpeakrightscale = (TH2D*)hPhi2Dpeak->Clone("AvgUSsubhPhi2Dpeakrightscale");
    AvgUSsubhPhi2Dpeakrightscale->Add(hPhiBGPeakRegion, -1.0*fullscaleright);
    AvgUSsubhPhi2Dpeakrightscale->GetXaxis()->SetRangeUser(-1.2 + epsilon, 1.2 - epsilon);
    AvgUSsubhPhi2Dpeakrightscale->Scale(weight*1.0/fitratioR);
    TH1D* AvgUSsubhPhi2Dpeakrightscale_deta = (TH1D*)AvgUSsubhPhi2Dpeakrightscale->ProjectionX("AvgUSsubhPhi2Dpeakrightscale_deta", 1, AvgUSsubhPhi2Dpeakrightscale->GetYaxis()->GetNbins());
    TH1D* AvgUSsubhPhi2Dpeakrightscale_dphi = (TH1D*)AvgUSsubhPhi2Dpeakrightscale->ProjectionY("AvgUSsubhPhi2Dpeakrightscale_dphi", AvgUSsubhPhi2Dpeakrightscale->GetXaxis()->FindBin(-1.2 + epsilon), AvgUSsubhPhi2Dpeakrightscale->GetXaxis()->FindBin(1.2 - epsilon));

    TH2D* AvgUSsubhPhi2Dpeakrightscaleplus = (TH2D*)hPhi2Dpeak->Clone("AvgUSsubhPhi2Dpeakrightscaleplus");
    AvgUSsubhPhi2Dpeakrightscaleplus->Add(hPhiBGPeakRegion, -1.0*(fullscaleright+fullerrorright));
    AvgUSsubhPhi2Dpeakrightscaleplus->GetXaxis()->SetRangeUser(-1.2 + epsilon, 1.2 - epsilon);
    AvgUSsubhPhi2Dpeakrightscaleplus->Scale(weight*1.0/fitratioR);
    TH1D* AvgUSsubhPhi2Dpeakrightscaleplus_deta = (TH1D*)AvgUSsubhPhi2Dpeakrightscaleplus->ProjectionX("AvgUSsubhPhi2Dpeakrightscaleplus_deta", 1, AvgUSsubhPhi2Dpeakrightscaleplus->GetYaxis()->GetNbins());
    TH1D* AvgUSsubhPhi2Dpeakrightscaleplus_dphi = (TH1D*)AvgUSsubhPhi2Dpeakrightscaleplus->ProjectionY("AvgUSsubhPhi2Dpeakrightscaleplus_dphi", AvgUSsubhPhi2Dpeakrightscaleplus->GetXaxis()->FindBin(-1.2 + epsilon), AvgUSsubhPhi2Dpeakrightscaleplus->GetXaxis()->FindBin(1.2 - epsilon));

    TH2D* AvgUSsubhPhi2Dpeakrightscaleminus = (TH2D*)hPhi2Dpeak->Clone("AvgUSsubhPhi2Dpeakrightscaleminus");
    AvgUSsubhPhi2Dpeakrightscaleminus->Add(hPhiBGPeakRegion, -1.0*(fullscaleright-fullerrorright));
    AvgUSsubhPhi2Dpeakrightscaleminus->GetXaxis()->SetRangeUser(-1.2 + epsilon, 1.2 - epsilon);
    AvgUSsubhPhi2Dpeakrightscaleminus->Scale(weight*1.0/fitratioR);
    TH1D* AvgUSsubhPhi2Dpeakrightscaleminus_deta = (TH1D*)AvgUSsubhPhi2Dpeakrightscaleminus->ProjectionX("AvgUSsubhPhi2Dpeakrightscaleminus_deta", 1, AvgUSsubhPhi2Dpeakrightscaleminus->GetYaxis()->GetNbins());
    TH1D* AvgUSsubhPhi2Dpeakrightscaleminus_dphi = (TH1D*)AvgUSsubhPhi2Dpeakrightscaleminus->ProjectionY("AvgUSsubhPhi2Dpeakrightscaleminus_dphi", AvgUSsubhPhi2Dpeakrightscaleminus->GetXaxis()->FindBin(-1.2 + epsilon), AvgUSsubhPhi2Dpeakrightscaleminus->GetXaxis()->FindBin(1.2 - epsilon));

    //Using left scale +- errors
    TH2D* AvgUSsubhPhi2Dpeakleftscale = (TH2D*)hPhi2Dpeak->Clone("AvgUSsubhPhi2Dpeakleftscale");
    AvgUSsubhPhi2Dpeakleftscale->Add(hPhiBGPeakRegion, -1.0*fullscaleleft);
    AvgUSsubhPhi2Dpeakleftscale->GetXaxis()->SetRangeUser(-1.2 + epsilon, 1.2 - epsilon);
    AvgUSsubhPhi2Dpeakleftscale->Scale(weight*1.0/fitratioR);
    TH1D* AvgUSsubhPhi2Dpeakleftscale_deta = (TH1D*)AvgUSsubhPhi2Dpeakleftscale->ProjectionX("AvgUSsubhPhi2Dpeakleftscale_deta", 1, AvgUSsubhPhi2Dpeakleftscale->GetYaxis()->GetNbins());
    TH1D* AvgUSsubhPhi2Dpeakleftscale_dphi = (TH1D*)AvgUSsubhPhi2Dpeakleftscale->ProjectionY("AvgUSsubhPhi2Dpeakleftscale_dphi", AvgUSsubhPhi2Dpeakleftscale->GetXaxis()->FindBin(-1.2 + epsilon), AvgUSsubhPhi2Dpeakleftscale->GetXaxis()->FindBin(1.2 - epsilon));

    TH2D* AvgUSsubhPhi2Dpeakleftscaleplus = (TH2D*)hPhi2Dpeak->Clone("AvgUSsubhPhi2Dpeakleftscaleplus");
    AvgUSsubhPhi2Dpeakleftscaleplus->Add(hPhiBGPeakRegion, -1.0*(fullscaleleft+fullerrorleft));
    AvgUSsubhPhi2Dpeakleftscaleplus->GetXaxis()->SetRangeUser(-1.2 + epsilon, 1.2 - epsilon);
    AvgUSsubhPhi2Dpeakleftscaleplus->Scale(weight*1.0/fitratioR);
    TH1D* AvgUSsubhPhi2Dpeakleftscaleplus_deta = (TH1D*)AvgUSsubhPhi2Dpeakleftscaleplus->ProjectionX("AvgUSsubhPhi2Dpeakleftscaleplus_deta", 1, AvgUSsubhPhi2Dpeakleftscaleplus->GetYaxis()->GetNbins());
    TH1D* AvgUSsubhPhi2Dpeakleftscaleplus_dphi = (TH1D*)AvgUSsubhPhi2Dpeakleftscaleplus->ProjectionY("AvgUSsubhPhi2Dpeakleftscaleplus_dphi", AvgUSsubhPhi2Dpeakleftscaleplus->GetXaxis()->FindBin(-1.2 + epsilon), AvgUSsubhPhi2Dpeakleftscaleplus->GetXaxis()->FindBin(1.2 - epsilon));

    TH2D* AvgUSsubhPhi2Dpeakleftscaleminus = (TH2D*)hPhi2Dpeak->Clone("AvgUSsubhPhi2Dpeakleftscaleminus");
    AvgUSsubhPhi2Dpeakleftscaleminus->Add(hPhiBGPeakRegion, -1.0*(fullscaleleft-fullerrorleft));
    AvgUSsubhPhi2Dpeakleftscaleminus->GetXaxis()->SetRangeUser(-1.2 + epsilon, 1.2 - epsilon);
    AvgUSsubhPhi2Dpeakleftscaleminus->Scale(weight*1.0/fitratioR);
    TH1D* AvgUSsubhPhi2Dpeakleftscaleminus_deta = (TH1D*)AvgUSsubhPhi2Dpeakleftscaleminus->ProjectionX("AvgUSsubhPhi2Dpeakleftscaleminus_deta", 1, AvgUSsubhPhi2Dpeakleftscaleminus->GetYaxis()->GetNbins());
    TH1D* AvgUSsubhPhi2Dpeakleftscaleminus_dphi = (TH1D*)AvgUSsubhPhi2Dpeakleftscaleminus->ProjectionY("AvgUSsubhPhi2Dpeakleftscaleminus_dphi", AvgUSsubhPhi2Dpeakleftscaleminus->GetXaxis()->FindBin(-1.2 + epsilon), AvgUSsubhPhi2Dpeakleftscaleminus->GetXaxis()->FindBin(1.2 - epsilon));

    //Using avg scale +- errors
    TH2D* AvgUSsubhPhi2Dpeakavgscale = (TH2D*)hPhi2Dpeak->Clone("AvgUSsubhPhi2Dpeakavgscale");
    AvgUSsubhPhi2Dpeakavgscale->Add(hPhiBGPeakRegion, -1.0*fullscaleavg);
    AvgUSsubhPhi2Dpeakavgscale->GetXaxis()->SetRangeUser(-1.2 + epsilon, 1.2 - epsilon);
    AvgUSsubhPhi2Dpeakavgscale->Scale(weight*1.0/fitratioR);
    TH1D* AvgUSsubhPhi2Dpeakavgscale_deta = (TH1D*)AvgUSsubhPhi2Dpeakavgscale->ProjectionX("AvgUSsubhPhi2Dpeakavgscale_deta", 1, AvgUSsubhPhi2Dpeakavgscale->GetYaxis()->GetNbins());
    TH1D* AvgUSsubhPhi2Dpeakavgscale_dphi = (TH1D*)AvgUSsubhPhi2Dpeakavgscale->ProjectionY("AvgUSsubhPhi2Dpeakavgscale_dphi", AvgUSsubhPhi2Dpeakavgscale->GetXaxis()->FindBin(-1.2 + epsilon), AvgUSsubhPhi2Dpeakavgscale->GetXaxis()->FindBin(1.2 - epsilon));

    TH2D* AvgUSsubhPhi2Dpeakavgscaleplus = (TH2D*)hPhi2Dpeak->Clone("AvgUSsubhPhi2Dpeakavgscaleplus");
    AvgUSsubhPhi2Dpeakavgscaleplus->Add(hPhiBGPeakRegion, -1.0*(fullscaleavg+fullerroravg));
    AvgUSsubhPhi2Dpeakavgscaleplus->GetXaxis()->SetRangeUser(-1.2 + epsilon, 1.2 - epsilon);
    AvgUSsubhPhi2Dpeakavgscaleplus->Scale(weight*1.0/fitratioR);
    TH1D* AvgUSsubhPhi2Dpeakavgscaleplus_deta = (TH1D*)AvgUSsubhPhi2Dpeakavgscaleplus->ProjectionX("AvgUSsubhPhi2Dpeakavgscaleplus_deta", 1, AvgUSsubhPhi2Dpeakavgscaleplus->GetYaxis()->GetNbins());
    TH1D* AvgUSsubhPhi2Dpeakavgscaleplus_dphi = (TH1D*)AvgUSsubhPhi2Dpeakavgscaleplus->ProjectionY("AvgUSsubhPhi2Dpeakavgscaleplus_dphi", AvgUSsubhPhi2Dpeakavgscaleplus->GetXaxis()->FindBin(-1.2 + epsilon), AvgUSsubhPhi2Dpeakavgscaleplus->GetXaxis()->FindBin(1.2 - epsilon));

    TH2D* AvgUSsubhPhi2Dpeakavgscaleminus = (TH2D*)hPhi2Dpeak->Clone("AvgUSsubhPhi2Dpeakavgscaleminus");
    AvgUSsubhPhi2Dpeakavgscaleminus->Add(hPhiBGPeakRegion, -1.0*(fullscaleavg-fullerroravg));
    AvgUSsubhPhi2Dpeakavgscaleminus->GetXaxis()->SetRangeUser(-1.2 + epsilon, 1.2 - epsilon);
    AvgUSsubhPhi2Dpeakavgscaleminus->Scale(weight*1.0/fitratioR);
    TH1D* AvgUSsubhPhi2Dpeakavgscaleminus_deta = (TH1D*)AvgUSsubhPhi2Dpeakavgscaleminus->ProjectionX("AvgUSsubhPhi2Dpeakavgscaleminus_deta", 1, AvgUSsubhPhi2Dpeakavgscaleminus->GetYaxis()->GetNbins());
    TH1D* AvgUSsubhPhi2Dpeakavgscaleminus_dphi = (TH1D*)AvgUSsubhPhi2Dpeakavgscaleminus->ProjectionY("AvgUSsubhPhi2Dpeakavgscaleminus_dphi", AvgUSsubhPhi2Dpeakavgscaleminus->GetXaxis()->FindBin(-1.2 + epsilon), AvgUSsubhPhi2Dpeakavgscaleminus->GetXaxis()->FindBin(1.2 - epsilon));

    

    //right side US sideband tests
    TH2D* RSUSsubhPhi2Dpeak = (TH2D*)hPhi2Dpeak->Clone("RSUSsubhPhi2Dpeak");
    RSUSsubhPhi2Dpeak->Add(hPhiBGPeakRegionR, -1.0*scaleUS);
    RSUSsubhPhi2Dpeak->GetXaxis()->SetRangeUser(-1.2 + epsilon, 1.2 - epsilon);
    TH1D* RSUSsubhPhi2Dpeak_deta = (TH1D*)RSUSsubhPhi2Dpeak->ProjectionX("RSUSsubhPhi2Dpeak_deta", 1, RSUSsubhPhi2Dpeak->GetYaxis()->GetNbins());
    TH1D* RSUSsubhPhi2Dpeak_dphi = (TH1D*)RSUSsubhPhi2Dpeak->ProjectionY("RSUSsubhPhi2Dpeak_dphi", RSUSsubhPhi2Dpeak->GetXaxis()->FindBin(-1.2 + epsilon), RSUSsubhPhi2Dpeak->GetXaxis()->FindBin(1.2 - epsilon));

    TH2D* RSUSsubhPhi2Dpeakleftscale = (TH2D*)hPhi2Dpeak->Clone("RSUSsubhPhi2Dpeakleftscale");
    RSUSsubhPhi2Dpeakleftscale->Add(hPhiBGPeakRegionR, -1.0*scaleUS*leftscale/rightscale);
    RSUSsubhPhi2Dpeakleftscale->GetXaxis()->SetRangeUser(-1.2 + epsilon, 1.2 - epsilon);
    TH1D* RSUSsubhPhi2Dpeakleftscale_deta = (TH1D*)RSUSsubhPhi2Dpeakleftscale->ProjectionX("RSUSsubhPhi2Dpeakleftscale_deta", 1, RSUSsubhPhi2Dpeakleftscale->GetYaxis()->GetNbins());
    TH1D* RSUSsubhPhi2Dpeakleftscale_dphi = (TH1D*)RSUSsubhPhi2Dpeakleftscale->ProjectionY("RSUSsubhPhi2Dpeakleftscale_dphi", RSUSsubhPhi2Dpeakleftscale->GetXaxis()->FindBin(-1.2 + epsilon), RSUSsubhPhi2Dpeakleftscale->GetXaxis()->FindBin(1.2 - epsilon));
    
    TH2D* RSUSsubhPhi2Dpeakavgscale = (TH2D*)hPhi2Dpeak->Clone("RSUSsubhPhi2Dpeakavgscale");
    RSUSsubhPhi2Dpeakavgscale->Add(hPhiBGPeakRegionR, -1.0*scaleUS*(leftscale+rightscale)/(2.0*rightscale));
    RSUSsubhPhi2Dpeakavgscale->GetXaxis()->SetRangeUser(-1.2 + epsilon, 1.2 - epsilon);
    TH1D* RSUSsubhPhi2Dpeakavgscale_deta = (TH1D*)RSUSsubhPhi2Dpeakavgscale->ProjectionX("RSUSsubhPhi2Dpeakavgscale_deta", 1, RSUSsubhPhi2Dpeakavgscale->GetYaxis()->GetNbins());
    TH1D* RSUSsubhPhi2Dpeakavgscale_dphi = (TH1D*)RSUSsubhPhi2Dpeakavgscale->ProjectionY("RSUSsubhPhi2Dpeakavgscale_dphi", RSUSsubhPhi2Dpeakavgscale->GetXaxis()->FindBin(-1.2 + epsilon), RSUSsubhPhi2Dpeakavgscale->GetXaxis()->FindBin(1.2 - epsilon));

    //left side US sideband tests
    TH2D* LSUSsubhPhi2Dpeak = (TH2D*)hPhi2Dpeak->Clone("LSUSsubhPhi2Dpeak");
    LSUSsubhPhi2Dpeak->Add(hPhiBGPeakRegionL, -1.0*scaleUS);
    LSUSsubhPhi2Dpeak->GetXaxis()->SetRangeUser(-1.2 + epsilon, 1.2 - epsilon);
    TH1D* LSUSsubhPhi2Dpeak_deta = (TH1D*)LSUSsubhPhi2Dpeak->ProjectionX("LSUSsubhPhi2Dpeak_deta", 1, LSUSsubhPhi2Dpeak->GetYaxis()->GetNbins());
    TH1D* LSUSsubhPhi2Dpeak_dphi = (TH1D*)LSUSsubhPhi2Dpeak->ProjectionY("LSUSsubhPhi2Dpeak_dphi", LSUSsubhPhi2Dpeak->GetXaxis()->FindBin(-1.2 + epsilon), LSUSsubhPhi2Dpeak->GetXaxis()->FindBin(1.2 - epsilon));

    TH2D* LSUSsubhPhi2Dpeakleftscale = (TH2D*)hPhi2Dpeak->Clone("LSUSsubhPhi2Dpeakleftscale");
    LSUSsubhPhi2Dpeakleftscale->Add(hPhiBGPeakRegionL, -1.0*scaleUS*leftscale/rightscale);
    LSUSsubhPhi2Dpeakleftscale->GetXaxis()->SetRangeUser(-1.2 + epsilon, 1.2 - epsilon);
    TH1D* LSUSsubhPhi2Dpeakleftscale_deta = (TH1D*)LSUSsubhPhi2Dpeakleftscale->ProjectionX("LSUSsubhPhi2Dpeakleftscale_deta", 1, LSUSsubhPhi2Dpeakleftscale->GetYaxis()->GetNbins());
    TH1D* LSUSsubhPhi2Dpeakleftscale_dphi = (TH1D*)LSUSsubhPhi2Dpeakleftscale->ProjectionY("LSUSsubhPhi2Dpeakleftscale_dphi", LSUSsubhPhi2Dpeakleftscale->GetXaxis()->FindBin(-1.2 + epsilon), LSUSsubhPhi2Dpeakleftscale->GetXaxis()->FindBin(1.2 - epsilon));
    
    TH2D* LSUSsubhPhi2Dpeakavgscale = (TH2D*)hPhi2Dpeak->Clone("LSUSsubhPhi2Dpeakavgscale");
    LSUSsubhPhi2Dpeakavgscale->Add(hPhiBGPeakRegionL, -1.0*scaleUS*(leftscale+rightscale)/(2.0*rightscale));
    LSUSsubhPhi2Dpeakavgscale->GetXaxis()->SetRangeUser(-1.2 + epsilon, 1.2 - epsilon);
    TH1D* LSUSsubhPhi2Dpeakavgscale_deta = (TH1D*)LSUSsubhPhi2Dpeakavgscale->ProjectionX("LSUSsubhPhi2Dpeakavgscale_deta", 1, LSUSsubhPhi2Dpeakavgscale->GetYaxis()->GetNbins());
    TH1D* LSUSsubhPhi2Dpeakavgscale_dphi = (TH1D*)LSUSsubhPhi2Dpeakavgscale->ProjectionY("LSUSsubhPhi2Dpeakavgscale_dphi", LSUSsubhPhi2Dpeakavgscale->GetXaxis()->FindBin(-1.2 + epsilon), LSUSsubhPhi2Dpeakavgscale->GetXaxis()->FindBin(1.2 - epsilon));


    TH2D* resUSvsLS = (TH2D*)AvgUSsubhPhi2Dpeak->Clone("resUSvsLS");
    resUSvsLS->Add(RLSsubhPhi2Dpeak, -1.0);
    resUSvsLS->Divide(AvgUSsubhPhi2Dpeak);
    resUSvsLS->GetXaxis()->SetRangeUser(-1.2 + epsilon, 1.2 - epsilon);
    TH1D* resUSvsLS_deta = (TH1D*)AvgUSsubhPhi2Dpeak_deta->Clone("resUSvsLS_deta");
    resUSvsLS_deta->Add(RLSsubhPhi2Dpeak_deta, -1.0);
    resUSvsLS_deta->Divide(AvgUSsubhPhi2Dpeak_deta);
    resUSvsLS_deta->GetXaxis()->SetRangeUser(-1.2 + epsilon, 1.2 - epsilon);
    TH1D* resUSvsLS_dphi = (TH1D*)AvgUSsubhPhi2Dpeak_dphi->Clone("resUSvsLS_dphi");
    resUSvsLS_dphi->Add(RLSsubhPhi2Dpeak_dphi, -1.0);
    resUSvsLS_dphi->Divide(AvgUSsubhPhi2Dpeak_dphi);



    TFile* output = new TFile(Form("US_syst_%s", inputFile.c_str()), "RECREATE");
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
    AvgUSsubhPhi2Dpeakrightscale->Write();
    AvgUSsubhPhi2Dpeakrightscale_deta->Write();
    AvgUSsubhPhi2Dpeakrightscale_dphi->Write();
    AvgUSsubhPhi2Dpeakrightscaleplus->Write();
    AvgUSsubhPhi2Dpeakrightscaleplus_deta->Write();
    AvgUSsubhPhi2Dpeakrightscaleplus_dphi->Write();
    AvgUSsubhPhi2Dpeakrightscaleminus->Write();
    AvgUSsubhPhi2Dpeakrightscaleminus_deta->Write();
    AvgUSsubhPhi2Dpeakrightscaleminus_dphi->Write();
    AvgUSsubhPhi2Dpeakleftscale->Write();
    AvgUSsubhPhi2Dpeakleftscale_deta->Write();
    AvgUSsubhPhi2Dpeakleftscale_dphi->Write();
    AvgUSsubhPhi2Dpeakleftscaleplus->Write();
    AvgUSsubhPhi2Dpeakleftscaleplus_deta->Write();
    AvgUSsubhPhi2Dpeakleftscaleplus_dphi->Write();
    AvgUSsubhPhi2Dpeakleftscaleminus->Write();
    AvgUSsubhPhi2Dpeakleftscaleminus_deta->Write();
    AvgUSsubhPhi2Dpeakleftscaleminus_dphi->Write();
    AvgUSsubhPhi2Dpeakavgscale->Write();
    AvgUSsubhPhi2Dpeakavgscale_deta->Write();
    AvgUSsubhPhi2Dpeakavgscale_dphi->Write();
    AvgUSsubhPhi2Dpeakavgscaleplus->Write();
    AvgUSsubhPhi2Dpeakavgscaleplus_deta->Write();
    AvgUSsubhPhi2Dpeakavgscaleplus_dphi->Write();
    AvgUSsubhPhi2Dpeakavgscaleminus->Write();
    AvgUSsubhPhi2Dpeakavgscaleminus_deta->Write();
    AvgUSsubhPhi2Dpeakavgscaleminus_dphi->Write();
    RSUSsubhPhi2Dpeak->Write();
    RSUSsubhPhi2Dpeak_deta->Write();
    RSUSsubhPhi2Dpeak_dphi->Write();
    RSUSsubhPhi2Dpeakleftscale->Write();
    RSUSsubhPhi2Dpeakleftscale_deta->Write();
    RSUSsubhPhi2Dpeakleftscale_dphi->Write();
    RSUSsubhPhi2Dpeakavgscale->Write();
    RSUSsubhPhi2Dpeakavgscale_deta->Write();
    RSUSsubhPhi2Dpeakavgscale_dphi->Write();
    LSUSsubhPhi2Dpeak->Write();
    LSUSsubhPhi2Dpeak_deta->Write();
    LSUSsubhPhi2Dpeak_dphi->Write();
    LSUSsubhPhi2Dpeakleftscale->Write();
    LSUSsubhPhi2Dpeakleftscale_deta->Write();
    LSUSsubhPhi2Dpeakleftscale_dphi->Write();
    LSUSsubhPhi2Dpeakavgscale->Write();
    LSUSsubhPhi2Dpeakavgscale_deta->Write();
    LSUSsubhPhi2Dpeakavgscale_dphi->Write();
    resUSvsLS->Write();
    resUSvsLS_deta->Write();
    resUSvsLS_dphi->Write();
    corrMass->Write();
}
