intRatioPlot(){
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
    corrFit->FixParameter(6, 0.3333*(hhdphi_0_20->GetBinContent(8)+hhdphi_0_20->GetBinContent(16)+hhdphi_0_20->GetBinContent(1)));
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
    corrFit2->FixParameter(6, 0.3333*(hPhidphi_0_20->GetBinContent(8)+hPhidphi_0_20->GetBinContent(16)+hPhidphi_0_20->GetBinContent(1)));
    //corrFit->SetParLimits(6, hPhidphi_0_20->GetBinContent(8)*0.9, hPhidphi_0_20->GetBinContent(8)*1.1);
    corrFit2->SetParameter(0, hPhidphi_0_20->GetBinContent(hPhidphi_0_20->GetXaxis()->FindBin(0)) - corrFit2->GetParameter(6));
    corrFit2->SetParLimits(0, corrFit2->GetParameter(0)*0.5, corrFit2->GetParameter(0)*1.5);
    corrFit2->SetParameter(1, 0.0);
    corrFit2->SetParLimits(1, -0.5, 0.5);
    corrFit2->SetParameter(2, 0.5);
    corrFit2->SetParLimits(2, 0.2, 0.9);
    corrFit2->SetParameter(3, hPhidphi_0_20->GetBinContent(hPhidphi_0_20->GetXaxis()->FindBin(3.14)) - corrFit2->GetParameter(6));
    corrFit2->SetParLimits(3, corrFit2->GetParameter(3)*0.5, corrFit2->GetParameter(3)*1.5);
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
    hhBG->SetParameter(0, 1.0*corrFit->GetParameter(6));
    hhBG->SetLineStyle(2);

    TF1 *hphiBG = new TF1("hphiBG", "pol0(0)", -1.4, 4.6);
    hphiBG->SetParLimits(0, 0.00001, 10000000.0);
    hphiBG->SetParameter(0, 1.0*corrFit2->GetParameter(6));
    hphiBG->SetLineStyle(2);

    float near0_20hPhiError = 0;
    float near0_20hhError = 0;
    float away0_20hPhiError = 0;
    float away0_20hhError = 0;
    float total0_20hPhiError = 0;
    float total0_20hhError = 0;
    float mid0_20hPhiError =(1.0/3.0)*TMath::Sqrt(TMath::Power(hPhidphi_0_20->GetBinError(8),2) + TMath::Power(hPhidphi_0_20->GetBinError(16),2) + TMath::Power(hPhidphi_0_20->GetBinError(1),2));
    float mid0_20hhError = (1.0/3.0)*TMath::Sqrt(TMath::Power(hhdphi_0_20->GetBinError(8),2) + TMath::Power(hhdphi_0_20->GetBinError(16),2) + TMath::Power(hhdphi_0_20->GetBinError(1),2));
    
    float near0_20hPhiYield = hPhidphi_0_20->IntegralAndError(2,7,near0_20hPhiError) - hphiBG->GetParameter(0)*6.0;
    float near0_20hPhiError = TMath::Sqrt(TMath::Power(near0_20hPhiError, 2) + TMath::Power(6.0*mid0_20hPhiError, 2));
    float near0_20hhYield = hhdphi_0_20->IntegralAndError(2,7,near0_20hhError) - hhBG->GetParameter(0)*6.0;
    float near0_20hhError = TMath::Sqrt(TMath::Power(near0_20hhError, 2) + TMath::Power(6.0*mid0_20hhError, 2));
    float away0_20hPhiYield = hPhidphi_0_20->IntegralAndError(9,16,away0_20hPhiError) - hphiBG->GetParameter(0)*8.0;
    float away0_20hPhiError = TMath::Sqrt(TMath::Power(away0_20hPhiError, 2) + TMath::Power(8.0*mid0_20hPhiError, 2));
    float away0_20hhYield = hhdphi_0_20->IntegralAndError(9,16,away0_20hPhiError)- hhBG->GetParameter(0)*8.0;
    float away0_20hhError = TMath::Sqrt(TMath::Power(away0_20hhError, 2) + TMath::Power(8.0*mid0_20hhError, 2));
    float mid0_20hPhiYield = hphiBG->GetParameter(0)*16.0;
    float mid0_20hPhiError = mid0_20hPhiError*16.0;
    float mid0_20hhYield = hhBG->GetParameter(0)*16.0;
    float mid0_20hhError = mid0_20hhError*16.0;
    float total0_20hPhiYield = hPhidphi_0_20->IntegralAndError(1,16,total0_20hPhiError);
    float total0_20hhYield = hhdphi_0_20->IntegralAndError(1,16,total0_20hhError);

    float near020 = near0_20hPhiYield/near0_20hhYield;
    float near020Er = near020*TMath::Sqrt(TMath::Power(near0_20hPhiError/near0_20hPhiYield, 2) + TMath::Power(near0_20hhError/near0_20hhYield, 2));
    float away020 = away0_20hPhiYield/away0_20hhYield;
    float away020Er = away020*TMath::Sqrt(TMath::Power(away0_20hPhiError/away0_20hPhiYield, 2) + TMath::Power(away0_20hhError/away0_20hhYield, 2));
    float mid020 = mid0_20hPhiYield/mid0_20hhYield;
    float mid020Er = mid020*TMath::Sqrt(TMath::Power(mid0_20hPhiError/mid0_20hPhiYield, 2) + TMath::Power(mid0_20hhError/mid0_20hhYield, 2));
    float total020 = total0_20hPhiYield/total0_20hhYield;
    float total020Er = total020*TMath::Sqrt(TMath::Power(total0_20hPhiError/total0_20hPhiYield, 2) + TMath::Power(total0_20hhError/total0_20hhYield, 2));

    TH1D *ratios020 = new TH1D("ratios020", "(h-#phi / h-h) Ratios", 4, 0, 4);
    ratios020->GetXaxis()->SetBinLabel(1, "near-side");
    ratios020->SetBinContent(1, near020);
    ratios020->SetBinError(1, near020Er);
    ratios020->GetXaxis()->SetBinLabel(2, "mid");
    ratios020->SetBinContent(2, mid020);
    ratios020->SetBinError(2, mid020Er);
    ratios020->GetXaxis()->SetBinLabel(3, "away-side");
    ratios020->SetBinContent(3, away020);
    ratios020->SetBinError(3, away020Er);
    ratios020->GetXaxis()->SetBinLabel(4, "total");
    ratios020->SetBinContent(4, total020);
    ratios020->SetBinError(4, total020Er);
    ratios020->SetMarkerStyle(22);
    ratios020->SetMarkerColor(kRed+2);
    ratios020->SetLineColor(kRed+2);
    ratios020->SetMarkerSize(2);
    ratios020->SetLineWidth(2);

   //20-50 section 
    TFile *hhFile = new TFile("~/phiStudies/LHC16q_FAST_hh_20_50/trig_4_8_assoc_2_4_hh_phiCorrelations_mult_20_50.root");
    TH2D* hh2D_20_50 = (TH2D*)hh2D->Clone("hhdphi");
    TH1D* hhdphi_20_50 = (TH1D*)hh2D_20_50->ProjectionY("hhdphi_50_100", hh2D_20_50->GetXaxis()->FindBin(-1.2), hh2D_20_50->GetXaxis()->FindBin(1.2));
    hhdphi_20_50->Rebin();
    hhdphi_20_50->SetLineWidth(4);
    hhdphi_20_50->SetLineColor(kBlue+2);
    hhdphi_20_50->SetMarkerColor(kBlue+2);
    hhdphi_20_50->SetMarkerStyle(21);
    hhdphi_20_50->SetMarkerSize(2);
    hhdphi_20_50->GetXaxis()->SetTitle("#Delta#varphi");
    hhdphi_20_50->SetTitle("");
    hhdphi_20_50->GetYaxis()->SetTitleOffset(1.60);
    hhdphi_20_50->GetXaxis()->SetTitleSize(0.05);
    hhdphi_20_50->GetXaxis()->SetTitleOffset(0.90);


    TFile* phiFile = new TFile("~/phiStudies/LHC16q_FAST_20_50_widecuts/LS_trig_4_8_assoc_2_4_mixcorr_phiCorrelations_mult_20_50.root");    
    TH2D* hPhi2D_20_50 = (TH2D*)RLSsubhPhi2Dpeak->Clone("hPhidphi");
    TH1D* hPhidphi_20_50 = (TH1D*)hPhi2D_20_50->ProjectionY("hPhidphi_50_100", hPhi2D_20_50->GetXaxis()->FindBin(-1.2), hPhi2D_20_50->GetXaxis()->FindBin(1.2));
    //hPhidphi_20_50->Rebin();
    hPhidphi_20_50->SetLineWidth(4);
    hPhidphi_20_50->SetLineColor(kRed+2);
    hPhidphi_20_50->SetMarkerColor(kRed+2);
    hPhidphi_20_50->SetMarkerStyle(22);
    hPhidphi_20_50->SetMarkerSize(2);
    hPhidphi_20_50->GetXaxis()->SetTitle("#Delta#varphi");
    hPhidphi_20_50->SetTitle("");
    hPhidphi_20_50->GetYaxis()->SetTitleOffset(1.20);
    hPhidphi_20_50->GetXaxis()->SetTitleSize(0.05);
    hPhidphi_20_50->GetXaxis()->SetTitleOffset(0.90);

    TF1 *corrFit2050 = new TF1("corrFit2050", "gaus(0) + gaus(3) + pol0(6)", -1.4, 4.6);
    corrFit2050->FixParameter(6, 0.3333*(hhdphi_20_50->GetBinContent(8)+hhdphi_20_50->GetBinContent(16)+hhdphi_20_50->GetBinContent(1)));
    //corrFit2050->SetParLimits(6, hhdphi_20_50->GetBinContent(8)*0.9, hhdphi_20_50->GetBinContent(8)*1.1);
    corrFit2050->SetParameter(0, hhdphi_20_50->GetBinContent(hhdphi_20_50->GetXaxis()->FindBin(0)) - corrFit2050->GetParameter(6));
    corrFit2050->SetParLimits(0, corrFit2050->GetParameter(0)*0.5, corrFit2050->GetParameter(0)*1.5);
    corrFit2050->SetParameter(1, 0.0);
    corrFit2050->SetParLimits(1, -0.5, 0.5);
    corrFit2050->SetParameter(2, 0.5);
    corrFit2050->SetParLimits(2, 0.2, 0.9);
    corrFit2050->SetParameter(3, hhdphi_20_50->GetBinContent(hhdphi_20_50->GetXaxis()->FindBin(3.14)) - corrFit2050->GetParameter(6));
    corrFit2050->SetParLimits(3, corrFit2050->GetParameter(3)*0.5, corrFit2050->GetParameter(3)*1.5);
    corrFit2050->SetParameter(4, 3.14);
    corrFit2050->SetParLimits(4, 3.0, 3.25);
    corrFit2050->SetParameter(5, 0.5);
    corrFit2050->SetParLimits(5, 0.2, 0.9);

    corrFit2050->SetParNames("peak1", "mean1", "sigma1", "peak2", "mean2", "sigma2", "BG");

    corrFit2050->SetLineColor(kBlue);
    corrFit2050->SetLineWidth(6);
    corrFit2050->SetLineStyle(7);

    TF1 *corrFit2_2050 = new TF1("corrFit2_2050", "gaus(0) + gaus(3) + pol0(6)", -1.4, 4.6);
    corrFit2_2050->FixParameter(6, 0.3333*(hPhidphi_20_50->GetBinContent(8)+hPhidphi_20_50->GetBinContent(16)+hPhidphi_20_50->GetBinContent(1)));
    //corrFit->SetParLimits(6, hPhidphi_20_50->GetBinContent(8)*0.9, hPhidphi_20_50->GetBinContent(8)*1.1);
    corrFit2_2050->SetParameter(0, hPhidphi_20_50->GetBinContent(hPhidphi_20_50->GetXaxis()->FindBin(0)) - corrFit2_2050->GetParameter(6));
    corrFit2_2050->SetParLimits(0, corrFit2_2050->GetParameter(0)*0.5, corrFit2_2050->GetParameter(0)*1.5);
    corrFit2_2050->SetParameter(1, 0.0);
    corrFit2_2050->SetParLimits(1, -0.5, 0.5);
    corrFit2_2050->SetParameter(2, 0.5);
    corrFit2_2050->SetParLimits(2, 0.2, 0.9);
    corrFit2_2050->SetParameter(3, hPhidphi_20_50->GetBinContent(hPhidphi_20_50->GetXaxis()->FindBin(3.14)) - corrFit2_2050->GetParameter(6));
    corrFit2_2050->SetParLimits(3, corrFit2_2050->GetParameter(3)*0.5, corrFit2_2050->GetParameter(3)*1.5);
    corrFit2_2050->SetParameter(4, 3.14);
    corrFit2_2050->SetParLimits(4, 3.0, 3.25);
    corrFit2_2050->SetParameter(5, 0.5);
    corrFit2_2050->SetParLimits(5, 0.1, 0.9);

    corrFit2_2050->SetParNames("peak1", "mean1", "sigma1", "peak2", "mean2", "sigma2", "BG");

    corrFit2_2050->SetLineColor(kRed);
    corrFit2_2050->SetLineWidth(6);
    corrFit2_2050->SetLineStyle(7);

    hhdphi_20_50->Fit("corrFit2050", "R0");
    hPhidphi_20_50->Fit("corrFit2_2050", "R0");

    TF1 *hhBG_20_50 = new TF1("hhBG_20_50", "pol0(0)", -1.4, 4.6);
    hhBG_20_50->SetParLimits(0, 0.00001, 10000000.0);
    hhBG_20_50->SetParameter(0, 1.0*corrFit2050->GetParameter(6));
    hhBG_20_50->SetLineStyle(2);

    TF1 *hphiBG_20_50 = new TF1("hphiBG_20_50", "pol0(0)", -1.4, 4.6);
    hphiBG_20_50->SetParLimits(0, 0.00001, 10000000.0);
    hphiBG_20_50->SetParameter(0, 1.0*corrFit2_2050->GetParameter(6));
    hphiBG_20_50->SetLineStyle(2);

    float near20_50hPhiError = 0;
    float near20_50hhError = 0;
    float away20_50hPhiError = 0;
    float away20_50hhError = 0;
    float mid20_50hPhiError =(1.0/3.0)*TMath::Sqrt(TMath::Power(hPhidphi_20_50->GetBinError(8),2) + TMath::Power(hPhidphi_20_50->GetBinError(16),2) + TMath::Power(hPhidphi_20_50->GetBinError(1),2));
    float mid20_50hhError = (1.0/3.0)*TMath::Sqrt(TMath::Power(hhdphi_20_50->GetBinError(8),2) + TMath::Power(hhdphi_20_50->GetBinError(16),2) + TMath::Power(hhdphi_20_50->GetBinError(1),2));
    float total20_50hPhiError = 0;
    float total20_50hhError = 0;

    float near20_50hPhiYield = hPhidphi_20_50->IntegralAndError(2,7,near20_50hPhiError) - hphiBG_20_50->GetParameter(0)*6.0;
    float near20_50hPhiError = TMath::Sqrt(TMath::Power(near20_50hPhiError, 2) + TMath::Power(6.0*mid20_50hPhiError, 2));
    float near20_50hhYield = hhdphi_20_50->IntegralAndError(2,7,near20_50hhError) - hhBG_20_50->GetParameter(0)*6.0;
    float near20_50hhError = TMath::Sqrt(TMath::Power(near20_50hhError, 2) + TMath::Power(6.0*mid20_50hhError, 2));
    float away20_50hPhiYield = hPhidphi_20_50->IntegralAndError(9,16,away20_50hPhiError) - hphiBG_20_50->GetParameter(0)*8.0;
    float away20_50hPhiError = TMath::Sqrt(TMath::Power(away20_50hPhiError, 2) + TMath::Power(8.0*mid20_50hPhiError, 2));
    float away20_50hhYield = hhdphi_20_50->IntegralAndError(9,16,away20_50hPhiError)- hhBG_20_50->GetParameter(0)*8.0;
    float away20_50hhError = TMath::Sqrt(TMath::Power(away20_50hhError, 2) + TMath::Power(8.0*mid20_50hhError, 2));
    float mid20_50hPhiYield = hphiBG_20_50->GetParameter(0)*16.0;
    float mid20_50hPhiError = mid20_50hPhiError*16.0;
    float mid20_50hhYield = hhBG_20_50->GetParameter(0)*16.0;
    float mid20_50hhError = mid20_50hhError*16.0;
    float total20_50hPhiYield = hPhidphi_20_50->IntegralAndError(1, 16,total20_50hPhiError);
    float total20_50hhYield = hhdphi_20_50->IntegralAndError(1, 16,total20_50hhError);

    float near2050 = near20_50hPhiYield/near20_50hhYield;
    float near2050Er = near2050*TMath::Sqrt(TMath::Power(near20_50hPhiError/near20_50hPhiYield, 2) + TMath::Power(near20_50hhError/near20_50hhYield, 2));
    float away2050 = away20_50hPhiYield/away20_50hhYield;
    float away2050Er = away2050*TMath::Sqrt(TMath::Power(away20_50hPhiError/away20_50hPhiYield, 2) + TMath::Power(away20_50hhError/away20_50hhYield, 2));
    float mid2050 = mid20_50hPhiYield/mid20_50hhYield;
    float mid2050Er = mid2050*TMath::Sqrt(TMath::Power(mid20_50hPhiError/mid20_50hPhiYield, 2) + TMath::Power(mid20_50hhError/mid20_50hhYield, 2));
    float total2050 = total20_50hPhiYield/total20_50hhYield;
    float total2050Er = total2050*TMath::Sqrt(TMath::Power(total20_50hPhiError/total20_50hPhiYield, 2) + TMath::Power(total20_50hhError/total20_50hhYield, 2));

    TH1D *ratios2050 = new TH1D("ratios2050", "(h-#phi / h-h) Ratios", 4, 0, 4);
    ratios2050->GetXaxis()->SetBinLabel(1, "near-side");
    ratios2050->SetBinContent(1, near2050);
    ratios2050->SetBinError(1, near2050Er);
    ratios2050->GetXaxis()->SetBinLabel(2, "mid");
    ratios2050->SetBinContent(2, mid2050);
    ratios2050->SetBinError(2, mid2050Er);
    ratios2050->GetXaxis()->SetBinLabel(3, "away-side");
    ratios2050->SetBinContent(3, away2050);
    ratios2050->SetBinError(3, away2050Er);
    ratios2050->GetXaxis()->SetBinLabel(4, "total");
    ratios2050->SetBinContent(4, total2050);
    ratios2050->SetBinError(4, total2050Er);
    ratios2050->SetMarkerStyle(21);
    ratios2050->SetMarkerColor(kBlue+2);
    ratios2050->SetLineColor(kBlue+2);
    ratios2050->SetMarkerSize(2);
    ratios2050->SetLineWidth(2);


    //50-100 section
    TFile *hhFile = new TFile("~/phiStudies/LHC16q_FAST_50_100_widecuts/trig_4_8_assoc_2_4_hh_phiCorrelations_mult_50_100.root");
    TH2D* hh2D_50_100 = (TH2D*)hh2D->Clone("hhdphi");
    TH1D* hhdphi_50_100 = (TH1D*)hh2D_50_100->ProjectionY("hhdphi_50_100", hh2D_50_100->GetXaxis()->FindBin(-1.2), hh2D_50_100->GetXaxis()->FindBin(1.2));
    hhdphi_50_100->Rebin();
    hhdphi_50_100->SetLineWidth(4);
    hhdphi_50_100->SetLineColor(kBlue+2);
    hhdphi_50_100->SetMarkerColor(kBlue+2);
    hhdphi_50_100->SetMarkerStyle(21);
    hhdphi_50_100->SetMarkerSize(2);
    hhdphi_50_100->GetXaxis()->SetTitle("#Delta#varphi");
    hhdphi_50_100->SetTitle("");
    hhdphi_50_100->GetYaxis()->SetTitleOffset(1.60);
    hhdphi_50_100->GetXaxis()->SetTitleSize(0.05);
    hhdphi_50_100->GetXaxis()->SetTitleOffset(0.90);


    TFile* phiFile = new TFile("~/phiStudies/LHC16q_FAST_50_100_widecuts/LS_trig_4_8_assoc_2_4_mixcorr_phiCorrelations_mult_50_100.root");    
    TH2D* hPhi2D_50_100 = (TH2D*)RLSsubhPhi2Dpeak->Clone("hPhidphi");
    TH1D* hPhidphi_50_100 = (TH1D*)hPhi2D_50_100->ProjectionY("hPhidphi_50_100", hPhi2D_50_100->GetXaxis()->FindBin(-1.2), hPhi2D_50_100->GetXaxis()->FindBin(1.2));
    //hPhidphi_50_100->Rebin();
    hPhidphi_50_100->SetLineWidth(4);
    hPhidphi_50_100->SetLineColor(kRed+2);
    hPhidphi_50_100->SetMarkerColor(kRed+2);
    hPhidphi_50_100->SetMarkerStyle(22);
    hPhidphi_50_100->SetMarkerSize(2);
    hPhidphi_50_100->GetXaxis()->SetTitle("#Delta#varphi");
    hPhidphi_50_100->SetTitle("");
    hPhidphi_50_100->GetYaxis()->SetTitleOffset(1.20);
    hPhidphi_50_100->GetXaxis()->SetTitleSize(0.05);
    hPhidphi_50_100->GetXaxis()->SetTitleOffset(0.90);

    TF1 *corrFit50100 = new TF1("corrFit50100", "gaus(0) + gaus(3) + pol0(6)", -1.4, 4.6);
    corrFit50100->FixParameter(6, 0.3333*(hhdphi_50_100->GetBinContent(8)+hhdphi_50_100->GetBinContent(16)+hhdphi_50_100->GetBinContent(1)));
    //corrFit50100->SetParLimits(6, hhdphi_50_100->GetBinContent(8)*0.9, hhdphi_50_100->GetBinContent(8)*1.1);
    corrFit50100->SetParameter(0, hhdphi_50_100->GetBinContent(hhdphi_50_100->GetXaxis()->FindBin(0)) - corrFit50100->GetParameter(6));
    corrFit50100->SetParLimits(0, corrFit50100->GetParameter(0)*0.5, corrFit50100->GetParameter(0)*1.5);
    corrFit50100->SetParameter(1, 0.0);
    corrFit50100->SetParLimits(1, -0.5, 0.5);
    corrFit50100->SetParameter(2, 0.5);
    corrFit50100->SetParLimits(2, 0.2, 0.9);
    corrFit50100->SetParameter(3, hhdphi_50_100->GetBinContent(hhdphi_50_100->GetXaxis()->FindBin(3.14)) - corrFit50100->GetParameter(6));
    corrFit50100->SetParLimits(3, corrFit50100->GetParameter(3)*0.5, corrFit50100->GetParameter(3)*1.5);
    corrFit50100->SetParameter(4, 3.14);
    corrFit50100->SetParLimits(4, 3.0, 3.25);
    corrFit50100->SetParameter(5, 0.5);
    corrFit50100->SetParLimits(5, 0.2, 0.9);

    corrFit50100->SetParNames("peak1", "mean1", "sigma1", "peak2", "mean2", "sigma2", "BG");

    corrFit50100->SetLineColor(kBlue);
    corrFit50100->SetLineWidth(6);
    corrFit50100->SetLineStyle(7);

    TF1 *corrFit2_50100 = new TF1("corrFit2_50100", "gaus(0) + gaus(3) + pol0(6)", -1.4, 4.6);
    corrFit2_50100->FixParameter(6, 0.3333*(hPhidphi_50_100->GetBinContent(8)+hPhidphi_50_100->GetBinContent(16)+hPhidphi_50_100->GetBinContent(1)));
    //corrFit->SetParLimits(6, hPhidphi_50_100->GetBinContent(8)*0.9, hPhidphi_50_100->GetBinContent(8)*1.1);
    corrFit2_50100->SetParameter(0, hPhidphi_50_100->GetBinContent(hPhidphi_50_100->GetXaxis()->FindBin(0)) - corrFit2_50100->GetParameter(6));
    corrFit2_50100->SetParLimits(0, corrFit2_50100->GetParameter(0)*0.5, corrFit2_50100->GetParameter(0)*1.5);
    corrFit2_50100->SetParameter(1, 0.0);
    corrFit2_50100->SetParLimits(1, -0.5, 0.5);
    corrFit2_50100->SetParameter(2, 0.5);
    corrFit2_50100->SetParLimits(2, 0.2, 0.9);
    corrFit2_50100->SetParameter(3, hPhidphi_50_100->GetBinContent(hPhidphi_50_100->GetXaxis()->FindBin(3.14)) - corrFit2_50100->GetParameter(6));
    corrFit2_50100->SetParLimits(3, corrFit2_50100->GetParameter(3)*0.5, corrFit2_50100->GetParameter(3)*1.5);
    corrFit2_50100->SetParameter(4, 3.14);
    corrFit2_50100->SetParLimits(4, 3.0, 3.25);
    corrFit2_50100->SetParameter(5, 0.5);
    corrFit2_50100->SetParLimits(5, 0.1, 0.9);

    corrFit2_50100->SetParNames("peak1", "mean1", "sigma1", "peak2", "mean2", "sigma2", "BG");

    corrFit2_50100->SetLineColor(kRed);
    corrFit2_50100->SetLineWidth(6);
    corrFit2_50100->SetLineStyle(7);

    hhdphi_50_100->Fit("corrFit50100", "R0");
    hPhidphi_50_100->Fit("corrFit2_50100", "R0");

    TF1 *hhBG_50_100 = new TF1("hhBG_50_100", "pol0(0)", -1.4, 4.6);
    hhBG_50_100->SetParLimits(0, 0.00001, 10000000.0);
    hhBG_50_100->SetParameter(0, 1.0*corrFit50100->GetParameter(6));
    hhBG_50_100->SetLineStyle(2);

    TF1 *hphiBG_50_100 = new TF1("hphiBG_50_100", "pol0(0)", -1.4, 4.6);
    hphiBG_50_100->SetParLimits(0, 0.00001, 10000000.0);
    hphiBG_50_100->SetParameter(0, 1.0*corrFit2_50100->GetParameter(6));
    hphiBG_50_100->SetLineStyle(2);

    float near50_100hPhiError = 0;
    float near50_100hhError = 0;
    float away50_100hPhiError = 0;
    float away50_100hhError = 0;
    float mid50_100hPhiError =(1.0/3.0)*TMath::Sqrt(TMath::Power(hPhidphi_50_100->GetBinError(8),2) + TMath::Power(hPhidphi_50_100->GetBinError(16),2) + TMath::Power(hPhidphi_50_100->GetBinError(1),2));
    float mid50_100hhError = (1.0/3.0)*TMath::Sqrt(TMath::Power(hhdphi_50_100->GetBinError(8),2) + TMath::Power(hhdphi_50_100->GetBinError(16),2) + TMath::Power(hhdphi_50_100->GetBinError(1),2));
    float total50_100hPhiError = 0;
    float total50_100hhError = 0;

    float near50_100hPhiYield = hPhidphi_50_100->IntegralAndError(2,7,near50_100hPhiError) - hphiBG_50_100->GetParameter(0)*6.0;
    float near50_100hPhiError = TMath::Sqrt(TMath::Power(near50_100hPhiError, 2) + TMath::Power(6.0*mid50_100hPhiError, 2));
    float near50_100hhYield = hhdphi_50_100->IntegralAndError(2,7,near50_100hhError) - hhBG_50_100->GetParameter(0)*6.0;
    float near50_100hhError = TMath::Sqrt(TMath::Power(near50_100hhError, 2) + TMath::Power(6.0*mid50_100hhError, 2));
    float away50_100hPhiYield = hPhidphi_50_100->IntegralAndError(9,16,away50_100hPhiError) - hphiBG_50_100->GetParameter(0)*8.0;
    float away50_100hPhiError = TMath::Sqrt(TMath::Power(away50_100hPhiError, 2) + TMath::Power(8.0*mid50_100hPhiError, 2));
    float away50_100hhYield = hhdphi_50_100->IntegralAndError(9,16,away50_100hPhiError)- hhBG_50_100->GetParameter(0)*8.0;
    float away50_100hhError = TMath::Sqrt(TMath::Power(away50_100hhError, 2) + TMath::Power(8.0*mid50_100hhError, 2));
    float mid50_100hPhiYield = hphiBG_50_100->GetParameter(0)*16.0;
    float mid50_100hPhiError = mid50_100hPhiError*16.0;
    float mid50_100hhYield = hhBG_50_100->GetParameter(0)*16.0;
    float mid50_100hhError = mid50_100hhError*16.0;
    float total50_100hPhiYield = hPhidphi_50_100->IntegralAndError(1, 16,total50_100hPhiError);
    float total50_100hhYield = hhdphi_50_100->IntegralAndError(1, 16,total50_100hhError);

    float near50100 = near50_100hPhiYield/near50_100hhYield;
    float near50100Er = near50100*TMath::Sqrt(TMath::Power(near50_100hPhiError/near50_100hPhiYield, 2) + TMath::Power(near50_100hhError/near50_100hhYield, 2));
    float away50100 = away50_100hPhiYield/away50_100hhYield;
    float away50100Er = away50100*TMath::Sqrt(TMath::Power(away50_100hPhiError/away50_100hPhiYield, 2) + TMath::Power(away50_100hhError/away50_100hhYield, 2));
    float mid50100 = mid50_100hPhiYield/mid50_100hhYield;
    float mid50100Er = mid50100*TMath::Sqrt(TMath::Power(mid50_100hPhiError/mid50_100hPhiYield, 2) + TMath::Power(mid50_100hhError/mid50_100hhYield, 2));
    float total50100 = total50_100hPhiYield/total50_100hhYield;
    float total50100Er = total50100*TMath::Sqrt(TMath::Power(total50_100hPhiError/total50_100hPhiYield, 2) + TMath::Power(total50_100hhError/total50_100hhYield, 2));

    TH1D *ratios50100 = new TH1D("ratios50100", "(h-#phi / h-h) Ratios", 4, 0, 4);
    ratios50100->GetXaxis()->SetBinLabel(1, "near-side");
    ratios50100->SetBinContent(1, near50100);
    ratios50100->SetBinError(1, near50100Er);
    ratios50100->GetXaxis()->SetBinLabel(2, "mid");
    ratios50100->SetBinContent(2, mid50100);
    ratios50100->SetBinError(2, mid50100Er);
    ratios50100->GetXaxis()->SetBinLabel(3, "away-side");
    ratios50100->SetBinContent(3, away50100);
    ratios50100->SetBinError(3, away50100Er);
    ratios50100->GetXaxis()->SetBinLabel(4, "total");
    ratios50100->SetBinContent(4, total50100);
    ratios50100->SetBinError(4, total50100Er);
    ratios50100->SetMarkerStyle(20);
    ratios50100->SetMarkerColor(kGreen+2);
    ratios50100->SetLineColor(kGreen+2);
    ratios50100->SetMarkerSize(2);
    ratios50100->SetLineWidth(2);

    float near0100 = (near0_20hPhiYield + near20_50hPhiYield + near50_100hPhiYield)/(near0_20hhYield + near20_50hhYield + near50_100hhYield);
    float near0100Er = near0100*TMath::Sqrt((TMath::Power(near0_20hPhiError,2) + TMath::Power(near20_50hPhiError,2) + TMath::Power(near50_100hPhiError,2))/TMath::Power((near0_20hPhiYield + near20_50hPhiYield + near50_100hPhiYield), 2) + (TMath::Power(near0_20hhError,2) + TMath::Power(near20_50hhError,2) + TMath::Power(near50_100hhError,2))/TMath::Power((near0_20hhYield + near20_50hhYield + near50_100hhYield), 2));
    float away0100 = (away0_20hPhiYield + away20_50hPhiYield + away50_100hPhiYield)/(away0_20hhYield + away20_50hhYield + away50_100hhYield);
    float away0100Er = away0100*TMath::Sqrt((TMath::Power(away0_20hPhiError,2) + TMath::Power(away20_50hPhiError,2) + TMath::Power(away50_100hPhiError,2))/TMath::Power((away0_20hPhiYield + away20_50hPhiYield + away50_100hPhiYield), 2) + (TMath::Power(away0_20hhError,2) + TMath::Power(away20_50hhError,2) + TMath::Power(away50_100hhError,2))/TMath::Power((away0_20hhYield + away20_50hhYield + away50_100hhYield), 2));
    float mid0100 = (mid0_20hPhiYield + mid20_50hPhiYield + mid50_100hPhiYield)/(mid0_20hhYield + mid20_50hhYield + mid50_100hhYield);
    float mid0100Er = mid0100*TMath::Sqrt((TMath::Power(mid0_20hPhiError,2) + TMath::Power(mid20_50hPhiError,2) + TMath::Power(mid50_100hPhiError,2))/TMath::Power((mid0_20hPhiYield + mid20_50hPhiYield + mid50_100hPhiYield), 2) + (TMath::Power(mid0_20hhError,2) + TMath::Power(mid20_50hhError,2) + TMath::Power(mid50_100hhError,2))/TMath::Power((mid0_20hhYield + mid20_50hhYield + mid50_100hhYield), 2));
    float total0100 = (total0_20hPhiYield + total20_50hPhiYield + total50_100hPhiYield)/(total0_20hhYield + total20_50hhYield + total50_100hhYield);
    float total0100Er = total0100*TMath::Sqrt((TMath::Power(total0_20hPhiError,2) + TMath::Power(total20_50hPhiError,2) + TMath::Power(total50_100hPhiError,2))/TMath::Power((total0_20hPhiYield + total20_50hPhiYield + total50_100hPhiYield), 2) + (TMath::Power(total0_20hhError,2) + TMath::Power(total20_50hhError,2) + TMath::Power(total50_100hhError,2))/TMath::Power((total0_20hhYield + total20_50hhYield + total50_100hhYield), 2));
    
    printf("near0100Er: %E\n", near0100Er);
    printf("total0100Er: %E\n", total0100Er);
    printf("total020Er: %E\n", total020Er);
    printf("total0_20hPhiError: %E\n", total0_20hPhiError);
    TH1D *ratios0100 = new TH1D("ratios0100", "(h-#phi / h-h) Ratios", 4, 0, 4);
    ratios0100->GetXaxis()->SetBinLabel(1, "near-side");
    ratios0100->SetBinContent(1, near0100);
    ratios0100->SetBinError(1, near0100Er);
    ratios0100->GetXaxis()->SetBinLabel(2, "mid");
    ratios0100->SetBinContent(2, mid0100);
    ratios0100->SetBinError(2, mid0100Er);
    ratios0100->GetXaxis()->SetBinLabel(3, "away-side");
    ratios0100->SetBinContent(3, away0100);
    ratios0100->SetBinError(3, away0100Er);
    ratios0100->GetXaxis()->SetBinLabel(4, "total");
    ratios0100->SetBinContent(4, total0100);
    ratios0100->SetBinError(4, total0100Er);
    ratios0100->SetMarkerStyle(29);
    ratios0100->SetMarkerColor(kViolet-1);
    ratios0100->SetLineColor(kViolet-1);
    ratios0100->SetLineWidth(2);
    ratios0100->SetMarkerSize(2);

    TLegend  *ratioslegend = new TLegend(0.3791, 0.1518, 0.8772, 0.2688);
    ratioslegend->SetMargin(0.35);
    ratioslegend->AddEntry(ratios020, "0-20%", "p");
    ratioslegend->AddEntry(ratios2050, "20-50%", "p");
    ratioslegend->AddEntry(ratios50100, "50-100%", "p");
    ratioslegend->AddEntry(ratios0100, "0-100%", "p");
    
    TLine *line = new TLine(3.0, 0.0, 3.0, 0.0040);
    line->SetLineStyle(7);
    line->SetLineWidth(2);
  
    TCanvas *testc = new TCanvas("test", "test",50, 50, 600, 600);
    testc->cd();
    ratios2050->GetYaxis()->SetRangeUser(0.0000, 0.0040);
    ratios2050->GetXaxis()->SetLabelSize(0.07);
    ratios2050->Draw("P SAME");
    ratios50100->Draw("P SAME");
    ratios020->Draw("P SAME");
    ratios0100->Draw("P SAME");
    line->Draw("SAME");
    ratioslegend->Draw("SAME");

    //Double Ratio plots:
    float doublenear020 = near020/mid020;
    float doublenear020Er = doublenear020*TMath::Sqrt(TMath::Power((near020Er)/(near020), 2) + TMath::Power((mid020Er)/(mid020), 2));
    float doubleaway020 = away020/mid020;
    float doubleaway020Er = doubleaway020*TMath::Sqrt(TMath::Power((away020Er)/(away020), 2) + TMath::Power((mid020Er)/(mid020), 2));

    float doublenear2050 = near2050/mid2050;
    float doublenear2050Er = doublenear2050*TMath::Sqrt(TMath::Power((near2050Er)/(near2050), 2) + TMath::Power((mid2050Er)/(mid2050), 2));
    float doubleaway2050 = away2050/mid2050;
    float doubleaway2050Er = doubleaway2050*TMath::Sqrt(TMath::Power((away2050Er)/(away2050), 2) + TMath::Power((mid2050Er)/(mid2050), 2));

    float doublenear50100 = near50100/mid50100;
    float doublenear50100Er = doublenear50100*TMath::Sqrt(TMath::Power((near50100Er)/(near50100), 2) + TMath::Power((mid50100Er)/(mid50100), 2));
    float doubleaway50100 = away50100/mid50100;
    float doubleaway50100Er = doubleaway50100*TMath::Sqrt(TMath::Power((away50100Er)/(away50100), 2) + TMath::Power((mid50100Er)/(mid50100), 2));

    float doublenear0100 = near0100/mid0100;
    float doublenear0100Er = doublenear0100*TMath::Sqrt(TMath::Power((near0100Er)/(near0100), 2) + TMath::Power((mid0100Er)/(mid0100), 2));
    float doubleaway0100 = away0100/mid0100;
    float doubleaway0100Er = doubleaway0100*TMath::Sqrt(TMath::Power((away0100Er)/(away0100), 2) + TMath::Power((mid0100Er)/(mid0100), 2));

    TH1D *doubleratios020 = new TH1D("doubleratios020", "(h-#phi / h-h) Double Ratios", 2, 0, 2);
    doubleratios020->GetXaxis()->SetBinLabel(1, "#frac{near-side}{mid}");
    doubleratios020->SetBinContent(1, doublenear020);
    doubleratios020->SetBinError(1, doublenear020Er);
    doubleratios020->GetXaxis()->SetBinLabel(2, "#frac{away-side}{mid}");
    doubleratios020->SetBinContent(2, doubleaway020);
    doubleratios020->SetBinError(2, doubleaway020Er);
    doubleratios020->SetMarkerStyle(22);
    doubleratios020->SetMarkerColor(kRed+2);
    doubleratios020->SetLineColor(kRed+2);
    doubleratios020->SetLineWidth(2);
    doubleratios020->SetMarkerSize(2);

    TH1D *doubleratios2050 = new TH1D("doubleratios2050", "(h-#phi / h-h) Double Ratios", 2, 0, 2);
    doubleratios2050->GetXaxis()->SetBinLabel(1, "#frac{near-side}{mid}");
    doubleratios2050->SetBinContent(1, doublenear2050);
    doubleratios2050->SetBinError(1, doublenear2050Er);
    doubleratios2050->GetXaxis()->SetBinLabel(2, "#frac{away-side}{mid}");
    doubleratios2050->SetBinContent(2, doubleaway2050);
    doubleratios2050->SetBinError(2, doubleaway2050Er);
    doubleratios2050->SetMarkerStyle(21);
    doubleratios2050->SetMarkerColor(kBlue+2);
    doubleratios2050->SetLineColor(kBlue+2);
    doubleratios2050->SetLineWidth(2);
    doubleratios2050->SetMarkerSize(2);

    TH1D *doubleratios50100 = new TH1D("doubleratios50100", "(h-#phi / h-h) Double Ratios", 2, 0, 2);
    doubleratios50100->GetXaxis()->SetBinLabel(1, "#frac{near-side}{mid}");
    doubleratios50100->SetBinContent(1, doublenear50100);
    doubleratios50100->SetBinError(1, doublenear50100Er);
    doubleratios50100->GetXaxis()->SetBinLabel(2, "#frac{away-side}{mid}");
    doubleratios50100->SetBinContent(2, doubleaway50100);
    doubleratios50100->SetBinError(2, doubleaway50100Er);
    doubleratios50100->SetMarkerStyle(20);
    doubleratios50100->SetMarkerColor(kGreen+2);
    doubleratios50100->SetLineColor(kGreen+2);
    doubleratios50100->SetLineWidth(2);
    doubleratios50100->SetMarkerSize(2);

    TH1D *doubleratios0100 = new TH1D("doubleratios0100", "(h-#phi / h-h) Double Ratios", 2, 0, 2);
    doubleratios0100->GetXaxis()->SetBinLabel(1, "#frac{near-side}{mid}");
    doubleratios0100->SetBinContent(1, doublenear0100);
    doubleratios0100->SetBinError(1, doublenear0100Er);
    doubleratios0100->GetXaxis()->SetBinLabel(2, "#frac{away-side}{mid}");
    doubleratios0100->SetBinContent(2, doubleaway0100);
    doubleratios0100->SetBinError(1, doubleaway0100Er);
    doubleratios0100->SetMarkerStyle(29);
    doubleratios0100->SetMarkerColor(kViolet-1);
    doubleratios0100->SetLineColor(kViolet-1);
    doubleratios0100->SetLineWidth(2);
    doubleratios0100->SetMarkerSize(2);


    TCanvas *testDoublec = new TCanvas("testdouble", "testdouble",50, 50, 600, 600);
    testDoublec->cd();
    //doubleratios2050->GetYaxis()->SetRangeUser(0.0000, 0.0040);
    doubleratios2050->GetXaxis()->SetLabelSize(0.07);
    doubleratios2050->Draw("P SAME");
    doubleratios50100->GetXaxis()->SetLimits(-0.1, 1.9);
    doubleratios50100->Draw("P SAME");
    doubleratios020->GetXaxis()->SetLimits(0.1, 2.1);
    doubleratios020->Draw("P SAME");
    doubleratios0100->GetXaxis()->SetLimits(0.2,2.2);
    doubleratios0100->Draw("P SAME");
    //line->Draw("SAME");
    ratioslegend->Draw("SAME");


/*
    TH1D *ratio = hPhidphi_0_20->Clone("ratio");
    ratio->Divide(hhdphi_0_20);
    

    TLegend  *legend = new TLegend(0.3791, 0.1518, 0.8772, 0.2688);
    legend->SetMargin(0.15);
    legend->AddEntry(corrFit2, "Hadron-#phi(1020) Correlation", "l");
    legend->AddEntry(corrFit, "Hadron-hadron Correlations", "l");
    
    TCanvas *c0_20 = new TCanvas("c0_20", "c0_20", 50, 50, 550, 600);
    c0_20->cd();
    c0_20->SetMargin(0.12, 0.05, 0.1, 0.05);
    hhdphi_0_20->GetYaxis()->SetRangeUser(0.001, 0.25);
    hhdphi_0_20->Draw("E0 X0");
    corrFit->Draw("SAME");
    corrFit2->Draw("SAME");
    hPhidphi_0_20->Draw("E0 X0 SAME");
    hhBG->Draw("SAME");
    hphiBG->Draw("SAME");
    legend->Draw();

    TCanvas *c20_50 = new TCanvas("c20_50", "c20_50", 50, 50, 550, 600);
    c20_50->cd();
    c20_50->SetMargin(0.12, 0.05, 0.1, 0.05);
    //hhdphi_0_20->GetYaxis()->SetRangeUser(0.03, 0.11);
    hhdphi_20_50->Draw("E0 X0");
    corrFit2050->Draw("SAME");
    corrFit2_2050->Draw("SAME");
    hPhidphi_20_50->Draw("E0 X0 SAME");
    hhBG_20_50->Draw("SAME");
    hphiBG_20_50->Draw("SAME");
    legend->Draw();

    TCanvas *c50_100 = new TCanvas("c50_100", "c50_100", 50, 50, 550, 600);
    c50_100->cd();
    c50_100->SetMargin(0.12, 0.05, 0.1, 0.05);
//    hhdphi_50_100->Draw("E0 X0");
    hPhidphi_50_100->Draw("E0 X0 SAME");
    corrFit50100->Draw("SAME");
    corrFit2_50100->Draw("SAME");
    hPhidphi_50_100->Draw("E0 X0 SAME");
    hhBG_50_100->Draw("SAME");
    hphiBG_50_100->Draw("SAME");
    legend->Draw();

    TCanvas *cratio = new TCanvas("cratio", "cratio", 50, 50, 550, 600);
    cratio->cd();
    cratio->SetMargin(0.12, 0.05, 0.1, 0.05);
    ratio->Draw("E0 X0");
*/

}

