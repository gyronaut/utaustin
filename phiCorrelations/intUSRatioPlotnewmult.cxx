void intUSRatioPlotnewmult(){
    gStyle->SetOptStat(0);
    gStyle->SetOptFit(0);
    gStyle->SetErrorX(0);

    TFile *hhFile = new TFile("~/phiStudies/LHC16q_FAST_newmult/trig_4_8_assoc_2_4_hh_0_10.root");
    TH2D* hh2D_0_10 = (TH2D*)hhFile->Get("hh2D");
    hh2D_0_10->SetName("hh2D_0_10");
    hh2D_0_10->Sumw2();
    TH1D* hhdphi_0_10 = (TH1D*)hh2D_0_10->ProjectionY("hhdphi_0_10", hh2D_0_10->GetXaxis()->FindBin(-1.2), hh2D_0_10->GetXaxis()->FindBin(1.2));
    hhdphi_0_10->Rebin();
    hhdphi_0_10->SetLineWidth(4);
    hhdphi_0_10->SetLineColor(kBlue+2);
    hhdphi_0_10->SetMarkerColor(kBlue+2);
    hhdphi_0_10->SetMarkerStyle(21);
    hhdphi_0_10->SetMarkerSize(2);
    hhdphi_0_10->GetXaxis()->SetTitle("#Delta#varphi");
    hhdphi_0_10->SetTitle("");
    hhdphi_0_10->GetYaxis()->SetTitleOffset(1.60);
    hhdphi_0_10->GetXaxis()->SetTitleSize(0.05);
    hhdphi_0_10->GetXaxis()->SetTitleOffset(0.90);


    TFile* phiFile = new TFile("~/phiStudies/LHC16q_FAST_newmult/US_trig_4_8_assoc_2_4_mixcorr_hPhi_0_10.root");    
    TH2D* hPhi2D_0_10 = (TH2D*)phiFile->Get("AvgUSsubhPhi2Dpeak");
    hPhi2D_0_10->Sumw2();
    TH1D* hPhidphi_0_10 = (TH1D*)hPhi2D_0_10->ProjectionY("hPhidphi_0_10", hPhi2D_0_10->GetXaxis()->FindBin(-1.2), hPhi2D_0_10->GetXaxis()->FindBin(1.2));
    //hPhidphi_0_10->Rebin();
    hPhidphi_0_10->SetLineWidth(4);
    hPhidphi_0_10->SetLineColor(kRed+2);
    hPhidphi_0_10->SetMarkerColor(kRed+2);
    hPhidphi_0_10->SetMarkerStyle(22);
    hPhidphi_0_10->SetMarkerSize(2);
    hPhidphi_0_10->GetXaxis()->SetTitle("#Delta#varphi");
    hPhidphi_0_10->SetTitle("");
    hPhidphi_0_10->GetYaxis()->SetTitleOffset(1.20);
    hPhidphi_0_10->GetXaxis()->SetTitleSize(0.05);
    hPhidphi_0_10->GetXaxis()->SetTitleOffset(0.90);

    TF1 *corrFit = new TF1("corrFit", "gaus(0) + gaus(3) + pol0(6)", -1.4, 4.6);
    corrFit->FixParameter(6, 0.3333*(hhdphi_0_10->GetBinContent(8)+hhdphi_0_10->GetBinContent(16)+hhdphi_0_10->GetBinContent(1)));
    //corrFit->SetParLimits(6, hhdphi_0_10->GetBinContent(8)*0.9, hhdphi_0_10->GetBinContent(8)*1.1);
    corrFit->SetParameter(0, hhdphi_0_10->GetBinContent(hhdphi_0_10->GetXaxis()->FindBin(0.0)) - corrFit->GetParameter(6));
    corrFit->SetParLimits(0, corrFit->GetParameter(0)*0.5, corrFit->GetParameter(0)*1.5);
    corrFit->SetParameter(1, 0.0);
    corrFit->SetParLimits(1, -0.5, 0.5);
    corrFit->SetParameter(2, 0.5);
    corrFit->SetParLimits(2, 0.2, 0.9);
    corrFit->SetParameter(3, hhdphi_0_10->GetBinContent(hhdphi_0_10->GetXaxis()->FindBin(3.14)) - corrFit->GetParameter(6));
    corrFit->SetParLimits(3, corrFit->GetParameter(3)*0.5, corrFit->GetParameter(3)*1.5);
    corrFit->SetParameter(4, 3.14);
    corrFit->SetParLimits(4, 3.0, 3.25);
    corrFit->SetParameter(5, 0.5);
    corrFit->SetParLimits(5, 0.2, 0.9);

    corrFit->SetParNames("peak1", "mean1", "sigma1", "peak2", "mean2", "sigma2", "BG");

    corrFit->SetLineColor(kBlue);
    corrFit->SetLineWidth(6);
    corrFit->SetLineStyle(7);

    TF1 *corrFit2 = new TF1("corrFit2", "gaus(0) + gaus(3) + ([3]/([5]))*exp(-(((x - [4] + 2.0*TMath::Pi())^2)/(2*[5]^2))) + pol0(6)", -1.4, 4.6);
    //corrFit2->FixParameter(6, 0.3333*(hPhidphi_0_10->GetBinContent(8)+hPhidphi_0_10->GetBinContent(16)+hPhidphi_0_10->GetBinContent(1)));
    corrFit2->SetParameter(6, 0.3333*(hPhidphi_0_10->GetBinContent(8)+hPhidphi_0_10->GetBinContent(16)+hPhidphi_0_10->GetBinContent(1)));
    corrFit2->SetParLimits(6, corrFit2->GetParameter(6)*0.75, corrFit2->GetParameter(6)*1.25);
    corrFit2->SetParameter(0, hPhidphi_0_10->GetBinContent(hPhidphi_0_10->GetXaxis()->FindBin(0.0)) - corrFit2->GetParameter(6));
    corrFit2->SetParLimits(0, corrFit2->GetParameter(0)*0.5, corrFit2->GetParameter(0)*1.5);
    corrFit2->SetParameter(1, 0.0);
    corrFit2->SetParLimits(1, -0.5, 0.5);
    corrFit2->SetParameter(2, 0.5);
    corrFit2->SetParLimits(2, 0.2, 0.9);
    corrFit2->SetParameter(3, hPhidphi_0_10->GetBinContent(hPhidphi_0_10->GetXaxis()->FindBin(3.14)) - corrFit2->GetParameter(6));
    corrFit2->SetParLimits(3, corrFit2->GetParameter(3)*0.5, corrFit2->GetParameter(3)*1.5);
    corrFit2->SetParameter(4, 3.14);
    corrFit2->SetParLimits(4, 3.1, 3.2);
    corrFit2->SetParameter(5, 0.5);
    corrFit2->SetParLimits(5, 0.1, 1.0);

    corrFit2->SetParNames("peak1", "mean1", "sigma1", "peak2", "mean2", "sigma2", "BG");

    corrFit2->SetLineColor(kRed);
    corrFit2->SetLineWidth(6);
    corrFit2->SetLineStyle(7);

    hhdphi_0_10->Fit("corrFit", "R0");
    hPhidphi_0_10->Fit("corrFit2", "R0");

    TF1 *hhBG = new TF1("hhBG", "pol0(0)", -1.4, 4.6);
    hhBG->SetParLimits(0, 0.00001, 10000000.0);
    hhBG->SetParameter(0, 1.0*corrFit->GetParameter(6));
    hhBG->SetLineStyle(2);

    TF1 *hphiBG = new TF1("hphiBG", "pol0(0)", -1.4, 4.6);
    hphiBG->SetParLimits(0, 0.00001, 10000000.0);
    hphiBG->SetParameter(0, 1.0*corrFit2->GetParameter(6));
    hphiBG->SetLineStyle(2);

    Double_t near0_10hPhiError = 0;
    Double_t near0_10hhError = 0;
    Double_t away0_10hPhiError = 0;
    Double_t away0_10hhError = 0;
    Double_t mid0_10hPhiError =(1.0/3.0)*TMath::Sqrt(TMath::Power(hPhidphi_0_10->GetBinError(8),2) + TMath::Power(hPhidphi_0_10->GetBinError(16),2) + TMath::Power(hPhidphi_0_10->GetBinError(1),2));
    Double_t mid0_10hhError = (1.0/3.0)*TMath::Sqrt(TMath::Power(hhdphi_0_10->GetBinError(8),2) + TMath::Power(hhdphi_0_10->GetBinError(16),2) + TMath::Power(hhdphi_0_10->GetBinError(1),2));
    Double_t total0_10hPhiError;
    Double_t total0_10hhError = 0;   

    Double_t near0_10hPhiYield = hPhidphi_0_10->IntegralAndError(2,7,near0_10hPhiError) - hphiBG->GetParameter(0)*6.0;
    near0_10hPhiError = TMath::Sqrt(TMath::Power(near0_10hPhiError, 2) + TMath::Power(6.0*mid0_10hPhiError, 2));
    Double_t near0_10hhYield = hhdphi_0_10->IntegralAndError(2,7,near0_10hhError) - hhBG->GetParameter(0)*6.0;
    near0_10hhError = TMath::Sqrt(TMath::Power(near0_10hhError, 2) + TMath::Power(6.0*mid0_10hhError, 2));
    Double_t away0_10hPhiYield = hPhidphi_0_10->IntegralAndError(9,16,away0_10hPhiError) - hphiBG->GetParameter(0)*8.0;
    away0_10hPhiError = TMath::Sqrt(TMath::Power(away0_10hPhiError, 2) + TMath::Power(8.0*mid0_10hPhiError, 2));
    Double_t away0_10hhYield = hhdphi_0_10->IntegralAndError(9,16,away0_10hhError)- hhBG->GetParameter(0)*8.0;
    away0_10hhError = TMath::Sqrt(TMath::Power(away0_10hhError, 2) + TMath::Power(8.0*mid0_10hhError, 2));
    Double_t total0_10hPhiYield = hPhidphi_0_10->IntegralAndError(1,16,total0_10hPhiError);
    Double_t total0_10hhYield = hhdphi_0_10->IntegralAndError(1,16,total0_10hhError);
    Double_t mid0_10hPhiYield = hphiBG->GetParameter(0)*16.0;
    mid0_10hPhiError = mid0_10hPhiError*16.0;
    Double_t mid0_10hhYield = hhBG->GetParameter(0)*16.0;
    mid0_10hhError = mid0_10hhError*16.0; 
       
    Double_t near010 = near0_10hPhiYield/near0_10hhYield;
    Double_t near010Er = near010*TMath::Sqrt(TMath::Power(near0_10hPhiError/near0_10hPhiYield, 2) + TMath::Power(near0_10hhError/near0_10hhYield, 2));
    Double_t away010 = away0_10hPhiYield/away0_10hhYield;
    Double_t away010Er = away010*TMath::Sqrt(TMath::Power(away0_10hPhiError/away0_10hPhiYield, 2) + TMath::Power(away0_10hhError/away0_10hhYield, 2));
    Double_t mid010 = mid0_10hPhiYield/mid0_10hhYield;
    Double_t mid010Er = mid010*TMath::Sqrt(TMath::Power(mid0_10hPhiError/mid0_10hPhiYield, 2) + TMath::Power(mid0_10hhError/mid0_10hhYield, 2));
    Double_t total010 = total0_10hPhiYield/total0_10hhYield;
    Double_t total010Er = total010*TMath::Sqrt(TMath::Power(total0_10hPhiError/total0_10hPhiYield, 2) + TMath::Power(total0_10hhError/total0_10hhYield, 2));

    TH1D *ratios010 = new TH1D("ratios010", "(h-#phi / h-h) Ratios", 4, 0, 4);
    ratios010->GetXaxis()->SetBinLabel(1, "near-side");
    ratios010->SetBinContent(1, near010);
    ratios010->SetBinError(1, near010Er);
    ratios010->GetXaxis()->SetBinLabel(2, "mid");
    ratios010->SetBinContent(2, mid010);
    ratios010->SetBinError(2, mid010Er);
    ratios010->GetXaxis()->SetBinLabel(3, "away-side");
    ratios010->SetBinContent(3, away010);
    ratios010->SetBinError(3, away010Er);
    ratios010->GetXaxis()->SetBinLabel(4, "total");
    ratios010->SetBinContent(4, total010);
    ratios010->SetBinError(4, total010Er);
    ratios010->SetMarkerStyle(22);
    ratios010->SetMarkerColor(kRed+2);
    ratios010->SetLineColor(kRed+2);
    ratios010->SetMarkerSize(2);
    ratios010->SetLineWidth(2);

    TFile *hhFile_10_20 = new TFile("~/phiStudies/LHC16q_FAST_newmult/trig_4_8_assoc_2_4_hh_10_20.root");
    TH2D* hh2D_10_20 = (TH2D*)hhFile_10_20->Get("hh2D");
    hh2D_10_20->SetName("hh2D_10_20");
    hh2D_10_20->Sumw2();
    TH1D* hhdphi_10_20 = (TH1D*)hh2D_10_20->ProjectionY("hhdphi_10_20", hh2D_10_20->GetXaxis()->FindBin(-1.2), hh2D_10_20->GetXaxis()->FindBin(1.2));
    hhdphi_10_20->Rebin();
    hhdphi_10_20->SetLineWidth(4);
    hhdphi_10_20->SetLineColor(kBlue+2);
    hhdphi_10_20->SetMarkerColor(kBlue+2);
    hhdphi_10_20->SetMarkerStyle(21);
    hhdphi_10_20->SetMarkerSize(2);
    hhdphi_10_20->GetXaxis()->SetTitle("#Delta#varphi");
    hhdphi_10_20->SetTitle("");
    hhdphi_10_20->GetYaxis()->SetTitleOffset(1.60);
    hhdphi_10_20->GetXaxis()->SetTitleSize(0.05);
    hhdphi_10_20->GetXaxis()->SetTitleOffset(0.90);


    //10-20 section
    TFile* phiFile_10_20 = new TFile("~/phiStudies/LHC16q_FAST_newmult/US_trig_4_8_assoc_2_4_mixcorr_hPhi_10_20.root");
    TH2D* hPhi2D_10_20 = (TH2D*)phiFile_10_20->Get("AvgUSsubhPhi2Dpeak");
    hPhi2D_10_20->Sumw2();
    TH1D* hPhidphi_10_20 = (TH1D*)hPhi2D_10_20->ProjectionY("hPhidphi_10_20", hPhi2D_10_20->GetXaxis()->FindBin(-1.2), hPhi2D_10_20->GetXaxis()->FindBin(1.2));
    //hPhidphi_10_20->Rebin();
    hPhidphi_10_20->SetLineWidth(4);
    hPhidphi_10_20->SetLineColor(kRed+2);
    hPhidphi_10_20->SetMarkerColor(kRed+2);
    hPhidphi_10_20->SetMarkerStyle(22);
    hPhidphi_10_20->SetMarkerSize(2);
    hPhidphi_10_20->GetXaxis()->SetTitle("#Delta#varphi");
    hPhidphi_10_20->SetTitle("");
    hPhidphi_10_20->GetYaxis()->SetTitleOffset(1.20);
    hPhidphi_10_20->GetXaxis()->SetTitleSize(0.05);
    hPhidphi_10_20->GetXaxis()->SetTitleOffset(0.90);

    TF1 *corrFit_10_20 = new TF1("corrFit_10_20", "gaus(0) + gaus(3) + pol0(6)", -1.4, 4.6);
    corrFit_10_20->FixParameter(6, 0.3333*(hhdphi_10_20->GetBinContent(8)+hhdphi_10_20->GetBinContent(16)+hhdphi_10_20->GetBinContent(1)));
    //corrFit_10_20->SetParLimits(6, hhdphi_10_20->GetBinContent(8)*0.9, hhdphi_10_20->GetBinContent(8)*1.1);
    corrFit_10_20->SetParameter(0, hhdphi_10_20->GetBinContent(hhdphi_10_20->GetXaxis()->FindBin(0.0)) - corrFit_10_20->GetParameter(6));
    corrFit_10_20->SetParLimits(0, corrFit_10_20->GetParameter(0)*0.5, corrFit_10_20->GetParameter(0)*1.5);
    corrFit_10_20->SetParameter(1, 0.0);
    corrFit_10_20->SetParLimits(1, -0.5, 0.5);
    corrFit_10_20->SetParameter(2, 0.5);
    corrFit_10_20->SetParLimits(2, 0.2, 0.9);
    corrFit_10_20->SetParameter(3, hhdphi_10_20->GetBinContent(hhdphi_10_20->GetXaxis()->FindBin(3.14)) - corrFit_10_20->GetParameter(6));
    corrFit_10_20->SetParLimits(3, corrFit_10_20->GetParameter(3)*0.5, corrFit_10_20->GetParameter(3)*1.5);
    corrFit_10_20->SetParameter(4, 3.14);
    corrFit_10_20->SetParLimits(4, 3.0, 3.25);
    corrFit_10_20->SetParameter(5, 0.5);
    corrFit_10_20->SetParLimits(5, 0.2, 0.9);

    corrFit_10_20->SetParNames("peak1", "mean1", "sigma1", "peak2", "mean2", "sigma2", "BG");

    corrFit_10_20->SetLineColor(kBlue);
    corrFit_10_20->SetLineWidth(6);
    corrFit_10_20->SetLineStyle(7);

    TF1 *corrFit2_10_20 = new TF1("corrFit2_10_20", "gaus(0) + gaus(3) + ([3]/([5]))*exp(-(((x - [4] + 2.0*TMath::Pi())^2)/(2*[5]^2))) + pol0(6)", -1.4, 4.6);
    //corrFit2_10_20->FixParameter(6, 0.3333*(hPhidphi_10_20->GetBinContent(8)+hPhidphi_10_20->GetBinContent(16)+hPhidphi_10_20->GetBinContent(1)));
    corrFit2_10_20->SetParameter(6, 0.3333*(hPhidphi_10_20->GetBinContent(8)+hPhidphi_10_20->GetBinContent(16)+hPhidphi_10_20->GetBinContent(1)));
    corrFit2_10_20->SetParLimits(6, corrFit2_10_20->GetParameter(6)*0.75, corrFit2_10_20->GetParameter(6)*1.25);
    corrFit2_10_20->SetParameter(0, hPhidphi_10_20->GetBinContent(hPhidphi_10_20->GetXaxis()->FindBin(0.0)) - corrFit2_10_20->GetParameter(6));
    corrFit2_10_20->SetParLimits(0, corrFit2_10_20->GetParameter(0)*0.5, corrFit2_10_20->GetParameter(0)*1.5);
    corrFit2_10_20->SetParameter(1, 0.0);
    corrFit2_10_20->SetParLimits(1, -0.5, 0.5);
    corrFit2_10_20->SetParameter(2, 0.5);
    corrFit2_10_20->SetParLimits(2, 0.2, 0.9);
    corrFit2_10_20->SetParameter(3, hPhidphi_10_20->GetBinContent(hPhidphi_10_20->GetXaxis()->FindBin(3.14)) - corrFit2_10_20->GetParameter(6));
    corrFit2_10_20->SetParLimits(3, corrFit2_10_20->GetParameter(3)*0.5, corrFit2_10_20->GetParameter(3)*1.5);
    corrFit2_10_20->SetParameter(4, 3.14);
    corrFit2_10_20->SetParLimits(4, 3.1, 3.2);
    corrFit2_10_20->SetParameter(5, 0.5);
    corrFit2_10_20->SetParLimits(5, 0.1, 1.0);

    corrFit2_10_20->SetParNames("peak1", "mean1", "sigma1", "peak2", "mean2", "sigma2", "BG");

    corrFit2_10_20->SetLineColor(kRed);
    corrFit2_10_20->SetLineWidth(6);
    corrFit2_10_20->SetLineStyle(7);

    hhdphi_10_20->Fit("corrFit_10_20", "R0");
    hPhidphi_10_20->Fit("corrFit2_10_20", "R0");

    TF1 *hhBG_10_20 = new TF1("hhBG_10_20", "pol0(0)", -1.4, 4.6);
    hhBG_10_20->SetParLimits(0, 0.00001, 10000000.0);
    hhBG_10_20->SetParameter(0, 1.0*corrFit_10_20->GetParameter(6));
    hhBG_10_20->SetLineStyle(2);

    TF1 *hphiBG_10_20 = new TF1("hphiBG_10_20", "pol0(0)", -1.4, 4.6);
    hphiBG_10_20->SetParLimits(0, 0.00001, 10000000.0);
    hphiBG_10_20->SetParameter(0, 1.0*corrFit2_10_20->GetParameter(6));
    hphiBG_10_20->SetLineStyle(2);

    Double_t near10_20hPhiError = 0;
    Double_t near10_20hhError = 0;
    Double_t away10_20hPhiError = 0;
    Double_t away10_20hhError = 0;
    Double_t mid10_20hPhiError =(1.0/3.0)*TMath::Sqrt(TMath::Power(hPhidphi_10_20->GetBinError(8),2) + TMath::Power(hPhidphi_10_20->GetBinError(16),2) + TMath::Power(hPhidphi_10_20->GetBinError(1),2));
    Double_t mid10_20hhError = (1.0/3.0)*TMath::Sqrt(TMath::Power(hhdphi_10_20->GetBinError(8),2) + TMath::Power(hhdphi_10_20->GetBinError(16),2) + TMath::Power(hhdphi_10_20->GetBinError(1),2));
    Double_t total10_20hPhiError;
    Double_t total10_20hhError = 0;   

    Double_t near10_20hPhiYield = hPhidphi_10_20->IntegralAndError(2,7,near10_20hPhiError) - hphiBG_10_20->GetParameter(0)*6.0;
    near10_20hPhiError = TMath::Sqrt(TMath::Power(near10_20hPhiError, 2) + TMath::Power(6.0*mid10_20hPhiError, 2));
    Double_t near10_20hhYield = hhdphi_10_20->IntegralAndError(2,7,near10_20hhError) - hhBG_10_20->GetParameter(0)*6.0;
    near10_20hhError = TMath::Sqrt(TMath::Power(near10_20hhError, 2) + TMath::Power(6.0*mid10_20hhError, 2));
    Double_t away10_20hPhiYield = hPhidphi_10_20->IntegralAndError(9,16,away10_20hPhiError) - hphiBG_10_20->GetParameter(0)*8.0;
    away10_20hPhiError = TMath::Sqrt(TMath::Power(away10_20hPhiError, 2) + TMath::Power(8.0*mid10_20hPhiError, 2));
    Double_t away10_20hhYield = hhdphi_10_20->IntegralAndError(9,16,away10_20hhError)- hhBG_10_20->GetParameter(0)*8.0;
    away10_20hhError = TMath::Sqrt(TMath::Power(away10_20hhError, 2) + TMath::Power(8.0*mid10_20hhError, 2));
    Double_t total10_20hPhiYield = hPhidphi_10_20->IntegralAndError(1,16,total10_20hPhiError);
    Double_t total10_20hhYield = hhdphi_10_20->IntegralAndError(1,16,total10_20hhError);
    Double_t mid10_20hPhiYield = hphiBG_10_20->GetParameter(0)*16.0;
    mid10_20hPhiError = mid10_20hPhiError*16.0;
    Double_t mid10_20hhYield = hhBG_10_20->GetParameter(0)*16.0;
    mid10_20hhError = mid10_20hhError*16.0; 
       
    Double_t near1020 = near10_20hPhiYield/near10_20hhYield;
    Double_t near1020Er = near1020*TMath::Sqrt(TMath::Power(near10_20hPhiError/near10_20hPhiYield, 2) + TMath::Power(near10_20hhError/near10_20hhYield, 2));
    Double_t away1020 = away10_20hPhiYield/away10_20hhYield;
    Double_t away1020Er = away1020*TMath::Sqrt(TMath::Power(away10_20hPhiError/away10_20hPhiYield, 2) + TMath::Power(away10_20hhError/away10_20hhYield, 2));
    Double_t mid1020 = mid10_20hPhiYield/mid10_20hhYield;
    Double_t mid1020Er = mid1020*TMath::Sqrt(TMath::Power(mid10_20hPhiError/mid10_20hPhiYield, 2) + TMath::Power(mid10_20hhError/mid10_20hhYield, 2));
    Double_t total1020 = total10_20hPhiYield/total10_20hhYield;
    Double_t total1020Er = total1020*TMath::Sqrt(TMath::Power(total10_20hPhiError/total10_20hPhiYield, 2) + TMath::Power(total10_20hhError/total10_20hhYield, 2));

    TH1D *ratios1020 = new TH1D("ratios1020", "(h-#phi / h-h) Ratios", 4, 0, 4);
    ratios1020->GetXaxis()->SetBinLabel(1, "near-side");
    ratios1020->SetBinContent(1, near1020);
    ratios1020->SetBinError(1, near1020Er);
    ratios1020->GetXaxis()->SetBinLabel(2, "mid");
    ratios1020->SetBinContent(2, mid1020);
    ratios1020->SetBinError(2, mid1020Er);
    ratios1020->GetXaxis()->SetBinLabel(3, "away-side");
    ratios1020->SetBinContent(3, away1020);
    ratios1020->SetBinError(3, away1020Er);
    ratios1020->GetXaxis()->SetBinLabel(4, "total");
    ratios1020->SetBinContent(4, total1020);
    ratios1020->SetBinError(4, total1020Er);
    ratios1020->SetMarkerStyle(22);
    ratios1020->SetMarkerColor(kRed+2);
    ratios1020->SetLineColor(kRed+2);
    ratios1020->SetMarkerSize(2);
    ratios1020->SetLineWidth(2);


   //20-40 section 
    TFile *hhFile_20_40 = new TFile("~/phiStudies/LHC16q_FAST_newmult/trig_4_8_assoc_2_4_hh_20_40.root");
    TH2D* hh2D_20_40 = (TH2D*)hhFile_20_40->Get("hh2D");
    hh2D_20_40->SetName("hh2D_20_40");
    TH1D* hhdphi_20_40 = (TH1D*)hh2D_20_40->ProjectionY("hhdphi_20_40", hh2D_20_40->GetXaxis()->FindBin(-1.2), hh2D_20_40->GetXaxis()->FindBin(1.2));
    hhdphi_20_40->Rebin();
    hhdphi_20_40->SetLineWidth(4);
    hhdphi_20_40->SetLineColor(kBlue+2);
    hhdphi_20_40->SetMarkerColor(kBlue+2);
    hhdphi_20_40->SetMarkerStyle(21);
    hhdphi_20_40->SetMarkerSize(2);
    hhdphi_20_40->GetXaxis()->SetTitle("#Delta#varphi");
    hhdphi_20_40->SetTitle("");
    hhdphi_20_40->GetYaxis()->SetTitleOffset(1.60);
    hhdphi_20_40->GetXaxis()->SetTitleSize(0.05);
    hhdphi_20_40->GetXaxis()->SetTitleOffset(0.90);


    TFile* phiFile_20_40 = new TFile("~/phiStudies/LHC16q_FAST_newmult/US_trig_4_8_assoc_2_4_mixcorr_hPhi_20_40.root");    
    TH2D* hPhi2D_20_40 = (TH2D*)phiFile_20_40->Get("AvgUSsubhPhi2Dpeak");
    hPhi2D_20_40->SetName("hPhi2D_20_40");
    TH1D* hPhidphi_20_40 = (TH1D*)hPhi2D_20_40->ProjectionY("hPhidphi_20_40", hPhi2D_20_40->GetXaxis()->FindBin(-1.2), hPhi2D_20_40->GetXaxis()->FindBin(1.2));
    //hPhidphi_20_40->Rebin();
    hPhidphi_20_40->SetLineWidth(4);
    hPhidphi_20_40->SetLineColor(kRed+2);
    hPhidphi_20_40->SetMarkerColor(kRed+2);
    hPhidphi_20_40->SetMarkerStyle(22);
    hPhidphi_20_40->SetMarkerSize(2);
    hPhidphi_20_40->GetXaxis()->SetTitle("#Delta#varphi");
    hPhidphi_20_40->SetTitle("");
    hPhidphi_20_40->GetYaxis()->SetTitleOffset(1.20);
    hPhidphi_20_40->GetXaxis()->SetTitleSize(0.05);
    hPhidphi_20_40->GetXaxis()->SetTitleOffset(0.90);

    TF1 *corrFit2040 = new TF1("corrFit2040", "gaus(0) + gaus(3) + ([3]/([5]))*exp(-(((x - [4] + 2.0*TMath::Pi())^2)/(2*[5]^2))) + pol0(6)", -1.4, 4.6);
    corrFit2040->FixParameter(6, 0.3333*(hhdphi_20_40->GetBinContent(8)+hhdphi_20_40->GetBinContent(16)+hhdphi_20_40->GetBinContent(1)));
    //corrFit2040->SetParLimits(6, hhdphi_20_40->GetBinContent(8)*0.9, hhdphi_20_40->GetBinContent(8)*1.1);
    corrFit2040->SetParameter(0, hhdphi_20_40->GetBinContent(hhdphi_20_40->GetXaxis()->FindBin(0.0)) - corrFit2040->GetParameter(6));
    corrFit2040->SetParLimits(0, corrFit2040->GetParameter(0)*0.5, corrFit2040->GetParameter(0)*1.5);
    corrFit2040->SetParameter(1, 0.0);
    corrFit2040->SetParLimits(1, -0.5, 0.5);
    corrFit2040->SetParameter(2, 0.5);
    corrFit2040->SetParLimits(2, 0.2, 0.9);
    corrFit2040->SetParameter(3, hhdphi_20_40->GetBinContent(hhdphi_20_40->GetXaxis()->FindBin(3.14)) - corrFit2040->GetParameter(6));
    corrFit2040->SetParLimits(3, corrFit2040->GetParameter(3)*0.5, corrFit2040->GetParameter(3)*1.5);
    corrFit2040->SetParameter(4, 3.14);
    corrFit2040->SetParLimits(4, 3.0, 3.25);
    corrFit2040->SetParameter(5, 0.5);
    corrFit2040->SetParLimits(5, 0.2, 0.9);

    corrFit2040->SetParNames("peak1", "mean1", "sigma1", "peak2", "mean2", "sigma2", "BG");

    corrFit2040->SetLineColor(kBlue);
    corrFit2040->SetLineWidth(6);
    corrFit2040->SetLineStyle(7);

    TF1 *corrFit2_2040 = new TF1("corrFit2_2040", "gaus(0) + gaus(3) + ([3]/([5]))*exp(-(((x - [4] + 2.0*TMath::Pi())^2)/(2*[5]^2))) + pol0(6)", -1.4, 4.6);
    //corrFit2_2040->FixParameter(6, 0.3333*(hPhidphi_20_40->GetBinContent(8)+hPhidphi_20_40->GetBinContent(16)+hPhidphi_20_40->GetBinContent(1)));
    corrFit2_2040->SetParameter(6, 0.3333*(hPhidphi_20_40->GetBinContent(8)+hPhidphi_20_40->GetBinContent(16)+hPhidphi_20_40->GetBinContent(1)));
    corrFit2_2040->SetParLimits(6, corrFit2_2040->GetParameter(6)*0.75, corrFit2_2040->GetParameter(6)*1.25);
    corrFit2_2040->SetParameter(0, hPhidphi_20_40->GetBinContent(hPhidphi_20_40->GetXaxis()->FindBin(0.0)) - corrFit2_2040->GetParameter(6));
    corrFit2_2040->SetParLimits(0, corrFit2_2040->GetParameter(0)*0.5, corrFit2_2040->GetParameter(0)*1.5);
    corrFit2_2040->SetParameter(1, 0.0);
    corrFit2_2040->SetParLimits(1, -0.5, 0.5);
    corrFit2_2040->SetParameter(2, 0.5);
    corrFit2_2040->SetParLimits(2, 0.2, 0.9);
    corrFit2_2040->SetParameter(3, hPhidphi_20_40->GetBinContent(hPhidphi_20_40->GetXaxis()->FindBin(3.14)) - corrFit2_2040->GetParameter(6));
    corrFit2_2040->SetParLimits(3, corrFit2_2040->GetParameter(3)*0.5, corrFit2_2040->GetParameter(3)*1.5);
    corrFit2_2040->SetParameter(4, 3.14);
    corrFit2_2040->SetParLimits(4, 3.1, 3.2);
    corrFit2_2040->SetParameter(5, 0.5);
    corrFit2_2040->SetParLimits(5, 0.1, 1.0);

    corrFit2_2040->SetParNames("peak1", "mean1", "sigma1", "peak2", "mean2", "sigma2", "BG");

    corrFit2_2040->SetLineColor(kRed);
    corrFit2_2040->SetLineWidth(6);
    corrFit2_2040->SetLineStyle(7);

    hhdphi_20_40->Fit("corrFit2040", "R0");
    hPhidphi_20_40->Fit("corrFit2_2040", "R0");

    TF1 *hhBG_20_40 = new TF1("hhBG_20_40", "pol0(0)", -1.4, 4.6);
    hhBG_20_40->SetParLimits(0, 0.00001, 10000000.0);
    hhBG_20_40->SetParameter(0, 1.0*corrFit2040->GetParameter(6));
    hhBG_20_40->SetLineStyle(2);

    TF1 *hphiBG_20_40 = new TF1("hphiBG_20_40", "pol0(0)", -1.4, 4.6);
    hphiBG_20_40->SetParLimits(0, 0.00001, 10000000.0);
    hphiBG_20_40->SetParameter(0, 1.0*corrFit2_2040->GetParameter(6));
    hphiBG_20_40->SetLineStyle(2);

    Double_t near20_40hPhiError = 0;
    Double_t near20_40hhError = 0;
    Double_t away20_40hPhiError = 0;
    Double_t away20_40hhError = 0;
    Double_t mid20_40hPhiError =(1.0/3.0)*TMath::Sqrt(TMath::Power(hPhidphi_20_40->GetBinError(8),2) + TMath::Power(hPhidphi_20_40->GetBinError(16),2) + TMath::Power(hPhidphi_20_40->GetBinError(1),2));
    Double_t mid20_40hhError = (1.0/3.0)*TMath::Sqrt(TMath::Power(hhdphi_20_40->GetBinError(8),2) + TMath::Power(hhdphi_20_40->GetBinError(16),2) + TMath::Power(hhdphi_20_40->GetBinError(1),2));
    Double_t total20_40hPhiError = 0;
    Double_t total20_40hhError = 0;

    Double_t near20_40hPhiYield = hPhidphi_20_40->IntegralAndError(2,7,near20_40hPhiError) - hphiBG_20_40->GetParameter(0)*6.0;
    near20_40hPhiError = TMath::Sqrt(TMath::Power(near20_40hPhiError, 2) + TMath::Power(6.0*mid20_40hPhiError, 2));
    Double_t near20_40hhYield = hhdphi_20_40->IntegralAndError(2,7,near20_40hhError) - hhBG_20_40->GetParameter(0)*6.0;
    near20_40hhError = TMath::Sqrt(TMath::Power(near20_40hhError, 2) + TMath::Power(6.0*mid20_40hhError, 2));
    Double_t away20_40hPhiYield = hPhidphi_20_40->IntegralAndError(9,16,away20_40hPhiError) - hphiBG_20_40->GetParameter(0)*8.0;
    away20_40hPhiError = TMath::Sqrt(TMath::Power(away20_40hPhiError, 2) + TMath::Power(8.0*mid20_40hPhiError, 2));
    Double_t away20_40hhYield = hhdphi_20_40->IntegralAndError(9,16,away20_40hhError)- hhBG_20_40->GetParameter(0)*8.0;
    away20_40hhError = TMath::Sqrt(TMath::Power(away20_40hhError, 2) + TMath::Power(8.0*mid20_40hhError, 2));
    Double_t mid20_40hPhiYield = hphiBG_20_40->GetParameter(0)*16.0;
    mid20_40hPhiError = mid20_40hPhiError*16.0;
    Double_t mid20_40hhYield = hhBG_20_40->GetParameter(0)*16.0;
    mid20_40hhError = mid20_40hhError*16.0;
    Double_t total20_40hPhiYield = hPhidphi_20_40->IntegralAndError(1, 16,total20_40hPhiError);
    Double_t total20_40hhYield = hhdphi_20_40->IntegralAndError(1, 16,total20_40hhError);

    Double_t near2040 = near20_40hPhiYield/near20_40hhYield;
    Double_t near2040Er = near2040*TMath::Sqrt(TMath::Power(near20_40hPhiError/near20_40hPhiYield, 2) + TMath::Power(near20_40hhError/near20_40hhYield, 2));
    Double_t away2040 = away20_40hPhiYield/away20_40hhYield;
    Double_t away2040Er = away2040*TMath::Sqrt(TMath::Power(away20_40hPhiError/away20_40hPhiYield, 2) + TMath::Power(away20_40hhError/away20_40hhYield, 2));
    Double_t mid2040 = mid20_40hPhiYield/mid20_40hhYield;
    Double_t mid2040Er = mid2040*TMath::Sqrt(TMath::Power(mid20_40hPhiError/mid20_40hPhiYield, 2) + TMath::Power(mid20_40hhError/mid20_40hhYield, 2));
    Double_t total2040 = total20_40hPhiYield/total20_40hhYield;
    Double_t total2040Er = total2040*TMath::Sqrt(TMath::Power(total20_40hPhiError/total20_40hPhiYield, 2) + TMath::Power(total20_40hhError/total20_40hhYield, 2));

    TH1D *ratios2040 = new TH1D("ratios2040", "(h-#phi / h-h) Ratios", 4, 0, 4);
    ratios2040->GetXaxis()->SetBinLabel(1, "near-side");
    ratios2040->SetBinContent(1, near2040);
    ratios2040->SetBinError(1, near2040Er);
    ratios2040->GetXaxis()->SetBinLabel(2, "mid");
    ratios2040->SetBinContent(2, mid2040);
    ratios2040->SetBinError(2, mid2040Er);
    ratios2040->GetXaxis()->SetBinLabel(3, "away-side");
    ratios2040->SetBinContent(3, away2040);
    ratios2040->SetBinError(3, away2040Er);
    ratios2040->GetXaxis()->SetBinLabel(4, "total");
    ratios2040->SetBinContent(4, total2040);
    ratios2040->SetBinError(4, total2040Er);
    ratios2040->SetMarkerStyle(21);
    ratios2040->SetMarkerColor(kBlue+2);
    ratios2040->SetLineColor(kBlue+2);
    ratios2040->SetMarkerSize(2);
    ratios2040->SetLineWidth(2);


    //40-90 section
    TFile *hhFile_40_90 = new TFile("~/phiStudies/LHC16q_FAST_newmult/trig_4_8_assoc_2_4_hh_40_90.root");
    TH2D* hh2D_40_90 = (TH2D*)hhFile_40_90->Get("hh2D");
    hh2D_40_90->SetName("hh2D_40_90");
    TH1D* hhdphi_40_90 = (TH1D*)hh2D_40_90->ProjectionY("hhdphi_40_90", hh2D_40_90->GetXaxis()->FindBin(-1.2), hh2D_40_90->GetXaxis()->FindBin(1.2));
    hhdphi_40_90->Rebin();
    hhdphi_40_90->SetLineWidth(4);
    hhdphi_40_90->SetLineColor(kBlue+2);
    hhdphi_40_90->SetMarkerColor(kBlue+2);
    hhdphi_40_90->SetMarkerStyle(21);
    hhdphi_40_90->SetMarkerSize(2);
    hhdphi_40_90->GetXaxis()->SetTitle("#Delta#varphi");
    hhdphi_40_90->SetTitle("");
    hhdphi_40_90->GetYaxis()->SetTitleOffset(1.60);
    hhdphi_40_90->GetXaxis()->SetTitleSize(0.05);
    hhdphi_40_90->GetXaxis()->SetTitleOffset(0.90);


    TFile* phiFile_40_90 = new TFile("~/phiStudies/LHC16q_FAST_newmult/US_trig_4_8_assoc_2_4_mixcorr_hPhi_40_90.root");    
    TH2D* hPhi2D_40_90 = (TH2D*)phiFile_40_90->Get("AvgUSsubhPhi2Dpeak");
    hPhi2D_40_90->SetName("hPhi2D_40_90");
    TH1D* hPhidphi_40_90 = (TH1D*)hPhi2D_40_90->ProjectionY("hPhidphi_40_90", hPhi2D_40_90->GetXaxis()->FindBin(-1.2), hPhi2D_40_90->GetXaxis()->FindBin(1.2));
    //hPhidphi_40_90->Rebin();
    hPhidphi_40_90->SetLineWidth(4);
    hPhidphi_40_90->SetLineColor(kRed+2);
    hPhidphi_40_90->SetMarkerColor(kRed+2);
    hPhidphi_40_90->SetMarkerStyle(22);
    hPhidphi_40_90->SetMarkerSize(2);
    hPhidphi_40_90->GetXaxis()->SetTitle("#Delta#varphi");
    hPhidphi_40_90->SetTitle("");
    hPhidphi_40_90->GetYaxis()->SetTitleOffset(1.20);
    hPhidphi_40_90->GetXaxis()->SetTitleSize(0.05);
    hPhidphi_40_90->GetXaxis()->SetTitleOffset(0.90);

    TF1 *corrFit4090 = new TF1("corrFit4090", "gaus(0) + gaus(3) + pol0(6)", -1.4, 4.6);
    corrFit4090->FixParameter(6, 0.3333*(hhdphi_40_90->GetBinContent(8)+hhdphi_40_90->GetBinContent(16)+hhdphi_40_90->GetBinContent(1)));
    //corrFit4090->SetParLimits(6, hhdphi_40_90->GetBinContent(8)*0.9, hhdphi_40_90->GetBinContent(8)*1.1);
    corrFit4090->SetParameter(0, hhdphi_40_90->GetBinContent(hhdphi_40_90->GetXaxis()->FindBin(0.0)) - corrFit4090->GetParameter(6));
    corrFit4090->SetParLimits(0, corrFit4090->GetParameter(0)*0.5, corrFit4090->GetParameter(0)*1.5);
    corrFit4090->SetParameter(1, 0.0);
    corrFit4090->SetParLimits(1, -0.5, 0.5);
    corrFit4090->SetParameter(2, 0.5);
    corrFit4090->SetParLimits(2, 0.2, 0.9);
    corrFit4090->SetParameter(3, hhdphi_40_90->GetBinContent(hhdphi_40_90->GetXaxis()->FindBin(3.14)) - corrFit4090->GetParameter(6));
    corrFit4090->SetParLimits(3, corrFit4090->GetParameter(3)*0.5, corrFit4090->GetParameter(3)*1.5);
    corrFit4090->SetParameter(4, 3.14);
    corrFit4090->SetParLimits(4, 3.0, 3.25);
    corrFit4090->SetParameter(5, 0.5);
    corrFit4090->SetParLimits(5, 0.1, 1.0);

    corrFit4090->SetParNames("peak1", "mean1", "sigma1", "peak2", "mean2", "sigma2", "BG");

    corrFit4090->SetLineColor(kBlue);
    corrFit4090->SetLineWidth(6);
    corrFit4090->SetLineStyle(7);

    TF1 *corrFit2_4090 = new TF1("corrFit2_4090", "gaus(0) + gaus(3) + ([3]/([5]))*exp(-(((x - [4] + 2.0*TMath::Pi())^2)/(2*[5]^2))) + pol0(6)", -1.4, 4.6);
    //corrFit2_4090->FixParameter(6, 0.3333*(hPhidphi_40_90->GetBinContent(8)+hPhidphi_40_90->GetBinContent(16)+hPhidphi_40_90->GetBinContent(1)));
    corrFit2_4090->SetParameter(6, 0.3333*(hPhidphi_40_90->GetBinContent(8)+hPhidphi_40_90->GetBinContent(16)+hPhidphi_40_90->GetBinContent(1)));
    corrFit2_4090->SetParLimits(6, corrFit2_4090->GetParameter(6)*0.75, corrFit2_4090->GetParameter(6)*1.25);
    corrFit2_4090->SetParameter(0, hPhidphi_40_90->GetBinContent(hPhidphi_40_90->GetXaxis()->FindBin(0.0)) - corrFit2_4090->GetParameter(6));
    corrFit2_4090->SetParLimits(0, corrFit2_4090->GetParameter(0)*0.5, corrFit2_4090->GetParameter(0)*1.5);
    corrFit2_4090->SetParameter(1, 0.0);
    corrFit2_4090->SetParLimits(1, -0.5, 0.5);
    corrFit2_4090->SetParameter(2, 0.5);
    corrFit2_4090->SetParLimits(2, 0.2, 0.9);
    corrFit2_4090->SetParameter(3, hPhidphi_40_90->GetBinContent(hPhidphi_40_90->GetXaxis()->FindBin(3.14)) - corrFit2_4090->GetParameter(6));
    corrFit2_4090->SetParLimits(3, corrFit2_4090->GetParameter(3)*0.5, corrFit2_4090->GetParameter(3)*1.5);
    corrFit2_4090->SetParameter(4, 3.14);
    corrFit2_4090->SetParLimits(4, 3.1, 3.2);
    corrFit2_4090->SetParameter(5, 0.5);
    corrFit2_4090->SetParLimits(5, 0.1, 0.9);

    corrFit2_4090->SetParNames("peak1", "mean1", "sigma1", "peak2", "mean2", "sigma2", "BG");

    corrFit2_4090->SetLineColor(kRed);
    corrFit2_4090->SetLineWidth(6);
    corrFit2_4090->SetLineStyle(7);

    hhdphi_40_90->Fit("corrFit4090", "R0");
    hPhidphi_40_90->Fit("corrFit2_4090", "R0");

    TF1 *hhBG_40_90 = new TF1("hhBG_40_90", "pol0(0)", -1.4, 4.6);
    hhBG_40_90->SetParLimits(0, 0.00001, 10000000.0);
    hhBG_40_90->SetParameter(0, 1.0*corrFit4090->GetParameter(6));
    hhBG_40_90->SetLineStyle(2);

    TF1 *hphiBG_40_90 = new TF1("hphiBG_40_90", "pol0(0)", -1.4, 4.6);
    hphiBG_40_90->SetParLimits(0, 0.00001, 10000000.0);
    hphiBG_40_90->SetParameter(0, 1.0*corrFit2_4090->GetParameter(6));
    hphiBG_40_90->SetLineStyle(2);

    Double_t near40_90hPhiError = 0;
    Double_t near40_90hhError = 0;
    Double_t away40_90hPhiError = 0;
    Double_t away40_90hhError = 0;
    Double_t mid40_90hPhiError =(1.0/3.0)*TMath::Sqrt(TMath::Power(hPhidphi_40_90->GetBinError(8),2) + TMath::Power(hPhidphi_40_90->GetBinError(16),2) + TMath::Power(hPhidphi_40_90->GetBinError(1),2));
    Double_t mid40_90hhError = (1.0/3.0)*TMath::Sqrt(TMath::Power(hhdphi_40_90->GetBinError(8),2) + TMath::Power(hhdphi_40_90->GetBinError(16),2) + TMath::Power(hhdphi_40_90->GetBinError(1),2));
    Double_t total40_90hPhiError = 0;
    Double_t total40_90hhError = 0;

    Double_t near40_90hPhiYield = hPhidphi_40_90->IntegralAndError(2,7,near40_90hPhiError) - hphiBG_40_90->GetParameter(0)*6.0;
    near40_90hPhiError = TMath::Sqrt(TMath::Power(near40_90hPhiError, 2) + TMath::Power(6.0*mid40_90hPhiError, 2));
    Double_t near40_90hhYield = hhdphi_40_90->IntegralAndError(2,7,near40_90hhError) - hhBG_40_90->GetParameter(0)*6.0;
    near40_90hhError = TMath::Sqrt(TMath::Power(near40_90hhError, 2) + TMath::Power(6.0*mid40_90hhError, 2));
    Double_t away40_90hPhiYield = hPhidphi_40_90->IntegralAndError(9,16,away40_90hPhiError) - hphiBG_40_90->GetParameter(0)*8.0;
    away40_90hPhiError = TMath::Sqrt(TMath::Power(away40_90hPhiError, 2) + TMath::Power(8.0*mid40_90hPhiError, 2));
    Double_t away40_90hhYield = hhdphi_40_90->IntegralAndError(9,16,away40_90hhError)- hhBG_40_90->GetParameter(0)*8.0;
    away40_90hhError = TMath::Sqrt(TMath::Power(away40_90hhError, 2) + TMath::Power(8.0*mid40_90hhError, 2));
    Double_t mid40_90hPhiYield = hphiBG_40_90->GetParameter(0)*16.0;
    mid40_90hPhiError = mid40_90hPhiError*16.0;
    Double_t mid40_90hhYield = hhBG_40_90->GetParameter(0)*16.0;
    mid40_90hhError = mid40_90hhError*16.0;
    Double_t total40_90hPhiYield = hPhidphi_40_90->IntegralAndError(1, 16,total40_90hPhiError);
    Double_t total40_90hhYield = hhdphi_40_90->IntegralAndError(1, 16,total40_90hhError);

    Double_t near4090 = near40_90hPhiYield/near40_90hhYield;
    Double_t near4090Er = near4090*TMath::Sqrt(TMath::Power(near40_90hPhiError/near40_90hPhiYield, 2) + TMath::Power(near40_90hhError/near40_90hhYield, 2));
    Double_t away4090 = away40_90hPhiYield/away40_90hhYield;
    Double_t away4090Er = away4090*TMath::Sqrt(TMath::Power(away40_90hPhiError/away40_90hPhiYield, 2) + TMath::Power(away40_90hhError/away40_90hhYield, 2));
    Double_t mid4090 = mid40_90hPhiYield/mid40_90hhYield;
    Double_t mid4090Er = mid4090*TMath::Sqrt(TMath::Power(mid40_90hPhiError/mid40_90hPhiYield, 2) + TMath::Power(mid40_90hhError/mid40_90hhYield, 2));
    Double_t total4090 = total40_90hPhiYield/total40_90hhYield;
    Double_t total4090Er = total4090*TMath::Sqrt(TMath::Power(total40_90hPhiError/total40_90hPhiYield, 2) + TMath::Power(total40_90hhError/total40_90hhYield, 2));

    TH1D *ratios4090 = new TH1D("ratios4090", "(h-#phi / h-h) Ratios", 4, 0, 4);
    ratios4090->GetXaxis()->SetBinLabel(1, "near-side");
    ratios4090->SetBinContent(1, near4090);
    ratios4090->SetBinError(1, near4090Er);
    ratios4090->GetXaxis()->SetBinLabel(2, "mid");
    ratios4090->SetBinContent(2, mid4090);
    ratios4090->SetBinError(2, mid4090Er);
    ratios4090->GetXaxis()->SetBinLabel(3, "away-side");
    ratios4090->SetBinContent(3, away4090);
    ratios4090->SetBinError(3, away4090Er);
    ratios4090->GetXaxis()->SetBinLabel(4, "total");
    ratios4090->SetBinContent(4, total4090);
    ratios4090->SetBinError(4, total4090Er);
    ratios4090->SetMarkerStyle(20);
    ratios4090->SetMarkerColor(kGreen+2);
    ratios4090->SetLineColor(kGreen+2);
    ratios4090->SetMarkerSize(2);
    ratios4090->SetLineWidth(2);

    Double_t near0_100hPhiYield = near0_10hPhiYield + near10_20hPhiYield + near20_40hPhiYield + near40_90hPhiYield;
    Double_t near0_100hPhiError = TMath::Sqrt((TMath::Power(near0_10hPhiError,2) + TMath::Power(near10_20hPhiError,2) + TMath::Power(near20_40hPhiError,2) + TMath::Power(near40_90hPhiError,2)));
    Double_t near0_100hhYield = near0_10hhYield + near10_20hhYield + near20_40hhYield + near40_90hhYield;
    Double_t near0_100hhError = TMath::Sqrt((TMath::Power(near0_10hhError,2) + TMath::Power(near10_20hhError,2) + TMath::Power(near20_40hhError,2) + TMath::Power(near40_90hhError,2)));
    Double_t mid0_100hPhiYield = mid0_10hPhiYield + mid10_20hPhiYield + mid20_40hPhiYield + mid40_90hPhiYield;
    Double_t mid0_100hPhiError = TMath::Sqrt((TMath::Power(mid0_10hPhiError,2) + TMath::Power(mid10_20hPhiError,2) + TMath::Power(mid20_40hPhiError,2) + TMath::Power(mid40_90hPhiError,2)));
    Double_t mid0_100hhYield = mid0_10hhYield + mid10_20hhYield + mid20_40hhYield + mid40_90hhYield;
    Double_t mid0_100hhError = TMath::Sqrt((TMath::Power(mid0_10hhError,2) + TMath::Power(mid10_20hhError,2) + TMath::Power(mid20_40hhError,2) + TMath::Power(mid40_90hhError,2)));
    Double_t away0_100hPhiYield = away0_10hPhiYield + away10_20hPhiYield + away20_40hPhiYield + away40_90hPhiYield;
    Double_t away0_100hPhiError = TMath::Sqrt((TMath::Power(away0_10hPhiError,2) + TMath::Power(away10_20hPhiError,2) + TMath::Power(away20_40hPhiError,2) + TMath::Power(away40_90hPhiError,2)));
    Double_t away0_100hhYield = away0_10hhYield + away10_20hhYield + away20_40hhYield + away40_90hhYield;
    Double_t away0_100hhError = TMath::Sqrt((TMath::Power(away0_10hhError,2) + TMath::Power(away10_20hhError,2) + TMath::Power(away20_40hhError,2) + TMath::Power(away40_90hhError,2)));
    Double_t total0_100hPhiYield = total0_10hPhiYield + total10_20hPhiYield + total20_40hPhiYield + total40_90hPhiYield;
    Double_t total0_100hPhiError = TMath::Sqrt((TMath::Power(total0_10hPhiError,2) + TMath::Power(total10_20hPhiError,2) + TMath::Power(total20_40hPhiError,2) + TMath::Power(total40_90hPhiError,2)));
    Double_t total0_100hhYield = total0_10hhYield + total10_20hhYield + total20_40hhYield + total40_90hhYield;
    Double_t total0_100hhError = TMath::Sqrt((TMath::Power(total0_10hhError,2) + TMath::Power(total10_20hhError,2) + TMath::Power(total20_40hhError,2) + TMath::Power(total40_90hhError,2)));
   


    Double_t near0100 = (near0_10hPhiYield + near20_40hPhiYield + near40_90hPhiYield)/(near0_10hhYield + near20_40hhYield + near40_90hhYield);
    Double_t near0100Er = near0100*TMath::Sqrt((TMath::Power(near0_10hPhiError,2) + TMath::Power(near20_40hPhiError,2) + TMath::Power(near40_90hPhiError,2))/TMath::Power((near0_10hPhiYield + near20_40hPhiYield + near40_90hPhiYield), 2) + (TMath::Power(near0_10hhError,2) + TMath::Power(near20_40hhError,2) + TMath::Power(near40_90hhError,2))/TMath::Power((near0_10hhYield + near20_40hhYield + near40_90hhYield), 2));
    Double_t away0100 = (away0_10hPhiYield + away20_40hPhiYield + away40_90hPhiYield)/(away0_10hhYield + away20_40hhYield + away40_90hhYield);
    Double_t away0100Er = away0100*TMath::Sqrt((TMath::Power(away0_10hPhiError,2) + TMath::Power(away20_40hPhiError,2) + TMath::Power(away40_90hPhiError,2))/TMath::Power((away0_10hPhiYield + away20_40hPhiYield + away40_90hPhiYield), 2) + (TMath::Power(away0_10hhError,2) + TMath::Power(away20_40hhError,2) + TMath::Power(away40_90hhError,2))/TMath::Power((away0_10hhYield + away20_40hhYield + away40_90hhYield), 2));
    Double_t mid0100 = (mid0_10hPhiYield + mid20_40hPhiYield + mid40_90hPhiYield)/(mid0_10hhYield + mid20_40hhYield + mid40_90hhYield);
    Double_t mid0100Er = mid0100*TMath::Sqrt((TMath::Power(mid0_10hPhiError,2) + TMath::Power(mid20_40hPhiError,2) + TMath::Power(mid40_90hPhiError,2))/TMath::Power((mid0_10hPhiYield + mid20_40hPhiYield + mid40_90hPhiYield), 2) + (TMath::Power(mid0_10hhError,2) + TMath::Power(mid20_40hhError,2) + TMath::Power(mid40_90hhError,2))/TMath::Power((mid0_10hhYield + mid20_40hhYield + mid40_90hhYield), 2));
    Double_t total0100 = (total0_10hPhiYield + total20_40hPhiYield + total40_90hPhiYield)/(total0_10hhYield + total20_40hhYield + total40_90hhYield);
    Double_t total0100Er = total0100*TMath::Sqrt((TMath::Power(total0_10hPhiError,2) + TMath::Power(total20_40hPhiError,2) + TMath::Power(total40_90hPhiError,2))/TMath::Power((total0_10hPhiYield + total20_40hPhiYield + total40_90hPhiYield), 2) + (TMath::Power(total0_10hhError,2) + TMath::Power(total20_40hhError,2) + TMath::Power(total40_90hhError,2))/TMath::Power((total0_10hhYield + total20_40hhYield + total40_90hhYield), 2));
    
    printf("near010Er: %E\n", near010Er);
    printf("mid010Er: %E\n", mid010Er);
    printf("total0100Er: %E\n", total0100Er);
    printf("total010Er: %E\n", total010Er);
    printf("total0_10hPhiError: %E\n", total0_10hPhiError);
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
    ratioslegend->AddEntry(ratios010, "0-20%", "p");
    ratioslegend->AddEntry(ratios2040, "20-50%", "p");
    ratioslegend->AddEntry(ratios4090, "50-100%", "p");
    ratioslegend->AddEntry(ratios0100, "0-100%", "p");
    
    TLine *line = new TLine(3.0, 0.0, 3.0, 0.0040);
    line->SetLineStyle(7);
    line->SetLineWidth(2);
  
    TCanvas *testc = new TCanvas("test", "test",50, 50, 600, 600);
    testc->cd();
    ratios2040->GetYaxis()->SetRangeUser(0.0000, 0.0040);
    ratios2040->GetXaxis()->SetLabelSize(0.07);
    ratios2040->Draw("P SAME");
    ratios4090->GetXaxis()->SetLimits(0.05, 4.05);
    ratios4090->Draw("P SAME");
    ratios1020->GetXaxis()->SetLimits(0.1, 4.1);
    ratios1020->Draw("P SAME");
    ratios010->GetXaxis()->SetLimits(-0.05, 3.95);
    ratios010->Draw("P SAME");
    ratios0100->GetXaxis()->SetLimits(-0.1, 3.9);
    ratios0100->Draw("P SAME");
    line->Draw("SAME");
    ratioslegend->Draw("SAME");

    //Plot ratio as a function of multiplicity for the different angular regions
    Double_t nearArray[4] = {ratios4090->GetBinContent(1), ratios2040->GetBinContent(1), ratios1020->GetBinContent(1), ratios010->GetBinContent(1)};
    Double_t nearArrayErr[4] = {ratios4090->GetBinError(1), ratios2040->GetBinError(1), ratios1020->GetBinError(1), ratios010->GetBinError(1)};
    Double_t awayArray[4] = {ratios4090->GetBinContent(3), ratios2040->GetBinContent(3), ratios1020->GetBinContent(3), ratios010->GetBinContent(3)};
    Double_t awayArrayErr[4] = {ratios4090->GetBinError(3), ratios2040->GetBinError(3), ratios1020->GetBinError(3), ratios010->GetBinError(3)};
    Double_t bulkArray[4] = {ratios4090->GetBinContent(2), ratios2040->GetBinContent(2), ratios1020->GetBinContent(2), ratios010->GetBinContent(2)};
    Double_t bulkArrayErr[4] = {ratios4090->GetBinError(2), ratios2040->GetBinError(2), ratios1020->GetBinError(2),  ratios010->GetBinError(2)};
    Double_t totalArray[4] = {ratios4090->GetBinContent(4), ratios2040->GetBinContent(4), ratios1020->GetBinContent(4), ratios010->GetBinContent(4)};
    Double_t totalArrayErr[4] = {ratios4090->GetBinError(4), ratios2040->GetBinError(4), ratios1020->GetBinError(4), ratios010->GetBinError(4)};
    Double_t multArray[4] = {35.0, 70.0, 85.0, 95.0};
    Double_t multArrayErr[4] = {25.0, 10.0, 5.0, 5.0};

    Double_t mult2Array[4] = {36.0, 71.0, 86.0, 96.0};
    Double_t mult2ArrayErr[4] = {25.0, 10.0, 5.0, 5.0};

    //trying instead with variable sized histograms:
    Double_t binwidths[6] = {0.0, 10.0, 60.0, 80.0, 90.0, 100.0};
    TH1D* ratioNearHist = new TH1D("ratioNearHist", "", 4, binwidths);
    for(int i =1; i<5; i++){
        ratioNearHist->SetBinContent(i+1, nearArray[i-1]);
        ratioNearHist->SetBinError(i+1, nearArrayErr[i-1]);
    }
    ratioNearHist->SetMarkerStyle(20);
    ratioNearHist->SetMarkerSize(3);
    ratioNearHist->SetMarkerColor(kRed+1);
    ratioNearHist->SetLineColor(kRed+2);
    ratioNearHist->SetLineWidth(2);
    ratioNearHist->GetXaxis()->SetTitle("Multiplicity Pctl.");
    ratioNearHist->GetXaxis()->SetTitleSize(0.05);
    ratioNearHist->GetXaxis()->SetLabelSize(0.04);
    ratioNearHist->GetXaxis()->SetTitleOffset(1.2);
    ratioNearHist->GetXaxis()->SetRangeUser(0.0, 100.0);
    ratioNearHist->GetYaxis()->SetTitle("Un-corrected #frac{h-#phi}{h-h} Ratio");
    ratioNearHist->GetYaxis()->SetTitleSize(0.04);
    ratioNearHist->GetYaxis()->SetTitleOffset(1.5); 
    ratioNearHist->GetYaxis()->SetRangeUser(0.0002, 0.0035);



    TGraphErrors* ratiosNear = new TGraphErrors(4, multArray, nearArray, multArrayErr, nearArrayErr);
    ratiosNear->SetMarkerStyle(20);
    ratiosNear->SetMarkerSize(3);
    ratiosNear->SetMarkerColor(kRed+1);
    ratiosNear->SetLineColor(kRed+2);
    ratiosNear->SetLineWidth(2);
    ratiosNear->GetXaxis()->SetTitle("Multiplicity Pctl.");
    ratiosNear->GetXaxis()->SetTitleSize(0.05);
    ratiosNear->GetXaxis()->SetLabelSize(0.04);
    ratiosNear->GetXaxis()->SetTitleOffset(0.9);
    ratiosNear->GetXaxis()->SetRangeUser(0.0, 100.0);
    ratiosNear->GetYaxis()->SetTitle("Un-corrected #frac{h-#phi}{h-h} Ratio");
    ratiosNear->GetYaxis()->SetTitleSize(0.04);
    ratiosNear->GetYaxis()->SetTitleOffset(1.5); 
    ratiosNear->GetYaxis()->SetRangeUser(0.0002, 0.0035);


    TGraphErrors* ratiosAway = new TGraphErrors(4, mult2Array, awayArray, mult2ArrayErr, awayArrayErr);
    ratiosAway->SetMarkerStyle(21);
    ratiosAway->SetMarkerSize(3);
    ratiosAway->SetMarkerColor(kBlue+1);
    ratiosAway->SetLineColor(kBlue+2);
    ratiosAway->SetLineWidth(2);

    TGraphErrors* ratiosBulk = new TGraphErrors(4, multArray, bulkArray, multArrayErr, bulkArrayErr);
    ratiosBulk->SetMarkerStyle(22);
    ratiosBulk->SetMarkerSize(3);
    ratiosBulk->SetMarkerColor(kGreen+2);
    ratiosBulk->SetLineColor(kGreen+3);
    ratiosBulk->SetLineWidth(2);
   
    TGraphErrors* ratiosTot = new TGraphErrors(4, multArray, totalArray, multArrayErr, totalArrayErr);
    ratiosTot->SetMarkerStyle(29);
    ratiosTot->SetMarkerSize(3);
    ratiosTot->SetMarkerColor(kMagenta+2);
    ratiosTot->SetLineColor(kMagenta+3);
    ratiosTot->SetLineWidth(2);
    ratiosTot->SetFillColor(kMagenta+1);
    ratiosTot->SetFillStyle(3144);

    TLegend  *ratiosMultlegend = new TLegend(0.183, 0.686, 0.461, 0.928);
    ratiosMultlegend->SetMargin(0.35);
    ratiosMultlegend->AddEntry(ratiosNear, "In Near-side Jet", "pl");
    ratiosMultlegend->AddEntry(ratiosAway, "In Away-side Jet", "pl");
    //ratiosMultlegend->AddEntry(ratiosBulk, "In U.E.", "pl");
    ratiosMultlegend->AddEntry(ratiosTot, "Total (Jet + UE)", "f");

    
    TCanvas* vsMultCanvas = new TCanvas("vsMultCanvas", "vsMultCanvas", 55, 55, 900, 600);
    vsMultCanvas->cd();
    vsMultCanvas->SetMargin(0.126, 0.05, 0.125, 0.05);
    //TH1F* hist = ratiosNear->GetHistogram();
    gStyle->SetErrorX(0.5);
    ratioNearHist->Draw("PE");

    ratioNearHist->GetXaxis()->SetLabelOffset(999);
    //ratioNearHist->GetXaxis()->SetTitleOffset(999);
    ratioNearHist->GetXaxis()->SetTickSize(0.0);

    //ratiosNear->Draw("P");
    gPad->Update();
    TGaxis *newaxis = new TGaxis(gPad->GetUxmax(),
            gPad->GetUymin(),
            gPad->GetUxmin(),
            gPad->GetUymin(),
            0,
            100,
            510,"-");
    newaxis->SetLabelOffset(-0.03);
    //newaxis->SetTitle("Multipliciy % (VOA)");
    //newaxis->SetTitleOffset(1.3);
    newaxis->Draw();   

    ratiosAway->Draw("P");
    
    ratiosBulk->Draw("P");
    
    ratiosTot->Draw("2");
    
    //ratiosTot->Draw("3");
    ratiosMultlegend->Draw();
    //ratiosNear->Draw("PL");
    //newaxis->Draw();
    //gPad->Update();

    //Double Ratio plots:
    Double_t doublenear010 = near010/mid010;
    Double_t doublenear010Er = doublenear010*TMath::Sqrt(TMath::Power((near010Er)/(near010), 2) + TMath::Power((mid010Er)/(mid010), 2));
    Double_t doubleaway010 = away010/mid010;
    Double_t doubleaway010Er = doubleaway010*TMath::Sqrt(TMath::Power((away010Er)/(away010), 2) + TMath::Power((mid010Er)/(mid010), 2));

    Double_t doublenear2040 = near2040/mid2040;
    Double_t doublenear2040Er = doublenear2040*TMath::Sqrt(TMath::Power((near2040Er)/(near2040), 2) + TMath::Power((mid2040Er)/(mid2040), 2));
    Double_t doubleaway2040 = away2040/mid2040;
    Double_t doubleaway2040Er = doubleaway2040*TMath::Sqrt(TMath::Power((away2040Er)/(away2040), 2) + TMath::Power((mid2040Er)/(mid2040), 2));

    Double_t doublenear4090 = near4090/mid4090;
    Double_t doublenear4090Er = doublenear4090*TMath::Sqrt(TMath::Power((near4090Er)/(near4090), 2) + TMath::Power((mid4090Er)/(mid4090), 2));
    Double_t doubleaway4090 = away4090/mid4090;
    Double_t doubleaway4090Er = doubleaway4090*TMath::Sqrt(TMath::Power((away4090Er)/(away4090), 2) + TMath::Power((mid4090Er)/(mid4090), 2));

    Double_t doublenear0100 = near0100/mid0100;
    Double_t doublenear0100Er = doublenear0100*TMath::Sqrt(TMath::Power((near0100Er)/(near0100), 2) + TMath::Power((mid0100Er)/(mid0100), 2));
    Double_t doubleaway0100 = away0100/mid0100;
    Double_t doubleaway0100Er = doubleaway0100*TMath::Sqrt(TMath::Power((away0100Er)/(away0100), 2) + TMath::Power((mid0100Er)/(mid0100), 2));

    TH1D *doubleratios010 = new TH1D("doubleratios010", "(h-#phi / h-h) Double Ratios", 2, 0, 2);
    doubleratios010->GetXaxis()->SetBinLabel(1, "#frac{near-side}{mid}");
    doubleratios010->SetBinContent(1, doublenear010);
    doubleratios010->SetBinError(1, doublenear010Er);
    doubleratios010->GetXaxis()->SetBinLabel(2, "#frac{away-side}{mid}");
    doubleratios010->SetBinContent(2, doubleaway010);
    doubleratios010->SetBinError(2, doubleaway010Er);
    doubleratios010->SetMarkerStyle(22);
    doubleratios010->SetMarkerColor(kRed+2);
    doubleratios010->SetLineColor(kRed+2);
    doubleratios010->SetLineWidth(2);
    doubleratios010->SetMarkerSize(2);

    TH1D *doubleratios2040 = new TH1D("doubleratios2040", "(h-#phi / h-h) Double Ratios", 2, 0, 2);
    doubleratios2040->GetXaxis()->SetBinLabel(1, "#frac{near-side}{mid}");
    doubleratios2040->SetBinContent(1, doublenear2040);
    doubleratios2040->SetBinError(1, doublenear2040Er);
    doubleratios2040->GetXaxis()->SetBinLabel(2, "#frac{away-side}{mid}");
    doubleratios2040->SetBinContent(2, doubleaway2040);
    doubleratios2040->SetBinError(2, doubleaway2040Er);
    doubleratios2040->SetMarkerStyle(21);
    doubleratios2040->SetMarkerColor(kBlue+2);
    doubleratios2040->SetLineColor(kBlue+2);
    doubleratios2040->SetLineWidth(2);
    doubleratios2040->SetMarkerSize(2);

    TH1D *doubleratios4090 = new TH1D("doubleratios4090", "(h-#phi / h-h) Double Ratios", 2, 0, 2);
    doubleratios4090->GetXaxis()->SetBinLabel(1, "#frac{near-side}{mid}");
    doubleratios4090->SetBinContent(1, doublenear4090);
    doubleratios4090->SetBinError(1, doublenear4090Er);
    doubleratios4090->GetXaxis()->SetBinLabel(2, "#frac{away-side}{mid}");
    doubleratios4090->SetBinContent(2, doubleaway4090);
    doubleratios4090->SetBinError(2, doubleaway4090Er);
    doubleratios4090->SetMarkerStyle(20);
    doubleratios4090->SetMarkerColor(kGreen+2);
    doubleratios4090->SetLineColor(kGreen+2);
    doubleratios4090->SetLineWidth(2);
    doubleratios4090->SetMarkerSize(2);

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
    //doubleratios2040->GetYaxis()->SetRangeUser(0.0000, 0.0040);
    doubleratios2040->GetXaxis()->SetLabelSize(0.07);
    doubleratios2040->Draw("P SAME");
    doubleratios4090->GetXaxis()->SetLimits(-0.05, 1.95);
    doubleratios4090->Draw("P SAME");
    doubleratios010->GetXaxis()->SetLimits(0.05, 2.05);
    doubleratios010->Draw("P SAME");
    doubleratios0100->GetXaxis()->SetLimits(-0.1,1.9);
    doubleratios0100->Draw("P SAME");
    //line->Draw("SAME");
    ratioslegend->Draw("SAME");

    printf("\n\n");
    printf(" h-Phi YIELDS  ||     0 - 20     ||    20 - 50     ||    50 - 100    ||     0 - 100    ||\n");
    printf("=========================================================================================\n");
    printf("      NEAR     ||  %E  ||  %E  ||  %E  ||  %E  ||\n", near0_10hPhiYield, near20_40hPhiYield, near40_90hPhiYield, near0_100hPhiYield);
    printf("      MID      ||  %E  ||  %E  ||  %E  ||  %E  ||\n", mid0_10hPhiYield, mid20_40hPhiYield, mid40_90hPhiYield, mid0_100hPhiYield);
    printf("      AWAY     ||  %E  ||  %E  ||  %E  ||  %E  ||\n", away0_10hPhiYield, away20_40hPhiYield, away40_90hPhiYield, away0_100hPhiYield);
    printf("      TOTAL    ||  %E  ||  %E  ||  %E  ||  %E  ||\n", total0_10hPhiYield, total20_40hPhiYield, total40_90hPhiYield, total0_100hPhiYield);

    printf("\n\n");
    printf(" h-Phi ERRORS  ||     0 - 20     ||    20 - 50     ||    50 - 100    ||     0 - 100    ||\n");
    printf("=========================================================================================\n");
    printf("      NEAR     ||  %E  ||  %E  ||  %E  ||  %E  ||\n", near0_10hPhiError, near20_40hPhiError, near40_90hPhiError, near0_100hPhiError);
    printf("      MID      ||  %E  ||  %E  ||  %E  ||  %E  ||\n", mid0_10hPhiError, mid20_40hPhiError, mid40_90hPhiError, mid0_100hPhiError);
    printf("      AWAY     ||  %E  ||  %E  ||  %E  ||  %E  ||\n", away0_10hPhiError, away20_40hPhiError, away40_90hPhiError, away0_100hPhiError);
    printf("      TOTAL    ||  %E  ||  %E  ||  %E  ||  %E  ||\n", total0_10hPhiError, total20_40hPhiError, total40_90hPhiError, total0_100hPhiError);

    printf("\n\n");
    printf(" h-Phi %%ERRORS ||     0 - 20     ||    20 - 50     ||    50 - 100    ||     0 - 100    ||\n");
    printf("=========================================================================================\n");
    printf("      NEAR     ||     %2.2f%%     ||     %2.2f%%     ||     %2.2f%%     ||     %2.2f%%     ||\n", 100*near0_10hPhiError/near0_10hPhiYield, 100*near20_40hPhiError/near20_40hPhiYield, 100*near40_90hPhiError/near40_90hPhiYield, 100*near0_100hPhiError/near0_100hPhiYield);
    printf("      MID      ||      %2.2f%%     ||      %2.2f%%     ||      %2.2f%%     ||      %2.2f%%     ||\n", 100*mid0_10hPhiError/mid0_10hPhiYield, 100*mid20_40hPhiError/mid20_40hPhiYield, 100*mid40_90hPhiError/mid40_90hPhiYield, 100*mid0_100hPhiError/mid0_100hPhiYield);
    printf("      AWAY     ||     %2.2f%%     ||     %2.2f%%     ||     %2.2f%%     ||     %2.2f%%     ||\n", 100*away0_10hPhiError/away0_10hPhiYield, 100*away20_40hPhiError/away20_40hPhiYield, 100*away40_90hPhiError/away40_90hPhiYield, 100*away0_100hPhiError/away0_100hPhiYield);
    printf("      TOTAL    ||      %2.2f%%     ||      %2.2f%%     ||      %2.2f%%     ||      %2.2f%%     ||\n", 100*total0_10hPhiError/total0_10hPhiYield, 100*total20_40hPhiError/total20_40hPhiYield, 100*total40_90hPhiError/total40_90hPhiYield, 100*total0_100hPhiError/total0_100hPhiYield);
   


    printf("\n\n");
    printf(" h-h   YIELDS  ||     0 - 20     ||    20 - 50     ||    50 - 100    ||     0 - 100    ||\n");
    printf("=========================================================================================\n");
    printf("      NEAR     ||  %E  ||  %E  ||  %E  ||  %E  ||\n", near0_10hhYield, near20_40hhYield, near40_90hhYield, near0_100hhYield);
    printf("      MID      ||  %E  ||  %E  ||  %E  ||  %E  ||\n", mid0_10hhYield, mid20_40hhYield, mid40_90hhYield, mid0_100hhYield);
    printf("      AWAY     ||  %E  ||  %E  ||  %E  ||  %E  ||\n", away0_10hhYield, away20_40hhYield, away40_90hhYield, away0_100hhYield);
    printf("      TOTAL    ||  %E  ||  %E  ||  %E  ||  %E  ||\n", total0_10hhYield, total20_40hhYield, total40_90hhYield, total0_100hhYield);

    printf("\n\n");
    printf(" h-h ERRORS  ||     0 - 20     ||    20 - 50     ||    50 - 100    ||     0 - 100    ||\n");
    printf("=========================================================================================\n");
    printf("      NEAR     ||  %E  ||  %E  ||  %E  ||  %E  ||\n", near0_10hhError, near20_40hhError, near40_90hhError, near0_100hhError);
    printf("      MID      ||  %E  ||  %E  ||  %E  ||  %E  ||\n", mid0_10hhError, mid20_40hhError, mid40_90hhError, mid0_100hhError);
    printf("      AWAY     ||  %E  ||  %E  ||  %E  ||  %E  ||\n", away0_10hhError, away20_40hhError, away40_90hhError, away0_100hhError);
    printf("      TOTAL    ||  %E  ||  %E  ||  %E  ||  %E  ||\n", total0_10hhError, total20_40hhError, total40_90hhError, total0_100hhError);
    printf("\n\n");


    TH1D *yields010hPhi = new TH1D("yields010hPhi", "h-#phi Per Trigger Yields", 4, 0, 4);
    yields010hPhi->GetXaxis()->SetBinLabel(1, "near-side");
    yields010hPhi->SetBinContent(1, near0_10hPhiYield);
    yields010hPhi->SetBinError(1, near0_10hPhiError);
    yields010hPhi->GetXaxis()->SetBinLabel(2, "mid");
    yields010hPhi->SetBinContent(2, mid0_10hPhiYield);
    yields010hPhi->SetBinError(2, mid0_10hPhiError);
    yields010hPhi->GetXaxis()->SetBinLabel(3, "away-side");
    yields010hPhi->SetBinContent(3, away0_10hPhiYield);
    yields010hPhi->SetBinError(3, away0_10hPhiError);
    yields010hPhi->GetXaxis()->SetBinLabel(4, "total");
    yields010hPhi->SetBinContent(4, total0_10hPhiYield);
    yields010hPhi->SetBinError(4, total0_10hPhiError);
    yields010hPhi->SetMarkerStyle(22);
    yields010hPhi->SetMarkerColor(kRed+2);
    yields010hPhi->SetLineColor(kRed+2);
    yields010hPhi->SetMarkerSize(2);
    yields010hPhi->SetLineWidth(2);

    TH1D *yields2040hPhi = new TH1D("yields2040hPhi", "h-#phi Per Trigger Yields", 4, 0, 4);
    yields2040hPhi->GetXaxis()->SetBinLabel(1, "near-side");
    yields2040hPhi->SetBinContent(1, near20_40hPhiYield);
    yields2040hPhi->SetBinError(1, near20_40hPhiError);
    yields2040hPhi->GetXaxis()->SetBinLabel(2, "mid");
    yields2040hPhi->SetBinContent(2, mid20_40hPhiYield);
    yields2040hPhi->SetBinError(2, mid20_40hPhiError);
    yields2040hPhi->GetXaxis()->SetBinLabel(3, "away-side");
    yields2040hPhi->SetBinContent(3, away20_40hPhiYield);
    yields2040hPhi->SetBinError(3, away20_40hPhiError);
    yields2040hPhi->GetXaxis()->SetBinLabel(4, "total");
    yields2040hPhi->SetBinContent(4, total20_40hPhiYield);
    yields2040hPhi->SetBinError(4, total20_40hPhiError);
    yields2040hPhi->SetMarkerStyle(21);
    yields2040hPhi->SetMarkerColor(kBlue+2);
    yields2040hPhi->SetLineColor(kBlue+2);
    yields2040hPhi->SetMarkerSize(2);
    yields2040hPhi->SetLineWidth(2);

    TH1D *yields4090hPhi = new TH1D("yields4090hPhi", "h-#phi Per Trigger Yields", 4, 0, 4);
    yields4090hPhi->GetXaxis()->SetBinLabel(1, "near-side");
    yields4090hPhi->SetBinContent(1, near40_90hPhiYield);
    yields4090hPhi->SetBinError(1, near40_90hPhiError);
    yields4090hPhi->GetXaxis()->SetBinLabel(2, "mid");
    yields4090hPhi->SetBinContent(2, mid40_90hPhiYield);
    yields4090hPhi->SetBinError(2, mid40_90hPhiError);
    yields4090hPhi->GetXaxis()->SetBinLabel(3, "away-side");
    yields4090hPhi->SetBinContent(3, away40_90hPhiYield);
    yields4090hPhi->SetBinError(3, away40_90hPhiError);
    yields4090hPhi->GetXaxis()->SetBinLabel(4, "total");
    yields4090hPhi->SetBinContent(4, total40_90hPhiYield);
    yields4090hPhi->SetBinError(4, total40_90hPhiError);
    yields4090hPhi->SetMarkerStyle(20);
    yields4090hPhi->SetMarkerColor(kGreen+2);
    yields4090hPhi->SetLineColor(kGreen+2);
    yields4090hPhi->SetMarkerSize(2);
    yields4090hPhi->SetLineWidth(2);

    TH1D *yields0100hPhi = new TH1D("yields0100hPhi", "h-#phi Per Trigger Yields", 4, 0, 4);
    yields0100hPhi->GetXaxis()->SetBinLabel(1, "near-side");
    yields0100hPhi->SetBinContent(1, near0_100hPhiYield);
    yields0100hPhi->SetBinError(1, near0_100hPhiError);
    yields0100hPhi->GetXaxis()->SetBinLabel(2, "mid");
    yields0100hPhi->SetBinContent(2, mid0_100hPhiYield);
    yields0100hPhi->SetBinError(2, mid0_100hPhiError);
    yields0100hPhi->GetXaxis()->SetBinLabel(3, "away-side");
    yields0100hPhi->SetBinContent(3, away0_100hPhiYield);
    yields0100hPhi->SetBinError(3, away0_100hPhiError);
    yields0100hPhi->GetXaxis()->SetBinLabel(4, "total");
    yields0100hPhi->SetBinContent(4, total0_100hPhiYield);
    yields0100hPhi->SetBinError(4, total0_100hPhiError);
    yields0100hPhi->SetMarkerStyle(29);
    yields0100hPhi->SetMarkerColor(kViolet-1);
    yields0100hPhi->SetLineColor(kViolet-1);
    yields0100hPhi->SetMarkerSize(2);
    yields0100hPhi->SetLineWidth(2);


    TCanvas *yieldhPhic = new TCanvas("yieldhPhi", "yieldhPhi",50, 50, 600, 600);
    yieldhPhic->cd();
    yields2040hPhi->GetYaxis()->SetRangeUser(0.0000, 0.0040);
    yields2040hPhi->GetXaxis()->SetLabelSize(0.07);
    yields2040hPhi->Draw("P SAME");
    yields4090hPhi->GetXaxis()->SetLimits(0.05, 4.05);
    yields4090hPhi->Draw("P SAME");
    yields010hPhi->GetXaxis()->SetLimits(-0.05, 3.95);
    yields010hPhi->Draw("P SAME");
    yields0100hPhi->GetXaxis()->SetLimits(-0.1, 3.9);
    yields0100hPhi->Draw("P SAME");
    line->Draw("SAME");
    ratioslegend->Draw("SAME");


    TH1D *yields010hh = new TH1D("yields010hh", "h-#phi Per Trigger Yields", 4, 0, 4);
    yields010hh->GetXaxis()->SetBinLabel(1, "near-side");
    yields010hh->SetBinContent(1, near0_10hhYield);
    yields010hh->SetBinError(1, near0_10hhError);
    yields010hh->GetXaxis()->SetBinLabel(2, "mid");
    yields010hh->SetBinContent(2, mid0_10hhYield);
    yields010hh->SetBinError(2, mid0_10hhError);
    yields010hh->GetXaxis()->SetBinLabel(3, "away-side");
    yields010hh->SetBinContent(3, away0_10hhYield);
    yields010hh->SetBinError(3, away0_10hhError);
    yields010hh->GetXaxis()->SetBinLabel(4, "total");
    yields010hh->SetBinContent(4, total0_10hhYield);
    yields010hh->SetBinError(4, total0_10hhError);
    yields010hh->SetMarkerStyle(22);
    yields010hh->SetMarkerColor(kRed+2);
    yields010hh->SetLineColor(kRed+2);
    yields010hh->SetMarkerSize(2);
    yields010hh->SetLineWidth(2);

    TH1D *yields2040hh = new TH1D("yields2040hh", "h-h Per Trigger Yields", 4, 0, 4);
    yields2040hh->GetXaxis()->SetBinLabel(1, "near-side");
    yields2040hh->SetBinContent(1, near20_40hhYield);
    yields2040hh->SetBinError(1, near20_40hhError);
    yields2040hh->GetXaxis()->SetBinLabel(2, "mid");
    yields2040hh->SetBinContent(2, mid20_40hhYield);
    yields2040hh->SetBinError(2, mid20_40hhError);
    yields2040hh->GetXaxis()->SetBinLabel(3, "away-side");
    yields2040hh->SetBinContent(3, away20_40hhYield);
    yields2040hh->SetBinError(3, away20_40hhError);
    yields2040hh->GetXaxis()->SetBinLabel(4, "total");
    yields2040hh->SetBinContent(4, total20_40hhYield);
    yields2040hh->SetBinError(4, total20_40hhError);
    yields2040hh->SetMarkerStyle(21);
    yields2040hh->SetMarkerColor(kBlue+2);
    yields2040hh->SetLineColor(kBlue+2);
    yields2040hh->SetMarkerSize(2);
    yields2040hh->SetLineWidth(2);

    TH1D *yields4090hh = new TH1D("yields4090hh", "h-#phi Per Trigger Yields", 4, 0, 4);
    yields4090hh->GetXaxis()->SetBinLabel(1, "near-side");
    yields4090hh->SetBinContent(1, near40_90hhYield);
    yields4090hh->SetBinError(1, near40_90hhError);
    yields4090hh->GetXaxis()->SetBinLabel(2, "mid");
    yields4090hh->SetBinContent(2, mid40_90hhYield);
    yields4090hh->SetBinError(2, mid40_90hhError);
    yields4090hh->GetXaxis()->SetBinLabel(3, "away-side");
    yields4090hh->SetBinContent(3, away40_90hhYield);
    yields4090hh->SetBinError(3, away40_90hhError);
    yields4090hh->GetXaxis()->SetBinLabel(4, "total");
    yields4090hh->SetBinContent(4, total40_90hhYield);
    yields4090hh->SetBinError(4, total40_90hhError);
    yields4090hh->SetMarkerStyle(20);
    yields4090hh->SetMarkerColor(kGreen+2);
    yields4090hh->SetLineColor(kGreen+2);
    yields4090hh->SetMarkerSize(2);
    yields4090hh->SetLineWidth(2);

    TH1D *yields0100hh = new TH1D("yields0100hh", "h-#phi Per Trigger Yields", 4, 0, 4);
    yields0100hh->GetXaxis()->SetBinLabel(1, "near-side");
    yields0100hh->SetBinContent(1, near0_100hhYield);
    yields0100hh->SetBinError(1, near0_100hhError);
    yields0100hh->GetXaxis()->SetBinLabel(2, "mid");
    yields0100hh->SetBinContent(2, mid0_100hhYield);
    yields0100hh->SetBinError(2, mid0_100hhError);
    yields0100hh->GetXaxis()->SetBinLabel(3, "away-side");
    yields0100hh->SetBinContent(3, away0_100hhYield);
    yields0100hh->SetBinError(3, away0_100hhError);
    yields0100hh->GetXaxis()->SetBinLabel(4, "total");
    yields0100hh->SetBinContent(4, total0_100hhYield);
    yields0100hh->SetBinError(4, total0_100hhError);
    yields0100hh->SetMarkerStyle(29);
    yields0100hh->SetMarkerColor(kViolet-1);
    yields0100hh->SetLineColor(kViolet-1);
    yields0100hh->SetMarkerSize(2);
    yields0100hh->SetLineWidth(2);

    TLine *linehh = new TLine(3.0, 0.0, 3.0, 10.000);
    linehh->SetLineStyle(7);
    linehh->SetLineWidth(2);

    TCanvas *yieldhhc = new TCanvas("yieldhh", "yieldhh",50, 50, 600, 600);
    yieldhhc->cd();
    yields2040hh->GetYaxis()->SetRangeUser(0.0000, 10.000);
    yields2040hh->GetXaxis()->SetLabelSize(0.07);
    yields2040hh->Draw("P SAME");
    yields4090hh->GetXaxis()->SetLimits(0.05, 4.05);
    yields4090hh->Draw("P SAME");
    yields010hh->GetXaxis()->SetLimits(-0.05, 3.95);
    yields010hh->Draw("P SAME");
    yields0100hh->GetXaxis()->SetLimits(-0.1, 3.9);
    yields0100hh->Draw("P SAME");
    linehh->Draw("SAME");
    ratioslegend->Draw("SAME");



    
    TH1D *ratio = (TH1D*)hPhidphi_0_10->Clone("ratio");
    ratio->Divide(hhdphi_0_10);
    

    TLegend  *legend = new TLegend(0.3791, 0.1518, 0.8772, 0.2688);
    legend->SetMargin(0.15);
    legend->AddEntry(corrFit2, "Hadron-#phi(1020) Correlation", "l");
    legend->AddEntry(corrFit, "Hadron-hadron Correlations", "l");

    TPaveText *text = new TPaveText(0.4815, 0.7056, 0.8658, 0.8551, "NDC");
    text->AddText("ALICE Work in Progress");
    text->AddText("p-Pb #sqrt{s_{NN}} = 5 TeV");
    text->AddText("0%-10% Multiplicity");
    text->SetTextSizePixels(20);
    text->SetBorderSize(0);
    text->SetFillColor(kWhite);

    TPaveText *text1020 = new TPaveText(0.4815, 0.7056, 0.8658, 0.8551, "NDC");
    text1020->AddText("ALICE Work in Progress");
    text1020->AddText("p-Pb #sqrt{s_{NN}} = 5 TeV");
    text1020->AddText("10%-20% Multiplicity");
    text1020->SetTextSizePixels(20);
    text1020->SetBorderSize(0);
    text1020->SetFillColor(kWhite);

    TPaveText *text2040 = new TPaveText(0.4815, 0.7056, 0.8658, 0.8551, "NDC");
    text2040->AddText("ALICE Work in Progress");
    text2040->AddText("p-Pb #sqrt{s_{NN}} = 5 TeV");
    text2040->AddText("20%-40% Multiplicity");
    text2040->SetTextSizePixels(20);
    text2040->SetBorderSize(0);
    text2040->SetFillColor(kWhite);

    TPaveText *text4090 = new TPaveText(0.4815, 0.7056, 0.8658, 0.8551, "NDC");
    text4090->AddText("ALICE Work in Progress");
    text4090->AddText("p-Pb #sqrt{s_{NN}} = 5 TeV");
    text4090->AddText("40%-90% Multiplicity");
    text4090->SetBorderSize(0);
    text4090->SetTextSizePixels(20);
    text4090->SetFillColor(kWhite);

    
    TPaveText *text2 = new TPaveText(0.6, 0.9, 0.85, 0.85, "NDC");
    text2->AddText("trigger: 4.0 < p_{T}^{h} < 8.0 GeV/c");
    text2->AddText("assoc: 2.0 < p_{T}^{#phi} < 4.0 GeV/c");
    text2->SetTextSizePixels(18);
    text2->SetFillColor(kWhite);
    text2->SetBorderSize(0);

    
    TCanvas *c0_10 = new TCanvas("c0_10", "c0_10", 50, 50, 550, 600);
    c0_10->cd();
    c0_10->SetMargin(0.12, 0.05, 0.1, 0.05);
    hhdphi_0_10->GetYaxis()->SetRangeUser(0.001, 0.25);
    //hhdphi_0_10->Draw("E0 X0");
    hPhidphi_0_10->Draw("E0 X0 SAME");
    //corrFit->Draw("SAME");
    corrFit2->Draw("SAME");
    //hhBG->Draw("SAME");
    hphiBG->Draw("SAME");
    text->Draw();
    text2->Draw();
    legend->Draw();

    TCanvas *c10_20 = new TCanvas("c10_20", "c10_20", 50, 50, 550, 600);
    c10_20->cd();
    c10_20->SetMargin(0.12, 0.05, 0.1, 0.05);
    hhdphi_0_10->GetYaxis()->SetRangeUser(0.03, 0.11);
    //hhdphi_10_20->Draw("E0 X0");
    hPhidphi_10_20->Draw("E0 X0 SAME");
    //corrFit1020->Draw("SAME");
    corrFit2_10_20->Draw("SAME");
    //hhBG_10_20->Draw("SAME");
    hphiBG_10_20->Draw("SAME");
    text1020->Draw();
    text2->Draw();
    legend->Draw();


    TCanvas *c20_40 = new TCanvas("c20_40", "c20_40", 50, 50, 550, 600);
    c20_40->cd();
    c20_40->SetMargin(0.12, 0.05, 0.1, 0.05);
    hhdphi_0_10->GetYaxis()->SetRangeUser(0.03, 0.11);
    //hhdphi_20_40->Draw("E0 X0");
    hPhidphi_20_40->Draw("E0 X0 SAME");
    //corrFit2040->Draw("SAME");
    corrFit2_2040->Draw("SAME");
    //hhBG_20_40->Draw("SAME");
    hphiBG_20_40->Draw("SAME");
    text2040->Draw();
    text2->Draw();
    legend->Draw();

    TCanvas *c40_90 = new TCanvas("c40_90", "c40_90", 50, 50, 550, 600);
    c40_90->cd();
    c40_90->SetMargin(0.12, 0.05, 0.1, 0.05);
    //hhdphi_40_90->Draw("E0 X0");
    hPhidphi_40_90->Draw("E0 X0 SAME");
    //corrFit4090->Draw("SAME");
    corrFit2_4090->Draw("SAME");
    //hhBG_40_90->Draw("SAME");
    hphiBG_40_90->Draw("SAME");
    text4090->Draw();
    text2->Draw();
    legend->Draw();
/*
    TCanvas *cratio = new TCanvas("cratio", "cratio", 50, 50, 550, 600);
    cratio->cd();
    cratio->SetMargin(0.12, 0.05, 0.1, 0.05);
    ratio->Draw("E0 X0");
*/

}

