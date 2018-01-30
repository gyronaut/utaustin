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

    float near = hPhidphi_0_20->Integral(2,7);
    near  = near - hphiBG->GetParameter(0)*6.0;
    near = near/(hhdphi_0_20->Integral(2,7) - hhBG->GetParameter(0)*6.0);

    float away = hPhidphi_0_20->Integral(9,16);
    away  = away - hphiBG->GetParameter(0)*8.0;
    away = away/(hhdphi_0_20->Integral(9,16) - hhBG->GetParameter(0)*8.0);

    float mid = hphiBG->GetParameter(0)/hhBG->GetParameter(0);

    TH1D *ratios = new TH1D("ratios", "(h-#phi / h-h) Ratios", 3, 0, 3);
    ratios->GetXaxis()->SetBinLabel(1, "near-side");
    ratios->SetBinContent(1, near);
    ratios->GetXaxis()->SetBinLabel(2, "mid");
    ratios->SetBinContent(2, mid);
    ratios->GetXaxis()->SetBinLabel(3, "away-side");
    ratios->SetBinContent(3, away);
    ratios->SetMarkerStyle(22);
    ratios->SetMarkerColor(kRed+2);
    ratios->SetMarkerSize(4);

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

    float near2050 = hPhidphi_20_50->Integral(2,7);
    near2050  = near2050 - hphiBG_20_50->GetParameter(0)*6.0;
    near2050 = near2050/(hhdphi_20_50->Integral(2,7) - hhBG_20_50->GetParameter(0)*6.0);

    float away2050 = hPhidphi_20_50->Integral(9,16);
    away2050  = away2050 - hphiBG_20_50->GetParameter(0)*8.0;
    away2050 = away2050/(hhdphi_20_50->Integral(9,16) - hhBG_20_50->GetParameter(0)*8.0);

    float mid2050 = hphiBG_20_50->GetParameter(0)/hhBG_20_50->GetParameter(0);

    TH1D *ratios2050 = new TH1D("ratios2050", "(h-#phi / h-h) Ratios", 3, 0, 3);
    ratios2050->GetXaxis()->SetBinLabel(1, "near-side");
    ratios2050->SetBinContent(1, near2050);
    ratios2050->GetXaxis()->SetBinLabel(2, "mid");
    ratios2050->SetBinContent(2, mid2050);
    ratios2050->GetXaxis()->SetBinLabel(3, "away-side");
    ratios2050->SetBinContent(3, away2050);
    ratios2050->SetMarkerStyle(21);
    ratios2050->SetMarkerColor(kBlue+2);
    ratios2050->SetMarkerSize(4);


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

    float near50100 = hPhidphi_50_100->Integral(2,7);
    near50100  = near50100 - hphiBG_50_100->GetParameter(0)*6.0;
    near50100 = near50100/(hhdphi_50_100->Integral(2,7) - hhBG_50_100->GetParameter(0)*6.0);

    float away50100 = hPhidphi_50_100->Integral(9,16);
    away50100  = away50100 - hphiBG_50_100->GetParameter(0)*8.0;
    away50100 = away50100/(hhdphi_50_100->Integral(9,16) - hhBG_50_100->GetParameter(0)*8.0);

    float mid50100 = hphiBG_50_100->GetParameter(0)/hhBG_50_100->GetParameter(0);

    TH1D *ratios50100 = new TH1D("ratios50100", "(h-#phi / h-h) Ratios", 3, 0, 3);
    ratios50100->GetXaxis()->SetBinLabel(1, "near-side");
    ratios50100->SetBinContent(1, near50100);
    ratios50100->GetXaxis()->SetBinLabel(2, "mid");
    ratios50100->SetBinContent(2, mid50100);
    ratios50100->GetXaxis()->SetBinLabel(3, "away-side");
    ratios50100->SetBinContent(3, away50100);
    ratios50100->SetMarkerStyle(20);
    ratios50100->SetMarkerColor(kGreen+2);
    ratios50100->SetMarkerSize(4);


    TLegend  *ratioslegend = new TLegend(0.3791, 0.1518, 0.8772, 0.2688);
    ratioslegend->SetMargin(0.35);
    ratioslegend->AddEntry(ratios, "0-20%", "p");
    ratioslegend->AddEntry(ratios2050, "20-50%", "p");
    ratioslegend->AddEntry(ratios50100, "50-100%", "p");
    
  
    TCanvas *testc = new TCanvas("test", "test",50, 50, 600, 600);
    testc->cd();
    ratios2050->GetYaxis()->SetRangeUser(0.0012, 0.0035);
    ratios2050->GetXaxis()->SetLabelSize(0.07);
    ratios2050->Draw("P SAME");
    ratios50100->Draw("P SAME");
    ratios->Draw("P SAME");
    ratioslegend->Draw("SAME");



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


}

