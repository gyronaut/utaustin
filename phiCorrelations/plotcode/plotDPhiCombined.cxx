void plotDPhiCombined(){
    gStyle->SetOptStat(0);
    gStyle->SetOptFit(0);


    TFile* file0_20_Combined = new TFile("~/phiStudies/results_onlineEff/Combined/US_syst_trig_4_8_assoc_2_4_effcorr_hPhi_0_20_Combined.root");    
    //string output = inputfile.substr(0,22);
    //output+= "dphi.pdf";
    TH2D* hPhi2D_0_20_Combined = (TH2D*)(file0_20_Combined->Get("AvgUSsubhPhi2Dpeakavgscale")->Clone("hphi2D_0_20_Combined"));
    TH1D* hPhidphi_0_20_Combined = (TH1D*)hPhi2D_0_20_Combined->ProjectionY("hPhidphi_0_20_Combined", hPhi2D_0_20_Combined->GetXaxis()->FindBin(-1.2), hPhi2D_0_20_Combined->GetXaxis()->FindBin(1.2));
    //hPhidphi_0_20_Combined->Rebin();
    hPhidphi_0_20_Combined->SetLineWidth(2);
    hPhidphi_0_20_Combined->SetLineColor(kViolet);
    hPhidphi_0_20_Combined->SetMarkerColor(kViolet);
    hPhidphi_0_20_Combined->SetMarkerStyle(21);
    hPhidphi_0_20_Combined->SetMarkerSize(2);
    hPhidphi_0_20_Combined->GetXaxis()->SetTitle("#Delta#varphi");
    hPhidphi_0_20_Combined->SetTitle("");
    //hPhidphi_0_20_Combined->Scale(1.0/(hPhidphi_0_20_Combined->Integral()));
    hPhidphi_0_20_Combined->GetYaxis()->SetTitle("Per Trigger Yield");
    hPhidphi_0_20_Combined->GetYaxis()->SetTitleOffset(1.70);
    hPhidphi_0_20_Combined->GetXaxis()->SetTitleSize(0.05);
    hPhidphi_0_20_Combined->GetXaxis()->SetTitleOffset(0.90);

    TF1 *corrFit_0_20_Combined = new TF1("corrFit_0_20_Combined", "gaus(0) + gaus(3) + pol0(6)", -1.4, 4.6);
    corrFit_0_20_Combined->FixParameter(6, 0.25*(hPhidphi_0_20_Combined->GetBinContent(8)+hPhidphi_0_20_Combined->GetBinContent(9)+hPhidphi_0_20_Combined->GetBinContent(16)+hPhidphi_0_20_Combined->GetBinContent(1)));
    //corrFit_0_20_Combined->SetParLimits(6, hPhidphi_0_20_Combined->GetBinContent(8)*0.9, hPhidphi_0_20_Combined->GetBinContent(8)*1.1);
    corrFit_0_20_Combined->SetParameter(0, hPhidphi_0_20_Combined->GetBinContent(hPhidphi_0_20_Combined->GetXaxis()->FindBin(0.0)) - corrFit_0_20_Combined->GetParameter(6));
    corrFit_0_20_Combined->SetParLimits(0, corrFit_0_20_Combined->GetParameter(0)*0.5, corrFit_0_20_Combined->GetParameter(0)*1.5);
    corrFit_0_20_Combined->SetParameter(1, 0.0);
    corrFit_0_20_Combined->SetParLimits(1, -0.5, 0.5);
    corrFit_0_20_Combined->SetParameter(2, 0.5);
    corrFit_0_20_Combined->SetParLimits(2, 0.2, 0.9);
    corrFit_0_20_Combined->SetParameter(3, hPhidphi_0_20_Combined->GetBinContent(hPhidphi_0_20_Combined->GetXaxis()->FindBin(3.14)) - corrFit_0_20_Combined->GetParameter(6));
    corrFit_0_20_Combined->SetParLimits(3, corrFit_0_20_Combined->GetParameter(3)*0.5, corrFit_0_20_Combined->GetParameter(3)*1.5);
    corrFit_0_20_Combined->SetParameter(4, 3.14);
    corrFit_0_20_Combined->SetParLimits(4, 3.0, 3.25);
    corrFit_0_20_Combined->SetParameter(5, 0.5);
    corrFit_0_20_Combined->SetParLimits(5, 0.2, 0.9);

    corrFit_0_20_Combined->SetParNames("peak1", "mean1", "sigma1", "peak2", "mean2", "sigma2", "BG");

    corrFit_0_20_Combined->SetLineColor(kViolet+1);
    corrFit_0_20_Combined->SetLineWidth(2);
    corrFit_0_20_Combined->SetLineStyle(7);



    TFile* file0_20_CENT = new TFile("~/phiStudies/results_onlineEff/CENTwoSDD/US_syst_trig_4_8_assoc_2_4_effcorr_hPhi_0_20_CENTwoSDD.root");    
    //string output = inputfile.substr(0,22);
    //output+= "dphi.pdf";
    TH2D* hPhi2D_0_20_CENT = (TH2D*)(file0_20_CENT->Get("AvgUSsubhPhi2Dpeakavgscale")->Clone("hphi2D_0_20_CENT"));
    TH1D* hPhidphi_0_20_CENT = (TH1D*)hPhi2D_0_20_CENT->ProjectionY("hPhidphi_0_20_CENT", hPhi2D_0_20_CENT->GetXaxis()->FindBin(-1.2), hPhi2D_0_20_CENT->GetXaxis()->FindBin(1.2));
    //hPhidphi_0_20_CENT->Rebin();
    hPhidphi_0_20_CENT->SetLineWidth(2);
    hPhidphi_0_20_CENT->SetLineColor(kBlue);
    hPhidphi_0_20_CENT->SetMarkerColor(kBlue);
    hPhidphi_0_20_CENT->SetMarkerStyle(22);
    hPhidphi_0_20_CENT->SetMarkerSize(2);
    hPhidphi_0_20_CENT->GetXaxis()->SetTitle("#Delta#varphi");
    hPhidphi_0_20_CENT->SetTitle("");
    //hPhidphi_0_20_CENT->Scale(1.0/(hPhidphi_0_20_CENT->Integral()));
    hPhidphi_0_20_CENT->GetYaxis()->SetTitle("Per Trigger Yield");
    hPhidphi_0_20_CENT->GetYaxis()->SetTitleOffset(1.70);
    hPhidphi_0_20_CENT->GetXaxis()->SetTitleSize(0.05);
    hPhidphi_0_20_CENT->GetXaxis()->SetTitleOffset(0.90);

    TF1 *corrFit_0_20_CENT = new TF1("corrFit_0_20_CENT", "gaus(0) + gaus(3) + pol0(6)", -1.4, 4.6);
    corrFit_0_20_CENT->FixParameter(6, 0.25*(hPhidphi_0_20_CENT->GetBinContent(8)+hPhidphi_0_20_CENT->GetBinContent(9)+hPhidphi_0_20_CENT->GetBinContent(16)+hPhidphi_0_20_CENT->GetBinContent(1)));
    //corrFit_0_20_CENT->SetParLimits(6, hPhidphi_0_20_CENT->GetBinContent(8)*0.9, hPhidphi_0_20_CENT->GetBinContent(8)*1.1);
    corrFit_0_20_CENT->SetParameter(0, hPhidphi_0_20_CENT->GetBinContent(hPhidphi_0_20_CENT->GetXaxis()->FindBin(0.0)) - corrFit_0_20_CENT->GetParameter(6));
    corrFit_0_20_CENT->SetParLimits(0, corrFit_0_20_CENT->GetParameter(0)*0.5, corrFit_0_20_CENT->GetParameter(0)*1.5);
    corrFit_0_20_CENT->SetParameter(1, 0.0);
    corrFit_0_20_CENT->SetParLimits(1, -0.5, 0.5);
    corrFit_0_20_CENT->SetParameter(2, 0.5);
    corrFit_0_20_CENT->SetParLimits(2, 0.2, 0.9);
    corrFit_0_20_CENT->SetParameter(3, hPhidphi_0_20_CENT->GetBinContent(hPhidphi_0_20_CENT->GetXaxis()->FindBin(3.14)) - corrFit_0_20_CENT->GetParameter(6));
    corrFit_0_20_CENT->SetParLimits(3, corrFit_0_20_CENT->GetParameter(3)*0.5, corrFit_0_20_CENT->GetParameter(3)*1.5);
    corrFit_0_20_CENT->SetParameter(4, 3.14);
    corrFit_0_20_CENT->SetParLimits(4, 3.0, 3.25);
    corrFit_0_20_CENT->SetParameter(5, 0.5);
    corrFit_0_20_CENT->SetParLimits(5, 0.2, 0.9);

    corrFit_0_20_CENT->SetParNames("peak1", "mean1", "sigma1", "peak2", "mean2", "sigma2", "BG");

    corrFit_0_20_CENT->SetLineColor(kBlue+1);
    corrFit_0_20_CENT->SetLineWidth(2);
    corrFit_0_20_CENT->SetLineStyle(7);

    TH1D* ratio_0_20_CENT = (TH1D*)hPhidphi_0_20_CENT->Clone("ratio_0_20_CENT");
    ratio_0_20_CENT->Divide(hPhidphi_0_20_Combined);
    ratio_0_20_CENT->GetYaxis()->SetTitle("Ratio to Combined");



    TFile* file0_20_FAST = new TFile("~/phiStudies/results_onlineEff/FAST/US_syst_trig_4_8_assoc_2_4_effcorr_hPhi_0_20_FAST.root");    
    TH2D* hPhi2D_0_20_FAST = (TH2D*)(file0_20_FAST->Get("AvgUSsubhPhi2Dpeakavgscale")->Clone("hphi2D_0_20_FAST"));
    TH1D* hPhidphi_0_20_FAST = (TH1D*)hPhi2D_0_20_FAST->ProjectionY("hPhidphi_0_20_FAST", hPhi2D_0_20_FAST->GetXaxis()->FindBin(-1.2), hPhi2D_0_20_FAST->GetXaxis()->FindBin(1.2));
    //hPhidphi_0_20_FAST->Rebin();
    hPhidphi_0_20_FAST->SetLineWidth(2);
    hPhidphi_0_20_FAST->SetLineColor(kRed);
    hPhidphi_0_20_FAST->SetMarkerColor(kRed);
    hPhidphi_0_20_FAST->SetMarkerStyle(23);
    hPhidphi_0_20_FAST->SetMarkerSize(2);
    hPhidphi_0_20_FAST->GetXaxis()->SetTitle("#Delta#varphi");
    hPhidphi_0_20_FAST->SetTitle("");
    //hPhidphi_0_20_FAST->Scale(1.0/(hPhidphi_0_20_FAST->Integral()));
    hPhidphi_0_20_FAST->GetYaxis()->SetTitle("Per Trigger Yield");
    hPhidphi_0_20_FAST->GetYaxis()->SetTitleOffset(1.70);
    hPhidphi_0_20_FAST->GetXaxis()->SetTitleSize(0.05);
    hPhidphi_0_20_FAST->GetXaxis()->SetTitleOffset(0.90);

    TF1 *corrFit_0_20_FAST = new TF1("corrFit_0_20_FAST", "gaus(0) + gaus(3) + pol0(6)", -1.4, 4.6);
    corrFit_0_20_FAST->FixParameter(6, 0.25*(hPhidphi_0_20_FAST->GetBinContent(8)+hPhidphi_0_20_FAST->GetBinContent(9)+hPhidphi_0_20_FAST->GetBinContent(16)+hPhidphi_0_20_FAST->GetBinContent(1)));
    //corrFit_0_20_FAST->SetParLimits(6, hPhidphi_0_20_FAST->GetBinContent(8)*0.9, hPhidphi_0_20_FAST->GetBinContent(8)*1.1);
    corrFit_0_20_FAST->SetParameter(0, hPhidphi_0_20_FAST->GetBinContent(hPhidphi_0_20_FAST->GetXaxis()->FindBin(0.0)) - corrFit_0_20_FAST->GetParameter(6));
    corrFit_0_20_FAST->SetParLimits(0, corrFit_0_20_FAST->GetParameter(0)*0.5, corrFit_0_20_FAST->GetParameter(0)*1.5);
    corrFit_0_20_FAST->SetParameter(1, 0.0);
    corrFit_0_20_FAST->SetParLimits(1, -0.5, 0.5);
    corrFit_0_20_FAST->SetParameter(2, 0.5);
    corrFit_0_20_FAST->SetParLimits(2, 0.2, 0.9);
    corrFit_0_20_FAST->SetParameter(3, hPhidphi_0_20_FAST->GetBinContent(hPhidphi_0_20_FAST->GetXaxis()->FindBin(3.14)) - corrFit_0_20_FAST->GetParameter(6));
    corrFit_0_20_FAST->SetParLimits(3, corrFit_0_20_FAST->GetParameter(3)*0.5, corrFit_0_20_FAST->GetParameter(3)*1.5);
    corrFit_0_20_FAST->SetParameter(4, 3.14);
    corrFit_0_20_FAST->SetParLimits(4, 3.0, 3.25);
    corrFit_0_20_FAST->SetParameter(5, 0.5);
    corrFit_0_20_FAST->SetParLimits(5, 0.2, 0.9);

    corrFit_0_20_FAST->SetParNames("peak1", "mean1", "sigma1", "peak2", "mean2", "sigma2", "BG");

    corrFit_0_20_FAST->SetLineColor(kRed+1);
    corrFit_0_20_FAST->SetLineWidth(2);
    corrFit_0_20_FAST->SetLineStyle(7);

    TH1D* ratio_0_20_FAST = (TH1D*)hPhidphi_0_20_FAST->Clone("ratio_0_20_FAST");
    ratio_0_20_FAST->Divide(hPhidphi_0_20_Combined);
    ratio_0_20_FAST->GetYaxis()->SetTitle("Ratio to Combined Data");
    ratio_0_20_FAST->GetYaxis()->SetRangeUser(0.5, 1.5);
    

    
    TPaveText *text = new TPaveText(0.4815, 0.7056, 0.8658, 0.8551, "NDC");
    //text->AddText("ALICE Work in Progress");
    //text->AddText("p-Pb #sqrt{s_{NN}} = 5 TeV");
    text->AddText("0%-20% Multiplicity");
    text->SetTextSizePixels(20);
    text->SetFillColor(kWhite);

    TPaveText *text2 = new TPaveText(0.6, 0.9, 0.85, 0.85, "NDC");
    text2->AddText("trigger: 4.0 < p_{T}^{h} < 8.0 GeV/c");
    text2->AddText("assoc: 2.0 < p_{T}^{#phi} < 4.0 GeV/c");
    text2->SetTextSizePixels(18);
    text2->SetFillColor(kWhite);

    TLegend  *legend = new TLegend(0.3791, 0.1518, 0.8772, 0.2688);
    legend->SetMargin(0.25);
    legend->AddEntry(hPhidphi_0_20_CENT, "CENT_wo_SDD", "pel");
    legend->AddEntry(hPhidphi_0_20_FAST, "FAST", "pel");
    legend->AddEntry(hPhidphi_0_20_Combined, "Combined", "pel");
   
    TCanvas *c0_20 = new TCanvas("c0_20", "c0_20", 50, 50, 550, 600);
    c0_20->cd();
    c0_20->SetMargin(0.12, 0.05, 0.1, 0.05);

    hPhidphi_0_20_FAST->Fit("corrFit_0_20_FAST", "R");
    hPhidphi_0_20_CENT->Fit("corrFit_0_20_CENT", "R");
    hPhidphi_0_20_Combined->Fit("corrFit_0_20_Combined", "R");

    hPhidphi_0_20_CENT->Draw();
    hPhidphi_0_20_FAST->Draw("SAME");
    hPhidphi_0_20_Combined->Draw("SAME");
    legend->Draw("SAME");
    text->Draw("SAME");
    text2->Draw("SAME");
    //c0_20->Print(output.c_str(),"pdf");
   
    TLegend  *legendRatio = new TLegend(0.3791, 0.1518, 0.8772, 0.2688);
    legendRatio->SetMargin(0.25);
    legendRatio->AddEntry(hPhidphi_0_20_CENT, "CENT/Combined", "pel");
    legendRatio->AddEntry(hPhidphi_0_20_FAST, "FAST/Combined", "pel");

    TCanvas* cratio0_20 = new TCanvas("cratio0_20", "cratio0_20", 50, 50, 550, 600);
    cratio0_20->cd();
    ratio_0_20_FAST->Draw();
    ratio_0_20_CENT->Draw("SAME");
    legendRatio->Draw();
    text->Draw();



   TFile* file20_50_Combined = new TFile("~/phiStudies/results_onlineeff/Combined/US_syst_trig_4_8_assoc_2_4_effcorr_hPhi_20_50_Combined.root");    
    //string output = inputfile.substr(0,22);
    //output+= "dphi.pdf";
    TH2D* hPhi2D_20_50_Combined = (TH2D*)(file20_50_Combined->Get("AvgUSsubhPhi2Dpeak")->Clone("hphi2D_20_50_Combined"));
    TH1D* hPhidphi_20_50_Combined = (TH1D*)hPhi2D_20_50_Combined->ProjectionY("hPhidphi_20_50_Combined", hPhi2D_20_50_Combined->GetXaxis()->FindBin(-1.2), hPhi2D_20_50_Combined->GetXaxis()->FindBin(1.2));
    //hPhidphi_20_50_Combined->Rebin();
    hPhidphi_20_50_Combined->SetLineWidth(2);
    hPhidphi_20_50_Combined->SetLineColor(kViolet);
    hPhidphi_20_50_Combined->SetMarkerColor(kViolet);
    hPhidphi_20_50_Combined->SetMarkerStyle(21);
    hPhidphi_20_50_Combined->SetMarkerSize(2);
    hPhidphi_20_50_Combined->GetXaxis()->SetTitle("#Delta#varphi");
    hPhidphi_20_50_Combined->SetTitle("");
    //hPhidphi_20_50_Combined->Scale(1.0/(hPhidphi_20_50_Combined->Integral()));
    hPhidphi_20_50_Combined->GetYaxis()->SetTitle("Per Trigger Yield");
    hPhidphi_20_50_Combined->GetYaxis()->SetTitleOffset(1.70);
    hPhidphi_20_50_Combined->GetXaxis()->SetTitleSize(0.05);
    hPhidphi_20_50_Combined->GetXaxis()->SetTitleOffset(0.90);

    TF1 *corrFit_20_50_Combined = new TF1("corrFit_20_50_Combined", "gaus(0) + gaus(3) + pol0(6)", -1.4, 4.6);
    corrFit_20_50_Combined->FixParameter(6, 0.25*(hPhidphi_20_50_Combined->GetBinContent(8)+hPhidphi_20_50_Combined->GetBinContent(9)+hPhidphi_20_50_Combined->GetBinContent(16)+hPhidphi_20_50_Combined->GetBinContent(1)));
    //corrFit_20_50_Combined->SetParLimits(6, hPhidphi_20_50_Combined->GetBinContent(8)*0.9, hPhidphi_20_50_Combined->GetBinContent(8)*1.1);
    corrFit_20_50_Combined->SetParameter(0, hPhidphi_20_50_Combined->GetBinContent(hPhidphi_20_50_Combined->GetXaxis()->FindBin(0.0)) - corrFit_20_50_Combined->GetParameter(6));
    corrFit_20_50_Combined->SetParLimits(0, corrFit_20_50_Combined->GetParameter(0)*0.5, corrFit_20_50_Combined->GetParameter(0)*1.5);
    corrFit_20_50_Combined->SetParameter(1, 0.0);
    corrFit_20_50_Combined->SetParLimits(1, -0.5, 0.5);
    corrFit_20_50_Combined->SetParameter(2, 0.5);
    corrFit_20_50_Combined->SetParLimits(2, 0.2, 0.9);
    corrFit_20_50_Combined->SetParameter(3, hPhidphi_20_50_Combined->GetBinContent(hPhidphi_20_50_Combined->GetXaxis()->FindBin(3.14)) - corrFit_20_50_Combined->GetParameter(6));
    corrFit_20_50_Combined->SetParLimits(3, corrFit_20_50_Combined->GetParameter(3)*0.5, corrFit_20_50_Combined->GetParameter(3)*1.5);
    corrFit_20_50_Combined->SetParameter(4, 3.14);
    corrFit_20_50_Combined->SetParLimits(4, 3.0, 3.25);
    corrFit_20_50_Combined->SetParameter(5, 0.5);
    corrFit_20_50_Combined->SetParLimits(5, 0.2, 0.9);

    corrFit_20_50_Combined->SetParNames("peak1", "mean1", "sigma1", "peak2", "mean2", "sigma2", "BG");

    corrFit_20_50_Combined->SetLineColor(kViolet+1);
    corrFit_20_50_Combined->SetLineWidth(2);
    corrFit_20_50_Combined->SetLineStyle(7);


    TFile* file20_50_CENT = new TFile("~/phiStudies/results_onlineeff/CENTwoSDD/US_syst_trig_4_8_assoc_2_4_effcorr_hPhi_20_50_CENTwoSDD.root");    
   
    TH2D* hPhi2D_20_50_CENT = (TH2D*)(file20_50_CENT->Get("AvgUSsubhPhi2Dpeak")->Clone("hphi2D_20_50_CENT"));
    TH1D* hPhidphi_20_50_CENT = (TH1D*)hPhi2D_20_50_CENT->ProjectionY("hPhidphi_20_50_CENT", hPhi2D_20_50_CENT->GetXaxis()->FindBin(-1.2), hPhi2D_20_50_CENT->GetXaxis()->FindBin(1.2));
    //hPhidphi_20_50_CENT->Rebin();
    hPhidphi_20_50_CENT->SetLineWidth(2);
    hPhidphi_20_50_CENT->SetLineColor(kBlue);
    hPhidphi_20_50_CENT->SetMarkerColor(kBlue);
    hPhidphi_20_50_CENT->SetMarkerStyle(22);
    hPhidphi_20_50_CENT->SetMarkerSize(2);
    hPhidphi_20_50_CENT->GetXaxis()->SetTitle("#Delta#varphi");
    hPhidphi_20_50_CENT->SetTitle("");
    //hPhidphi_20_50_CENT->Scale(1.0/(hPhidphi_20_50_CENT->Integral()));
    hPhidphi_20_50_CENT->GetYaxis()->SetTitle("Per Trigger Yield");
    hPhidphi_20_50_CENT->GetYaxis()->SetTitleOffset(1.70);
    hPhidphi_20_50_CENT->GetXaxis()->SetTitleSize(0.05);
    hPhidphi_20_50_CENT->GetXaxis()->SetTitleOffset(0.90);

    TF1 *corrFit_20_50_CENT = new TF1("corrFit_20_50_CENT", "gaus(0) + gaus(3) + pol0(6)", -1.4, 4.6);
    corrFit_20_50_CENT->FixParameter(6, 0.25*(hPhidphi_20_50_CENT->GetBinContent(8)+hPhidphi_20_50_CENT->GetBinContent(9)+hPhidphi_20_50_CENT->GetBinContent(16)+hPhidphi_20_50_CENT->GetBinContent(1)));
    //corrFit_20_50_CENT->SetParLimits(6, hPhidphi_20_50_CENT->GetBinContent(8)*0.9, hPhidphi_20_50_CENT->GetBinContent(8)*1.1);
    corrFit_20_50_CENT->SetParameter(0, hPhidphi_20_50_CENT->GetBinContent(hPhidphi_20_50_CENT->GetXaxis()->FindBin(0.0)) - corrFit_20_50_CENT->GetParameter(6));
    corrFit_20_50_CENT->SetParLimits(0, corrFit_20_50_CENT->GetParameter(0)*0.5, corrFit_20_50_CENT->GetParameter(0)*1.5);
    corrFit_20_50_CENT->SetParameter(1, 0.0);
    corrFit_20_50_CENT->SetParLimits(1, -0.5, 0.5);
    corrFit_20_50_CENT->SetParameter(2, 0.5);
    corrFit_20_50_CENT->SetParLimits(2, 0.2, 0.9);
    corrFit_20_50_CENT->SetParameter(3, hPhidphi_20_50_CENT->GetBinContent(hPhidphi_20_50_CENT->GetXaxis()->FindBin(3.14)) - corrFit_20_50_CENT->GetParameter(6));
    corrFit_20_50_CENT->SetParLimits(3, corrFit_20_50_CENT->GetParameter(3)*0.5, corrFit_20_50_CENT->GetParameter(3)*1.5);
    corrFit_20_50_CENT->SetParameter(4, 3.14);
    corrFit_20_50_CENT->SetParLimits(4, 3.0, 3.25);
    corrFit_20_50_CENT->SetParameter(5, 0.5);
    corrFit_20_50_CENT->SetParLimits(5, 0.2, 0.9);

    corrFit_20_50_CENT->SetParNames("peak1", "mean1", "sigma1", "peak2", "mean2", "sigma2", "BG");

    corrFit_20_50_CENT->SetLineColor(kBlue+1);
    corrFit_20_50_CENT->SetLineWidth(2);
    corrFit_20_50_CENT->SetLineStyle(7);

    TH1D* ratio_20_50_CENT = (TH1D*)hPhidphi_20_50_CENT->Clone("ratio_20_50_CENT");
    ratio_20_50_CENT->Divide(hPhidphi_20_50_Combined);
    ratio_20_50_CENT->GetYaxis()->SetTitle("Ratio to Combined");



    TFile* file20_50_FAST = new TFile("~/phiStudies/results_onlineeff/FAST/US_syst_trig_4_8_assoc_2_4_effcorr_hPhi_20_50_FAST.root");    
    TH2D* hPhi2D_20_50_FAST = (TH2D*)(file20_50_FAST->Get("AvgUSsubhPhi2Dpeak")->Clone("hphi2D_20_50_FAST"));
    TH1D* hPhidphi_20_50_FAST = (TH1D*)hPhi2D_20_50_FAST->ProjectionY("hPhidphi_20_50_FAST", hPhi2D_20_50_FAST->GetXaxis()->FindBin(-1.2), hPhi2D_20_50_FAST->GetXaxis()->FindBin(1.2));
    //hPhidphi_20_50_FAST->Rebin();
    hPhidphi_20_50_FAST->SetLineWidth(2);
    hPhidphi_20_50_FAST->SetLineColor(kRed);
    hPhidphi_20_50_FAST->SetMarkerColor(kRed);
    hPhidphi_20_50_FAST->SetMarkerStyle(23);
    hPhidphi_20_50_FAST->SetMarkerSize(2);
    hPhidphi_20_50_FAST->GetXaxis()->SetTitle("#Delta#varphi");
    hPhidphi_20_50_FAST->SetTitle("");
    //hPhidphi_20_50_FAST->Scale(1.0/(hPhidphi_20_50_FAST->Integral()));
    hPhidphi_20_50_FAST->GetYaxis()->SetTitle("Per Trigger Yield");
    hPhidphi_20_50_FAST->GetYaxis()->SetTitleOffset(1.70);
    hPhidphi_20_50_FAST->GetXaxis()->SetTitleSize(0.05);
    hPhidphi_20_50_FAST->GetXaxis()->SetTitleOffset(0.90);

    TF1 *corrFit_20_50_FAST = new TF1("corrFit_20_50_FAST", "gaus(0) + gaus(3) + pol0(6)", -1.4, 4.6);
    corrFit_20_50_FAST->FixParameter(6, 0.25*(hPhidphi_20_50_FAST->GetBinContent(8)+hPhidphi_20_50_FAST->GetBinContent(9)+hPhidphi_20_50_FAST->GetBinContent(16)+hPhidphi_20_50_FAST->GetBinContent(1)));
    //corrFit_20_50_FAST->SetParLimits(6, hPhidphi_20_50_FAST->GetBinContent(8)*0.9, hPhidphi_20_50_FAST->GetBinContent(8)*1.1);
    corrFit_20_50_FAST->SetParameter(0, hPhidphi_20_50_FAST->GetBinContent(hPhidphi_20_50_FAST->GetXaxis()->FindBin(0.0)) - corrFit_20_50_FAST->GetParameter(6));
    corrFit_20_50_FAST->SetParLimits(0, corrFit_20_50_FAST->GetParameter(0)*0.5, corrFit_20_50_FAST->GetParameter(0)*1.5);
    corrFit_20_50_FAST->SetParameter(1, 0.0);
    corrFit_20_50_FAST->SetParLimits(1, -0.5, 0.5);
    corrFit_20_50_FAST->SetParameter(2, 0.5);
    corrFit_20_50_FAST->SetParLimits(2, 0.2, 0.9);
    corrFit_20_50_FAST->SetParameter(3, hPhidphi_20_50_FAST->GetBinContent(hPhidphi_20_50_FAST->GetXaxis()->FindBin(3.14)) - corrFit_20_50_FAST->GetParameter(6));
    corrFit_20_50_FAST->SetParLimits(3, corrFit_20_50_FAST->GetParameter(3)*0.5, corrFit_20_50_FAST->GetParameter(3)*1.5);
    corrFit_20_50_FAST->SetParameter(4, 3.14);
    corrFit_20_50_FAST->SetParLimits(4, 3.0, 3.25);
    corrFit_20_50_FAST->SetParameter(5, 0.5);
    corrFit_20_50_FAST->SetParLimits(5, 0.2, 0.9);

    corrFit_20_50_FAST->SetParNames("peak1", "mean1", "sigma1", "peak2", "mean2", "sigma2", "BG");

    corrFit_20_50_FAST->SetLineColor(kRed+1);
    corrFit_20_50_FAST->SetLineWidth(2);
    corrFit_20_50_FAST->SetLineStyle(7);

    TH1D* ratio_20_50_FAST = (TH1D*)hPhidphi_20_50_FAST->Clone("ratio_20_50_FAST");
    ratio_20_50_FAST->Divide(hPhidphi_20_50_Combined);
    ratio_20_50_FAST->GetYaxis()->SetTitle("Ratio to Combined Data");
    ratio_20_50_FAST->GetYaxis()->SetRangeUser(0.5, 1.5);
 
    TPaveText *text_20_50 = new TPaveText(0.4815, 0.7056, 0.8658, 0.8551, "NDC");
    text_20_50->AddText("20%-50% Multiplicity");
    text_20_50->SetTextSizePixels(20);
    text_20_50->SetFillColor(kWhite);
   
   
    TCanvas *c20_50 = new TCanvas("c20_50", "c20_50", 50, 50, 550, 600);
    c20_50->cd();
    c20_50->SetMargin(0.12, 0.05, 0.1, 0.05);

    hPhidphi_20_50_FAST->Fit("corrFit_20_50_FAST", "R");
    hPhidphi_20_50_CENT->Fit("corrFit_20_50_CENT", "R");
    hPhidphi_20_50_Combined->Fit("corrFit_20_50_Combined", "R");

    hPhidphi_20_50_CENT->Draw();
    hPhidphi_20_50_FAST->Draw("SAME");
    hPhidphi_20_50_Combined->Draw("SAME");
    legend->Draw("SAME");
    text_20_50->Draw("SAME");
    text2->Draw("SAME");
   
    TCanvas* cratio20_50 = new TCanvas("cratio20_50", "cratio20_50", 50, 50, 550, 600);
    cratio20_50->cd();
    ratio_20_50_FAST->Draw();
    ratio_20_50_CENT->Draw("SAME");
    legendRatio->Draw();
    text_20_50->Draw();


   TFile* file50_80_Combined = new TFile("~/phiStudies/results_onlineeff/Combined/US_syst_trig_4_8_assoc_2_4_effcorr_hPhi_50_80_Combined.root");    
    //string output = inputfile.substr(0,22);
    //output+= "dphi.pdf";
    TH2D* hPhi2D_50_80_Combined = (TH2D*)(file50_80_Combined->Get("AvgUSsubhPhi2Dpeak")->Clone("hphi2D_50_80_Combined"));
    TH1D* hPhidphi_50_80_Combined = (TH1D*)hPhi2D_50_80_Combined->ProjectionY("hPhidphi_50_80_Combined", hPhi2D_50_80_Combined->GetXaxis()->FindBin(-1.2), hPhi2D_50_80_Combined->GetXaxis()->FindBin(1.2));
    //hPhidphi_50_80_Combined->Rebin();
    hPhidphi_50_80_Combined->SetLineWidth(2);
    hPhidphi_50_80_Combined->SetLineColor(kViolet);
    hPhidphi_50_80_Combined->SetMarkerColor(kViolet);
    hPhidphi_50_80_Combined->SetMarkerStyle(21);
    hPhidphi_50_80_Combined->SetMarkerSize(2);
    hPhidphi_50_80_Combined->GetXaxis()->SetTitle("#Delta#varphi");
    hPhidphi_50_80_Combined->SetTitle("");
    //hPhidphi_50_80_Combined->Scale(1.0/(hPhidphi_50_80_Combined->Integral()));
    hPhidphi_50_80_Combined->GetYaxis()->SetTitle("Per Trigger Yield");
    hPhidphi_50_80_Combined->GetYaxis()->SetTitleOffset(1.70);
    hPhidphi_50_80_Combined->GetXaxis()->SetTitleSize(0.05);
    hPhidphi_50_80_Combined->GetXaxis()->SetTitleOffset(0.90);

    TF1 *corrFit_50_80_Combined = new TF1("corrFit_50_80_Combined", "gaus(0) + gaus(3) + pol0(6)", -1.4, 4.6);
    corrFit_50_80_Combined->FixParameter(6, 0.25*(hPhidphi_50_80_Combined->GetBinContent(8)+hPhidphi_50_80_Combined->GetBinContent(9)+hPhidphi_50_80_Combined->GetBinContent(16)+hPhidphi_50_80_Combined->GetBinContent(1)));
    //corrFit_50_80_Combined->SetParLimits(6, hPhidphi_50_80_Combined->GetBinContent(8)*0.9, hPhidphi_50_80_Combined->GetBinContent(8)*1.1);
    corrFit_50_80_Combined->SetParameter(0, hPhidphi_50_80_Combined->GetBinContent(hPhidphi_50_80_Combined->GetXaxis()->FindBin(0.0)) - corrFit_50_80_Combined->GetParameter(6));
    corrFit_50_80_Combined->SetParLimits(0, corrFit_50_80_Combined->GetParameter(0)*0.5, corrFit_50_80_Combined->GetParameter(0)*1.5);
    corrFit_50_80_Combined->SetParameter(1, 0.0);
    corrFit_50_80_Combined->SetParLimits(1, -0.5, 0.5);
    corrFit_50_80_Combined->SetParameter(2, 0.5);
    corrFit_50_80_Combined->SetParLimits(2, 0.2, 0.9);
    corrFit_50_80_Combined->SetParameter(3, hPhidphi_50_80_Combined->GetBinContent(hPhidphi_50_80_Combined->GetXaxis()->FindBin(3.14)) - corrFit_50_80_Combined->GetParameter(6));
    corrFit_50_80_Combined->SetParLimits(3, corrFit_50_80_Combined->GetParameter(3)*0.5, corrFit_50_80_Combined->GetParameter(3)*1.5);
    corrFit_50_80_Combined->SetParameter(4, 3.14);
    corrFit_50_80_Combined->SetParLimits(4, 3.0, 3.25);
    corrFit_50_80_Combined->SetParameter(5, 0.5);
    corrFit_50_80_Combined->SetParLimits(5, 0.2, 0.9);

    corrFit_50_80_Combined->SetParNames("peak1", "mean1", "sigma1", "peak2", "mean2", "sigma2", "BG");

    corrFit_50_80_Combined->SetLineColor(kViolet+1);
    corrFit_50_80_Combined->SetLineWidth(2);
    corrFit_50_80_Combined->SetLineStyle(7); 


    TFile* file50_80_CENT = new TFile("~/phiStudies/results_onlineeff/CENTwoSDD/US_syst_trig_4_8_assoc_2_4_effcorr_hPhi_50_80_CENTwoSDD.root");    
    //string output = inputfile.substr(0,22);
    //output+= "dphi.pdf";
    TH2D* hPhi2D_50_80_CENT = (TH2D*)(file50_80_CENT->Get("AvgUSsubhPhi2Dpeak")->Clone("hphi2D_50_80_CENT"));
    TH1D* hPhidphi_50_80_CENT = (TH1D*)hPhi2D_50_80_CENT->ProjectionY("hPhidphi_50_80_CENT", hPhi2D_50_80_CENT->GetXaxis()->FindBin(-1.2), hPhi2D_50_80_CENT->GetXaxis()->FindBin(1.2));
    //hPhidphi_50_80_CENT->Rebin();
    hPhidphi_50_80_CENT->SetLineWidth(2);
    hPhidphi_50_80_CENT->SetLineColor(kBlue);
    hPhidphi_50_80_CENT->SetMarkerColor(kBlue);
    hPhidphi_50_80_CENT->SetMarkerStyle(22);
    hPhidphi_50_80_CENT->SetMarkerSize(2);
    hPhidphi_50_80_CENT->GetXaxis()->SetTitle("#Delta#varphi");
    hPhidphi_50_80_CENT->SetTitle("");
    //hPhidphi_50_80_CENT->Scale(1.0/(hPhidphi_50_80_CENT->Integral()));
    hPhidphi_50_80_CENT->GetYaxis()->SetTitle("Per Trigger Yield");
    hPhidphi_50_80_CENT->GetYaxis()->SetTitleOffset(1.70);
    hPhidphi_50_80_CENT->GetXaxis()->SetTitleSize(0.05);
    hPhidphi_50_80_CENT->GetXaxis()->SetTitleOffset(0.90);

    TF1 *corrFit_50_80_CENT = new TF1("corrFit_50_80_CENT", "gaus(0) + gaus(3) + pol0(6)", -1.4, 4.6);
    corrFit_50_80_CENT->FixParameter(6, 0.25*(hPhidphi_50_80_CENT->GetBinContent(8)+hPhidphi_50_80_CENT->GetBinContent(9)+hPhidphi_50_80_CENT->GetBinContent(16)+hPhidphi_50_80_CENT->GetBinContent(1)));
    //corrFit_50_80_CENT->SetParLimits(6, hPhidphi_50_80_CENT->GetBinContent(8)*0.9, hPhidphi_50_80_CENT->GetBinContent(8)*1.1);
    corrFit_50_80_CENT->SetParameter(0, hPhidphi_50_80_CENT->GetBinContent(hPhidphi_50_80_CENT->GetXaxis()->FindBin(0.0)) - corrFit_50_80_CENT->GetParameter(6));
    corrFit_50_80_CENT->SetParLimits(0, corrFit_50_80_CENT->GetParameter(0)*0.5, corrFit_50_80_CENT->GetParameter(0)*1.5);
    corrFit_50_80_CENT->SetParameter(1, 0.0);
    corrFit_50_80_CENT->SetParLimits(1, -0.5, 0.5);
    corrFit_50_80_CENT->SetParameter(2, 0.5);
    corrFit_50_80_CENT->SetParLimits(2, 0.2, 0.9);
    corrFit_50_80_CENT->SetParameter(3, hPhidphi_50_80_CENT->GetBinContent(hPhidphi_50_80_CENT->GetXaxis()->FindBin(3.14)) - corrFit_50_80_CENT->GetParameter(6));
    corrFit_50_80_CENT->SetParLimits(3, corrFit_50_80_CENT->GetParameter(3)*0.5, corrFit_50_80_CENT->GetParameter(3)*1.5);
    corrFit_50_80_CENT->SetParameter(4, 3.14);
    corrFit_50_80_CENT->SetParLimits(4, 3.0, 3.25);
    corrFit_50_80_CENT->SetParameter(5, 0.5);
    corrFit_50_80_CENT->SetParLimits(5, 0.2, 0.9);

    corrFit_50_80_CENT->SetParNames("peak1", "mean1", "sigma1", "peak2", "mean2", "sigma2", "BG");

    corrFit_50_80_CENT->SetLineColor(kBlue+1);
    corrFit_50_80_CENT->SetLineWidth(2);
    corrFit_50_80_CENT->SetLineStyle(7);

    TH1D* ratio_50_80_CENT = (TH1D*)hPhidphi_50_80_CENT->Clone("ratio_50_80_CENT");
    ratio_50_80_CENT->Divide(hPhidphi_50_80_Combined);
    ratio_50_80_CENT->GetYaxis()->SetTitle("Ratio to Combined");


    TFile* file50_80_FAST = new TFile("~/phiStudies/results_onlineeff/FAST/US_syst_trig_4_8_assoc_2_4_effcorr_hPhi_50_80_FAST.root");    
    TH2D* hPhi2D_50_80_FAST = (TH2D*)(file50_80_FAST->Get("AvgUSsubhPhi2Dpeak")->Clone("hphi2D_50_80_FAST"));
    TH1D* hPhidphi_50_80_FAST = (TH1D*)hPhi2D_50_80_FAST->ProjectionY("hPhidphi_50_80_FAST", hPhi2D_50_80_FAST->GetXaxis()->FindBin(-1.2), hPhi2D_50_80_FAST->GetXaxis()->FindBin(1.2));
    //hPhidphi_50_80_FAST->Rebin();
    hPhidphi_50_80_FAST->SetLineWidth(2);
    hPhidphi_50_80_FAST->SetLineColor(kRed);
    hPhidphi_50_80_FAST->SetMarkerColor(kRed);
    hPhidphi_50_80_FAST->SetMarkerStyle(23);
    hPhidphi_50_80_FAST->SetMarkerSize(2);
    hPhidphi_50_80_FAST->GetXaxis()->SetTitle("#Delta#varphi");
    hPhidphi_50_80_FAST->SetTitle("");
    //hPhidphi_50_80_FAST->Scale(1.0/(hPhidphi_50_80_FAST->Integral()));
    hPhidphi_50_80_FAST->GetYaxis()->SetTitle("Per Trigger Yield");
    hPhidphi_50_80_FAST->GetYaxis()->SetTitleOffset(1.70);
    hPhidphi_50_80_FAST->GetXaxis()->SetTitleSize(0.05);
    hPhidphi_50_80_FAST->GetXaxis()->SetTitleOffset(0.90);

    TF1 *corrFit_50_80_FAST = new TF1("corrFit_50_80_FAST", "gaus(0) + gaus(3) + pol0(6)", -1.4, 4.6);
    corrFit_50_80_FAST->FixParameter(6, 0.25*(hPhidphi_50_80_FAST->GetBinContent(8)+hPhidphi_50_80_FAST->GetBinContent(9)+hPhidphi_50_80_FAST->GetBinContent(16)+hPhidphi_50_80_FAST->GetBinContent(1)));
    //corrFit_50_80_FAST->SetParLimits(6, hPhidphi_50_80_FAST->GetBinContent(8)*0.9, hPhidphi_50_80_FAST->GetBinContent(8)*1.1);
    corrFit_50_80_FAST->SetParameter(0, hPhidphi_50_80_FAST->GetBinContent(hPhidphi_50_80_FAST->GetXaxis()->FindBin(0.0)) - corrFit_50_80_FAST->GetParameter(6));
    corrFit_50_80_FAST->SetParLimits(0, corrFit_50_80_FAST->GetParameter(0)*0.5, corrFit_50_80_FAST->GetParameter(0)*1.5);
    corrFit_50_80_FAST->SetParameter(1, 0.0);
    corrFit_50_80_FAST->SetParLimits(1, -0.5, 0.5);
    corrFit_50_80_FAST->SetParameter(2, 0.5);
    corrFit_50_80_FAST->SetParLimits(2, 0.2, 0.9);
    corrFit_50_80_FAST->SetParameter(3, hPhidphi_50_80_FAST->GetBinContent(hPhidphi_50_80_FAST->GetXaxis()->FindBin(3.14)) - corrFit_50_80_FAST->GetParameter(6));
    corrFit_50_80_FAST->SetParLimits(3, corrFit_50_80_FAST->GetParameter(3)*0.5, corrFit_50_80_FAST->GetParameter(3)*1.5);
    corrFit_50_80_FAST->SetParameter(4, 3.14);
    corrFit_50_80_FAST->SetParLimits(4, 3.0, 3.25);
    corrFit_50_80_FAST->SetParameter(5, 0.5);
    corrFit_50_80_FAST->SetParLimits(5, 0.2, 0.9);

    corrFit_50_80_FAST->SetParNames("peak1", "mean1", "sigma1", "peak2", "mean2", "sigma2", "BG");

    corrFit_50_80_FAST->SetLineColor(kRed+1);
    corrFit_50_80_FAST->SetLineWidth(2);
    corrFit_50_80_FAST->SetLineStyle(7);

    TH1D* ratio_50_80_FAST = (TH1D*)hPhidphi_50_80_FAST->Clone("ratio_50_80_FAST");
    ratio_50_80_FAST->Divide(hPhidphi_50_80_Combined);
    ratio_50_80_FAST->GetYaxis()->SetTitle("Ratio to Combined Data");
    ratio_50_80_FAST->GetYaxis()->SetRangeUser(0.5, 1.5);
    
    TPaveText *text_50_80 = new TPaveText(0.4815, 0.7056, 0.8658, 0.8551, "NDC");
    text_50_80->AddText("50%-80% Multiplicity");
    text_50_80->SetTextSizePixels(20);
    text_50_80->SetFillColor(kWhite);
   
       
    TCanvas *c50_80 = new TCanvas("c50_80", "c50_80", 50, 50, 550, 600);
    c50_80->cd();
    c50_80->SetMargin(0.12, 0.05, 0.1, 0.05);

    hPhidphi_50_80_FAST->Fit("corrFit_50_80_FAST", "R");
    hPhidphi_50_80_CENT->Fit("corrFit_50_80_CENT", "R");
    hPhidphi_50_80_Combined->Fit("corrFit_50_80_Combined", "R");

    hPhidphi_50_80_CENT->Draw();
    hPhidphi_50_80_FAST->Draw("SAME");
    hPhidphi_50_80_Combined->Draw("SAME");
    legend->Draw("SAME");
    text_50_80->Draw("SAME");
    text2->Draw("SAME");
    //c50_80->Print(output.c_str(),"pdf"); 
    
    TCanvas* cratio50_80 = new TCanvas("cratio50_80", "cratio50_80", 50, 50, 550, 600);
    cratio50_80->cd();
    ratio_50_80_FAST->Draw();
    ratio_50_80_CENT->Draw("SAME");
    legendRatio->Draw();
    text_50_80->Draw();


 }
