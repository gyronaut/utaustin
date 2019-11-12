void plotHHDphi(){
    gStyle->SetOptStat(0);
    gStyle->SetOptFit(0);

    Float_t epsilon = 0.0001;

    TFile* filehh_0_20 = new TFile("~/phiStudies/results_onlineEff/Combined/trig_2_4_assoc_1_2_effcorr_hh_0_20.root"); 
    //string output = inputfile.substr(0,22);
    //output+= "dphi.pdf";
    TH2D* hh2D_0_20 = (TH2D*)(filehh_0_20->Get("hh2D")->Clone("hh2D_0_20"));
    TH1D* hhdphi_0_20 = (TH1D*)hh2D_0_20->ProjectionY("hhdphi_0_20", hh2D_0_20->GetXaxis()->FindBin(-1.6 + epsilon), hh2D_0_20->GetXaxis()->FindBin(1.6 - epsilon));
    hhdphi_0_20->SetLineWidth(4);
    hhdphi_0_20->SetLineColor(kBlue+2);
    hhdphi_0_20->SetMarkerColor(kBlue+2);
    hhdphi_0_20->SetMarkerStyle(21);
    hhdphi_0_20->SetMarkerSize(2);
    hhdphi_0_20->GetXaxis()->SetTitle("#Delta#varphi");
    hhdphi_0_20->SetTitle("");
    //hhdphi_0_20->Scale(1.0/(hhdphi_0_20->Integral()));
    hhdphi_0_20->GetYaxis()->SetTitle("1/N_{Trig} dN/d#Delta#varphi per #Delta#eta - const.");
    hhdphi_0_20->GetYaxis()->SetTitleOffset(1.60);
    hhdphi_0_20->GetXaxis()->SetTitleSize(0.05);
    hhdphi_0_20->GetXaxis()->SetTitleOffset(0.90);

    TF1 *corrFit_0_20 = new TF1("corrFit_0_20", "gaus(0) + gaus(3) + pol0(6)", -1.4, 4.6);
    corrFit_0_20->FixParameter(6, 0.25*(hhdphi_0_20->GetBinContent(8)+hhdphi_0_20->GetBinContent(9)+hhdphi_0_20->GetBinContent(16)+hhdphi_0_20->GetBinContent(1)));
    //corrFit_0_20->SetParLimits(6, hhdphi_0_20->GetBinContent(8)*0.9, hhdphi_0_20->GetBinContent(8)*1.1);
    corrFit_0_20->SetParameter(0, hhdphi_0_20->GetBinContent(hhdphi_0_20->GetXaxis()->FindBin(0.0)) - corrFit_0_20->GetParameter(6));
    corrFit_0_20->SetParLimits(0, corrFit_0_20->GetParameter(0)*0.5, corrFit_0_20->GetParameter(0)*2.0);
    corrFit_0_20->SetParameter(1, 0.0);
    corrFit_0_20->SetParLimits(1, -0.5, 0.5);
    corrFit_0_20->SetParameter(2, 0.5);
    corrFit_0_20->SetParLimits(2, 0.2, 0.9);
    corrFit_0_20->SetParameter(3, hhdphi_0_20->GetBinContent(hhdphi_0_20->GetXaxis()->FindBin(3.14)) - corrFit_0_20->GetParameter(6));
    corrFit_0_20->SetParLimits(3, corrFit_0_20->GetParameter(3)*0.5, corrFit_0_20->GetParameter(3)*1.5);
    corrFit_0_20->SetParameter(4, 3.14);
    corrFit_0_20->SetParLimits(4, 3.0, 3.25);
    corrFit_0_20->SetParameter(5, 0.5);
    corrFit_0_20->SetParLimits(5, 0.2, 0.9);

    corrFit_0_20->SetParNames("peak1", "mean1", "sigma1", "peak2", "mean2", "sigma2", "BG");

    corrFit_0_20->SetLineColor(kBlue);
    corrFit_0_20->SetLineWidth(6);
    corrFit_0_20->SetLineStyle(7);

    hhdphi_0_20->Fit("corrFit_0_20", "R0");
    
    TF1 *gaus1 = new TF1("gaus1", "gaus(0)", -1.4, 4.6);
    TF1 *gaus2 = new TF1("gaus2", "gaus(0)", -1.4, 4.6);
    TF1 *bg = new TF1("bg", "pol0(0)", -1.4, 4.6);
    gaus1->SetParameter(0, corrFit_0_20->GetParameter(0));
    gaus1->SetParameter(1, corrFit_0_20->GetParameter(1));
    gaus1->SetParameter(2, corrFit_0_20->GetParameter(2));
    gaus2->SetParameter(0, corrFit_0_20->GetParameter(3));
    gaus2->SetParameter(1, corrFit_0_20->GetParameter(4));
    gaus2->SetParameter(2, corrFit_0_20->GetParameter(5));
    
    //bg->SetParameter(0, corrFit->GetParameter(6));    
    //bg->FixParameter(0, 0.25*(hhdphi_0_20->GetBinContent(1) + hhdphi_0_20->GetBinContent(8) +hhdphi_0_20->GetBinContent(9) +hhdphi_0_20->GetBinContent(16)));
    bg->FixParameter(0, 0.5*(hhdphi_0_20->GetBinContent(1) + hhdphi_0_20->GetBinContent(8)));

    TH1D* hhdphi_0_20_corr = (TH1D*)hhdphi_0_20->Clone("hhdphi_0_20_corr");
    hhdphi_0_20_corr->Add(bg, -1.0);
    hhdphi_0_20_corr->Scale(1.0/3.2); //correct for larger delta eta range "per Delta eta"
    for(int i = 1; i <= hhdphi_0_20_corr->GetXaxis()->GetNbins(); i++){
        hhdphi_0_20_corr->SetBinContent(i, hhdphi_0_20_corr->GetBinContent(i)/hhdphi_0_20_corr->GetXaxis()->GetBinWidth(i));
    }

   
    TFile* filehh_20_50 = new TFile("~/phiStudies/results_onlineEff/Combined/trig_2_4_assoc_1_2_effcorr_hh_20_50.root"); 
    //string output = inputfile.substr(0,22);
    //output+= "dphi.pdf";
    TH2D* hh2D_20_50 = (TH2D*)(filehh_20_50->Get("hh2D")->Clone("hh2D_20_50"));
    TH1D* hhdphi_20_50 = (TH1D*)hh2D_20_50->ProjectionY("hhdphi_20_50", hh2D_20_50->GetXaxis()->FindBin(-1.6 + epsilon), hh2D_20_50->GetXaxis()->FindBin(1.6 - epsilon));
    hhdphi_20_50->SetLineWidth(4);
    hhdphi_20_50->SetLineColor(kBlue+2);
    hhdphi_20_50->SetMarkerColor(kBlue+2);
    hhdphi_20_50->SetMarkerStyle(21);
    hhdphi_20_50->SetMarkerSize(2);
    hhdphi_20_50->GetXaxis()->SetTitle("#Delta#varphi");
    hhdphi_20_50->SetTitle("");
    //hhdphi_20_50->Scale(1.0/(hhdphi_20_50->Integral()));
    hhdphi_20_50->GetYaxis()->SetTitle("1/N_{Trig} dN/d#Delta#varphi per #Delta#eta - const.");
    hhdphi_20_50->GetYaxis()->SetTitleOffset(1.60);
    hhdphi_20_50->GetXaxis()->SetTitleSize(0.05);
    hhdphi_20_50->GetXaxis()->SetTitleOffset(0.90);

    TF1 *corrFit_20_50 = new TF1("corrFit_20_50", "gaus(0) + gaus(3) + pol0(6)", -1.4, 4.6);
    corrFit_20_50->FixParameter(6, 0.25*(hhdphi_20_50->GetBinContent(8)+hhdphi_20_50->GetBinContent(9)+hhdphi_20_50->GetBinContent(16)+hhdphi_20_50->GetBinContent(1)));
    //corrFit_20_50->SetParLimits(6, hhdphi_20_50->GetBinContent(8)*0.9, hhdphi_20_50->GetBinContent(8)*1.1);
    corrFit_20_50->SetParameter(0, hhdphi_20_50->GetBinContent(hhdphi_20_50->GetXaxis()->FindBin(0.0)) - corrFit_20_50->GetParameter(6));
    corrFit_20_50->SetParLimits(0, corrFit_20_50->GetParameter(0)*0.5, corrFit_20_50->GetParameter(0)*2.0);
    corrFit_20_50->SetParameter(1, 0.0);
    corrFit_20_50->SetParLimits(1, -0.5, 0.5);
    corrFit_20_50->SetParameter(2, 0.5);
    corrFit_20_50->SetParLimits(2, 0.2, 0.9);
    corrFit_20_50->SetParameter(3, hhdphi_20_50->GetBinContent(hhdphi_20_50->GetXaxis()->FindBin(3.14)) - corrFit_20_50->GetParameter(6));
    corrFit_20_50->SetParLimits(3, corrFit_20_50->GetParameter(3)*0.5, corrFit_20_50->GetParameter(3)*1.5);
    corrFit_20_50->SetParameter(4, 3.14);
    corrFit_20_50->SetParLimits(4, 3.0, 3.25);
    corrFit_20_50->SetParameter(5, 0.5);
    corrFit_20_50->SetParLimits(5, 0.2, 0.9);

    corrFit_20_50->SetParNames("peak1", "mean1", "sigma1", "peak2", "mean2", "sigma2", "BG");

    corrFit_20_50->SetLineColor(kBlue);
    corrFit_20_50->SetLineWidth(6);
    corrFit_20_50->SetLineStyle(7);

    hhdphi_20_50->Fit("corrFit_20_50", "R0");
    
    TF1 *gaus1_20_50 = new TF1("gaus1_20_50", "gaus(0)", -1.4, 4.6);
    TF1 *gaus2_20_50 = new TF1("gaus2_20_50", "gaus(0)", -1.4, 4.6);
    TF1 *bg_20_50 = new TF1("bg_20_50", "pol0(0)", -1.4, 4.6);
    gaus1_20_50->SetParameter(0, corrFit_20_50->GetParameter(0));
    gaus1_20_50->SetParameter(1, corrFit_20_50->GetParameter(1));
    gaus1_20_50->SetParameter(2, corrFit_20_50->GetParameter(2));
    gaus2_20_50->SetParameter(0, corrFit_20_50->GetParameter(3));
    gaus2_20_50->SetParameter(1, corrFit_20_50->GetParameter(4));
    gaus2_20_50->SetParameter(2, corrFit_20_50->GetParameter(5));
    
    //bg->SetParameter(0, corrFit->GetParameter(6));    
    //bg->FixParameter(0, 0.25*(hhdphi_20_50->GetBinContent(1) + hhdphi_20_50->GetBinContent(8) +hhdphi_20_50->GetBinContent(9) +hhdphi_20_50->GetBinContent(16)));
    bg_20_50->FixParameter(0, 0.5*(hhdphi_20_50->GetBinContent(1) + hhdphi_20_50->GetBinContent(8)));

    TH1D* hhdphi_20_50_corr = (TH1D*)hhdphi_20_50->Clone("hhdphi_20_50_corr");
    hhdphi_20_50_corr->Add(bg_20_50, -1.0);
    hhdphi_20_50_corr->Scale(1.0/3.2); //correct for larger delta eta range "per Delta eta"
    for(int i = 1; i <= hhdphi_20_50_corr->GetXaxis()->GetNbins(); i++){
        hhdphi_20_50_corr->SetBinContent(i, hhdphi_20_50_corr->GetBinContent(i)/hhdphi_20_50_corr->GetXaxis()->GetBinWidth(i));
    }

      TFile* filehh_50_80 = new TFile("~/phiStudies/results_onlineEff/Combined/trig_2_4_assoc_1_2_effcorr_hh_50_80.root"); 
    //string output = inputfile.substr(0,22);
    //output+= "dphi.pdf";
    TH2D* hh2D_50_80 = (TH2D*)(filehh_50_80->Get("hh2D")->Clone("hh2D_50_80"));
    TH1D* hhdphi_50_80 = (TH1D*)hh2D_50_80->ProjectionY("hhdphi_50_80", hh2D_50_80->GetXaxis()->FindBin(-1.6 + epsilon), hh2D_50_80->GetXaxis()->FindBin(1.6 - epsilon));
    hhdphi_50_80->SetLineWidth(4);
    hhdphi_50_80->SetLineColor(kBlue+2);
    hhdphi_50_80->SetMarkerColor(kBlue+2);
    hhdphi_50_80->SetMarkerStyle(21);
    hhdphi_50_80->SetMarkerSize(2);
    hhdphi_50_80->GetXaxis()->SetTitle("#Delta#varphi");
    hhdphi_50_80->SetTitle("");
    //hhdphi_50_80->Scale(1.0/(hhdphi_50_80->Integral()));
    hhdphi_50_80->GetYaxis()->SetTitle("1/N_{Trig} dN/d#Delta#varphi per #Delta#eta - const.");
    hhdphi_50_80->GetYaxis()->SetTitleOffset(1.60);
    hhdphi_50_80->GetXaxis()->SetTitleSize(0.05);
    hhdphi_50_80->GetXaxis()->SetTitleOffset(0.90);

    TF1 *corrFit_50_80 = new TF1("corrFit_50_80", "gaus(0) + gaus(3) + pol0(6)", -1.4, 4.6);
    corrFit_50_80->FixParameter(6, 0.25*(hhdphi_50_80->GetBinContent(8)+hhdphi_50_80->GetBinContent(9)+hhdphi_50_80->GetBinContent(16)+hhdphi_50_80->GetBinContent(1)));
    //corrFit_50_80->SetParLimits(6, hhdphi_50_80->GetBinContent(8)*0.9, hhdphi_50_80->GetBinContent(8)*1.1);
    corrFit_50_80->SetParameter(0, hhdphi_50_80->GetBinContent(hhdphi_50_80->GetXaxis()->FindBin(0.0)) - corrFit_50_80->GetParameter(6));
    corrFit_50_80->SetParLimits(0, corrFit_50_80->GetParameter(0)*0.5, corrFit_50_80->GetParameter(0)*2.0);
    corrFit_50_80->SetParameter(1, 0.0);
    corrFit_50_80->SetParLimits(1, -0.5, 0.5);
    corrFit_50_80->SetParameter(2, 0.5);
    corrFit_50_80->SetParLimits(2, 0.2, 0.9);
    corrFit_50_80->SetParameter(3, hhdphi_50_80->GetBinContent(hhdphi_50_80->GetXaxis()->FindBin(3.14)) - corrFit_50_80->GetParameter(6));
    corrFit_50_80->SetParLimits(3, corrFit_50_80->GetParameter(3)*0.5, corrFit_50_80->GetParameter(3)*1.5);
    corrFit_50_80->SetParameter(4, 3.14);
    corrFit_50_80->SetParLimits(4, 3.0, 3.25);
    corrFit_50_80->SetParameter(5, 0.5);
    corrFit_50_80->SetParLimits(5, 0.2, 0.9);

    corrFit_50_80->SetParNames("peak1", "mean1", "sigma1", "peak2", "mean2", "sigma2", "BG");

    corrFit_50_80->SetLineColor(kBlue);
    corrFit_50_80->SetLineWidth(6);
    corrFit_50_80->SetLineStyle(7);

    hhdphi_50_80->Fit("corrFit_50_80", "R0");
    
    TF1 *gaus1_50_80 = new TF1("gaus1_50_80", "gaus(0)", -1.4, 4.6);
    TF1 *gaus2_50_80 = new TF1("gaus2_50_80", "gaus(0)", -1.4, 4.6);
    TF1 *bg_50_80 = new TF1("bg_50_80", "pol0(0)", -1.4, 4.6);
    gaus1_50_80->SetParameter(0, corrFit_50_80->GetParameter(0));
    gaus1_50_80->SetParameter(1, corrFit_50_80->GetParameter(1));
    gaus1_50_80->SetParameter(2, corrFit_50_80->GetParameter(2));
    gaus2_50_80->SetParameter(0, corrFit_50_80->GetParameter(3));
    gaus2_50_80->SetParameter(1, corrFit_50_80->GetParameter(4));
    gaus2_50_80->SetParameter(2, corrFit_50_80->GetParameter(5));
    
    //bg->SetParameter(0, corrFit->GetParameter(6));    
    //bg->FixParameter(0, 0.25*(hhdphi_50_80->GetBinContent(1) + hhdphi_50_80->GetBinContent(8) +hhdphi_50_80->GetBinContent(9) +hhdphi_50_80->GetBinContent(16)));
    bg_50_80->FixParameter(0, 0.5*(hhdphi_50_80->GetBinContent(1) + hhdphi_50_80->GetBinContent(8)));

    TH1D* hhdphi_50_80_corr = (TH1D*)hhdphi_50_80->Clone("hhdphi_50_80_corr");
    hhdphi_50_80_corr->Add(bg_50_80, -1.0);
    hhdphi_50_80_corr->Scale(1.0/3.2); //correct for larger delta eta range "per Delta eta"
    for(int i = 1; i <= hhdphi_50_80_corr->GetXaxis()->GetNbins(); i++){
        hhdphi_50_80_corr->SetBinContent(i, hhdphi_50_80_corr->GetBinContent(i)/hhdphi_50_80_corr->GetXaxis()->GetBinWidth(i));
    }


    TPaveText *text = new TPaveText(0.4815, 0.7056, 0.8658, 0.8551, "NDC");
    text->AddText("ALICE Work in Progress");
    text->AddText("p-Pb #sqrt{#it{s}_{NN}} = 5 TeV");
    text->AddText("0%-20% Multiplicity");
    text->SetTextSizePixels(20);
    text->SetFillColor(kWhite);

    TPaveText *text2 = new TPaveText(0.6, 0.9, 0.85, 0.85, "NDC");
    text2->AddText("trigger: 4.0 < #it{p}_{T}^{h} < 8.0 GeV/c");
    text2->AddText("assoc: 2.0 < #it{p}_{T}^{a} < 4.0 GeV/c");
    text2->SetTextSizePixels(18);
    text2->SetFillColor(kWhite);

    TLegend  *legend = new TLegend(0.3791, 0.1518, 0.8772, 0.2688);
    legend->SetMargin(0.15);
    //legend->AddEntry(corrFit2, "Hadron-#phi(1020) Correlation", "l");
    legend->AddEntry(corrFit_0_20, "Hadron-hadron Correlations", "l");
    TCanvas *c0_20 = new TCanvas("c0_20", "c0_20", 50, 50, 720, 522);
    c0_20->cd();
    c0_20->SetMargin(0.12, 0.05, 0.1, 0.05);
    //hhdphi_0_20->GetYaxis()->SetRangeUser(0.03, 0.11);
    hhdphi_0_20_corr->GetYaxis()->SetRangeUser(-0.02, 0.2);
    hhdphi_0_20_corr->SetMarkerColor(kRed);
    hhdphi_0_20_corr->Draw("E0 X0");
    hhdphi_20_50_corr->SetMarkerColor(kOrange+8);
    hhdphi_20_50_corr->Draw("E0 X0 SAME");
    hhdphi_50_80_corr->Draw("E0 X0 SAME");
    
}
