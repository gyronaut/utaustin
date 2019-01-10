void compareDPhi(){
    TFile* filenarrow_0_20 = new TFile("~/phiStudies/results_newmult/US_syst_trig_4_8_assoc_2_4_smallmass12_hPhi_0_20.root");
    TH2D* narrow2D_0_20 = (TH2D*)filenarrow_0_20->Get("AvgUSsubhPhi2Dpeak");
    TH1D* narrow1D_0_20 = (TH1D*)narrow2D_0_20->ProjectionY("narrow1D_0_20", narrow2D_0_20->GetXaxis()->FindBin(-1.2), narrow2D_0_20->GetXaxis()->FindBin(1.2));
    TFile* filewide_0_20 = new TFile("~/phiStudies/results_newmult/US_syst_trig_4_8_assoc_2_4_mixcorr_hPhi_0_20.root");
    TH2D* wide2D_0_20 = (TH2D*)filewide_0_20->Get("AvgUSsubhPhi2Dpeak");
    TH1D* wide1D_0_20 = (TH1D*)wide2D_0_20->ProjectionY("wide1D", wide2D_0_20->GetXaxis()->FindBin(-1.2), wide2D_0_20->GetXaxis()->FindBin(1.2));

    TFile* filenarroweff_0_20 = new TFile("~/phiStudies/results_newmult/US_syst_trig_4_8_assoc_2_4_effcorr_smallmass12_hPhi_0_20.root");
    TH2D* narroweff2D_0_20 = (TH2D*)filenarroweff_0_20->Get("AvgUSsubhPhi2Dpeak");
    TH1D* narroweff1D_0_20 = (TH1D*)narroweff2D_0_20->ProjectionY("narroweff1D", narroweff2D_0_20->GetXaxis()->FindBin(-1.2), narroweff2D_0_20->GetXaxis()->FindBin(1.2));

    TFile* filewideeff_0_20 = new TFile("~/phiStudies/results_newmult/US_syst_trig_4_8_assoc_2_4_effcorr_hPhi_0_20.root");
    TH2D* wideeff2D_0_20 = (TH2D*)filewideeff_0_20->Get("AvgUSsubhPhi2Dpeak");
    TH1D* wideeff1D_0_20 = (TH1D*)wideeff2D_0_20->ProjectionY("wideeff1D", wideeff2D_0_20->GetXaxis()->FindBin(-1.2), wideeff2D_0_20->GetXaxis()->FindBin(1.2));

    //scale for invariant mass cut
    narrow1D_0_20->Scale(1./.803);
    narrow1D_0_20->SetMarkerColor(kRed+1);
    narrow1D_0_20->SetLineColor(kRed+1);
    narrow1D_0_20->SetLineWidth(2);
    narrow1D_0_20->SetMarkerSize(2);
    narrow1D_0_20->SetMarkerStyle(20);
    narrow1D_0_20->SetStats(kFALSE);
    narroweff1D_0_20->Scale(1./.803);
    narroweff1D_0_20->SetMarkerColor(kRed+1);
    narroweff1D_0_20->SetLineColor(kRed+1);
    narroweff1D_0_20->SetLineWidth(2);
    narroweff1D_0_20->SetMarkerSize(2);
    narroweff1D_0_20->SetMarkerStyle(22);
    narroweff1D_0_20->SetStats(kFALSE);

    wide1D_0_20->Scale(1./.897);
    wide1D_0_20->SetMarkerColor(kBlue+1);
    wide1D_0_20->SetLineColor(kBlue+1);
    wide1D_0_20->SetLineWidth(2);
    wide1D_0_20->SetMarkerSize(2);
    wide1D_0_20->SetMarkerStyle(21);
    wideeff1D_0_20->Scale(1./.897);
    wideeff1D_0_20->SetMarkerColor(kBlue+1);
    wideeff1D_0_20->SetLineColor(kBlue+1);
    wideeff1D_0_20->SetLineWidth(2);
    wideeff1D_0_20->SetMarkerSize(2);
    wideeff1D_0_20->SetMarkerStyle(23);

    TLegend *legend_0_20 = new TLegend(0.49, 0.68, 0.88, 0.88);
    legend_0_20->AddEntry(narrow1D_0_20, "Narrow (1.014-1.026)", "le");
    legend_0_20->AddEntry(wide1D_0_20, "Wide (1.010-1.030)", "le");

    TCanvas* cnoeff_0_20 = new TCanvas("cnoeff_0_20", "cnoeff_0_20", 50, 50, 600, 600);
    cnoeff_0_20->cd();
    narrow1D_0_20->Draw("P E");
    wide1D_0_20->Draw("P E SAME");
    legend_0_20->Draw();
    
    TCanvas* ceff_0_20 = new TCanvas("ceff_0_20", "ceff_0_20", 50, 50, 600, 600);
    ceff_0_20->cd();
    narroweff1D_0_20->Draw("P E");
    wideeff1D_0_20->Draw("P E SAME");
    legend_0_20->Draw();


    //comparison of all together, scaled to integral of 1
    TH1D* narrowNorm_0_20 = (TH1D*)narrow1D_0_20->Clone("narrowNorm_0_20");
    narrowNorm_0_20->Scale(1.0/narrowNorm_0_20->Integral());
    TH1D* wideNorm_0_20 = (TH1D*)wide1D_0_20->Clone("wideNorm_0_20");
    wideNorm_0_20->Scale(1.0/wideNorm_0_20->Integral());
    TH1D* narrowNormEff_0_20 = (TH1D*)narroweff1D_0_20->Clone("narrowNormEff_0_20");
    narrowNormEff_0_20->Scale(1.0/narrowNormEff_0_20->Integral());
    TH1D* wideNormEff_0_20 = (TH1D*)wideeff1D_0_20->Clone("wideNormEff_0_20");
    wideNormEff_0_20->Scale(1.0/wideNormEff_0_20->Integral());

    TLegend *legendAll_0_20 = new TLegend(0.49, 0.68, 0.88, 0.88);
    legendAll_0_20->AddEntry(narrowNorm_0_20, "Narrow, No Eff.", "lep");
    legendAll_0_20->AddEntry(wideNorm_0_20, "Wide, No Eff.", "lep");
    legendAll_0_20->AddEntry(narrowNormEff_0_20, "Narrow, With Eff.", "lep");
    legendAll_0_20->AddEntry(wideNormEff_0_20, "Wide, With Eff.", "lep");

    TCanvas *call_0_20 = new TCanvas("call_0_20", "call_0_20", 50, 50, 600, 600);
    call_0_20->cd();
    narrowNorm_0_20->Draw("P E");
    narrowNormEff_0_20->Draw("SAME P E");
    wideNorm_0_20->Draw("SAME P E");
    wideNormEff_0_20->Draw("SAME P E");
    legendAll_0_20->Draw();

    TFile* filenarrow_20_50 = new TFile("~/phiStudies/results_newmult/US_syst_trig_4_8_assoc_2_4_smallmass12_hPhi_20_50.root");
    TH2D* narrow2D_20_50 = (TH2D*)filenarrow_20_50->Get("AvgUSsubhPhi2Dpeak");
    TH1D* narrow1D_20_50 = (TH1D*)narrow2D_20_50->ProjectionY("narrow1D_20_50", narrow2D_20_50->GetXaxis()->FindBin(-1.2), narrow2D_20_50->GetXaxis()->FindBin(1.2));
    TFile* filewide_20_50 = new TFile("~/phiStudies/results_newmult/US_syst_trig_4_8_assoc_2_4_mixcorr_hPhi_20_50.root");
    TH2D* wide2D_20_50 = (TH2D*)filewide_20_50->Get("AvgUSsubhPhi2Dpeak");
    TH1D* wide1D_20_50 = (TH1D*)wide2D_20_50->ProjectionY("wide1D", wide2D_20_50->GetXaxis()->FindBin(-1.2), wide2D_20_50->GetXaxis()->FindBin(1.2));

    TFile* filenarroweff_20_50 = new TFile("~/phiStudies/results_newmult/US_syst_trig_4_8_assoc_2_4_effcorr_smallmass12_hPhi_20_50.root");
    TH2D* narroweff2D_20_50 = (TH2D*)filenarroweff_20_50->Get("AvgUSsubhPhi2Dpeak");
    TH1D* narroweff1D_20_50 = (TH1D*)narroweff2D_20_50->ProjectionY("narroweff1D", narroweff2D_20_50->GetXaxis()->FindBin(-1.2), narroweff2D_20_50->GetXaxis()->FindBin(1.2));

    TFile* filewideeff_20_50 = new TFile("~/phiStudies/results_newmult/US_syst_trig_4_8_assoc_2_4_effcorr_hPhi_20_50.root");
    TH2D* wideeff2D_20_50 = (TH2D*)filewideeff_20_50->Get("AvgUSsubhPhi2Dpeak");
    TH1D* wideeff1D_20_50 = (TH1D*)wideeff2D_20_50->ProjectionY("wideeff1D", wideeff2D_20_50->GetXaxis()->FindBin(-1.2), wideeff2D_20_50->GetXaxis()->FindBin(1.2));

    //scale for invariant mass cut
    narrow1D_20_50->Scale(1./.803);
    narrow1D_20_50->SetMarkerColor(kRed+1);
    narrow1D_20_50->SetLineColor(kRed+1);
    narrow1D_20_50->SetLineWidth(2);
    narrow1D_20_50->SetMarkerSize(2);
    narrow1D_20_50->SetMarkerStyle(20);
    narrow1D_20_50->SetStats(kFALSE);
    narroweff1D_20_50->Scale(1./.803);
    narroweff1D_20_50->SetMarkerColor(kRed+1);
    narroweff1D_20_50->SetLineColor(kRed+1);
    narroweff1D_20_50->SetLineWidth(2);
    narroweff1D_20_50->SetMarkerSize(2);
    narroweff1D_20_50->SetMarkerStyle(22);
    narroweff1D_20_50->SetStats(kFALSE);

    wide1D_20_50->Scale(1./.897);
    wide1D_20_50->SetMarkerColor(kBlue+1);
    wide1D_20_50->SetLineColor(kBlue+1);
    wide1D_20_50->SetLineWidth(2);
    wide1D_20_50->SetMarkerSize(2);
    wide1D_20_50->SetMarkerStyle(21);
    wideeff1D_20_50->Scale(1./.897);
    wideeff1D_20_50->SetMarkerColor(kBlue+1);
    wideeff1D_20_50->SetLineColor(kBlue+1);
    wideeff1D_20_50->SetLineWidth(2);
    wideeff1D_20_50->SetMarkerSize(2);
    wideeff1D_20_50->SetMarkerStyle(23);

    TLegend *legend_20_50 = new TLegend(0.49, 0.68, 0.88, 0.88);
    legend_20_50->AddEntry(narrow1D_20_50, "Narrow (1.014-1.026)", "le");
    legend_20_50->AddEntry(wide1D_20_50, "Wide (1.010-1.030)", "le");

    TCanvas* cnoeff_20_50 = new TCanvas("cnoeff_20_50", "cnoeff_20_50", 50, 50, 600, 600);
    cnoeff_20_50->cd();
    narrow1D_20_50->Draw("P E");
    wide1D_20_50->Draw("P E SAME");
    legend_20_50->Draw();
    
    TCanvas* ceff_20_50 = new TCanvas("ceff_20_50", "ceff_20_50", 50, 50, 600, 600);
    ceff_20_50->cd();
    narroweff1D_20_50->Draw("P E");
    wideeff1D_20_50->Draw("P E SAME");
    legend_20_50->Draw();


    //comparison of all together, scaled to integral of 1
    TH1D* narrowNorm_20_50 = (TH1D*)narrow1D_20_50->Clone("narrowNorm_20_50");
    narrowNorm_20_50->Scale(1.0/narrowNorm_20_50->Integral());
    TH1D* wideNorm_20_50 = (TH1D*)wide1D_20_50->Clone("wideNorm_20_50");
    wideNorm_20_50->Scale(1.0/wideNorm_20_50->Integral());
    TH1D* narrowNormEff_20_50 = (TH1D*)narroweff1D_20_50->Clone("narrowNormEff_20_50");
    narrowNormEff_20_50->Scale(1.0/narrowNormEff_20_50->Integral());
    TH1D* wideNormEff_20_50 = (TH1D*)wideeff1D_20_50->Clone("wideNormEff_20_50");
    wideNormEff_20_50->Scale(1.0/wideNormEff_20_50->Integral());

    TLegend *legendAll_20_50 = new TLegend(0.49, 0.68, 0.88, 0.88);
    legendAll_20_50->AddEntry(narrowNorm_20_50, "Narrow, No Eff.", "lep");
    legendAll_20_50->AddEntry(wideNorm_20_50, "Wide, No Eff.", "lep");
    legendAll_20_50->AddEntry(narrowNormEff_20_50, "Narrow, With Eff.", "lep");
    legendAll_20_50->AddEntry(wideNormEff_20_50, "Wide, With Eff.", "lep");

    TCanvas *call_20_50 = new TCanvas("call_20_50", "call_20_50", 50, 50, 600, 600);
    call_20_50->cd();
    narrowNorm_20_50->Draw("P E");
    narrowNormEff_20_50->Draw("SAME P E");
    wideNorm_20_50->Draw("SAME P E");
    wideNormEff_20_50->Draw("SAME P E");
    legendAll_20_50->Draw();

    TFile* filenarrow_50_80 = new TFile("~/phiStudies/results_newmult/US_syst_trig_4_8_assoc_2_4_smallmass12_hPhi_50_80.root");
    TH2D* narrow2D_50_80 = (TH2D*)filenarrow_50_80->Get("AvgUSsubhPhi2Dpeak");
    TH1D* narrow1D_50_80 = (TH1D*)narrow2D_50_80->ProjectionY("narrow1D_50_80", narrow2D_50_80->GetXaxis()->FindBin(-1.2), narrow2D_50_80->GetXaxis()->FindBin(1.2));
    TFile* filewide_50_80 = new TFile("~/phiStudies/results_newmult/US_syst_trig_4_8_assoc_2_4_mixcorr_hPhi_50_80.root");
    TH2D* wide2D_50_80 = (TH2D*)filewide_50_80->Get("AvgUSsubhPhi2Dpeak");
    TH1D* wide1D_50_80 = (TH1D*)wide2D_50_80->ProjectionY("wide1D", wide2D_50_80->GetXaxis()->FindBin(-1.2), wide2D_50_80->GetXaxis()->FindBin(1.2));

    TFile* filenarroweff_50_80 = new TFile("~/phiStudies/results_newmult/US_syst_trig_4_8_assoc_2_4_effcorr_smallmass12_hPhi_50_80.root");
    TH2D* narroweff2D_50_80 = (TH2D*)filenarroweff_50_80->Get("AvgUSsubhPhi2Dpeak");
    TH1D* narroweff1D_50_80 = (TH1D*)narroweff2D_50_80->ProjectionY("narroweff1D", narroweff2D_50_80->GetXaxis()->FindBin(-1.2), narroweff2D_50_80->GetXaxis()->FindBin(1.2));

    TFile* filewideeff_50_80 = new TFile("~/phiStudies/results_newmult/US_syst_trig_4_8_assoc_2_4_effcorr_hPhi_50_80.root");
    TH2D* wideeff2D_50_80 = (TH2D*)filewideeff_50_80->Get("AvgUSsubhPhi2Dpeak");
    TH1D* wideeff1D_50_80 = (TH1D*)wideeff2D_50_80->ProjectionY("wideeff1D", wideeff2D_50_80->GetXaxis()->FindBin(-1.2), wideeff2D_50_80->GetXaxis()->FindBin(1.2));

    //scale for invariant mass cut
    narrow1D_50_80->Scale(1./.803);
    narrow1D_50_80->SetMarkerColor(kRed+1);
    narrow1D_50_80->SetLineColor(kRed+1);
    narrow1D_50_80->SetLineWidth(2);
    narrow1D_50_80->SetMarkerSize(2);
    narrow1D_50_80->SetMarkerStyle(20);
    narrow1D_50_80->SetStats(kFALSE);
    narroweff1D_50_80->Scale(1./.803);
    narroweff1D_50_80->SetMarkerColor(kRed+1);
    narroweff1D_50_80->SetLineColor(kRed+1);
    narroweff1D_50_80->SetLineWidth(2);
    narroweff1D_50_80->SetMarkerSize(2);
    narroweff1D_50_80->SetMarkerStyle(22);
    narroweff1D_50_80->SetStats(kFALSE);

    wide1D_50_80->Scale(1./.897);
    wide1D_50_80->SetMarkerColor(kBlue+1);
    wide1D_50_80->SetLineColor(kBlue+1);
    wide1D_50_80->SetLineWidth(2);
    wide1D_50_80->SetMarkerSize(2);
    wide1D_50_80->SetMarkerStyle(21);
    wideeff1D_50_80->Scale(1./.897);
    wideeff1D_50_80->SetMarkerColor(kBlue+1);
    wideeff1D_50_80->SetLineColor(kBlue+1);
    wideeff1D_50_80->SetLineWidth(2);
    wideeff1D_50_80->SetMarkerSize(2);
    wideeff1D_50_80->SetMarkerStyle(23);

    TLegend *legend_50_80 = new TLegend(0.49, 0.68, 0.88, 0.88);
    legend_50_80->AddEntry(narrow1D_50_80, "Narrow (1.014-1.026)", "le");
    legend_50_80->AddEntry(wide1D_50_80, "Wide (1.010-1.030)", "le");

    TCanvas* cnoeff_50_80 = new TCanvas("cnoeff_50_80", "cnoeff_50_80", 50, 50, 600, 600);
    cnoeff_50_80->cd();
    narrow1D_50_80->Draw("P E");
    wide1D_50_80->Draw("P E SAME");
    legend_50_80->Draw();
    
    TCanvas* ceff_50_80 = new TCanvas("ceff_50_80", "ceff_50_80", 50, 50, 600, 600);
    ceff_50_80->cd();
    narroweff1D_50_80->Draw("P E");
    wideeff1D_50_80->Draw("P E SAME");
    legend_50_80->Draw();


    //comparison of all together, scaled to integral of 1
    TH1D* narrowNorm_50_80 = (TH1D*)narrow1D_50_80->Clone("narrowNorm_50_80");
    narrowNorm_50_80->Scale(1.0/narrowNorm_50_80->Integral());
    TH1D* wideNorm_50_80 = (TH1D*)wide1D_50_80->Clone("wideNorm_50_80");
    wideNorm_50_80->Scale(1.0/wideNorm_50_80->Integral());
    TH1D* narrowNormEff_50_80 = (TH1D*)narroweff1D_50_80->Clone("narrowNormEff_50_80");
    narrowNormEff_50_80->Scale(1.0/narrowNormEff_50_80->Integral());
    TH1D* wideNormEff_50_80 = (TH1D*)wideeff1D_50_80->Clone("wideNormEff_50_80");
    wideNormEff_50_80->Scale(1.0/wideNormEff_50_80->Integral());

    TLegend *legendAll_50_80 = new TLegend(0.49, 0.68, 0.88, 0.88);
    legendAll_50_80->AddEntry(narrowNorm_50_80, "Narrow, No Eff.", "lep");
    legendAll_50_80->AddEntry(wideNorm_50_80, "Wide, No Eff.", "lep");
    legendAll_50_80->AddEntry(narrowNormEff_50_80, "Narrow, With Eff.", "lep");
    legendAll_50_80->AddEntry(wideNormEff_50_80, "Wide, With Eff.", "lep");

    TCanvas *call_50_80 = new TCanvas("call_50_80", "call_50_80", 50, 50, 600, 600);
    call_50_80->cd();
    narrowNorm_50_80->Draw("P E");
    narrowNormEff_50_80->Draw("SAME P E");
    wideNorm_50_80->Draw("SAME P E");
    wideNormEff_50_80->Draw("SAME P E");
    legendAll_50_80->Draw();


}
