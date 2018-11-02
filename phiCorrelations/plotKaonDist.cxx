void plotKaonDist(TString inputname){
    TFile* infile = new TFile(inputname.Data());
    TList* list = (TList*)infile->Get("phiCorr_mult_0_20");

    THnSparseF* USdist = (THnSparseF*)list->FindObject("fkkUSDist");
    USdist->GetAxis(0)->SetRangeUser(2.0, 4.0);
    TH1D* invmass = USdist->Projection(1);
    invmass->SetTitle("KK Invariant mass for 2.0 < p_{T} < 4.0");
    invmass->SetStats(kFALSE);

    THnSparseF* LSdist = (THnSparseF*)list->FindObject("fkkLSDist");
    LSdist->GetAxis(0)->SetRangeUser(2.0, 4.0);
    TH1D* lsmass = LSdist->Projection(1);

    TH3D* kaonDist = (TH3D*)list->FindObject("fKaonDist");
    TH3D* kaonPID = (TH3D*)list->FindObject("fKaonPID");

    TH2D* kaonTPC = (TH2D*)kaonPID->Project3D("yx");
    kaonTPC->SetTitle("n#sigma_{TPC} vs p_{T}");

    TH2D* kaonTOF = (TH2D*)kaonPID->Project3D("zx");
    kaonTOF->SetTitle("n#sigma_{TOF} vs p_{T}");

    kaonPID->GetXaxis()->SetRangeUser(1.0, 5.0);
    TH2D* kaonTPCTOF = (TH2D*)kaonPID->Project3D("yz");
    kaonTPCTOF->SetTitle("n#sigma_{TPC} vs. n#sigma_{TOF} for 1.0 < p_{T} < 5.0 GeV/c");

    TH1D* kaonPT = (TH1D*)kaonDist->ProjectionX("Kaon Candidate p_{T}");
    kaonPT->SetTitle("Kaon Candidate p_{T}");
    TH1D* kaonPhi = (TH1D*)kaonDist->ProjectionY("Kaon Candidate #varphi");
    kaonPhi->SetTitle("Kaon Candidate #varphi");
    TH1D* kaonEta = (TH1D*)kaonDist->ProjectionZ("Kaon Candidate #eta");
    kaonEta->SetTitle("Kaon Candidate #eta");

    TCanvas* c1 = new TCanvas("c1", "c1", 50, 50, 600, 600);
    c1->cd()->SetLogz();
    kaonTPC->Draw("colz");

    TCanvas* c2 = new TCanvas("c2", "c2", 50, 50, 600, 600);
    c2->cd()->SetLogz();
    kaonTOF->Draw("colz");

    TCanvas* c3 = new TCanvas("c3", "c3", 50, 50, 600, 600);
    c3->cd()->SetLogz();
    kaonTPCTOF->Draw("colz");

    TCanvas* c4 = new TCanvas("c4", "c4", 50, 50, 600, 600);
    c4->cd();
    kaonPT->Draw();

    TCanvas* c5 = new TCanvas("c5", "c5", 50, 50, 600, 600);
    c5->cd();
    kaonPhi->Draw();

    TCanvas* c6 = new TCanvas("c6", "c6", 50, 50, 600, 600);
    c6->cd();
    kaonEta->Draw();

    //invmass fitting and plotting
    TF1* massfit = new TF1("massfit", "[0]*TMath::Voigt(x - [1], [2], [3], 4) + pol2(4)", 1.0, 1.04);
    massfit->SetParLimits(0, 1.0, 100000);
    massfit->SetParameter(1, 1.0196);
    massfit->SetParameter(2, 0.0002);
    massfit->FixParameter(3, 0.00426);
    //massfit->SetParLimits(5, -1000, 0.0);
    massfit->SetParameter(4, -100000);
    massfit->SetParameter(5, 100000);
    massfit->SetParameter(6, -100000);
    invmass->Fit(massfit);

    TF1* bgfit = new TF1("bgfit", "pol2(0)", 1.0, 1.04);
    bgfit->SetParameters(massfit->GetParameter(4), massfit->GetParameter(5), massfit->GetParameter(6));
    bgfit->SetLineColor(kBlue+1);

    float scale = invmass->Integral(invmass->GetXaxis()->FindBin(1.040), invmass->GetXaxis()->FindBin(1.06))/lsmass->Integral(lsmass->GetXaxis()->FindBin(1.040), lsmass->GetXaxis()->FindBin(1.060));
    
    lsmass->Scale(scale);
    lsmass->SetLineStyle(4);
    lsmass->SetLineColor(kMagenta+1);

    float bg = bgfit->Integral(1.01, 1.03);
    float signal = invmass->Integral(invmass->GetXaxis()->FindBin(1.01), invmass->GetXaxis()->FindBin(1.03), "width") - bg;

    float purity = signal/(signal+bg);
    float significance = (signal)/TMath::Sqrt(signal + bg);

    TPaveText* text = new TPaveText(0.51, 0.69, 0.86, 0.86, "NDC");
    text->AddText(Form("signal: %.4e", signal));
    text->AddText(Form("BG: %.4e", bg));
    text->AddText(Form("purity: %f", purity));
    text->SetFillColor(kWhite);

    TCanvas* c7 = new TCanvas("c7", "c7", 50, 50, 600, 600);
    c7->cd();
    invmass->Draw();
    bgfit->Draw("SAME");
    lsmass->Draw("SAME");
    text->Draw();

}
