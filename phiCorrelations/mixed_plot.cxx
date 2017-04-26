void mixed_plot(){
    TFile *histofile = new TFile("~/phiStudies/phiCorrelations_LHC16q_FAST_mixed_20170416.root");
    histofile->cd("PhiReconstruction");
    THnSparseF *ffMixedUS = (THnSparseF*) InvMass->FindObject("fDphiHPhiMixed");
   
    TFile* output = new TFile("mixed.root", "RECREATE");

    ffMixedUS->GetAxis(0)->SetRangeUser(4.0, 8.0);
    ffMixedUS->GetAxis(1)->SetRangeUser(2.0, 4.0);
    TH2D *twoCorr = ffMixedUS->Projection(2,3);
    twoCorr->SetName("twocorr_4trigger8_2phi4");
    twoCorr->Sumw2();

    Double_t centerPoint = 0.25*(twoCorr->GetBinContent(32, 16) + twoCorr->GetBinContent(32, 17) + twoCorr->GetBinContent(33, 16) + twoCorr->GetBinContent(33, 17));
    twoCorr->Scale(1.0/centerPoint);

    twoCorr->Write();

    ffMixedUS->GetAxis(1)->SetRangeUser(4.0, 8.0);
    TH2D *twoCorrHigh = ffMixedUS->Projection(2,3);
    twoCorrHigh->SetName("twocorr_4trigger8_4phi8");
    twoCorrHigh->Sumw2();
 
    Double_t centerPointHigh = 0.25*(twoCorrHigh->GetBinContent(32, 16) + twoCorrHigh->GetBinContent(32, 17) + twoCorrHigh->GetBinContent(33, 16) + twoCorrHigh->GetBinContent(33, 17));
    twoCorrHigh->Scale(1.0/centerPointHigh);

    twoCorrHigh->Write();

    TH1D *phi = twoCorr->ProjectionY("2phi4_phiproj");
    phi->Scale(2.0/(phi->GetBinContent(16) + phi->GetBinContent(17)));
    phi->Write();

    TH1D *phiHigh = twoCorrHigh->ProjectionY("4phi8_phiproj");
    phiHigh->Scale(2.0/(phiHigh->GetBinContent(16) + phiHigh->GetBinContent(17)));
    phiHigh->Write();

    ffMixedUS->GetAxis(0)->SetRangeUser(4.0, 8.0);
    ffMixedUS->GetAxis(1)->SetRangeUser(2.0, 4.0);
    ffMixedUS->GetAxis(4)->SetRangeUser(1.01, 1.03);
    TH2D *twoCorrPeak = ffMixedUS->Projection(2,3);
    twoCorrPeak->Sumw2();
    twoCorrPeak->SetName("twocorr_peak");
    twoCorrPeak->Scale(1.0/(twoCorrPeak->Integral()));
    twoCorrPeak->Write();

    ffMixedUS->GetAxis(4)->SetRangeUser(1.04, 1.06);
    TH2D *twoCorrRSB = ffMixedUS->Projection(2,3);
    twoCorrRSB->Sumw2();
    twoCorrRSB->SetName("twocorr_RSB");
    twoCorrRSB->Scale(1.0/(twoCorrRSB->Integral()));
    twoCorrRSB->Write();

    ffMixedUS->GetAxis(4)->SetRangeUser(0.995, 1.005);
    TH2D *twoCorrLSB = ffMixedUS->Projection(2,3);
    twoCorrLSB->Sumw2();
    twoCorrLSB->SetName("twocorr_LSB");
    twoCorrLSB->Scale(1.0/(twoCorrLSB->Integral()));
    twoCorrLSB->Write();

    TH2D *ratioR = twoCorrRSB->Clone("ratioR");
    ratioR->Divide(twoCorrPeak, twoCorrRSB);
    ratioR->Write();

    TH2D *ratioL = twoCorrLSB->Clone("ratioL");
    ratioL->Divide(twoCorrPeak, twoCorrLSB);
    ratioL->Write();

    TH2D *ratioSB = ratioL->Clone("ratioSB");
    ratioSB->Divide(ratioL, ratioR);
    ratioSB->Write();

    TH1D *dEta[3];
    TH1D *dPhiRestricted[3];
    TH1D *dPhi[3];

    dEta[0] = twoCorrPeak->ProjectionX();
    dEta[0]->SetName("dEta_peak");
    dEta[0]->Write();

    dEta[1] = twoCorrRSB->ProjectionX();
    dEta[1]->SetName("dEta_RSB");
    dEta[1]->Write();

    dEta[2] = twoCorrLSB->ProjectionX();
    dEta[2]->SetName("dEta_LSB");
    dEta[2]->Write();

    dPhi[0] = twoCorrPeak->ProjectionY();
    dPhi[0]->SetName("dPhi_peak");
    dPhi[0]->Write();

    dPhi[1] = twoCorrRSB->ProjectionY();
    dPhi[1]->SetName("dPhi_RSB");
    dPhi[1]->Write();

    dPhi[2] = twoCorrLSB->ProjectionY();
    dPhi[2]->SetName("dPhi_LSB");
    dPhi[2]->Write();

    TCanvas *cdEta = new TCanvas("cdEta", "cdEta", 50, 50, 600, 600);
    cdEta->cd();
    dEta[0]->SetLineColor(1);
    dEta[0]->Draw("SAME H");
    dEta[1]->SetLineColor(2);
    dEta[1]->Draw("SAME H");
    dEta[2]->SetLineColor(4);
    dEta[2]->Draw("SAME H");

    TCanvas *cdPhi = new TCanvas("cdPhi", "cdPhi", 50, 50, 600, 600);
    cdPhi->cd();
    dPhi[0]->SetLineColor(1);
    dPhi[0]->Draw("SAME H");
    dPhi[1]->SetLineColor(2);
    dPhi[1]->Draw("SAME H");
    dPhi[2]->SetLineColor(4);
    dPhi[2]->Draw("SAME H");


}
