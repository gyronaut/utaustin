TH2D* makeCorrections(THnSparse* same, THnSparse* mixed, Float_t lowmass, Float_t highmass, float totalTrigSame){
    same->GetAxis(3)->SetRangeUser(lowmass, highmass);
    mixed->GetAxis(3)->SetRangeUser(lowmass, highmass);
    same->GetAxis(2)->SetRangeUser(-10., 10.);
    mixed->GetAxis(2)->SetRangeUser(-10., 10.);
    TH2D* same2D = same->Projection(0, 1);
    same2D->Sumw2();
    TH2D* mix2D = mixed->Projection(0, 1);
    mix2D->Sumw2();

    same2D->RebinX(4);
    mix2D->RebinX(4);
    same2D->RebinY(4);
    mix2D->RebinY(4);

    same2D->GetXaxis()->SetRange(1,same2D->GetXaxis()->GetNbins());
    mix2D->GetXaxis()->SetRange(1,mix2D->GetXaxis()->GetNbins());

    Float_t scale = 0.0;

    scale = 0.5*(float)(mix2D->GetBinContent(mix2D->GetXaxis()->FindBin(0.01), mix2D->GetYaxis()->FindBin(0)) + mix2D->GetBinContent(mix2D->GetXaxis()->FindBin(-0.01), mix2D->GetYaxis()->FindBin(0.0)));
    printf("scale: %e \n", scale);
    same2D->Divide(mix2D);
    same2D->Scale(scale);
    same2D->Scale(1.0/totalTrigSame);

    same->GetAxis(3)->SetRange(0,0);
    mixed->GetAxis(3)->SetRange(0,0);
    return same2D;
}

makeSimpleMixCorrections(string inputName, float trigPTLow, float trigPTHigh, float assocPTLow, float assocPTHigh){
    TFile *histoFile = new TFile(inputName.c_str()); 
    string mult = inputName.substr(inputName.find("_", inputName.find("_")+1), inputName.find(".") - inputName.find("_", inputName.find("_")+1));
    TList* list = (TList*) histoFile->Get(Form("phiCorr_mult%s", mult.c_str()));

    THnSparseF *dphiHPhi = (THnSparseF *)list->FindObject("fDphiHPhi");
    THnSparseF *dphiHKK = (THnSparseF *)list->FindObject("fDphiHKK");
    THnSparseF *dphiHPhiMixed = (THnSparseF *)list->FindObject("fDphiHPhiMixed");
    THnSparseF *dphiHKKMixed = (THnSparseF *)list->FindObject("fDphiHKKMixed");

    TH2D *trigSameUSDist = (TH2D*)list->FindObject("fTrigSameUSDist");
    TH2D *trigSameLSDist = (TH2D*)list->FindObject("fTrigSameLSDist");

    float totalTrigSameUS = (float)trigSameUSDist->Integral(trigSameUSDist->GetXaxis()->FindBin(trigPTLow), trigSameUSDist->GetXaxis()->FindBin(trigPTHigh), 1, trigSameUSDist->GetYaxis()->GetNbins());
    float totalTrigSameLS = (float)trigSameUSDist->Integral(trigSameLSDist->GetXaxis()->FindBin(trigPTLow), trigSameLSDist->GetXaxis()->FindBin(trigPTHigh), 1, trigSameLSDist->GetYaxis()->GetNbins());


    dphiHPhi->GetAxis(0)->SetRangeUser(trigPTLow, trigPTHigh); 
    dphiHPhi->GetAxis(1)->SetRangeUser(assocPTLow, assocPTHigh); 
    dphiHKK->GetAxis(0)->SetRangeUser(trigPTLow, trigPTHigh);
    dphiHKK->GetAxis(1)->SetRangeUser(assocPTLow, assocPTHigh);
   
    dphiHPhiMixed->GetAxis(0)->SetRangeUser(trigPTLow, trigPTHigh); 
    dphiHPhiMixed->GetAxis(1)->SetRangeUser(assocPTLow, assocPTHigh); 
    dphiHKKMixed->GetAxis(0)->SetRangeUser(trigPTLow, trigPTHigh);
    dphiHKKMixed->GetAxis(1)->SetRangeUser(assocPTLow, assocPTHigh);

    dphiHPhi->GetAxis(5)->SetRange(1,dphiHPhi->GetAxis(5)->GetNbins());
    dphiHPhi->GetAxis(4)->SetRangeUser(-10.0, 10.0);

    dphiHKK->GetAxis(5)->SetRange(1,dphiHKK->GetAxis(5)->GetNbins());
    dphiHKK->GetAxis(4)->SetRangeUser(-10.0, 10.0);

    dphiHPhiMixed->GetAxis(5)->SetRange(1,dphiHPhiMixed->GetAxis(5)->GetNbins());
    dphiHPhiMixed->GetAxis(4)->SetRangeUser(-10.0, 10.0);

    dphiHKKMixed->GetAxis(5)->SetRange(1,dphiHKKMixed->GetAxis(5)->GetNbins());
    dphiHKKMixed->GetAxis(4)->SetRangeUser(-10.0, 10.0);

    Int_t axes[] = {2,3,4,5};

    THnSparseF* hPhi = dphiHPhi->Projection(4, axes);
    THnSparseF* hKK = dphiHKK->Projection(4, axes);
    THnSparseF* hPhiMixed = dphiHPhiMixed->Projection(4, axes);
    THnSparseF* hKKMixed = dphiHKKMixed->Projection(4, axes);


    TH2D* hPhi2Dpeak = makeCorrections(hPhi, hPhiMixed, 1.010, 1.030, totalTrigSameUS);
    hPhi2Dpeak->SetName("hPhi2Dpeak");
    TH2D* hKK2Dpeak = makeCorrections(hKK, hKKMixed, 1.010, 1.030, totalTrigSameLS);
    hKK2Dpeak->SetName("hKK2Dpeak");
    TH2D* hPhi2DRside = makeCorrections(hPhi, hPhiMixed, 1.040, 1.060, totalTrigSameUS);
    hPhi2DRside->SetName("hPhi2DRside");
    TH2D* hKK2DRside = makeCorrections(hKK, hKKMixed, 1.040, 1.060, totalTrigSameLS);
    hKK2DRside->SetName("hKK2DRside");
    TH2D* hPhi2DLside = makeCorrections(hPhi, hPhiMixed, 0.995, 1.005, totalTrigSameUS);
    hPhi2DLside->SetName("hPhi2DLside");
    TH2D* hKK2DLside = makeCorrections(hKK, hKKMixed, 0.995, 1.005, totalTrigSameLS);
    hKK2DLside->SetName("hKK2DLside");

    //Create some uncorrected same/mixed event 2D histos
    hPhi->GetAxis(3)->SetRangeUser(1.010, 1.030);
    hPhiMixed->GetAxis(3)->SetRangeUser(1.010, 1.030);
    hKK->GetAxis(3)->SetRangeUser(1.010, 1.030);
    hKKMixed->GetAxis(3)->SetRangeUser(1.010, 1.030);
    TH2D* uncorrhPhi2Dpeak = hPhi->Projection(0,1);
    uncorrhPhi2Dpeak->Sumw2();
    uncorrhPhi2Dpeak->SetName("uncorrhPhi2Dpeak");
    TH2D* uncorrhKK2Dpeak = hKK->Projection(0,1);
    uncorrhKK2Dpeak->Sumw2();
    uncorrhKK2Dpeak->SetName("uncorrhKK2Dpeak");
    TH2D* uncorrhPhiMixed2Dpeak = hPhiMixed->Projection(0,1);
    uncorrhPhiMixed2Dpeak->Sumw2();
    uncorrhPhiMixed2Dpeak->SetName("uncorrhPhiMixed2Dpeak");
    TH2D* uncorrhKKMixed2Dpeak = hKKMixed->Projection(0,1);
    uncorrhKKMixed2Dpeak->Sumw2();
    uncorrhKKMixed2Dpeak->SetName("uncorrhKKMixed2Dpeak");

    hPhi->GetAxis(3)->SetRangeUser(1.040, 1.060);
    hPhiMixed->GetAxis(3)->SetRangeUser(1.040, 1.060);
    hKK->GetAxis(3)->SetRangeUser(1.040, 1.060);
    hKKMixed->GetAxis(3)->SetRangeUser(1.040, 1.060);
    TH2D* uncorrhPhi2DRside = hPhi->Projection(0,1);
    uncorrhPhi2DRside->Sumw2();
    uncorrhPhi2DRside->SetName("uncorrhPhi2DRside");
    TH2D* uncorrhKK2DRside = hKK->Projection(0,1);
    uncorrhKK2DRside->Sumw2();
    uncorrhKK2DRside->SetName("uncorrhKK2DRside");
    TH2D* uncorrhPhiMixed2DRside = hPhiMixed->Projection(0,1);
    uncorrhPhiMixed2DRside->Sumw2();
    uncorrhPhiMixed2DRside->SetName("uncorrhPhiMixed2DRside");
    TH2D* uncorrhKKMixed2DRside = hKKMixed->Projection(0,1);
    uncorrhKKMixed2DRside->Sumw2();
    uncorrhKKMixed2DRside->SetName("uncorrhKKMixed2DRside");

    hPhi->GetAxis(3)->SetRangeUser(0.995, 1.005);
    hPhiMixed->GetAxis(3)->SetRangeUser(0.995, 1.005);
    hKK->GetAxis(3)->SetRangeUser(0.995, 1.005);
    hKKMixed->GetAxis(3)->SetRangeUser(0.995, 1.005);
    TH2D* uncorrhPhi2DLside = hPhi->Projection(0,1);
    uncorrhPhi2DLside->Sumw2();
    uncorrhPhi2DLside->SetName("uncorrhPhi2DLside");
    TH2D* uncorrhKK2DLside = hKK->Projection(0,1);
    uncorrhKK2DLside->Sumw2();
    uncorrhKK2DLside->SetName("uncorrhKK2DLside");
    TH2D* uncorrhPhiMixed2DLside = hPhiMixed->Projection(0,1);
    uncorrhPhiMixed2DLside->Sumw2();
    uncorrhPhiMixed2DLside->SetName("uncorrhPhiMixed2DLside");
    TH2D* uncorrhKKMixed2DLside = hKKMixed->Projection(0,1);
    uncorrhKKMixed2DLside->Sumw2();
    uncorrhKKMixed2DLside->SetName("uncorrhKKMixed2DLside");

    //Make some ratio plots of mixed event US over LS for sidbeand and peak regions
    //also make projections onto delta eta and delta phi (on restricted delta eta range)
    TH2D* mixedratioRSB = uncorrhPhiMixed2DRside->Clone("mixedratioRSB");
    mixedratioRSB->Divide(uncorrhKKMixed2DRside);
    TH1D* mixedratioRSBdeta = mixedratioRSB->ProjectionX("mixedratioRSBdeta", 1, mixedratioRSB->GetYaxis()->GetNbins());
    TH1D* mixedratioRSBdphi = mixedratioRSB->ProjectionY("mixedratioRSBdphi", mixedratioRSB->GetXaxis()->FindBin(-1.2), mixedratioRSB->GetXaxis()->FindBin(1.2));
    mixedratioRSB->GetXaxis()->SetRange(0,0);

    TH2D* mixedratioLSB = uncorrhPhiMixed2DLside->Clone("mixedratioLSB");
    mixedratioLSB->Divide(uncorrhKKMixed2DLside);
    TH1D* mixedratioLSBdeta = mixedratioLSB->ProjectionX("mixedratioLSBdeta", 1, mixedratioLSB->GetYaxis()->GetNbins());
    TH1D* mixedratioLSBdphi = mixedratioLSB->ProjectionY("mixedratioLSBdphi", mixedratioLSB->GetXaxis()->FindBin(-1.2), mixedratioLSB->GetXaxis()->FindBin(1.2));
    mixedratioLSB->GetXaxis()->SetRange(0,0);

    TH2D* mixedratioPeak = uncorrhPhiMixed2Dpeak->Clone("mixedratioPeak");
    mixedratioPeak->Divide(uncorrhKKMixed2Dpeak);
    TH1D* mixedratioPeakdeta = mixedratioPeak->ProjectionX("mixedratioPeakdeta", 1, mixedratioPeak->GetYaxis()->GetNbins());
    TH1D* mixedratioPeakdphi = mixedratioPeak->ProjectionY("mixedratioPeakdphi", mixedratioPeak->GetXaxis()->FindBin(-1.2), mixedratioPeak->GetXaxis()->FindBin(1.2));
    mixedratioPeak->GetXaxis()->SetRange(0,0);


    TFile* output = new TFile(Form("trig_%i_%i_assoc_%i_%i_simplecorr_%s", (int)trigPTLow, (int)trigPTHigh, (int)assocPTLow, (int)assocPTHigh, inputName.c_str()), "RECREATE");
    hPhi2Dpeak->Write();
    hKK2Dpeak->Write();
    hPhi2DRside->Write();
    hKK2DRside->Write();
    hPhi2DLside->Write();
    hKK2DLside->Write();
    uncorrhPhi2Dpeak->Write();
    uncorrhKK2Dpeak->Write();
    uncorrhPhiMixed2Dpeak->Write();
    uncorrhKKMixed2Dpeak->Write();
    uncorrhPhi2DRside->Write();
    uncorrhKK2DRside->Write();
    uncorrhPhiMixed2DRside->Write();
    uncorrhKKMixed2DRside->Write();
    uncorrhPhi2DLside->Write();
    uncorrhKK2DLside->Write();
    uncorrhPhiMixed2DLside->Write();
    uncorrhKKMixed2DLside->Write();

    mixedratioRSB->Write();
    mixedratioRSBdeta->Write();
    mixedratioRSBdphi->Write();
    mixedratioLSB->Write();
    mixedratioLSBdeta->Write();
    mixedratioLSBdphi->Write();
    mixedratioPeak->Write();
    mixedratioPeakdeta->Write();
    mixedratioPeakdphi->Write();

}
