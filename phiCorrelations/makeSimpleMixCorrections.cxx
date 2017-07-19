TH2D* makeCorrections(THnSparse* same, THnSparse* mixed, Float_t lowmass, Float_t highmass){
    same->GetAxis(3)->SetRangeUser(lowmass, highmass);
    mixed->GetAxis(3)->SetRangeUser(lowmass, highmass);
    same->GetAxis(2)->SetRangeUser(-10., 10.);
    mixed->GetAxis(2)->SetRangeUser(-10., 10.);
    TH2D* same2D = same->Projection(0, 1);
    same2D->Sumw2();
    TH2D* mix2D = mixed->Projection(0, 1);
    mix2D->Sumw2();

    same2D->RebinX();
    mix2D->RebinX();
    same2D->RebinY();
    mix2D->RebinY();

    same2D->GetXaxis()->SetRange(1,same2D->GetXaxis()->GetNbins());
    mix2D->GetXaxis()->SetRange(1,mix2D->GetXaxis()->GetNbins());

    Float_t scale = 0.0;

    scale = 0.5*(float)(mix2D->GetBinContent(mix2D->GetXaxis()->FindBin(0.01), mix2D->GetYaxis()->FindBin(0)) + mix2D->GetBinContent(mix2D->GetXaxis()->FindBin(-0.01), mix2D->GetYaxis()->FindBin(0.0)));
    printf("scale: %e \n", scale);
    same2D->Divide(mix2D);
    same2D->Scale(scale);

    same->GetAxis(3)->SetRange(0,0);
    mixed->GetAxis(3)->SetRange(0,0);
    return same2D;
}

makeSimpleMixCorrections(string inputName){
    TFile *histoFile = new TFile(inputName.c_str()); 
    TList* list_0_20 = (TList*) histoFile->Get("phiCorr_mult_0_20");
    TList* list_20_50 = (TList*) histoFile->Get("phiCorr_mult_20_50");
    TList* list_50_100 = (TList*) histoFile->Get("phiCorr_mult_50_100");
    //histoFile->cd("PhiReconstruction");
/*
    THnSparseF *fkkUSDist = (THnSparseF *)InvMass->FindObject("fkkUSDist");
    THnSparseF *fkkLSDist = (THnSparseF *)InvMass->FindObject("fkkLSDist");
    THnSparseF *fTrigDist = (THnSparseF *)InvMass->FindObject("fTrigDist");
    THnSparseF *dphiHPhi = (THnSparseF *)InvMass->FindObject("fDphiHPhi");
    THnSparseF *dphiHKK = (THnSparseF *)InvMass->FindObject("fDphiHKK");
    THnSparseF *dphiHPhiMixed = (THnSparseF *)InvMass->FindObject("fDphiHPhiMixed");
    THnSparseF *dphiHKKMixed = (THnSparseF *)InvMass->FindObject("fDphiHKKMixed");
    TH1D* zVtx = InvMass->FindObject("fVtxZ");
*/
    THnSparseF *dphiHPhi = (THnSparseF *)list_0_20->FindObject("fDphiHPhi");
    THnSparseF *dphiHKK = (THnSparseF *)list_0_20->FindObject("fDphiHKK");
    THnSparseF *dphiHPhiMixed = (THnSparseF *)list_0_20->FindObject("fDphiHPhiMixed");
    THnSparseF *dphiHKKMixed = (THnSparseF *)list_0_20->FindObject("fDphiHKKMixed");

    dphiHPhi->Add((THnSparseF *)list_20_50->FindObject("fDphiHPhi"));
    dphiHKK->Add((THnSparseF *)list_20_50->FindObject("fDphiHKK"));
    dphiHPhiMixed->Add((THnSparseF *)list_20_50->FindObject("fDphiHPhiMixed"));
    dphiHKKMixed->Add((THnSparseF *)list_20_50->FindObject("fDphiHKKMixed"));

    dphiHPhi->Add((THnSparseF *)list_50_100->FindObject("fDphiHPhi"));
    dphiHKK->Add((THnSparseF *)list_50_100->FindObject("fDphiHKK"));
    dphiHPhiMixed->Add((THnSparseF *)list_50_100->FindObject("fDphiHPhiMixed"));
    dphiHKKMixed->Add((THnSparseF *)list_50_100->FindObject("fDphiHKKMixed"));


    dphiHPhi->GetAxis(0)->SetRangeUser(4.0, 8.0); 
    dphiHPhi->GetAxis(1)->SetRangeUser(2.0,4.0); 
    dphiHKK->GetAxis(0)->SetRangeUser(4.0, 8.0);
    dphiHKK->GetAxis(1)->SetRangeUser(2.0,4.0);
   
    dphiHPhiMixed->GetAxis(0)->SetRangeUser(4.0, 8.0); 
    dphiHPhiMixed->GetAxis(1)->SetRangeUser(2.0,4.0); 
    dphiHKKMixed->GetAxis(0)->SetRangeUser(4.0, 8.0);
    dphiHKKMixed->GetAxis(1)->SetRangeUser(2.0,4.0);

    dphiHPhi->GetAxis(5)->SetRange(0,0);
    dphiHPhi->GetAxis(4)->SetRangeUser(-10.0, 10.0);

    Int_t axes[] = {2,3,4,5};

    THnSparseF* hPhi = dphiHPhi->Projection(4, axes);
    THnSparseF* hKK = dphiHKK->Projection(4, axes);
    THnSparseF* hPhiMixed = dphiHPhiMixed->Projection(4, axes);
    THnSparseF* hKKMixed = dphiHKKMixed->Projection(4, axes);



    TH2D* hPhi2Dpeak = makeCorrections(hPhi, hPhiMixed, 1.010, 1.030);
    hPhi2Dpeak->SetName("hPhi2Dpeak");
    TH2D* hKK2Dpeak = makeCorrections(hKK, hKKMixed, 1.010, 1.030);
    hKK2Dpeak->SetName("hKK2Dpeak");
    TH2D* hPhi2DRside = makeCorrections(hPhi, hPhiMixed, 1.040, 1.060);
    hPhi2DRside->SetName("hPhi2DRside");
    TH2D* hKK2DRside = makeCorrections(hKK, hKKMixed, 1.040, 1.060);
    hKK2DRside->SetName("hKK2DRside");
    TH2D* hPhi2DLside = makeCorrections(hPhi, hPhiMixed, 0.995, 1.005);
    hPhi2DLside->SetName("hPhi2DLside");
    TH2D* hKK2DLside = makeCorrections(hKK, hKKMixed, 0.995, 1.005);
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


    TFile* output = new TFile(Form("eta20_simplecorrected_%s", inputName.c_str()), "RECREATE");
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
}
