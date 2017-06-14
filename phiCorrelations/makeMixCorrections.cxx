TH2D* makeCorrections(THnSparse* same, THnSparse* mixed, Float_t lowmass, Float_t highmass){
    same->GetAxis(3)->SetRangeUser(lowmass, highmass);
    mixed->GetAxis(3)->SetRangeUser(lowmass, highmass);
    TH3D* same3D = same->Projection(0, 1, 2);
    TH3D* mix3D = mixed->Projection(0, 1, 2);

    same3D->RebinX();
    mix3D->RebinX();
    same3D->RebinY();
    mix3D->RebinY();

    same3D->GetYaxis()->SetRangeUser(-1.2, 1.2);
    mix3D->GetYaxis()->SetRangeUser(-1.2, 1.2);

    Float_t scale = 0.0;
    TH2D* same2D[10];
    TH2D* mix2D[10];
    TH2D* same2DTotal;

    for(int zbin = 0; zbin < 10; zbin++){
        same3D->GetZaxis()->SetRange(zbin+1, zbin+1);
        same2D[zbin] = (TH2D*)same3D->Project3D("xye");
        same2D[zbin]->SetName(Form("2dproj_zbin%i", zbin));
        
        mix3D->GetZaxis()->SetRange(zbin+1, zbin+1);
        mix2D[zbin] = (TH2D*)mix3D->Project3D("xye");
        mix2D[zbin]->SetName(Form("mix2droj_zbin%i", zbin));

        scale = 0.5*(float)(mix2D[zbin]->GetBinContent(mix2D[zbin]->GetXaxis()->FindBin(0), mix2D[zbin]->GetYaxis()->FindBin(0.01)) + mix2D[zbin]->GetBinContent(mix2D[zbin]->GetXaxis()->FindBin(0), mix2D[zbin]->GetYaxis()->FindBin(-0.01)));
        printf("scale: %e \n", scale);
        same2D[zbin]->Divide(mix2D[zbin]);
        same2D[zbin]->Scale(scale);
        if(zbin==0){
            same2DTotal = (TH2D*)same2D[zbin]->Clone("2dproj_total");
        }else{
            same2DTotal->Add(same2D[zbin]);
        }
        
        same3D->GetZaxis()->SetRange(0,0);
        mix3D->GetZaxis()->SetRange(0,0);
    }

    same->GetAxis(3)->SetRange(0,0);
    mixed->GetAxis(3)->SetRange(0,0);
    return same2DTotal;
}

makeMixCorrections(string inputName){
    TFile *histoFile = new TFile(inputName.c_str());
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
    THnSparseF *dphiHPhi = (THnSparseF *)phiCorr_mult_50_100->FindObject("fDphiHPhi");
    THnSparseF *dphiHKK = (THnSparseF *)phiCorr_mult_50_100->FindObject("fDphiHKK");
    THnSparseF *dphiHPhiMixed = (THnSparseF *)phiCorr_mult_50_100->FindObject("fDphiHPhiMixed");
    THnSparseF *dphiHKKMixed = (THnSparseF *)phiCorr_mult_50_100->FindObject("fDphiHKKMixed");

    //make 4D THnProjections projection to do mixed event corrections
    
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

    TFile* output = new TFile(Form("corrected_%s", inputName.c_str()), "RECREATE");
    hPhi2Dpeak->Write();
    hKK2Dpeak->Write();
    hPhi2DRside->Write();
    hKK2DRside->Write();
    hPhi2DLside->Write();
    hKK2DLside->Write();

/*    
    hPhi->Write();
    hPhiMixed->Write();

    hPhiMixed->GetAxis(3)->SetRangeUser(1.01, 1.03);
    TH3D* mix3D = hPhiMixed->Projection(0, 1, 2);
    mix3D->Write();

    hPhi->GetAxis(3)->SetRangeUser(1.01, 1.03);
    TH3D* same3D = hPhi->Projection(0,1,2);
    same3D->Write();

    same3D->GetZaxis()->SetRange(1,1);
    TH2D* same2D = (TH2D*)same3D->Project3D("xye");
    same2D->Write();

    TH1D* dphi = hPhi->Projection(1);
    dphi->Write();
*/
}
