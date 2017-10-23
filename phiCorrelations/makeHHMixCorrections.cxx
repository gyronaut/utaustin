TH2D* makeCorrections(THnSparse* same, THnSparse* mixed, Float_t lowmass, Float_t highmass, TH1D** sameEta, TH1D** mixedEta, float* trigMixScale, float totalTrigSame){
    same->GetAxis(3)->SetRangeUser(lowmass, highmass);
    mixed->GetAxis(3)->SetRangeUser(lowmass, highmass);
    TH3D* same3D = same->Projection(0, 1, 2);
    same3D->Sumw2();
    TH3D* mix3D = mixed->Projection(0, 1, 2);
    mix3D->Sumw2();

    same3D->RebinX(4);
    mix3D->RebinX(4);
    same3D->RebinY(4);
    mix3D->RebinY(4);

    //same3D->GetYaxis()->SetRangeUser(-1.2, 1.2);
    //mix3D->GetYaxis()->SetRangeUser(-1.2, 1.2);
    same3D->GetYaxis()->SetRange(1,same3D->GetYaxis()->GetNbins());
    mix3D->GetYaxis()->SetRange(1, mix3D->GetYaxis()->GetNbins());

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

        //get d-eta 1D plots for same and mixed distributions for each zbin
        sameEta[zbin] = same2D[zbin]->ProjectionX(Form("sameEta_zvtx_%i", zbin), 1, same2D[zbin]->GetYaxis()->GetNbins());
        mixedEta[zbin] = mix2D[zbin]->ProjectionX(Form("mixedEta_zvtx_%i", zbin), 1, mix2D[zbin]->GetYaxis()->GetNbins());

        //scale mixed by number of triggers in zvtx bin
        //mix2D[zbin]->Scale(1.0/trigMixScale[zbin]);

        scale = 0.5*(float)(mix2D[zbin]->GetBinContent(mix2D[zbin]->GetXaxis()->FindBin(0.01), mix2D[zbin]->GetYaxis()->FindBin(0.0)) + mix2D[zbin]->GetBinContent(mix2D[zbin]->GetXaxis()->FindBin(-0.01), mix2D[zbin]->GetYaxis()->FindBin(0.0)));
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
    same2DTotal->Scale(1.0/totalTrigSame);
    return same2DTotal;
}

//---------------------------------------------------------------------------------------------
TH2D* makehhCorrections(TH3D* same3D, TH3D* mix3D){
    same3D->RebinX();
    mix3D->RebinX();
    same3D->RebinY();
    mix3D->RebinY();

    //same3D->GetYaxis()->SetRangeUser(-1.2, 1.2);
    //mix3D->GetYaxis()->SetRangeUser(-1.2, 1.2);
    same3D->GetYaxis()->SetRange(1,same3D->GetYaxis()->GetNbins());
    mix3D->GetYaxis()->SetRange(1, mix3D->GetYaxis()->GetNbins());

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

        scale = 0.5*(float)(mix2D[zbin]->GetBinContent(mix2D[zbin]->GetXaxis()->FindBin(0.01), mix2D[zbin]->GetYaxis()->FindBin(0.0)) + mix2D[zbin]->GetBinContent(mix2D[zbin]->GetXaxis()->FindBin(-0.01), mix2D[zbin]->GetYaxis()->FindBin(0.0)));
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

    return same2DTotal;
}

//--------------------------------------------------------------------------------------------
makeHHMixCorrections(string inputName, float trigPTLow, float trigPTHigh, float assocPTLow, float assocPTHigh){
    TFile *histoFile = new TFile(inputName.c_str());
    string mult = inputName.substr(inputName.find("_", inputName.find("_")+1), inputName.find(".") - inputName.find("_", inputName.find("_")+1));
    TList* list = (TList*) histoFile->Get(Form("phiCorr_mult%s", mult.c_str()));
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

    TH2D *trigHHDist = (TH2D*)list->FindObject("fTrigHHDist");

    float trigMixScales[10] = {};
    for(int i = 0; i < 10; i++){
        trigMixScales[i] = (float) trigHHDist->Integral(trigHHDist->GetXaxis()->FindBin(trigPTLow), trigHHDist->GetXaxis()->FindBin(trigPTHigh), i+1, i+1);
    }

    float totalTrig = (float)trigHHDist->Integral(trigHHDist->GetXaxis()->FindBin(trigPTLow), trigHHDist->GetXaxis()->FindBin(trigPTHigh), 1, trigHHDist->GetYaxis()->GetNbins());
      
    THnSparseF *dphiHH = (THnSparseF*)list->FindObject("fDphiHH");
    THnSparseF *dphiHHMixed = (THnSparseF*)list->FindObject("fDPhiHHMixed");

    //make 4D THnProjections projection to do mixed event corrections    
    dphiHH->GetAxis(0)->SetRangeUser(trigPTLow, trigPTHigh); 
    dphiHH->GetAxis(1)->SetRangeUser(assocPTLow,assocPTHigh); 
   
    dphiHHMixed->GetAxis(0)->SetRangeUser(trigPTLow, trigPTHigh); 
    dphiHHMixed->GetAxis(1)->SetRangeUser(assocPTLow,assocPTHigh); 

    dphiHH->GetAxis(4)->SetRangeUser(-10.0, 10.0);
    dphiHHMixed->GetAxis(4)->SetRangeUser(-10.0, 10.0);

    Int_t axes[] = {2,3,4,5};

    TH3D* hh = dphiHH->Projection(2,3,4);
    hh->Sumw2();
    TH3D* hhMixed = dphiHHMixed->Projection(2,3,4);
    hhMixed->Sumw2();

    TH1D* sameUSPeakEta[10];
    TH1D* mixedUSPeakEta[10];
    TH1D* sameLSPeakEta[10];
    TH1D* mixedLSPeakEta[10];
    TH1D* sameUSRsideEta[10];
    TH1D* mixedUSRsideEta[10];
    TH1D* sameLSRsideEta[10];
    TH1D* mixedLSRsideEta[10];
    TH1D* sameUSLsideEta[10];
    TH1D* mixedUSLsideEta[10];
    TH1D* sameLSLsideEta[10];
    TH1D* mixedLSLsideEta[10];

    
    TH2D* hh2D = makehhCorrections(hh, hhMixed);
    hh2D->SetName("hh2D");
    hh2D->Scale(1.0/(hh2D->Integral(hh2D->GetXaxis()->FindBin(-1.2), hh2D->GetXaxis()->FindBin(1.2), 1, hh2D->GetYaxis()->GetNbins())));

    TH1D* hhdphi = hh2D->ProjectionY("hhdphi", hh2D->GetXaxis()->FindBin(-1.2), hh2D->GetXaxis()->FindBin(1.2));
    hhdphi->Scale(1.0/(hhdphi->Integral()));

    hh->GetZaxis()->SetRangeUser(-10.0, 10.0);
    TH2D* uncorrhh2D = hh->Project3D("xye");
    uncorrhh2D->SetName("uncorrhh2D");
    uncorrhh2D->Scale(1.0/(uncorrhh2D->Integral(uncorrhh2D->GetXaxis()->FindBin(-1.2), uncorrhh2D->GetXaxis()->FindBin(1.2), 1, uncorrhh2D->GetYaxis()->GetNbins())));
    
    TFile* output = new TFile(Form("trig_%i_%i_assoc_%i_%i_mixcorr_%s", (int)trigPTLow, (int)trigPTHigh, (int)assocPTLow, (int)assocPTHigh, inputName.c_str()), "RECREATE");
    hh2D->Write();
    hhdphi->Write();
    uncorrhh2D->Write(); 
}
