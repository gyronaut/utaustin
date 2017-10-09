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
makeMixCorrections(string inputName, float trigPTLow, float trigPTHigh, float assocPTLow, float assocPTHigh){
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

    TH2D *trigSameUSDist = (TH2D*)list->FindObject("fTrigSameUSDist");
    TH2D *trigSameLSDist = (TH2D*)list->FindObject("fTrigSameLSDist");

    float trigMixScalesUS[10] = {};
    float trigMixScalesLS[10] = {};
    for(int i = 0; i < 10; i++){
        trigMixScalesUS[i] = (float) trigSameUSDist->Integral(trigSameUSDist->GetXaxis()->FindBin(trigPTLow), trigSameUSDist->GetXaxis()->FindBin(trigPTHigh), i+1, i+1);
        trigMixScalesLS[i] = (float) trigSameLSDist->Integral(trigSameLSDist->GetXaxis()->FindBin(trigPTLow), trigSameLSDist->GetXaxis()->FindBin(trigPTHigh), i+1, i+1);
    }

    float totalTrigSameUS = (float)trigSameUSDist->Integral(trigSameUSDist->GetXaxis()->FindBin(trigPTLow), trigSameUSDist->GetXaxis()->FindBin(trigPTHigh), 1, trigSameUSDist->GetYaxis()->GetNbins());
    float totalTrigSameLS = (float)trigSameUSDist->Integral(trigSameLSDist->GetXaxis()->FindBin(trigPTLow), trigSameLSDist->GetXaxis()->FindBin(trigPTHigh), 1, trigSameLSDist->GetYaxis()->GetNbins());
   

    THnSparseF *dphiHPhi = (THnSparseF *)list->FindObject("fDphiHPhi");
    THnSparseF *dphiHKK = (THnSparseF *)list->FindObject("fDphiHKK");
    THnSparseF *dphiHPhiMixed = (THnSparseF *)list->FindObject("fDphiHPhiMixed");
    THnSparseF *dphiHKKMixed = (THnSparseF *)list->FindObject("fDphiHKKMixed");

    THnSparseF *dphiHH = (THnSparseF*)list->FindObject("fDphiHH");
    THnSparseF *dphiHHMixed = (THnSparseF*)list->FindObject("fDPhiHHMixed");

    //make 4D THnProjections projection to do mixed event corrections    
    dphiHPhi->GetAxis(0)->SetRangeUser(trigPTLow, trigPTHigh); 
    dphiHPhi->GetAxis(1)->SetRangeUser(assocPTLow,assocPTHigh); 
    dphiHKK->GetAxis(0)->SetRangeUser(trigPTLow, trigPTHigh);
    dphiHKK->GetAxis(1)->SetRangeUser(assocPTLow,assocPTHigh);
    dphiHH->GetAxis(0)->SetRangeUser(trigPTLow,trigPTHigh);
    dphiHH->GetAxis(1)->SetRangeUser(assocPTLow,assocPTHigh);
   
    dphiHPhiMixed->GetAxis(0)->SetRangeUser(trigPTLow, trigPTHigh); 
    dphiHPhiMixed->GetAxis(1)->SetRangeUser(assocPTLow,assocPTHigh); 
    dphiHKKMixed->GetAxis(0)->SetRangeUser(trigPTLow, trigPTHigh);
    dphiHKKMixed->GetAxis(1)->SetRangeUser(assocPTLow,assocPTHigh);
    dphiHHMixed->GetAxis(0)->SetRangeUser(trigPTLow,trigPTHigh);
    dphiHHMixed->GetAxis(1)->SetRangeUser(assocPTLow,assocPTHigh);

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
    TH3D* hh = dphiHH->Projection(2,3,4);
    hh->Sumw2();
    THnSparseF* hPhiMixed = dphiHPhiMixed->Projection(4, axes);
    THnSparseF* hKKMixed = dphiHKKMixed->Projection(4, axes);
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

    TH2D* hPhi2Dpeak = makeCorrections(hPhi, hPhiMixed, 1.010, 1.030, sameUSPeakEta, mixedUSPeakEta, trigMixScalesUS, totalTrigSameUS);
    hPhi2Dpeak->SetName("hPhi2Dpeak");
    for(int i = 0; i<10; i++){
        sameUSPeakEta[i]->SetName(Form("sameUSPeakEta_zvtx_%i", i));
        mixedUSPeakEta[i]->SetName(Form("mixedUSPeakEta_zvtx_%i", i));
    }
    TH2D* hKK2Dpeak = makeCorrections(hKK, hKKMixed, 1.010, 1.030, sameLSPeakEta, mixedLSPeakEta, trigMixScalesLS, totalTrigSameLS);
    hKK2Dpeak->SetName("hKK2Dpeak");
    for(int i = 0; i<10; i++){
        sameLSPeakEta[i]->SetName(Form("sameLSPeakEta_zvtx_%i", i));
        mixedLSPeakEta[i]->SetName(Form("mixedLSPeakEta_zvtx_%i", i));
    } 
    TH2D* hPhi2DRside = makeCorrections(hPhi, hPhiMixed, 1.040, 1.060, sameUSRsideEta, mixedUSRsideEta, trigMixScalesUS, totalTrigSameUS);
    hPhi2DRside->SetName("hPhi2DRside");
    for(int i = 0; i<10; i++){
        sameUSRsideEta[i]->SetName(Form("sameUSRsideEta_zvtx_%i", i));
        mixedUSRsideEta[i]->SetName(Form("mixedUSRsideEta_zvtx_%i", i));
    }
    TH2D* hKK2DRside = makeCorrections(hKK, hKKMixed, 1.040, 1.060, sameLSRsideEta, mixedLSRsideEta, trigMixScalesLS, totalTrigSameLS);
    hKK2DRside->SetName("hKK2DRside");
    for(int i = 0; i<10; i++){
        sameLSRsideEta[i]->SetName(Form("sameLSRsideEta_zvtx_%i", i));
        mixedLSRsideEta[i]->SetName(Form("mixedLSRsideEta_zvtx_%i", i));
    }
    TH2D* hPhi2DLside = makeCorrections(hPhi, hPhiMixed, 0.995, 1.005, sameUSLsideEta, mixedUSLsideEta, trigMixScalesUS, totalTrigSameUS);
    hPhi2DLside->SetName("hPhi2DLside");
    for(int i = 0; i<10; i++){
        sameUSLsideEta[i]->SetName(Form("sameUSLsideEta_zvtx_%i", i));
        mixedUSLsideEta[i]->SetName(Form("mixedUSLsideEta_zvtx_%i", i));
    }
    TH2D* hKK2DLside = makeCorrections(hKK, hKKMixed, 0.995, 1.005, sameLSLsideEta, mixedLSLsideEta, trigMixScalesLS, totalTrigSameLS);
    hKK2DLside->SetName("hKK2DLside");
    for(int i = 0; i<10; i++){
        sameLSLsideEta[i]->SetName(Form("sameLSLsideEta_zvtx_%i", i));
        mixedLSLsideEta[i]->SetName(Form("mixedLSLsideEta_zvtx_%i", i));
    }
    /*
    TH2D* hh2D = makehhCorrections(hh, hhMixed);
    hh2D->SetName("hh2D");

    hh->GetZaxis()->SetRangeUser(-10.0, 10.0);
    TH2D* uncorrhh2D = hh->Project3D("xye");
    uncorrhh2D->SetName("uncorrhh2D");
    */
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

    //make ratio plot of just 1 zvtx bin as a check:
    hPhiMixed->GetAxis(3)->SetRangeUser(1.01, 1.03);
    hPhiMixed->GetAxis(2)->SetRange(6,6);
    TH2D* mixedratioPeakZ2 = hPhiMixed->Projection(0,1);
    hKKMixed->GetAxis(3)->SetRangeUser(1.01, 1.03);
    hKKMixed->GetAxis(2)->SetRange(6,6);
    TH2D* hist = hKKMixed->Projection(0,1);
    mixedratioPeakZ2->Divide(hist);
    TH1D* mixedratioPeakZ2deta = mixedratioPeakZ2->ProjectionX("mixedratioPeakZ2deta");


    //project just zvtx distributions for hPhiMixed points and hKKMixed points in peak region:
    hPhiMixed->GetAxis(2)->SetRange(0,0);
    hKKMixed->GetAxis(2)->SetRange(0,0);
    TH1D* mixedUSzvtx = hPhiMixed->Projection(2);
    mixedUSzvtx->SetName("mixedUSzvtx");
    TH1D* mixedLSzvtx = hKKMixed->Projection(2);
    mixedLSzvtx->SetName("mixedLSzvtx");

    TFile* output = new TFile(Form("trig_%i_%i_assoc_%i_%i_mixcorr_%s", (int)trigPTLow, (int)trigPTHigh, (int)assocPTLow, (int)assocPTHigh, inputName.c_str()), "RECREATE");
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

    mixedratioPeakZ2->Write();
    mixedratioPeakZ2deta->Write();

    mixedUSzvtx->Write();
    mixedLSzvtx->Write();

    for(int i=0; i< 10; i++){
        sameUSPeakEta[i]->Write();
        mixedUSPeakEta[i]->Write();
        sameLSPeakEta[i]->Write();
        mixedLSPeakEta[i]->Write();
        sameUSRsideEta[i]->Write();
        mixedUSRsideEta[i]->Write();
        sameLSRsideEta[i]->Write();
        mixedLSRsideEta[i]->Write();
        sameUSLsideEta[i]->Write();
        mixedUSLsideEta[i]->Write();
        sameLSLsideEta[i]->Write();
        mixedLSLsideEta[i]->Write();
    }

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
