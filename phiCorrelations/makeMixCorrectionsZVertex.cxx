TH2D* makeCorrections(THnSparse* same, THnSparse* mixed, Float_t lowmass, Float_t highmass, TH1D* sameEta, TH1D* mixedEta, Int_t zbin){
    same->GetAxis(2)->SetRangeUser(lowmass, highmass);
    mixed->GetAxis(2)->SetRangeUser(lowmass, highmass);
    TH2D* same2D = same->Projection(0, 1);
    same2D->Sumw2();
    TH2D* mix2D = mixed->Projection(0, 1);
    mix2D->Sumw2();

    //same2D->RebinX(4);
    //mix2D->RebinX(4);
    //same2D->RebinY(4);
    //mix2D->RebinY(4);

    //same3D->GetYaxis()->SetRangeUser(-1.2, 1.2);
    //mix3D->GetYaxis()->SetRangeUser(-1.2, 1.2);
    same2D->GetYaxis()->SetRange(1,same2D->GetYaxis()->GetNbins());
    mix2D->GetYaxis()->SetRange(1, mix2D->GetYaxis()->GetNbins());

    Float_t scale = 0.0;
    TH2D* same2DTotal;


    //get d-eta 1D plots for same and mixed distributions for each zbin
    sameEta = same2D->ProjectionX(Form("sameEta_zvtx_%i", zbin), 1, same2D->GetYaxis()->GetNbins());
    mixedEta = mix2D->ProjectionX(Form("mixedEta_zvtx_%i", zbin), 1, mix2D->GetYaxis()->GetNbins());


    scale = 0.5*(float)(mix2D->GetBinContent(mix2D->GetXaxis()->FindBin(0.01), mix2D->GetYaxis()->FindBin(0.0)) + mix2D->GetBinContent(mix2D->GetXaxis()->FindBin(-0.01), mix2D->GetYaxis()->FindBin(0.0)));
    printf("scale: %e \n", scale);
    same2D->Divide(mix2D);
    same2D->Scale(scale);

    same->GetAxis(2)->SetRange(0,0);
    mixed->GetAxis(2)->SetRange(0,0);
    return same2D;
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
void makeMixCorrectionsZVertex(string inputName, int multLow, int multHigh, float trigPTLow, float trigPTHigh, float assocPTLow, float assocPTHigh){
    TFile *histoFile = new TFile(inputName.c_str());
    //string mult = inputName.substr(inputName.find("_", inputName.find("_")+1), inputName.find(".") - inputName.find("_", inputName.find("_")+1));
    string mult = "_" + std::to_string(multLow) + "_" + std::to_string(multHigh);
    TList* list = (TList*) histoFile->Get(Form("phiCorr_mult%s", mult.c_str()));

    TH2D *trigSameUSDist = (TH2D*)list->FindObject("fTrigSameUSDist");
    TH2D *trigSameLSDist = (TH2D*)list->FindObject("fTrigSameLSDist");

    float totalTrigSameUS = (float)trigSameUSDist->Integral(trigSameUSDist->GetXaxis()->FindBin(trigPTLow), trigSameUSDist->GetXaxis()->FindBin(trigPTHigh), 1, trigSameUSDist->GetYaxis()->GetNbins());
    float totalTrigSameLS = (float)trigSameLSDist->Integral(trigSameLSDist->GetXaxis()->FindBin(trigPTLow), trigSameLSDist->GetXaxis()->FindBin(trigPTHigh), 1, trigSameLSDist->GetYaxis()->GetNbins());
   
    //TH1D* vtxZmixbins = (TH1D*)list->FindObject("fVtxZmixbins");
    //Int_t numbinsZvtx = vtxZmixbins->GetXaxis()->GetNbins();
 
    Int_t numbinsZvtx = 8;

    TH1D* sameUSPeakEta[numbinsZvtx];
    TH1D* mixedUSPeakEta[numbinsZvtx];
    TH1D* sameLSPeakEta[numbinsZvtx];
    TH1D* mixedLSPeakEta[numbinsZvtx];
    TH1D* sameUSRsideEta[numbinsZvtx];
    TH1D* mixedUSRsideEta[numbinsZvtx];
    TH1D* sameLSRsideEta[numbinsZvtx];
    TH1D* mixedLSRsideEta[numbinsZvtx];
    TH1D* sameUSLsideEta[numbinsZvtx];
    TH1D* mixedUSLsideEta[numbinsZvtx];
    TH1D* sameLSLsideEta[numbinsZvtx];
    TH1D* mixedLSLsideEta[numbinsZvtx];

    THnSparseF *dphiHPhi[numbinsZvtx];
    THnSparseF *dphiHKK[numbinsZvtx];
    THnSparseF *dphiHPhiMixed[numbinsZvtx];
    THnSparseF *dphiHKKMixed[numbinsZvtx];

    THnSparseF *hPhi[numbinsZvtx];
    THnSparseF *hKK[numbinsZvtx];
    THnSparseF *hPhiMixed[numbinsZvtx];
    THnSparseF *hKKMixed[numbinsZvtx];

    TH2D* hPhi2Dpeak[numbinsZvtx];
    TH2D* hKK2Dpeak[numbinsZvtx];
    TH2D* hPhi2DRside[numbinsZvtx];
    TH2D* hKK2DRside[numbinsZvtx];
    TH2D* hPhi2DLside[numbinsZvtx];
    TH2D* hKK2DLside[numbinsZvtx];

    TH2D* hPhi2DpeakTotal;
    TH2D* hKK2DpeakTotal;
    TH2D* hPhi2DRsideTotal;
    TH2D* hKK2DRsideTotal;
    TH2D* hPhi2DLsideTotal;
    TH2D* hKK2DLsideTotal;

    THnSparseF* hPhiTotal;
    THnSparseF* hPhiMixedTotal;
    THnSparseF* hKKTotal;
    THnSparseF* hKKMixedTotal;

    for(int izvtx = 0; izvtx < numbinsZvtx; izvtx++){
        dphiHPhi[izvtx] = (THnSparseF *)list->FindObject(Form("fDphiHPhiz%i", izvtx));
        dphiHKK[izvtx] = (THnSparseF *)list->FindObject(Form("fDphiHKKz%i", izvtx));
        dphiHPhiMixed[izvtx] = (THnSparseF *)list->FindObject(Form("fDphiHPhiMixedz%i", izvtx));
        dphiHKKMixed[izvtx] = (THnSparseF *)list->FindObject(Form("fDphiHKKMixedz%i", izvtx));


        //make 4D THnProjections projection to do mixed event corrections    
        dphiHPhi[izvtx]->GetAxis(0)->SetRangeUser(trigPTLow, trigPTHigh); 
        dphiHPhi[izvtx]->GetAxis(1)->SetRangeUser(assocPTLow,assocPTHigh); 
        dphiHKK[izvtx]->GetAxis(0)->SetRangeUser(trigPTLow, trigPTHigh);
        dphiHKK[izvtx]->GetAxis(1)->SetRangeUser(assocPTLow,assocPTHigh);

        dphiHPhiMixed[izvtx]->GetAxis(0)->SetRangeUser(trigPTLow, trigPTHigh); 
        dphiHPhiMixed[izvtx]->GetAxis(1)->SetRangeUser(assocPTLow,assocPTHigh); 
        dphiHKKMixed[izvtx]->GetAxis(0)->SetRangeUser(trigPTLow, trigPTHigh);
        dphiHKKMixed[izvtx]->GetAxis(1)->SetRangeUser(assocPTLow,assocPTHigh);

        dphiHPhi[izvtx]->GetAxis(4)->SetRange(1,dphiHPhi[izvtx]->GetAxis(4)->GetNbins()); 
        dphiHKK[izvtx]->GetAxis(4)->SetRange(1,dphiHKK[izvtx]->GetAxis(4)->GetNbins());
        dphiHPhiMixed[izvtx]->GetAxis(4)->SetRange(1,dphiHPhiMixed[izvtx]->GetAxis(4)->GetNbins());
        dphiHKKMixed[izvtx]->GetAxis(4)->SetRange(1,dphiHKKMixed[izvtx]->GetAxis(4)->GetNbins());

        Int_t axes[] = {2,3,4};

        hPhi[izvtx] = (THnSparseF*)dphiHPhi[izvtx]->Projection(3, axes);
        hKK[izvtx] = (THnSparseF*)dphiHKK[izvtx]->Projection(3, axes);
        hPhiMixed[izvtx] = (THnSparseF*)dphiHPhiMixed[izvtx]->Projection(3, axes);
        hKKMixed[izvtx] = (THnSparseF*)dphiHKKMixed[izvtx]->Projection(3, axes);

        hPhi2Dpeak[izvtx] = makeCorrections(hPhi[izvtx], hPhiMixed[izvtx], 1.010, 1.030, sameUSPeakEta[izvtx], mixedUSPeakEta[izvtx], izvtx);
        hPhi2Dpeak[izvtx]->SetName(Form("hPhi2Dpeakz%i", izvtx));
        //sameUSPeakEta[izvtx]->SetName(Form("sameUSPeakEta_zvtx_%i", izvtx));
        //mixedUSPeakEta[izvtx]->SetName(Form("mixedUSPeakEta_zvtx_%i", izvtx));

        hKK2Dpeak[izvtx] = makeCorrections(hKK[izvtx], hKKMixed[izvtx], 1.010, 1.030, sameLSPeakEta[izvtx], mixedLSPeakEta[izvtx], izvtx);
        hKK2Dpeak[izvtx]->SetName(Form("hKK2Dpeakz%i", izvtx));
        //sameLSPeakEta[izvtx]->SetName(Form("sameLSPeakEta_zvtx_%i", izvtx));
        //mixedLSPeakEta[izvtx]->SetName(Form("mixedLSPeakEta_zvtx_%i", izvtx));

        hPhi2DRside[izvtx] = makeCorrections(hPhi[izvtx], hPhiMixed[izvtx], 1.040, 1.060, sameUSRsideEta[izvtx], mixedUSRsideEta[izvtx], izvtx);
        hPhi2DRside[izvtx]->SetName(Form("hPhi2DRsidez%i", izvtx));
        //sameUSRsideEta[izvtx]->SetName(Form("sameUSRsideEta_zvtx_%i", izvtx));
        //mixedUSRsideEta[izvtx]->SetName(Form("mixedUSRsideEta_zvtx_%i", izvtx));

        hKK2DRside[izvtx] = makeCorrections(hKK[izvtx], hKKMixed[izvtx], 1.040, 1.060, sameLSRsideEta[izvtx], mixedLSRsideEta[izvtx], izvtx);
        hKK2DRside[izvtx]->SetName(Form("hKK2DRsidez%i", izvtx));
        //sameLSRsideEta[izvtx]->SetName(Form("sameLSRsideEta_zvtx_%i", izvtx));
        //mixedLSRsideEta[izvtx]->SetName(Form("mixedLSRsideEta_zvtx_%i", izvtx));

        hPhi2DLside[izvtx] = makeCorrections(hPhi[izvtx], hPhiMixed[izvtx], 0.995, 1.005, sameUSLsideEta[izvtx], mixedUSLsideEta[izvtx], izvtx);
        hPhi2DLside[izvtx]->SetName(Form("hPhi2DLsidez%i", izvtx));
        //sameUSLsideEta[izvtx]->SetName(Form("sameUSLsideEta_zvtx_%i", izvtx));
        //mixedUSLsideEta[izvtx]->SetName(Form("mixedUSLsideEta_zvtx_%i", izvtx));

        hKK2DLside[izvtx] = makeCorrections(hKK[izvtx], hKKMixed[izvtx], 0.995, 1.005, sameLSLsideEta[izvtx], mixedLSLsideEta[izvtx], izvtx);
        hKK2DLside[izvtx]->SetName(Form("hKK2DLsidez%i", izvtx));
        //sameLSLsideEta[izvtx]->SetName(Form("sameLSLsideEta_zvtx_%i", izvtx));
        //mixedLSLsideEta[izvtx]->SetName(Form("mixedLSLsideEta_zvtx_%i", izvtx));

        if(izvtx==0){
            hPhi2DpeakTotal = (TH2D*)hPhi2Dpeak[izvtx]->Clone("hPhi2Dpeak");
            hPhi2DLsideTotal = (TH2D*)hPhi2DLside[izvtx]->Clone("hPhi2DLside");
            hPhi2DRsideTotal = (TH2D*)hPhi2DRside[izvtx]->Clone("hPhi2DRside");
            hKK2DpeakTotal = (TH2D*)hKK2Dpeak[izvtx]->Clone("hKK2Dpeak");
            hKK2DLsideTotal = (TH2D*)hKK2DLside[izvtx]->Clone("hKK2DLside");
            hKK2DRsideTotal = (TH2D*)hKK2DRside[izvtx]->Clone("hKK2DRside");

            hPhiTotal = (THnSparseF*)hPhi[izvtx]->Clone("hPhiTotal");
            hPhiMixedTotal = (THnSparseF*)hPhiMixed[izvtx]->Clone("hPhiMixedTotal");
            hKKTotal = (THnSparseF*)hKK[izvtx]->Clone("hKKTotal");
            hKKMixedTotal = (THnSparseF*)hKKMixed[izvtx]->Clone("hKKMixedTotal");
       }else{
            hPhi2DpeakTotal->Add(hPhi2Dpeak[izvtx]);
            hPhi2DRsideTotal->Add(hPhi2DRside[izvtx]);
            hPhi2DLsideTotal->Add(hPhi2DLside[izvtx]);
            hKK2DpeakTotal->Add(hKK2Dpeak[izvtx]);
            hKK2DRsideTotal->Add(hKK2DRside[izvtx]);
            hKK2DLsideTotal->Add(hKK2DLside[izvtx]);

            hPhiTotal->Add(hPhi[izvtx]);
            hPhiMixedTotal->Add(hPhiMixed[izvtx]);
            hKKTotal->Add(hKK[izvtx]);
            hKKMixedTotal->Add(hKKMixed[izvtx]);
       }
    }

    hPhi2DpeakTotal->Scale(1.0/totalTrigSameUS);
    hPhi2DRsideTotal->Scale(1.0/totalTrigSameUS);
    hPhi2DLsideTotal->Scale(1.0/totalTrigSameUS);
    hKK2DpeakTotal->Scale(1.0/totalTrigSameLS);
    hKK2DRsideTotal->Scale(1.0/totalTrigSameLS);
    hKK2DLsideTotal->Scale(1.0/totalTrigSameLS);

    //Create some uncorrected same/mixed event 2D histos
    hPhiTotal->GetAxis(2)->SetRangeUser(1.010, 1.030);
    hPhiMixedTotal->GetAxis(2)->SetRangeUser(1.010, 1.030);
    hKKTotal->GetAxis(2)->SetRangeUser(1.010, 1.030);
    hKKMixedTotal->GetAxis(2)->SetRangeUser(1.010, 1.030);
    TH2D* uncorrhPhi2Dpeak = hPhiTotal->Projection(0,1);
    uncorrhPhi2Dpeak->Sumw2();
    uncorrhPhi2Dpeak->SetName("uncorrhPhi2Dpeak");
    TH2D* uncorrhKK2Dpeak = hKKTotal->Projection(0,1);
    uncorrhKK2Dpeak->Sumw2();
    uncorrhKK2Dpeak->SetName("uncorrhKK2Dpeak");
    TH2D* uncorrhPhiMixed2Dpeak = hPhiMixedTotal->Projection(0,1);
    uncorrhPhiMixed2Dpeak->Sumw2();
    uncorrhPhiMixed2Dpeak->SetName("uncorrhPhiMixed2Dpeak");
    TH2D* uncorrhKKMixed2Dpeak = hKKMixedTotal->Projection(0,1);
    uncorrhKKMixed2Dpeak->Sumw2();
    uncorrhKKMixed2Dpeak->SetName("uncorrhKKMixed2Dpeak");

    hPhiTotal->GetAxis(2)->SetRangeUser(1.040, 1.060);
    hPhiMixedTotal->GetAxis(2)->SetRangeUser(1.040, 1.060);
    hKKTotal->GetAxis(2)->SetRangeUser(1.040, 1.060);
    hKKMixedTotal->GetAxis(2)->SetRangeUser(1.040, 1.060);
    TH2D* uncorrhPhi2DRside = hPhiTotal->Projection(0,1);
    uncorrhPhi2DRside->Sumw2();
    uncorrhPhi2DRside->SetName("uncorrhPhi2DRside");
    TH2D* uncorrhKK2DRside = hKKTotal->Projection(0,1);
    uncorrhKK2DRside->Sumw2();
    uncorrhKK2DRside->SetName("uncorrhKK2DRside");
    TH2D* uncorrhPhiMixed2DRside = hPhiMixedTotal->Projection(0,1);
    uncorrhPhiMixed2DRside->Sumw2();
    uncorrhPhiMixed2DRside->SetName("uncorrhPhiMixed2DRside");
    TH2D* uncorrhKKMixed2DRside = hKKMixedTotal->Projection(0,1);
    uncorrhKKMixed2DRside->Sumw2();
    uncorrhKKMixed2DRside->SetName("uncorrhKKMixed2DRside");

    hPhiTotal->GetAxis(2)->SetRangeUser(0.995, 1.005);
    hPhiMixedTotal->GetAxis(2)->SetRangeUser(0.995, 1.005);
    hKKTotal->GetAxis(2)->SetRangeUser(0.995, 1.005);
    hKKMixedTotal->GetAxis(2)->SetRangeUser(0.995, 1.005);
    TH2D* uncorrhPhi2DLside = hPhiTotal->Projection(0,1);
    uncorrhPhi2DLside->Sumw2();
    uncorrhPhi2DLside->SetName("uncorrhPhi2DLside");
    TH2D* uncorrhKK2DLside = hKKTotal->Projection(0,1);
    uncorrhKK2DLside->Sumw2();
    uncorrhKK2DLside->SetName("uncorrhKK2DLside");
    TH2D* uncorrhPhiMixed2DLside = hPhiMixedTotal->Projection(0,1);
    uncorrhPhiMixed2DLside->Sumw2();
    uncorrhPhiMixed2DLside->SetName("uncorrhPhiMixed2DLside");
    TH2D* uncorrhKKMixed2DLside = hKKMixedTotal->Projection(0,1);
    uncorrhKKMixed2DLside->Sumw2();
    uncorrhKKMixed2DLside->SetName("uncorrhKKMixed2DLside");

    //Make some ratio plots of mixed event US over LS for sidbeand and peak regions
    //also make projections onto delta eta and delta phi (on restricted delta eta range)
    TH2D* mixedratioRSB = (TH2D*)uncorrhPhiMixed2DRside->Clone("mixedratioRSB");
    mixedratioRSB->Divide(uncorrhKKMixed2DRside);
    TH1D* mixedratioRSBdeta = mixedratioRSB->ProjectionX("mixedratioRSBdeta", 1, mixedratioRSB->GetYaxis()->GetNbins());
    TH1D* mixedratioRSBdphi = mixedratioRSB->ProjectionY("mixedratioRSBdphi", mixedratioRSB->GetXaxis()->FindBin(-1.2), mixedratioRSB->GetXaxis()->FindBin(1.2));
    mixedratioRSB->GetXaxis()->SetRange(0,0);

    TH2D* mixedratioLSB = (TH2D*)uncorrhPhiMixed2DLside->Clone("mixedratioLSB");
    mixedratioLSB->Divide(uncorrhKKMixed2DLside);
    TH1D* mixedratioLSBdeta = mixedratioLSB->ProjectionX("mixedratioLSBdeta", 1, mixedratioLSB->GetYaxis()->GetNbins());
    TH1D* mixedratioLSBdphi = mixedratioLSB->ProjectionY("mixedratioLSBdphi", mixedratioLSB->GetXaxis()->FindBin(-1.2), mixedratioLSB->GetXaxis()->FindBin(1.2));
    mixedratioLSB->GetXaxis()->SetRange(0,0);

    TH2D* mixedratioPeak = (TH2D*)uncorrhPhiMixed2Dpeak->Clone("mixedratioPeak");
    mixedratioPeak->Divide(uncorrhKKMixed2Dpeak);
    TH1D* mixedratioPeakdeta = mixedratioPeak->ProjectionX("mixedratioPeakdeta", 1, mixedratioPeak->GetYaxis()->GetNbins());
    TH1D* mixedratioPeakdphi = mixedratioPeak->ProjectionY("mixedratioPeakdphi", mixedratioPeak->GetXaxis()->FindBin(-1.2), mixedratioPeak->GetXaxis()->FindBin(1.2));
    mixedratioPeak->GetXaxis()->SetRange(0,0);

    //make ratio plot of just 1 zvtx bin as a check:
    hPhiMixed[6]->GetAxis(2)->SetRangeUser(1.01, 1.03);
    TH2D* mixedratioPeakZ2 = hPhiMixed[6]->Projection(0,1);
    hKKMixed[6]->GetAxis(2)->SetRangeUser(1.01, 1.03);
    TH2D* hist = hKKMixed[6]->Projection(0,1);
    mixedratioPeakZ2->Divide(hist);
    TH1D* mixedratioPeakZ2deta = mixedratioPeakZ2->ProjectionX("mixedratioPeakZ2deta");


    
    TFile* output = new TFile(Form("trig_%i_%i_assoc_%i_%i_mixcorr_hPhi%s.root", (int)trigPTLow, (int)trigPTHigh, (int)assocPTLow, (int)assocPTHigh, mult.c_str()), "RECREATE");
    hPhi2DpeakTotal->Write();
    hKK2DpeakTotal->Write();
    hPhi2DRsideTotal->Write();
    hKK2DRsideTotal->Write();
    hPhi2DLsideTotal->Write();
    hKK2DLsideTotal->Write();
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
/*
    for(int i=0; i< numbinsZvtx; i++){
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
*/
    trigSameUSDist->Write();
    trigSameLSDist->Write();

}
