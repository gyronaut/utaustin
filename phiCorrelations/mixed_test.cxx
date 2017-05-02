void mixed_test(){
    TFile* file = new TFile("~/phiStudies/phiCorrelations_LHC16q_FAST_mixedtest.root");
    if(!file) return;
    file->cd("PhiReconstruction");
    
    THnSparse *fUSDist = (THnSparse*) InvMass->FindObject("fkkUSDist");
    TH1D *phidist = fUSDist->Projection(2);
    THnSparse *fTrigDist1 = (THnSparse*) InvMass->FindObject("fTrigDist");
    fTrigDist1->SetName("fTrigDist1");
    TH1D *trig1 = fTrigDist1->Projection(1);


    if(!fUSDist){
        printf("didn't load histogram!!\n");
        return;
    }

    //file->Close();

    TFile* trigfile = new TFile("~/phiStudies/phiCorrelations_LHC16q_CENT_wSDD_mixedtest.root");
    trigfile->cd("PhiReconstruction");

    THnSparse *fTrigDist2 = (THnSparse*) InvMass->FindObject("fTrigDist");
    fTrigDist2->SetName("fTrigDist2");
    TH1D *trig2 = fTrigDist2->Projection(1);

    //trigfile->Close();

    if(!trig2){
        printf("closing file killed the projection!\n");
        return;
    }
    /*TH1D *trig = new TH1D("trig", "trig", 64, 0.0, 6.28);
    for(int i = 0; i < 64; i++){
        flat->SetBinContent(i+1, 1.0);
    }*/

    TH1D *dphi = new TH1D("dphi", "dphi", 64, -1.57, 4.71);
    TH1D *dphiMix = new TH1D("dphiMix", "dphiMix", 64, -1.57, 4.71);
    TH2D *phiVphi = new TH2D("phiVphi", "#varphi_{h} vs. #varphi_{K^{+}K^{-}}; #varphi_{h}; #varphi_{K^{+}K^{-}}", 64, 0.0, 6.28, 64, 0.0, 6.28);

    for(int j = 0; j < 1000000000; j++){
        Double_t randomTrig = trig1->GetRandom();
        Double_t randomTrigMix = trig2->GetRandom();
        Double_t randomPhi = (phidist->GetRandom() - TMath::Pi());
        if(randomPhi < 0){
            randomPhi += 2.0*TMath::Pi();
        }

        Double_t random = randomTrig - randomPhi;
        if(random < -TMath::Pi()/2.0){
            random += 2.0*TMath::Pi();
        }else if(random > 3.0*TMath::Pi()/2.0){
            random -= 2.0*TMath::Pi();
        }
        
        Double_t randomMix = randomTrigMix - randomPhi;
        if(randomMix < -TMath::Pi()/2.0){
            randomMix += 2.0*TMath::Pi();
        }else if(randomMix > 3.0*TMath::Pi()/2.0){
            randomMix -= 2.0*TMath::Pi();
        }

        dphi->Fill(random);
        dphiMix->Fill(randomMix);
        phiVphi->Fill(randomTrigMix, randomPhi);
    }

    TH1D *ratio = dphi->Clone("ratio");
    ratio->Divide(dphi, dphiMix);

    TFile* out = new TFile("truemix_testoutput.root", "RECREATE");

    dphi->Write();
    dphiMix->Write();
    ratio->Write();
    phiVphi->Write();
    trig1->Write();
    trig2->Write();
    phidist->Write();
}
