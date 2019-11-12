void checkMassVsDEta(string inputName, int multLow, int multHigh, float trigPTLow, float trigPTHigh, float assocPTLow, float assocPTHigh){

    TFile *histoFile = new TFile(inputName.c_str());
    string mult = "_" + std::to_string(multLow) + "_" + std::to_string(multHigh);
    TList* list = (TList*) histoFile->Get(Form("truePhiCorr_mult%s", mult.c_str()));
    list->SetName("histoList");

    TFile *trueFile = new TFile("trueResults_etacut_02_07.root");
    TList* trueList = (TList*)trueFile->Get(Form("truePhiCorr_mult%s", mult.c_str()));
   
    THnSparseF* trigDist = (THnSparseF*)list->FindObject("fTrigDist");
    TH1D* trigDist1D = (TH1D*)trigDist->Projection(0);
    float totalTrigSameUS = (float)trigDist1D->Integral(trigDist1D->GetXaxis()->FindBin(trigPTLow), trigDist1D->GetXaxis()->FindBin(trigPTHigh));

    trigDist = (THnSparseF*)trueList->FindObject("fTrigDist");
    trigDist1D = (TH1D*)trigDist->Projection(0);
    float totalTrigTrue = (float)trigDist1D->Integral(trigDist1D->GetXaxis()->FindBin(trigPTLow), trigDist1D->GetXaxis()->FindBin(trigPTHigh));

    TH1D* vtxZmixbins = (TH1D*)list->FindObject("fVtxZmixbins");
    Int_t numbinsZvtx = vtxZmixbins->GetXaxis()->GetNbins();


    THnSparseF *dphiHPhi[numbinsZvtx];
    THnSparseF *dphiHKK[numbinsZvtx];
    THnSparseF *dphiHPhiMixed[numbinsZvtx];
    THnSparseF *dphiHKKMixed[numbinsZvtx];

    THnSparseF *dphiHPhiTrue[numbinsZvtx];

    TH3D *hPhi[numbinsZvtx];
    TH3D *hPhiTrue[numbinsZvtx];
    TH3D *hKK[numbinsZvtx];

    TH1D* USmassLowDEta[numbinsZvtx];
    TH1D* USmassHighDEta[numbinsZvtx];
    TH1D* LSmassLowDEta[numbinsZvtx];
    TH1D* LSmassHighDEta[numbinsZvtx];
    TH1D* USmass[numbinsZvtx];
    TH1D* LSmass[numbinsZvtx];
    TH1D* trueMass[numbinsZvtx];

    TH1D* USmassLowDEtaTotal;
    TH1D* USmassHighDEtaTotal;
    TH1D* USmassTotal;
    TH1D* LSmassLowDEtaTotal;
    TH1D* LSmassHighDEtaTotal;
    TH1D* LSmassTotal;
    TH1D* trueMassTotal;

    for(int izvtx = 0; izvtx < numbinsZvtx; izvtx++){

        dphiHPhi[izvtx] = (THnSparseF *)list->FindObject(Form("fDphiHPhiz%i", izvtx));
        dphiHKK[izvtx] = (THnSparseF *)list->FindObject(Form("fDphiHKKz%i", izvtx));
        dphiHPhiTrue[izvtx] = (THnSparseF *)trueList->FindObject(Form("fDphiTrueAcceptanceHPhiz%i", izvtx));

        dphiHPhi[izvtx]->GetAxis(0)->SetRangeUser(trigPTLow, trigPTHigh); 
        dphiHPhi[izvtx]->GetAxis(1)->SetRangeUser(assocPTLow,assocPTHigh); 
        dphiHKK[izvtx]->GetAxis(0)->SetRangeUser(trigPTLow, trigPTHigh);
        dphiHKK[izvtx]->GetAxis(1)->SetRangeUser(assocPTLow,assocPTHigh);
        dphiHPhiTrue[izvtx]->GetAxis(0)->SetRangeUser(trigPTLow, trigPTHigh); 
        dphiHPhiTrue[izvtx]->GetAxis(1)->SetRangeUser(assocPTLow,assocPTHigh);

        dphiHPhi[izvtx]->GetAxis(4)->SetRange(1,dphiHPhi[izvtx]->GetAxis(4)->GetNbins()); 
        dphiHPhiTrue[izvtx]->GetAxis(4)->SetRange(1,dphiHPhi[izvtx]->GetAxis(4)->GetNbins());
        dphiHKK[izvtx]->GetAxis(4)->SetRange(1,dphiHKK[izvtx]->GetAxis(4)->GetNbins());
        
        hPhi[izvtx] = (TH3D*)dphiHPhi[izvtx]->Projection(2,3,4);
        hPhiTrue[izvtx] = (TH3D*)dphiHPhiTrue[izvtx]->Projection(2,3,4);
        hKK[izvtx] = (TH3D*)dphiHKK[izvtx]->Projection(2,3,4);

        hPhi[izvtx]->GetXaxis()->SetRangeUser(-0.499, 0.499);
        USmassLowDEta[izvtx] = (TH1D*)hPhi[izvtx]->Project3D("ze");
        USmassLowDEta[izvtx]->SetName(Form("USmassLowDEta%d", izvtx));
        hPhi[izvtx]->GetXaxis()->SetRangeUser(-2.0, -0.5);
        USmassHighDEta[izvtx] = (TH1D*)hPhi[izvtx]->Project3D("ze");
        USmassHighDEta[izvtx]->SetName(Form("USmassHighDEta%d", izvtx));
        hPhi[izvtx]->GetXaxis()->SetRangeUser(0.5, 2.0);
        TH1D* USbuffer = (TH1D*)hPhi[izvtx]->Project3D("ze");
        USbuffer->SetName("USbuffer");
        USmassHighDEta[izvtx]->Add(USbuffer);
        hPhi[izvtx]->GetXaxis()->SetRangeUser(-1.2001, 1.1999);
        USmass[izvtx] = (TH1D*)hPhi[izvtx]->Project3D("ze");
        USmass[izvtx]->SetName(Form("USmass%d", izvtx));

        hKK[izvtx]->GetXaxis()->SetRangeUser(-0.499, 0.499);
        LSmassLowDEta[izvtx] = (TH1D*)hKK[izvtx]->Project3D("ze");
        LSmassLowDEta[izvtx]->SetName(Form("LSmassLowDEta%d", izvtx));
        hKK[izvtx]->GetXaxis()->SetRangeUser(-2.0, -0.5);
        LSmassHighDEta[izvtx] = (TH1D*)hKK[izvtx]->Project3D("ze");
        LSmassHighDEta[izvtx]->SetName(Form("LSmassHighDEta%d", izvtx));
        hKK[izvtx]->GetXaxis()->SetRangeUser(0.5, 2.0);
        TH1D* LSbuffer = (TH1D*)hKK[izvtx]->Project3D("ze");
        LSbuffer->SetName("LSbuffer");
        LSmassHighDEta[izvtx]->Add(LSbuffer);
        hKK[izvtx]->GetXaxis()->SetRangeUser(-1.2001, 1.1999);
        LSmass[izvtx] = (TH1D*)hKK[izvtx]->Project3D("ze");
        LSmass[izvtx]->SetName(Form("LSmass%d", izvtx));

        hPhiTrue[izvtx]->GetXaxis()->SetRangeUser(-1.2001, 1.1999);
        trueMass[izvtx] = (TH1D*)hPhiTrue[izvtx]->Project3D("ze");
        trueMass[izvtx]->SetName(Form("trueMass%d", izvtx));

        if(izvtx == 0){
            USmassLowDEtaTotal = (TH1D*)USmassLowDEta[izvtx]->Clone("USmassLowDEta");
            USmassHighDEtaTotal = (TH1D*)USmassHighDEta[izvtx]->Clone("USmassHighDEta");
            USmassTotal = (TH1D*)USmass[izvtx]->Clone("USmassTotal");
            LSmassLowDEtaTotal = (TH1D*)LSmassLowDEta[izvtx]->Clone("LSmassLowDEta");
            LSmassHighDEtaTotal = (TH1D*)LSmassHighDEta[izvtx]->Clone("LSmassHighDEta");
            LSmassTotal = (TH1D*)LSmass[izvtx]->Clone("LSmassTotal");
            trueMassTotal = (TH1D*)trueMass[izvtx]->Clone("trueMassTotal");
        }else{
            USmassLowDEtaTotal->Add(USmassLowDEta[izvtx]);
            USmassHighDEtaTotal->Add(USmassHighDEta[izvtx]);
            USmassTotal->Add(USmass[izvtx]);
            LSmassLowDEtaTotal->Add(LSmassLowDEta[izvtx]);
            LSmassHighDEtaTotal->Add(LSmassHighDEta[izvtx]);
            LSmassTotal->Add(LSmass[izvtx]);
            trueMassTotal->Add(trueMass[izvtx]);
        }
    }
    
    USmassLowDEtaTotal->Scale(1.0/USmassLowDEtaTotal->Integral(1, USmassLowDEtaTotal->GetXaxis()->GetNbins()));
    USmassHighDEtaTotal->Scale(1.0/USmassHighDEtaTotal->Integral(1, USmassHighDEtaTotal->GetXaxis()->GetNbins()));
    //USmassTotal->Scale(1.0/USmassTotal->Integral(1, USmassTotal->GetXaxis()->GetNbins()));
    USmassTotal->Scale(1.0/totalTrigSameUS);

    LSmassLowDEtaTotal->Scale(1.0/LSmassLowDEtaTotal->Integral(1, LSmassLowDEtaTotal->GetXaxis()->GetNbins()));
    LSmassHighDEtaTotal->Scale(1.0/LSmassHighDEtaTotal->Integral(1, LSmassHighDEtaTotal->GetXaxis()->GetNbins()));
    //LSmassTotal->Scale(1.0/LSmassTotal->Integral(1, LSmassTotal->GetXaxis()->GetNbins()));
    
    trueMassTotal->Scale(1.0/totalTrigTrue);

    printf("kaon trig: %f\ntrue trig: %f\n", totalTrigSameUS, totalTrigTrue);

    USmassLowDEtaTotal->SetLineColor(kBlue);
    USmassHighDEtaTotal->SetLineColor(kBlue+2);
    USmassTotal->SetLineColor(kCyan+1);
    LSmassLowDEtaTotal->SetLineColor(kRed);
    LSmassHighDEtaTotal->SetLineColor(kRed+2);
    LSmassTotal->SetLineColor(kMagenta+1);

    TCanvas *cUS = new TCanvas("cUS", "cUS", 50, 50, 600, 600);
    cUS->cd();
    USmassLowDEtaTotal->Draw();
    USmassHighDEtaTotal->Draw("SAME");
    USmassTotal->Draw("SAME");

    TCanvas *cLS = new TCanvas("cLS", "cLS", 50, 50, 600, 600);
    cLS->cd();
    LSmassLowDEtaTotal->Draw();
    LSmassHighDEtaTotal->Draw("SAME");
    LSmassTotal->Draw("SAME");

    TH1D* scaledLSmassLowDEta = (TH1D*)LSmassLowDEtaTotal->Clone("scaledLSmassLowDEta");
    scaledLSmassLowDEta->Scale(USmassLowDEtaTotal->Integral(USmassLowDEtaTotal->GetXaxis()->FindBin(1.04001), USmassLowDEtaTotal->GetXaxis()->FindBin(1.05999))/LSmassLowDEtaTotal->Integral(LSmassLowDEtaTotal->GetXaxis()->FindBin(1.04001), LSmassLowDEtaTotal->GetXaxis()->FindBin(1.05999)));

    TH1D* scaledLSmassHighDEta = (TH1D*)LSmassHighDEtaTotal->Clone("scaledLSmassHighDEta");
    scaledLSmassHighDEta->Scale(USmassHighDEtaTotal->Integral(USmassHighDEtaTotal->GetXaxis()->FindBin(1.04001), USmassHighDEtaTotal->GetXaxis()->FindBin(1.05999))/LSmassHighDEtaTotal->Integral(LSmassHighDEtaTotal->GetXaxis()->FindBin(1.04001), LSmassHighDEtaTotal->GetXaxis()->FindBin(1.05999)));


    TH1D* scaledLSmassTotal = (TH1D*)LSmassTotal->Clone("scaledLSmassTotal");
    scaledLSmassTotal->Scale(USmassTotal->Integral(USmassTotal->GetXaxis()->FindBin(1.04001), USmassTotal->GetXaxis()->FindBin(1.05999))/LSmassTotal->Integral(LSmassTotal->GetXaxis()->FindBin(1.04001), LSmassTotal->GetXaxis()->FindBin(1.05999)));

    TCanvas *cLow = new TCanvas("cLow", "cLow", 60, 60, 600, 600);
    cLow->cd();
    USmassLowDEtaTotal->Draw();
    scaledLSmassLowDEta->Draw("SAME");

    TCanvas *cHigh = new TCanvas("cHigh", "cHigh", 60, 60, 600, 600);
    cHigh->cd();
    USmassHighDEtaTotal->Draw();
    scaledLSmassHighDEta->Draw("SAME");

    float highYield = USmassHighDEtaTotal->Integral(USmassHighDEtaTotal->GetXaxis()->FindBin(1.014001), USmassHighDEtaTotal->GetXaxis()->FindBin(1.025999)) - scaledLSmassHighDEta->Integral(scaledLSmassHighDEta->GetXaxis()->FindBin(1.014001), scaledLSmassHighDEta->GetXaxis()->FindBin(1.025999));
    float lowYield = USmassLowDEtaTotal->Integral(USmassLowDEtaTotal->GetXaxis()->FindBin(1.014001), USmassLowDEtaTotal->GetXaxis()->FindBin(1.025999)) - scaledLSmassLowDEta->Integral(scaledLSmassLowDEta->GetXaxis()->FindBin(1.014001), scaledLSmassLowDEta->GetXaxis()->FindBin(1.025999));
    float ratio = highYield/lowYield;

    printf("High Yield: %f\nLow  Yield: %f\nRatio     : %f\n", highYield, lowYield, ratio);

    TH1D* otherscaledLSmassLowDEta = (TH1D*)LSmassLowDEtaTotal->Clone("otherscaledLSmassLowDEta");
    otherscaledLSmassLowDEta->Scale(USmassLowDEtaTotal->Integral(USmassLowDEtaTotal->GetXaxis()->FindBin(0.99501), USmassLowDEtaTotal->GetXaxis()->FindBin(1.00499))/LSmassLowDEtaTotal->Integral(LSmassLowDEtaTotal->GetXaxis()->FindBin(0.99501), LSmassLowDEtaTotal->GetXaxis()->FindBin(1.00499)));

    TH1D* otherscaledLSmassHighDEta = (TH1D*)LSmassHighDEtaTotal->Clone("otherscaledLSmassHighDEta");
    otherscaledLSmassHighDEta->Scale(USmassHighDEtaTotal->Integral(USmassHighDEtaTotal->GetXaxis()->FindBin(0.99501), USmassHighDEtaTotal->GetXaxis()->FindBin(1.00499))/LSmassHighDEtaTotal->Integral(LSmassHighDEtaTotal->GetXaxis()->FindBin(0.99501), LSmassHighDEtaTotal->GetXaxis()->FindBin(1.00499)));


    TCanvas *cotherLow = new TCanvas("cotherLow", "cotherLow", 60, 60, 600, 600);
    cotherLow->cd();
    USmassLowDEtaTotal->Draw();
    otherscaledLSmassLowDEta->Draw("SAME");

    TCanvas *cotherHigh = new TCanvas("cotherHigh", "cotherHigh", 60, 60, 600, 600);
    cotherHigh->cd();
    USmassHighDEtaTotal->Draw();
    otherscaledLSmassHighDEta->Draw("SAME");

    float otherhighYield = USmassHighDEtaTotal->Integral(USmassHighDEtaTotal->GetXaxis()->FindBin(1.014001), USmassHighDEtaTotal->GetXaxis()->FindBin(1.025999)) - otherscaledLSmassHighDEta->Integral(scaledLSmassHighDEta->GetXaxis()->FindBin(1.014001), scaledLSmassHighDEta->GetXaxis()->FindBin(1.025999));
    float otherlowYield = USmassLowDEtaTotal->Integral(USmassLowDEtaTotal->GetXaxis()->FindBin(1.014001), USmassLowDEtaTotal->GetXaxis()->FindBin(1.025999)) - otherscaledLSmassLowDEta->Integral(scaledLSmassLowDEta->GetXaxis()->FindBin(1.014001), scaledLSmassLowDEta->GetXaxis()->FindBin(1.025999));
    float otherratio = otherhighYield/otherlowYield;

    printf("other High Yield: %f\nother Low  Yield: %f\nother Ratio     : %f\n", otherhighYield, otherlowYield, otherratio);


    TH1D* corrMass = (TH1D*)USmassTotal->Clone("corrMass");
    corrMass->Add(scaledLSmassTotal, -1.0);
    trueMassTotal->SetLineColor(kGreen+1);

    TCanvas *ctrue = new TCanvas("ctrue", "ctrue", 60, 60, 600, 600);
    ctrue->cd();
    //corrMass->Draw();
    trueMassTotal->SetTitle("Invariant mass of #phi_{MC}");
    trueMassTotal->GetXaxis()->SetTitle("M_{#phi} (GeV/c^{2})");
    trueMassTotal->GetYaxis()->SetTitle("Counts");
    trueMassTotal->Draw("HIST E SAME");

    float trueyield = trueMassTotal->Integral(trueMassTotal->GetXaxis()->FindBin(1.014001), trueMassTotal->GetXaxis()->FindBin(1.025999));
    float trueyieldTotal = trueMassTotal->Integral(trueMassTotal->GetXaxis()->FindBin(1.00001), trueMassTotal->GetXaxis()->FindBin(1.039999));
    float corryield = corrMass->Integral(corrMass->GetXaxis()->FindBin(1.014001), corrMass->GetXaxis()->FindBin(1.025999));
    float corryieldTotal = corrMass->Integral(corrMass->GetXaxis()->FindBin(1.0001), corrMass->GetXaxis()->FindBin(1.039999));

    float trueratio = corryield/trueyield;
    float recoratio = corryield/corryieldTotal;
    float truerecoratio = trueyield/trueyieldTotal;

    printf("true yield: %f\ncorr yield: %f\nratio    : %f\n", trueyield, corryield, trueratio);

    printf("recoratio: %f\ntrueratio: %f\n", recoratio, truerecoratio);
    
    
}
