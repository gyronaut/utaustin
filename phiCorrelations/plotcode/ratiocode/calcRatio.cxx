void calcRatio(TString filename){
    TFile* input = new TFile(filename.Data());
    TList* list = (TList*) input->Get("phiCorr_mult_50_80_");

    THnSparseF* hadronDist = (THnSparseF*)list->FindObject("fTrigDist");
    THnSparseF* KKUSDist = (THnSparseF*)list->FindObject("fkkUSTrigDist")->Clone("KKUSDist");
    THnSparseF* KKLSDist = (THnSparseF*)list->FindObject("fkkLSTrigDist")->Clone("KKLSDist");

    Float_t eps = 0.000001;
    hadronDist->GetAxis(2)->SetRangeUser(-0.8+eps, 0.8-eps);
    TH1D* hadronPT = (TH1D*)hadronDist->Projection(0); 
    //TH1D* hadronPT = (TH1D*)list->FindObject("fHadronTrigPT");

    TFile* efffile = new TFile("~/alidock/alirepos/utaustin/efficiency/fits_17f2b_secondarytest.root");
    TF1* phiEff = (TF1*)efffile->Get("phiFit")->Clone("phiEff");

    KKUSDist->GetAxis(3)->SetRangeUser(-0.8+eps, 0.8-eps);
    TH2D* KKUS = (TH2D*)KKUSDist->Projection(1,0);
    KKUS->RebinX();

    KKLSDist->GetAxis(3)->SetRangeUser(-0.8+eps, 0.8-eps);
    TH2D* KKLS = (TH2D*)KKLSDist->Projection(1,0);
    KKLS->RebinX();

    Int_t first_bin = KKUS->GetXaxis()->FindBin(2.0+eps);
    Int_t last_bin = KKUS->GetXaxis()->FindBin(4.0-eps);
    
    printf("PhiEff test: %f\n", phiEff->Eval(2.25));

    printf("first bin: %d, last bin: %d\n", first_bin, last_bin);

    TCanvas* c2d = new TCanvas("c2d", "c2d", 50, 50, 600, 600);
    c2d->cd();
    KKUS->Draw("colz");


    TH1D* USmass[4];
    TH1D* LSmass[4];
    TH1D* mass[4];
    Float_t phis = 0.0;
    Float_t bgsubphis = 0.0;
    for(int istep = 0; istep < (last_bin - first_bin)+1; istep++){

        KKUS->GetXaxis()->SetRange(istep+first_bin, istep+first_bin);
        USmass[istep] = (TH1D*)KKUS->ProjectionY(Form("USmass_%d", istep));
        printf("test2: %f\n", USmass[istep]->GetXaxis()->GetBinCenter(1));

        //USmass[istep]->Scale(1.0/phiEff->Eval(KKUS->GetXaxis()->GetBinCenter(istep+first_bin)));
        USmass[istep]->Scale(1.0/0.49);
        KKLS->GetXaxis()->SetRange(istep+first_bin, istep+first_bin);
        LSmass[istep] = (TH1D*)KKLS->ProjectionY(Form("LSmass_%d", istep));

        Float_t scalefactor = USmass[istep]->Integral(USmass[istep]->GetXaxis()->FindBin(1.04 + eps), USmass[istep]->GetXaxis()->FindBin(1.06 - eps));
        scalefactor = scalefactor/(LSmass[istep]->Integral(LSmass[istep]->GetXaxis()->FindBin(1.04 + eps), LSmass[istep]->GetXaxis()->FindBin(1.06 - eps)));

        LSmass[istep]->Scale(scalefactor);
        
        mass[istep] = (TH1D*)USmass[istep]->Clone(Form("mass_%d", istep));
        mass[istep]->Add(LSmass[istep], -1.0);

        TF1* fit = new TF1("fit",  "[0]*TMath::Voigt(x - [1], [2], [3], 4) + pol2(4)",0.99, 1.07);
        fit->SetParameter(1, 1.020);
        fit->SetParameter(2, 0.0002);
        fit->SetParameter(0, 600);
        fit->FixParameter(3, 0.00426);
        fit->SetParLimits(1, 1.010, 1.030);

        mass[istep]->Fit(fit, "R");
        TF1* bgfit = new TF1("bgfit", "pol2(0)", 0.99, 1.07);
        bgfit->SetParameters(fit->GetParameter(4), fit->GetParameter(5), fit->GetParameter(6));
        
        Float_t histoint =  mass[istep]->Integral(mass[istep]->GetXaxis()->FindBin(1.005 + eps), mass[istep]->GetXaxis()->FindBin(1.035 - eps));
        Float_t bgint = bgfit->Integral(1.005, 1.035);
        


        phis += histoint;
        bgsubphis += histoint - bgint;
        printf("phis %d: %f, bgsubphis: %f, diff: %f\n", istep, phis, bgsubphis, (phis-bgsubphis)/phis);
    }
    
    TCanvas* c1 = new TCanvas("c1", "c1", 50, 50, 600, 600);
    c1->cd();
    USmass[1]->Draw();
    LSmass[1]->Draw();

    Float_t hadrons = hadronPT->Integral(hadronPT->GetXaxis()->FindBin(2.0+eps), hadronPT->GetXaxis()->FindBin(4.0-eps));

    printf("phis: %f, hadrons: %f, ratio: %f, bgsubratio: %f\n", phis, hadrons, phis/hadrons, bgsubphis/hadrons);
}

