void plotMultEff(THnSparse* reco, THnSparse* real, TH1D** recoVar, TH1D** realVar, TH1D** eff, TH1D** ratio, Float_t* mult, Int_t numMultBins, Int_t axis, TCanvas* ceff, TCanvas* cratio){

    TString var = "";
    TString label = "";
    Int_t rebin = 1;

    Int_t colors[] = {kGreen+4,kGreen-2,kSpring+5,kRed+1};

    printf("num mult bins: %i\n", numMultBins);
    switch(axis){
        case 0: var = "PT";
                label = "p_{T}";
                rebin = 5;
                break;
        case 1: var = "Phi";
                label = "#varphi";
                rebin = 4;
                break;
        case 2: var = "Eta";
                label = "#eta";
                rebin = 2;
                break;
        case 3: var = "Y";
                label = "y";
                rebin = 2;
                break;
        case 4: var = "Z";
                label = "Z Vertex (cm);";
                break;
        default: var = "Var";
                 label = "var";
                 break;
    }

    for(int imult = 0; imult < numMultBins; imult++){
        reco->GetAxis(6)->SetRangeUser(mult[imult], mult[imult+1]);
        real->GetAxis(6)->SetRangeUser(mult[imult], mult[imult+1]);
        recoVar[imult] = reco->Projection(axis);
        recoVar[imult]->SetName(Form("reco%s_mult%i%i", var.Data(), int(mult[imult]), int(mult[imult+1])));
        recoVar[imult]->Rebin(rebin);
        realVar[imult] = real->Projection(axis);
        realVar[imult]->SetName(Form("reco%s_mult%i%i", var.Data(), int(mult[imult]), int(mult[imult+1])));
        realVar[imult]->Rebin(rebin);
        eff[imult] = (TH1D*)recoVar[imult]->Clone(Form("eff%s_mult%i%i", var.Data(), int(mult[imult]), int(mult[imult+1])));
        eff[imult]->Divide(realVar[imult]);
        eff[imult]->SetTitle(Form("Efficiency vs. %s for Mult. Bins;%s", label.Data(), label.Data()));
        eff[imult]->SetLineColor(colors[imult]);
        eff[imult]->SetMarkerColor(colors[imult]);
        eff[imult]->SetMarkerSize(1);
        eff[imult]->SetMarkerStyle(22);
    }
 
    reco->GetAxis(6)->SetRangeUser(0.0, 90.0);
    real->GetAxis(6)->SetRangeUser(0.0, 90.0);
    recoVar[numMultBins] = reco->Projection(axis);
    recoVar[numMultBins]->Rebin(rebin);
    recoVar[numMultBins]->SetName(Form("reco%s_mult090", var.Data()));
    realVar[numMultBins] = real->Projection(axis);
    realVar[numMultBins]->SetName(Form("real%s_mult090", var.Data()));
    realVar[numMultBins]->Rebin(rebin);
    eff[numMultBins] = (TH1D*)recoVar[numMultBins]->Clone(Form("eff%s_mult090", var.Data()));
    eff[numMultBins]->Divide(realVar[numMultBins]);
    eff[numMultBins]->SetTitle(Form("Efficiency vs. %s for Mult. Bins", label.Data()));
    eff[numMultBins]->SetLineColor(kRed+1);
    eff[numMultBins]->SetMarkerColor(kRed+1);
    eff[numMultBins]->SetMarkerSize(1);
    eff[numMultBins]->SetMarkerStyle(22);

    cratio->cd();
    for(int i=0; i<numMultBins; i++){
        ratio[i] = (TH1D*)eff[i]->Clone(Form("ratioPT_mult%i", i));
        ratio[i]->Divide(eff[3]);
        ratio[i]->Draw("P SAME");
    } 

    ceff->cd();
    for(int j=0; j<=numMultBins; j++){
        eff[j]->Draw("P SAME");
    }

}; 

void plotEfficiencyBetter(TString input){

    TFile* infile = new TFile(input.Data());
    TList* list = (TList*)infile->Get("phiEff_mult_0_100");

    THnSparseF* recoPhi = (THnSparseF*)list->FindObject("fRecoPhiDist");
    THnSparseF* realPhi = (THnSparseF*)list->FindObject("fRealPhiDist");

    recoPhi->GetAxis(2)->SetRangeUser(-0.8, 0.8);
    realPhi->GetAxis(2)->SetRangeUser(-0.8, 0.8);
    recoPhi->GetAxis(5)->SetRangeUser(1.01, 1.03);
    realPhi->GetAxis(5)->SetRangeUser(1.01, 1.03);
    recoPhi->GetAxis(4)->SetRangeUser(-10.0, 10.0);
    realPhi->GetAxis(4)->SetRangeUser(-10.0, 10.0);

    recoPhi->Sumw2();
    realPhi->Sumw2();

    // pT vs mult
    TH1D* recoPhi_PT_mult[4];
    TH1D* realPhi_PT_mult[4];
    TH1D* effPT_mult[4];
    TH1D* ratioPT_mult[3];
    Float_t mult[4] = {0.0, 20.0, 50.0, 90.0};
    TCanvas* cratioPT = new TCanvas("cratioPT", "cratioPT", 55, 55, 600, 600);
    TCanvas* ceffPT = new TCanvas("ceffPT", "ceffPTmult", 55, 55, 600, 600);

    plotMultEff(recoPhi, realPhi, recoPhi_PT_mult, realPhi_PT_mult, effPT_mult, ratioPT_mult, mult, 3, 0, ceffPT, cratioPT);


    recoPhi->GetAxis(0)->SetRangeUser(1.0, 10.0);
    realPhi->GetAxis(0)->SetRangeUser(1.0, 10.0);
    //phi vs. Mult
    TH1D* recoPhi_Phi_mult[4];
    TH1D* realPhi_Phi_mult[4];
    TH1D* effPhi_mult[4];
    TH1D* ratioPhi_mult[3];
    TCanvas* ceffphi = new TCanvas("ceffphi", "ceffphimult", 55, 55, 600, 600);
    TCanvas* cratioPhi = new TCanvas("cratioPhi", "cratioPhi", 55, 55, 600, 600);

    plotMultEff(recoPhi, realPhi, recoPhi_Phi_mult, realPhi_Phi_mult, effPhi_mult, ratioPhi_mult, mult, 3, 1, ceffphi, cratioPhi);

    
    //eta vs. Mult
    TH1D* recoPhi_Eta_mult[4];
    TH1D* realPhi_Eta_mult[4];
    TH1D* effEta_mult[4];
    TH1D* ratioEta_mult[4];
    TCanvas* ceffeta = new TCanvas("ceffeta", "ceffetamult", 55, 55, 600, 600);
    TCanvas* cratioEta = new TCanvas("cratioEta", "cratioEta", 55, 55, 600, 600);

    plotMultEff(recoPhi, realPhi, recoPhi_Eta_mult, realPhi_Eta_mult, effEta_mult, ratioEta_mult, mult, 3, 2, ceffeta, cratioEta);


    //Z vs. Mult
    TH1D* recoPhi_Z_mult[4];
    TH1D* realPhi_Z_mult[4];
    TH1D* effZ_mult[4];
    TH1D* ratioZ_mult[4];
    TCanvas* ceffZ = new TCanvas("ceffZ", "ceffZmult", 55, 55, 600, 600);
    TCanvas* cratioZ = new TCanvas("cratioZ", "cratioZ", 55, 55, 600, 600);
    
    plotMultEff(recoPhi, realPhi, recoPhi_Z_mult, realPhi_Z_mult, effZ_mult, ratioZ_mult, mult, 3, 4, ceffZ, cratioZ);

    // 0-100% mult with Z vertex dependency 
    realPhi->GetAxis(6)->SetRangeUser(0.0, 100.0);
    recoPhi->GetAxis(6)->SetRangeUser(0.0, 100.0);

    TH1D* recoPhi_PT = recoPhi->Projection(0);
    recoPhi_PT->SetName("recoPT");
    recoPhi_PT->Rebin(5);
    TH1D* realPhi_PT = realPhi->Projection(0);
    realPhi_PT->SetName("realPT");
    realPhi_PT->Rebin(5);

    TH1D* effPT = (TH1D*)recoPhi_PT->Clone("Efficiency vs. p_{T};p_{T}");
    effPT->Divide(realPhi_PT);
    effPT->SetTitle("Efficiency vs. p_{T};p_{T} (GeV/c)");
   
    TCanvas* c1 = new TCanvas("c1", "recon pT", 50, 50, 600, 600);
    c1->cd();
    recoPhi_PT->Draw();

    TCanvas* c2 = new TCanvas("c2", "real pT", 50, 50, 600, 600);
    c2->cd();
    realPhi_PT->Draw();

    TCanvas* c3 = new TCanvas ("c3", "eff PT", 50, 50, 600, 600);
    c3->cd();
    effPT->Draw();

    //recoPhi->GetAxis(0)->SetRangeUser(2.0, 4.0);
    //realPhi->GetAxis(0)->SetRangeUser(2.0, 4.0);

    TH1D* recoPhi_Phi = recoPhi->Projection(1);
    recoPhi_Phi->SetName("recoPhiPhi");

    TH1D* realPhi_Phi = realPhi->Projection(1);
    realPhi_Phi->SetName("realPhiPhi");


    TH1D* effPhi = (TH1D*)recoPhi_Phi->Clone("Efficiency vs. #varphi for 2.0 < p^{#phi}_{T} < 4.0; #phi");
    effPhi->Divide(realPhi_Phi);
    effPhi->SetTitle("Efficiency vs. #varphi for 2.0 < p^{#phi}_{T} < 4.0 GeV/c; #phi");

    TCanvas* c4 = new TCanvas("c4" , "recon phi", 50, 50, 600, 600);
    c4->cd();
    recoPhi_Phi->Draw();

    TCanvas* c5 = new TCanvas("c5", "real phi", 50, 50, 600, 600);
    c5->cd();
    realPhi_Phi->Draw();

    TCanvas* c6 = new TCanvas("c6", "eff phi", 55, 55, 600, 600);
    c6->cd();
    effPhi->Draw();

    TH1D* recoPhi_Eta = recoPhi->Projection(2);
    recoPhi_Eta->SetName("recoPhiEta");
    TH1D* realPhi_Eta = realPhi->Projection(2);
    realPhi_Eta->SetName("realPhiEta");

    TH1D* effEta = (TH1D*)recoPhi_Eta->Clone("Efficiency vs. #eta for 2.0 < p^{#phi}_{T} < 4.0 GeV/c; #eta");
    effEta->Divide(realPhi_Eta);
    effEta->SetTitle("Efficiency vs. #eta for 2.0 < p^{#phi}_{T} < 4.0; #eta");

    TCanvas* c7 = new TCanvas("c7" , "recon eta", 50, 50, 600, 600);
    c7->cd();
    recoPhi_Eta->Draw();

    TCanvas* c8 = new TCanvas("c8", "real eta", 50, 50, 600, 600);
    c8->cd();
    realPhi_Eta->Draw();

    TCanvas* c9 = new TCanvas("c9", "eff eta", 55, 55, 600, 600);
    c9->cd();
    effEta->Draw();

    
    TH1D* recoPhi_Z = recoPhi->Projection(4);
    recoPhi_Z->SetName("recoPhiZ");
    TH1D* realPhi_Z = realPhi->Projection(4);
    realPhi_Z->SetName("realPhiZ");

    TH1D* effZ = (TH1D*)recoPhi_Z->Clone("Efficiency vs. ZVertex for 2.0 < p^{#phi}_{T} < 4.0 GeV/c; Z Vertex (cm)");
    effZ->Divide(realPhi_Z);
    effZ->SetTitle("Efficiency vs. Zvertex for 2.0 < p^{#phi}_{T} < 4.0; Z Vertex (cm)");

    TCanvas* c10 = new TCanvas("c10" , "recon Z", 50, 50, 600, 600);
    c10->cd();
    recoPhi_Z->Draw();

    TCanvas* c11 = new TCanvas("c11", "real Z", 50, 50, 600, 600);
    c11->cd();
    realPhi_Z->Draw();

    TCanvas* c12 = new TCanvas("c12", "eff Z", 55, 55, 600, 600);
    c12->cd();
    effZ->Draw();

    


    recoPhi->GetAxis(0)->SetRange(0,0);
    realPhi->GetAxis(0)->SetRange(0,0);
    TH1D* recoPhi_PT_Z[10];
    TH1D* realPhi_PT_Z[10];
    TH1D* effPT_Z[10];

    //Int_t colors[10] = {kRed+1, kOrange+2, kOrange-2, kSpring+9, kGreen+2, kTeal+2, kCyan+2, kBlue+1, kViolet+7, kViolet-3};
    Int_t colors[10] = {kRed-10, kRed-9, kRed-7, kRed-3, kRed-4, kRed, kRed+1, kRed+2, kRed+3, kRed+4};


    TCanvas* cptz = new TCanvas("cptz", "eff pt z", 60, 60, 600, 600);
    cptz->cd();
    for(int iz = 0; iz < 10; iz++){
        recoPhi->GetAxis(4)->SetRange(iz+1, iz+1);
        realPhi->GetAxis(4)->SetRange(iz+1, iz+1);

        recoPhi_PT_Z[iz] = recoPhi->Projection(0);
        recoPhi_PT_Z[iz]->SetName(Form("reco_PT_Z_%i", iz));
        recoPhi_PT_Z[iz]->Rebin(5);

        realPhi_PT_Z[iz] = realPhi->Projection(0);
        realPhi_PT_Z[iz]->SetName(Form("real_PT_Z_%i", iz));
        realPhi_PT_Z[iz]->Rebin(5);

        effPT_Z[iz] = (TH1D*)recoPhi_PT_Z[iz]->Clone(Form("effPT_Z_%i", iz));
        effPT_Z[iz]->Divide(realPhi_PT_Z[iz]);
        effPT_Z[iz]->SetTitle("Efficiency vs. p_{T} for ZVertex bins:p_{T} (GeV/c)");
        effPT_Z[iz]->SetLineColor(colors[iz]);
        effPT_Z[iz]->SetMarkerColor(colors[iz]);
        effPT_Z[iz]->SetMarkerSize(2);
        effPT_Z[iz]->SetMarkerStyle(22);
        //effPT_Z[iz]->SetLineWidth(2);
        effPT_Z[iz]->Draw("P SAME");
    }

    //recoPhi->GetAxis(0)->SetRangeUser(2.0,4.0);
    //realPhi->GetAxis(0)->SetRange(2.0,4.0);
    TH1D* recoPhi_Phi_Z[10];
    TH1D* realPhi_Phi_Z[10];
    TH1D* effPhi_Z[10];

    TCanvas* cphiz = new TCanvas("cphiz", "eff phi z", 60, 60, 600, 600);
    cphiz->cd();
    for(int iz = 0; iz < 10; iz++){
        recoPhi->GetAxis(4)->SetRange(iz+1, iz+1);
        realPhi->GetAxis(4)->SetRange(iz+1, iz+1);

        recoPhi_Phi_Z[iz] = recoPhi->Projection(1);
        recoPhi_Phi_Z[iz]->SetName(Form("reco_Phi_Z_%i", iz));
        recoPhi_Phi_Z[iz]->Rebin(4);

        realPhi_Phi_Z[iz] = realPhi->Projection(1);
        realPhi_Phi_Z[iz]->SetName(Form("real_Phi_Z_%i", iz));
        realPhi_Phi_Z[iz]->Rebin(4);

        effPhi_Z[iz] = (TH1D*)recoPhi_Phi_Z[iz]->Clone(Form("effPhi_Z_%i", iz));
        effPhi_Z[iz]->Divide(realPhi_Phi_Z[iz]);
        effPhi_Z[iz]->SetTitle("Efficiency vs. #varphi for ZVertex bins;#varphi");
        effPhi_Z[iz]->SetLineColor(colors[iz]);
        effPhi_Z[iz]->SetMarkerColor(colors[iz]);
        effPhi_Z[iz]->SetMarkerStyle(23);
        effPhi_Z[iz]->SetMarkerSize(2);
        //effPhi_Z[iz]->SetLineWidth(2);
        effPhi_Z[iz]->Draw("P SAME");
    }

    recoPhi->GetAxis(2)->SetRange(0,0);
    realPhi->GetAxis(2)->SetRange(0,0);
    TH1D* recoPhi_Eta_Z[10];
    TH1D* realPhi_Eta_Z[10];
    TH1D* effEta_Z[10];

    TCanvas* cetaz = new TCanvas("cetaz", "eff eta z", 60, 60, 600, 600);
    cetaz->cd();
    for(int iz = 0; iz < 10; iz++){
        recoPhi->GetAxis(4)->SetRange(iz+1, iz+1);
        realPhi->GetAxis(4)->SetRange(iz+1, iz+1);

        recoPhi_Eta_Z[iz] = recoPhi->Projection(2);
        recoPhi_Eta_Z[iz]->SetName(Form("reco_Eta_Z_%i", iz));
        recoPhi_Eta_Z[iz]->Rebin(4);

        realPhi_Eta_Z[iz] = realPhi->Projection(2);
        realPhi_Eta_Z[iz]->SetName(Form("real_Eta_Z_%i", iz));
        realPhi_Eta_Z[iz]->Rebin(4);

        effEta_Z[iz] = (TH1D*)recoPhi_Eta_Z[iz]->Clone(Form("effEta_Z_%i", iz));
        effEta_Z[iz]->Divide(realPhi_Eta_Z[iz]);
        effEta_Z[iz]->SetTitle("Efficiency vs. #eta for ZVertex bins;#eta");
        effEta_Z[iz]->SetLineColor(colors[iz]);
        effEta_Z[iz]->SetMarkerColor(colors[iz]);
        effEta_Z[iz]->SetMarkerStyle(23);
        effEta_Z[iz]->SetMarkerSize(2);
        //effEta_Z[iz]->SetLineWidth(2);
        effEta_Z[iz]->Draw("P SAME");
    }


}
