void plotEfficiency(TString input){

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

    recoPhi->GetAxis(6)->SetRangeUser(0.0, 20.0);
    realPhi->GetAxis(6)->SetRangeUser(0.0, 20.0);
    recoPhi_PT_mult[0] = recoPhi->Projection(0);
    recoPhi_PT_mult[0]->SetName("recoPT_mult020");
    recoPhi_PT_mult[0]->Rebin(5);
    realPhi_PT_mult[0] = realPhi->Projection(0);
    realPhi_PT_mult[0]->SetName("realPT_mult020");
    realPhi_PT_mult[0]->Rebin(5);
    effPT_mult[0] = (TH1D*)recoPhi_PT_mult[0]->Clone("effPT_mult020");
    effPT_mult[0]->Divide(realPhi_PT_mult[0]);
    effPT_mult[0]->SetTitle("Efficiency vs. p_{T} for Mult. Bins;p_{T} (GeV/c)");
    effPT_mult[0]->SetLineColor(kGreen+4);
    effPT_mult[0]->SetMarkerColor(kGreen+4);
    effPT_mult[0]->SetMarkerSize(1);
    effPT_mult[0]->SetMarkerStyle(22);

    recoPhi->GetAxis(6)->SetRangeUser(20.0, 50.0);
    realPhi->GetAxis(6)->SetRangeUser(20.0, 50.0);
    recoPhi_PT_mult[1] = recoPhi->Projection(0);
    recoPhi_PT_mult[1]->SetName("recoPT_mult2050");
    recoPhi_PT_mult[1]->Rebin(5);
    realPhi_PT_mult[1] = realPhi->Projection(0);
    realPhi_PT_mult[1]->SetName("realPT_mult2050");
    realPhi_PT_mult[1]->Rebin(5);
    effPT_mult[1] = (TH1D*)recoPhi_PT_mult[1]->Clone("effPT_mult2050");
    effPT_mult[1]->Divide(realPhi_PT_mult[1]);
    effPT_mult[1]->SetTitle("Efficiency vs. p_{T} for Mult. Bins");
    effPT_mult[1]->SetLineColor(kGreen-2);
    effPT_mult[1]->SetMarkerColor(kGreen-2);
    effPT_mult[1]->SetMarkerSize(1);
    effPT_mult[1]->SetMarkerStyle(22);

    recoPhi->GetAxis(6)->SetRangeUser(50.0, 90.0);
    realPhi->GetAxis(6)->SetRangeUser(50.0, 90.0);
    recoPhi_PT_mult[2] = recoPhi->Projection(0);
    recoPhi_PT_mult[2]->Rebin(5);
    recoPhi_PT_mult[2]->SetName("recoPT_mult5090");
    realPhi_PT_mult[2] = realPhi->Projection(0);
    realPhi_PT_mult[2]->SetName("realPT_mult5090");
    realPhi_PT_mult[2]->Rebin(5);
    effPT_mult[2] = (TH1D*)recoPhi_PT_mult[2]->Clone("effPT_mult5090");
    effPT_mult[2]->Divide(realPhi_PT_mult[2]);
    effPT_mult[2]->SetTitle("Efficiency vs. p_{T} for Mult. Bins");
    effPT_mult[2]->SetLineColor(kSpring+5);
    effPT_mult[2]->SetMarkerColor(kSpring+5);
    effPT_mult[2]->SetMarkerSize(1);
    effPT_mult[2]->SetMarkerStyle(22);

    recoPhi->GetAxis(6)->SetRangeUser(0.0, 90.0);
    realPhi->GetAxis(6)->SetRangeUser(0.0, 90.0);
    recoPhi_PT_mult[3] = recoPhi->Projection(0);
    recoPhi_PT_mult[3]->Rebin(5);
    recoPhi_PT_mult[3]->SetName("recoPT_mult090");
    realPhi_PT_mult[3] = realPhi->Projection(0);
    realPhi_PT_mult[3]->SetName("realPT_mult090");
    realPhi_PT_mult[3]->Rebin(5);
    effPT_mult[3] = (TH1D*)recoPhi_PT_mult[3]->Clone("effPT_mult090");
    effPT_mult[3]->Divide(realPhi_PT_mult[3]);
    effPT_mult[3]->SetTitle("Efficiency vs. p_{T} for Mult. Bins");
    effPT_mult[3]->SetLineColor(kRed+1);
    effPT_mult[3]->SetMarkerColor(kRed+1);
    effPT_mult[3]->SetMarkerSize(1);
    effPT_mult[3]->SetMarkerStyle(22);

    for(int i=0; i<3; i++){
        ratioPT_mult[i] = (TH1D*)effPT_mult[i]->Clone(Form("ratioPT_mult%i", i));
        ratioPT_mult[i]->Divide(effPT_mult[3]);
    }
    TCanvas* cratioPT = new TCanvas("cratioPT", "cratioPT", 55, 55, 600, 600);
    cratioPT->cd();
    ratioPT_mult[0]->Draw("P");
    ratioPT_mult[1]->Draw("P SAME");
    ratioPT_mult[2]->Draw("P SAME");

    TCanvas* ceffPT = new TCanvas("ceffPT", "ceffPTmult", 55, 55, 600, 600);
    ceffPT->cd();
    effPT_mult[0]->Draw("P");
    effPT_mult[1]->Draw("P SAME");
    effPT_mult[2]->Draw("P SAME");
    effPT_mult[3]->Draw("P SAME");

    //phi vs. Mult
    recoPhi->GetAxis(0)->SetRangeUser(1.0, 10.0);
    realPhi->GetAxis(0)->SetRangeUser(1.0, 10.0);


    TH1D* recoPhi_Phi_mult[4];
    TH1D* realPhi_Phi_mult[4];
    TH1D* effPhi_mult[4];
    TH1D* ratioPhi_mult[3];

    recoPhi->GetAxis(6)->SetRangeUser(0.0, 20.0);
    realPhi->GetAxis(6)->SetRangeUser(0.0, 20.0);
    recoPhi_Phi_mult[0] = recoPhi->Projection(1);
    recoPhi_Phi_mult[0]->SetName("recoPhi_mult020");
    recoPhi_Phi_mult[0]->Rebin(4);
    realPhi_Phi_mult[0] = realPhi->Projection(1);
    realPhi_Phi_mult[0]->SetName("realPhi_mult020");
    realPhi_Phi_mult[0]->Rebin(4);
    effPhi_mult[0] = (TH1D*)recoPhi_Phi_mult[0]->Clone("effPhi_mult020");
    effPhi_mult[0]->Divide(realPhi_Phi_mult[0]);
    effPhi_mult[0]->SetTitle("Efficiency vs. #varphi for Mult. Bins;#varphi");
    effPhi_mult[0]->SetLineColor(kGreen+4);
    effPhi_mult[0]->SetMarkerColor(kGreen+4);
    effPhi_mult[0]->SetMarkerSize(1);
    effPhi_mult[0]->SetMarkerStyle(22);

    recoPhi->GetAxis(6)->SetRangeUser(20.0, 50.0);
    realPhi->GetAxis(6)->SetRangeUser(20.0, 50.0);
    recoPhi_Phi_mult[1] = recoPhi->Projection(1);
    recoPhi_Phi_mult[1]->SetName("recoPhi_mult2050");
    recoPhi_Phi_mult[1]->Rebin(4);
    realPhi_Phi_mult[1] = realPhi->Projection(1);
    realPhi_Phi_mult[1]->SetName("realPhi_mult2050");
    realPhi_Phi_mult[1]->Rebin(4);
    effPhi_mult[1] = (TH1D*)recoPhi_Phi_mult[1]->Clone("effPhi_mult2050");
    effPhi_mult[1]->Divide(realPhi_Phi_mult[1]);
    effPhi_mult[1]->SetTitle("Efficiency vs. #varphi for Mult. Bins");
    effPhi_mult[1]->SetLineColor(kGreen-2);
    effPhi_mult[1]->SetMarkerColor(kGreen-2);
    effPhi_mult[1]->SetMarkerSize(1);
    effPhi_mult[1]->SetMarkerStyle(22);

    recoPhi->GetAxis(6)->SetRangeUser(50.0, 90.0);
    realPhi->GetAxis(6)->SetRangeUser(50.0, 90.0);
    recoPhi_Phi_mult[2] = recoPhi->Projection(1);
    recoPhi_Phi_mult[2]->Rebin(4);
    recoPhi_Phi_mult[2]->SetName("recoPhi_mult5090");
    realPhi_Phi_mult[2] = realPhi->Projection(1);
    realPhi_Phi_mult[2]->SetName("realPhi_mult5090");
    realPhi_Phi_mult[2]->Rebin(4);
    effPhi_mult[2] = (TH1D*)recoPhi_Phi_mult[2]->Clone("effPhi_mult5090");
    effPhi_mult[2]->Divide(realPhi_Phi_mult[2]);
    effPhi_mult[2]->SetTitle("Efficiency vs. #varphi for Mult. Bins");
    effPhi_mult[2]->SetLineColor(kSpring+5);
    effPhi_mult[2]->SetMarkerColor(kSpring+5);
    effPhi_mult[2]->SetMarkerSize(1);
    effPhi_mult[2]->SetMarkerStyle(22);

    recoPhi->GetAxis(6)->SetRangeUser(0.0, 90.0);
    realPhi->GetAxis(6)->SetRangeUser(0.0, 90.0);
    recoPhi_Phi_mult[3] = recoPhi->Projection(1);
    recoPhi_Phi_mult[3]->Rebin(4);
    recoPhi_Phi_mult[3]->SetName("recoPhi_mult090");
    realPhi_Phi_mult[3] = realPhi->Projection(1);
    realPhi_Phi_mult[3]->SetName("realPhi_mult090");
    realPhi_Phi_mult[3]->Rebin(4);
    effPhi_mult[3] = (TH1D*)recoPhi_Phi_mult[3]->Clone("effPhi_mult090");
    effPhi_mult[3]->Divide(realPhi_Phi_mult[3]);
    effPhi_mult[3]->SetTitle("Efficiency vs. #varphi for Mult. Bins");
    effPhi_mult[3]->SetLineColor(kRed+1);
    effPhi_mult[3]->SetMarkerColor(kRed+1);
    effPhi_mult[3]->SetMarkerSize(1);
    effPhi_mult[3]->SetMarkerStyle(22);


    TCanvas* ceffphi = new TCanvas("ceffphi", "ceffphimult", 55, 55, 600, 600);
    ceffphi->cd();
    effPhi_mult[0]->Draw("P");
    effPhi_mult[1]->Draw("P SAME");
    effPhi_mult[2]->Draw("P SAME");
    effPhi_mult[2]->Draw("P SAME");

    for(int i=0; i<3; i++){
        ratioPhi_mult[i] = (TH1D*)effPhi_mult[i]->Clone(Form("ratioPhi_mult%i", i));
        ratioPhi_mult[i]->Divide(effPhi_mult[3]);
    }
    TCanvas* cratioPhi = new TCanvas("cratioPhi", "cratioPhi", 55, 55, 600, 600);
    cratioPhi->cd();
    ratioPhi_mult[0]->Draw("P");
    ratioPhi_mult[1]->Draw("P SAME");
    ratioPhi_mult[2]->Draw("P SAME");

    
    //eta vs. Mult
    TH1D* recoPhi_Eta_mult[4];
    TH1D* realPhi_Eta_mult[4];
    TH1D* effEta_mult[4];
    TH1D* ratioEta_mult[4];

    recoPhi->GetAxis(6)->SetRangeUser(0.0, 20.0);
    realPhi->GetAxis(6)->SetRangeUser(0.0, 20.0);
    recoPhi_Eta_mult[0] = recoPhi->Projection(2);
    recoPhi_Eta_mult[0]->SetName("recoEta_mult020");
    recoPhi_Eta_mult[0]->Rebin(2);
    realPhi_Eta_mult[0] = realPhi->Projection(2);
    realPhi_Eta_mult[0]->SetName("realEta_mult020");
    realPhi_Eta_mult[0]->Rebin(2);
    effEta_mult[0] = (TH1D*)recoPhi_Eta_mult[0]->Clone("effEta_mult020");
    effEta_mult[0]->Divide(realPhi_Eta_mult[0]);
    effEta_mult[0]->SetTitle("Efficiency vs. #eta for Mult. Bins;#eta");
    effEta_mult[0]->SetLineColor(kGreen+4);
    effEta_mult[0]->SetMarkerColor(kGreen+4);
    effEta_mult[0]->SetMarkerSize(1);
    effEta_mult[0]->SetMarkerStyle(22);

    recoPhi->GetAxis(6)->SetRangeUser(20.0, 50.0);
    realPhi->GetAxis(6)->SetRangeUser(20.0, 50.0);
    recoPhi_Eta_mult[1] = recoPhi->Projection(2);
    recoPhi_Eta_mult[1]->SetName("recoEta_mult2050");
    recoPhi_Eta_mult[1]->Rebin(2);
    realPhi_Eta_mult[1] = realPhi->Projection(2);
    realPhi_Eta_mult[1]->SetName("realEta_mult2050");
    realPhi_Eta_mult[1]->Rebin(2);
    effEta_mult[1] = (TH1D*)recoPhi_Eta_mult[1]->Clone("effEta_mult2050");
    effEta_mult[1]->Divide(realPhi_Eta_mult[1]);
    effEta_mult[1]->SetTitle("Efficiency vs. #eta for Mult. Bins");
    effEta_mult[1]->SetLineColor(kGreen-2);
    effEta_mult[1]->SetMarkerColor(kGreen-2);
    effEta_mult[1]->SetMarkerSize(1);
    effEta_mult[1]->SetMarkerStyle(22);

    recoPhi->GetAxis(6)->SetRangeUser(50.0, 90.0);
    realPhi->GetAxis(6)->SetRangeUser(50.0, 90.0);
    recoPhi_Eta_mult[2] = recoPhi->Projection(2);
    recoPhi_Eta_mult[2]->Rebin(2);
    recoPhi_Eta_mult[2]->SetName("recoEta_mult5090");
    realPhi_Eta_mult[2] = realPhi->Projection(2);
    realPhi_Eta_mult[2]->SetName("realEta_mult5090");
    realPhi_Eta_mult[2]->Rebin(2);
    effEta_mult[2] = (TH1D*)recoPhi_Eta_mult[2]->Clone("effEta_mult5090");
    effEta_mult[2]->Divide(realPhi_Eta_mult[2]);
    effEta_mult[2]->SetTitle("Efficiency vs. #eta for Mult. Bins");
    effEta_mult[2]->SetLineColor(kSpring+5);
    effEta_mult[2]->SetMarkerColor(kSpring+5);
    effEta_mult[2]->SetMarkerSize(1);
    effEta_mult[2]->SetMarkerStyle(22);

    recoPhi->GetAxis(6)->SetRangeUser(0.0, 90.0);
    realPhi->GetAxis(6)->SetRangeUser(0.0, 90.0);
    recoPhi_Eta_mult[3] = recoPhi->Projection(2);
    recoPhi_Eta_mult[3]->Rebin(2);
    recoPhi_Eta_mult[3]->SetName("recoEta_mult090");
    realPhi_Eta_mult[3] = realPhi->Projection(2);
    realPhi_Eta_mult[3]->SetName("realEta_mult090");
    realPhi_Eta_mult[3]->Rebin(2);
    effEta_mult[3] = (TH1D*)recoPhi_Eta_mult[3]->Clone("effEta_mult090");
    effEta_mult[3]->Divide(realPhi_Eta_mult[3]);
    effEta_mult[3]->SetTitle("Efficiency vs. #eta for Mult. Bins");
    effEta_mult[3]->SetLineColor(kRed+1);
    effEta_mult[3]->SetMarkerColor(kRed+1);
    effEta_mult[3]->SetMarkerSize(1);
    effEta_mult[3]->SetMarkerStyle(22);


    TCanvas* ceffeta = new TCanvas("ceffeta", "ceffetamult", 55, 55, 600, 600);
    ceffeta->cd();
    effEta_mult[0]->Draw("P");
    effEta_mult[1]->Draw("P SAME");
    effEta_mult[2]->Draw("P SAME");
    effEta_mult[3]->Draw("P SAME");
      
    for(int i=0; i<3; i++){
        ratioEta_mult[i] = (TH1D*)effEta_mult[i]->Clone(Form("ratioEta_mult%i", i));
        ratioEta_mult[i]->Divide(effEta_mult[3]);
    }
    TCanvas* cratioEta = new TCanvas("cratioEta", "cratioEta", 55, 55, 600, 600);
    cratioEta->cd();
    ratioEta_mult[0]->Draw("P");
    ratioEta_mult[1]->Draw("P SAME");
    ratioEta_mult[2]->Draw("P SAME");

    //Z vs. Mult
    TH1D* recoPhi_Z_mult[4];
    TH1D* realPhi_Z_mult[4];
    TH1D* effZ_mult[4];
    TH1D* ratioZ_mult[4];

    recoPhi->GetAxis(6)->SetRangeUser(0.0, 20.0);
    realPhi->GetAxis(6)->SetRangeUser(0.0, 20.0);
    recoPhi_Z_mult[0] = recoPhi->Projection(4);
    recoPhi_Z_mult[0]->SetName("recoZ_mult020");
    realPhi_Z_mult[0] = realPhi->Projection(4);
    realPhi_Z_mult[0]->SetName("realZ_mult020");
    effZ_mult[0] = (TH1D*)recoPhi_Z_mult[0]->Clone("effZ_mult020");
    effZ_mult[0]->Divide(realPhi_Z_mult[0]);
    effZ_mult[0]->SetTitle("Efficiency vs. Z Vertex for Mult. Bins;Z Vertex (cm)");
    effZ_mult[0]->SetLineColor(kGreen+4);
    effZ_mult[0]->SetMarkerColor(kGreen+4);
    effZ_mult[0]->SetMarkerSize(1);
    effZ_mult[0]->SetMarkerStyle(22);

    recoPhi->GetAxis(6)->SetRangeUser(20.0, 50.0);
    realPhi->GetAxis(6)->SetRangeUser(20.0, 50.0);
    recoPhi_Z_mult[1] = recoPhi->Projection(4);
    recoPhi_Z_mult[1]->SetName("recoZ_mult2050");
    realPhi_Z_mult[1] = realPhi->Projection(4);
    realPhi_Z_mult[1]->SetName("realZ_mult2050");
    effZ_mult[1] = (TH1D*)recoPhi_Z_mult[1]->Clone("effZ_mult2050");
    effZ_mult[1]->Divide(realPhi_Z_mult[1]);
    effZ_mult[1]->SetTitle("Efficiency vs. Z Vertex for Mult. Bins");
    effZ_mult[1]->SetLineColor(kGreen-2);
    effZ_mult[1]->SetMarkerColor(kGreen-2);
    effZ_mult[1]->SetMarkerSize(1);
    effZ_mult[1]->SetMarkerStyle(22);

    recoPhi->GetAxis(6)->SetRangeUser(50.0, 90.0);
    realPhi->GetAxis(6)->SetRangeUser(50.0, 90.0);
    recoPhi_Z_mult[2] = recoPhi->Projection(4);
    recoPhi_Z_mult[2]->SetName("recoZ_mult5090");
    realPhi_Z_mult[2] = realPhi->Projection(4);
    realPhi_Z_mult[2]->SetName("realZ_mult5090");
    effZ_mult[2] = (TH1D*)recoPhi_Z_mult[2]->Clone("effZ_mult5090");
    effZ_mult[2]->Divide(realPhi_Z_mult[2]);
    effZ_mult[2]->SetTitle("Efficiency vs. Z Vertex for Mult. Bins");
    effZ_mult[2]->SetLineColor(kSpring+5);
    effZ_mult[2]->SetMarkerColor(kSpring+5);
    effZ_mult[2]->SetMarkerSize(1);
    effZ_mult[2]->SetMarkerStyle(22);

    recoPhi->GetAxis(6)->SetRangeUser(0.0, 90.0);
    realPhi->GetAxis(6)->SetRangeUser(0.0, 90.0);
    recoPhi_Z_mult[3] = recoPhi->Projection(4);
    recoPhi_Z_mult[3]->SetName("recoZ_mult090");
    realPhi_Z_mult[3] = realPhi->Projection(4);
    realPhi_Z_mult[3]->SetName("realZ_mult090");
    effZ_mult[3] = (TH1D*)recoPhi_Z_mult[3]->Clone("effZ_mult090");
    effZ_mult[3]->Divide(realPhi_Z_mult[3]);
    effZ_mult[3]->SetTitle("Efficiency vs. Z Vertex for Mult. Bins");
    effZ_mult[3]->SetLineColor(kRed+1);
    effZ_mult[3]->SetMarkerColor(kRed+1);
    effZ_mult[3]->SetMarkerSize(1);
    effZ_mult[3]->SetMarkerStyle(22);

    TCanvas* ceffZ = new TCanvas("ceffZ", "ceffZmult", 55, 55, 600, 600);
    ceffZ->cd();
    effZ_mult[0]->Draw("P");
    effZ_mult[1]->Draw("P SAME");
    effZ_mult[2]->Draw("P SAME");
    effZ_mult[3]->Draw("P SAME");
       
    for(int i=0; i<3; i++){
        ratioZ_mult[i] = (TH1D*)effZ_mult[i]->Clone(Form("ratioZ_mult%i", i));
        ratioZ_mult[i]->Divide(effZ_mult[3]);
    }
    TCanvas* cratioZ = new TCanvas("cratioZ", "cratioZ", 55, 55, 600, 600);
    cratioZ->cd();
    ratioZ_mult[0]->Draw("P");
    ratioZ_mult[1]->Draw("P SAME");
    ratioZ_mult[2]->Draw("P SAME");

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
