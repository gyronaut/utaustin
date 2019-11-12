void checkMCClosureEtacut(){
    // true h-phi results //

    gStyle->SetOptStat(0);

    TFile* fileTrue = new TFile("~/phiStudies/MCclosure/trueResults_etacut_02_07.root");
    TList* listTrue = (TList*)fileTrue->Get("truePhiCorr_mult_0_100");

    THnSparseF* trueSparse = 0x0;
    for(int i =0; i< 10; i++){
        THnSparseF* bufSparse = (THnSparseF*)listTrue->FindObject(Form("fDphiTrueHPhiz%i", i));
        if(i == 0){
            trueSparse = (THnSparseF*)bufSparse->Clone("trueSparse");
        }else{
            trueSparse->Add(bufSparse);
        }
    }

    THnSparseF* trigDistSparse = (THnSparseF*)listTrue->FindObject("fTrigDist");
    TH3D* trigDist = (TH3D*)trigDistSparse->Projection(0, 1, 2);
    float totalTrig = (float)trigDist->Integral(trigDist->GetXaxis()->FindBin(4.0), trigDist->GetXaxis()->FindBin(8.0), 1, trigDist->GetYaxis()->GetNbins(), 1, trigDist->GetZaxis()->GetNbins());

    trueSparse->GetAxis(0)->SetRangeUser(4.0, 8.0);
    trueSparse->GetAxis(1)->SetRangeUser(2.0, 4.0);
    trueSparse->GetAxis(4)->SetRangeUser(1, trueSparse->GetAxis(4)->GetNbins());

    TH3D* true3D = (TH3D*)trueSparse->Projection(2, 3, 4);
    true3D->Sumw2();
    //true3D->GetZaxis()->SetRangeUser(1.014, 1.026);
    TH2D* true2D = (TH2D*)true3D->Project3D("xye");
    true2D->Sumw2();
    true2D->Scale(1.0/totalTrig);
    true2D->GetXaxis()->SetRangeUser(-1.0, 1.0);

    TH1D* trueDphi = (TH1D*)true2D->ProjectionY("trueDphi", true2D->GetXaxis()->FindBin(-1.0), true2D->GetXaxis()->FindBin(1.0));
   
    //TH1D* trueDeta = (TH1D*)true2D->ProjectionX("trueDeta", 1, true2D->GetXaxis()->GetNbins());
    TH1D* trueDeta = (TH1D*)true2D->ProjectionX("trueDeta", 1, 8);


    // true h-phi with acceptance //
    
    TFile* fileAccpt = new TFile("~/phiStudies/MCclosure/trig_4_8_assoc_2_4_MC_hPhi_0_100_02_07.root");
    TH2D* accpt2D = (TH2D*)fileAccpt->Get("MChPhi2D");
    TH1D* accptDphi = (TH1D*)accpt2D->ProjectionY("kaonDphi", accpt2D->GetXaxis()->FindBin(-1.0), accpt2D->GetXaxis()->FindBin(1.0));
    accptDphi->GetXaxis()->SetTitle("#Delta#varphi");
    //TH1D* accptDeta = (TH1D*)accpt2D->ProjectionX("accptDeta", 1, accpt2D->GetXaxis()->GetNbins());
    TH1D* accptDeta = (TH1D*)accpt2D->ProjectionX("accptDeta", 1, 8);



    // correlation with real kaons (no efficiency) //
    TFile* fileMC = new TFile("~/phiStudies/MCclosure/US_syst_trig_4_8_assoc_2_4_kaons_hPhi_0_100_02_07.root");

    TH2D* MC2D = (TH2D*)fileMC->Get("AvgUSsubhPhi2Dpeak");
    //scale for invariant mass range
    //MC2D->Scale(1.0/0.803); // percent from Data
    MC2D->Scale(1.0/0.8470); // percent from MC
    TH1D* MCDphi = (TH1D*)MC2D->ProjectionY("trueDphi", MC2D->GetXaxis()->FindBin(-1.0), MC2D->GetXaxis()->FindBin(1.0));
    MCDphi->GetXaxis()->SetTitle("#Delta#varphi");
    //TH1D* MCDeta = (TH1D*)MC2D->ProjectionX("MCDeta", 1, MC2D->GetXaxis()->GetNbins());
    TH1D* MCDeta = (TH1D*)MC2D->ProjectionX("MCDeta", 1, 8);


    // full correlation (efficiency included) //
    //TFile* fileMCfull = new TFile("~/phiStudies/MCclosure/US_syst_trig_4_8_assoc_2_4_etacut_hPhi_0_100_02_06.root");
    TFile* fileMCfull = new TFile("~/phiStudies/MCclosure/US_syst_trig_4_8_assoc_2_4_tpc80realeff_hPhi_0_100_02_08.root");
    //TFile* fileMCfull = new TFile("~/phiStudies/MCclosure/US_syst_trig_4_8_assoc_2_4_tpc80eff_hPhi_0_100_02_20.root");
    //TFile* fileMCfull = new TFile("~/phiStudies/MCclosure/US_syst_trig_4_8_assoc_2_4_ktrack_hPhi_0_100_02_11.root");
    //TFile* fileMCfull = new TFile("~/phiStudies/MCclosure/US_syst_trig_4_8_assoc_2_4_effcorr_hPhi_0_100_fulleff_02_26.root");
    //TFile* fileMCfull = new TFile("~/phiStudies/MCclosure/US_syst_trig_4_8_assoc_2_4_effcorr_hPhi_0_100_CENT_02_28.root");

    TH2D* MCfull2D = (TH2D*)fileMCfull->Get("AvgUSsubhPhi2Dpeak");
    //correct for wrong MC PID that was included accidentally:
    MCfull2D->Scale(1.0/0.86558);
    //scale for invariant mass range
    //MCfull2D->Scale(1.0/0.803); // percent from data
    MCfull2D->Scale(1.0/0.8470); // percent from MC
    //MCfull2D->Scale(1.081) //account for difference of CENT and FAST?
    TH1D* MCfullDphi = (TH1D*)MCfull2D->ProjectionY("MCDphi", MCfull2D->GetXaxis()->FindBin(-1.0), MCfull2D->GetXaxis()->FindBin(1.0));
    MCfullDphi->GetXaxis()->SetTitle("#Delta#varphi");

    //TFile* fileMCnoeff = new TFile("~/phiStudies/MCclosure/US_syst_trig_4_8_assoc_2_4_tpc80noeff_hPhi_0_100_02_08.root");
    TFile* fileMCnoeff = new TFile("~/phiStudies/MCclosure/US_syst_trig_4_8_assoc_2_4_effcorr_hPhi_0_100_fulleff_02_26.root");
    TH2D* MCnoeff2D = (TH2D*)fileMCnoeff->Get("AvgUSsubhPhi2Dpeak");
    //correct for wrong MC PID that was included accidentally:
    MCnoeff2D->Scale(1.0/0.86558);
    //scale for invariant mass range
    //MCfull2D->Scale(1.0/0.803); // percent from data
    //MCnoeff2D->Scale(1.0/0.8470); // percent from MC
    TH1D* MCnoeffDphi = (TH1D*)MCnoeff2D->ProjectionY("MCnoeffDphi", MCfull2D->GetXaxis()->FindBin(-1.0), MCfull2D->GetXaxis()->FindBin(1.0));
    MCnoeffDphi->GetXaxis()->SetTitle("#Delta#varphi");



    //TFile* fileMCCENT = new TFile("~/phiStudies/MCclosure/US_syst_trig_4_8_assoc_2_4_effcorr_hPhi_0_100_CENT_02_28.root");
    TFile* fileMCCENT = new TFile("~/phiStudies/MCclosure/US_syst_trig_4_8_assoc_2_4_effcorr_hPhi_0_100_CENT_03_07.root");
    TH2D* MCCENT2D = (TH2D*)fileMCCENT->Get("AvgUSsubhPhi2Dpeak");
    //correct for wrong MC PID that was included accidentally:
    MCCENT2D->Scale(1.0/0.86558);
    //MCCENT2D->Scale(1/1.081);
    TH1D* MCCENTDphi = (TH1D*)MCnoeff2D->ProjectionY("MCCENTDphi", MCfull2D->GetXaxis()->FindBin(-1.0), MCfull2D->GetXaxis()->FindBin(1.0));
    MCCENTDphi->GetXaxis()->SetTitle("#Delta#varphi");   

    TFile* fileCENT = new TFile("~/phiStudies/MCclosure/trig_4_8_assoc_2_4_MC_hPhi_0_100_CENT_02_28.root");
    TH2D* CENT2D = (TH2D*)fileCENT->Get("MChPhi2D");
    TH1D* CENTDphi = (TH1D*)CENT2D->ProjectionY("CENTDphi", CENT2D->GetXaxis()->FindBin(-1.0), CENT2D->GetXaxis()->FindBin(1.0));
    CENTDphi->GetXaxis()->SetTitle("#Delta#varphi");
    //TH1D* CENTDeta = (TH1D*)CENT2D->ProjectionX("CENTDeta", 1, CENT2D->GetXaxis()->GetNbins());
    TH1D* CENTDeta = (TH1D*)CENT2D->ProjectionX("CENTDeta", 1, 8);

    TFile* fileCombined = new TFile("~/phiStudies/MCclosure/true_centfast_02_28.root");
    TH2D* Combined2D = (TH2D*)fileCombined->Get("CombinedhPhi2D");
    TH1D* CombinedDphi = (TH1D*)Combined2D->ProjectionY("CombinedDphi", Combined2D->GetXaxis()->FindBin(-1.0), Combined2D->GetXaxis()->FindBin(1.0));
    CombinedDphi->GetXaxis()->SetTitle("#Delta#varphi");
    //TH1D* CombinedDeta = (TH1D*)Combined2D->ProjectionX("CombinedDeta", 1, Combined2D->GetXaxis()->GetNbins());
    TH1D* CombinedDeta = (TH1D*)Combined2D->ProjectionX("CombinedDeta", 1, 8);

    TFile* fileMCCombined = new TFile("~/phiStudies/MCclosure/US_syst_trig_4_8_assoc_2_4_effcorr_hPhi_0_100_FASTCENT_03_07.root");
    TH2D* MCCombined2D = (TH2D*)fileMCCombined->Get("AvgUSsubhPhi2Dpeak");
    //correct for wrong MC PID that was included accidentally:
    MCCombined2D->Scale(1.0/0.86558);
    //MCCombined2D->Scale(1/1.081);
    TH1D* MCCombinedDphi = (TH1D*)MCnoeff2D->ProjectionY("MCCombinedDphi", MCfull2D->GetXaxis()->FindBin(-1.0), MCfull2D->GetXaxis()->FindBin(1.0));
    MCCombinedDphi->GetXaxis()->SetTitle("#Delta#varphi");



    TCanvas* ctrue = new TCanvas("ctrue", "ctrue", 50, 50, 600, 600);
    ctrue->cd();
    true2D->Draw("surf1");

    TCanvas* caccpt = new TCanvas("caccpt", "caccpt", 50, 50, 600, 600);
    caccpt->cd();
    accpt2D->Draw("surf1");

    TCanvas* cMC = new TCanvas("cMC", "cMC", 50, 50, 600, 600);
    cMC->cd();
    MC2D->Draw("surf1");
    //accpt2D->Draw("surf1");

    TCanvas* c1D = new TCanvas("c1D", "c1D", 50, 50, 1000, 400);
    c1D->Divide(4,1);
    c1D->cd(1);
    //trueDphi->GetYaxis()->SetRangeUser(0.0005, 0.0030);
    trueDphi->Draw();
    c1D->cd(2);
    //accptDphi->GetYaxis()->SetRangeUser(0.0005, 0.0030);
    accptDphi->Draw();
    c1D->cd(3);
    //MCDphi->GetYaxis()->SetRangeUser(0.0005, 0.0030);
    MCDphi->Draw();
    c1D->cd(4);
    //MCfullDphi->GetYaxis()->SetRangeUser(0.0005, 0.0030);
    MCfullDphi->Draw();

    TCanvas* c1Deta = new TCanvas("c1Deta", "c1Deta", 50, 50, 1000, 400);
    c1Deta->Divide(3,1);
    c1Deta->cd(1);
    trueDeta->Draw();
    c1Deta->cd(2);
    accptDeta->Draw();
    c1Deta->cd(3);
    MCDeta->Draw();


    TH1D* ratiofullDphi = (TH1D*)MCfullDphi->Clone("ratiofullDphi");
    ratiofullDphi->Divide(MCDphi);
    TCanvas* cfullratio = new TCanvas ("cratioFinal2Kaons", "cratioFinal2Kaons", 50, 50, 600, 600);
    cfullratio->cd();
    ratiofullDphi->Draw();

    TH1D* ratioDphi = (TH1D*)MCDphi->Clone("ratioDphi");
    ratioDphi->Divide(accptDphi);
    ratioDphi->SetTitle("Ratio of (h-KK_{MC})/(h-#phi_{MC})");
    TF1* finalFit2 = new TF1("finalFit2", "pol0(0)", -0.5*TMath::Pi(), 1.5*TMath::Pi());
    ratioDphi->Fit(finalFit2);
    TCanvas* cratio = new TCanvas ("cratioKaons2AccptPhi", "cratioKaons2AccptPhi", 50, 50, 600, 600);
    cratio->cd();
    ratioDphi->Draw();

    TH1D *hint2 = new TH1D("hint2","hint2", 100, -0.5*TMath::Pi(), 1.5*TMath::Pi());
    (TVirtualFitter::GetFitter())->GetConfidenceIntervals(hint2);
    hint2->SetFillColor(2);
    hint2->SetFillStyle(3001);
    hint2->Draw("e3 same");


    TH2D* ratio2D = (TH2D*)MC2D->Clone("ratio2D");
    ratio2D->Divide(accpt2D);
    TCanvas *c2Dratio = new TCanvas("c2Dratio", "c2Dratio", 50, 50, 600, 600);
    c2Dratio->cd();
    ratio2D->Draw("colz");

    TH1D* ratioTrueDphi = (TH1D*)accptDphi->Clone("ratioTrueDphi");
    ratioTrueDphi->Divide(trueDphi);
    TCanvas* ctrueratio = new TCanvas ("cratioTrue2Accpt", "cratioTrue2Accpt", 50, 50, 600, 600);
    ctrueratio->cd();
    ratioTrueDphi->Draw();

    TH2D* ratioTrue2D = (TH2D*)accpt2D->Clone("ratioTrue2D");
    ratioTrue2D->Divide(true2D);
    TCanvas *ctrue2Dratio = new TCanvas("ctrue2Dratio", "ctrue2Dratio", 50, 50, 600, 600);
    ctrue2Dratio->cd();
    ratioTrue2D->Draw("colz");

    TH1D* ratioFinal = (TH1D*)MCfullDphi->Clone("ratioFinalDphi");
    ratioFinal->Divide(trueDphi);
    TCanvas *cfinalratio = new TCanvas("cratioFinal2True", "cratioFinal2True", 50, 50, 600, 600);
    cfinalratio->cd();
    ratioFinal->Draw();

    TH1D* ratioFinal2 = (TH1D*)MCfullDphi->Clone("ratioFinal2Dphi");
    ratioFinal2->Divide(accptDphi);
    ratioFinal2->SetTitle("Ratio of (h-KK_{track})/(h-#phi_{MC})");
    TF1* finalFit = new TF1("finalFit", "pol0(0)", -0.5*TMath::Pi(), 1.5*TMath::Pi());
    ratioFinal2->Fit(finalFit);

    TCanvas *cfinal2ratio = new TCanvas("cratioFinal2Accpt", "cratioFinal2Accpt", 50, 50, 600, 600);
    cfinal2ratio->cd();
    ratioFinal2->Draw();
    TH1D *hint = new TH1D("hint","hint", 100, -0.5*TMath::Pi(), 1.5*TMath::Pi());
    (TVirtualFitter::GetFitter())->GetConfidenceIntervals(hint);
    hint->SetFillColor(2);
    hint->SetFillStyle(3001);
    hint->Draw("e3 same");


    TCanvas *contop = new TCanvas("cFinalKaonsOnTop", "cFinalKaonsOnTop", 50, 50, 600, 600);
    contop->cd();
    MCfullDphi->SetLineColor(kBlue+1);
    MCfullDphi->SetLineWidth(2);
    MCfullDphi->Draw("P E");
    MCDphi->SetLineColor(kRed+1);
    MCDphi->SetLineWidth(2);
    MCDphi->Draw("P E SAME");
 
    TCanvas *contop2 = new TCanvas("cKaonsAccptOnTop", "cKaonsAccptOnTop", 50, 50, 600, 600);
    contop2->cd();
    accptDphi->SetLineColor(kGreen+2);
    accptDphi->SetLineWidth(2);
    accptDphi->Draw("P E");
    MCDphi->Draw("P E SAME");   

    TCanvas *contop3 = new TCanvas("cFinalAccptOnTop", "cFinalAccptOnTop", 50, 50, 600, 600);
    contop3->cd();
    accptDphi->Draw("P E");
    MCfullDphi->Draw("P E SAME");

    TCanvas *contop4 = new TCanvas("cFinalNoEffOnTop", "cFinalNoeffOnTop", 50, 50, 600, 600);
    contop4->cd();
    MCnoeffDphi->Draw("P E");
    MCfullDphi->Draw("P E SAME");

    TH1D* rationoeff = (TH1D*)MCfullDphi->Clone("rationoeff");
    rationoeff->Divide(MCnoeffDphi);
    TCanvas* crationoeff = new TCanvas("crationoeff", "crationoeff", 60, 60, 600, 600);
    crationoeff->cd();
    rationoeff->Draw();

    TCanvas *contop5 = new TCanvas("cCENTFinalAccptOnTop", "cCENTFinalAccptOnTop", 50, 50, 600, 600);
    contop5->cd();
    CENTDphi->Draw("P E");
    MCCENTDphi->Draw("P E SAME");

    TH1D* ratioFinalCENT = (TH1D*)MCCENTDphi->Clone("ratioFinalCENTDphi");
    ratioFinalCENT->Divide(CENTDphi);
    TF1* finalCENT = new TF1("finalCENT", "pol0(0)", -0.5*TMath::Pi(), 1.5*TMath::Pi());

    TCanvas *cfinalCENTratio = new TCanvas("cratioCENTFinal2True", "cratioCENTFinal2True", 50, 50, 600, 600);
    cfinalCENTratio->cd();

    ratioFinalCENT->Fit(finalCENT);
    TH1D *hintCENT = new TH1D("hintCENT","hintCENT", 100, -0.5*TMath::Pi(), 1.5*TMath::Pi());
    (TVirtualFitter::GetFitter())->GetConfidenceIntervals(hintCENT);
    hintCENT->SetFillColor(2);
    hintCENT->SetFillStyle(3001);
    hintCENT->Draw("e3 same");
    ratioFinalCENT->Draw("SAME");

    TF1* centfit = new TF1("centfit", "pol0(0)", -0.5*TMath::Pi(), 1.5*TMath::Pi());
    TH1D* ratioFASTCENT = (TH1D*)CENTDphi->Clone("ratioFASTCENTDphi");
    ratioFASTCENT->Divide(accptDphi);
    TCanvas *cFASTCENTratio = new TCanvas("cratioCENT2FAST", "cratioCENT2FAST", 50, 50, 600, 600);
    cFASTCENTratio->cd();
    ratioFASTCENT->Fit(centfit);
    ratioFASTCENT->Draw();

    TCanvas *contop6 = new TCanvas("cCENTFASTOnTop", "cCENTFASTOnTop", 50, 50, 600, 600);
    contop6->cd();
    CENTDphi->SetLineColor(kMagenta+1);
    CENTDphi->Draw("P E");
    accptDphi->Draw("SAME P E");
    CombinedDphi->SetLineColor(kViolet+2);
    CombinedDphi->Draw("P E SAME");

    TCanvas *contop7 = new TCanvas("cCombinedOnTop", "cCombinedOnTop", 50, 50, 600, 600);
    contop7->cd();
    CombinedDphi->SetLineColor(kGreen+2);
    CombinedDphi->SetLineWidth(2);
    CombinedDphi->SetTitle("(h-#phi_{MC}) vs. (h-KK_{track})");
    CombinedDphi->Draw("P E SAME");
    MCCombinedDphi->SetLineColor(kBlue+1);
    MCCombinedDphi->SetLineWidth(2);
    MCCombinedDphi->Draw("P E SAME");

    TH1D* ratioFinalCombined = (TH1D*)MCCombinedDphi->Clone("ratioFinalCombinedDphi");
    ratioFinalCombined->Divide(CombinedDphi);
    TF1* finalCombined = new TF1("finalCombined", "pol0(0)", -0.5*TMath::Pi(), 1.5*TMath::Pi());

    TCanvas *cCombinedRatio = new TCanvas("cratioCombined", "cratioCombined", 50, 50, 600, 600);
    cCombinedRatio->cd();
    ratioFinalCombined->Fit(finalCombined);
    ratioFinalCombined->SetTitle("Ratio of (h-KK_{track})/(h-#phi_{MC})");
    ratioFinalCombined->Draw();
    TH1D *hintcomb = new TH1D("hintcomb","hintcomb", 100, -0.5*TMath::Pi(), 1.5*TMath::Pi());
    (TVirtualFitter::GetFitter())->GetConfidenceIntervals(hintcomb);
    hintcomb->SetFillColor(2);
    hintcomb->SetFillStyle(3001);
    hintcomb->Draw("e3 same");





}
