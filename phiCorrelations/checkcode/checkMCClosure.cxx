void checkMCClosure(){
    // true h-phi results //

    TFile* fileTrue = new TFile("~/phiStudies/MCclosure/trueResults.root");
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
    true2D->GetXaxis()->SetRangeUser(-1.2, 1.2);

    TH1D* trueDphi = (TH1D*)true2D->ProjectionY("trueDphi", true2D->GetXaxis()->FindBin(-1.2), true2D->GetXaxis()->FindBin(1.2));
   
    //TH1D* trueDeta = (TH1D*)true2D->ProjectionX("trueDeta", 1, true2D->GetXaxis()->GetNbins());
    TH1D* trueDeta = (TH1D*)true2D->ProjectionX("trueDeta", 1, 8);


    // true h-phi with acceptance
    
    TFile* fileAccpt = new TFile("~/phiStudies/MCclosure/trig_4_8_assoc_2_4_true_hPhi_0_100.root");
    TH2D* accpt2D = (TH2D*)fileAccpt->Get("hPhi2Dpeak");
    TH1D* accptDphi = (TH1D*)accpt2D->ProjectionY("kaonDphi", accpt2D->GetXaxis()->FindBin(-1.2), accpt2D->GetXaxis()->FindBin(1.2));
    //TH1D* accptDeta = (TH1D*)accpt2D->ProjectionX("accptDeta", 1, accpt2D->GetXaxis()->GetNbins());
    TH1D* accptDeta = (TH1D*)accpt2D->ProjectionX("accptDeta", 1, 8);



    // correlation with real kaons (no efficiency) //
    TFile* fileMC = new TFile("~/phiStudies/MCclosure/US_syst_trig_4_8_assoc_2_4_kaons_hPhi_0_100_02_03.root");

    TH2D* MC2D = (TH2D*)fileMC->Get("AvgUSsubhPhi2Dpeak");
    TH1D* MCDphi = (TH1D*)MC2D->ProjectionY("trueDphi", MC2D->GetXaxis()->FindBin(-1.2), MC2D->GetXaxis()->FindBin(1.2));
    //TH1D* MCDeta = (TH1D*)MC2D->ProjectionX("MCDeta", 1, MC2D->GetXaxis()->GetNbins());
    TH1D* MCDeta = (TH1D*)MC2D->ProjectionX("MCDeta", 1, 8);


    // full correlation (efficiency included) //
    TFile* fileMCfull = new TFile("~/phiStudies/MCclosure/US_syst_trig_4_8_assoc_2_4_effcorr_hPhi_0_100_02_04.root");

    TH2D* MCfull2D = (TH2D*)fileMCfull->Get("AvgUSsubhPhi2Dpeak");
    //correct for wrong MC PID that was included accidentally:
    MCfull2D->Scale(1.0/0.86558);
    //scale for invariant mass range
    MCfull2D->Scale(1.0/0.803);
    TH1D* MCfullDphi = (TH1D*)MCfull2D->ProjectionY("MCDphi", MCfull2D->GetXaxis()->FindBin(-1.2), MCfull2D->GetXaxis()->FindBin(1.2));


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
    trueDphi->Draw();
    c1D->cd(2);
    accptDphi->Draw();
    c1D->cd(3);
    MCDphi->Draw();
    c1D->cd(4);
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
    TCanvas* cfullratio = new TCanvas ("cfullratio", "cfullratio", 50, 50, 600, 600);
    cfullratio->cd();
    ratiofullDphi->Draw();

    TH1D* ratioDphi = (TH1D*)MCDphi->Clone("ratioDphi");
    ratioDphi->Divide(accptDphi);
    TCanvas* cratio = new TCanvas ("cratio", "cratio", 50, 50, 600, 600);
    cratio->cd();
    ratioDphi->Draw();

    TH2D* ratio2D = (TH2D*)MC2D->Clone("ratio2D");
    ratio2D->Divide(accpt2D);
    TCanvas *c2Dratio = new TCanvas("c2Dratio", "c2Dratio", 50, 50, 600, 600);
    c2Dratio->cd();
    ratio2D->Draw("colz");

    TH1D* ratioTrueDphi = (TH1D*)accptDphi->Clone("ratioTrueDphi");
    ratioTrueDphi->Divide(trueDphi);
    TCanvas* ctrueratio = new TCanvas ("ctrueratio", "ctrueratio", 50, 50, 600, 600);
    ctrueratio->cd();
    ratioTrueDphi->Draw();

    TH2D* ratioTrue2D = (TH2D*)accpt2D->Clone("ratioTrue2D");
    ratioTrue2D->Divide(true2D);
    TCanvas *ctrue2Dratio = new TCanvas("ctrue2Dratio", "ctrue2Dratio", 50, 50, 600, 600);
    ctrue2Dratio->cd();
    ratioTrue2D->Draw("colz");


}
