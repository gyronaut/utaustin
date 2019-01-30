void checkMCClosure(){
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


    TFile* fileMC = new TFile("~/phiStudies/MCclosure/US_syst_trig_4_8_assoc_2_4_effcorr_hPhi_0_100.root");

    TH2D* MC2D = (TH2D*)fileMC->Get("AvgUSsubhPhi2Dpeak");
    //correct for wrong MC PID that was included accidentally:
    MC2D->Scale(1.0/0.86558);
    //scale for invariant mass range
    MC2D->Scale(1.0/0.803);
    TH1D* MCDphi = (TH1D*)MC2D->ProjectionY("MCDphi", MC2D->GetXaxis()->FindBin(-1.2), MC2D->GetXaxis()->FindBin(1.2));

    TCanvas* ctrue = new TCanvas("ctrue", "ctrue", 50, 50, 600, 600);
    ctrue->cd();
    true2D->Draw("surf1");

    TCanvas* cMC = new TCanvas("cMC", "cMC", 50, 50, 600, 600);
    cMC->cd();
    MC2D->Draw("surf1");

    TCanvas* c1D = new TCanvas("c1D", "c1D", 50, 50, 1000, 600);
    c1D->Divide(2,1);
    c1D->cd(1);
    trueDphi->Draw();
    c1D->cd(2);
    MCDphi->Draw();


    TH1D* ratioDphi = (TH1D*)MCDphi->Clone("ratioDphi");
    ratioDphi->Divide(trueDphi);
    TCanvas* cratio = new TCanvas ("cratio", "cratio", 50, 50, 600, 600);
    cratio->cd();
    ratioDphi->Draw();

    TH2D* ratio2D = (TH2D*)MC2D->Clone("ratio2D");
    ratio2D->Divide(true2D);
    TCanvas *c2Dratio = new TCanvas("c2Dratio", "c2Dratio", 50, 50, 600, 600);
    c2Dratio->cd();
    ratio2D->Draw("colz");


}
