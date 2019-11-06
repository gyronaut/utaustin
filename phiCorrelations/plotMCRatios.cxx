


void plotMCRatios(TString inputname){
    TFile *inputfile = new TFile(inputname.Data());

    TList* inlist = (TList*)inputfile->Get("hhCorr_mult_0_100_");

    THnSparseF* hdist = (THnSparseF*)inlist->FindObject("fTrueHDist");
    hdist->Sumw2();
    THnSparseF* primhdist = (THnSparseF*)inlist->FindObject("fTruePrimHDist");
    THnSparseF* sechdist = (THnSparseF*)inlist->FindObject("fTrueSecHDist");

    THnSparseF* phidist = (THnSparseF*)inlist->FindObject("fTruePhiDist");
    phidist->Sumw2();

    //set hadron axis ranges and project onto pT
    //hdist->GetAxis(0)->SetRange(hdist->GetAxis(0)->FindBin(2.001), hdist->GetAxis(0)->FindBin(3.999));
    hdist->GetAxis(1)->SetRange(-1, hdist->GetAxis(1)->GetNbins()+1);
    hdist->GetAxis(2)->SetRange(hdist->GetAxis(2)->FindBin(-0.79999), hdist->GetAxis(2)->FindBin(0.79999));
    
    TH1D* hpT = (TH1D*)hdist->Projection(0);
    hpT->SetLineColor(kBlue+1);
    hpT->SetLineWidth(2);

    primhdist->GetAxis(1)->SetRange(-1, primhdist->GetAxis(1)->GetNbins()+1);
    primhdist->GetAxis(2)->SetRange(primhdist->GetAxis(2)->FindBin(-0.79999), primhdist->GetAxis(2)->FindBin(0.79999));
    
    TH1D* primhpT = (TH1D*)primhdist->Projection(0);
    primhpT->SetLineColor(kRed+1);
    primhpT->SetLineWidth(2);

    sechdist->GetAxis(1)->SetRange(-1, sechdist->GetAxis(1)->GetNbins()+1);
    sechdist->GetAxis(2)->SetRange(sechdist->GetAxis(2)->FindBin(-0.79999), sechdist->GetAxis(2)->FindBin(0.79999));
    
    TH1D* sechpT = (TH1D*)sechdist->Projection(0);
    sechpT->SetLineColor(kGreen+2);
    sechpT->SetLineWidth(2);

    //set phi axis ranges and project onto pT
    phidist->GetAxis(1)->SetRange(-1, phidist->GetAxis(1)->GetNbins()+1);
    phidist->GetAxis(2)->SetRange(-1, phidist->GetAxis(2)->GetNbins()+1);
    phidist->GetAxis(3)->SetRange(phidist->GetAxis(3)->FindBin(-0.79999), phidist->GetAxis(3)->FindBin(0.79999));
    phidist->GetAxis(4)->SetRange(-1, phidist->GetAxis(4)->GetNbins()+1);

    TH1D* phipT = (TH1D*)phidist->Projection(0);
    phipT->Scale(1.0/0.49);

    //Create ratio histogram as function of pT
    TH1D* pTratio = (TH1D*)phipT->Clone("pTratio");
    pTratio->SetTitle("Ratio of #phi/h vs. p_{T}");
    pTratio->Sumw2();
    pTratio->GetXaxis()->SetRangeUser(0.1, 10.0);
    hpT->GetXaxis()->SetRangeUser(0.1, 10.0);
    Bool_t check = pTratio->Divide(hpT);

    printf("division returned %d\n", check);

    //Calculate total ratio and ratio in 2-4 GeV/c pT region
    
    Double_t totalh = hpT->Integral(-1, hpT->GetXaxis()->GetNbins()+1);
    Double_t totalphi = phipT->Integral(-1, phipT->GetXaxis()->GetNbins()+1);
    Double_t midh = hpT->Integral(hpT->GetXaxis()->FindBin(2.00001), hpT->GetXaxis()->FindBin(3.9999));
    Double_t midprimh = primhpT->Integral(primhpT->GetXaxis()->FindBin(2.00001), primhpT->GetXaxis()->FindBin(3.9999));
    Double_t midsech = sechpT->Integral(sechpT->GetXaxis()->FindBin(2.00001), primhpT->GetXaxis()->FindBin(3.9999));
    Double_t midphi = phipT->Integral(phipT->GetXaxis()->FindBin(2.00001), phipT->GetXaxis()->FindBin(3.9999));

    Double_t totratio = totalphi/totalh;
    Double_t midratio = midphi/midh;

    printf("total h/phi ratio: %f, 2-4 GeV/c ratio: %f\n", totratio, midratio);
    printf("ratio of prim to tot in 2-4 GeV/c range: %f\n", midprimh/midh);
    printf("ratio of sec to tot in 2-4 GeV/c range: %f\n", midsech/midh);

    //plot the pT spectra
    TCanvas* chpT = new TCanvas("chpT", "chpT", 50, 50, 600, 600);
    chpT->cd();
    hpT->Draw("HIST");
    primhpT->Draw("SAME HIST");
    sechpT->Draw("SAME HIST");

    TCanvas* cphipT = new TCanvas("cphipT", "cphipT", 50, 50, 600, 600);
    cphipT->cd();
    phipT->Draw("HIST");

    //plot the ratio vs. pT
    TCanvas* cratiopT = new TCanvas("cratiopT", "cratiopT", 50, 50, 600, 600);
    cratiopT->cd();
    pTratio->Draw("HIST"); 

}
