TH1D* projectPT(THnSparseF* hist, TString outname, Int_t particle = -1, Bool_t withRapCut = kFALSE){
    Int_t ndims = hist->GetNdimensions();
    for(int i = 0; i < ndims; i++){
        hist->GetAxis(i)->SetRange(-1, hist->GetAxis(i)->GetNbins()+1);
    }
    Int_t etaaxis = 2;
    if(particle == -99) etaaxis=3;
    hist->GetAxis(etaaxis)->SetRange(hist->GetAxis(etaaxis)->FindBin(-0.79999), hist->GetAxis(etaaxis)->FindBin(0.79999));

    if(withRapCut){
        Int_t rapaxis = 3;
        if(particle == -99) rapaxis = 4;
        hist->GetAxis(rapaxis)->SetRange(hist->GetAxis(rapaxis)->FindBin(-0.4999), hist->GetAxis(rapaxis)->FindBin(0.4999));
    }

    if(particle > 0){
        hist->GetAxis(4)->SetRange(particle, particle);
    }

    TH1D* ptdist = (TH1D*)hist->Projection(0);
    ptdist->SetName(outname.Data());
    ptdist->SetTitle(outname.Data());
    ptdist->SetLineWidth(2);
    return ptdist;
}

void plotMCRatios(TString inputname){
    TFile *inputfile = new TFile(inputname.Data());

    TList* inlist = (TList*)inputfile->Get("hhCorr_mult_0_100_");

    THnSparseF* hdist = (THnSparseF*)inlist->FindObject("fTrueHDist");
    THnSparseF* primhdist = (THnSparseF*)inlist->FindObject("fTruePrimHDist");
    THnSparseF* notprimhdist = (THnSparseF*)inlist->FindObject("fTrueNotPrimHDist");
    THnSparseF* sechdist = (THnSparseF*)inlist->FindObject("fTrueSecHDist");
    THnSparseF* phidist = (THnSparseF*)inlist->FindObject("fTruePhiDist");

    THnSparseF* trighdist = (THnSparseF*)inlist->FindObject("fTriggeredTrueHDist");
    THnSparseF* trigprimhdist = (THnSparseF*)inlist->FindObject("fTriggeredTruePrimHDist");
    THnSparseF* trignotprimhdist = (THnSparseF*)inlist->FindObject("fTriggeredTrueNotPrimHDist");
    THnSparseF* trigsechdist = (THnSparseF*)inlist->FindObject("fTriggeredTrueSecHDist");
    THnSparseF* trigphidist = (THnSparseF*)inlist->FindObject("fTriggeredTruePhiDist");

    TH1D* hpT = projectPT(hdist, "hadronPT");
    hpT->SetLineColor(kBlue+1);

    TH1D* trigprimhpT = projectPT(trigprimhdist, "trigPrimHadronPT");
    trigprimhpT->SetLineColor(kAzure+1);
    
    TH1D* primhpT = projectPT(primhdist, "primaryPT");
    primhpT->SetLineColor(kRed+1);

    TH1D* primRaphpT = projectPT(primhdist, "primaryRapPT", -1, kTRUE);
    primRaphpT->SetLineColor(kRed+2);

    TH1D* primpipT = projectPT(primhdist, "primaryPionPT", 1);
    primpipT->SetLineColor(kViolet+1);

    TH1D* primKpT = projectPT(primhdist, "primaryKPT", 2);
    primKpT->SetLineColor(kBlue+1);

    TH1D* primppT = projectPT(primhdist, "primaryPPT", 3);
    primppT->SetLineColor(kOrange+2);

    TH1D* primepT = projectPT(primhdist, "primaryePT", 4);
    primepT->SetLineColor(kGreen-4);

    TH1D* primmupT = projectPT(primhdist, "primarymuPT", 5);
    primmupT->SetLineColor(kSpring+1);

    TH1D* sechpT = projectPT(sechdist, "sechPT");
    sechpT->SetLineColor(kGreen+2);

    // project phi distribution and scale for branching ratio
    TH1D* phipT = projectPT(phidist, "phiPT", -99, kFALSE);
    phipT->Scale(1.0/0.49);
    phipT->SetLineColor(kRed-4);

    TH1D* phiRappT = projectPT(phidist, "phiRapPT", -99, kTRUE);
    phiRappT->Scale(1.0/0.49);
    phiRappT->SetLineColor(kRed-3);

    TH1D* trigphipT = projectPT(trigphidist, "trigphiPT");
    trigphipT->Scale(1.0/0.49);
    trigphipT->SetLineColor(kAzure+2);

    //Create ratio histogram as function of pT
    TH1D* pTratio = (TH1D*)phipT->Clone("pTratio");
    pTratio->SetTitle("Ratio of #phi/h vs. p_{T}");
    pTratio->Sumw2();
    pTratio->GetXaxis()->SetRangeUser(0.1, 10.0);
    primhpT->GetXaxis()->SetRangeUser(0.1, 10.0);
    Bool_t check = pTratio->Divide(primhpT);

    //Do same ratio histogram for triggered case
    TH1D* trigpTratio = (TH1D*)trigphipT->Clone("trigpTratio");
    trigpTratio->SetTitle("Ratio of triggered #phi/h vs. p_{T}");
    trigpTratio->Sumw2();
    trigpTratio->GetXaxis()->SetRangeUser(0.1, 10.0);
    trigprimhpT->GetXaxis()->SetRangeUser(0.1, 10.0);
    trigpTratio->Divide(trigprimhpT);
    trigpTratio->SetLineColor(kViolet+2);

    printf("division returned %d\n", check);

    //Calculate total ratio and ratio in 2-4 GeV/c pT region
    
    Double_t totalh = primhpT->Integral(-1, primhpT->GetXaxis()->GetNbins()+1);
    Double_t totalphi = phipT->Integral(-1, phipT->GetXaxis()->GetNbins()+1);
    Double_t totpi = primpipT->Integral(-1, primpipT->GetXaxis()->GetNbins()+1);
    Double_t totK = primKpT->Integral(-1, primKpT->GetXaxis()->GetNbins()+1);
    Double_t totp = primppT->Integral(-1, primppT->GetXaxis()->GetNbins()+1);
    Double_t tote = primepT->Integral(-1, primepT->GetXaxis()->GetNbins()+1);
    Double_t totmu = primmupT->Integral(-1, primmupT->GetXaxis()->GetNbins()+1);

    Double_t midh = hpT->Integral(hpT->GetXaxis()->FindBin(2.00001), hpT->GetXaxis()->FindBin(3.9999));
    Double_t midprimh = primhpT->Integral(primhpT->GetXaxis()->FindBin(2.00001), primhpT->GetXaxis()->FindBin(3.9999));
    Double_t midrapprimh = primRaphpT->Integral(primRaphpT->GetXaxis()->FindBin(2.00001), primRaphpT->GetXaxis()->FindBin(3.9999));
    Double_t midprimpi = primpipT->Integral(primpipT->GetXaxis()->FindBin(2.00001), primpipT->GetXaxis()->FindBin(3.9999));
    Double_t midprimK = primKpT->Integral(primKpT->GetXaxis()->FindBin(2.00001), primKpT->GetXaxis()->FindBin(3.9999));
    Double_t midprimp = primppT->Integral(primppT->GetXaxis()->FindBin(2.00001), primppT->GetXaxis()->FindBin(3.9999));
    Double_t midprime = primepT->Integral(primepT->GetXaxis()->FindBin(2.00001), primepT->GetXaxis()->FindBin(3.9999));
    Double_t midprimmu = primmupT->Integral(primmupT->GetXaxis()->FindBin(2.00001), primmupT->GetXaxis()->FindBin(3.9999));
    Double_t midsech = sechpT->Integral(sechpT->GetXaxis()->FindBin(2.00001), primhpT->GetXaxis()->FindBin(3.9999));
    Double_t midphi = phipT->Integral(phipT->GetXaxis()->FindBin(2.00001), phipT->GetXaxis()->FindBin(3.9999));
    Double_t midrapphi =phiRappT->Integral(phiRappT->GetXaxis()->FindBin(2.00001), phiRappT->GetXaxis()->FindBin(3.9999));

    Double_t midtrigprimh = trigprimhpT->Integral(trigprimhpT->GetXaxis()->FindBin(2.0001), trigprimhpT->GetXaxis()->FindBin(3.999));
    Double_t midtrigphi = trigphipT->Integral(trigphipT->GetXaxis()->FindBin(2.001), trigphipT->GetXaxis()->FindBin(3.999));


    Double_t totratio = totalphi/totalh;
    Double_t midratio = midphi/midprimh;
    Double_t midrapratio = midrapphi/midrapprimh;
    Double_t midtrigratio = midtrigphi/midtrigprimh;

    printf("total phi/h ratio: %f, 2-4 GeV/c ratio: %f\n", totratio, midratio);
    //printf("ratio of prim to tot in 2-4 GeV/c range: %f\n", midprimh/midh);
    printf("ratio of triggered phi/h in 2-4 GeV/c: %f\n", midtrigratio);
    //printf("ratio of sec to tot in 2-4 GeV/c range: %f\n", midsech/midh);
    printf("ratio of h/phi in mid rapidity, 2-4GeV/c: %f\n", midrapratio);
    printf("pion: %2.2f%%, kaon: %2.2f%%, proton: %2.2f%%, e: %2.2f%%, muon: %2.2f%%\n", 100.0*totpi/totalh, 100.0*totK/totalh, 100.0*totp/totalh, 100.0*tote/totalh, 100.0*totmu/totalh);
    printf("pion: %2.2f%%, kaon: %2.2f%%, proton: %2.2f%%, e: %2.2f%%, muon: %2.2f%%\n", 100.0*midprimpi/midprimh, 100.0*midprimK/midprimh, 100.0*midprimp/midprimh, 100.0*midprime/midprimh, 100.0*midprimmu/midprimh);

    //plot the pT spectra
    TCanvas* chpT = new TCanvas("chpT", "chpT", 50, 50, 600, 600);
    chpT->cd();
    //hpT->Draw("HIST");
    primhpT->Draw("HIST");
    trigprimhpT->Draw("SAME HIST");
    //primpipT->Draw("SAME HIST");
    //primKpT->Draw("SAME HIST");
    //primppT->Draw("SAME HIST");
    //primepT->Draw("SAME HIST");
    //primmupT->Draw("SAME HIST");
    //sechpT->Draw("SAME HIST");

    TCanvas* cphipT = new TCanvas("cphipT", "cphipT", 50, 50, 600, 600);
    cphipT->cd();
    phipT->Draw("HIST");
    trigphipT->Draw("HIST SAME");

    //plot the ratio vs. pT
    TCanvas* cratiopT = new TCanvas("cratiopT", "cratiopT", 50, 50, 600, 600);
    cratiopT->cd();
    pTratio->Draw("HIST"); 

    //plot triggered ratio vs. pT
    TCanvas* ctrigratio = new TCanvas("ctrigratio", "ctrigratio", 50, 50, 600, 600);
    ctrigratio->cd();
    pTratio->Draw("HIST");
    trigpTratio->Draw("SAME HIST");

}
