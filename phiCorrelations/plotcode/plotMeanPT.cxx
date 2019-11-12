float MeanPT(TH1D* hist, float low, float high){
    float mean = 0.0;
    float tot = 0.0;
    int lowbin = hist->GetXaxis()->FindBin(low + 0.001);
    int highbin = hist->GetXaxis()->FindBin(high - 0.001);
    tot = hist->Integral(lowbin, highbin);
    for(int ibin = lowbin; ibin <= highbin; ibin++){
        mean+= hist->GetXaxis()->GetBinCenter(ibin)*hist->GetBinContent(ibin);
    }
    mean = mean/tot;

    return mean;
}

TH1D* Normalize(TH1D* hist, TString name, float low, float high){
    float eps = 0.001;
    TH1D* norm = (TH1D*)hist->Clone(name.Data());
    norm->GetXaxis()->SetRangeUser(low+eps, high-eps);
    norm->Scale(1.0/norm->Integral());
    norm->SetLineWidth(2);

    return norm;
}

TH2D* ProjectKKDist(THnSparseF* hist, TString name){

    hist->GetAxis(0)->SetRangeUser(2.01, 3.99);
    hist->GetAxis(1)->SetRangeUser(0, 0);
    hist->GetAxis(2)->SetRange(0, 0);
    hist->GetAxis(3)->SetRange(0, 0);
    hist->GetAxis(4)->SetRange(0, 0);

    TH2D* proj = (TH2D*)hist->Projection(1, 0);
    proj->SetName(name.Data());
    proj->SetTitle(name.Data());

    return proj;
}

TH1D* ProjectMassToPT(TH2D* hist, TString name, Float_t lowmass, Float_t highmass){
    Float_t eps = 0.0001;
    hist->GetYaxis()->SetRangeUser(lowmass+eps, highmass-eps);
    TH1D* proj = (TH1D*)hist->ProjectionX(name.Data());
    proj->SetName(name.Data());
    proj->SetTitle(name.Data());
    
    return proj;
}

TH1D* ReconPhiPT(TH2D* UShist, TH2D* LShist){
    TH1D* pthist;

    return pthist;

}

void plotMeanPT(TString inputname, int multlow, int multhigh){

    TFile *input = new TFile(inputname.Data());

    TString liststr(Form("phiCorr_mult_%d_%d_", multlow, multhigh));
    TList* list = (TList*)input->Get(liststr.Data());

    TH1D* hadronPT = (TH1D*)list->FindObject("fHadronPT");
    TH1D* hadronTrigPT = (TH1D*)list->FindObject("fHadronTrigPT");
    TH1D* hadronTrigPhiPT = (TH1D*)list->FindObject("fHadronTrigPhiPT");
    THnSparseF* kkUSDist= (THnSparseF*)list->FindObject("fkkUSDist");
    THnSparseF* kkLSDist= (THnSparseF*)list->FindObject("fkkLSDist");
    THnSparseF* kkUSTrigDist= (THnSparseF*)list->FindObject("fkkUSTrigDist");
    THnSparseF* kkLSTrigDist= (THnSparseF*)list->FindObject("fkkLSTrigDist");


    /* find mean pT in 2-4 GeV range for hadrons */
    Float_t avg24 = MeanPT(hadronPT, 2.0, 4.0);
    Float_t avgtrig24 = MeanPT(hadronTrigPT, 2.0, 4.0);
    Float_t avgtrigphi24 = MeanPT(hadronTrigPhiPT, 2.0, 4.0);
    printf("Avg in 2-4 range, minbias: %f, trigger:%f, trigger+phi:%f\n", avg24, avgtrig24, avgtrigphi24);
    printf("%% increase, trigger: %2.2f%%, trigger+phi: %2.2f%%\n", (avgtrig24/avg24-1.0)*100.0, (avgtrigphi24/avg24 - 1.0)*100.0);

    TH1D* hPTnorm = Normalize(hadronPT, "hPTnorm", 2.0, 4.0);
    hPTnorm->SetLineColor(kBlack);
/*
    TH1D* hPTnorm = (TH1D*)hadronPT->Clone("hPTnorm");
    hPTnorm->GetXaxis()->SetRangeUser(2.01, 3.99);
    hPTnorm->Scale(1.0/hPTnorm->Integral());
    hPTnorm->SetLineColor(kBlack);
    hPTnorm->SetLineWidth(2);
*/
    TH1D* hTrigPTnorm = Normalize(hadronTrigPT, "hTrigPTnorm", 2.0, 4.0);
    hTrigPTnorm->SetLineColor(kBlue+1);
   
   /* 
    TH1D* hTrigPTnorm = (TH1D*)hadronTrigPT->Clone("hTrigPTnorm");
    hTrigPTnorm->GetXaxis()->SetRangeUser(2.01, 3.99);
    hTrigPTnorm->Scale(1.0/hTrigPTnorm->Integral());
    hTrigPTnorm->SetLineColor(kBlue+1);
    hTrigPTnorm->SetLineWidth(2);
   */

    TH1D* hTrigPhiPTnorm = Normalize(hadronTrigPhiPT, "hTrigPhiPTnorm", 2.0, 4.0);
    hTrigPhiPTnorm->SetLineColor(kRed+1);
    
   /* TH1D* hTrigPhiPTnorm = (TH1D*)hadronTrigPhiPT->Clone("hTrigPhiPTnorm");
    hTrigPhiPTnorm->GetXaxis()->SetRangeUser(2.01, 3.99);
    hTrigPhiPTnorm->Scale(1.0/hTrigPhiPTnorm->Integral());
    hTrigPhiPTnorm->SetLineColor(kRed+1);
    hTrigPhiPTnorm->SetLineWidth(2);
  */
    TCanvas* c1 = new TCanvas("c1", "c1", 50, 50, 600, 600);
    c1->cd();
    hPTnorm->Draw("H");
    hTrigPTnorm->Draw("H SAME");
    hTrigPhiPTnorm->Draw("H SAME");

    /* find mean pT in 2-4 GeV range for phi mesons */
    TH2D* kkUSPT = ProjectKKDist(kkUSDist, "kkUSPT");
    TH2D* kkLSPT = ProjectKKDist(kkLSDist, "kkLSPT");
    TH2D* kkUSTrigPT = ProjectKKDist(kkUSTrigDist, "kkUSTrigPT");
    TH2D* kkLSTrigPT = ProjectKKDist(kkLSTrigDist, "kkLSTrigPT");

    TH1D* kkUSPeakPT = ProjectMassToPT(kkUSPT, "kkUSPeakPT", 1.014, 1.026);
    TH1D* kkLSPeakPT = ProjectMassToPT(kkLSPT, "kkLSPeakPT", 1.014, 1.026);
    TH1D* kkUSTrigPeakPT = ProjectMassToPT(kkUSTrigPT, "kkUSTrigPeakPT", 1.014, 1.026);
    TH1D* kkLSTrigPeakPT = ProjectMassToPT(kkLSTrigPT, "kkUSTrigPeakPT", 1.014, 1.026);

    TH1D* kkUSLSBPT = ProjectMassToPT(kkUSPT, "kkUSLSBPT", 0.995, 1.010);
    TH1D* kkLSLSBPT = ProjectMassToPT(kkLSPT, "kkLSLSBPT", 0.995, 1.010);
    TH1D* kkUSTrigLSBPT = ProjectMassToPT(kkUSTrigPT, "kkUSTrigLSBPT", 0.995, 1.010);
    TH1D* kkLSTrigLSBPT = ProjectMassToPT(kkLSTrigPT, "kkUSTrigLSBPT", 0.995, 1.010);

    TH1D* kkUSRSBPT = ProjectMassToPT(kkUSPT, "kkUSRSBPT", 1.020, 1.040);
    TH1D* kkLSRSBPT = ProjectMassToPT(kkLSPT, "kkLSRSBPT", 1.020, 1.040);
    TH1D* kkUSTrigRSBPT = ProjectMassToPT(kkUSTrigPT, "kkUSTrigRSBPT", 1.020, 1.040);
    TH1D* kkLSTrigRSBPT = ProjectMassToPT(kkLSTrigPT, "kkUSTrigRSBPT", 1.020, 1.040);

    Float_t avgphiPT = MeanPT(kkUSPeakPT, 2.0, 4.0);
    Float_t avgphiTrigPT = MeanPT(kkUSTrigPeakPT, 2.0, 4.0);

    printf("Avg in 2-4 range for US KK, minbias: %f, trig: %f\n", avgphiPT, avgphiTrigPT);
    printf("%% increase, trigger: %2.2f\n%%", 100.0*(avgphiTrigPT/avgphiPT - 1.0));

    TH1D* kkUSPeakPTnorm = Normalize(kkUSPeakPT, "kkUSPeakPTnorm", 2.0, 4.0);
    TH1D* kkUSTrigPeakPTnorm = Normalize(kkUSTrigPeakPT, "kkUSTrigPeakPTnorm", 2.0, 4.0);
    TH1D* kkLSPeakPTnorm = Normalize(kkLSPeakPT, "kkLSPeakPTnorm", 2.0, 4.0);
    TH1D* kkLSTrigPeakPTnorm = Normalize(kkLSTrigPeakPT, "kkLSTrigPeakPTnorm", 2.0, 4.0);
    TH1D* kkUSRSBPTnorm = Normalize(kkUSRSBPT, "kkUSRSBPTnorm", 2.0, 4.0);
    TH1D* kkUSTrigRSBPTnorm = Normalize(kkUSTrigRSBPT, "kkUSTrigRSBPTnorm", 2.0, 4.0);
    TH1D* kkLSRSBPTnorm = Normalize(kkLSRSBPT, "kkLSRSBPTnorm", 2.0, 4.0);
    TH1D* kkLSTrigRSBPTnorm = Normalize(kkLSTrigRSBPT, "kkLSTrigRSBPTnorm", 2.0, 4.0);

    TCanvas* c2 = new TCanvas("c2", "c2", 50, 50, 600, 600);
    c2->cd();
    kkUSPeakPT->SetLineColor(kBlack);
    kkUSPeakPT->Draw("H");
    kkUSTrigPeakPT->SetLineColor(kRed+1);
    kkUSTrigPeakPT->Draw("H SAME");

    TCanvas* c3 = new TCanvas("c3", "c3", 50, 50, 600, 600);
    c3->cd();
    kkUSPeakPTnorm->Draw("H");
    kkUSRSBPTnorm->SetLineColor(kGreen+2);
    kkUSRSBPTnorm->Draw("H SAME");
    kkLSPeakPTnorm->SetLineColor(kAzure);
    kkLSPeakPTnorm->Draw("H SAME");

}
