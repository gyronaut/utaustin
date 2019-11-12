void checkWeights(){
    TFile* nwfile = new TFile("/Users/jtblair/Downloads/AnalysisResults_noweight.root");
    TList* nwlist = (TList*)nwfile->Get("truePhiCorr_mult_0_100");
    THnSparseF* nwhist = (THnSparseF*)nwlist->FindObject("fDphiHPhiz5");
    nwhist->GetAxis(0)->SetRangeUser(4.0, 8.0);
    nwhist->GetAxis(1)->SetRangeUser(2.0, 4.0);
    TH1D* nwhist1D = (TH1D*)nwhist->Projection(2);
    nwhist1D->SetName("noweights");

    TFile* wfile = new TFile("/Users/jtblair/Downloads/AnalysisResults_weights_root6.root");
    TList* wlist = (TList*)wfile->Get("truePhiCorr_mult_0_100");
    THnSparseF* whist = (THnSparseF*)wlist->FindObject("fDphiHPhiz5");
    whist->GetAxis(0)->SetRangeUser(4.0, 8.0);
    whist->GetAxis(1)->SetRangeUser(2.0, 4.0);
    TH1D* whist1D = (TH1D*)whist->Projection(2);
    whist1D->SetName("weights");
    whist1D->SetLineColor(kGreen+2);

    TH1D* ratio = (TH1D*)whist1D->Clone("ratio");
    ratio->Divide(nwhist1D);

    TCanvas* c1 = new TCanvas("c1", "c1", 50, 50, 600, 600);
    c1->cd();
    ratio->Draw();

    TCanvas *c2 = new TCanvas("c2", "c2", 50, 50, 600, 600);
    c2->cd();
    nwhist1D->Draw();
    whist1D->Draw("SAME");

}
