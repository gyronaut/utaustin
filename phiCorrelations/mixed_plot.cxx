void mixed_plot(){
    TFile *histofile = new TFile("AnalysisResults_0001.root");
    histofile->cd("PhiReconstruction");
    

    THnSparseF *ffMixedUS = (THnSparseF*) InvMass->FindObject("fDphiHPhiMixed");
    
    ffMixedUS->GetAxis(0)->SetRangeUser(4.0, 8.0);
    ffMixedUS->GetAxis(1)->SetRangeUser(2.0, 4.0);
    TH2D *twocorr = ffMixedUS->Projection(2,3);
    twocorr->SetName("twocorr");

    TCanvas *c1 = new TCanvas("c1", "c1", 50, 50, 600, 600);
    c1->cd();
    twocorr->Draw("Surf2");

    TFile* output = new TFile("mixed.root", "RECREATE");
    twocorr->Write();
}
