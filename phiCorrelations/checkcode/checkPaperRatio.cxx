void checkPaperRatio(){
    TFile* hfile = new TFile("~/alirepos/utaustin/phiCorrelations/results/paperresults/HEPData-ins1684320-v1-Table_8.root");
    TDirectoryFile* hdir = (TDirectoryFile*)hfile->GetDirectory("Table 8");
    TH1D* hhist = (TH1D*)hdir->Get("Hist1D_y1")->Clone("hhist");
    
    TFile* phifile = new TFile("~/alirepos/utaustin/phiCorrelations/results/paperresults/HEPData-ins1684320-v1-Table_56.root");
    TDirectoryFile* phidir = (TDirectoryFile*)phifile->GetDirectory("Table 56");
    TH1D* phihist = (TH1D*)phidir->Get("Hist1D_y1")->Clone("phihist");
 
    TCanvas* c1 = new TCanvas("c1", "c1", 50, 50, 600, 600);
    c1->cd();
    hhist->Draw("HIST");
    phihist->Draw("SAME HIST");

    Double_t hyield = hhist->Integral(hhist->GetXaxis()->FindBin(2.001), hhist->GetXaxis()->FindBin(3.999), "binwidth");
    Double_t phiyield = phihist->Integral(phihist->GetXaxis()->FindBin(2.001), phihist->GetXaxis()->FindBin(3.999), "binwidth");
    printf("ratio: %f\n", phiyield/hyield);
}
