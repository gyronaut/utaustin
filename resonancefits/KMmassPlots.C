void KMmassPlots(){

    TLegend* legend = new TLegend(0.4614,0.1658,0.8826,0.3944);

    TFile *file1 = new TFile("KMoutput_SMinvm_masswidth_pf160_error05.root");
    TH1D* mass1 = file1->Get("kstar0mass");
    mass1->SetName("mwError05");
    mass1->SetTitle("Fit Mass Peak for K*^{-}");
    mass1->SetMarkerStyle(20);
    mass1->SetMarkerSize(1.5);
    mass1->SetMarkerColor(1);
    TH1D* width1 = file1->Get("kstar0width");
    width1->SetName("mwWidth05");
    width1->SetTitle("Fit Width for K*^{-}");
    width1->SetMarkerStyle(20);
    width1->SetMarkerSize(1.5);
    width1->SetMarkerColor(1);
    TH1D* single1 = file1->Get("ptbin05particle4");
    TF1* fit1 = file1->Get("fitPTbin500particle4");
    fit1->SetLineColor(1);
    fit1->SetLineWidth(3);

    TFile *file3 = new TFile("KMoutput_SMinvm_simplewidth_pf160_error05.root");
    TH1D* mass3 = file3->Get("kstar0mass");
    mass3->SetName("swError05");
    mass3->SetTitle("simple width, error: 5%");
    mass3->SetMarkerStyle(21);
    mass3->SetMarkerSize(1.5);
    mass3->SetMarkerColor(4);
    mass3->SetLineColor(4);
    TH1D* width3 = file3->Get("kstar0collWidth");
    width3->SetName("swWidth05");
    width3->SetTitle("Fit Width for K*^{-}");
    width3->SetMarkerStyle(21);
    width3->SetMarkerSize(1.5);
    width3->SetMarkerColor(4);
    width3->SetLineColor(4);
    TF1* fit3 = file3->Get("fitPTbin500particle4");
    fit3->SetLineColor(4);
    fit3->SetLineStyle(7);
    fit3->SetLineWidth(3);   


    TCanvas *cMass = new TCanvas("cMass", "cMass", 50, 50, 600, 600);

    TF1 *pdg = new TF1("pdg", "[0]", 0.0, 4.0);
    pdg->SetParameter(0, 0.8912);
    pdg->SetLineStyle(7);
    pdg->SetLineColor(1);
    pdg->SetLineWidth(2);

    cMass->cd();
    legend->AddEntry(mass1, "Mass Width, 5% Error", "lpe");
    legend->AddEntry(mass3, "Simple Width, 5% Error", "lpe");
    //legend->AddEntry(mass5, "Fixed Width (50 MeV/c^{2}), 5% Error", "lpe");
    //legend->AddEntry(mass7, "Fixed Width (70 MeV/c^{2}), 5% Error", "lpe");

    mass1->GetYaxis()->SetRangeUser(0.885, 0.91);
    mass1->GetYaxis()->SetLabelSize(0.03);
    mass1->GetYaxis()->SetTitleOffset(1.4);
    mass1->Draw("P E1");
    //mass2->Draw("SAME");
    mass3->Draw("SAME P E1");
    mass1->Draw("SAME P E1");
    //mass4->Draw("SAME");
    //mass5->Draw("SAME P E1");
    pdg->Draw("SAME");
    //mass6->Draw("SAME");
    //mass7->Draw("SAME P E1");

    legend->Draw();

    TLegend *widthLegend = new TLegend(0.4614,0.1658,0.8826,0.3944);
    widthLegend->AddEntry(width1, "Mass Width, 5% Error", "lpe");
    widthLegend->AddEntry(width3, "Simple Width, 5% Error", "lpe");
    //widthLegend->AddEntry(width5, "Fixed Width (50 MeV/c^{2}), 5% Error", "lpe");
    //widthLegend->AddEntry(width7, "Fixed Width (70 MeV/c^{2}), 5% Error", "lpw");
    
    TCanvas *cWidth = new TCanvas("cWidth", "cWidth", 60, 60, 600, 600);
    cWidth->cd();
    width1->GetYaxis()->SetLabelSize(0.03);
    width1->GetYaxis()->SetTitleOffset(1.4);
    width1->GetYaxis()->SetRangeUser(0.03, 0.12);
    width1->Draw("P E1");
    width3->Draw("SAME P E1");
    width1->Draw("SAME P E1");
    //width5->Draw("SAME P E1");
    //width7->Draw("SAME P E1");
    widthLegend->Draw();

    TLegend *singleLegend = new TLegend(0.4614, 0.1658, 0.8826, 0.3944);
    singleLegend->AddEntry(fit1, "Mass Width Fit", "l");
    singleLegend->AddEntry(fit3, "Simple Width Fit", "l");
    //singleLegend->AddEntry(fit5, "Fixed Width (50 MeV/c^{2})", "l");
    //singleLegend->AddEntry(fit7, "Fixed Width (70 MeV/c^{2})", "l");
    TCanvas *cSingle = new TCanvas("single", "single", 70, 70, 600, 600);
    cSingle->cd();
    single1->SetStats(kFALSE);
    single1->GetYaxis()->SetTitleOffset(1.40);
    single1->GetYaxis()->SetLabelSize(0.03);
    single1->Draw("H");
    fit1->Draw("SAME");
    fit3->Draw("SAME");
    //fit5->Draw("SAME");
    //fit7->Draw("SAME");
    singleLegend->Draw();

}
