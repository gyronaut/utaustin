void massTPlots(){

    TLegend* legend = new TLegend(0.4614,0.1658,0.8826,0.3944);

    TFile *file1 = new TFile("20170524_Kstar0bar_recon_masswidth_pf70_scaled_error05.root");
    TH1D* mass1 = file1->Get("kstar0mass");
    mass1->SetName("mwScaled");
    mass1->SetTitle("Fit Mass Peak for #bar{K}*^{0}");
    mass1->SetMarkerStyle(20);
    mass1->SetMarkerSize(1.5);
    mass1->SetMarkerColor(1);
    TH1D* width1 = file1->Get("kstar0width");
    width1->SetName("mwWidthScaled");
    width1->SetTitle("Fit Width for #bar{K}*^{0}");
    width1->SetMarkerStyle(20);
    width1->SetMarkerSize(1.5);
    width1->SetMarkerColor(1);
    TH1D* single1 = file1->Get("ptbin21particle5");
    single1->SetLineColor(1);
    single1->SetLineWidth(2);
    TH1D* othersingle1 = file1->Get("ptbin05particle5");
    othersingle1->SetLineColor(1);
    othersingle1->SetLineWidth(2);
    TF1* fit1 = single1->GetFunction("fitPTbin2100particle5");
    TF1* otherfit1 = othersingle1->GetFunction("fitPTbin500particle5");
    fit1->SetBit(TF1::kNotDraw);
    fit1->SetLineColor(kGray+2);
    fit1->SetLineWidth(2);
    otherfit1->SetBit(TF1::kNotDraw);
    otherfit1->SetLineColor(kGray+2);
    otherfit1->SetLineWidth(2);

    TFile *file2 = new TFile("20170524_Kstar0bar_recon_masswidth_pf80_scaled_error05.root");
    TH1D* mass2 = file2->Get("kstar0mass");
    mass2->SetName("mwError05");
    mass2->SetTitle("Fit Mass Peak for K*^{-}");
    mass2->SetMarkerStyle(33);
    mass2->SetMarkerSize(1.5);
    mass2->SetLineColor(kGreen-2);
    mass2->SetMarkerColor(kGreen-2);
    TH1D* width2 = file2->Get("kstar0width");
    width2->SetName("mwWidth05");
    width2->SetTitle("Fit Width for K*^{0}");
    width2->SetMarkerStyle(33);
    width2->SetMarkerSize(1.5);
    width2->SetMarkerColor(kGreen-2);
    width2->SetLineColor(kGreen-2);
    TF1* fit2 = file2->Get("fitPTbin2100particle5");
    TF1* otherfit2 = file2->Get("fitPTbin500particle5");
    fit2->SetLineColor(kGreen-2);
    fit2->SetLineStyle(5);
    fit2->SetLineWidth(4);
    otherfit2->SetLineColor(kGreen-2);
    otherfit2->SetLineWidth(4);
    otherfit2->SetLineStyle(5);



    TFile *file3 = new TFile("20170524_Kstar0bar_recon_masswidth_pf90_scaled_error05.root");
    TH1D* mass3 = file3->Get("kstar0mass");
    mass3->SetName("swError05");
    mass3->SetTitle("simple width, error: 5%");
    mass3->SetMarkerStyle(21);
    mass3->SetMarkerSize(1.5);
    mass3->SetMarkerColor(4);
    mass3->SetLineColor(4);
    TH1D* width3 = file3->Get("kstar0width");
    width3->SetName("swWidth05");
    width3->SetTitle("Fit Width for (K*^{0} + K*^{0})");
    width3->SetMarkerStyle(21);
    width3->SetMarkerSize(1.5);
    width3->SetMarkerColor(4);
    width3->SetLineColor(4);
    TF1* fit3 = file3->Get("fitPTbin2100particle5");
    fit3->SetLineColor(4);
    fit3->SetLineStyle(7);
    fit3->SetLineWidth(4);
    TF1* otherfit3 = file3->Get("fitPTbin500particle5");
    otherfit3->SetLineColor(4);
    otherfit3->SetLineStyle(7);
    otherfit3->SetLineWidth(4);

/*
    TFile *file4 = new TFile("20170522_KKbar_reconsimplewidth_pf160_error01.root");
    TH1D* mass4 = file4->Get("kstar0mass");
    mass4->SetName("swError10");
    mass4->SetTitle("simple width, error: 10%");
    mass4->SetMarkerStyle(25);
    mass4->SetMarkerSize(1.5);
    mass4->SetMarkerColor(3);
*/
    TFile *file5 = new TFile("20170524_Kstar0bar_recon_masswidth_pf100_scaled_error05.root");
    TH1D* mass5 = file5->Get("kstar0mass");
    mass5->SetName("fwError05");
    mass5->SetTitle("fixed width, error: 5%");
    mass5->SetMarkerStyle(22);
    mass5->SetMarkerSize(1.5);
    mass5->SetMarkerColor(2);
    mass5->SetLineColor(2);
    TH1D* width5 = file5->Get("kstar0width");
    width5->SetName("fwWidth05");
    width5->SetTitle("Fit Width for K*^0");
    width5->SetMarkerStyle(22);
    width5->SetMarkerSize(1.5);
    width5->SetMarkerColor(2);
    width5->SetLineColor(2);
    TF1* fit5 = file5->Get("fitPTbin2100particle5");
    fit5->SetLineColor(2);
    fit5->SetLineStyle(3);
    fit5->SetLineWidth(4);
    TF1* otherfit5 = file5->Get("fitPTbin500particle5");
    otherfit5->SetLineColor(2);
    otherfit5->SetLineStyle(3);
    otherfit5->SetLineWidth(4);
    /*
    TFile *file6 = new TFile("20170522_KKbar_reconfixedwidth_pf160_error10.root");
    TH1D* mass6 = file6->Get("kstar0mass");
    mass6->SetName("fwError05");
    mass6->SetTitle("fixed width, error: 10%");
    mass6->SetMarkerStyle(26);
    mass6->SetMarkerSize(1.5);
    mass6->SetMarkerColor(3);

    TFile *file7 = new TFile("20170522_KKbar_reconfixedwidth70_pf160_scaled.root");
    TH1D* mass7 = file7->Get("kstar0mass");
    mass7->SetName("fw70Error05");
    mass7->SetMarkerStyle(22);
    mass7->SetMarkerSize(1.5);
    mass7->SetMarkerColor(kRed+3);
    TH1D* width7 = file7->Get("kstar0collWidth");
    width7->SetName("fw70Width05");
    width7->SetTitle("Fit Width for K*^{0}");
    width7->SetMarkerStyle(22);
    width7->SetMarkerSize(1.5);
    width7->SetMarkerColor(kRed+3);
    width7->SetLineColor(2);
    TF1* fit7 = file7->Get("fitPTbin2100particle5");
    fit7->SetLineColor(kRed+3);
    fit7->SetLineStyle(3);
    fit7->SetLineWidth(3);
*/
   


    TCanvas *cMass = new TCanvas("cMass", "cMass", 50, 50, 600, 600);

    TF1 *pdg = new TF1("pdg", "[0]", 0.0, 4.0);
    pdg->SetParameter(0, 0.8958);
    pdg->SetLineStyle(7);
    pdg->SetLineColor(1);
    pdg->SetLineWidth(2);

    cMass->cd();
    legend->AddEntry(mass1, "Mass Dep. Width, T = 70 MeV", "lpe");
    legend->AddEntry(mass2, "Mass Dep. Width, T = 80 MeV", "lpe");
    legend->AddEntry(mass3, "Mass Dep. Width, T = 90 MeV", "lpe");
    legend->AddEntry(mass5, "Mass Dep. Width, T = 100 MeV", "lpe");
    //legend->AddEntry(mass7, "Fixed Width (70 MeV/c^{2}), 5% Error", "lpe");

    mass1->GetYaxis()->SetRangeUser(0.885, 0.90);
    mass1->GetYaxis()->SetLabelSize(0.03);
    mass1->GetYaxis()->SetTitleOffset(1.4);
    mass1->Draw("P E1");
    mass3->Draw("SAME P E1");
    //mass4->Draw("SAME");
    mass5->Draw("SAME P E1");
    mass2->Draw("SAME P E1");
    pdg->Draw("SAME");
    //mass6->Draw("SAME");
    //mass7->Draw("SAME P E1");

    legend->Draw();

    TLegend *widthLegend = new TLegend(0.4614,0.1658,0.8826,0.3944);
    widthLegend->AddEntry(width1, "Mass Dep. Width, T = 70 MeV", "lpe");
    widthLegend->AddEntry(width2, "Mass Dep. Width, T = 80 MeV", "lpe");
    widthLegend->AddEntry(width3, "Mass Dep. Width, T = 90 MeV", "lpe");
    widthLegend->AddEntry(width5, "Mass Dep. Width, T = 100 MeV", "lpe");
    //widthLegend->AddEntry(width7, "Fixed Width (70 MeV/c^{2}), 5% Error", "lpw");
    
    TCanvas *cWidth = new TCanvas("cWidth", "cWidth", 60, 60, 600, 600);
    cWidth->cd();
    width1->GetYaxis()->SetLabelSize(0.03);
    width1->GetYaxis()->SetTitleOffset(1.4);
    width1->GetYaxis()->SetRangeUser(0.03, 0.12);
    width1->Draw("P E1");
    width3->Draw("SAME P E1");
    width1->Draw("SAME P E1");
    width5->Draw("SAME P E1");
    width2->Draw("SAME P E1");
   //width7->Draw("SAME P E1");
    widthLegend->Draw();

    TLegend *singleLegend = new TLegend(0.4614, 0.1658, 0.8826, 0.3944);
    singleLegend->AddEntry(fit1, "Mass Dep. Width, T = 70 MeV", "l");
    singleLegend->AddEntry(fit2, "Mass Dep. Width, T = 80 MeV", "l");
    singleLegend->AddEntry(fit3, "Mass Dep. Width, T = 90 MeV", "l");
    singleLegend->AddEntry(fit5, "Mass Dep. Width, T = 100 MeV", "l");
    //singleLegend->AddEntry(fit7, "Fixed Width (70 MeV/c^{2})", "l");
    TCanvas *cSingle = new TCanvas("single", "single", 70, 70, 600, 600);
    cSingle->cd();
    single1->SetStats(kFALSE);
    single1->GetYaxis()->SetTitleOffset(1.40);
    single1->GetYaxis()->SetLabelSize(0.03);
    single1->GetXaxis()->SetRangeUser(0.71, 1.03);
    single1->Draw("H");
    fit1->Draw("SAME");
    fit2->Draw("SAME");
    fit3->Draw("SAME");
    fit5->Draw("SAME");
    single1->Draw("H SAME");
    //fit7->Draw("SAME");
    singleLegend->Draw();

    TLegend *othersingleLegend = new TLegend(0.4614, 0.1658, 0.8826, 0.3944);
    othersingleLegend->AddEntry(otherfit1, "Mass Dep. Width, T = 70 MeV", "l");
    othersingleLegend->AddEntry(otherfit2, "Mass Dep. Width, T = 80 MeV", "l");
    othersingleLegend->AddEntry(otherfit3, "Mass Dep. Width, T = 90 MeV", "l");
    othersingleLegend->AddEntry(otherfit5, "Mass Dep. Width, T = 100 MeV", "l");
    //singleLegend->AddEntry(fit7, "Fixed Width (70 MeV/c^{2})", "l");
    TCanvas *cotherSingle = new TCanvas("othersingle", "othersingle", 70, 70, 600, 600);
    cotherSingle->cd();
    othersingle1->SetStats(kFALSE);
    othersingle1->GetYaxis()->SetTitleOffset(1.40);
    othersingle1->GetYaxis()->SetLabelSize(0.03);
    othersingle1->GetXaxis()->SetRangeUser(0.71, 1.03);
    othersingle1->Draw("H");
    otherfit1->Draw("SAME");
    otherfit2->Draw("SAME");
    otherfit3->Draw("SAME");
    otherfit5->Draw("SAME");
    othersingle1->Draw("H SAME");
    //fit7->Draw("SAME");
    othersingleLegend->Draw();   
}
