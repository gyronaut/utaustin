void massPlots(){

    TFile *file1 = new TFile("20170721_Kstar0_masswidth_pf160_scaled.root");
    TH1D* mass1 = (TH1D*)file1->Get("kstar0mass");
    mass1->SetName("mwScaled");
    mass1->SetTitle("Fit Mass Peak for K*^{0}");
    mass1->SetMarkerStyle(20);
    mass1->SetMarkerSize(1.5);
    mass1->SetMarkerColor(1);
    TH1D* width1 = (TH1D*)file1->Get("kstar0width");
    width1->SetName("mwWidthScaled");
    width1->SetTitle("Fit Width for K*^{0}");
    width1->SetMarkerStyle(20);
    width1->SetMarkerSize(1.5);
    width1->SetMarkerColor(1);
    TH1D* single1 = (TH1D*)file1->Get("ptbin21particle4");
    single1->SetLineColor(1);
    single1->SetMarkerStyle(21);
    single1->SetMarkerSize(0.3);
    single1->SetLineWidth(3);
    TH1D* othersingle1 = (TH1D*)file1->Get("ptbin05particle4");
    othersingle1->SetLineColor(1);
    othersingle1->SetMarkerStyle(21);
    othersingle1->SetMarkerSize(0.3);
    othersingle1->SetLineWidth(3);
    TF1* fit1 = (TF1*)single1->GetFunction("fitPTbin2100particle4");
    TF1* otherfit1 = (TF1*)othersingle1->GetFunction("fitPTbin500particle4");
    fit1->SetBit(TF1::kNotDraw);
    fit1->SetLineColor(kGray+3);
    fit1->SetLineWidth(4);
    otherfit1->SetBit(TF1::kNotDraw);
    otherfit1->SetLineColor(kGray+3);
    otherfit1->SetLineWidth(4);

    TFile *file2 = new TFile("20170721_Kstar0_masswidth_pf160_scaled_error05.root");
    TH1D* mass2 = (TH1D*)file2->Get("kstar0mass");
    mass2->SetName("mwError05");
    mass2->SetTitle("Fit Mass Peak for K*^{-}");
    mass2->SetMarkerStyle(33);
    mass2->SetMarkerSize(1.5);
    mass2->SetLineColor(kGreen-2);
    mass2->SetMarkerColor(kGreen-2);
    TH1D* width2 = (TH1D*)file2->Get("kstar0width");
    width2->SetName("mwWidth05");
    width2->SetTitle("Fit Width for K*^{0}");
    width2->SetMarkerStyle(33);
    width2->SetMarkerSize(1.5);
    width2->SetMarkerColor(kGreen-2);
    width2->SetLineColor(kGreen-2);
    TH1D* single2 = (TH1D*)file2->Get("ptbin21particle4");
    single2->SetLineColor(17);
    single2->SetFillColor(17);
    single2->SetLineWidth(3);
    TH1D* othersingle2 = (TH1D*)file2->Get("ptbin05particle4");
    othersingle2->SetLineColor(17);
    othersingle2->SetFillColor(17);
    othersingle2->SetLineWidth(3);
    TF1* fit2 = (TF1*)single2->GetFunction("fitPTbin2100particle4");
    TF1* otherfit2 = (TF1*)othersingle2->GetFunction("fitPTbin500particle4");
    fit2->SetLineColor(kGreen-2);
    fit2->SetLineStyle(5);
    fit2->SetLineWidth(4);
    otherfit2->SetLineColor(kGreen-2);
    otherfit2->SetLineWidth(4);
    otherfit2->SetLineStyle(5);



    TFile *file3 = new TFile("20170721_Kstar0_simplewidth_pf160_scaled_error05.root");
    TH1D* mass3 = (TH1D*)file3->Get("kstar0mass");
    mass3->SetName("swError05");
    mass3->SetTitle("simple width, error: 5%");
    mass3->SetMarkerStyle(21);
    mass3->SetMarkerSize(1.5);
    mass3->SetMarkerColor(4);
    mass3->SetLineColor(4);
    TH1D* width3 = (TH1D*)file3->Get("kstar0collWidth");
    width3->SetName("swWidth05");
    width3->SetTitle("Fit Width for (K*^{0} + K*^{0})");
    width3->SetMarkerStyle(21);
    width3->SetMarkerSize(1.5);
    width3->SetMarkerColor(4);
    width3->SetLineColor(4);
    TF1* fit3 = (TF1*)file3->Get("fitPTbin2100particle4");
    fit3->SetLineColor(4);
    fit3->SetLineStyle(7);
    fit3->SetLineWidth(4);
    TF1* otherfit3 = (TF1*)file3->Get("fitPTbin500particle4");
    otherfit3->SetLineColor(4);
    otherfit3->SetLineStyle(7);
    otherfit3->SetLineWidth(4);

/*
    TFile *file4 = new TFile("20170522_Kstar0_reconsimplewidth_pf160_error01.root");
    TH1D* mass4 = file4->Get("kstar0mass");
    mass4->SetName("swError10");
    mass4->SetTitle("simple width, error: 10%");
    mass4->SetMarkerStyle(25);
    mass4->SetMarkerSize(1.5);
    mass4->SetMarkerColor(3);
*/
    TFile *file5 = new TFile("20170721_Kstar0_fixedwidth42_pf160_scaled_error05.root");
    TH1D* mass5 = (TH1D*)file5->Get("kstar0mass");
    mass5->SetName("fwError05");
    mass5->SetTitle("fixed width, error: 5%");
    mass5->SetMarkerStyle(22);
    mass5->SetMarkerSize(1.5);
    mass5->SetMarkerColor(2);
    mass5->SetLineColor(2);
    TH1D* width5 = (TH1D*)file5->Get("kstar0collWidth");
    width5->SetName("fwWidth05");
    width5->SetTitle("Fit Width for K*^0");
    width5->SetMarkerStyle(22);
    width5->SetMarkerSize(1.5);
    width5->SetMarkerColor(2);
    width5->SetLineColor(2);
    TF1* fit5 = (TF1*)file5->Get("fitPTbin2100particle4");
    fit5->SetLineColor(2);
    fit5->SetLineStyle(3);
    fit5->SetLineWidth(4);
    TF1* otherfit5 = (TF1*)file5->Get("fitPTbin500particle4");
    otherfit5->SetLineColor(2);
    otherfit5->SetLineStyle(3);
    otherfit5->SetLineWidth(4);
    /*
    TFile *file6 = new TFile("20170522_Kstar0_reconfixedwidth_pf160_error10.root");
    TH1D* mass6 = file6->Get("kstar0mass");
    mass6->SetName("fwError05");
    mass6->SetTitle("fixed width, error: 10%");
    mass6->SetMarkerStyle(26);
    mass6->SetMarkerSize(1.5);
    mass6->SetMarkerColor(3);

    TFile *file7 = new TFile("20170522_Kstar0_reconfixedwidth70_pf160_scaled.root");
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
    TF1* fit7 = file7->Get("fitPTbin2100particle4");
    fit7->SetLineColor(kRed+3);
    fit7->SetLineStyle(3);
    fit7->SetLineWidth(3);
*/
   


    TCanvas *cMass = new TCanvas("cMass", "cMass", 50, 50, 600, 600);

    TF1 *pdg = new TF1("pdg", "[0]", 0.0, 4.0);
    pdg->SetParameter(0, 0.892);
    pdg->SetLineStyle(7);
    pdg->SetLineColor(1);
    pdg->SetLineWidth(4);

    TPaveText *pdgtext = new TPaveText(0.5520, 0.5637, 0.8322, 0.6946, "NDC");
    pdgtext->AddText("PHSD vacuum mass");
    pdgtext->SetTextSizePixels(32);
    pdgtext->SetBorderSize(0);
    pdgtext->SetTextFont(42);
    pdgtext->SetFillStyle(0);

    cMass->cd();
    TLegend* legend = new TLegend(0.1611,0.6370,0.4832,0.8935);
    legend->AddEntry(mass1, "Mass Dep. Width", "lpe");
    legend->AddEntry(mass2, "Mass Dep. Width (+5% Error)", "lpe");
    legend->AddEntry(mass3, "Simple Width (+5% Error)", "lpe");
    legend->AddEntry(mass5, "Fixed Vacuum Width (+5% Error)", "lpe");
    legend->SetMargin(0.1);
    legend->SetTextSizePixels(28);
    //legend->AddEntry(mass7, "Fixed Width (70 MeV/c^{2}), 5% Error", "lpe");

    mass1->GetYaxis()->SetRangeUser(0.87, 0.91);
    mass1->GetYaxis()->SetLabelSize(0.06);
    mass1->GetYaxis()->SetTitleOffset(1.3);
    mass1->GetYaxis()->SetTitleSize(0.07);
    mass1->GetYaxis()->SetTitleFont(42);
    mass1->GetYaxis()->SetLabelFont(42);
    mass1->GetXaxis()->SetLabelSize(0.06);
    mass1->GetXaxis()->SetTitleSize(0.07);
    mass1->GetXaxis()->SetLabelFont(42);
    mass1->GetXaxis()->SetTitleFont(42);
    mass1->SetTitle("");
    mass1->Draw("P E1");
    mass3->Draw("SAME P E1");
    //mass4->Draw("SAME");
    mass5->Draw("SAME P E1");
    mass2->Draw("SAME P E1");
    pdg->Draw("SAME");
    pdgtext->Draw("SAME");
    //mass6->Draw("SAME");
    //mass7->Draw("SAME P E1");
    TPaveText *masstext = new TPaveText(0.5520, 0.7208, 0.8389, 0.8517, "NDC");
    masstext->AddText("K*^{0} Mass");
    masstext->AddText("At Decay Point");
    //masstext->AddText("Reconstructed");
    masstext->SetBorderSize(0);
    masstext->SetFillStyle(0);
    masstext->SetTextSizePixels(36);
    masstext->SetTextFont(42);
    masstext->Draw();

    legend->Draw();

    TLegend *widthLegend = new TLegend(0.1913,0.6108,0.5134,0.8656);
    widthLegend->AddEntry(width1, "Mass Dep. Width", "lpe");
    widthLegend->AddEntry(width2, "Mass Dep. Width (+5% Error)", "lpe");
    widthLegend->AddEntry(width3, "Simple Width (+5% Error)", "lpe");
    widthLegend->AddEntry(width5, "Fixed Vacuum Width (+5% Error)", "lpe");
    widthLegend->SetTextSizePixels(20);
    //widthLegend->AddEntry(width7, "Fixed Width (70 MeV/c^{2}), 5% Error", "lpw");
    
    TCanvas *cWidth = new TCanvas("cWidth", "cWidth", 60, 60, 600, 600);
    cWidth->cd();
    width1->GetYaxis()->SetLabelSize(0.06);
    width1->GetYaxis()->SetTitleOffset(1.4);
    width1->GetYaxis()->SetRangeUser(0.02, 0.12);
    width1->GetYaxis()->SetTitleSize(0.07);
    width1->GetYaxis()->SetTitleFont(42);
    width1->GetYaxis()->SetLabelFont(42);
    width1->GetXaxis()->SetTitleSize(0.07);
    width1->GetXaxis()->SetLabelSize(0.06);
    width1->GetXaxis()->SetTitleFont(42);
    width1->GetXaxis()->SetLabelFont(42);
    width1->GetYaxis()->SetTitle("Width (GeV/c^{2})");
    width1->SetTitle("");
    width1->Draw("P E1");
    width3->Draw("SAME P E1");
    width1->Draw("SAME P E1");
    width5->Draw("SAME P E1");
    width2->Draw("SAME P E1");
   //width7->Draw("SAME P E1");
    widthLegend->Draw();

    TPaveText *widthtext = new TPaveText(0.5822, 0.7024, 0.8675, 0.8333, "NDC");
    widthtext->AddText("K*^{0} Width");
    widthtext->AddText("At Decay Point");
    //widthtext->AddText("Reconstructed");
    widthtext->SetBorderSize(0);
    widthtext->SetFillStyle(0);
    widthtext->SetTextSizePixels(36);
    widthtext->SetTextFont(42);
    widthtext->Draw();

    TExec *exec1 = new TExec("exec1", "gStyle->SetErrorX(0)");
    TExec *exec2 = new TExec("exec2", "gStyle->SetErrorX(0.5)");

    TLegend *singleLegend = new TLegend(0.1376, 0.5585, 0.5201, 0.8778);
    singleLegend->AddEntry(fit1, "Mass Dep. Width", "l");
    singleLegend->AddEntry(fit2, "#splitline{Mass Dep. Width}{+5% Error}", "l");
    singleLegend->AddEntry(fit3, "#splitline{Simple Width}{+5% Error}", "l");
    singleLegend->AddEntry(fit5, "#splitline{Fixed #Gamma = 50 MeV/c^{2}}{+5% Error}", "l");
    singleLegend->SetTextSizePixels(20);
    //singleLegend->AddEntry(fit7, "Fixed Width (70 MeV/c^{2})", "l");
    TCanvas *cSingle = new TCanvas("single", "single", 70, 70, 600, 600);
    cSingle->cd();
    single1->SetStats(kFALSE);
    single1->SetTitle("");
    single1->GetYaxis()->SetTitleOffset(1.50);
    single1->GetYaxis()->SetLabelSize(0.03);
    single1->GetYaxis()->SetTitleSize(0.04);
    single1->GetYaxis()->SetTitle("Counts / 8 MeV/c^{2}");
    single1->GetXaxis()->SetRangeUser(0.61, 1.09);
    single1->GetXaxis()->SetLabelSize(0.03);
    single1->GetXaxis()->SetTitleSize(0.04);
    single1->GetXaxis()->SetTitle("K^{+}#pi^{-} invariant mass (GeV/c^{2})");
    single1->Draw("E");
    single2->SetStats(kFALSE);
    single2->GetXaxis()->SetRangeUser(0.61, 1.09);
    single2->Draw("SAME E2");
    fit1->Draw("SAME");
    fit2->Draw("SAME");
    fit3->Draw("SAME");
    fit5->Draw("SAME");
    exec1->Draw();
    single1->Draw("E SAME");
    exec2->Draw();
    //fit7->Draw("SAME");
    TLegend *singleLegend2 = new TLegend(0.1695, 0.4738, 0.4211, 0.5279);
    singleLegend2->AddEntry(single2, "#splitline{Added Error}{(5% of peak bin)}", "f");
    singleLegend2->SetFillStyle(0);
    singleLegend2->SetBorderSize(0);
    singleLegend2->SetTextSizePixels(20);
 
    TPaveText *text = new TPaveText(0.5822, 0.7024, 0.8675, 0.8333, "NDC");
    text->AddText("Reconstructed K*^{0}");
    text->AddText("2.0 < p_{T} < 2.2 GeV/c");
    text->SetBorderSize(0);
    text->SetFillStyle(0);
    text->GetLine(1)->SetTextSizePixels(20);
    text->GetLine(0)->SetTextSizePixels(22);
  
    singleLegend->Draw();
    singleLegend2->Draw();
    text->Draw();

    TLegend *othersingleLegend = new TLegend(0.1376, 0.5585, 0.5201, 0.8778);
    othersingleLegend->AddEntry(otherfit1, "Mass Dep. Width", "l");
    othersingleLegend->AddEntry(otherfit2, "#splitline{Mass Dep. Width}{+5% Error}", "l");
    othersingleLegend->AddEntry(otherfit3, "#splitline{Simple Width}{+5% Error}", "l");
    othersingleLegend->AddEntry(otherfit5, "#splitline{Fixed #Gamma = 50 MeV/c^{2}}{+5% Error}", "l");
    //othersingleLegend->SetBorderSize(0);
    //othersingleLegend->SetFillStyle(0);
    othersingleLegend->SetTextSizePixels(20);
    TCanvas *cotherSingle = new TCanvas("othersingle", "othersingle", 70, 70, 600, 600);
    cotherSingle->cd();
    othersingle1->SetStats(kFALSE);
    othersingle1->SetTitle("");
    othersingle1->GetYaxis()->SetTitleOffset(1.50);
    othersingle1->GetYaxis()->SetTitleSize(0.06);
    othersingle1->GetYaxis()->SetLabelSize(0.07);
    othersingle1->GetYaxis()->SetTitleFont(42);
    othersingle1->GetYaxis()->SetLabelFont(42);
    othersingle1->GetYaxis()->SetTitle("Counts / 8 MeV/c^{2}");
    othersingle1->GetXaxis()->SetRangeUser(0.61, 1.09);
    othersingle1->GetXaxis()->SetLabelSize(0.06);
    othersingle1->GetXaxis()->SetTitleSize(0.07);
    othersingle1->GetXaxis()->SetLabelFont(42);
    othersingle1->GetXaxis()->SetTitleFont(42);
    othersingle1->GetXaxis()->SetTitle("K^{+}#pi^{-} invariant mass (GeV/c^{2})");
    othersingle1->Draw("E");
    othersingle2->SetStats(kFALSE);
    othersingle2->GetXaxis()->SetRangeUser(0.61, 1.09);
    othersingle2->Draw("SAME E2");
    otherfit1->Draw("SAME");
    otherfit2->Draw("SAME");
    otherfit3->Draw("SAME");
    otherfit5->Draw("SAME");
    exec1->Draw();
    othersingle1->Draw("E SAME");
    exec2->Draw();
    othersingle1->GetXaxis()->Draw();
    //fit7->Draw("SAME");
    TLegend *othersingleLegend2 = new TLegend(0.1695, 0.4738, 0.4211, 0.5279);
   //othersingleLegend2->AddEntry(othersingle1, "Invariant Mass at Decay Point", "l");
    othersingleLegend2->AddEntry(othersingle2, "#splitline{Added Error}{(5% of peak bin)}", "f");
    othersingleLegend2->SetTextSizePixels(20);
    othersingleLegend2->SetBorderSize(0);
    othersingleLegend2->SetFillStyle(0);
   
    TPaveText *othertext = new TPaveText(0.5882, 0.7024, 0.8675, 0.8333, "NDC");
    othertext->AddText("Reconstructed K*^{0}");
    othertext->AddText("0.4 < p_{T} < 0.6 GeV/c");
    othertext->GetLine(0)->SetTextSizePixels(22);
    othertext->GetLine(1)->SetTextSizePixels(20);
    othertext->SetBorderSize(0);
    othertext->SetFillStyle(0);

    othersingleLegend->Draw();
    othersingleLegend2->Draw();
    othertext->Draw();  

    //testing drawing mass and width on same plot
    TCanvas* cMassWidth = new TCanvas("cmasswidth", "cmasswidth", 50, 50, 1000, 500);
    cMassWidth->SetMargin(0.0, 0.0, 0.0, 0.0);

    cMassWidth->Divide(2, 1, 0.0);
    
    cMassWidth->cd(1)->SetMargin(0.1827, 0.0, 0.1543, 0.0994);
    cMassWidth->cd(1)->SetTicks(0,1);
   
    mass1->Draw("P E1");
    mass3->Draw("SAME P E1");
    //mass4->Draw("SAME");
    mass5->Draw("SAME P E1");
    mass2->Draw("SAME P E1");
    pdg->Draw("SAME");
    pdgtext->Draw("SAME");
    masstext->Draw();
    legend->Draw();

    cMassWidth->cd(2)->SetMargin(0.0, 0.1968, 0.1543, 0.0994);
    cMassWidth->cd(2)->SetTicks(0,1);

    width1->Draw("P E1 Y+");
    width3->Draw("SAME P E1");
    width1->Draw("SAME P E1");
    width5->Draw("SAME P E1");
    width2->Draw("SAME P E1");
    //widthLegend->Draw("SAME");
    widthtext->Draw();


}
