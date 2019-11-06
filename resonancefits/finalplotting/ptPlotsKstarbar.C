void ptPlotsKstarbar(){
    //TFile *file1 = new TFile("20170721_Kstar0bar_masswidth_recon_pf100_scaled.root");
    TFile *file1 = new TFile("20170527_Kstar0bar_recon_masswidth_pf100_wide_scaled.root");
    TH1D* single1 = (TH1D*)file1->Get("ptbin21particle5");
    single1->SetLineColor(1);
    single1->SetMarkerStyle(21);
    single1->SetMarkerSize(0.5);
    single1->SetLineWidth(3);
    TH1D* othersingle1 = (TH1D*)file1->Get("ptbin05particle5");
    othersingle1->SetLineColor(1);
    othersingle1->SetMarkerStyle(21);
    othersingle1->SetMarkerSize(0.5);
    othersingle1->SetLineWidth(3);
    TF1* fit1 = (TF1*)single1->GetFunction("fitPTbin2100particle5");
    TF1* otherfit1 = (TF1*)othersingle1->GetFunction("fitPTbin500particle5");
    fit1->SetBit(TF1::kNotDraw);
    fit1->SetLineColor(kGray+3);
    fit1->SetLineWidth(3);
    otherfit1->SetBit(TF1::kNotDraw);
    otherfit1->SetLineColor(kGray+3);
    otherfit1->SetLineWidth(3);

    //TFile *file2 = new TFile("20170721_Kstar0bar_masswidth_recon_pf100_scaled_error05.root");
    TFile *file2 = new TFile("20170527_Kstar0bar_recon_masswidth_pf100_wide_scaled_error05.root");
    TH1D* single2 = (TH1D*)file2->Get("ptbin21particle5");
    single2->SetLineColor(17);
    single2->SetFillColor(17);
    single2->SetLineWidth(3);
    TH1D* othersingle2 = (TH1D*)file2->Get("ptbin05particle5");
    othersingle2->SetLineColor(17);
    othersingle2->SetFillColor(17);
    othersingle2->SetLineWidth(3);
    TF1* fit2 = (TF1*)single2->GetFunction("fitPTbin2100particle5");
    TF1* otherfit2 = (TF1*)othersingle2->GetFunction("fitPTbin500particle5");
    fit2->SetLineColor(kGreen-2);
    fit2->SetLineStyle(5);
    fit2->SetLineWidth(5);
    otherfit2->SetLineColor(kGreen-2);
    otherfit2->SetLineWidth(5);
    otherfit2->SetLineStyle(5);

    //TFile *file3 = new TFile("20170721_Kstar0bar_simplewidth_recon_pf100_scaled_error05.root");
    TFile *file3 = new TFile("20170527_Kstar0bar_recon_simplewidth_pf100_wide_scaled_error05.root");
    TF1* fit3 = (TF1*)file3->Get("fitPTbin2100particle5");
    fit3->SetLineColor(4);
    fit3->SetLineStyle(7);
    fit3->SetLineWidth(5);
    TF1* otherfit3 = (TF1*)file3->Get("fitPTbin500particle5");
    otherfit3->SetLineColor(4);
    otherfit3->SetLineStyle(7);
    otherfit3->SetLineWidth(5);

    //TFile *file5 = new TFile("20170721_Kstar0bar_fixedwidth42_recon_pf100_scaled_error05.root");
    TFile *file5 = new TFile("20170527_Kstar0bar_recon_fixedwidth_pf100_wide_scaled_error05.root");
    TF1* fit5 = (TF1*)file5->Get("fitPTbin2100particle5");
    fit5->SetLineColor(2);
    fit5->SetLineStyle(3);
    fit5->SetLineWidth(5);
    TF1* otherfit5 = (TF1*)file5->Get("fitPTbin500particle5");
    otherfit5->SetLineColor(2);
    otherfit5->SetLineStyle(3);
    otherfit5->SetLineWidth(5);

    //KStar0 AT DECAY //
    //TFile *decayfile1 = new TFile("20170721_Kstar0bar_masswidth_pf160_scaled.root");
    TFile *decayfile1 = new TFile("20170527_Kstar0bar_masswidth_pf160_scaled.root");
    TH1D* decaysingle1 = (TH1D*)decayfile1->Get("ptbin21particle5");
    decaysingle1->SetLineColor(1);
    decaysingle1->SetMarkerStyle(21);
    decaysingle1->SetMarkerSize(0.5);
    decaysingle1->SetLineWidth(3);
    TH1D* otherdecaysingle1 = (TH1D*)decayfile1->Get("ptbin05particle5");
    otherdecaysingle1->SetLineColor(1);
    otherdecaysingle1->SetMarkerStyle(21);
    otherdecaysingle1->SetMarkerSize(0.5);
    otherdecaysingle1->SetLineWidth(3);
    TF1* decayfit1 = (TF1*)decaysingle1->GetFunction("fitPTbin2100particle5");
    TF1* otherdecayfit1 = (TF1*)otherdecaysingle1->GetFunction("fitPTbin500particle5");
    decayfit1->SetBit(TF1::kNotDraw);
    decayfit1->SetLineColor(kGray+3);
    decayfit1->SetLineWidth(3);
    otherdecayfit1->SetBit(TF1::kNotDraw);
    otherdecayfit1->SetLineColor(kGray+3);
    otherdecayfit1->SetLineWidth(3);

    //TFile *decayfile2 = new TFile("20170721_Kstar0bar_masswidth_pf160_scaled_error05.root");
    TFile *decayfile2 = new TFile("20170527_Kstar0bar_masswidth_pf160_scaled_error05.root");
    TH1D* decaysingle2 = (TH1D*)decayfile2->Get("ptbin21particle5");
    decaysingle2->SetLineColor(17);
    decaysingle2->SetFillColor(17);
    decaysingle2->SetLineWidth(3);
    TH1D* otherdecaysingle2 = (TH1D*)decayfile2->Get("ptbin05particle5");
    otherdecaysingle2->SetLineColor(17);
    otherdecaysingle2->SetFillColor(17);
    otherdecaysingle2->SetLineWidth(3);
    TF1* decayfit2 = (TF1*)decaysingle2->GetFunction("fitPTbin2100particle5");
    TF1* otherdecayfit2 = (TF1*)otherdecaysingle2->GetFunction("fitPTbin500particle5");
    decayfit2->SetLineColor(kGreen-2);
    decayfit2->SetLineStyle(5);
    decayfit2->SetLineWidth(5);
    otherdecayfit2->SetLineColor(kGreen-2);
    otherdecayfit2->SetLineWidth(4);
    otherdecayfit2->SetLineStyle(5);

    //TFile *decayfile3 = new TFile("20170721_Kstar0bar_simplewidth_pf160_scaled_error05.root");
    TFile *decayfile3 = new TFile("20170527_Kstar0bar_simplewidth_pf160_scaled_error05.root");
    TF1* decayfit3 = (TF1*)decayfile3->Get("fitPTbin2100particle5");
    decayfit3->SetLineColor(4);
    decayfit3->SetLineStyle(7);
    decayfit3->SetLineWidth(5);
    TF1* otherdecayfit3 = (TF1*)decayfile3->Get("fitPTbin500particle5");
    otherdecayfit3->SetLineColor(4);
    otherdecayfit3->SetLineStyle(7);
    otherdecayfit3->SetLineWidth(5);

    //TFile *decayfile5 = new TFile("20170721_Kstar0bar_fixedwidth42_pf160_scaled_error05.root");
    TFile *decayfile5 = new TFile("20170527_Kstar0bar_fixedwidth_pf160_scaled_error05.root");
    TF1* decayfit5 = (TF1*)decayfile5->Get("fitPTbin2100particle5");
    decayfit5->SetLineColor(2);
    decayfit5->SetLineStyle(3);
    decayfit5->SetLineWidth(5);
    TF1* otherdecayfit5 = (TF1*)decayfile5->Get("fitPTbin500particle5");
    otherdecayfit5->SetLineColor(2);
    otherdecayfit5->SetLineStyle(3);
    otherdecayfit5->SetLineWidth(5);


    //Do plots for decay and recon side by side, for high PT////////////////////////////////////////////////
    TExec *exec1 = new TExec("exec1", "gStyle->SetErrorX(0)");
    TExec *exec2 = new TExec("exec2", "gStyle->SetErrorX(0.5)");

    TCanvas *cSingle = new TCanvas("single", "single", 70, 70, 1000, 600);
    cSingle->SetMargin(0.0, 0.0, 0.0, 0.0);
    cSingle->Divide(2,1,0.0);
    cSingle->cd(2)->SetMargin(0.0, 0.1867, 0.1326, 0.0977);
    cSingle->cd(2)->SetTicks(0,1);

    single1->SetStats(kFALSE);
    single1->SetTitle("");
    single1->GetYaxis()->SetTitleOffset(1.50);
    single1->GetYaxis()->SetLabelSize(0.05);
    single1->GetYaxis()->SetTitleSize(0.06);
    single1->GetYaxis()->SetTitleFont(42);
    single1->GetYaxis()->SetLabelFont(42);
    single1->GetYaxis()->SetTitle("Counts / 8 MeV/c^{2}");
    single1->GetXaxis()->SetRangeUser(0.61, 1.09);
    single1->GetXaxis()->SetLabelSize(0.05);
    single1->GetXaxis()->SetTitleSize(0.06);
    single1->GetXaxis()->SetLabelFont(42);
    single1->GetXaxis()->SetTitleFont(42);
    single1->GetXaxis()->SetTitle("K^{-}#pi^{+} invariant mass (GeV/c^{2})");
    single1->Draw("E Y+");
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

    TPaveText *text = new TPaveText(0.3715, 0.7592, 0.6586, 0.8901, "NDC");
    text->AddText("Reconstructed #bar{K}*^{0}");
    text->AddText("2.0 < p_{T} < 2.2 GeV/c");
    text->SetBorderSize(0);
    text->SetFillStyle(0);
    text->GetLine(1)->SetTextSizePixels(28);
    text->GetLine(0)->SetTextSizePixels(32);
  
    text->Draw();

    gPad->RedrawAxis();

    //Do DECAY part
    cSingle->cd(1)->SetMargin(0.1727, 0.0, 0.1326, 0.0977);
    cSingle->cd(1)->SetTicks(0,1);
    TLegend *singleLegend = new TLegend(0.2048, 0.4223, 0.5884, 0.7400);
    singleLegend->AddEntry(fit1, "Mass Dep. Width", "l");
    singleLegend->AddEntry(fit2, "Mass Dep. Width (+5% Error)", "l");
    singleLegend->AddEntry(fit3, "Simple Width (+5% Error)", "l");
    singleLegend->AddEntry(fit5, "Fixed Vacuum Width (+5% Error)", "l");
    singleLegend->SetTextSizePixels(20);

    TLegend *singleLegend2 = new TLegend(0.2430, 0.3351, 0.4940, 0.3892);
    singleLegend2->AddEntry(single2, "#splitline{Added Error}{(5% of peak bin)}", "f");
    singleLegend2->SetFillStyle(0);
    singleLegend2->SetBorderSize(0);
    singleLegend2->SetTextSizePixels(20); 

    decaysingle1->SetStats(kFALSE);
    decaysingle1->SetTitle("");
    decaysingle1->GetYaxis()->SetTitleOffset(1.50);
    decaysingle1->GetYaxis()->SetLabelSize(0.05);
    decaysingle1->GetYaxis()->SetTitleSize(0.06);
    decaysingle1->GetYaxis()->SetTitleFont(42);
    decaysingle1->GetYaxis()->SetLabelFont(42);
    decaysingle1->GetYaxis()->SetTitle("Counts / 8 MeV/c^{2}");
    decaysingle1->GetXaxis()->SetRangeUser(0.61, 1.09);
    decaysingle1->GetXaxis()->SetLabelSize(0.05);
    decaysingle1->GetXaxis()->SetTitleSize(0.06);
    decaysingle1->GetXaxis()->SetTitleFont(42);
    decaysingle1->GetXaxis()->SetLabelFont(42);
    decaysingle1->GetXaxis()->SetTitle("K^{-}#pi^{+} invariant mass (GeV/c^{2})");
    decaysingle1->Draw("E");
    decaysingle2->SetStats(kFALSE);
    decaysingle2->GetXaxis()->SetRangeUser(0.61, 1.09);
    decaysingle2->Draw("SAME E2");
    decayfit1->Draw("SAME");
    decayfit2->Draw("SAME");
    decayfit3->Draw("SAME");
    decayfit5->Draw("SAME");
    exec1->Draw();
    decaysingle1->Draw("E SAME");
    exec2->Draw();
    //fit7->Draw("SAME");

    TPaveText *text2 = new TPaveText(0.5964, 0.7539, 0.8835, 0.8848, "NDC");
    text2->AddText("#bar{K}*^{0} at Decay Point");
    text2->AddText("2.0 < p_{T} < 2.2 GeV/c");
    text2->SetBorderSize(0);
    text2->SetFillStyle(0);
    text2->GetLine(1)->SetTextSizePixels(28);
    text2->GetLine(0)->SetTextSizePixels(32);
  
    singleLegend->Draw();
    singleLegend2->Draw();
    text2->Draw();

    gPad->RedrawAxis();

    //DO LOW MOMENTUM/////////////////////////////////////////////////////////////

    TCanvas *cotherSingle = new TCanvas("othersingle", "othersingle", 70, 70, 1000, 800);
    cotherSingle->SetMargin(0.0, 0.0, 0.0, 0.0);
    cotherSingle->Divide(2,1,0.0);
    cotherSingle->cd(2)->SetMargin(0.0, 0.1567, 0.1326, 0.0577);
    cotherSingle->cd(2)->SetTicks(0,1);

    TGaxis::SetMaxDigits(3);
    othersingle1->SetStats(kFALSE);
    othersingle1->SetTitle("");
    othersingle1->GetYaxis()->SetTitleOffset(1.25);
    othersingle1->GetYaxis()->SetLabelSize(0.05);
    othersingle1->GetYaxis()->SetTitleSize(0.06);
    othersingle1->GetYaxis()->SetTitleFont(42);
    othersingle1->GetYaxis()->SetLabelFont(42);
    othersingle1->GetYaxis()->SetTitle("Counts / (8 MeV/c^{2})");
    othersingle1->GetXaxis()->SetRangeUser(0.61, 1.09);
    othersingle1->GetXaxis()->SetNdivisions(509);
    othersingle1->GetXaxis()->SetLabelSize(0.05);
    othersingle1->GetXaxis()->SetTitleSize(0.06);
    othersingle1->GetXaxis()->SetLabelFont(42);
    othersingle1->GetXaxis()->SetTitleFont(42);
    othersingle1->GetXaxis()->SetTitle("K^{-}#pi^{+} invariant mass (GeV/c^{2})");
    othersingle1->Draw("E Y+");
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
    //fit7->Draw("SAME");

    TPaveText *othertext = new TPaveText(0.4309, 0.7869, 0.7174, 0.9168, "NDC");
    othertext->AddText("Reconstructed #bar{K}*^{0}");
    othertext->AddText("0.4 < p_{T} < 0.6 GeV/c");
    othertext->SetBorderSize(0);
    othertext->SetFillStyle(0);
    othertext->GetLine(1)->SetTextSizePixels(32);
    othertext->GetLine(0)->SetTextSizePixels(36);
    othertext->SetTextFont(42);
  
    othertext->Draw();
 
    gPad->RedrawAxis();

    //Do DECAY part
    cotherSingle->cd(1)->SetMargin(0.1427, 0.0, 0.1326, 0.0577);
    cotherSingle->cd(1)->SetTicks(0,1);

    otherdecaysingle1->SetStats(kFALSE);
    otherdecaysingle1->SetTitle("");
    otherdecaysingle1->GetYaxis()->SetTitleOffset(1.15);
    otherdecaysingle1->GetYaxis()->SetLabelSize(0.05);
    otherdecaysingle1->GetYaxis()->SetTitleSize(0.06);
    otherdecaysingle1->GetYaxis()->SetTitleFont(42);
    otherdecaysingle1->GetYaxis()->SetLabelFont(42);
    otherdecaysingle1->GetYaxis()->SetTitle("Counts / (8 MeV/c^{2})");
    otherdecaysingle1->GetXaxis()->SetRangeUser(0.61, 1.09);
    otherdecaysingle1->SetNdivisions(509);
    otherdecaysingle1->GetXaxis()->SetLabelSize(0.05);
    otherdecaysingle1->GetXaxis()->SetTitleSize(0.06);
    otherdecaysingle1->GetXaxis()->SetTitleFont(42);
    otherdecaysingle1->GetXaxis()->SetLabelFont(42);
    otherdecaysingle1->GetXaxis()->SetTitle("K^{-}#pi^{+} invariant mass (GeV/c^{2})");
    otherdecaysingle1->Draw("E");
    otherdecaysingle2->SetStats(kFALSE);
    otherdecaysingle2->GetXaxis()->SetRangeUser(0.61, 1.09);
    otherdecaysingle2->Draw("SAME E2");
    otherdecayfit1->Draw("SAME");
    otherdecayfit2->Draw("SAME");
    otherdecayfit3->Draw("SAME");
    otherdecayfit5->Draw("SAME");
    exec1->Draw();
    otherdecaysingle1->Draw("E SAME");
    exec2->Draw();
    //fit7->Draw("SAME");

    TPaveText *othertext2 = new TPaveText(0.2445, 0.5460, 0.5311, 0.6774, "NDC");
    othertext2->AddText("#bar{K}*^{0} at Decay Point");
    othertext2->AddText("0.4 < p_{T} < 0.6 GeV/c");
    othertext2->SetBorderSize(0);
    othertext2->SetFillStyle(0);
    othertext2->GetLine(1)->SetTextSizePixels(30);
    othertext2->GetLine(0)->SetTextSizePixels(34);
    othertext2->SetTextFont(42);

    othertext2->Draw();

    TLegend *othersingleLegend = new TLegend(0.1683, 0.6993, 0.9739, 0.9358);
    othersingleLegend->AddEntry(fit1, "Mass Dep. Width", "l");
    othersingleLegend->AddEntry(fit2, "Mass Dep. Width (+5% Error)", "l");
    othersingleLegend->AddEntry(fit3, "Simple Width (+5% Error)", "l");
    othersingleLegend->AddEntry(fit5, "Fixed Vacuum Width (+5% Error)", "l");
    othersingleLegend->SetMargin(0.15);
    othersingleLegend->SetTextSizePixels(34);

    TLegend *othersingleLegend2 = new TLegend(0.1723, 0.4000, 0.5210, 0.4555);
    othersingleLegend2->AddEntry(othersingle2, "#splitline{Added Error}{(5% of peak bin)}", "f");
    othersingleLegend2->SetFillStyle(0);
    othersingleLegend2->SetBorderSize(0);
    othersingleLegend2->SetMargin(0.20);
    othersingleLegend2->SetTextSizePixels(34); 
    othersingleLegend->Draw();
    //othersingleLegend2->Draw();

    cotherSingle->cd(2);
    othersingleLegend2->Draw();

    gPad->RedrawAxis();

}
