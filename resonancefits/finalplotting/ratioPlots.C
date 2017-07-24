void ratioPlots(){
    TH1D* totHist1;
    TH1D* totHist2;
    TH1D* totBarHist1;
    TH1D* totBarHist2;
    TFile* kstarRecon  = new TFile("20170721_Kstar0_masswidth_recon_pf100_scaled.root");
    TFile* kstarDecay = new TFile("20170721_Kstar0_masswidth_pf160_scaled.root");
    TFile* kstarbarRecon = new TFile("20170721_Kstar0bar_masswidth_recon_pf100_scaled.root");
    TFile* kstarbarDecay = new TFile("20170721_Kstar0bar_masswidth_pf160_scaled.root");

    for(int i = 0; i < 20; i++){
        if(i <5){  
            TH1D* hist1 = kstarRecon->Get(Form("ptbin0%iparticle3", (i*2)+1));
            TH1D* hist2 = kstarDecay->Get(Form("ptbin0%iparticle4", (i*2)+1));
            TH1D* barhist1 = kstarbarRecon->Get(Form("ptbin0%iparticle5", (i*2)+1));
            TH1D* barhist2 = kstarbarDecay->Get(Form("ptbin0%iparticle5", (i*2)+1));
        }else{
            TH1D* hist1 = kstarRecon->Get(Form("ptbin%iparticle3", (i*2)+1));
            TH1D* hist2 = kstarDecay->Get(Form("ptbin%iparticle4", (i*2)+1));
            TH1D* barhist1 = kstarbarRecon->Get(Form("ptbin%iparticle5", (i*2)+1));
            TH1D* barhist2 = kstarbarDecay->Get(Form("ptbin%iparticle5", (i*2)+1));
        }
        hist2->SetName("kstardecay");
        hist2->GetXaxis()->SetRangeUser(0.6, 1.2);
        hist2->Sumw2();
        hist1->SetName("kstarrecon");
        hist1->GetXaxis()->SetRangeUser(0.6, 1.2);
        hist1->Sumw2();
        barhist1->SetName("kstarbarrecon");
        barhist1->GetXaxis()->SetRangeUser(0.6, 1.2);
        barhist1->Sumw2();
        barhist2->SetName("kstarbardecay");
        barhist2->GetXaxis()->SetRangeUser(0.6, 1.2);
        barhist2->Sumw2();

        if(i==0){
            hist1->GetFunction("fitPTbin100particle3")->Delete();
            hist2->GetFunction("fitPTbin100particle4")->Delete();
            barhist1->GetFunction("fitPTbin100particle5")->Delete();
            barhist2->GetFunction("fitPTbin100particle5")->Delete();
            totHist1 = (TH1D*)hist1->Clone();
            totHist2 = (TH1D*)hist2->Clone();
            totBarHist1 = (TH1D*)barhist1->Clone();
            totBarHist2 = (TH1D*)barhist2->Clone();
       }else{
            totHist1->Add(hist1);
            totHist2->Add(hist2);
            totBarHist1->Add(barhist1);
            totBarHist2->Add(barhist2);
       }
    }

    //settting the pole mass bin:
    //totBarHist1->SetBinContent(totBarHist1->GetXaxis()->FindBin(0.892), 0.5*(totBarHist1->GetBinContent(totBarHist1->GetXaxis()->FindBin(0.892)-1) + totBarHist1->GetBinContent(totBarHist1->GetXaxis()->FindBin(0.892)+1)));
    //totHist1->SetBinContent(totHist1->GetXaxis()->FindBin(0.892), 0.5*(totHist1->GetBinContent(totHist1->GetXaxis()->FindBin(0.892)-1) + totHist1->GetBinContent(totHist1->GetXaxis()->FindBin(0.892)+1)));


    TH1D* ratio = (TH1D*)totHist1->Clone("ratio");
    ratio->Divide(totHist2);
    ratio->SetTitle("");
    ratio->GetXaxis()->SetTitle("Invariant Mass (GeV/c^{2})");
    ratio->GetXaxis()->SetTitleFont(42);
    ratio->GetXaxis()->SetTitleSize(0.06);
    ratio->GetXaxis()->SetLabelSize(0.05);
    ratio->GetXaxis()->SetTitleOffset(1.05);
    ratio->GetXaxis()->SetRangeUser(0.6, 1.2);
    ratio->GetYaxis()->SetTitle("#frac{Reconstructed}{At Decay Point}");
    ratio->GetYaxis()->SetTitleFont(42);
    ratio->GetYaxis()->SetTitleSize(0.06);
    ratio->GetYaxis()->SetLabelSize(0.05);
    ratio->GetYaxis()->SetTitleOffset(1.2);
    ratio->GetYaxis()->SetRangeUser(0.05, 100.0);
    ratio->SetMarkerColor(4);
    ratio->SetMarkerSize(1.6);
    ratio->SetMarkerStyle(22);
    ratio->SetLineColor(4);
    ratio->SetLineWidth(3);

    ratio->SetStats(kFALSE);

    TH1D* testratio = (TH1D*)totHist1->Clone("testratio");
    testratio->Divide(totBarHist2);
    testratio->SetTitle("");
    testratio->GetXaxis()->SetTitle("Invariant Mass (GeV/c^{2})");
    testratio->GetXaxis()->SetTitleFont(42);
    testratio->GetXaxis()->SetTitleSize(0.06);
    testratio->GetXaxis()->SetLabelSize(0.05);
    testratio->GetXaxis()->SetTitleOffset(1.05);
    testratio->GetXaxis()->SetRangeUser(0.6, 1.2);
    testratio->GetYaxis()->SetTitle("#frac{Reconstructed}{At Decay Point}");
    testratio->GetYaxis()->SetTitleFont(42);
    testratio->GetYaxis()->SetTitleSize(0.06);
    testratio->GetYaxis()->SetLabelSize(0.05);
    testratio->GetYaxis()->SetTitleOffset(1.2);
    testratio->GetYaxis()->SetRangeUser(0.15, 11);
    testratio->SetMarkerColor(2);
    testratio->SetMarkerSize(1.6);
    testratio->SetMarkerStyle(22);
    testratio->SetStats(kFALSE);

    TH1D* barratio = (TH1D*)totBarHist1->Clone("barratio");
    barratio->Divide(totBarHist2);
    barratio->SetTitle("");
    barratio->GetXaxis()->SetTitle("Invariant Mass (GeV/c^{2})");
    barratio->GetXaxis()->SetTitleFont(42);
    barratio->GetXaxis()->SetTitleSize(0.06);
    barratio->GetXaxis()->SetLabelSize(0.05);
    barratio->GetXaxis()->SetTitleOffset(1.05);
    barratio->GetXaxis()->SetRangeUser(0.6, 1.2);
    barratio->GetYaxis()->SetTitle("#frac{Reconstructed}{At Decay Point}");
    barratio->GetYaxis()->SetTitleFont(42);
    barratio->GetYaxis()->SetTitleSize(0.06);
    barratio->GetYaxis()->SetLabelSize(0.05);
    barratio->GetYaxis()->SetTitleOffset(1.2);
    barratio->SetMarkerColor(kRed-4);
    barratio->SetMarkerSize(1.6);
    barratio->SetMarkerStyle(34);
    barratio->SetLineColor(kRed-4);
    barratio->SetLineWidth(2.0);
    barratio->SetStats(kFALSE);

    TH1D* testbarratio = (TH1D*)totBarHist1->Clone("testbarratio");
    testbarratio->Divide(totHist2);
    testbarratio->SetTitle("");
    testbarratio->GetXaxis()->SetTitle("Invariant Mass (GeV/c^{2})");
    testbarratio->GetXaxis()->SetTitleFont(42);
    testbarratio->GetXaxis()->SetTitleSize(0.06);
    testbarratio->GetXaxis()->SetLabelSize(0.05);
    testbarratio->GetXaxis()->SetTitleOffset(1.05);
    testbarratio->GetXaxis()->SetRangeUser(0.6, 1.2);
    testbarratio->GetYaxis()->SetTitle("#frac{Reconstructed}{At Decay Point}");
    testbarratio->GetYaxis()->SetTitleFont(42);
    testbarratio->GetYaxis()->SetTitleSize(0.06);
    testbarratio->GetYaxis()->SetLabelSize(0.05);
    testbarratio->GetYaxis()->SetTitleOffset(1.2);
    testbarratio->SetMarkerColor(8);
    testbarratio->SetMarkerSize(1.6);
    testbarratio->SetMarkerStyle(34);
    testbarratio->SetStats(kFALSE);
    
    TF1* line = new TF1("line", "[0]", 0.6, 1.2);
    line->SetParameter(0, 1.0);
    line->SetLineColor(1);
    line->SetLineWidth(3);
    line->SetLineStyle(7);

    TLegend *legend = new TLegend(0.3758, 0.6981, 0.9362, 0.8743);
    legend->AddEntry(ratio, "K*^{0} Reconstructed/Decay", "lp");
    legend->AddEntry(barratio, "#bar{K}*^{0} Reconstructed/Decay", "lp");
    //legend->AddEntry(testratio, "(K*^{0} Reconstructed)/(#bar{K}*^{0} at Decay)", "p");
    legend->SetBorderSize(0);
    legend->SetTextSizePixels(29);
    legend->SetMargin(0.12);

    TGraph* oldgratio = new TGraph(ratio);
    Double_t* xgratio = oldgratio->GetX();
    Double_t* ygratio = oldgratio->GetY();
    TGraph* gratio = new TGraph(70, &xgratio[80], &ygratio[80]);
    gratio->SetTitle("");
    gratio->GetXaxis()->SetTitle("Invariant Mass (GeV/c^{2})");
    gratio->GetXaxis()->SetTitleFont(42);
    gratio->GetXaxis()->SetTitleSize(0.06);
    gratio->GetXaxis()->SetLabelSize(0.05);
    gratio->GetXaxis()->SetTitleOffset(1.05);
    gratio->GetXaxis()->SetRangeUser(0.6, 1.2);
    gratio->GetYaxis()->SetTitle("#frac{Reconstructed}{At Decay Point}");
    gratio->GetYaxis()->SetTitleFont(42);
    gratio->GetYaxis()->SetTitleSize(0.06);
    gratio->GetYaxis()->SetLabelSize(0.05);
    gratio->GetYaxis()->SetTitleOffset(1.0);
    gratio->GetYaxis()->SetAxisColor(kBlue);
    gratio->GetYaxis()->SetLabelColor(kBlue);
    gratio->GetYaxis()->SetRangeUser(0.1, 4000.0);
    gratio->SetMarkerSize(1.6);
    gratio->SetMarkerColor(4);
    gratio->SetMarkerStyle(22);
    gratio->SetLineColor(4);
    gratio->SetLineWidth(4);
    TGraph* oldgbarratio = new TGraph(barratio);
    Double_t* xgbarratio = oldgbarratio->GetX();
    Double_t* ygbarratio = oldgbarratio->GetY();
    TGraph* gbarratio = new TGraph(71, &xgbarratio[80], &ygbarratio[80]);
    gbarratio->SetMarkerSize(1.6);
    gbarratio->SetMarkerStyle(34);
    gbarratio->SetMarkerColor(kRed-4);
    gbarratio->SetLineColor(kRed-4);
    gbarratio->SetLineStyle(1);
    gbarratio->SetLineWidth(4);

    TCanvas *cratio = new TCanvas("cratio", "cratio",50, 50, 800, 600);
    //line->Draw("SAME"); 
    //testratio->Draw("P SAME");
    //testbarratio->Draw("P SAME");
    //legend->Draw("SAME");
    
    //TGaxis *axis = new TGaxis(0.6, 0.01 , 0.6, 100.0, 0.01, 100, 510, "G");
    //axis->SetLineColor(kBlue);
    //axis->SetTextColor(kBlue);
    //axis->Draw("SAME");

    
    Float_t histmax[4] = {1.3*totHist1->GetMaximum(), 1.3*totHist2->GetMaximum(), 1.3*totBarHist1->GetMaximum(), 1.3*totBarHist2->GetMaximum()};
    Float_t max = TMath::MaxElement(4, histmax);
    Float_t scale = 1000.0/max;
    

    TPad* overlay = new TPad("overlay", "", 0,0,1,1);
    overlay->SetFillColor(0);
    overlay->SetFillStyle(4000);
    overlay->SetFrameFillStyle(4000);
    overlay->SetMargin(0.1403, 0.1122, 0.1497, 0.06);
    overlay->SetTickx(0);
    overlay->Draw();
    //overlay->SetLogy();
    overlay->cd();
    totHist1->SetTitle("");
    totHist1->SetLineColor(kOrange-7);
    totHist1->SetLineWidth(3);
    totHist1->SetStats(kFALSE);
    totHist1->GetYaxis()->SetRangeUser(0.0, 125000);
    totHist1->GetYaxis()->SetTitleFont(42);
    totHist1->GetYaxis()->SetLabelFont(42);
    totHist1->GetYaxis()->SetLabelSize(0.05);
    totHist1->GetYaxis()->SetLabelOffset(0.01);
    totHist1->GetYaxis()->SetTitleSize(0.06);
    totHist1->GetXaxis()->SetLabelOffset(99);
    totHist1->GetXaxis()->SetTitle("");
    totHist1->GetXaxis()->SetLabelSize(0.05);
    totHist1->GetXaxis()->SetRangeUser(0.6, 1.2);
    totHist1->GetXaxis()->SetTitleOffset(99);
    totHist1->GetXaxis()->SetTickSize(0.0);
    //totHist1->GetXaxis()->Set
    totHist2->SetLineColor(kBlack);
    totHist2->SetLineWidth(3);
    totHist2->SetStats(kFALSE);
    //totHist1->Scale(scale);
    //totHist2->Scale(scale);
    totHist1->Draw("Y+ H");
    totHist2->Draw("H SAME");

    totBarHist2->SetLineColor(kSpring-7);
    totBarHist2->SetLineWidth(4);
    totBarHist2->SetLineStyle(2);
    totBarHist2->SetStats(kFALSE);
    totBarHist1->SetLineColor(kGreen-3);
    totBarHist1->SetLineWidth(4);
    totBarHist1->SetLineStyle(2);
    totBarHist1->SetStats(kFALSE);
    totBarHist1->Draw("H SAME");
    totBarHist2->Draw("H SAME");

    TLegend* legend2 = new TLegend(0.6144, 0.6667, 0.8613, 0.9267);
    legend2->AddEntry(totHist2, "K*^{0} decay", "l");
    legend2->AddEntry(totHist1, "K*^{0} recon.", "l");
    legend2->AddEntry(totBarHist2, "#bar{K}*^{0} decay", "l");
    legend2->AddEntry(totBarHist1, "#bar{K}*^{0} recon.", "l");
    legend2->SetMargin(0.25);
    legend2->SetTextSizePixels(26);
    legend2->SetTextAlign(22);
    legend2->Draw("SAME");

    TLegend* legend3 = new TLegend(0.2852, 0.7356, 0.4636, 0.8909);
    legend3->AddEntry(gratio, "K*^{0} ratio", "lp");
    legend3->AddEntry(gbarratio, "#bar{K}*^{0} ratio", "lp");
    legend3->SetMargin(0.33);
    legend3->SetTextSizePixels(28);
    legend3->SetBorderSize(0);
    legend3->SetTextAlign(22);
    legend3->Draw("SAME");

    TPad* pad = new TPad("pad", "", 0,0,1,1);
    pad->SetFillStyle(4000);
    pad->SetFrameFillStyle(4000);
    pad->SetMargin(0.1403, 0.1122, 0.1497, 0.06);
    pad->Draw();
    pad->cd();
    pad->SetLogy(kTRUE);
    gratio->Draw("ALP SAME");
    //ratio->Draw("P SAME");
    TLine* vertline1 = new TLine(0.7, gPad->GetUymin(), 0.7, 4000.0);
    vertline1->SetLineWidth(2);
    vertline1->SetLineStyle(7);
    vertline1->SetLineColor(16);
    vertline1->Draw("SAME");
    TLine* vertline2 = new TLine(1.1, gPad->GetUymin(), 1.1, 4000.0);
    vertline2->SetLineWidth(2);
    vertline2->SetLineStyle(7);
    vertline2->SetLineColor(16);
    vertline2->Draw("SAME"); 
    legend2->Draw("SAME");
    
    pad->cd();
    pad->Draw();
    gratio->Draw("LP SAME");
    gbarratio->Draw("LP SAME");
  
    printf("kstar decay: %f, kstarbar decay: %f, percent: %f\n", totHist2->Integral(), totBarHist2->Integral(), totHist2->Integral()/totBarHist2->Integral() - 1.0);
    printf("kstar recon: %f, kstarbar recon: %f, percent:%f\n", totHist1->Integral(), totBarHist1->Integral(), totHist1->Integral()/totBarHist1->Integral() - 1.0);

/*   
    TH1D* invratio = (TH1D*)totHist2->Clone("invratio");
    invratio->Divide(totHist1);
    invratio->SetTitle("");
    invratio->GetXaxis()->SetTitle("Invariant Mass (GeV/c^{2})");
    invratio->GetXaxis()->SetTitleFont(42);
    invratio->GetXaxis()->SetTitleSize(0.06);
    invratio->GetXaxis()->SetLabelSize(0.05);
    invratio->GetXaxis()->SetTitleOffset(1.05);
    invratio->GetXaxis()->SetRangeUser(0.6, 1.35);
    invratio->GetYaxis()->SetTitle("#frac{At Decay Point}{Reconstructed}");
    invratio->GetYaxis()->SetTitleFont(42);
    invratio->GetYaxis()->SetTitleSize(0.06);
    invratio->GetYaxis()->SetLabelSize(0.05);
    invratio->GetYaxis()->SetTitleOffset(1.2);
    invratio->GetYaxis()->SetRangeUser(0.005, 100.0);
    invratio->SetMarkerColor(4);
    invratio->SetMarkerSize(1.6);
    invratio->SetMarkerStyle(22);
    invratio->SetStats(kFALSE);

    TH1D* invtestratio = (TH1D*)totBarHist2->Clone("invtestratio");
    invtestratio->Divide(totHist1);
    invtestratio->SetTitle("");
    invtestratio->GetXaxis()->SetTitle("Invariant Mass (GeV/c^{2})");
    invtestratio->GetXaxis()->SetTitleFont(42);
    invtestratio->GetXaxis()->SetTitleSize(0.06);
    invtestratio->GetXaxis()->SetLabelSize(0.05);
    invtestratio->GetXaxis()->SetTitleOffset(1.05);
    invtestratio->GetXaxis()->SetRangeUser(0.6, 1.35);
    invtestratio->GetYaxis()->SetTitle("#frac{Reconstructed}{At Decay Point}");
    invtestratio->GetYaxis()->SetTitleFont(42);
    invtestratio->GetYaxis()->SetTitleSize(0.06);
    invtestratio->GetYaxis()->SetLabelSize(0.05);
    invtestratio->GetYaxis()->SetTitleOffset(1.2);
    invtestratio->GetYaxis()->SetRangeUser(0.15, 11);
    invtestratio->SetMarkerColor(2);
    invtestratio->SetMarkerSize(1.6);
    invtestratio->SetMarkerStyle(22);
    invtestratio->SetStats(kFALSE);

    TH1D* invbarratio = (TH1D*)totBarHist2->Clone("invbarratio");
    invbarratio->Divide(totBarHist1);
    invbarratio->SetTitle("");
    invbarratio->GetXaxis()->SetTitle("Invariant Mass (GeV/c^{2})");
    invbarratio->GetXaxis()->SetTitleFont(42);
    invbarratio->GetXaxis()->SetTitleSize(0.06);
    invbarratio->GetXaxis()->SetLabelSize(0.05);
    invbarratio->GetXaxis()->SetTitleOffset(1.05);
    invbarratio->GetXaxis()->SetRangeUser(0.6, 1.35);
    invbarratio->GetYaxis()->SetTitle("#frac{Reconstructed}{At Decay Point}");
    invbarratio->GetYaxis()->SetTitleFont(42);
    invbarratio->GetYaxis()->SetTitleSize(0.06);
    invbarratio->GetYaxis()->SetLabelSize(0.05);
    invbarratio->GetYaxis()->SetTitleOffset(1.2);
    invbarratio->SetMarkerColor(kMagenta-6);
    invbarratio->SetMarkerSize(1.6);
    invbarratio->SetMarkerStyle(34);
    invbarratio->SetStats(kFALSE);

    TLegend *invlegend = new TLegend(0.3758, 0.6981, 0.9362, 0.8743);
    invlegend->AddEntry(invratio, "K*^{0} Decay/Reconstructed", "p");
    invlegend->AddEntry(invbarratio, "#bar{K}*^{0} Decay/Reconstructed", "p");
    invlegend->AddEntry(invtestratio, "(#bar{K}*^{0} at Decay)/(K*^{0} Reconstructed)", "p");
    invlegend->SetBorderSize(0);
    invlegend->SetTextSizePixels(29);
    invlegend->SetMargin(0.12);


    TCanvas *cinvratio = new TCanvas("cinvratio", "cinvratio",50,50, 600, 600);
    cinvratio->SetMargin(0.1695, 0.0369, 0.1501, 0.0471);
    gPad->SetLogy(kTRUE);
    invratio->Draw("P");
    line->Draw("SAME");
    invbarratio->Draw("P SAME");
    invtestratio->Draw("P SAME");
    //testbarratio->Draw("P SAME");
    invlegend->Draw("SAME");

*/

    TFile* fratio = new TFile("test_kstar_ratio.root", "RECREATE");
    ratio->Write();
    gratio->Write();
    gbarratio->Write();
    cratio->Write();
    //cinvratio->Write();
    barratio->Write();
    testratio->Write();
    testbarratio->Write();
    totHist1->Write();
    totHist2->Write();
    totBarHist1->Write();
    totBarHist2->Write();
    
}
