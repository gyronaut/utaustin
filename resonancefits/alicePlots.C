void alicePlots(){
    
    TFile* alice = new TFile("~/Downloads/HEPData-ins1288320-v1-root.root");
    alice->cd("Table 16");
    TGraph* aliceData = Graph1D_y1;
    TH1F* hist = Hist1D_y1;
    TH1F* stat = Hist1D_y1_e1;
    TH1F* syst = Hist1D_y1_e2;
    TGraphAsymmErrors* graph2 = (TGraphAsymmErrors*)aliceData->Clone("graph2");

    Int_t numPts = aliceData->GetN();
    Double_t x, y;
    for(int i = 0; i<numPts; i++){
        aliceData->GetPoint(i, x, y);
        aliceData->SetPoint(i, x, (y - 0.89581));
        graph2->SetPoint(i, x, (y- 0.89581));
        hist->SetBinContent(i+1, hist->GetBinContent(i+1) - 0.89581);
        hist->SetBinError(i+1, stat->GetBinContent(i+1));
        graph2->SetPointEXhigh(i, 0.1);
        graph2->SetPointEXlow(i, 0.1);
    }

    graph2->SetLineColor(kBlue-10);
    graph2->SetLineWidth(2);
    graph2->SetMarkerColor(kBlue-10);
    graph2->SetFillColor(kBlue-10);
    hist->SetLineColor(kBlue-2);
    hist->SetLineWidth(2);

    aliceData->SetTitle("");
    aliceData->GetYaxis()->SetTitle("Mass - Vacuum Mass (GeV/c^{2})");
    aliceData->GetYaxis()->SetTitleSize(0.06);
    aliceData->GetYaxis()->SetLabelSize(0.04);
    aliceData->GetYaxis()->SetTitleOffset(1.65);
    aliceData->GetYaxis()->SetTitleFont(42);
    aliceData->GetYaxis()->SetLabelFont(42);
    aliceData->GetXaxis()->SetTitle("p_{T} (GeV/c)");
    aliceData->GetXaxis()->SetTitleSize(0.06);
    aliceData->GetXaxis()->SetLabelSize(0.05);
    aliceData->GetXaxis()->SetTitleFont(42);
    aliceData->GetXaxis()->SetLabelFont(42);
    aliceData->SetMarkerStyle(29);
    aliceData->SetMarkerSize(2.5);
    aliceData->SetMarkerColor(kBlue-2);
    aliceData->SetLineColor(kBlue-2);

    aliceData->GetYaxis()->SetRangeUser(-0.02, 0.015);
    aliceData->GetXaxis()->SetRangeUser(0, 5);

    TFile* phsd = new TFile("~/utaustin/resonancefits/finalplotting/20170721_KKbarAdded2_fixedwidth42_recon_pf100_scaled_error05.root");
    TH1D* mass = phsd->Get("kstar0mass");
    mass->SetName("mass");
    mass->SetMarkerStyle(22);
    mass->SetMarkerSize(2.5);
    mass->SetMarkerColor(2);
    mass->SetLineColor(2);

    TF1* line = new TF1("line", "[0]", 0.0, 5.0);
    line->SetParameter(0, 0.0);
    line->SetLineColor(1);
    line->SetLineStyle(7);
    line->SetLineWidth(3);

    for(int j = 0; j<mass->GetNbinsX(); j++){
        mass->SetBinContent(j+1, (mass->GetBinContent(j+1) - 0.892));
    }

    TExec *exec1 = new TExec("exec1", "gStyle->SetErrorX(0.1)");
    TExec *exec2 = new TExec("exec2", "gStyle->SetErrorX(0.5)");
    
    TCanvas *c = new TCanvas ("c", "c", 50, 50, 650, 600);
    c->cd()->SetMargin(0.1997, 0.0369, 0.1396, 0.0681);
    aliceData->Draw("APX");
    //exec1->Draw();
    graph2->Draw("SAME P2");
    //exec2->Draw();
    hist->Draw("SAME E1");
    line->Draw("SAME");
    mass->Draw("SAME P E1");
    aliceData->Draw("SAME PX");
    

    TLegend* legend = new TLegend(0.5836, 0.1815, 0.9489, 0.3438);
    legend->SetMargin(0.2);
    legend->SetTextSizePixels(20);
    legend->AddEntry(aliceData, "ALICE data, 0-20%", "p");
    legend->AddEntry(mass, "Fit to PHSD with Fixed Vacuum #Gamma", "p");
    legend->Draw("SAME"); 
  
    TPaveText* text = new TPaveText(0.2554, 0.7243, 0.6006, 0.9162, "NDC");
    text->AddText("(K*^{0} + #bar{K}*^{0})");
    text->AddText("Pb-Pb #sqrt{s_{NN}} = 2.76 TeV");
    text->GetLine(0)->SetTextSizePixels(36);
    text->GetLine(1)->SetTextSizePixels(24);
    text->SetTextFont(42);
    text->SetBorderSize(0);
    text->SetFillStyle(0);
    text->Draw();
}
