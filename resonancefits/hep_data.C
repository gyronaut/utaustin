{
//=========Macro generated from canvas: Canvas_1/Canvas_1
//=========  (Tue Jun 13 18:50:56 2017) by ROOT version5.34/30
   TCanvas *Canvas_1 = new TCanvas("Canvas_1", "Canvas_1",260,94,536,324);
   Canvas_1->SetHighLightColor(2);
   Canvas_1->Range(-0.2875,-0.0009390626,5.5875,0.005951562);
   Canvas_1->SetFillColor(0);
   Canvas_1->SetBorderMode(0);
   Canvas_1->SetBorderSize(2);
   Canvas_1->SetFrameBorderMode(0);
   Canvas_1->SetFrameBorderMode(0);
   Double_t xAxis1[9] = {0.3, 0.8, 1.2, 1.6, 2, 2.5, 3, 4, 5}; 
   
   TH1F *Hist1D_y1_e1 = new TH1F("Hist1D_y1_e1","doi:10.17182/hepdata.66630.v1/t16",8, xAxis1);
   Hist1D_y1_e1->SetBinContent(1,0.0023);
   Hist1D_y1_e1->SetBinContent(2,0.0025);
   Hist1D_y1_e1->SetBinContent(3,0.0017);
   Hist1D_y1_e1->SetBinContent(4,0.0015);
   Hist1D_y1_e1->SetBinContent(5,0.0012);
   Hist1D_y1_e1->SetBinContent(6,0.0016);
   Hist1D_y1_e1->SetBinContent(7,0.0013);
   Hist1D_y1_e1->SetBinContent(8,0.0025);
   Hist1D_y1_e1->SetBinError(1,0.0023);
   Hist1D_y1_e1->SetBinError(2,0.0025);
   Hist1D_y1_e1->SetBinError(3,0.0017);
   Hist1D_y1_e1->SetBinError(4,0.0015);
   Hist1D_y1_e1->SetBinError(5,0.0012);
   Hist1D_y1_e1->SetBinError(6,0.0016);
   Hist1D_y1_e1->SetBinError(7,0.0013);
   Hist1D_y1_e1->SetBinError(8,0.0025);
   Hist1D_y1_e1->SetEntries(8);
   
   TPaveStats *ptstats = new TPaveStats(0.78,0.775,0.98,0.935,"brNDC");
   ptstats->SetName("stats");
   ptstats->SetBorderSize(1);
   ptstats->SetFillColor(0);
   ptstats->SetTextAlign(12);
   ptstats->SetTextFont(42);
   TText *text = ptstats->AddText("Hist1D_y1_e1");
   text->SetTextSize(0.0368);
   text = ptstats->AddText("Entries = 8      ");
   text = ptstats->AddText("Mean  =  2.174");
   text = ptstats->AddText("RMS   =  1.362");
   ptstats->SetOptStat(1111);
   ptstats->SetOptFit(0);
   ptstats->Draw();
   Hist1D_y1_e1->GetListOfFunctions()->Add(ptstats);
   ptstats->SetParent(Hist1D_y1_e1);

   Int_t ci;      // for color index setting
   TColor *color; // for color definition with alpha
   ci = TColor::GetColor("#000099");
   Hist1D_y1_e1->SetLineColor(ci);
   Hist1D_y1_e1->GetXaxis()->SetTitle("PT [GEV]");
   Hist1D_y1_e1->GetXaxis()->SetLabelFont(42);
   Hist1D_y1_e1->GetXaxis()->SetLabelSize(0.035);
   Hist1D_y1_e1->GetXaxis()->SetTitleSize(0.035);
   Hist1D_y1_e1->GetXaxis()->SetTitleFont(42);
   Hist1D_y1_e1->GetYaxis()->SetTitle("stat");
   Hist1D_y1_e1->GetYaxis()->SetLabelFont(42);
   Hist1D_y1_e1->GetYaxis()->SetLabelSize(0.035);
   Hist1D_y1_e1->GetYaxis()->SetTitleSize(0.035);
   Hist1D_y1_e1->GetYaxis()->SetTitleFont(42);
   Hist1D_y1_e1->GetZaxis()->SetLabelFont(42);
   Hist1D_y1_e1->GetZaxis()->SetLabelSize(0.035);
   Hist1D_y1_e1->GetZaxis()->SetTitleSize(0.035);
   Hist1D_y1_e1->GetZaxis()->SetTitleFont(42);
   Hist1D_y1_e1->Draw("");
   
   TPaveText *pt = new TPaveText(0.2456391,0.935,0.7543609,0.995,"blNDC");
   pt->SetName("title");
   pt->SetBorderSize(0);
   pt->SetFillColor(0);
   pt->SetFillStyle(0);
   pt->SetTextFont(42);
   text = pt->AddText("doi:10.17182/hepdata.66630.v1/t16");
   pt->Draw();
   Canvas_1->Modified();
   Canvas_1->cd();
   Canvas_1->SetSelected(Canvas_1);
}
