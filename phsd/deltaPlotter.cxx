void deltaPlotter(string inputname){

    TFile* file = new TFile(inputname.c_str());

    TH2D* massallDelta = (TH2D*)file->Get("massallDelta");
    TH2D* massallantiDelta = (TH2D*)file->Get("massallantiDelta");
    TH2D* massTotalDelta = (TH2D*)file->Get("massTotalDelta");
    TH1D* massallDeltaPT[4];
    TH1D* massallantiDeltaPT[4];
    TH1D* massTotalDeltaPT[4];

    TPaveText* label[4];
    TLegend* legend = new TLegend(0.1253, 0.3700, 0.5344, 0.5433);

    for(int ipt = 0; ipt < 4; ipt++){
        massallDeltaPT[ipt] = massallDelta->ProjectionX(Form("allDeltaptbin%d", ipt+1), ipt+1, ipt+1);
        massallDeltaPT[ipt]->SetStats(kFALSE);
        massallDeltaPT[ipt]->SetLineColor(kBlue+1);
        massallDeltaPT[ipt]->SetLineWidth(2);
        massallantiDeltaPT[ipt] = massallantiDelta->ProjectionX(Form("allantiDeltaptbin%d", ipt+1), ipt+1, ipt+1);
        massallantiDeltaPT[ipt]->SetStats(kFALSE);
        massallantiDeltaPT[ipt]->SetLineColor(kRed+1);
        massallantiDeltaPT[ipt]->SetLineWidth(2);
        massTotalDeltaPT[ipt] = massTotalDelta->ProjectionX(Form("TotalDeltaptbin%d", ipt+1), ipt+1, ipt+1);
        massTotalDeltaPT[ipt]->SetStats(kFALSE);
        massTotalDeltaPT[ipt]->SetLineColor(kBlack);
        massTotalDeltaPT[ipt]->SetLineWidth(2);
        label[ipt] = new TPaveText(0.1545, 0.6984, 0.4050, 0.8505, "NDC");
        label[ipt]->AddText(Form("%2.1f < p_{T}^{#Delta} < %2.1f", (float)(ipt+1)*0.5, (float)(ipt+2)*0.5));
        label[ipt]->SetBorderSize(0);
        label[ipt]->SetFillStyle(0);
        //label[ipt]->GetLine(0)->SetTextSizePixels(20);
        if(ipt == 0){
            legend->AddEntry(massTotalDeltaPT[ipt], "Delta baryons + antibaryons", "l");
            legend->AddEntry(massallDeltaPT[ipt], "Delta baryons", "l");
            legend->AddEntry(massallantiDeltaPT[ipt], "Delta antibaryons", "l");
        }
    }

    TCanvas* c1 = new TCanvas("c1", "c1", 50, 50, 1000, 800);
    c1->Divide(2,2);
    for(int ipt = 0; ipt<4; ipt++){
        c1->cd(ipt+1);
        massTotalDeltaPT[ipt]->Draw("H");
        massallDeltaPT[ipt]->Draw("H SAME");
        massallantiDeltaPT[ipt]->Draw("H SAME");
        label[ipt]->Draw();
        if(ipt ==1) legend->Draw("SAME");
    }
}
