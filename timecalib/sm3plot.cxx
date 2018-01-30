void sm3plot(){

    //TFile* file = new TFile("~/TimeCalibration/LHC16ktest/withbc/257139");
    TFile* file = new TFile("looseE_results.root");
    TList* list = file->Get("chistolist");
    TH3F* histo = list->FindObject("SupMod3Detail");
    TH2D* sections[12];

    for(int i = 0; i < 12; i++){
        int low = i*100+3456;
        int high = (i+1)*100+3456;
        if(high > 4607) high = 4607;
        histo->GetZaxis()->SetRangeUser(low, high);
        //histo->GetYaxis()->SetRangeUser(0.0, 5.0);
        sections[i] = (TH2D*)histo->Project3D("yx");
        sections[i]->SetName(Form("projection_%i", i));
        sections[i]->SetTitle(Form("Cell Range: %d - %d", low, high));
        //sections[i]->SetTitle(Form("Cell: %d", 842+low));
        sections[i]->GetXaxis()->SetTitle("E (GeV)");
        sections[i]->GetYaxis()->SetTitle("Corrected Time (ns)");
        sections[i]->SetStats(0);
        histo->GetZaxis()->SetRange(0,1152);
    }
     
    TCanvas* c = new TCanvas("c", "c", 50, 50, 600, 600);
    c->Divide(3,4);
    for(int i =0; i<12; i++){
        c->cd(i+1)->SetLogz();
        sections[i]->Draw("colz");
    }   

    int smLow = 4300;
    int smHigh = 900;
    histo->GetZaxis()->SetRangeUser(smLow-1+.01, smLow-1+.01);
    TH2D* smallsection = (TH2D*)histo->Project3D("yx");
    smallsection->SetName("small");
//    smallsection->SetTitle(Form("Cell Range %d - %d", 3456+smLow, 3456+smHigh));
    smallsection->SetTitle(Form("Cell: %d", smLow-1));
    smallsection->GetXaxis()->SetTitle("E (GeV)");
    smallsection->GetYaxis()->SetTitle("Corrected Time (ns)");
   TCanvas* c1 = new TCanvas("c1", "c1", 60, 60, 600, 600);
    c1->cd()->SetLogz();
    smallsection->Draw("colz");

    TCanvas *cenergy = new TCanvas("ce", "ce", 50, 50, 600, 600);
    cenergy->cd()->SetLogy();
    TH1D *energy = smallsection->ProjectionX();
    energy->Draw();


    histo->GetZaxis()->SetRange(0, 1152);
    TH2D* test = (TH2D*)histo->Project3D("yz");
    TCanvas* c2 = new TCanvas("c2", "c2", 70, 70, 600, 600);
    c2->cd()->SetLogz();
    test->Draw("colz");


    TH2D* all = (TH2D*)histo->Project3D("yx");
    all->SetTitle("Time Vs. Energy for SM3, all BC merged");
    all->GetXaxis()->SetTitle("E (GeV)");
    all->GetYaxis()->SetTitle("Corrected Time (ns)");
    all->SetStats(0);
    TCanvas *call = new TCanvas("call", "call", 50, 50, 600, 600);
    call->cd()->SetLogz();
    all->Draw("colz");


}
