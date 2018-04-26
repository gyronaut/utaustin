void makeLSCorrections(string inputFile){
    TFile* input = new TFile(inputFile.c_str());
    TH2D* hPhi2Dpeak = input->Get("hPhi2Dpeak");
    TH2D* hPhi2DLside = input->Get("hPhi2DLside");
    TH2D* hPhi2DRside = input->Get("hPhi2DRside");
    TH2D* hKK2Dpeak = input->Get("hKK2Dpeak");
    TH2D* hKK2DLside = input->Get("hKK2DLside");
    TH2D* hKK2DRside = input->Get("hKK2DRside");

    TH2D* trigDistSameUS = input->Get("fTrigSameUSDist");
    TH2D* trigDistSameLS = input->Get("fTrigSameLSDist");

    hPhi2Dpeak->SetName("uncorrectedhPhi2Dpeak");

    Float_t totalTrigUS = trigDistSameUS->Integral(trigDistSameUS->GetXaxis()->FindBin(4.0), trigDistSameUS->GetXaxis()->FindBin(8.0));
    Float_t totalTrigLS = trigDistSameLS->Integral(trigDistSameLS->GetXaxis()->FindBin(4.0), trigDistSameLS->GetXaxis()->FindBin(8.0));
    hPhi2Dpeak->Scale(1.0/totalTrigUS);
    hPhi2DLside->Scale(1.0/totalTrigUS);
    hPhi2DRside->Scale(1.0/totalTrigUS);
    hKK2Dpeak->Scale(1.0/totalTrigLS);
    hKK2DLside->Scale(1.0/totalTrigLS);
    hKK2DRside->Scale(1.0/totalTrigLS);



    //Float_t leftscale = hPhi2DLside->Integral(1, hPhi2DLside->GetXaxis()->GetNbins(), 1, hPhi2DLside->GetYaxis()->GetNbins())/hKK2DLside->Integral(1, hPhi2DLside->GetXaxis()->GetNbins(), 1, hPhi2DLside->GetYaxis()->GetNbins());
    Float_t leftscale = hPhi2DLside->Integral(hPhi2DLside->GetXaxis()->FindBin(-1.2), hPhi2DLside->GetXaxis()->FindBin(1.2), 1, hPhi2DLside->GetYaxis()->GetNbins())/hKK2DLside->Integral(hPhi2DLside->GetXaxis()->FindBin(-1.2), hPhi2DLside->GetXaxis()->FindBin(1.2), 1, hPhi2DLside->GetYaxis()->GetNbins());
    TH2D* LLSsubhPhi2DLside = hPhi2DLside->Clone("LLSsubhPhi2DLside");
    TH2D* LLSsubhPhi2Dpeak = hPhi2Dpeak->Clone("LLSsubhPhi2Dpeak");
    LLSsubhPhi2DLside->Add(hKK2DLside, -1.0*leftscale);
    //LLSsubhPhi2DLside->Divide(hKK2DLside);
    //LLSsubhPhi2DLside->Scale(1.0/leftscale);
    LLSsubhPhi2Dpeak->Add(hKK2Dpeak, -1.0*leftscale);

    TH1D* LLSsubhPhi2DLside_deta = LLSsubhPhi2DLside->ProjectionX("LLSsubhPhi2DLside_deta", 1, LLSsubhPhi2DLside->GetYaxis()->GetNbins());
    TH1D* LLSsubhPhi2DLside_dphi = LLSsubhPhi2DLside->ProjectionY("LLSsubhPhi2DLside_dphi", LLSsubhPhi2DLside->GetXaxis()->FindBin(-1.2), LLSsubhPhi2DLside->GetXaxis()->FindBin(1.2));

    //Float_t rightscale = hPhi2DRside->Integral(1, hPhi2DRside->GetXaxis()->GetNbins(), 1, hPhi2DRside->GetYaxis()->GetNbins())/hKK2DRside->Integral(1, hPhi2DRside->GetXaxis()->GetNbins(), 1, hPhi2DRside->GetYaxis()->GetNbins());
    Float_t rightscale = hPhi2DRside->Integral(hPhi2DRside->GetXaxis()->FindBin(-1.2), hPhi2DRside->GetXaxis()->FindBin(1.2), 1, hPhi2DRside->GetYaxis()->GetNbins())/hKK2DRside->Integral(hPhi2DRside->GetXaxis()->FindBin(-1.2), hPhi2DRside->GetXaxis()->FindBin(1.2), 1, hPhi2DRside->GetYaxis()->GetNbins());
    TH2D* RLSsubhPhi2DRside = hPhi2DRside->Clone("RLSsubhPhi2DRside");
    TH2D* RLSsubhPhi2Dpeak = hPhi2Dpeak->Clone("RLSsubhPhi2Dpeak");
    RLSsubhPhi2DRside->Add(hKK2DRside, -1.0*rightscale);
    //RLSsubhPhi2DRside->Divide(hKK2DRside);
    //RLSsubhPhi2DRside->Scale(1.0/rightscale);
    RLSsubhPhi2Dpeak->Add(hKK2Dpeak, -1.0*rightscale);

    TH1D* RLSsubhPhi2DRside_deta = RLSsubhPhi2DRside->ProjectionX("RLSsubhPhi2DRside_deta", 1, RLSsubhPhi2DRside->GetYaxis()->GetNbins());
    TH1D* RLSsubhPhi2DRside_dphi = RLSsubhPhi2DRside->ProjectionY("RLSsubhPhi2DRside_dphi", RLSsubhPhi2DRside->GetXaxis()->FindBin(-1.2), RLSsubhPhi2DRside->GetXaxis()->FindBin(1.2));

    TH1D* RLSsubhPhi2Dpeak_deta = RLSsubhPhi2Dpeak->ProjectionX("RLSsubhPhi2Dpeak_deta", 1, RLSsubhPhi2Dpeak->GetYaxis()->GetNbins());
    TH1D* RLSsubhPhi2Dpeak_dphi = RLSsubhPhi2Dpeak->ProjectionY("RLSsubhPhi2Dpeak_dphi", RLSsubhPhi2Dpeak->GetXaxis()->FindBin(-1.2), RLSsubhPhi2Dpeak->GetXaxis()->FindBin(1.2));


    TH1D* scales = new TH1D("scales", "scales", 2, -1, 1);
    scales->SetBinContent(1, leftscale);
    scales->SetBinContent(2, rightscale);

    TH2D* rebinRLSsubhPhi2Dpeak = RLSsubhPhi2Dpeak->Clone("rebinRLSsubhPhi2Dpeak");
    rebinRLSsubhPhi2Dpeak->Rebin2D(2, 2);

    TFile* output = new TFile(Form("LS_%s", inputFile.c_str()), "RECREATE");
    LLSsubhPhi2DLside->Write();
    LLSsubhPhi2DLside_deta->Write();
    LLSsubhPhi2DLside_dphi->Write();
    LLSsubhPhi2Dpeak->Write();
    RLSsubhPhi2DRside->Write();
    RLSsubhPhi2DRside_deta->Write();
    RLSsubhPhi2DRside_dphi->Write();
    RLSsubhPhi2Dpeak->Write();
    RLSsubhPhi2Dpeak_deta->Write();
    RLSsubhPhi2Dpeak_dphi->Write();
    rebinRLSsubhPhi2Dpeak->Write();
    scales->Write();
    hPhi2Dpeak->Write();

}
