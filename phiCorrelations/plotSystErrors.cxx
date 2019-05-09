void plotSystErrors(){
    Double_t dphistat[3] = {2.2, 3.2, 7.3};
    Double_t dphisyst_masspeak[3] = {0.8, 1.3, 2.7};
    Double_t dphisyst_massrsb[3] = {1.1, 0.9, 2.6};
    Double_t dphisyst_masslsb[3] = {0.2, 0.3, 1.2};
    Double_t dphisyst_kls[3] = {0.9, 0.6, 3.0};
    Double_t dphisyst_ksignal[3] = {0.8, 0.8, 0.8};
    Double_t dphisyst_pid[3] = {1.0, 1.0, 1.0};
    Double_t dphisyst_trackingeff[3] = {5.2, 5.2, 5.2};
    
    Double_t dphisyst_total[3] = {0.0, 0.0, 0.0};

    for(int i = 0; i < 3; i++){
        dphisyst_total[i] += TMath::Power(dphisyst_masspeak[i], 2.0);
        dphisyst_total[i] += TMath::Power(dphisyst_massrsb[i], 2.0);
        dphisyst_total[i] += TMath::Power(dphisyst_masslsb[i], 2.0);
        dphisyst_total[i] += TMath::Power(dphisyst_kls[i], 2.0);
        dphisyst_total[i] += TMath::Power(dphisyst_ksignal[i], 2.0);
        dphisyst_total[i] += TMath::Power(dphisyst_pid[i], 2.0);
        dphisyst_total[i] += TMath::Power(dphisyst_trackingeff[i], 2.0);

        dphisyst_total[i] = TMath::Sqrt(dphisyst_total[i]);
        printf("total syst %d: %4.2f\n", i, dphisyst_total[i]);
    }
    
    
    TH1D *hdphistat = new TH1D("hdphistat", "hdphistat", 3, 0, 3);
    hdphistat->GetXaxis()->SetBinLabel(1, "50-80%");
    hdphistat->GetXaxis()->SetBinLabel(2, "20-50%");
    hdphistat->GetXaxis()->SetBinLabel(3, "0-20%");

    TH1D* hdphisyst_masspeak = (TH1D*)hdphistat->Clone("hdphisyst_masspeak");
    TH1D* hdphisyst_massrsb = (TH1D*)hdphistat->Clone("hdphisyst_massrsb");
    TH1D* hdphisyst_masslsb = (TH1D*)hdphistat->Clone("hdphisyst_masslsb");
    TH1D* hdphisyst_kls = (TH1D*)hdphistat->Clone("hdphisyst_kls");
    TH1D* hdphisyst_ksignal = (TH1D*)hdphistat->Clone("hdphisyst_ksignal");
    TH1D* hdphisyst_pid = (TH1D*)hdphistat->Clone("hdphisyst_pid");
    TH1D* hdphisyst_trackingeff = (TH1D*)hdphistat->Clone("hdphisyst_trackingeff");
    TH1D* hdphisyst_total = (TH1D*)hdphistat->Clone("hdphisyst_total");
    
    for(int i = 1; i<=3; i++){
        hdphistat->SetBinContent(i, dphistat[i-1]);
        hdphisyst_masspeak->SetBinContent(i, dphisyst_masspeak[i-1]);
        hdphisyst_massrsb->SetBinContent(i, dphisyst_massrsb[i-1]);
        hdphisyst_masslsb->SetBinContent(i, dphisyst_masslsb[i-1]);
        hdphisyst_kls->SetBinContent(i, dphisyst_kls[i-1]);
        hdphisyst_ksignal->SetBinContent(i, dphisyst_ksignal[i-1]);
        hdphisyst_pid->SetBinContent(i, dphisyst_pid[i-1]);
        hdphisyst_trackingeff->SetBinContent(i, dphisyst_trackingeff[i-1]);
        hdphisyst_total->SetBinContent(i, dphisyst_total[i-1]);
    }


}
