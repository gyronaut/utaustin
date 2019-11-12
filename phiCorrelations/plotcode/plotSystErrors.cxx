void plotSystErrors(){
    //systematic erros for dphi distribution
    Double_t dphistat[3] = {2.2, 3.2, 7.3};
    Double_t dphisyst_masspeak[3] = {0.8, 1.3, 2.7};
    Double_t dphisyst_massrsb[3] = {1.1, 0.9, 2.6};
    Double_t dphisyst_masslsb[3] = {0.2, 0.3, 1.2};
    Double_t dphisyst_kls[3] = {0.9, 0.6, 3.0};
    Double_t dphisyst_ksignal[3] = {0.8, 0.8, 0.8};
    Double_t dphisyst_pid[3] = {1.0, 1.0, 1.0};
    Double_t dphisyst_TOFpid[3] = {1.4, 1.4, 1.4};
    Double_t dphisyst_trackingeff[3] = {4.2, 4.2, 4.2};
    Double_t dphisyst_total[3] = {0.0, 0.0, 0.0};

    //systematic errors for near-side yields
    Double_t nearstat[3] = {13.2, 11.55, 26.4};
    Double_t nearsyst_masspeak[3] = {4.3, 4.3, 3.5};
    Double_t nearsyst_massrsb[3] = {3.6, 2.0, 7.6};
    Double_t nearsyst_masslsb[3] = {1.4, 1.4, 2.2};
    Double_t nearsyst_kls[3] = {1.5, 0.8, 6.1};
    Double_t nearsyst_ksignal[3] = {0.8, 0.8, 0.8};
    Double_t nearsyst_pid[3] = {1.0, 1.0, 1.0};
    Double_t nearsyst_TOFpid[3] = {1.4, 1.4, 1.4};
    Double_t nearsyst_trackingeff[3] = {4.2, 4.2, 4.2};
    Double_t nearsyst_ue[3] = {3.2, 5.8, 8.0};
    Double_t nearsyst_total[3] = {0.0, 0.0, 0.0};

    //systematic errors for away-side yields
    Double_t awaystat[3] = {10.6, 8.9, 23.2};
    Double_t awaysyst_masspeak[3] = {3.2, 3.3, 5.0};
    Double_t awaysyst_massrsb[3] = {2.1, 1.0, 4.2};
    Double_t awaysyst_masslsb[3] = {0.9, 0.4, 2.1};
    Double_t awaysyst_kls[3] = {1.2, 0.5, 3.9};
    Double_t awaysyst_ksignal[3] = {0.8, 0.8, 0.8};
    Double_t awaysyst_pid[3] = {1.0, 1.0, 1.0};
    Double_t awaysyst_TOFpid[3] = {1.4, 1.4, 1.4};
    Double_t awaysyst_trackingeff[3] = {4.2, 4.2, 4.2};
    Double_t awaysyst_ue[3] = {8.8, 4.9, 7.7};
    Double_t awaysyst_total[3] = {0.0, 0.0, 0.0};

    //systematic errors for U.E. yields
    Double_t uestat[3] = {1.6, 2.3, 3.3};
    Double_t uesyst_masspeak[3] = {0.5, 0.3, 1.3};
    Double_t uesyst_massrsb[3] = {0.1, 0.1, 0.4};
    Double_t uesyst_masslsb[3] = {0.8, 0.3, 0.8};
    Double_t uesyst_kls[3] = {0.9, 0.6, 2.6};
    Double_t uesyst_ksignal[3] = {0.8, 0.8, 0.8};
    Double_t uesyst_pid[3] = {1.0, 1.0, 1.0};
    Double_t uesyst_TOFpid[3] = {1.4, 1.4, 1.4};
    Double_t uesyst_trackingeff[3] = {4.2, 4.2, 4.2};
    Double_t uesyst_ue[3] = {0.7, 0.7, 2.0};
    Double_t uesyst_total[3] = {0.0, 0.0, 0.0};

    //systematic errors for total yields
    Double_t totalstat[3] = {0.6, 0.8, 1.9};
    Double_t totalsyst_masspeak[3] = {0.4, 0.6, 1.6};
    Double_t totalsyst_massrsb[3] = {0.1, 0.1, 0.5};
    Double_t totalsyst_masslsb[3] = {0.7, 0.4, 0.6};
    Double_t totalsyst_kls[3] = {0.8, 0.6, 2.6};
    Double_t totalsyst_ksignal[3] = {0.8, 0.8, 0.8};
    Double_t totalsyst_pid[3] = {1.0, 1.0, 1.0};
    Double_t totalsyst_TOFpid[3] = {1.4, 1.4, 1.4};
    Double_t totalsyst_trackingeff[3] = {4.2, 4.2, 4.2};
    Double_t totalsyst_ue[3] = {0.0, 0.0, 0.0};
    Double_t totalsyst_total[3] = {0.0, 0.0, 0.0};

    for(int i = 0; i < 3; i++){
        dphisyst_total[i] += TMath::Power(dphisyst_masspeak[i], 2.0);
        dphisyst_total[i] += TMath::Power(dphisyst_massrsb[i], 2.0);
        dphisyst_total[i] += TMath::Power(dphisyst_masslsb[i], 2.0);
        dphisyst_total[i] += TMath::Power(dphisyst_kls[i], 2.0);
        dphisyst_total[i] += TMath::Power(dphisyst_ksignal[i], 2.0);
        dphisyst_total[i] += TMath::Power(dphisyst_pid[i], 2.0);
        dphisyst_total[i] += TMath::Power(dphisyst_TOFpid[i], 2.0);
        dphisyst_total[i] += TMath::Power(dphisyst_trackingeff[i], 2.0);

        dphisyst_total[i] = TMath::Sqrt(dphisyst_total[i]);
        printf("total syst %d: %4.2f\n", i, dphisyst_total[i]);

        nearsyst_total[i] += TMath::Power(nearsyst_masspeak[i], 2.0);
        nearsyst_total[i] += TMath::Power(nearsyst_massrsb[i], 2.0);
        nearsyst_total[i] += TMath::Power(nearsyst_masslsb[i], 2.0);
        nearsyst_total[i] += TMath::Power(nearsyst_kls[i], 2.0);
        nearsyst_total[i] += TMath::Power(nearsyst_ksignal[i], 2.0);
        nearsyst_total[i] += TMath::Power(nearsyst_pid[i], 2.0);
        nearsyst_total[i] += TMath::Power(nearsyst_TOFpid[i], 2.0);
        nearsyst_total[i] += TMath::Power(nearsyst_ue[i], 2.0);
        nearsyst_total[i] += TMath::Power(nearsyst_trackingeff[i], 2.0);

        nearsyst_total[i] = TMath::Sqrt(nearsyst_total[i]);
        printf("total near syst %d: %4.2f\n", i, nearsyst_total[i]);

        awaysyst_total[i] += TMath::Power(awaysyst_masspeak[i], 2.0);
        awaysyst_total[i] += TMath::Power(awaysyst_massrsb[i], 2.0);
        awaysyst_total[i] += TMath::Power(awaysyst_masslsb[i], 2.0);
        awaysyst_total[i] += TMath::Power(awaysyst_kls[i], 2.0);
        awaysyst_total[i] += TMath::Power(awaysyst_ksignal[i], 2.0);
        awaysyst_total[i] += TMath::Power(awaysyst_pid[i], 2.0);
        awaysyst_total[i] += TMath::Power(awaysyst_TOFpid[i], 2.0);
        awaysyst_total[i] += TMath::Power(awaysyst_ue[i], 2.0);
        awaysyst_total[i] += TMath::Power(awaysyst_trackingeff[i], 2.0);

        awaysyst_total[i] = TMath::Sqrt(awaysyst_total[i]);
        printf("total away syst %d: %4.2f\n", i, awaysyst_total[i]);

        uesyst_total[i] += TMath::Power(uesyst_masspeak[i], 2.0);
        uesyst_total[i] += TMath::Power(uesyst_massrsb[i], 2.0);
        uesyst_total[i] += TMath::Power(uesyst_masslsb[i], 2.0);
        uesyst_total[i] += TMath::Power(uesyst_kls[i], 2.0);
        uesyst_total[i] += TMath::Power(uesyst_ksignal[i], 2.0);
        uesyst_total[i] += TMath::Power(uesyst_pid[i], 2.0);
        uesyst_total[i] += TMath::Power(uesyst_TOFpid[i], 2.0);
        uesyst_total[i] += TMath::Power(uesyst_ue[i], 2.0);
        uesyst_total[i] += TMath::Power(uesyst_trackingeff[i], 2.0);

        uesyst_total[i] = TMath::Sqrt(uesyst_total[i]);
        printf("total ue syst %d: %4.2f\n", i, uesyst_total[i]);

        totalsyst_total[i] += TMath::Power(totalsyst_masspeak[i], 2.0);
        totalsyst_total[i] += TMath::Power(totalsyst_massrsb[i], 2.0);
        totalsyst_total[i] += TMath::Power(totalsyst_masslsb[i], 2.0);
        totalsyst_total[i] += TMath::Power(totalsyst_kls[i], 2.0);
        totalsyst_total[i] += TMath::Power(totalsyst_ksignal[i], 2.0);
        totalsyst_total[i] += TMath::Power(totalsyst_pid[i], 2.0);
        totalsyst_total[i] += TMath::Power(totalsyst_TOFpid[i], 2.0);
        totalsyst_total[i] += TMath::Power(totalsyst_ue[i], 2.0);
        totalsyst_total[i] += TMath::Power(totalsyst_trackingeff[i], 2.0);

        totalsyst_total[i] = TMath::Sqrt(totalsyst_total[i]);
        printf("total total syst %d: %4.2f\n", i, totalsyst_total[i]);

    }
    
    
    //plotting dphi systematics
    TH1D *hdphistat = new TH1D("hdphistat", "hdphistat", 3, 0, 3);
    hdphistat->GetXaxis()->SetBinLabel(1, "50-80%");
    hdphistat->GetXaxis()->SetBinLabel(2, "20-50%");
    hdphistat->GetXaxis()->SetBinLabel(3, "0-20%");
    hdphistat->GetXaxis()->SetLabelSize(0.06);
    hdphistat->GetXaxis()->SetTitle("Multiplicity");
    hdphistat->GetXaxis()->SetTitleOffset(1.2);
    hdphistat->GetYaxis()->SetTitle("Relative Uncertainty (for each #Delta#varphi bin)");
    hdphistat->SetTitle("Systematic Error Sources for #Delta#varphi (h-#phi) Correlation");
    
    TH1D* hdphisyst_masspeak = (TH1D*)hdphistat->Clone("hdphisyst_masspeak");
    hdphisyst_masspeak->SetLineColor(kRed+1);
    TH1D* hdphisyst_massrsb = (TH1D*)hdphistat->Clone("hdphisyst_massrsb");
    hdphisyst_massrsb->SetLineColor(kOrange+1);
    TH1D* hdphisyst_masslsb = (TH1D*)hdphistat->Clone("hdphisyst_masslsb");
    hdphisyst_masslsb->SetLineColor(kOrange+2);
    TH1D* hdphisyst_kls = (TH1D*)hdphistat->Clone("hdphisyst_kls");
    hdphisyst_kls->SetLineColor(kSpring+5);
    TH1D* hdphisyst_ksignal = (TH1D*)hdphistat->Clone("hdphisyst_ksignal");
    hdphisyst_ksignal->SetLineColor(kGreen+1);
    TH1D* hdphisyst_pid = (TH1D*)hdphistat->Clone("hdphisyst_pid");
    hdphisyst_pid->SetLineColor(kCyan+1);
    TH1D* hdphisyst_TOFpid = (TH1D*)hdphistat->Clone("hdphisyst_TOFpid");
    hdphisyst_TOFpid->SetLineColor(kAzure+1);
    TH1D* hdphisyst_trackingeff = (TH1D*)hdphistat->Clone("hdphisyst_trackingeff");
    hdphisyst_trackingeff->SetLineColor(kViolet+10);
    TH1D* hdphisyst_total = (TH1D*)hdphistat->Clone("hdphisyst_total");
    hdphisyst_total->SetLineColor(kBlack);
    hdphisyst_total->SetLineWidth(2);

    hdphistat->SetLineColor(kBlue+2);
    hdphistat->SetLineStyle(7);
    hdphistat->SetLineWidth(2);
    hdphistat->SetStats(kFALSE);

    for(int i = 0; i<3; i++){
        hdphistat->SetBinContent(3-i, dphistat[i]/100.0);
        hdphisyst_masspeak->SetBinContent(3-i, dphisyst_masspeak[i]/100.0);
        hdphisyst_massrsb->SetBinContent(3-i, dphisyst_massrsb[i]/100.0);
        hdphisyst_masslsb->SetBinContent(3-i, dphisyst_masslsb[i]/100.0);
        hdphisyst_kls->SetBinContent(3-i, dphisyst_kls[i]/100.0);
        hdphisyst_ksignal->SetBinContent(3-i, dphisyst_ksignal[i]/100.0);
        hdphisyst_pid->SetBinContent(3-i, dphisyst_pid[i]/100.0);
        hdphisyst_TOFpid->SetBinContent(3-i, dphisyst_TOFpid[i]/100.0);
        hdphisyst_trackingeff->SetBinContent(3-i, dphisyst_trackingeff[i]/100.0);
        hdphisyst_total->SetBinContent(3-i, dphisyst_total[i]/100.0);
    }

    TLegend* dphileg = new TLegend(0.14, 0.47, 0.34, 0.87);
    dphileg->AddEntry(hdphistat, "Avg. Statistical Error", "l");
    dphileg->AddEntry(hdphisyst_masspeak, "Mass peak region", "l");
    dphileg->AddEntry(hdphisyst_massrsb, "right SB", "l");
    dphileg->AddEntry(hdphisyst_masslsb, "left SB", "l");
    dphileg->AddEntry(hdphisyst_kls, "LS scale factor", "l");
    dphileg->AddEntry(hdphisyst_ksignal, "signal scale factor", "l");
    dphileg->AddEntry(hdphisyst_pid, "TPC PID cut", "l");
    dphileg->AddEntry(hdphisyst_TOFpid, "TOF PID cut", "l");
    dphileg->AddEntry(hdphisyst_trackingeff, "tracking eff.", "l");
    dphileg->AddEntry(hdphisyst_total, "Total Syst. Error", "l");

    TCanvas* cdphi = new TCanvas("cdphi", "cdphi", 50, 50, 700, 500);
    cdphi->cd();
    hdphistat->Draw("HIST");
    hdphisyst_masspeak->Draw("HIST SAME");
    hdphisyst_massrsb->Draw("HIST SAME");
    hdphisyst_masslsb->Draw("HIST SAME");
    hdphisyst_kls->Draw("HIST SAME");
    hdphisyst_ksignal->Draw("HIST SAME");
    hdphisyst_pid->Draw("HIST SAME");
    hdphisyst_TOFpid->Draw("HIST SAME");
    hdphisyst_trackingeff->Draw("HIST SAME");
    hdphisyst_total->Draw("HIST SAME");
    dphileg->Draw();


    //plotting near systematics
    TH1D *hnearstat = new TH1D("hnearstat", "hnearstat", 3, 0, 3);
    hnearstat->GetXaxis()->SetBinLabel(1, "50-80%");
    hnearstat->GetXaxis()->SetBinLabel(2, "20-50%");
    hnearstat->GetXaxis()->SetBinLabel(3, "0-20%");
    hnearstat->GetXaxis()->SetLabelSize(0.06);
    hnearstat->GetXaxis()->SetTitle("Multiplicity");
    hnearstat->GetXaxis()->SetTitleOffset(1.2);
    hnearstat->GetYaxis()->SetTitle("Relative Uncertainty");
    hnearstat->SetTitle("Systematic Error Sources for Near-side (h-#phi) Yield");
    
    TH1D* hnearsyst_masspeak = (TH1D*)hnearstat->Clone("hnearsyst_masspeak");
    hnearsyst_masspeak->SetLineColor(kRed+1);
    TH1D* hnearsyst_massrsb = (TH1D*)hnearstat->Clone("hnearsyst_massrsb");
    hnearsyst_massrsb->SetLineColor(kOrange+1);
    TH1D* hnearsyst_masslsb = (TH1D*)hnearstat->Clone("hnearsyst_masslsb");
    hnearsyst_masslsb->SetLineColor(kOrange+2);
    TH1D* hnearsyst_kls = (TH1D*)hnearstat->Clone("hnearsyst_kls");
    hnearsyst_kls->SetLineColor(kSpring+5);
    TH1D* hnearsyst_ksignal = (TH1D*)hnearstat->Clone("hnearsyst_ksignal");
    hnearsyst_ksignal->SetLineColor(kGreen+1);
    TH1D* hnearsyst_pid = (TH1D*)hnearstat->Clone("hnearsyst_pid");
    hnearsyst_pid->SetLineColor(kCyan+1);
    TH1D* hnearsyst_TOFpid = (TH1D*)hnearstat->Clone("hnearsyst_TOFpid");
    hnearsyst_TOFpid->SetLineColor(kAzure+1);
    TH1D* hnearsyst_trackingeff = (TH1D*)hnearstat->Clone("hnearsyst_trackingeff");
    hnearsyst_trackingeff->SetLineColor(kViolet+10);
    TH1D* hnearsyst_ue = (TH1D*)hnearstat->Clone("hnearsyst_ue");
    hnearsyst_ue->SetLineColor(kViolet+1);
    TH1D* hnearsyst_total = (TH1D*)hnearstat->Clone("hnearsyst_total");
    hnearsyst_total->SetLineColor(kBlack);
    hnearsyst_total->SetLineWidth(2);

    hnearstat->SetLineColor(kBlue+2);
    hnearstat->SetLineStyle(7);
    hnearstat->SetLineWidth(2);
    hnearstat->SetStats(kFALSE);

    for(int i = 0; i<3; i++){
        hnearstat->SetBinContent(3-i, nearstat[i]/100.0);
        hnearsyst_masspeak->SetBinContent(3-i, nearsyst_masspeak[i]/100.0);
        hnearsyst_massrsb->SetBinContent(3-i, nearsyst_massrsb[i]/100.0);
        hnearsyst_masslsb->SetBinContent(3-i, nearsyst_masslsb[i]/100.0);
        hnearsyst_kls->SetBinContent(3-i, nearsyst_kls[i]/100.0);
        hnearsyst_ksignal->SetBinContent(3-i, nearsyst_ksignal[i]/100.0);
        hnearsyst_pid->SetBinContent(3-i, nearsyst_pid[i]/100.0);
        hnearsyst_TOFpid->SetBinContent(3-i, nearsyst_TOFpid[i]/100.0);
        hnearsyst_trackingeff->SetBinContent(3-i, nearsyst_trackingeff[i]/100.0);
        hnearsyst_ue->SetBinContent(3-i, nearsyst_ue[i]/100.0);
        hnearsyst_total->SetBinContent(3-i, nearsyst_total[i]/100.0);
    }

    TLegend* nearleg = new TLegend(0.14, 0.47, 0.34, 0.87);
    nearleg->AddEntry(hnearstat, "Statistical Error", "l");
    nearleg->AddEntry(hnearsyst_masspeak, "Mass peak region", "l");
    nearleg->AddEntry(hnearsyst_massrsb, "right SB", "l");
    nearleg->AddEntry(hnearsyst_masslsb, "left SB", "l");
    nearleg->AddEntry(hnearsyst_kls, "LS scale factor", "l");
    nearleg->AddEntry(hnearsyst_ksignal, "signal scale factor", "l");
    nearleg->AddEntry(hnearsyst_pid, "TPC PID cut", "l");
    nearleg->AddEntry(hnearsyst_TOFpid, "TOF PID cut", "l");
    nearleg->AddEntry(hnearsyst_trackingeff, "tracking eff.", "l");
    nearleg->AddEntry(hnearsyst_ue, "U.E estimation", "l");
    nearleg->AddEntry(hnearsyst_total, "Total Syst. Error", "l");

    TCanvas* cnear = new TCanvas("cnear", "cnear", 50, 50, 700, 500);
    cnear->cd();
    hnearstat->Draw("HIST");
    hnearsyst_masspeak->Draw("HIST SAME");
    hnearsyst_massrsb->Draw("HIST SAME");
    hnearsyst_masslsb->Draw("HIST SAME");
    hnearsyst_kls->Draw("HIST SAME");
    hnearsyst_ksignal->Draw("HIST SAME");
    hnearsyst_pid->Draw("HIST SAME");
    hnearsyst_TOFpid->Draw("HIST SAME");
    hnearsyst_trackingeff->Draw("HIST SAME");
    hnearsyst_ue->Draw("HIST SAME");
    hnearsyst_total->Draw("HIST SAME");
    nearleg->Draw();

    //plotting away systematic
    TH1D *hawaystat = new TH1D("hawaystat", "hawaystat", 3, 0, 3);
    hawaystat->GetXaxis()->SetBinLabel(1, "50-80%");
    hawaystat->GetXaxis()->SetBinLabel(2, "20-50%");
    hawaystat->GetXaxis()->SetBinLabel(3, "0-20%");
    hawaystat->GetXaxis()->SetLabelSize(0.06);
    hawaystat->GetXaxis()->SetTitle("Multiplicity");
    hawaystat->GetXaxis()->SetTitleOffset(1.2);
    hawaystat->GetYaxis()->SetTitle("Relative Uncertainty");
    hawaystat->SetTitle("Systematic Error Sources for Away-side (h-#phi) Yield");
    
    TH1D* hawaysyst_masspeak = (TH1D*)hawaystat->Clone("hawaysyst_masspeak");
    hawaysyst_masspeak->SetLineColor(kRed+1);
    TH1D* hawaysyst_massrsb = (TH1D*)hawaystat->Clone("hawaysyst_massrsb");
    hawaysyst_massrsb->SetLineColor(kOrange+1);
    TH1D* hawaysyst_masslsb = (TH1D*)hawaystat->Clone("hawaysyst_masslsb");
    hawaysyst_masslsb->SetLineColor(kOrange+2);
    TH1D* hawaysyst_kls = (TH1D*)hawaystat->Clone("hawaysyst_kls");
    hawaysyst_kls->SetLineColor(kSpring+5);
    TH1D* hawaysyst_ksignal = (TH1D*)hawaystat->Clone("hawaysyst_ksignal");
    hawaysyst_ksignal->SetLineColor(kGreen+1);
    TH1D* hawaysyst_pid = (TH1D*)hawaystat->Clone("hawaysyst_pid");
    hawaysyst_pid->SetLineColor(kCyan+1);
    TH1D* hawaysyst_TOFpid = (TH1D*)hawaystat->Clone("hawaysyst_TOFpid");
    hawaysyst_TOFpid->SetLineColor(kAzure+1);
    TH1D* hawaysyst_trackingeff = (TH1D*)hawaystat->Clone("hawaysyst_trackingeff");
    hawaysyst_trackingeff->SetLineColor(kViolet+10);
    TH1D* hawaysyst_ue = (TH1D*)hawaystat->Clone("hawaysyst_ue");
    hawaysyst_ue->SetLineColor(kViolet+1);
    TH1D* hawaysyst_total = (TH1D*)hawaystat->Clone("hawaysyst_total");
    hawaysyst_total->SetLineColor(kBlack);
    hawaysyst_total->SetLineWidth(2);

    hawaystat->SetLineColor(kBlue+2);
    hawaystat->SetLineStyle(7);
    hawaystat->SetLineWidth(2);
    hawaystat->SetStats(kFALSE);

    for(int i = 0; i<3; i++){
        hawaystat->SetBinContent(3-i, awaystat[i]/100.0);
        hawaysyst_masspeak->SetBinContent(3-i, awaysyst_masspeak[i]/100.0);
        hawaysyst_massrsb->SetBinContent(3-i, awaysyst_massrsb[i]/100.0);
        hawaysyst_masslsb->SetBinContent(3-i, awaysyst_masslsb[i]/100.0);
        hawaysyst_kls->SetBinContent(3-i, awaysyst_kls[i]/100.0);
        hawaysyst_ksignal->SetBinContent(3-i, awaysyst_ksignal[i]/100.0);
        hawaysyst_pid->SetBinContent(3-i, awaysyst_pid[i]/100.0);
        hawaysyst_TOFpid->SetBinContent(3-i, awaysyst_TOFpid[i]/100.0);
        hawaysyst_trackingeff->SetBinContent(3-i, awaysyst_trackingeff[i]/100.0);
        hawaysyst_ue->SetBinContent(3-i, awaysyst_ue[i]/100.0);
        hawaysyst_total->SetBinContent(3-i, awaysyst_total[i]/100.0);
    }

    TLegend* awayleg = new TLegend(0.14, 0.47, 0.34, 0.87);
    awayleg->AddEntry(hawaystat, "Statistical Error", "l");
    awayleg->AddEntry(hawaysyst_masspeak, "Mass peak region", "l");
    awayleg->AddEntry(hawaysyst_massrsb, "right SB", "l");
    awayleg->AddEntry(hawaysyst_masslsb, "left SB", "l");
    awayleg->AddEntry(hawaysyst_kls, "LS scale factor", "l");
    awayleg->AddEntry(hawaysyst_ksignal, "signal scale factor", "l");
    awayleg->AddEntry(hawaysyst_pid, "TPC PID cut", "l");
    awayleg->AddEntry(hawaysyst_TOFpid, "TOF PID cut", "l");
    awayleg->AddEntry(hawaysyst_trackingeff, "tracking eff.", "l");
    awayleg->AddEntry(hawaysyst_ue, "U.E estimation", "l");
    awayleg->AddEntry(hawaysyst_total, "Total Syst. Error", "l");

    TCanvas* caway = new TCanvas("caway", "caway", 50, 50, 700, 500);
    caway->cd();
    hawaystat->Draw("HIST");
    hawaysyst_masspeak->Draw("HIST SAME");
    hawaysyst_massrsb->Draw("HIST SAME");
    hawaysyst_masslsb->Draw("HIST SAME");
    hawaysyst_kls->Draw("HIST SAME");
    hawaysyst_ksignal->Draw("HIST SAME");
    hawaysyst_pid->Draw("HIST SAME");
    hawaysyst_TOFpid->Draw("HIST SAME");
    hawaysyst_trackingeff->Draw("HIST SAME");
    hawaysyst_ue->Draw("HIST SAME");
    hawaysyst_total->Draw("HIST SAME");
    awayleg->Draw();

    //plotting ue systematics
    TH1D *huestat = new TH1D("huestat", "huestat", 3, 0, 3);
    huestat->GetXaxis()->SetBinLabel(1, "50-80%");
    huestat->GetXaxis()->SetBinLabel(2, "20-50%");
    huestat->GetXaxis()->SetBinLabel(3, "0-20%");
    huestat->GetXaxis()->SetLabelSize(0.06);
    huestat->GetXaxis()->SetTitle("Multiplicity");
    huestat->GetXaxis()->SetTitleOffset(1.2);
    huestat->GetYaxis()->SetTitle("Relative Uncertainty");
    huestat->SetTitle("Systematic Error Sources for Underlying Event (h-#phi) Yield");
    
    TH1D* huesyst_masspeak = (TH1D*)huestat->Clone("huesyst_masspeak");
    huesyst_masspeak->SetLineColor(kRed+1);
    TH1D* huesyst_massrsb = (TH1D*)huestat->Clone("huesyst_massrsb");
    huesyst_massrsb->SetLineColor(kOrange+1);
    TH1D* huesyst_masslsb = (TH1D*)huestat->Clone("huesyst_masslsb");
    huesyst_masslsb->SetLineColor(kOrange+2);
    TH1D* huesyst_kls = (TH1D*)huestat->Clone("huesyst_kls");
    huesyst_kls->SetLineColor(kSpring+5);
    TH1D* huesyst_ksignal = (TH1D*)huestat->Clone("huesyst_ksignal");
    huesyst_ksignal->SetLineColor(kGreen+1);
    TH1D* huesyst_pid = (TH1D*)huestat->Clone("huesyst_pid");
    huesyst_pid->SetLineColor(kCyan+1);
    TH1D* huesyst_TOFpid = (TH1D*)huestat->Clone("huesyst_TOFpid");
    huesyst_TOFpid->SetLineColor(kAzure+1);
    TH1D* huesyst_trackingeff = (TH1D*)huestat->Clone("huesyst_trackingeff");
    huesyst_trackingeff->SetLineColor(kViolet+10);
    TH1D* huesyst_ue = (TH1D*)huestat->Clone("huesyst_ue");
    huesyst_ue->SetLineColor(kViolet+1);
    TH1D* huesyst_total = (TH1D*)huestat->Clone("huesyst_total");
    huesyst_total->SetLineColor(kBlack);
    huesyst_total->SetLineWidth(2);

    huestat->SetLineColor(kBlue+2);
    huestat->SetLineStyle(7);
    huestat->SetLineWidth(2);
    huestat->SetStats(kFALSE);

    for(int i = 0; i<3; i++){
        huestat->SetBinContent(3-i, uestat[i]/100.0);
        huesyst_masspeak->SetBinContent(3-i, uesyst_masspeak[i]/100.0);
        huesyst_massrsb->SetBinContent(3-i, uesyst_massrsb[i]/100.0);
        huesyst_masslsb->SetBinContent(3-i, uesyst_masslsb[i]/100.0);
        huesyst_kls->SetBinContent(3-i, uesyst_kls[i]/100.0);
        huesyst_ksignal->SetBinContent(3-i, uesyst_ksignal[i]/100.0);
        huesyst_pid->SetBinContent(3-i, uesyst_pid[i]/100.0);
        huesyst_TOFpid->SetBinContent(3-i, uesyst_TOFpid[i]/100.0);
        huesyst_trackingeff->SetBinContent(3-i, uesyst_trackingeff[i]/100.0);
        huesyst_ue->SetBinContent(3-i, uesyst_ue[i]/100.0);
        huesyst_total->SetBinContent(3-i, uesyst_total[i]/100.0);
    }

    TLegend* ueleg = new TLegend(0.14, 0.47, 0.34, 0.87);
    ueleg->AddEntry(huestat, "Statistical Error", "l");
    ueleg->AddEntry(huesyst_masspeak, "Mass peak region", "l");
    ueleg->AddEntry(huesyst_massrsb, "right SB", "l");
    ueleg->AddEntry(huesyst_masslsb, "left SB", "l");
    ueleg->AddEntry(huesyst_kls, "LS scale factor", "l");
    ueleg->AddEntry(huesyst_ksignal, "signal scale factor", "l");
    ueleg->AddEntry(huesyst_pid, "TPC PID cut", "l");
    ueleg->AddEntry(huesyst_TOFpid, "TOF PID cut", "l");
    ueleg->AddEntry(huesyst_trackingeff, "tracking eff.", "l");
    ueleg->AddEntry(huesyst_ue, "U.E estimation", "l");
    ueleg->AddEntry(huesyst_total, "Total Syst. Error", "l");

    TCanvas* cue = new TCanvas("cue", "cue", 50, 50, 700, 500);
    cue->cd();
    huestat->Draw("HIST");
    huesyst_masspeak->Draw("HIST SAME");
    huesyst_massrsb->Draw("HIST SAME");
    huesyst_masslsb->Draw("HIST SAME");
    huesyst_kls->Draw("HIST SAME");
    huesyst_ksignal->Draw("HIST SAME");
    huesyst_pid->Draw("HIST SAME");
    huesyst_TOFpid->Draw("HIST SAME");
    huesyst_trackingeff->Draw("HIST SAME");
    huesyst_ue->Draw("HIST SAME");
    huesyst_total->Draw("HIST SAME");
    ueleg->Draw();

    //plotting total systematics
    TH1D *htotalstat = new TH1D("htotalstat", "htotalstat", 3, 0, 3);
    htotalstat->GetXaxis()->SetBinLabel(1, "50-80%");
    htotalstat->GetXaxis()->SetBinLabel(2, "20-50%");
    htotalstat->GetXaxis()->SetBinLabel(3, "0-20%");
    htotalstat->GetXaxis()->SetLabelSize(0.06);
    htotalstat->GetXaxis()->SetTitle("Multiplicity");
    htotalstat->GetXaxis()->SetTitleOffset(1.2);
    htotalstat->GetYaxis()->SetTitle("Relative Uncertainty");
    htotalstat->SetTitle("Systematic Error Sources for Total (h-#phi) Yield");
    
    TH1D* htotalsyst_masspeak = (TH1D*)htotalstat->Clone("htotalsyst_masspeak");
    htotalsyst_masspeak->SetLineColor(kRed+1);
    TH1D* htotalsyst_massrsb = (TH1D*)htotalstat->Clone("htotalsyst_massrsb");
    htotalsyst_massrsb->SetLineColor(kOrange+1);
    TH1D* htotalsyst_masslsb = (TH1D*)htotalstat->Clone("htotalsyst_masslsb");
    htotalsyst_masslsb->SetLineColor(kOrange+2);
    TH1D* htotalsyst_kls = (TH1D*)htotalstat->Clone("htotalsyst_kls");
    htotalsyst_kls->SetLineColor(kSpring+5);
    TH1D* htotalsyst_ksignal = (TH1D*)htotalstat->Clone("htotalsyst_ksignal");
    htotalsyst_ksignal->SetLineColor(kGreen+1);
    TH1D* htotalsyst_pid = (TH1D*)htotalstat->Clone("htotalsyst_pid");
    htotalsyst_pid->SetLineColor(kCyan+1);
    TH1D* htotalsyst_TOFpid = (TH1D*)htotalstat->Clone("htotalsyst_TOFpid");
    htotalsyst_TOFpid->SetLineColor(kAzure+1);
    TH1D* htotalsyst_trackingeff = (TH1D*)htotalstat->Clone("htotalsyst_trackingeff");
    htotalsyst_trackingeff->SetLineColor(kViolet+10);
    TH1D* htotalsyst_ue = (TH1D*)htotalstat->Clone("htotalsyst_ue");
    htotalsyst_ue->SetLineColor(kViolet+1);
    TH1D* htotalsyst_total = (TH1D*)htotalstat->Clone("htotalsyst_total");
    htotalsyst_total->SetLineColor(kBlack);
    htotalsyst_total->SetLineWidth(2);

    htotalstat->SetLineColor(kBlue+2);
    htotalstat->SetLineStyle(7);
    htotalstat->SetLineWidth(2);
    htotalstat->SetStats(kFALSE);

    for(int i = 0; i<3; i++){
        htotalstat->SetBinContent(3-i, totalstat[i]/100.0);
        htotalsyst_masspeak->SetBinContent(3-i, totalsyst_masspeak[i]/100.0);
        htotalsyst_massrsb->SetBinContent(3-i, totalsyst_massrsb[i]/100.0);
        htotalsyst_masslsb->SetBinContent(3-i, totalsyst_masslsb[i]/100.0);
        htotalsyst_kls->SetBinContent(3-i, totalsyst_kls[i]/100.0);
        htotalsyst_ksignal->SetBinContent(3-i, totalsyst_ksignal[i]/100.0);
        htotalsyst_pid->SetBinContent(3-i, totalsyst_pid[i]/100.0);
        htotalsyst_TOFpid->SetBinContent(3-i, totalsyst_TOFpid[i]/100.0);
        htotalsyst_trackingeff->SetBinContent(3-i, totalsyst_trackingeff[i]/100.0);
        htotalsyst_ue->SetBinContent(3-i, totalsyst_ue[i]/100.0);
        htotalsyst_total->SetBinContent(3-i, totalsyst_total[i]/100.0);
    }

    TLegend* totalleg = new TLegend(0.14, 0.47, 0.34, 0.87);
    totalleg->AddEntry(htotalstat, "Statistical Error", "l");
    totalleg->AddEntry(htotalsyst_masspeak, "Mass peak region", "l");
    totalleg->AddEntry(htotalsyst_massrsb, "right SB", "l");
    totalleg->AddEntry(htotalsyst_masslsb, "left SB", "l");
    totalleg->AddEntry(htotalsyst_kls, "LS scale factor", "l");
    totalleg->AddEntry(htotalsyst_ksignal, "signal scale factor", "l");
    totalleg->AddEntry(htotalsyst_pid, "TPC PID cut", "l");
    totalleg->AddEntry(htotalsyst_TOFpid, "TOF PID cut", "l");
    totalleg->AddEntry(htotalsyst_trackingeff, "tracking eff.", "l");
    totalleg->AddEntry(htotalsyst_total, "Total Syst. Error", "l");

    TCanvas* ctotal = new TCanvas("ctotal", "ctotal", 50, 50, 700, 500);
    ctotal->cd();
    htotalstat->Draw("HIST");
    htotalsyst_masspeak->Draw("HIST SAME");
    htotalsyst_massrsb->Draw("HIST SAME");
    htotalsyst_masslsb->Draw("HIST SAME");
    htotalsyst_kls->Draw("HIST SAME");
    htotalsyst_ksignal->Draw("HIST SAME");
    htotalsyst_pid->Draw("HIST SAME");
    htotalsyst_TOFpid->Draw("HIST SAME");
    htotalsyst_trackingeff->Draw("HIST SAME");
    htotalsyst_total->Draw("HIST SAME");
    totalleg->Draw();

}
