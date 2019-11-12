Float_t quotientErrors(Float_t a, Float_t aErr, Float_t b, Float_t bErr){
    return (a/b)*TMath::Sqrt(TMath::Power(aErr/a, 2.0) + TMath::Power(bErr/b, 2.0));
}

Float_t mySimpleErrInt(TH1D* hist, Int_t lowbin, Int_t highbin){
    Float_t totInt = 0.0;
    for(int i = lowbin; i<=highbin; i++){
        totInt += TMath::Power(hist->GetBinContent(i)*hist->GetXaxis()->GetBinWidth(i), 2.0);
    }
    return TMath::Sqrt(totInt);
}

Float_t myIntegral(TH1D* hist, Int_t lowbin, Int_t highbin){
    Float_t totInt = 0.0;
    for(int i = lowbin; i<=highbin; i++){
        totInt += (hist->GetBinContent(i))*(hist->GetXaxis()->GetBinWidth(i))*(hist->GetXaxis()->GetBinCenter(i))*2.0*TMath::Pi();
    }
    return totInt;
}

Float_t myIntegralErr(TH1D* hist, Int_t lowbin, Int_t highbin){
    Float_t totInt = 0.0;
    for(int i = lowbin; i<=highbin; i++){
        totInt += TMath::Power((hist->GetBinContent(i))*(hist->GetXaxis()->GetBinWidth(i))*(hist->GetXaxis()->GetBinCenter(i))*2.0*TMath::Pi(), 2.0);
    }
    return TMath::Sqrt(totInt);
}

Float_t addInQuad(TH1D* hist, Int_t lowbin, Int_t highbin){
    Float_t totInt = 0.0;
    for(int i = lowbin; i <= highbin; i++){
        totInt += TMath::Power(hist->GetBinContent(i), 2.0);
    }
    return TMath::Sqrt(totInt);
}

void plotPaperPhiRatios(){
    TFile* piFile = new TFile("~/phiStudies/pubRatioCheck/HEPData-ins1415274-v1-Table_1.root");
    TDirectory* pidir = (TDirectory*) piFile->GetDirectory("Table 1");
    TFile* kFile = new TFile("~/phiStudies/pubRatioCheck/HEPData-ins1415274-v1-Table_3.root");
    TDirectory* kdir = (TDirectory*) kFile->GetDirectory("Table 3");
    TFile* protonFile = new TFile("~/phiStudies/pubRatioCheck/HEPData-ins1415274-v1-Table_5.root");
    TDirectory* protondir = (TDirectory*) protonFile->GetDirectory("Table 5");

    TFile* chargeFile = new TFile("~/phiStudies/pubRatioCheck/HEPData-ins1657384-v1-Table_3.root");
    TDirectory* chargeDir = (TDirectory*)chargeFile->GetDirectory("Table 3");
    TH1D* chargeMB = (TH1D*)(chargeDir->Get("Hist1D_y1")->Clone("chargeMB"));

    TH1D* pihist[7];
    TH1D* pihistStat[7];
    TH1D* pihistSyst[7];
    TH1D* pihistSystUncorr[7];
    TH1D* khist[7];
    TH1D* khistSyst[7];
    TH1D* khistStat[7];
    TH1D* khistSystUncorr[7];
    TH1D* protonhist[7];
    TH1D* protonhistStat[7];
    TH1D* protonhistSyst[7];
    TH1D* protonhistSystUncorr[7];

    //setup histograms (and plot) for different particles and multiplicities
    TCanvas* cfilepi = new TCanvas("cfilepi", "cfilepi", 50, 50, 600, 600);
    TCanvas* cfilek = new TCanvas("cfilek", "cfilek", 50, 50, 600, 600);
    TCanvas* cfileproton = new TCanvas("cfileproton", "cfileproton", 50, 50, 600, 600);

    for(int i = 0; i < 7; i++){
        pihist[i] = (TH1D*)(pidir->Get(Form("Hist1D_y%d", i+1))->Clone(Form("pihist%d", i+1)));
        for(int ibin = 1; ibin <= pihist[i]->GetXaxis()->GetNbins(); ibin++){
            pihist[i]->SetBinContent(ibin, pihist[i]->GetBinContent(ibin)*2.0*TMath::Pi()*pihist[i]->GetXaxis()->GetBinCenter(ibin)*pihist[i]->GetXaxis()->GetBinWidth(ibin));         
            //pihist[i]->SetBinContent(ibin, pihist[i]->GetBinContent(ibin)*2.0*TMath::Pi()*pihist[i]->GetXaxis()->GetBinWidth(ibin));         
        }
        pihistStat[i] = (TH1D*)(pidir->Get(Form("Hist1D_y%d_e1", i+1))->Clone(Form("pihistStat%d", i+1)));
        pihistSyst[i] = (TH1D*)(pidir->Get(Form("Hist1D_y%d_e2", i+1))->Clone(Form("pihistSyst%d", i+1)));
        pihistSystUncorr[i] = (TH1D*)(pidir->Get(Form("Hist1D_y%d_e3", i+1))->Clone(Form("pihistSystUncorr%d", i+1)));
        cfilepi->cd();
        pihist[i]->Draw("HIST SAME");
    }
    for(int i = 0; i < 7; i++){
        khist[i] = (TH1D*)(kdir->Get(Form("Hist1D_y%d", i+1))->Clone(Form("khist%d", i+1)));
        for(int ibin = 1; ibin <= khist[i]->GetXaxis()->GetNbins(); ibin++){
            khist[i]->SetBinContent(ibin, khist[i]->GetBinContent(ibin)*2.0*TMath::Pi()*khist[i]->GetXaxis()->GetBinCenter(ibin)*khist[i]->GetXaxis()->GetBinWidth(ibin));
            //khist[i]->SetBinContent(ibin, khist[i]->GetBinContent(ibin)*2.0*TMath::Pi()*khist[i]->GetXaxis()->GetBinWidth(ibin));
        }

        khistStat[i] = (TH1D*)(kdir->Get(Form("Hist1D_y%d_e1", i+1))->Clone(Form("khistStat%d", i+1)));
        khistSyst[i] = (TH1D*)(kdir->Get(Form("Hist1D_y%d_e2", i+1))->Clone(Form("khistSyst%d", i+1)));
        khistSystUncorr[i] = (TH1D*)(kdir->Get(Form("Hist1D_y%d_e3", i+1))->Clone(Form("khistSystUncorr%d", i+1)));
        cfilek->cd();
        khist[i]->Draw("HIST SAME");
    }
    for(int i = 0; i < 7; i++){
        protonhist[i] = (TH1D*)(protondir->Get(Form("Hist1D_y%d", i+1))->Clone(Form("protonhist%d", i+1)));
        for(int ibin = 1; ibin <= protonhist[i]->GetXaxis()->GetNbins(); ibin++){
            protonhist[i]->SetBinContent(ibin, protonhist[i]->GetBinContent(ibin)*2.0*TMath::Pi()*protonhist[i]->GetXaxis()->GetBinCenter(ibin)*protonhist[i]->GetXaxis()->GetBinWidth(ibin));
            //protonhist[i]->SetBinContent(ibin, protonhist[i]->GetBinContent(ibin)*2.0*TMath::Pi()*protonhist[i]->GetXaxis()->GetBinWidth(ibin));
        }
        protonhistStat[i] = (TH1D*)(protondir->Get(Form("Hist1D_y%d_e1", i+1))->Clone(Form("protonhistStat%d", i+1)));
        protonhistSyst[i] = (TH1D*)(protondir->Get(Form("Hist1D_y%d_e2", i+1))->Clone(Form("protonhistSyst%d", i+1)));
        protonhistSystUncorr[i] = (TH1D*)(protondir->Get(Form("Hist1D_y%d_e3", i+1))->Clone(Form("protonhistSystUncorr%d", i+1)));
        cfileproton->cd();
        protonhist[i]->Draw("HIST SAME");
    }
    TCanvas* cspeciescomp = new TCanvas("cspeciescomp", "cspeciescomp", 50, 50, 600, 600);
    cspeciescomp->cd();
    pihist[0]->Draw("SAME HIST");
    khist[0]->Draw("SAME HIST");
    protonhist[0]->Draw("SAME HIST");
    pihist[6]->Draw("SAME HIST");
    khist[6]->Draw("SAME HIST");
    protonhist[6]->Draw("SAME HIST");

    printf("past species compare\n");
    TH1D* phihist[7];
    TH1D* phihistStat[7];
    TH1D* phihistSyst[7];
    TH1D* phihistSystUncorr[7];

    TFile* phiFilemult1 = new TFile("~/phiStudies/pubRatioCheck/HEPData-ins1418181-v1-Table_9_mult_0_5.root");
    TDirectory* phidirmult1 = (TDirectory*) phiFilemult1->GetDirectory("Table 9");
    phihist[0] = (TH1D*)(phidirmult1->Get("Hist1D_y1")->Clone("phihist1"));
    phihistStat[0] = (TH1D*)(phidirmult1->Get("Hist1D_y1_e1")->Clone("phihistStat1"));
    phihistSyst[0] = (TH1D*)(phidirmult1->Get("Hist1D_y1_e2")->Clone("phihistSyst1"));
    TFile* phiFilemult2 = new TFile("~/phiStudies/pubRatioCheck/HEPData-ins1418181-v1-Table_10_mult_5_10.root");
    TDirectory* phidirmult2 = (TDirectory*) phiFilemult2->GetDirectory("Table 10");
    phihist[1] = (TH1D*)(phidirmult2->Get("Hist1D_y1")->Clone("phihist2"));
    phihistStat[1] = (TH1D*)(phidirmult2->Get("Hist1D_y1_e1")->Clone("phihistStat2"));
    phihistSyst[1] = (TH1D*)(phidirmult2->Get("Hist1D_y1_e2")->Clone("phihistSyst2"));
    TFile* phiFilemult3 = new TFile("~/phiStudies/pubRatioCheck/HEPData-ins1418181-v1-Table_11_mult_10_20.root");
    TDirectory* phidirmult3 = (TDirectory*) phiFilemult3->GetDirectory("Table 11");
    phihist[2] = (TH1D*)(phidirmult3->Get("Hist1D_y1")->Clone("phihist3"));
    phihistStat[2] = (TH1D*)(phidirmult3->Get("Hist1D_y1_e1")->Clone("phihistStat3"));
    phihistSyst[2] = (TH1D*)(phidirmult3->Get("Hist1D_y1_e2")->Clone("phihistSyst3"));
    TFile* phiFilemult4 = new TFile("~/phiStudies/pubRatioCheck/HEPData-ins1418181-v1-Table_12_mult_20_40.root");
    TDirectory* phidirmult4 = (TDirectory*) phiFilemult4->GetDirectory("Table 12");
    phihist[3] = (TH1D*)(phidirmult4->Get("Hist1D_y1")->Clone("phihist4"));
    phihistStat[3] = (TH1D*)(phidirmult4->Get("Hist1D_y1_e1")->Clone("phihistStat4"));
    phihistSyst[3] = (TH1D*)(phidirmult4->Get("Hist1D_y1_e2")->Clone("phihistSyst4"));
    TFile* phiFilemult5 = new TFile("~/phiStudies/pubRatioCheck/HEPData-ins1418181-v1-Table_13_mult_40_60.root");
    TDirectory* phidirmult5 = (TDirectory*) phiFilemult5->GetDirectory("Table 13");
    phihist[4] = (TH1D*)(phidirmult5->Get("Hist1D_y1")->Clone("phihist5"));
    phihistStat[4] = (TH1D*)(phidirmult5->Get("Hist1D_y1_e1")->Clone("phihistStat5"));
    phihistSyst[4] = (TH1D*)(phidirmult5->Get("Hist1D_y1_e2")->Clone("phihistSyst5"));
    TFile* phiFilemult6 = new TFile("~/phiStudies/pubRatioCheck/HEPData-ins1418181-v1-Table_14_mult_60_80.root");
    TDirectory* phidirmult6 = (TDirectory*) phiFilemult6->GetDirectory("Table 14");
    phihist[5] = (TH1D*)(phidirmult6->Get("Hist1D_y1")->Clone("phihist6"));
    phihistStat[5] = (TH1D*)(phidirmult6->Get("Hist1D_y1_e1")->Clone("phihistStat6"));
    phihistSyst[5] = (TH1D*)(phidirmult6->Get("Hist1D_y1_e2")->Clone("phihistSyst6"));
    TFile* phiFilemult7 = new TFile("~/phiStudies/pubRatioCheck/HEPData-ins1418181-v1-Table_15_mult_80_100.root");
    TDirectory* phidirmult7 = (TDirectory*) phiFilemult7->GetDirectory("Table 15");
    phihist[6] = (TH1D*)(phidirmult7->Get("Hist1D_y1")->Clone("phihist7"));
    phihistStat[6] = (TH1D*)(phidirmult7->Get("Hist1D_y1_e1")->Clone("phihistStat7"));
    phihistSyst[6] = (TH1D*)(phidirmult7->Get("Hist1D_y1_e2")->Clone("phihistSyst7"));

    TH1D* phiMB = (TH1D*)phihist[6]->Clone("phiMB");
    for(int i = 1; i <= phiMB->GetXaxis()->GetNbins(); i++){
        phiMB->SetBinContent(i, .05*phihist[0]->GetBinContent(i) + .05*phihist[1]->GetBinContent(i) + .1*phihist[2]->GetBinContent(i) + 0.2*phihist[3]->GetBinContent(i) + 0.2*phihist[4]->GetBinContent(i) + 0.2*phihist[5]->GetBinContent(i) + 0.2*phihist[6]->GetBinContent(i));
    }

    for(int i = 0; i<7; i++){
        for(int ibin = 1; ibin <= phihist[i]->GetXaxis()->GetNbins(); ibin++){
            phihist[i]->SetBinContent(ibin, phihist[i]->GetBinContent(ibin)*phihist[i]->GetXaxis()->GetBinWidth(ibin));
        }
    }

    printf("past phi setup\n");

    Float_t multbins[] = {2.5, 7.5, 15.0, 30.0, 50.0, 70.0, 90.0};
    Float_t revmultbins[] = { 10, 30, 50, 70, 85, 92.5, 97.5};
    Float_t revmultbins6[] = {30, 50, 70, 85, 92.5, 97.5};
    Float_t revmultwidths[] = {10, 10, 10, 10, 5, 2.5, 2.5};
    Float_t revmultwidths6[] = {10, 10, 10, 5, 2.5, 2.5};
    Float_t histbins[] = {0.0, 20.0, 40.0, 60.0, 80.0, 90.0, 95.0, 100.0};
    
    Float_t charged[7];
    Float_t chargedErr[7];
    Float_t chargedSystErr[7];
    Float_t totcharged[7];
    Float_t pions[7];
    Float_t totpions[7];
    Float_t kaons[7];
    Float_t totkaons[7];
    Float_t protons[7];
    Float_t totprotons[7];
    Float_t phi[7];
    Float_t phiErr[7];
    Float_t phiSystErr[7];
    Float_t totphi[7];
    Float_t ratio[7];
    Float_t ratioErr[7];
    Float_t ratioSystErr[7];
    Float_t totratio[7];
    Float_t ratio2pi[7];
    Float_t totratio2pi[7];
    Float_t ratio2k[7];
    Float_t totratio2k[7];
    Float_t ratio2proton[7];
    Float_t totratio2proton[7];
    Float_t ratiop2phi[7];

    Float_t eps = 0.001;
    Float_t lowpt = 2.0;
    Float_t highpt = 4.0;

    //changing pt binning of proton hist
    Double_t ptlow[] = {0.4, 0.6, 0.8, 1.0, 1.2, 1.4, 1.6, 1.8, 2.0, 4.0};
    Double_t pthigh[] = {0.6, 0.8, 1.0, 1.2, 1.4, 1.6, 1.8, 2.0, 4.0};
    
    
    TH1D* protonRebin[7];
    TH1D* kaonRebin[7];
    TH1D* pionRebin[7];
    TH1D* phiRebin[7];
    TH1D* chargedRebin[7];
    TH1D* chargedRebinMB = new TH1D("chargedRebinMB", "chargedRebinMB", 10, ptlow);
    
    TH1D* protonRebinError[7];
    TH1D* kaonRebinError[7];
    TH1D* pionRebinError[7];
    TH1D* phiRebinError[7];
    TH1D* chargedRebinError[7];

    TH1D* protonRebinSystErr[7];
    TH1D* kaonRebinSystErr[7];
    TH1D* pionRebinSystErr[7];
    TH1D* phiRebinSystErr[7];
    TH1D* chargedRebinSystErr[7];



    TH1D* ptratiop2phi[7];

    TH1D* ptratiophi2p[7];
    TH1D* ptratiophi2k[7];
    TH1D* ptratiophi2pi[7];
    TH1D* ptratiophi2charged[7];
    
    TH1D* pi2charged[7];
    TH1D* k2charged[7];
    TH1D* proton2charged[7];
    TH1D* proton2pion[7];

    for(int imult = 0; imult<7; imult++){
        protonRebin[imult] = new TH1D(Form("protonRebin%d", imult), Form("protonRebin%d", imult), 9, ptlow);
        kaonRebin[imult] = new TH1D(Form("kaonRebin%d", imult), Form("kaonRebin%d", imult), 9, ptlow);
        pionRebin[imult] = new TH1D(Form("pionRebin%d", imult), Form("pionRebin%d", imult), 9, ptlow);
        phiRebin[imult] = new TH1D(Form("phiRebin%d", imult), Form("phiRebin%d", imult), 9, ptlow);
        chargedRebin[imult] = new TH1D(Form("chargedRebin%d", imult), Form("chargedRebin%d", imult), 9, ptlow);

        protonRebinError[imult] = new TH1D(Form("protonRebinError%d", imult), Form("protonRebinError%d", imult), 9, ptlow);
        kaonRebinError[imult] = new TH1D(Form("kaonRebinError%d", imult), Form("kaonRebinError%d", imult), 9, ptlow);
        pionRebinError[imult] = new TH1D(Form("pionRebinError%d", imult), Form("pionRebinError%d", imult), 9, ptlow);
        phiRebinError[imult] = new TH1D(Form("phiRebinError%d", imult), Form("phiRebinError%d", imult), 9, ptlow);
        chargedRebinError[imult] = new TH1D(Form("chargedRebinError%d", imult), Form("chargedRebinError%d", imult), 9, ptlow);

        protonRebinSystErr[imult] = new TH1D(Form("protonRebinSystErr%d", imult), Form("protonRebinSystErr%d", imult), 9, ptlow);
        kaonRebinSystErr[imult] = new TH1D(Form("kaonRebinSystErr%d", imult), Form("kaonRebinSystErr%d", imult), 9, ptlow);
        pionRebinSystErr[imult] = new TH1D(Form("pionRebinSystErr%d", imult), Form("pionRebinSystErr%d", imult), 9, ptlow);
        phiRebinSystErr[imult] = new TH1D(Form("phiRebinSystErr%d", imult), Form("phiRebinSystErr%d", imult), 9, ptlow);
        chargedRebinSystErr[imult] = new TH1D(Form("chargedRebinSystErr%d", imult), Form("chargedRebinSystErr%d", imult), 9, ptlow);


        ptratiop2phi[imult] = new TH1D(Form("ptratiop2phi%d", imult), Form("ptratiop2phi%d", imult), 9, ptlow);
        ptratiophi2p[imult] = new TH1D(Form("ptratiophi2p%d", imult), Form("ptratiophi2p%d", imult), 9, ptlow);
        ptratiophi2k[imult] = new TH1D(Form("ptratiophi2k%d", imult), Form("ptratiophi2k%d", imult), 9, ptlow);
        ptratiophi2pi[imult] = new TH1D(Form("ptratiophi2pi%d", imult), Form("ptratiophi2pi%d", imult), 9, ptlow);
        ptratiophi2charged[imult] = new TH1D(Form("ptratiophi2charged%d", imult), Form("ptratiophi2charged%d", imult), 9, ptlow);

        pi2charged[imult] = new TH1D(Form("pi2charged%d", imult), Form("pi2charged%d", imult), 9, ptlow);
        k2charged[imult] = new TH1D(Form("k2charged%d", imult), Form("k2charged%d", imult), 9, ptlow);
        proton2charged[imult] = new TH1D(Form("proton2charged%d", imult), Form("proton2charged%d", imult), 9, ptlow);
        proton2pion[imult] = new TH1D(Form("proton2pion%d", imult), Form("proton2pion%d", imult), 9, ptlow);

        pionRebin[imult]->SetLineWidth(2);
        pionRebin[imult]->SetLineColor(kRed+1);
        kaonRebin[imult]->SetLineWidth(2);
        kaonRebin[imult]->SetLineColor(kGreen+2);
        protonRebin[imult]->SetLineWidth(2);
        protonRebin[imult]->SetLineColor(kMagenta+2);

        pi2charged[imult]->SetLineWidth(2);
        pi2charged[imult]->SetLineColor(kRed+1);
        k2charged[imult]->SetLineWidth(2);
        k2charged[imult]->SetLineColor(kGreen+2);
        proton2charged[imult]->SetLineWidth(2);
        proton2charged[imult]->SetLineColor(kMagenta+2);

        for(int i = 1; i <= 9; i++){
            Int_t protonlowbin = protonhist[imult]->GetXaxis()->FindBin(ptlow[i-1] + eps);
            Int_t protonhighbin = protonhist[imult]->GetXaxis()->FindBin(pthigh[i-1] - eps);
            Int_t kaonlowbin = khist[imult]->GetXaxis()->FindBin(ptlow[i-1] + eps);
            Int_t kaonhighbin = khist[imult]->GetXaxis()->FindBin(pthigh[i-1] - eps);
            Int_t pionlowbin = pihist[imult]->GetXaxis()->FindBin(ptlow[i-1] + eps);
            Int_t pionhighbin = pihist[imult]->GetXaxis()->FindBin(pthigh[i-1] - eps);
            Int_t philowbin = phihist[imult]->GetXaxis()->FindBin(ptlow[i-1] + eps);
            Int_t phihighbin = phihist[imult]->GetXaxis()->FindBin(pthigh[i-1] - eps);
          //dealing with re-binning histograms after they've already been re-scaled to just yields
            protonRebin[imult]->SetBinContent(i, protonhist[imult]->Integral(protonlowbin, protonhighbin));
            kaonRebin[imult]->SetBinContent(i, khist[imult]->Integral(kaonlowbin, kaonhighbin));
            pionRebin[imult]->SetBinContent(i, pihist[imult]->Integral(pionlowbin, pionhighbin));
            phiRebin[imult]->SetBinContent(i, phihist[imult]->Integral(philowbin, phihighbin));

          //dealing with re-binning histograms in their original forms
          /* 
            protonRebin[imult]->SetBinContent(i, myIntegral(protonhist[imult], protonlowbin, protonhighbin));
            kaonRebin[imult]->SetBinContent(i, myIntegral(khist[imult], kaonlowbin, kaonhighbin));
            pionRebin[imult]->SetBinContent(i, myIntegral(pihist[imult], pionlowbin, pionhighbin));
            
            protonRebinError[imult]->SetBinContent(i, myIntegralErr(protonhistStat[imult], protonlowbin, protonhighbin));
            kaonRebinError[imult]->SetBinContent(i, myIntegralErr(khistStat[imult], kaonlowbin, kaonhighbin));
            pionRebinError[imult]->SetBinContent(i, myIntegralErr(pihistStat[imult], pionlowbin, pionhighbin));
            protonRebin[imult]->SetBinError(i, protonRebinError[imult]->GetBinContent(i));
            kaonRebin[imult]->SetBinError(i, kaonRebinError[imult]->GetBinContent(i));
            pionRebin[imult]->SetBinError(i, pionRebinError[imult]->GetBinContent(i));

            protonRebinSystErr[imult]->SetBinContent(i, myIntegralErr(protonhistSystUncorr[imult], protonlowbin, protonhighbin));
            kaonRebinSystErr[imult]->SetBinContent(i, myIntegralErr(khistSystUncorr[imult], kaonlowbin, kaonhighbin));
            pionRebinSystErr[imult]->SetBinContent(i, myIntegralErr(pihistSystUncorr[imult], pionlowbin, pionhighbin));

          
            phiRebin[imult]->SetBinContent(i, phihist[imult]->Integral(philowbin, phihighbin, "width"));
            
            phiRebinError[imult]->SetBinContent(i, mySimpleErrInt(phihistStat[imult], philowbin, phihighbin));
            phiRebinSystErr[imult]->SetBinContent(i, mySimpleErrInt(phihistSyst[imult], philowbin, phihighbin));
                
            phiRebin[imult]->SetBinError(i, phiRebinError[imult]->GetBinContent(i));
          */

            chargedRebin[imult]->SetBinContent(i, protonRebin[imult]->GetBinContent(i) + kaonRebin[imult]->GetBinContent(i) + pionRebin[imult]->GetBinContent(i));
            //chargedRebinError[imult]->SetBinContent(i, TMath::Sqrt(TMath::Power(protonRebinError[imult]->GetBinContent(i), 2.0) + TMath::Power(kaonRebinError[imult]->GetBinContent(i), 2.0) + TMath::Power(pionRebinError[imult]->GetBinContent(i), 2.0)));
            //chargedRebinSystErr[imult]->SetBinContent(i, TMath::Sqrt(TMath::Power(protonRebinSystErr[imult]->GetBinContent(i), 2.0) + TMath::Power(kaonRebinSystErr[imult]->GetBinContent(i), 2.0) + TMath::Power(pionRebinSystErr[imult]->GetBinContent(i), 2.0)));

            //chargedRebin[imult]->SetBinError(i, TMath::Sqrt(TMath::Power(protonRebinError[imult]->GetBinContent(i), 2.0) + TMath::Power(kaonRebinError[imult]->GetBinContent(i), 2.0) + TMath::Power(pionRebinError[imult]->GetBinContent(i), 2.0)));

            ptratiop2phi[imult]->SetBinContent(i, protonRebin[imult]->GetBinContent(i)/(phiRebin[imult]->GetBinContent(i)));
            ptratiophi2p[imult]->SetBinContent(i, phiRebin[imult]->GetBinContent(i)/(protonRebin[imult]->GetBinContent(i)));
            ptratiophi2k[imult]->SetBinContent(i, phiRebin[imult]->GetBinContent(i)/(kaonRebin[imult]->GetBinContent(i)));
            ptratiophi2pi[imult]->SetBinContent(i, phiRebin[imult]->GetBinContent(i)/(pionRebin[imult]->GetBinContent(i)));
            ptratiophi2charged[imult]->SetBinContent(i, phiRebin[imult]->GetBinContent(i)/chargedRebin[imult]->GetBinContent(i));

            pi2charged[imult]->SetBinContent(i, pionRebin[imult]->GetBinContent(i)/chargedRebin[imult]->GetBinContent(i));
            k2charged[imult]->SetBinContent(i, kaonRebin[imult]->GetBinContent(i)/chargedRebin[imult]->GetBinContent(i));
            proton2charged[imult]->SetBinContent(i, protonRebin[imult]->GetBinContent(i)/chargedRebin[imult]->GetBinContent(i));
            proton2pion[imult]->SetBinContent(i, protonRebin[imult]->GetBinContent(i)/pionRebin[imult]->GetBinContent(i));

        }
        printf("past setup of all rebin for imult %d\n", imult);

        pions[6-imult] = pionRebin[imult]->Integral(pionRebin[imult]->GetXaxis()->FindBin(lowpt+eps), pionRebin[imult]->GetXaxis()->FindBin(highpt-eps));
        kaons[6-imult] = kaonRebin[imult]->Integral(kaonRebin[imult]->GetXaxis()->FindBin(lowpt+eps), kaonRebin[imult]->GetXaxis()->FindBin(highpt-eps));
        protons[6-imult] = protonRebin[imult]->Integral(protonRebin[imult]->GetXaxis()->FindBin(lowpt+eps), protonRebin[imult]->GetXaxis()->FindBin(highpt-eps));
        phi[6-imult] = phiRebin[imult]->Integral(phiRebin[imult]->GetXaxis()->FindBin(lowpt+eps), phiRebin[imult]->GetXaxis()->FindBin(highpt-eps));
        //phiErr[6-imult] = addInQuad(phiRebinError[imult], phiRebinError[imult]->GetXaxis()->FindBin(lowpt+eps), phiRebinError[imult]->GetXaxis()->FindBin(highpt-eps));
        //phiSystErr[6-imult] = addInQuad(phiRebinSystErr[imult], phiRebinSystErr[imult]->GetXaxis()->FindBin(lowpt+eps), phiRebinSystErr[imult]->GetXaxis()->FindBin(highpt-eps));
        charged[6-imult] = chargedRebin[imult]->Integral(chargedRebin[imult]->GetXaxis()->FindBin(lowpt+eps), chargedRebin[imult]->GetXaxis()->FindBin(highpt-eps));
        //chargedErr[6-imult] = addInQuad(chargedRebinError[imult], chargedRebinError[imult]->GetXaxis()->FindBin(lowpt+eps), chargedRebinError[imult]->GetXaxis()->FindBin(highpt-eps));
        //chargedSystErr[6-imult] = addInQuad(chargedRebinSystErr[imult], chargedRebinSystErr[imult]->GetXaxis()->FindBin(lowpt+eps), chargedRebinSystErr[imult]->GetXaxis()->FindBin(highpt-eps));
        ratio[6-imult] = phi[6-imult]/charged[6-imult];
        //ratioErr[6-imult] = quotientErrors(phi[6-imult], phiErr[6-imult], charged[6-imult], chargedErr[6-imult]);
        //ratioSystErr[6-imult] = quotientErrors(phi[6-imult], phiSystErr[6-imult], charged[6-imult], chargedSystErr[6-imult]);
        ratio2proton[6-imult] = phi[6-imult]/protons[6-imult];
        ratio2k[6-imult] = phi[6-imult]/kaons[6-imult];
        ratio2pi[6-imult] = phi[6-imult]/pions[6-imult];


        totpions[6-imult] = pionRebin[imult]->Integral(1, 10);
        totkaons[6-imult] = kaonRebin[imult]->Integral(1, 10);
        totprotons[6-imult] = protonRebin[imult]->Integral(1, 10);
        totphi[6-imult] = phiRebin[imult]->Integral(1, 10);
        totcharged[6-imult] = chargedRebin[imult]->Integral(1, 10);
        totratio[6-imult] = totphi[6-imult]/totcharged[6-imult];
        totratio2pi[6-imult] = totphi[6-imult]/totpions[6-imult];
        totratio2k[6-imult] = totphi[6-imult]/totkaons[6-imult];
        totratio2proton[6-imult] = totphi[6-imult]/totprotons[6-imult];


        printf("mult: %d, ratio2charged: %f\n", imult, totratio[6-imult]);
    }

    for(int i = 1; i <=10; i++){
        chargedRebinMB->SetBinContent(i, 0.05*chargedRebin[0]->GetBinContent(i) + 0.05*chargedRebin[1]->GetBinContent(i) + 0.1*chargedRebin[2]->GetBinContent(i) + 0.2*chargedRebin[3]->GetBinContent(i) + 0.2*chargedRebin[4]->GetBinContent(i) + 0.2*chargedRebin[5]->GetBinContent(i) + 0.2*chargedRebin[6]->GetBinContent(i));
        chargedRebinMB->SetBinContent(i, chargedRebinMB->GetBinContent(i)/chargedRebinMB->GetXaxis()->GetBinWidth(i));
    }
    Float_t ratio6[6];
    Float_t ratioErr6[6];
    Float_t ratioSystErr6[6];
    for(int i = 0; i < 6; i++){
        ratio6[i] = ratio[i+1];
        //ratioErr6[i] = ratioErr[i+1];
        //ratioSystErr6[i] = ratioSystErr[i+1];
    }

    TCanvas* cchpt= new TCanvas("cchpt", "cchpt", 50, 50, 600, 600);
    chargedRebin[0]->Draw();

    TGraph* ratiograph = new TGraph(6, revmultbins6, ratio6);
    //TGraphErrors* ratiographEr = new TGraphErrors(6, revmultbins6, ratio6, revmultwidths6, ratioErr6);
    //TGraphErrors* ratiographSyst = new TGraphErrors(6, revmultbins6, ratio6, revmultwidths6, ratioSystErr6);
    TGraph* totratiograph = new TGraph(7, revmultbins, totratio);
    TGraph* ratio2protongraph = new TGraph(7, revmultbins, ratio2proton);
    ratio2protongraph->SetMarkerStyle(23);
    ratio2protongraph->SetMarkerSize(2);
    ratio2protongraph->SetMarkerColor(kBlue);
    ratio2protongraph->SetLineColor(kBlue);
    
    TGraph* ratio2kgraph = new TGraph(7, revmultbins, ratio2k);
    ratio2kgraph->SetMarkerStyle(23);
    ratio2kgraph->SetMarkerSize(2);
    ratio2kgraph->SetMarkerColor(kRed);
    ratio2kgraph->SetLineColor(kRed);
    
    TGraph* ratio2pigraph = new TGraph(7, revmultbins, ratio2pi);
    ratio2pigraph->SetMarkerStyle(23);
    ratio2pigraph->SetMarkerSize(2);
    ratio2pigraph->SetMarkerColor(kGreen+2);
    ratio2pigraph->SetLineColor(kGreen+2);

    //setup histogram for axis drawing
    TH1D* ratiohist = new TH1D("ratiohist", "", 7, histbins);
    for(int i =1; i<=6; i++){
        ratiohist->SetBinContent(i+1, ratio6[i-1]);
    }
    
    //total ratio graphs
    TGraph* totratio2pigraph = new TGraph(7, revmultbins, totratio2pi);
    totratio2pigraph->SetMarkerStyle(23);
    totratio2pigraph->SetMarkerSize(2);
    totratio2pigraph->SetMarkerColor(kGray);
    totratio2pigraph->SetLineColor(kGray);
    TGraph* totratio2kgraph = new TGraph(7, revmultbins, totratio2k);
    totratio2kgraph->SetMarkerStyle(23);
    totratio2kgraph->SetMarkerSize(2);
    totratio2kgraph->SetMarkerColor(kRed);
    totratio2kgraph->SetLineColor(kRed);
    TGraph* totratio2protongraph = new TGraph(7, revmultbins, totratio2proton);
    totratio2protongraph->SetMarkerStyle(23);
    totratio2protongraph->SetMarkerSize(2);
    totratio2protongraph->SetMarkerColor(kBlue);
    totratio2protongraph->SetLineColor(kBlue);

    //individual particle yield graphs
    TGraph* piongraph = new TGraph(7, revmultbins, pions);
    piongraph->SetMarkerStyle(23);
    piongraph->SetMarkerSize(2);
    piongraph->SetMarkerColor(kBlue+1);
    piongraph->SetLineColor(kBlue+1);
    TGraph* kaongraph = new TGraph(7, revmultbins, kaons);
    kaongraph->SetMarkerStyle(23);
    kaongraph->SetMarkerSize(2);
    kaongraph->SetMarkerColor(kGreen+1);
    kaongraph->SetLineColor(kGreen+1);
    TGraph* protongraph = new TGraph(7, revmultbins, protons);
    protongraph->SetMarkerStyle(23);
    protongraph->SetMarkerSize(2);
    protongraph->SetMarkerColor(kRed+1);
    protongraph->SetLineColor(kRed+1);
    TGraph* chargedgraph = new TGraph(7, revmultbins, charged);
    chargedgraph->SetMarkerStyle(23);
    chargedgraph->SetMarkerSize(2);
    chargedgraph->SetMarkerColor(kBlack);
    chargedgraph->SetLineColor(kBlack);
    TGraph* phigraph = new TGraph(7, revmultbins, phi);
    phigraph->SetMarkerStyle(23);
    phigraph->SetMarkerSize(2);
    phigraph->SetMarkerColor(kViolet);
    phigraph->SetLineColor(kViolet);
  

    TCanvas* cyields = new TCanvas("cyields", "cyields", 50, 50, 600, 600);
    cyields->cd();
    //protonRebin[0]->Draw();
    //pionRebin[0]->Draw("SAME");
    //chargedgraph->Draw("ALP");
    //piongraph->Draw("LP");
    //kaongraph->Draw("LP");
    //protongraph->Draw("LP");
    //phigraph->Draw("LP");
    
    //ratiograph->Draw("ALP");
    //totratio2protongraph->Draw("ALP");
    //totratio2kgraph->Draw("LP");
    totratio2pigraph->Draw("ALP");

    //Check on (p+p)/phi pt dependent ratio from the paper

    TH1D* paperratio[7];

    TFile* ratioFilemult1 = new TFile("~/phiStudies/pubRatioCheck/HEPData-ins1418181-v1-Table_22.root");
    TDirectory* ratiodirmult1 = (TDirectory*) ratioFilemult1->GetDirectory("Table 22");
    paperratio[0] = (TH1D*)(ratiodirmult1->Get("Hist1D_y1")->Clone("paperratio1"));
    TFile* ratioFilemult2 = new TFile("~/phiStudies/pubRatioCheck/HEPData-ins1418181-v1-Table_23.root");
    TDirectory* ratiodirmult2 = (TDirectory*) ratioFilemult2->GetDirectory("Table 23");
    paperratio[1] = (TH1D*)(ratiodirmult2->Get("Hist1D_y1")->Clone("paperratio2"));
    TFile* ratioFilemult3 = new TFile("~/phiStudies/pubRatioCheck/HEPData-ins1418181-v1-Table_24.root");
    TDirectory* ratiodirmult3 = (TDirectory*) ratioFilemult3->GetDirectory("Table 24");
    paperratio[2] = (TH1D*)(ratiodirmult3->Get("Hist1D_y1")->Clone("paperratio3"));
    TFile* ratioFilemult4 = new TFile("~/phiStudies/pubRatioCheck/HEPData-ins1418181-v1-Table_25.root");
    TDirectory* ratiodirmult4 = (TDirectory*) ratioFilemult4->GetDirectory("Table 25");
    paperratio[3] = (TH1D*)(ratiodirmult4->Get("Hist1D_y1")->Clone("paperratio4"));
    TFile* ratioFilemult5 = new TFile("~/phiStudies/pubRatioCheck/HEPData-ins1418181-v1-Table_26.root");
    TDirectory* ratiodirmult5 = (TDirectory*) ratioFilemult5->GetDirectory("Table 26");
    paperratio[4] = (TH1D*)(ratiodirmult5->Get("Hist1D_y1")->Clone("paperratio5"));
    TFile* ratioFilemult6 = new TFile("~/phiStudies/pubRatioCheck/HEPData-ins1418181-v1-Table_27.root");
    TDirectory* ratiodirmult6 = (TDirectory*) ratioFilemult6->GetDirectory("Table 27");
    paperratio[5] = (TH1D*)(ratiodirmult6->Get("Hist1D_y1")->Clone("paperratio6"));
    TFile* ratioFilemult7 = new TFile("~/phiStudies/pubRatioCheck/HEPData-ins1418181-v1-Table_28.root");
    TDirectory* ratiodirmult7 = (TDirectory*) ratioFilemult7->GetDirectory("Table 28");
    paperratio[6] = (TH1D*)(ratiodirmult7->Get("Hist1D_y1")->Clone("paperratio7"));

    
    TCanvas* cpi = new TCanvas("cpi", "cpi", 50, 50, 600, 600);
    cpi->cd();
    Int_t colors[] = {kRed+2, kOrange+9, kOrange+7, kOrange, kSpring+9, kGreen+2, kAzure-3};
    for(int i = 0; i < 7; i++){
        ptratiophi2pi[i]->SetLineColor(colors[i]);
        ptratiophi2pi[i]->SetLineWidth(2);
        ptratiophi2pi[i]->Draw("SAME");
    }

    TCanvas* ck = new TCanvas("ck", "ck", 50, 50, 600, 600);
    ck->cd();
    for(int i = 0; i < 7; i++){
        ptratiophi2k[i]->SetLineColor(colors[i]);
        ptratiophi2k[i]->SetLineWidth(2);
        ptratiophi2k[i]->Draw("SAME");
    }

    TCanvas* cp = new TCanvas("cp", "cp", 50, 50, 600, 600);
    cp->cd();
    for(int i = 0; i < 7; i++){
        ptratiophi2p[i]->SetLineColor(colors[i]);
        ptratiophi2p[i]->SetLineWidth(2);
        ptratiophi2p[i]->Draw("SAME");
    }

    TCanvas* cp2[7];   
    for(int i = 0; i < 7; i++){
        cp2[i] = new TCanvas(Form("cp2_%d", i), Form("cp2_%d", i), 50, 50, 600, 600);
        cp2[i]->cd();
        ptratiop2phi[i]->SetLineColor(colors[i]);
        ptratiop2phi[i]->SetLineWidth(2);
        ptratiop2phi[i]->Draw("SAME");
        paperratio[i]->Draw("SAME");
    }


    //draw ratio of ratios
    TH1D* paperratioratio[7];

    for(int i = 0; i<7; i++){
        paperratioratio[i] = (TH1D*)ptratiop2phi[i]->Clone(Form("paperratioratio%d",i));
        //paperratioratio[i]->Divide(paperratio[i]);
    }


    TPaveText* ratiotext = new TPaveText(0.4, 0.6, 0.5, 0.7, "NDC");
    ratiotext->SetFillColor(kWhite);
    ratiotext->SetBorderSize(0);
    ratiotext->AddText("#phi/h Ratio from published spectra");
    ratiotext->AddText("2.0 < p_{T} < 4.0 GeV/c");

    TCanvas* ccheck = new TCanvas("ccheck", "ccheck", 50, 50, 600, 600);
    ccheck->cd();
    ptratiop2phi[0]->Draw();
    ptratiop2phi[6]->Draw("SAME"); 

    TCanvas* cratio = new TCanvas("cratio", "cratio", 50, 50, 600, 600);
    cratio->cd();

    gStyle->SetErrorX(0.5);
    ratiohist->GetYaxis()->SetRangeUser(0.035, 0.05);
    ratiohist->GetYaxis()->SetMaxDigits(3);
    ratiohist->GetYaxis()->SetTitle("#phi/h");
    ratiohist->SetStats(kFALSE);
    ratiohist->Draw("AXIS");

    ratiohist->GetXaxis()->SetLabelOffset(999);
    //ratioNearHist->GetXaxis()->SetTitleOffset(999);
    ratiohist->GetXaxis()->SetTickSize(0.0);

    //ratiosNear->Draw("P");
    gPad->Update();
    TGaxis* newaxis = new TGaxis(gPad->GetUxmax(),
            gPad->GetUymin(),
            gPad->GetUxmin(),
            gPad->GetUymin(),
            ratiohist->GetXaxis()->GetXmin(),
            ratiohist->GetXaxis()->GetXmax(),
            510,"-");
    newaxis->SetLabelOffset(-0.03);
    newaxis->SetTitle("Multipliciy % (VOA)");
    //newaxis->SetTitleOffset(1.3);
    newaxis->Draw(); 

    //ratiographEr->SetMarkerStyle(23);
    //ratiographEr->SetMarkerSize(2);
    //ratiographSyst->Draw("5");
    //ratiographEr->Draw("LP");
    ratiotext->Draw();
    ratio2protongraph->Draw("LP");
    ratio2kgraph->Draw("LP");
    ratio2pigraph->Draw("LP");
    //ptratiophi2k[6]->Draw("SAME");
    //ptratiophi2pi[6]->Draw("SAME");
    
    //protonhist[6]->Draw();
    //protonRebin[6]->SetLineColor(kRed+1);
    //protonRebin[6]->Draw("SAME");
    

    TCanvas* cratio2charge = new TCanvas("cratio2charge", "cratio2charge", 50, 50, 600, 600);
    cratio2charge->cd();
    pi2charged[0]->Draw("");
    k2charged[0]->Draw("SAME");
    proton2charged[0]->Draw("SAME");
    //ptratiophi2charged[0]->Draw("SAME");
    /*for(int i = 0; i<6; i++){
        proton2pion[i]->SetLineColor(colors[i]);
        proton2pion[i]->SetLineWidth(2);
        proton2pion[i]->Draw("SAME");
    }*/

    TCanvas* cphi2charge = new TCanvas("cphi2charge", "cphi2charge", 50, 50, 600, 600);
    cphi2charge->cd();
    for(int i=0; i< 6; i++){
        ptratiophi2charged[i]->SetLineColor(colors[i]);
        ptratiophi2charged[i]->Draw("SAME");
    }

    TCanvas* cphi2pi = new TCanvas("cphi2pi", "cphi2pi", 50, 50, 600, 600);
    cphi2pi->cd();
    for(int i=0; i< 6; i++){
        ptratiophi2pi[i]->SetLineColor(colors[i]);
        ptratiophi2pi[i]->Draw("SAME");
    }
    
    TCanvas* cphipt = new TCanvas("cphipt", "cphipt", 50, 50, 600, 600);
    cphipt->cd();
    phihist[0]->Draw("HIST");
    phihist[5]->Draw("HIST SAME");
    phiMB->SetLineColor(kRed);
    phiMB->Draw("HIST SAME");
    chargeMB->SetLineColor(kGreen+1);
    chargeMB->Draw("HIST P SAME");
    chargedRebinMB->SetLineColor(kMagenta+1);
    chargedRebinMB->Draw("HIST P SAME");

    Double_t phiMBInt = phiMB->Integral(phiMB->GetXaxis()->FindBin(2.0 + 0.0001), phiMB->GetXaxis()->FindBin(4.0 - 0.0001), "width");
    Double_t chargeMBInt = chargeMB->Integral(chargeMB->GetXaxis()->FindBin(2.0 + 0.0001), chargeMB->GetXaxis()->FindBin(4.0 - 0.0001), "width");
    Double_t chargedRebinMBInt = chargedRebinMB->Integral(chargedRebinMB->GetXaxis()->FindBin(2.0 + 0.0001), chargedRebinMB->GetXaxis()->FindBin(4.0 - 0.0001), "width");


    printf("phi over MB paper ratio: %f\n", phiMBInt/chargeMBInt);
    //printf("charged MB ratio: %f\n", chargedRebinMBInt/chargeMBInt);
    printf("phi over MB calc ratio: %f\n", phiMBInt/chargedRebinMBInt);


    TFile* file = new TFile("~/phiStudies/results_onlineEff/FAST/hphi_0_20_rapiditytest.root");
    TList* list = (TList*) file->Get("phiCorr_mult_0_20_");
    THnSparseF* trigDist = (THnSparseF*)list->FindObject("fTrigDist");
    trigDist->GetAxis(0)->SetRange(0, 0);
    trigDist->GetAxis(2)->SetRangeUser(-0.8 + 0.0001, 0.8 - 0.00001);
    TH1D* mychargedpt = (TH1D*)trigDist->Projection(0);

    /*
    for(int ibin = 1; ibin < mychargept->GetXaxis()->GetNbins(); ibin++){
        mychargept->SetBinContent(ibin, mychargept->GetBinContent(ibin)/(2.0*TMath::Pi()*mychargept->GetXaxis()->GetBinWidth(ibin)8mychargept->GetXaxis()->GetBinCenter(ibin)));
    }
    */
    TH1D* charged020 = (TH1D*)chargedRebin[0]->Clone("charged020");
    for(int ibin = 1; ibin < charged020->GetXaxis()->GetNbins(); ibin++){
        
    }

}
