#include "TLegend.h"
#include <string>
#include <sstream>

void plot_phi_histo(string inputName){
    gROOT->Reset();
    gROOT->SetStyle("Plain");
    gStyle->SetPalette(1);
    gStyle->SetCanvasColor(kWhite);
    gStyle->SetCanvasBorderMode(0);
    gStyle->SetPadBorderMode(0);
    gStyle->SetTitleBorderSize(0);
    gStyle->SetOptStat(0);
    gStyle->SetErrorX(0);
    gStyle->SetTitleW(0.9);
    gStyle->SetTitleSize(0.05, "xyz");
    gStyle->SetTitleSize(0.06, "h");

    TFile *histoFile = new TFile(inputName.c_str());
    histoFile->cd("PhiReconstruction");
    THnSparseF *phiInvMass = (THnSparseF *)InvMass->FindObject("fPhiInvMass");
    THnSparseF *likeSignInvMass = (THnSparseF *)InvMass->FindObject("fPhiLikeSignInvMass");
    TH1D *corrInvMass[16];

    Double_t sideband = 0.0;
    Double_t likeSignSideBand = 0.0;
    Double_t scaleFactor = 0.0;

    if(phiInvMass && likeSignInvMass){
        phiInvMass->GetAxis(1)->SetRange(320,400);
        phiInvMass->GetAxis(1)->SetTitle("Inv Mass (GeV/c^2)");
        likeSignInvMass->GetAxis(1)->SetRange(320,400);
        likeSignInvMass->SetTitle("Inv Mass (GeV/c^2)");

        phiInvMass->GetAxis(0)->SetRange(4,6);
        TH1D *phiInvMass_04_06 = (TH1D*)phiInvMass->Projection(1);
        phiInvMass_04_06->SetTitle("0.4 < p_{T} < 0.6 GeV/c");
        sideband = phiInvMass_04_06->Integral(40, 53);
        likeSignInvMass->GetAxis(0)->SetRange(4,6);
        TH1D *likeSignInvMass_04_06 = (TH1D*)likeSignInvMass->Projection(1);
        likeSignInvMass_04_06->SetTitle("0.4 < p_{T} < 0.6 GeV/c");
        likeSignSideBand = likeSignInvMass_04_06->Integral(40, 53);
//      scaleFactor = sideband/likeSignSideBand;
//      printf("sideband: %d, likesignsideband: %d, scalefactor: %d", sideband, likeSignSideBand, scaleFactor);
        likeSignInvMass_04_06->Scale(sideband/likeSignSideBand);
        likeSignInvMass_04_06->SetLineColor(2);
        corrInvMass[0] = (TH1D*)phiInvMass_04_06->Clone("corrInvMass_04_06");
        corrInvMass[0]->Add(likeSignInvMass_04_06, -1);

        phiInvMass->GetAxis(0)->SetRange(6,8);
        TH1D *phiInvMass_06_08 = (TH1D*)phiInvMass->Projection(1);
        phiInvMass_06_08->SetTitle("0.6 < p_{T} < 0.8 GeV/c");
        sideband = phiInvMass_06_08->Integral(40,53);
        likeSignInvMass->GetAxis(0)->SetRange(6,8);
        TH1D *likeSignInvMass_06_08 = (TH1D*)likeSignInvMass->Projection(1);
        likeSignInvMass_06_08->SetTitle("0.6 < p_{T} < 0.8 GeV/c");
        likeSignSideBand = likeSignInvMass_06_08->Integral(40,53);
        likeSignInvMass_06_08->Scale(sideband/likeSignSideBand);
        likeSignInvMass_06_08->SetLineColor(2);
        corrInvMass[1] = (TH1D*)phiInvMass_06_08->Clone("corrInvMass_06_08");
        corrInvMass[1]->Add(likeSignInvMass_06_08, -1);


        phiInvMass->GetAxis(0)->SetRange(8,10);
        TH1D *phiInvMass_08_10 = (TH1D*)phiInvMass->Projection(1);
        phiInvMass_08_10->SetTitle("0.8 < p_{T} < 1 GeV/c");
        sideband = phiInvMass_08_10->Integral(40,53);
        likeSignInvMass->GetAxis(0)->SetRange(8,10);
        TH1D *likeSignInvMass_08_10 = (TH1D*)likeSignInvMass->Projection(1);
        likeSignInvMass_08_10->SetTitle("0.8 < p_{T} < 1 GeV/c");
        likeSignSideBand = likeSignInvMass_08_10->Integral(40,53);
        likeSignInvMass_08_10->Scale(sideband/likeSignSideBand);
        likeSignInvMass_08_10->SetLineColor(2);
        corrInvMass[2] = (TH1D*)phiInvMass_08_10->Clone("corrInvMass_08_10");
        corrInvMass[2]->Add(likeSignInvMass_08_10, -1);


        phiInvMass->GetAxis(0)->SetRange(10,12);
        TH1D *phiInvMass_10_12 = (TH1D*)phiInvMass->Projection(1);
        phiInvMass_10_12->SetTitle("1 < p_{T} < 1.2 GeV/c");
        sideband = phiInvMass_10_12->Integral(40,53);
        likeSignInvMass->GetAxis(0)->SetRange(10,12);
        TH1D *likeSignInvMass_10_12 = (TH1D*)likeSignInvMass->Projection(1);
        likeSignInvMass_10_12->SetTitle("1 < p_{T} < 1.2 GeV/c");
        likeSignSideBand = likeSignInvMass_10_12->Integral(40,53);
        likeSignInvMass_10_12->Scale(sideband/likeSignSideBand);
        likeSignInvMass_10_12->SetLineColor(2);
        corrInvMass[3] = (TH1D*)phiInvMass_10_12->Clone("corrInvMass_10_12");
        corrInvMass[3]->Add(likeSignInvMass_10_12, -1);


        phiInvMass->GetAxis(0)->SetRange(12,14);
        TH1D *phiInvMass_12_14 = (TH1D*)phiInvMass->Projection(1);
        phiInvMass_12_14->SetTitle("1.2 < p_{T} < 1.4 GeV/c");
        sideband = phiInvMass_12_14->Integral(40,53);
        likeSignInvMass->GetAxis(0)->SetRange(12,14);
        TH1D *likeSignInvMass_12_14 = (TH1D*)likeSignInvMass->Projection(1);
        likeSignInvMass_12_14->SetTitle("1.2 < p_{T} < 1.4 GeV/c");
        likeSignSideBand = likeSignInvMass_12_14->Integral(40,53);
        likeSignInvMass_12_14->Scale(sideband/likeSignSideBand);
        likeSignInvMass_12_14->SetLineColor(2);
        corrInvMass[4] = (TH1D*)phiInvMass_12_14->Clone("corrInvMass_12_14");
        corrInvMass[4]->Add(likeSignInvMass_12_14, -1);


        phiInvMass->GetAxis(0)->SetRange(14,16);
        TH1D *phiInvMass_14_16 = (TH1D*)phiInvMass->Projection(1);
        phiInvMass_14_16->SetTitle("1.4 < p_{T} < 1.6 GeV/c");
        sideband = phiInvMass_14_16->Integral(40,53);
        likeSignInvMass->GetAxis(0)->SetRange(14,16);
        TH1D *likeSignInvMass_14_16 = (TH1D*)likeSignInvMass->Projection(1);
        likeSignInvMass_14_16->SetTitle("1.4 < p_{T} < 1.6 GeV/c");
        likeSignSideBand = likeSignInvMass_14_16->Integral(40,53);
        likeSignInvMass_14_16->Scale(sideband/likeSignSideBand);
        likeSignInvMass_14_16->SetLineColor(2);
        corrInvMass[5] = (TH1D*)phiInvMass_14_16->Clone("corrInvMass_14_16");
        corrInvMass[5]->Add(likeSignInvMass_14_16, -1);


        phiInvMass->GetAxis(0)->SetRange(16,18);
        TH1D *phiInvMass_16_18 = (TH1D*)phiInvMass->Projection(1);
        phiInvMass_16_18->SetTitle("1.6 < p_{T} < 1.8 GeV/c");
        sideband = phiInvMass_16_18->Integral(40,53);
        likeSignInvMass->GetAxis(0)->SetRange(16,18);
        TH1D *likeSignInvMass_16_18 = (TH1D*)likeSignInvMass->Projection(1);
        likeSignInvMass_16_18->SetTitle("1.6 < p_{T} < 1.8 GeV/c");
        likeSignSideBand = likeSignInvMass_16_18->Integral(40,53);
        likeSignInvMass_16_18->Scale(sideband/likeSignSideBand);
        likeSignInvMass_16_18->SetLineColor(2);
        corrInvMass[6] = (TH1D*)phiInvMass_16_18->Clone("corrInvMass_16_18");
        corrInvMass[6]->Add(likeSignInvMass_16_18, -1);


        phiInvMass->GetAxis(0)->SetRange(18,20);
        TH1D *phiInvMass_18_20 = (TH1D*)phiInvMass->Projection(1);
        phiInvMass_18_20->SetTitle("1.8 < p_{T} < 2 GeV/c");
        sideband = phiInvMass_18_20->Integral(40,53);
        likeSignInvMass->GetAxis(0)->SetRange(18,20);
        TH1D *likeSignInvMass_18_20 = (TH1D*)likeSignInvMass->Projection(1);
        likeSignInvMass_18_20->SetTitle("1.8 < p_{T} < 2 GeV/c");
        likeSignSideBand = likeSignInvMass_18_20->Integral(40,53);
        likeSignInvMass_18_20->Scale(sideband/likeSignSideBand);
        likeSignInvMass_18_20->SetLineColor(2);
        corrInvMass[7] = (TH1D*)phiInvMass_18_20->Clone("corrInvMass_18_20");
        corrInvMass[7]->Add(likeSignInvMass_18_20, -1);


        phiInvMass->GetAxis(0)->SetRange(20,25);
        TH1D *phiInvMass_20_25 = (TH1D*)phiInvMass->Projection(1);
        phiInvMass_20_25->SetTitle("2 < p_{T} < 2.5 GeV/c");
        sideband = phiInvMass_20_25->Integral(40,53);
        likeSignInvMass->GetAxis(0)->SetRange(20,25);
        TH1D *likeSignInvMass_20_25 = (TH1D*)likeSignInvMass->Projection(1);
        likeSignInvMass_20_25->SetTitle("2 < p_{T} < 2.5 GeV/c");
        likeSignSideBand = likeSignInvMass_20_25->Integral(40,53);
        likeSignInvMass_20_25->Scale(sideband/likeSignSideBand);
        likeSignInvMass_20_25->SetLineColor(2);


        phiInvMass->GetAxis(0)->SetRange(25,30);
        TH1D *phiInvMass_25_30 = (TH1D*)phiInvMass->Projection(1);
        phiInvMass_25_30->SetTitle("2.5 < p_{T} < 3 GeV/c");
        sideband = phiInvMass_25_30->Integral(40,53);
        likeSignInvMass->GetAxis(0)->SetRange(25,30);
        TH1D *likeSignInvMass_25_30 = (TH1D*)likeSignInvMass->Projection(1);
        likeSignInvMass_25_30->SetTitle("2.5 < p_{T} < 3 GeV/c");
        likeSignSideBand = likeSignInvMass_25_30->Integral(40,53);
        likeSignInvMass_25_30->Scale(sideband/likeSignSideBand);
        likeSignInvMass_25_30->SetLineColor(2);


        phiInvMass->GetAxis(0)->SetRange(30,40);
        TH1D *phiInvMass_30_40 = (TH1D*)phiInvMass->Projection(1);
        phiInvMass_30_40->SetTitle("3 < p_{T} < 4 GeV/c");
        sideband = phiInvMass_30_40->Integral(40,53);
        likeSignInvMass->GetAxis(0)->SetRange(30,40);
        TH1D *likeSignInvMass_30_40 = (TH1D*)likeSignInvMass->Projection(1);
        likeSignInvMass_30_40->SetTitle("3 < p_{T} < 4 GeV/c");
        likeSignSideBand = likeSignInvMass_30_40->Integral(40,53);
        likeSignInvMass_30_40->Scale(sideband/likeSignSideBand);
        likeSignInvMass_30_40->SetLineColor(2);


        phiInvMass->GetAxis(0)->SetRange(40,50);
        TH1D *phiInvMass_40_50 = (TH1D*)phiInvMass->Projection(1);
        phiInvMass_40_50->SetTitle("4 < p_{T} < 5 GeV/c");
        sideband = phiInvMass_40_50->Integral(40,53);
        likeSignInvMass->GetAxis(0)->SetRange(40,50);
        TH1D *likeSignInvMass_40_50 = (TH1D*)likeSignInvMass->Projection(1);
        likeSignInvMass_40_50->SetTitle("4 < p_{T} < 5 GeV/c");
        likeSignSideBand = likeSignInvMass_40_50->Integral(40,53);
        likeSignInvMass_40_50->Scale(sideband/likeSignSideBand);
        likeSignInvMass_40_50->SetLineColor(2);


        phiInvMass->GetAxis(0)->SetRange(50,60);
        TH1D *phiInvMass_50_60 = (TH1D*)phiInvMass->Projection(1);
        phiInvMass_50_60->SetTitle("5 < p_{T} < 6 GeV/c");
        sideband = phiInvMass_50_60->Integral(40,53);
        likeSignInvMass->GetAxis(0)->SetRange(50,60);
        TH1D *likeSignInvMass_50_60 = (TH1D*)likeSignInvMass->Projection(1);
        likeSignInvMass_50_60->SetTitle("5 < p_{T} < 6 GeV/c");
        likeSignSideBand = likeSignInvMass_50_60->Integral(40,53);
        likeSignInvMass_50_60->Scale(sideband/likeSignSideBand);
        likeSignInvMass_50_60->SetLineColor(2);


        phiInvMass->GetAxis(0)->SetRange(60,70);
        TH1D *phiInvMass_60_70 = (TH1D*)phiInvMass->Projection(1);
        phiInvMass_60_70->SetTitle("6 < p_{T} < 7 GeV/c");
        sideband = phiInvMass_60_70->Integral(40,53);
        likeSignInvMass->GetAxis(0)->SetRange(60,70);
        TH1D *likeSignInvMass_60_70 = (TH1D*)likeSignInvMass->Projection(1);
        likeSignInvMass_60_70->SetTitle("6 < p_{T} < 7 GeV/c");
        likeSignSideBand = likeSignInvMass_60_70->Integral(40,53);
        likeSignInvMass_60_70->Scale(sideband/likeSignSideBand);
        likeSignInvMass_60_70->SetLineColor(2);


        phiInvMass->GetAxis(0)->SetRange(70,100);
        TH1D *phiInvMass_70_100 = (TH1D*)phiInvMass->Projection(1);
        phiInvMass_70_100->SetTitle("7 < p_{T} < 10 GeV/c");
        sideband = phiInvMass_70_100->Integral(40,53);
        likeSignInvMass->GetAxis(0)->SetRange(70,100);
        TH1D *likeSignInvMass_70_100 = (TH1D*)likeSignInvMass->Projection(1);
        likeSignInvMass_70_100->SetTitle("7 < p_{T} < 10 GeV/c");
        likeSignSideBand = likeSignInvMass_70_100->Integral(40,53);
        likeSignInvMass_70_100->Scale(sideband/likeSignSideBand);
        likeSignInvMass_70_100->SetLineColor(2);

    }   


    TCanvas* c0 = new TCanvas("cInvMassPtBins", "cInvMassPtBins", 50,50, 1000, 1000);
    c0->Divide(4,4);
    c0->cd(1);
    phiInvMass_04_06->Draw("PE");
    likeSignInvMass_04_06->Draw("SAME");
    c0->cd(2);
    phiInvMass_06_08->Draw("E");
    likeSignInvMass_06_08->Draw("SAME");
    c0->cd(3);
    phiInvMass_08_10->Draw("E");
    likeSignInvMass_08_10->Draw("SAME");
    c0->cd(4);
    phiInvMass_10_12->Draw("E");
    likeSignInvMass_10_12->Draw("SAME");
    c0->cd(5);
    phiInvMass_12_14->Draw("E");
    likeSignInvMass_12_14->Draw("SAME");
    c0->cd(6);
    phiInvMass_14_16->Draw("E");
    likeSignInvMass_14_16->Draw("SAME");
    c0->cd(7);
    phiInvMass_16_18->Draw("E");
    likeSignInvMass_16_18->Draw("SAME");
    c0->cd(8);
    phiInvMass_18_20->Draw("E");
    likeSignInvMass_18_20->Draw("SAME");
    c0->cd(9);
    phiInvMass_20_25->Draw("E");
    likeSignInvMass_20_25->Draw("SAME");
    c0->cd(10);
    phiInvMass_25_30->Draw("E");
    likeSignInvMass_25_30->Draw("SAME");
    c0->cd(11);
    phiInvMass_30_40->Draw("E");
    likeSignInvMass_30_40->Draw("SAME");
    c0->cd(12);
    phiInvMass_40_50->Draw("E");
    likeSignInvMass_40_50->Draw("SAME");
    c0->cd(13);
    phiInvMass_50_60->Draw("E");
    likeSignInvMass_50_60->Draw("SAME");
    c0->cd(14);
    phiInvMass_60_70->Draw("E");
    likeSignInvMass_60_70->Draw("SAME");
    c0->cd(15);
    phiInvMass_70_100->Draw("E");
    likeSignInvMass_70_100->Draw("SAME");

    TCanvas *c1 = new TCanvas("cCorrInvMass", "cCorrInvMass", 50, 50, 800, 800);
    c1->Divide(4,4);

    for(int k =0; k < 8; k++){
        c1->cd(k+1);
        corrInvMass[k]->Draw();
    }
    /* 
    THnSparseF *dphiHPhi = (THnSparseF *)InvMass->FindObject("fDphiHPhi");
    THnSparseF *dphiHKstar = (THnSparseF *)InvMass->FindObject("fDphiHKstar");
    THnSparseF *dphiHK = (THnSparseF *)InvMass->FindObject("fDphiHK");
    THnSparseF *dphiHPi = (THnSparseF *)InvMass->FindObject("fDphiHPi");
    THnSparseF *dphiHp = (THnSparseF *)InvMass->FindObject("fDphiHp");
    THnSparseF *dphiHK0 = (THnSparseF *)InvMass->FindObject("fDphiHK0");

    TH3D *HPhiDphi = dphiHPhi->Projection(0,1,2);
    TH3D *HKstarDphi = dphiHKstar->Projection(0,1,2);
    TH3D *HKDphi = dphiHK->Projection(0,1,2);
    TH3D *HPiDphi = dphiHPi->Projection(0,1,2);
    TH3D *HpDphi = dphiHp->Projection(0,1,2);
    TH3D *HK0Dphi = dphiHK0->Projection(0,1,2);
*/
    // Graph of BG corrected inv-mass distribution
/*
    if(phiInvMass){
        TH1F *bgCorPhiMass = phiInvMass->Projection(1);
    }else{
        printf("No Histograms set up!");
        break;
    }

    //Doing Background scaling/subtraction
    Double_t num_unlike = phiInvMass->GetEntries();
    Double_t sideband = 0.0;
    Double_t bin_counter = 0.0;
    for(int i = 0; i < bgCorPhiMass->GetNbinsX(); i++){
        bgCorPhiMass->SetBinContent(i, bgCorPhiMass->GetBinContent(i) - (likeSignInvMass->GetBinContent(i))*num_unlike/(2.0*TMath::Sqrt(likeSignCounter->GetBinContent(likeSignCounter->FindBin(-1))*likeSignCounter->GetBinContent(likeSignCounter->FindBin(1)))));
        if(i<=bgCorPhiMass->FindBin(1.1) && i>=bgCorPhiMass->FindBin(1.05)){
            sideband+= bgCorPhiMass->GetBinContent(i);
            bin_counter++;
        }
    }
    sideband = sideband/(bin_counter);
    printf("sideband: %f\n", sideband);
    //Trying additional side-band subtraction?
//    for(int i = 0; i < bgCorPhiMass->GetNbinsX(); i++){
//        bgCorPhiMass->SetBinContent(i, bgCorPhiMass->GetBinContent(i) - sideband);
//    }

    TCanvas *cInv = new TCanvas("bgCorInvMass", "Like-Sign BG corrected InvMass distribution", 50, 50, 1000, 1000);
    cInv->cd();
    bgCorPhiMass->SetLineColor(1);
    bgCorPhiMass->SetFillColor(16);
    bgCorPhiMass->DrawCopy();
    bgCorPhiMass->Sumw2();
    bgCorPhiMass->SetMarkerStyle(21);
    bgCorPhiMass->SetMarkerSize(1);
    bgCorPhiMass->Draw("SAME E1");
*/
/*
    TH1D *dphiHPhiArray[5];
    dphiHPhiArray[0] = HPhiDphi->ProjectionZ("ptH4_ptPhiInc", 20, 100, 0, 100);
    dphiHPhiArray[0]->SetTitle("p_{T}^{H} > 4 GeV/c, p_{T}^{#Phi} Inclusive");
    dphiHPhiArray[1] = HPhiDphi->ProjectionZ("ptH4_1ptPhi2", 20, 100, 5, 10);
    dphiHPhiArray[1]->SetTitle("p_{T}^{H} > 4 GeV/c, 1 < p_{T}^{#Phi} < 2 GeV/c");
    dphiHPhiArray[2] = HPhiDphi->ProjectionZ("ptH4_2ptPhi4", 20, 100, 10, 20);
    dphiHPhiArray[2]->SetTitle("p_{T}^{H} > 4 GeV/c, 2 < p_{T}^{#Phi} < 4 GeV/c");
    dphiHPhiArray[3] = HPhiDphi->ProjectionZ("2ptH_ptPhi2", 0, 10, 0, 10);
    dphiHPhiArray[3]->SetTitle("p_{T}^{H} < 2 GeV/c, p_{T}^{#Phi} < 2 GeV/c");
    dphiHPhiArray[4] = HPhiDphi->ProjectionZ("ptH4_2ptPhi", 20, 100, 10, 100);
    dphiHPhiArray[4]->SetTitle("p_{T}^{H} > 4 GeV/c, p_{T}^{#Phi} > 2 GeV/c");

    TCanvas *cDphiHPhi = new TCanvas("cDphiHPhi", "cDphiHPhi", 50, 50, 600, 600);
    cDphiHPhi->Divide(2,2);
    for(int i =0; i<4; i++){
        cDphiHPhi->cd(i+1);
        dphiHPhiArray[i]->Draw("SAME");
    }  

    TH1D *dphiHKstarArray[5];
    dphiHKstarArray[0] = HKstarDphi->ProjectionZ("ptH4_ptKstarInc", 20, 100, 0, 100);
    dphiHKstarArray[0]->SetTitle("p_{T}^{H} > 4 GeV/c, p_{T}^{K*} Inclusive");
    dphiHKstarArray[1] = HKstarDphi->ProjectionZ("ptH4_1ptKstar2", 20, 100, 5, 10);
    dphiHKstarArray[1]->SetTitle("p_{T}^{H} > 4 GeV/c, 1 < p_{T}^{K*} < 2 GeV/c");
    dphiHKstarArray[2] = HKstarDphi->ProjectionZ("ptH4_2ptKstar4", 20, 100, 10, 20);
    dphiHKstarArray[2]->SetTitle("p_{T}^{H} > 4 GeV/c, 2 < p_{T}^{K*} < 4 GeV/c");
    dphiHKstarArray[3] = HKstarDphi->ProjectionZ("2ptH_ptKstar2", 0, 10, 0, 10);
    dphiHKstarArray[3]->SetTitle("p_{T}^{H} < 2 GeV/c, p_{T}^{K*} < 2 GeV/c");
    dphiHKstarArray[4] = HKstarDphi->ProjectionZ("ptH4_2ptKstar", 20, 100, 10, 100);
    dphiHKstarArray[4]->SetTitle("p_{T}^{H} > 4 GeV/c, p_{T}^{K*} > 2 GeV/c");

  
    TCanvas *cDphiHKstar = new TCanvas("cDphiHKstar", "cDphiHKstar", 50, 50, 600, 600);
    cDphiHKstar->Divide(2,2);
    for(int i =0; i<4; i++){
        cDphiHKstar->cd(i+1);
        dphiHKstarArray[i]->Draw("SAME");
    }  

    TH1D *dphiHKArray[5];
    dphiHKArray[0] = HKDphi->ProjectionZ("ptH4_ptKInc", 20, 100, 0, 100);
    dphiHKArray[0]->SetTitle("p_{T}^{H} > 4 GeV/c, p_{T}^{K} Inclusive");
    dphiHKArray[1] = HKDphi->ProjectionZ("ptH4_1ptK2", 20, 100, 5, 10);
    dphiHKArray[1]->SetTitle("p_{T}^{H} > 4 GeV/c, 1 < p_{T}^{K} < 2 GeV/c");
    dphiHKArray[2] = HKDphi->ProjectionZ("ptH4_2ptK4", 20, 100, 10, 20);
    dphiHKArray[2]->SetTitle("p_{T}^{H} > 4 GeV/c, 2 < p_{T}^{K} < 4 GeV/c");
    dphiHKArray[3] = HKDphi->ProjectionZ("2ptH_ptK2", 0, 10, 0, 10);
    dphiHKArray[3]->SetTitle("p_{T}^{H} < 2 GeV/c, p_{T}^{K} < 2 GeV/c");
    dphiHKArray[4] = HKDphi->ProjectionZ("ptH4_2ptK", 20, 100, 10, 100);
    dphiHKArray[4]->SetTitle("p_{T}^{H} > 4 GeV/c, p_{T}^{K} > 2 GeV/c");


    TCanvas *cDphiHK = new TCanvas("cDphiHK", "cDphiHK", 50, 50, 600, 600);
    cDphiHK->Divide(2,2);
    for(int i =0; i<4; i++){
        cDphiHK->cd(i+1);
        dphiHKArray[i]->Draw("SAME");
    }  

    TH1D *dphiHPiArray[5];
    dphiHPiArray[0] = HPiDphi->ProjectionZ("ptH4_ptPiInc", 20, 100, 0, 100);
    dphiHPiArray[0]->SetTitle("p_{T}^{H} > 4 GeV/c, p_{T}^{#pi} Inclusive");
    dphiHPiArray[1] = HPiDphi->ProjectionZ("ptH4_1ptPi2", 20, 100, 5, 10);
    dphiHPiArray[1]->SetTitle("p_{T}^{H} > 4 GeV/c, 1 < p_{T}^{#pi} < 2 GeV/c");
    dphiHPiArray[2] = HPiDphi->ProjectionZ("ptH4_2ptPi4", 20, 100, 10, 20);
    dphiHPiArray[2]->SetTitle("p_{T}^{H} > 4 GeV/c, 2 < p_{T}^{#pi} < 4 GeV/c");
    dphiHPiArray[3] = HPiDphi->ProjectionZ("2ptH_ptPi2", 0, 5, 0, 10);
    dphiHPiArray[3]->SetTitle("p_{T}^{H} < 1 GeV/c, p_{T}^{#pi} < 2 GeV/c");
    dphiHPiArray[4] = HPiDphi->ProjectionZ("ptH4_2ptPi", 20, 100, 10, 100);
    dphiHPiArray[4]->SetTitle("p_{T}^{H} > 4 GeV/c, p_{T}^{#pi} > 2 GeV/c");


    TCanvas *cDphiHPi = new TCanvas("cDphiHPi", "cDphiHPi", 50, 50, 600, 600);
    cDphiHPi->Divide(2,2);
    for(int i =0; i<4; i++){
        cDphiHPi->cd(i+1);
        dphiHPiArray[i]->Draw("SAME");
    }  

    TH1D *dphiHpArray[5];
    dphiHpArray[0] = HpDphi->ProjectionZ("ptH4_ptpInc", 20, 100, 0, 100);
    dphiHpArray[0]->SetTitle("p_{T}^{H} > 4 GeV/c, p_{T}^{p} Inclusive");
    dphiHpArray[1] = HpDphi->ProjectionZ("ptH4_1ptp2", 20, 100, 5, 10);
    dphiHpArray[1]->SetTitle("p_{T}^{H} > 4 GeV/c, 1 < p_{T}^{p} < 2 GeV/c");
    dphiHpArray[2] = HpDphi->ProjectionZ("ptH4_2ptp4", 20, 100, 10, 20);
    dphiHpArray[2]->SetTitle("p_{T}^{H} > 4 GeV/c, 2 < p_{T}^{p} < 4 GeV/c");
    dphiHpArray[3] = HpDphi->ProjectionZ("2ptH_ptp2", 0, 5, 0, 10);
    dphiHpArray[3]->SetTitle("p_{T}^{H} < 1 GeV/c, p_{T}^{p} < 2 GeV/c");
    dphiHpArray[4] = HpDphi->ProjectionZ("ptH4_2ptp", 20, 100, 10, 100);
    dphiHpArray[4]->SetTitle("p_{T}^{H} > 4 GeV/c, p_{T}^{p} > 2 GeV/c");

  
    TCanvas *cDphiHp = new TCanvas("cDphiHp", "cDphiHp", 50, 50, 600, 600);
    cDphiHp->Divide(2,2);
    for(int i =0; i<4; i++){
        cDphiHp->cd(i+1);
        dphiHpArray[i]->Draw("SAME");
    }

    TH1D *dphiHK0Array[5];
    dphiHK0Array[0] = HK0Dphi->ProjectionZ("ptH4_ptK0Inc", 20, 100, 0, 100);
    dphiHK0Array[0]->SetTitle("p_{T}^{H} > 4 GeV/c, p_{T}^{K0} Inclusive");
    dphiHK0Array[1] = HK0Dphi->ProjectionZ("ptH4_1ptK02", 20, 100, 5, 10);
    dphiHK0Array[1]->SetTitle("p_{T}^{H} > 4 GeV/c, 1 < p_{T}^{K0} < 2 GeV/c");
    dphiHK0Array[2] = HK0Dphi->ProjectionZ("ptH4_2ptK04", 20, 100, 10, 20);
    dphiHK0Array[2]->SetTitle("p_{T}^{H} > 4 GeV/c, 2 < p_{T}^{K0} < 4 GeV/c");
    dphiHK0Array[3] = HK0Dphi->ProjectionZ("2ptH_ptK02", 0, 10, 0, 10);
    dphiHK0Array[3]->SetTitle("p_{T}^{H} < 2 GeV/c, p_{T}^{K0} < 2 GeV/c");
    dphiHK0Array[4] = HK0Dphi->ProjectionZ("ptH4_2ptK0", 20, 100, 10, 100);
    dphiHK0Array[4]->SetTitle("p_{T}^{H} > 4 GeV/c, p_{T}^{K0} > 2 GeV/c");


    TCanvas *cDphiHK0 = new TCanvas("cDphiHK0", "cDphiHK0", 50, 50, 600, 600);
    cDphiHK0->Divide(2,2);
    for(int i =0; i<4; i++){
        cDphiHK0->cd(i+1);
        dphiHK0Array[i]->Draw("SAME");
    }  

//Putting several different Species on one canvas for comparisons
    TCanvas *cDphiAll = new TCanvas("cDphiAll", "cDphiAll", 50, 50, 600, 600);
    cDphiAll->Divide(2,2);
    cDphiAll->cd(1);
    dphiHK0Array[4]->Draw("SAME");
    cDphiAll->cd(2);
    dphiHpArray[2]->Draw("SAME");
    cDphiAll->cd(3);
    dphiHKstarArray[4]->Draw("SAME");
    cDphiAll->cd(4);
    dphiHPhiArray[4]->Draw("SAME");

//Setting up ratio between h-phi and h-K0
    TH1D* dphiHPhiRebin = dphiHPhiArray[4]->Rebin(16, "dphiHPhiRebin");
    TH1D* dphiHK0Rebin = dphiHK0Array[4]->Rebin(16, "dphiHK0Rebin");

    TH1D* dphiRatio = dphiHK0Rebin->Clone("dphiRatio");
    dphiRatio->Divide(dphiHPhiRebin);
    dphiRatio->SetTitle("Ratio of K_{S}^{0} / #Phi(1020)");

    TCanvas *cHPhiRebin = new TCanvas("cHPhiRebin", "cHPhiRebin", 50, 50, 600, 600);
    cHPhiRebin->cd();
    dphiHPhiRebin->Draw();
    dphiHPhiRebin->Draw("E SAME");

    TCanvas *cHK0Rebin = new TCanvas("cHK0Rebin", "cHK0Rebin", 50, 50, 600, 600);
    cHK0Rebin->cd();
    dphiHK0Rebin->Draw();
    dphiHK0Rebin->Draw("E SAME");

    TCanvas *cRatio = new TCanvas("cRatio", "cRatio", 50,50, 600,600);
    cRatio->cd();
    dphiRatio->Draw();
    */
}
