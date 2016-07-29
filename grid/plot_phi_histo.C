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
    TH1D *phiInvMassBinned[16];
    TH1D *likeSignInvMassBinned[16];
    TF1 *fits[16];
    for(int n = 0; n < 16; n++){
        fits[n] = new TF1(Form("f%i", n), "gaus", 1.005,1.035);
    } 

    Double_t pt_bounds[] = {0.4, 0.6, 0.8, 1.0, 1.2, 1.4, 1.6, 1.8, 2.0, 2.5, 3.0, 4.0, 5.0, 6.0, 7.0, 10.0};
    Int_t pt_bins[] = {4, 6, 8, 10, 12, 14, 16, 18, 20, 25, 30, 40, 50, 60, 70, 100};
    Double_t sigmas[16];    

    Double_t sideband = 0.0;
    Double_t likeSignSideBand = 0.0;
    Double_t scaleFactor = 0.0;

    if(phiInvMass && likeSignInvMass){
        phiInvMass->GetAxis(1)->SetRange(320,400);
        phiInvMass->GetAxis(1)->SetTitle("Inv Mass (GeV/c^2)");
        likeSignInvMass->GetAxis(1)->SetRange(320,400);
        likeSignInvMass->SetTitle("Inv Mass (GeV/c^2)");

        for(int i =0; i<15; i++){
            phiInvMass->GetAxis(0)->SetRange(pt_bins[i], pt_bins[i+1]);
            phiInvMassBinned[i] = (TH1D*)phiInvMass->Projection(1);
            phiInvMassBinned[i]->SetTitle(Form("%.1f < p_{T} < %.1f GeV/c", pt_bounds[i], pt_bounds[i+1]));
            sideband = phiInvMassBinned[i]->Integral(61,81); //integrating from mass of 1.04 to ~1.06 for the sideband scaling
            likeSignInvMass->GetAxis(0)->SetRange(pt_bins[i], pt_bins[i+1]);
            likeSignInvMassBinned[i] = (TH1D*)likeSignInvMass->Projection(1);
            likeSignInvMassBinned[i]->SetTitle(Form("%.1f < p_{T} < %.1f GeV/c", pt_bounds[i], pt_bounds[i+1]));
            likeSignSideBand = likeSignInvMassBinned[i]->Integral(61,81);
            likeSignInvMassBinned[i]->Scale(sideband/likeSignSideBand);
            likeSignInvMassBinned[i]->SetLineColor(2);
            corrInvMass[i] = (TH1D*)phiInvMassBinned[i]->Clone(Form("corrInvMass_%i_%i", pt_bins[i], pt_bins[i]));
            corrInvMass[i]->Add(likeSignInvMassBinned[i], -1.0);
            corrInvMass[i]->Fit(Form("f%i", i), "R");
            sigmas[i] = fits[i]->GetParameter("Sigma");
        }
    }   

    TCanvas *c0 = new TCanvas("cTest", "cTest", 50, 50, 800, 800);
    c0->Divide(4,4);

    for(int j = 0; j< 15; j++){
        c0->cd(j+1);
        phiInvMassBinned[j]->Draw("E");
        likeSignInvMassBinned[j]->Draw("SAME");
    }


    TCanvas *c1 = new TCanvas("cCorrInvMass", "cCorrInvMass", 50, 50, 800, 800);
    c1->Divide(4,4);

    for(int k =0; k < 15; k++){
        c1->cd(k+1);
        printf("sigma: %E\n", sigmas[k]);
        corrInvMass[k]->Draw();
        fits[k]->Draw("SAME");
    }
     
    THnSparseF *dphiHPhi = (THnSparseF *)InvMass->FindObject("fDphiHPhi");
    THnSparseF *dphiHKK = (THnSparseF *)InvMass->FindObject("fDphiHKK");

    /*
    THnSparseF *dphiHKstar = (THnSparseF *)InvMass->FindObject("fDphiHKstar");
    THnSparseF *dphiHK = (THnSparseF *)InvMass->FindObject("fDphiHK");
    THnSparseF *dphiHPi = (THnSparseF *)InvMass->FindObject("fDphiHPi");
    THnSparseF *dphiHp = (THnSparseF *)InvMass->FindObject("fDphiHp");
    THnSparseF *dphiHK0 = (THnSparseF *)InvMass->FindObject("fDphiHK0");
    */
    TH3D *HPhiDphi = dphiHPhi->Projection(0,1,2);
    
    dphiHPhi->GetAxis(0)->SetRange(20, 100); //cutting on trigger pt greater than 4 GeV/c
    dphiHPhi->GetAxis(1)->SetRange(5,15); //cutting on phi pt between 1 GeV/c and 3 GeV/c
    
    //performing same cuts for KK pairs:
    dphiHKK->GetAxis(0)->SetRange(20, 100);
    dphiHKK->GetAxis(1)->SetRange(5,15);

    TH1D *phiInvMassPerDPhi[16];
    TH1D *likesignInvMassPerDPhi[16]; 
    TH1D *corrPhiInvMassPerDPhi[16];
    TCanvas *cInvDPhi = new TCanvas("cInvDPhi", "cIncDPhi", 50,50,800,800);
    cInvDPhi->Divide(4,4);
    for(int i = 0; i < 16; i++){
        scaleFactor=0.0;
        cInvDPhi->cd(i+1);
        dphiHPhi->GetAxis(2)->SetRange((4*i)+1, 4*(i+1));
        dphiHKK->GetAxis(2)->SetRange((4*i)+1, 4*(i+1));
        
        phiInvMassPerDPhi[i] = dphiHPhi->Projection(4);
        likesignInvMassPerDPhi[i] = dphiHKK->Projection(4);
        sideband = phiInvMassPerDPhi[i]->Integral(61, 81);
        likeSignSideBand = likesignInvMassPerDPhi[i]->Integral(61,81);

        phiInvMassPerDPhi[i]->SetTitle(Form("Inv Mass dist. for %.2f < #Delta#Phi < %.2f", (-1.57 + i*(0.3925)), (-1.57+(i+1)*(0.3925))));
        phiInvMassPerDPhi[i]->Draw();

        scaleFactor = sideband/likeSignSideBand;
        likesignInvMassPerDPhi[i]->Scale(scaleFactor);
        likeSignSideBand = likesignInvMassPerDPhi[i]->Integral(61,81);
        printf("%d. sideband: %f, likeSignSideBand: %f, difference: %f\n", i, sideband, likeSignSideBand, sideband-likeSignSideBand);
        likesignInvMassPerDPhi[i]->SetLineColor(2);
        likesignInvMassPerDPhi[i]->Draw("SAME");
        corrPhiInvMassPerDPhi[i] = (TH1D*)phiInvMassPerDPhi[i]->Clone();
        corrPhiInvMassPerDPhi[i]->Add(likesignInvMassPerDPhi[i], -1.0);
    }

    TCanvas *cCorrInvDPhi = new TCanvas("cCorrInvDPhi", "cCorrInvDPhi", 40, 40, 800, 800);
    cCorrInvDPhi->Divide(4,4);
    TF1 *dPhiFits[16];
    double numPhi[16];
    TH1D* corrDPhi = new TH1D("corrDPhi", "BG-Corrected Hadron-Phi #Delta#Phi distribution", 16, -1.57, 4.71);
    for(int i = 0; i<16; i++){
        cCorrInvDPhi->cd(i+1);
        //dPhiFits[i] = new TF1(Form("dphif%i", i), "gaus", 1.005,1.035); 
        //dPhiFits[i]->SetParameter(1, 1.020);
        //dPhiFits[i]->SetParameter(2, 0.0045);
        //dPhiFits[i]->SetParameter(0, 200);
        //dPhiFits[i]->SetParLimits(1, 1.015, 1.025);
        //dPhiFits[i]->SetParLimits(2, 0.0010, 0.007);
        //dPhiFits[i]->SetParLimits(0, 100, 10000);
        corrPhiInvMassPerDPhi[i]->Draw();
        //corrPhiInvMassPerDPhi[i]->Fit(Form("dphif%i", i), "R");
        //printf("%i. Total: %f\n", i, dPhiFits[i]->Integral(1.005, 1.035));
        numPhi[i] = corrPhiInvMassPerDPhi[i]->Integral(31, 51);
        corrDPhi->SetBinContent(i+1, numPhi[i]);
    }

    TCanvas *cCorrDPhi = new TCanvas("cCorrDPhi", "cCorrDPhi", 30, 30, 600, 600);
    cCorrDPhi->cd();
    corrDPhi->Draw();

    dphiHPhi->GetAxis(2)->SetRange(0,0);

    TH1D *dphiPhiPeak;
    dphiHPhi->GetAxis(4)->SetRangeUser(1.010, 1.030);
    dphiPhiPeak = dphiHPhi->Projection(2);
    dphiPhiPeak->SetTitle("Hadron-Phi #Delta#Phi for Phi Peak"); 

    TH1D *dphiPhiSideband;
    dphiHPhi->GetAxis(4)->SetRangeUser(1.040, 1.060);
    dphiPhiSideband = dphiHPhi->Projection(2);
    dphiPhiSideband->SetTitle("Hadron-Phi #Delta#Phi for Phi Sideband");
    
    /*
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
    
    TCanvas *cDphiHPhi2 = new TCanvas("cDphiHPhi2", "cDphiHPhi2", 50, 50, 600, 600);
    cDphiHPhi2->cd();
    dphiPhiPeak->Draw("SAME");
    
    TCanvas *cDphiHPhi3 = new TCanvas("cDphiHPhi3", "cDphiHPhi3", 50, 50, 600, 600);
    dphiPhiSideband->Draw("SAME");

    /*
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
