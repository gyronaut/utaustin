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
            scaleFactor = sideband/likeSignSideBand;
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
        //printf("sigma: %E\n", sigmas[k]);
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
        phiInvMassPerDPhi[i]->Sumw2();
        sideband = phiInvMassPerDPhi[i]->Integral(61, 81);
        likeSignSideBand = likesignInvMassPerDPhi[i]->Integral(61,81);

        phiInvMassPerDPhi[i]->SetTitle(Form("Inv Mass dist. for %.2f < #Delta#Phi < %.2f", (-1.57 + i*(0.3925)), (-1.57+(i+1)*(0.3925))));
        phiInvMassPerDPhi[i]->Draw("HIST E");

        scaleFactor = sideband/likeSignSideBand;
        likesignInvMassPerDPhi[i]->Scale(scaleFactor);
        likesignInvMassPerDPhi[i]->Sumw2();
        likesignInvMassPerDPhi[i]->SetLineColor(2);
        likesignInvMassPerDPhi[i]->Draw("SAME HIST E");
        
        //For corrected plots, need to rebin the distributions to get rid of some of the noise
        phiInvMassPerDPhi[i]->Rebin(4);
        phiInvMassPerDPhi[i]->Sumw2();
        likesignInvMassPerDPhi[i]->Rebin(4);
        likesignInvMassPerDPhi[i]->Sumw2();
        corrPhiInvMassPerDPhi[i] = (TH1D*)phiInvMassPerDPhi[i]->Clone();
        corrPhiInvMassPerDPhi[i]->Add(likesignInvMassPerDPhi[i], -1.0);
    }

    TCanvas *cCorrInvDPhi = new TCanvas("cCorrInvDPhi", "cCorrInvDPhi", 40, 40, 800, 800);
    cCorrInvDPhi->Divide(4,4);
    TF1 *dPhiFits[16];
    Double_t numPhi[16];
    TH1D* corrDPhi = new TH1D("corrDPhi", "BG-Corrected Hadron-Phi #Delta#Phi distribution", 16, -1.57, 4.71);
    for(int i = 0; i<16; i++){
        cCorrInvDPhi->cd(i+1);
        dPhiFits[i] = new TF1(Form("dphif%i", i), "gaus(0) + pol0(3)",0.98,1.1); 
        dPhiFits[i]->SetParameter(1, 1.020);
        dPhiFits[i]->SetParameter(2, 0.0045);
        dPhiFits[i]->SetParameter(0, 200);
        //dPhiFits[i]->SetParameter(3, 0);
        dPhiFits[i]->SetParLimits(1, 1.015, 1.025);
        dPhiFits[i]->SetParLimits(2, 0.0010, 0.007);
        dPhiFits[i]->SetParLimits(0, 100, 10000);
        //dPhiFits[i]->SetParLimits(3, -100, 100);
        corrPhiInvMassPerDPhi[i]->Fit(Form("dphif%i", i), "R");
        corrPhiInvMassPerDPhi[i]->Draw("SAME HIST E FUNC");
        //printf("%i. Total: %f\n", i, dPhiFits[i]->Integral(1.005, 1.035));
        printf("correction: %f", dPhiFits[i]->GetParameter(3)*3.0);
        numPhi[i] = corrPhiInvMassPerDPhi[i]->Integral(9, 12) - dPhiFits[i]->GetParameter(3)*3.0;
        printf("Corrected num Phis: %f\n", numPhi[i]);
        corrDPhi->SetBinContent(i+1, numPhi[i]);
    }

    TCanvas *cCorrDPhi = new TCanvas("cCorrDPhi", "cCorrDPhi", 30, 30, 600, 600);
    cCorrDPhi->cd();
    corrDPhi->Draw();

    //Reset the Delta-phi axis range after the dphi binned projections are done being created above.
    dphiHPhi->GetAxis(2)->SetRange(0,0);
    dphiHKK->GetAxis(2)->SetRange(0,0);

    TH1D *dphiPhiPeak;
    dphiHPhi->GetAxis(4)->SetRangeUser(1.010, 1.030);
    dphiPhiPeak = dphiHPhi->Projection(2);
    dphiPhiPeak->SetTitle("Hadron-Phi #Delta#Phi for Phi Peak"); 

    TH1D *dphiPhiSideband;
    dphiHPhi->GetAxis(4)->SetRangeUser(1.040, 1.060);
    dphiPhiSideband = dphiHPhi->Projection(2);
    dphiPhiSideband->SetTitle("Hadron-Phi #Delta#Phi for Phi Sideband");
 
// Setting up different trigger/assoc momentum cuts for delta-phi distribution
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

// Plotting delta-phi distribution for different trigger/assoc. momentum cuts
    TCanvas *cDphiHPhi = new TCanvas("cDphiHPhi", "cDphiHPhi", 50, 50, 600, 600);
    cDphiHPhi->Divide(2,2);
    for(int i =0; i<4; i++){
        cDphiHPhi->cd(i+1);
        dphiHPhiArray[i]->Draw("SAME");
    }  
    
// Setting & plotting Delta Phi distribution for Unlike sign and Likesign pairs in the
// sideband region (scaled by the ratio of the integral of the sideband region)
    dphiHPhi->GetAxis(4)->SetRangeUser(1.040,1.060);
    TH1D* dphiUSSideband = dphiHPhi->Projection(2);
    dphiHKK->GetAxis(4)->SetRangeUser(1.040,1.060);
    TH1D* dphiLSSideband = dphiHKK->Projection(2);
    sideband = dphiHPhi->Projection(4)->Integral(15, 20);
    likeSignSideBand = dphiHKK->Projection(4)->Integral(15, 20);
    scaleFactor = sideband/likeSignSideBand;
    //printf("sideband: %f, likeSignSideBand: %f, scaleFactor: %f \n", sideband, likeSignSideBand, scaleFactor);
    dphiLSSideband->Scale(scaleFactor);

    TCanvas *cDphiHPhi3 = new TCanvas("cDphiHPhi3", "cDphiHPhi3", 50, 50, 600, 600);
    cDphiHPhi3->cd();
    dphiUSSideband->SetLineWidth(2);
    dphiUSSideband->SetTitle("Hadron-Phi #Delta#Phi for Sideband Region");
    dphiUSSideband->Draw("SAME");
    dphiLSSideband->SetLineColor(2);
    dphiLSSideband->Draw("SAME");


// Setting & plotting delta-phi distribution for Unlike and Like sign pairs in the phi mass peak
// region (scaled by the ratio of the integral of the sideband region).
    dphiHPhi->GetAxis(4)->SetRangeUser(1.010,1.030);
    TH1D* dphiUSPeak = dphiHPhi->Projection(2);
    dphiHKK->GetAxis(4)->SetRangeUser(1.010,1.030);
    TH1D* dphiLSPeak = dphiHKK->Projection(2);
    dphiLSPeak->Scale(scaleFactor);

    TCanvas *cDphiHPhi2 = new TCanvas("cDphiHPhi2", "cDphiHPhi2", 50, 50, 600, 600);
    cDphiHPhi2->cd();
    dphiUSPeak->SetLineWidth(2);
    dphiUSPeak->SetTitle("Hadron-Phi #Delta#Phi for Peak Region");
    dphiUSPeak->Draw("SAME");
    dphiLSPeak->SetLineColor(2);
    dphiLSPeak->Draw("SAME");
}
