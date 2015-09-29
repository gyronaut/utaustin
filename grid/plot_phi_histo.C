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

    TFile *histoFile = new TFile(inputName.c_str());
    histoFile->cd("PhiReconstruction");
    THnSparseF *phiInvMass = (THnSparseF *)InvMass->FindObject("fPhiInvMass");
    THnSparseF *likeSignInvMass = (THnSparseF *)InvMass->FindObject("fPhiLikeSignInvMass");
   
    THnSparseF *dphiHPhi = (THnSparseF *)InvMass->FindObject("fDphiHPhi");
    THnSparseF *dphiHKstar = (THnSparseF *)InvMass->FindObject("fDphiHKstar");
    THnSparseF *dphiHK = (THnSparseF *)InvMass->FindObject("fDphiHK");
    THnSparseF *dphiHPi = (THnSparseF *)InvMass->FindObject("fDphiHPi");
    THnSparseF *dphiHp = (THnSparseF *)InvMass->FindObject("fDphiHp");

    TH3D *HPhiDphi = dphiHPhi->Projection(0,1,2);
    TH3D *HKstarDphi = dphiHKstar->Projection(0,1,2);
    TH3D *HKDphi = dphiHK->Projection(0,1,2);
    TH3D *HPiDphi = dphiHPi->Projection(0,1,2);
    TH3D *HpDphi = dphiHp->Projection(0,1,2);

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
    TH1D *dphiHPhiArray[3];
    dphiHPhiArray[0] = HPhiDphi->ProjectionZ("ptH4_ptPhiInc", 20, 100, 0, 100);
    dphiHPhiArray[0]->SetTitle("p_{T}^{H} > 4 GeV/c, p_{T}^{#Phi} Inclusive");
    dphiHPhiArray[1] = HPhiDphi->ProjectionZ("ptH4_1ptPhi2", 20, 100, 5, 10);
    dphiHPhiArray[1]->SetTitle("p_{T}^{H} > 4 GeV/c, 1 < p_{T}^{#Phi} < 2 GeV/c");
    dphiHPhiArray[2] = HPhiDphi->ProjectionZ("ptH4_2ptPhi4", 20, 100, 10, 20);
    dphiHPhiArray[2]->SetTitle("p_{T}^{H} > 4 GeV/c, 2 < p_{T}^{#Phi} < 4 GeV/c");

    TCanvas *cDphiHPhi = new TCanvas("cDphiHPhi", "cDphiHPhi", 50, 50, 600, 600);
    cDphiHPhi->Divide(2,2);
    for(int i =0; i<3; i++){
        cDphiHPhi->cd(i+1);
        dphiHPhiArray[i]->Draw("SAME");
    }  

    TH1D *dphiHKstarArray[3];
    dphiHKstarArray[0] = HKstarDphi->ProjectionZ("ptH4_ptKstarInc", 20, 100, 0, 100);
    dphiHKstarArray[0]->SetTitle("p_{T}^{H} > 4 GeV/c, p_{T}^{K*} Inclusive");
    dphiHKstarArray[1] = HKstarDphi->ProjectionZ("ptH4_1ptKstar2", 20, 100, 5, 10);
    dphiHKstarArray[1]->SetTitle("p_{T}^{H} > 4 GeV/c, 1 < p_{T}^{K*} < 2 GeV/c");
    dphiHKstarArray[2] = HKstarDphi->ProjectionZ("ptH4_2ptKstar4", 20, 100, 10, 20);
    dphiHKstarArray[2]->SetTitle("p_{T}^{H} > 4 GeV/c, 2 < p_{T}^{K*} < 4 GeV/c");
    
    TCanvas *cDphiHKstar = new TCanvas("cDphiHKstar", "cDphiHKstar", 50, 50, 600, 600);
    cDphiHKstar->Divide(2,2);
    for(int i =0; i<3; i++){
        cDphiHKstar->cd(i+1);
        dphiHKstarArray[i]->Draw("SAME");
    }  

    TH1D *dphiHKArray[3];
    dphiHKArray[0] = HKDphi->ProjectionZ("ptH4_ptKInc", 20, 100, 0, 100);
    dphiHKArray[0]->SetTitle("p_{T}^{H} > 4 GeV/c, p_{T}^{K} Inclusive");
    dphiHKArray[1] = HKDphi->ProjectionZ("ptH4_1ptK2", 20, 100, 5, 10);
    dphiHKArray[1]->SetTitle("p_{T}^{H} > 4 GeV/c, 1 < p_{T}^{K} < 2 GeV/c");
    dphiHKArray[2] = HKDphi->ProjectionZ("ptH4_2ptK4", 20, 100, 10, 20);
    dphiHKArray[2]->SetTitle("p_{T}^{H} > 4 GeV/c, 2 < p_{T}^{K} < 4 GeV/c");

    TCanvas *cDphiHK = new TCanvas("cDphiHK", "cDphiHK", 50, 50, 600, 600);
    cDphiHK->Divide(2,2);
    for(int i =0; i<3; i++){
        cDphiHK->cd(i+1);
        dphiHKArray[i]->Draw("SAME");
    }  

    TH1D *dphiHPiArray[3];
    dphiHPiArray[0] = HPiDphi->ProjectionZ("ptH4_ptPiInc", 20, 100, 0, 100);
    dphiHPiArray[0]->SetTitle("p_{T}^{H} > 4 GeV/c, p_{T}^{#pi} Inclusive");
    dphiHPiArray[1] = HPiDphi->ProjectionZ("ptH4_1ptPi2", 20, 100, 5, 10);
    dphiHPiArray[1]->SetTitle("p_{T}^{H} > 4 GeV/c, 1 < p_{T}^{#pi} < 2 GeV/c");
    dphiHPiArray[2] = HPiDphi->ProjectionZ("ptH4_2ptPi4", 20, 100, 10, 20);
    dphiHPiArray[2]->SetTitle("p_{T}^{H} > 4 GeV/c, 2 < p_{T}^{#pi} < 4 GeV/c");

    TCanvas *cDphiHPi = new TCanvas("cDphiHPi", "cDphiHPi", 50, 50, 600, 600);
    cDphiHPi->Divide(2,2);
    for(int i =0; i<3; i++){
        cDphiHPi->cd(i+1);
        dphiHPiArray[i]->Draw("SAME");
    }  

    TH1D *dphiHpArray[3];
    dphiHpArray[0] = HpDphi->ProjectionZ("ptH4_ptpInc", 20, 100, 0, 100);
    dphiHpArray[0]->SetTitle("p_{T}^{H} > 4 GeV/c, p_{T}^{p} Inclusive");
    dphiHpArray[1] = HpDphi->ProjectionZ("ptH4_1ptp2", 20, 100, 5, 10);
    dphiHpArray[1]->SetTitle("p_{T}^{H} > 4 GeV/c, 1 < p_{T}^{p} < 2 GeV/c");
    dphiHpArray[2] = HpDphi->ProjectionZ("ptH4_2ptp4", 20, 100, 10, 20);
    dphiHpArray[2]->SetTitle("p_{T}^{H} > 4 GeV/c, 2 < p_{T}^{p} < 4 GeV/c");
    
    TCanvas *cDphiHp = new TCanvas("cDphiHp", "cDphiHp", 50, 50, 600, 600);
    cDphiHp->Divide(2,2);
    for(int i =0; i<3; i++){
        cDphiHp->cd(i+1);
        dphiHpArray[i]->Draw("SAME");
    }  

}
