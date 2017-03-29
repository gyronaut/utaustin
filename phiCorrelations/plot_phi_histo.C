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
    THnSparseF *fkkUSDist = (THnSparseF *)InvMass->FindObject("fkkUSDist");
    THnSparseF *fkkLSDist = (THnSparseF *)InvMass->FindObject("fkkLSDist");
    THnSparseF *fTrigDist = (THnSparseF *)InvMass->FindObject("fTrigDist");
    THnSparseF *dphiHPhi = (THnSparseF *)InvMass->FindObject("fDphiHPhi");
    THnSparseF *dphiHKK = (THnSparseF *)InvMass->FindObject("fDphiHKK");

    TH1D *phiPTSpectrum;
    TH1D *corrInvMass[16];
    TH1D *phiInvMassBinned[16];
    TH1D *likeSignInvMassBinned[16];
    TF1 *fits[16];
    for(int n = 0; n < 16; n++){
        fits[n] = new TF1(Form("f%i", n),  "[0]*TMath::Voigt(x - [1], [2], [3], 4) + pol0(4)",0.99, 1.1);
        fits[n]->SetParameter(1, 1.020);
        fits[n]->SetParameter(2, 0.0002);
        fits[n]->SetParameter(0, 600);
        fits[n]->FixParameter(3, 0.00426);
        fits[n]->SetParLimits(1, 1.010, 1.030);
        //fits[n]->SetParLimits(0, 225, 400);
        //fits[n]->SetParLimits(2, 0.0010, 0.007);
        //fits[n]->SetParLimits(3, 0.0010, 0.010);

    } 

    Double_t pt_bounds[] = {0.4, 0.6, 0.8, 1.0, 1.2, 1.4, 1.6, 1.8, 2.0, 2.5, 3.0, 4.0, 5.0, 6.0, 7.0, 10.0};
    Int_t pt_bins[] = {4, 6, 8, 10, 12, 14, 16, 18, 20, 25, 30, 40, 50, 60, 70, 100};
    Double_t sigmas[16];    

    Double_t error;
    Double_t sideband = 0.0;
    Double_t likeSignSideBand = 0.0;
    Double_t scaleFactor = 0.0;
    Double_t invMin = 0.99;
    Double_t invMax = 1.1;
    Int_t numBins = 200 - (int)(((invMin-0.98)+(1.1-invMax))*200./(1.1-0.98));

    /* Setting up the invariant mass dist for LS and US Kaon pairs (scaled to the sideband region),
     * and a "BG corrected" invmass per pt bin*/
    phiPTSpectrum = new TH1D("phiPTspectrum", "p_{T}^{#Phi} Spectrum Corrected with LS BG", 15, pt_bounds);
    if(fkkUSDist && fkkLSDist){
        double bins_per_mass = numBins/(invMax-invMin);

        fkkUSDist->GetAxis(1)->SetRange(200-numBins,200);
        fkkUSDist->GetAxis(1)->SetTitle("Inv Mass (GeV/c^2)");
        //fkkUSDist->Sumw2();
        fkkLSDist->GetAxis(1)->SetRange(200-numBins,200);
        fkkLSDist->SetTitle("Inv Mass (GeV/c^2)");
        //fkkLSDist->Sumw2();
        
        for(int i =0; i<15; i++){
            fkkUSDist->GetAxis(0)->SetRange(pt_bins[i], pt_bins[i+1]);
            phiInvMassBinned[i] = (TH1D*)fkkUSDist->Projection(1);
            //phiInvMassBinned[i]->Sumw2();
            phiInvMassBinned[i]->SetTitle(Form("%.1f < p_{T} < %.1f GeV/c", pt_bounds[i], pt_bounds[i+1]));
            sideband = phiInvMassBinned[i]->Integral((int)(bins_per_mass*(1.04-invMin)),(int)(bins_per_mass*(1.06-invMin))); //integrating from mass of 1.04 to ~1.06 for the sideband scaling
            fkkLSDist->GetAxis(0)->SetRange(pt_bins[i], pt_bins[i+1]);
            likeSignInvMassBinned[i] = (TH1D*)fkkLSDist->Projection(1);
            //likeSignInvMassBinned[i]->Sumw2();
            likeSignInvMassBinned[i]->SetTitle(Form("%.1f < p_{T} < %.1f GeV/c", pt_bounds[i], pt_bounds[i+1]));
            likeSignSideBand = likeSignInvMassBinned[i]->Integral((int)(bins_per_mass*(1.04-invMin)),(int)(bins_per_mass*(1.06-invMin)));
            likeSignInvMassBinned[i]->Scale(sideband/likeSignSideBand);
            scaleFactor = sideband/likeSignSideBand;
            likeSignInvMassBinned[i]->SetLineColor(2);
            corrInvMass[i] = (TH1D*)phiInvMassBinned[i]->Clone(Form("corrInvMass_%i_%i", pt_bins[i], pt_bins[i]));
            corrInvMass[i]->Add(likeSignInvMassBinned[i], -1.0);
            //corrInvMass[i]->Sumw2();
            corrInvMass[i]->Fit(Form("f%i", i), "R");
            phiPTSpectrum->SetBinContent(i+1, (corrInvMass[i]->IntegralAndError((Int_t)(bins_per_mass*(1.01-invMin)),(Int_t)(bins_per_mass*(1.03-invMin)), error, "width"))/((pt_bounds[i+1] - pt_bounds[i])*1.6));
            //phiPTSpectrum->SetBinError(i, error);
            phiPTSpectrum->SetBinError(i+1, TMath::Sqrt(phiPTSpectrum->GetBinContent(i)));
        }
    }else{
      printf("couldn't open kk distributions!\n");
      return 0;
    }

    //Plot the Like-sign background corrected phi(1020) pT distribution:
    TCanvas *ptCanvas = new TCanvas("ptCanvas", "ptCanvas", 50, 50, 800, 800);
    ptCanvas->cd();
    phiPTSpectrum->GetXaxis()->SetTitle("p_{T} (GeV/c^{2})");
    phiPTSpectrum->GetYaxis()->SetTitle("d^{2}N/dp_{T}d#eta");
    phiPTSpectrum->SetMarkerStyle(3);
    phiPTSpectrum->SetMarkerSize(0.8);
    phiPTSpectrum->Draw();

    // Plotting distributions for KK pairs (US and LS) and trigger particles

    TCanvas *cLSDist = new TCanvas("cLSDist", "cLSDist", 50, 50, 800, 800);
    cLSDist->Divide(2,2);
    TCanvas *cUSDist = new TCanvas("cUSDist", "cUSDist", 50, 50, 800, 800);
    cUSDist->Divide(2,2);
    TCanvas *cTrigDist = new TCanvas("cTrigDist", "cTrigDist", 50, 50, 800, 800);
    cTrigDist->Divide(2,2);

    fkkLSDist->GetAxis(0)->SetRangeUser(2.0, 4.0);
    fkkLSDist->GetAxis(0)->SetTitle("p_{T} (GeV/c)");
    fkkUSDist->GetAxis(0)->SetRangeUser(2.0, 4.0);
    fkkUSDist->GetAxis(0)->SetTitle("p_{T} (GeV/c)");
   //fkkUSDist->GetAxis(1)->SetRangeUser(1.01, 1.03);
    fTrigDist->GetAxis(0)->SetRangeUser(4.0, 8.0);
    fTrigDist->GetAxis(0)->SetTitle("p_{T} (GeV/c)");
    
    fkkLSDist->GetAxis(2)->SetTitle("#varphi");
    fkkUSDist->GetAxis(2)->SetTitle("#varphi");

    fkkLSDist->GetAxis(3)->SetTitle("#eta");
    fkkUSDist->GetAxis(3)->SetTitle("#eta");
    fTrigDist->GetAxis(2)->SetTitle("#eta");

    TH1D* hLSDist[3];
    TH1D* hUSDist[3];
    TH1D* hTrigDist[2];
    for(int i=0; i<3; i++){
        hLSDist[i] = fkkLSDist->Projection(i+1);
        hUSDist[i] = fkkUSDist->Projection(i+1);
        if(i<2){
            hTrigDist[i] = fTrigDist->Projection(i+1);
        }
    }
    fkkLSDist->GetAxis(0)->SetRangeUser(0.1, 10.1);
    fkkUSDist->GetAxis(0)->SetRangeUser(0.1, 10.1);
    fkkUSDist->GetAxis(1)->SetRangeUser(0.98, 1.1);
    fTrigDist->GetAxis(0)->SetRangeUser(0.1, 10.1);

    for(int i=0; i<4; i++){
        if(i==0){
            cLSDist->cd(i+1);
            fkkLSDist->Projection(i)->Draw("SAME");
            cUSDist->cd(i+1);
            fkkUSDist->Projection(i)->Draw("SAME");
            cTrigDist->cd(i+1);
            fTrigDist->Projection(i)->Draw("SAME");
        }else if(i==1){
            cLSDist->cd(i+1);
            hLSDist[i-1]->Draw("SAME");
            cUSDist->cd(i+1);
            hUSDist[i-1]->Draw("SAME");            
        }else{
            cLSDist->cd(i+1);
            hLSDist[i-1]->Draw("SAME");
            cUSDist->cd(i+1);
            hUSDist[i-1]->Draw("SAME");
            cTrigDist->cd(i+1);
            hTrigDist[i-2]->Draw("SAME");
        }
    }

    // Plotting the US and LS invariant mass per pT bin on the same plot
    TCanvas *c0 = new TCanvas("cTest", "cTest", 50, 50, 800, 800);
/*    c0->Divide(4,4);

    for(int j = 0; j< 15; j++){
        c0->cd(j+1);
        phiInvMassBinned[j]->Draw("E");
        likeSignInvMassBinned[j]->Draw("SAME");
    }
*/  
    c0->cd();
    phiInvMassBinned[0]->Draw("HIST E");
    likeSignInvMassBinned[0]->Draw("SAME");

    // Plotting the 'corrected' invariant mass per pT bin
    TCanvas *c1 = new TCanvas("cCorrInvMass", "cCorrInvMass", 50, 50, 800, 800);
    c1->Divide(4,4);

    for(int k =0; k < 15; k++){
        c1->cd(k+1);
        corrInvMass[k]->Draw();
        fits[k]->Draw("SAME");
    }
/*
    c1->cd();
    corrInvMass[0]->Draw();
    fits[0]->Draw("SAME");    
*/
   
    dphiHPhi->GetAxis(0)->SetRangeUser(4.0, 8.0); //cutting on trigger pt between 3 GeV/c and 6 GeV/c
    dphiHPhi->GetAxis(1)->SetRangeUser(2.0,4.0); //cutting on phi pt between 1 GeV/c and 3 GeV/c
  
    //performing same cuts for KK pairs:
    dphiHKK->GetAxis(0)->SetRangeUser(4.0, 8.0);
    dphiHKK->GetAxis(1)->SetRangeUser(2.0,4.0);
 
    //3D delta eta, delta phi, inv mass dist (all candidates)
    TH3D *dEtadPhiDist = dphiHPhi->Projection(2, 3, 4);
    dEtadPhiDist->Rebin3D(4,2,1);
    TH3D *dEtadPhiLSDist = dphiHKK->Projection(2, 3, 4);
    dEtadPhiLSDist->Rebin3D(4,2,1);

    dEtadPhiDist->GetZaxis()->SetRangeUser(1.04, 1.06);
    dEtadPhiLSDist->GetZaxis()->SetRangeUser(1.04, 1.06);

    TH1D *scaleProj = dEtadPhiDist->Project3D("x");
    TH1D *LSscaleProj = dEtadPhiLSDist->Project3D("x");
    Double_t testScale = (scaleProj->Integral())/(LSscaleProj->Integral());

    dEtadPhiDist->GetZaxis()->SetRange(0,0);
    dEtadPhiLSDist->GetZaxis()->SetRange(0,0);

    dEtadPhiLSDist->Scale(testScale);
    dEtadPhiDist->Add(dEtadPhiLSDist, -1.0);

    dEtadPhiDist->GetZaxis()->SetRangeUser(1.01, 1.03);
    
    TH2D *twoCorr = dEtadPhiDist->Project3D("xy");

    twoCorr->GetYaxis()->SetRange(7,9);
    TH1D *dEta = twoCorr->ProjectionX();
    twoCorr->GetYaxis()->SetRange(0,0);
    TF1 *dEtaFit = new TF1("dEtaFit", "[0] - [1]*TMath::Abs(x)", -1.6, 1.6);
    dEtaFit->SetParameters(1000.0, 500.0);
   
    TCanvas *dEtaCanvas = new TCanvas("cdeta", "cdeta", 50, 50, 800, 800);
    dEtaCanvas->cd();
    dEta->Fit("dEtaFit", "R");
    dEta->Draw();

    TH2D *correctedTwoCorr = twoCorr->Clone("correctedTwoCorr");
    for(int iphi = 1; iphi < 16+1; iphi++){
        for(int ieta = 4; ieta < 10; ieta++){
            Double_t value = correctedTwoCorr->GetBinContent(ieta, iphi);
            value = 2.0*value/(dEtaFit->Eval(-2.75+(0.5*(ieta-1))));
            correctedTwoCorr->SetBinContent(ieta, iphi, value);
        }
    }
    TCanvas *twocanvas = new TCanvas("2dcanvas", "2dcanvas", 50,50,800,800);
    twocanvas->cd();
    twoCorr->GetXaxis()->SetRangeUser(-1.5, 1.5);
    twoCorr->Draw("SURF1");

    TCanvas *twoCorrcanvas = new TCanvas("2dCorrcanvas", "2dCorrcanvas", 50,50,800,800);
    twoCorrcanvas->cd();
    correctedTwoCorr->GetXaxis()->SetRangeUser(-1.5, 1.5);
    correctedTwoCorr->Draw("SURF1");



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
        numBins = (Double_t)phiInvMassPerDPhi[i]->GetNbinsX();
        sideband = phiInvMassPerDPhi[i]->Integral((int)((1.04-invMin)*numBins/(invMax-invMin)), (int)((1.06-invMin)*numBins/(invMax-invMin)));
        likeSignSideBand = likesignInvMassPerDPhi[i]->Integral((int)((1.04-invMin)*numBins/(invMax-invMin)), (int)((1.06-invMin)*numBins/(invMax-invMin)));

        phiInvMassPerDPhi[i]->SetTitle(Form("Inv Mass dist. for %.2f < #Delta#varphi < %.2f", (-1.57 + i*(0.3925)), (-1.57+(i+1)*(0.3925))));
        phiInvMassPerDPhi[i]->Draw("HIST E");

        scaleFactor = sideband/likeSignSideBand;
        likesignInvMassPerDPhi[i]->Scale(scaleFactor);
        likesignInvMassPerDPhi[i]->Sumw2();
        likesignInvMassPerDPhi[i]->SetLineColor(2);
        likesignInvMassPerDPhi[i]->Draw("SAME HIST E");
        
        //For corrected plots, need to rebin the distributions to get rid of some of the noise
        phiInvMassPerDPhi[i]->Rebin();
        phiInvMassPerDPhi[i]->Sumw2();
        likesignInvMassPerDPhi[i]->Rebin();
        likesignInvMassPerDPhi[i]->Sumw2();
        corrPhiInvMassPerDPhi[i] = (TH1D*)phiInvMassPerDPhi[i]->Clone();
        corrPhiInvMassPerDPhi[i]->Add(likesignInvMassPerDPhi[i], -1.0);
    }

    TCanvas *cCorrInvDPhi = new TCanvas("cCorrInvDPhi", "cCorrInvDPhi", 40, 40, 800, 800);
    cCorrInvDPhi->Divide(4,4);
    TF1 *dPhiFits[16];
    Double_t numPhi[16];
    TH1D* corrDPhi = new TH1D("corrDPhi", "BG-Corrected #Delta#varphi for Hadron-#Phi(1020)", 16, -1.57, 4.71);
    for(int i = 0; i<16; i++){
        cCorrInvDPhi->cd(i+1);
        //dPhiFits[i] = new TF1(Form("dphif%i", i), "gaus(0) + pol0(3)",0.98,1.1); 
        dPhiFits[i] = new TF1(Form("dphif%i", i), "[0]*TMath::Voigt(x - [1], [2], [3], 3) + pol0(4)",0.98, 1.1); 
        dPhiFits[i]->SetParameter(1, 1.020);
        dPhiFits[i]->SetParameter(2, 0.002);
        dPhiFits[i]->SetParameter(0, 1);
        dPhiFits[i]->FixParameter(3, 0.00426);
        dPhiFits[i]->SetParLimits(1, 1.015, 1.025);
        //dPhiFits[i]->SetParLimits(2, 0.0010, 0.007);
        //dPhiFits[i]->SetParLimits(3, 0.0010, 0.010);
        //dPhiFits[i]->SetParLimits(0, 100, 10000);
        //dPhiFits[i]->SetParLimits(3, -200, 200);
        corrPhiInvMassPerDPhi[i]->Fit(Form("dphif%i", i), "R");
        corrPhiInvMassPerDPhi[i]->Draw("SAME HIST E FUNC");
        //printf("%i. Total: %f\n", i, dPhiFits[i]->Integral(1.005, 1.035));
        numBins = (Double_t)corrPhiInvMassPerDPhi[i]->GetNbinsX();
        printf("correction: %f", dPhiFits[i]->GetParameter(3)*(int)((1.03-1.01)*numBins/(invMax-invMin)));
        numPhi[i] = corrPhiInvMassPerDPhi[i]->Integral((int)((1.01-invMin)*numBins/(invMax-invMin)), (int)((1.03-invMin)*numBins/(invMax-invMin))) - dPhiFits[i]->GetParameter(4)*(int)((1.03-1.01)*numBins/(invMax-invMin));
        printf("Corrected num Phis: %f\n", numPhi[i]);
        corrDPhi->SetBinContent(i+1, numPhi[i]);
    }

    TCanvas *cCorrDPhi = new TCanvas("cCorrDPhi", "cCorrDPhi", 30, 30, 600, 600);
    cCorrDPhi->cd();
    TF1 *corrFit = new TF1("corrFit", "gaus(0) + gaus(3) + pol0(6)", -1.4, 4.6);
    corrFit->SetParameter(6, 500);
    corrFit->SetParameter(0, 100);
    corrFit->SetParameter(1, 0.0);
    corrFit->SetParameter(2, 1.0);
    corrFit->SetParameter(3, 50);
    corrFit->SetParameter(4, 3.14);
    corrFit->SetParLimits(4, 3.0, 3.25);
    corrFit->SetParameter(5, 1.5);

    corrDPhi->Fit("corrFit", "R");
    corrDPhi->Draw();
    //corrFit->Draw("SAME");

    //Reset the Delta-phi axis range after the dphi binned projections are done being created above.
    dphiHPhi->GetAxis(2)->SetRange(0,0);
    dphiHKK->GetAxis(2)->SetRange(0,0);

    TH1D *dphiPhiPeak;
    dphiHPhi->GetAxis(4)->SetRangeUser(1.010, 1.030);
    dphiPhiPeak = dphiHPhi->Projection(2);
    dphiPhiPeak->SetTitle("#Delta#varphi for Hadron-#Phi(1020) in Peak"); 

    TH1D *dphiPhiSideband;
    dphiHPhi->GetAxis(4)->SetRangeUser(1.040, 1.060);
    dphiPhiSideband = dphiHPhi->Projection(2);
    dphiPhiSideband->SetTitle("#Delta#varphi for Hadron-#Phi(1020) in Sideband");
    
// Setting & plotting Delta Phi distribution for Unlike sign and Likesign pairs in the
// sideband region (scaled by the ratio of the integral of the sideband region)
    dphiHPhi->GetAxis(4)->SetRangeUser(1.040,1.060);
    TH1D* dphiUSSideband = dphiHPhi->Projection(2);
    dphiHKK->GetAxis(4)->SetRangeUser(1.040,1.060);
    TH1D* dphiLSSideband = dphiHKK->Projection(2);
    dphiHKK->GetAxis(4)->SetRangeUser(0.98, 1.1);
    dphiHPhi->GetAxis(4)->SetRangeUser(0.98, 1.1);
    numBins = (Double_t)dphiHPhi->GetAxis(4)->GetNbins();
    sideband = dphiHPhi->Projection(4)->Integral((int)((1.040-invMin)*numBins/(invMax-invMin)), (int)((1.060-invMin)*numBins/(invMax-invMin)));
    likeSignSideBand = dphiHKK->Projection(4)->Integral((int)((1.040-invMin)*numBins/(invMax-invMin)), (int)((1.060-invMin)*numBins/(invMax-invMin)));
    scaleFactor = sideband/likeSignSideBand;
    printf("sideband: %f, likeSignSideBand: %f, scaleFactor: %f \n", sideband, likeSignSideBand, scaleFactor);
    dphiLSSideband->Scale(scaleFactor);

    TCanvas *cDphiHPhi3 = new TCanvas("cDphiHPhi3", "cDphiHPhi3", 50, 50, 600, 600);
    cDphiHPhi3->cd();
    dphiUSSideband->SetLineWidth(2);
    dphiUSSideband->SetTitle("#Delta#varphi for Hadron-#Phi(1020) in Sideband Region");
    dphiUSSideband->GetXaxis()->SetTitle("#Delta#varphi");
    dphiUSSideband->Draw("SAME HIST E");
    dphiLSSideband->SetLineColor(2);
    dphiLSSideband->Draw("SAME HIST E");


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
    dphiUSPeak->SetTitle("#Delta#varphi for Hadron-#Phi(1020) in Peak Region");
    dphiUSPeak->GetXaxis()->SetTitle("#Delta#varphi");
    dphiUSPeak->Draw("SAME HIST E");
    dphiLSPeak->SetLineColor(2);
    dphiLSPeak->Draw("SAME HIST E");

//Calculating and plotting the difference between US and LS pairs for the phi mass peak region (second method)
    dphiUSPeak->Sumw2();
    dphiLSPeak->Sumw2();
    TH1D *corrDPhiPeakRegion = (TH1D*)dphiUSPeak->Clone();
    corrDPhiPeakRegion->Add(dphiLSPeak, -1.0);
    

    Double_t xbins[18]= {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};

    double low = corrDPhiPeakRegion->GetXaxis()->GetBinLowEdge(1);
    printf("low %i: %f\n", 0, low);
    xbins[0] = low;
    for(int i=1;i<17; i++){
       low  = corrDPhiPeakRegion->GetXaxis()->GetBinLowEdge(3+4*(i-1));
       printf("low %i: %f\n", i, low);
       xbins[i] = low;
    }
    xbins[17] = 4.71; 

    corrDPhiPeakRegion->Rebin(17, "hnew", xbins);
    
    TCanvas *cDphiHPhi4 = new TCanvas("cDphiHPhi4", "cDphiHPhi4", 50, 50, 600, 600);
    cDphiHPhi4->cd();

    TF1* dphifit = new TF1("dphifit", "gaus(0) + gaus(3) + pol0(6)", -1.3, 4.5);
    dphifit->SetParameter(6, hnew->GetBinContent(9));
    dphifit->SetParameter(0, 200);
    dphifit->SetParLimits(0, 100.0, 10000.0); 
    dphifit->SetParameter(1, 0.0);
    dphifit->SetParameter(2, 1.0);
    dphifit->SetParameter(3, 100);
    dphifit->SetParLimits(3, 1.0, 10000.0);
    dphifit->SetParameter(4, 3.14);
    dphifit->SetParLimits(4, 3.1, 3.2);
    dphifit->SetParameter(5, 1.5);

    hnew->Fit("dphifit", "R");


    hnew->SetTitle("Corrected #Delta#varphi for Hadron-#Phi(1020) in Peak Region");
    hnew->GetXaxis()->SetTitle("#Delta#varphi");
    hnew->GetXaxis()->SetRangeUser(-1.3, 4.4);
    
    hnew->Draw("SAME HIST E");


}


