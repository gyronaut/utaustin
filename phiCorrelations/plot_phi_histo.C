#include "TLegend.h"
#include <string>
#include <sstream>

/* function to plot the Raw PT spectrum of US and LS Kaon pairs, as well as the
 * (LS BG) corrected PT spectrum of Phi mesons
 */
void plotRawAndCorrectedPT(THnSparseF *kkUSDist, THnSparseF *kkLSDist){
    Double_t error;
    Double_t sideband = 0.0;
    Double_t likeSignSideBand = 0.0;
    Double_t scaleFactor = 0.0;
    
    Int_t numInvMassBins = kkUSDist->GetAxis(1)->GetNbins();
    Double_t invMin = kkUSDist->GetAxis(1)->GetBinLowEdge(1);
    Double_t invMax = kkUSDist->GetAxis(1)->GetBinUpEdge(numInvMassBins);
    Double_t binsPerMass = (Double_t)(numInvMassBins)/(invMax-invMin);
    Double_t peakBounds[2] = {1.01, 1.03};
    Double_t sidebandBounds[2] = {1.05, 1.07};

    //Int_t numPTBins = (Int_t)((sizeof(ptbounds)/sizeof(ptbounds[0]))); 
    //printf("sizeof(ptbounds): %i, sizeof(ptbounds[0]): %i, numPTBins: %i\n", sizeof(ptbounds), sizeof(ptbounds[0]), numPTBins);
    const Int_t numPTBins = 17;
    Double_t ptbounds[] = {0.4, 0.6, 0.8, 1.0, 1.2, 1.4, 1.6, 1.8, 2.0, 2.5, 3.0, 4.0, 5.0, 6.0, 7.0, 10.0, 2.0, 4.0};
    Int_t ptbins[] = {4, 6, 8, 10, 12, 14, 16, 18, 20, 25, 30, 40, 50, 60, 70, 100, 20, 40};

    TH1D *phiPTSpectrum = new TH1D("phiPTspectrum", "p_{T}^{#Phi} Spectrum Corrected with LS BG", 15, ptbounds);
    TH1D *corrInvMass[numPTBins];
    TH1D *phiInvMassBinned[numPTBins];
    TH1D *likeSignInvMassBinned[numPTBins];

    TF1 *fits[numPTBins];
    TF1 *bgFits[numPTBins];
    for(int n = 0; n < numPTBins; n++){
        fits[n] = new TF1(Form("f%i", n),  "[0]*TMath::Voigt(x - [1], [2], [3], 4) + pol2(4)",0.99, 1.1);
        fits[n]->SetParameter(1, 1.020);
        fits[n]->SetParameter(2, 0.0002);
        fits[n]->SetParameter(0, 600);
        fits[n]->FixParameter(3, 0.00426);
        fits[n]->SetParLimits(1, 1.010, 1.030);
        //fits[n]->SetParLimits(0, 225, 400);
        //fits[n]->SetParLimits(2, 0.0010, 0.007);
        //fits[n]->SetParLimits(3, 0.0010, 0.010);
        bgFits[n] = new TF1(Form("bgfit%i", n), "pol2(0)", 0.99, 1.1);

    } 

    Double_t numPhiPT[15];
    for(int i =0; i<numPTBins; i++){
        //Set up the InvMass projections for each PT Bin, for US and LS distributions
        kkUSDist->GetAxis(0)->SetRange(ptbins[i], ptbins[i+1]);
        phiInvMassBinned[i] = (TH1D*)kkUSDist->Projection(1);
        phiInvMassBinned[i]->SetTitle(Form("%.1f < p_{T} < %.1f GeV/c", ptbounds[i], ptbounds[i+1]));
        phiInvMassBinned[i]->Sumw2();
        kkLSDist->GetAxis(0)->SetRange(ptbins[i], ptbins[i+1]);
        likeSignInvMassBinned[i] = (TH1D*)kkLSDist->Projection(1);
        likeSignInvMassBinned[i]->SetTitle(Form("%.1f < p_{T} < %.1f GeV/c", ptbounds[i], ptbounds[i+1]));
        likeSignInvMassBinned[i]->Sumw2();

        //Integrate both US and LS distributions over the sideband region to calculate a scalefactor.
        sideband = phiInvMassBinned[i]->Integral((int)(binsPerMass*(sidebandBounds[0]-invMin)),(int)(binsPerMass*(sidebandBounds[1]-invMin)));
        likeSignSideBand = likeSignInvMassBinned[i]->Integral((int)(binsPerMass*(sidebandBounds[0]-invMin)),(int)(binsPerMass*(sidebandBounds[1]-invMin)));
        scaleFactor = sideband/likeSignSideBand;

        //Create the "corrected" inv mass distribution by subtracting the scaled LS distribution from 
        //the US distribution, then fit with a Voigt + 0th order polynomial
        corrInvMass[i] = (TH1D*)phiInvMassBinned[i]->Clone(Form("corrInvMass_%i_%i", ptbins[i], ptbins[i]));
        likeSignInvMassBinned[i]->Scale(scaleFactor);
        corrInvMass[i]->Add(likeSignInvMassBinned[i], -1.0);
        if(i<15){
            corrInvMass[i]->Fit(Form("f%i", i), "R");
        }else{
            corrInvMass[i]->Fit(Form("f%i", i), "N R");
        }
        //Get the number of phi mesons in each pT bin by integrating the corrected distribution in the
        //mass peak range, subtracting the 2nd order polynomial contribution in that range, and 
        //dividing by the bin width
        bgFits[i]->SetParameters(fits[i]->GetParameter(4), fits[i]->GetParameter(5), fits[i]->GetParameter(6));

        if(i < 15){
            numPhiPT[i] = (corrInvMass[i]->IntegralAndError((Int_t)(binsPerMass*(peakBounds[0]-invMin)),(Int_t)(binsPerMass*(peakBounds[1]-invMin)), error) - bgFits[i]->Integral(peakBounds[0], peakBounds[1]));

            printf("lowInvMassBin: %i, lowInvMassBinCenter: %.4f, upBin: %i, upBinCenter: %.4f\n", (Int_t)(binsPerMass*(peakBounds[0]-invMin)), corrInvMass[i]->GetBinCenter((Int_t)(binsPerMass*(peakBounds[0]-invMin))), (Int_t)(binsPerMass*(peakBounds[1]-invMin)), corrInvMass[i]->GetBinCenter((Int_t)(binsPerMass*(peakBounds[1]-invMin))));
            printf("ptbin: (%.1f, %.1f), binwidth: %f, numphiPT: %i\n", ptbounds[i], ptbounds[i+1], phiPTSpectrum->GetBinWidth(i+1), numPhiPT[i]);
            phiPTSpectrum->SetBinContent(i+1, (double)(numPhiPT[i])/(phiPTSpectrum->GetBinWidth(i+1)));
            phiPTSpectrum->SetBinError(i+1, (TMath::Sqrt((double)numPhiPT[i]))/phiPTSpectrum->GetBinWidth(i+1));
        }
        likeSignInvMassBinned[i]->SetLineColor(2); 
    }

    // Plotting the US and LS invariant mass per pT bin on the same plot
    TCanvas *cInvMassPTBins = new TCanvas("cInvMassPtBins", "cInvMassPtBins", 50, 50, 800, 800);
    cInvMassPTBins->Divide(4,4);

    for(int j = 1; j< 15; j++){
        cInvMassPTBins->cd(j+1);
        phiInvMassBinned[j]->Draw("E");
        likeSignInvMassBinned[j]->Draw("SAME");
    }

    TCanvas *cSinglePTBin = new TCanvas("cSinglePTBin", "cSinglePTBin", 50, 50, 800, 800);
    cSinglePTBin->cd();
    phiInvMassBinned[16]->GetXaxis()->SetRangeUser(0.99, 1.08);
    phiInvMassBinned[16]->GetXaxis()->SetTitle("Inv Mass (GeV/c^{2})");
    phiInvMassBinned[16]->SetLineWidth(2);
    phiInvMassBinned[16]->SetTitle("Invariant Mass for Kaon Pairs");
    phiInvMassBinned[16]->Draw("E H");
    likeSignInvMassBinned[16]->SetLineWidth(2);
    likeSignInvMassBinned[16]->Draw("SAME E H");

    // Plotting the 'corrected' invariant mass per pT bin
    TCanvas *cCorrInvMassPtBins = new TCanvas("cCorrInvMassPtBins", "cCorrInvMassPtBins", 50, 50, 800, 800);
    cCorrInvMassPtBins->Divide(4,4);

    TPaveText *corrText;
    for(int k =1; k < 15; k++){
        corrText = new TPaveText(0.5, 0.68, 0.9, 0.85, "NDC");
        corrText->AddText(Form("Num #phi: %i", (Int_t)numPhiPT[k]));
        cCorrInvMassPtBins->cd(k+1);
        corrInvMass[k]->Draw("E H");
        fits[k]->Draw("SAME");
        bgFits[k]->SetLineColor(2);
        bgFits[k]->Draw("SAME");
        corrText->Draw();
    }

    TCanvas *cCorrSinglePT = new TCanvas("cCorrSinglePT", "cCorrSinglePT", 50, 50, 800, 800);
    cCorrSinglePT->cd();
    corrInvMass[16]->GetXaxis()->SetRangeUser(0.99, 1.08);
    corrInvMass[16]->GetXaxis()->SetTitle("Inv Mass (GeV/c^{2})");
    corrInvMass[16]->SetTitle("Background Corrected Invariant Mass");
    corrInvMass[16]->SetLineWidth(2);
    corrInvMass[16]->Draw("E H");
    fits[16]->SetLineColor(4);
    fits[16]->SetLineStyle(2);
    fits[16]->SetLineWidth(4);
    fits[16]->Draw("SAME");
    bgFits[16]->SetLineColor(8);
    //bgFits[16]->Draw("SAME");

    //Plot the Like-sign background corrected phi(1020) pT distribution:
    TCanvas *cPhiPt = new TCanvas("cPhiPt", "cPhiPt", 50, 50, 800, 800);
    cPhiPt->cd();
    gPad->SetLogy();
    phiPTSpectrum->GetXaxis()->SetTitle("p_{T} (GeV/c^{2})");
    phiPTSpectrum->GetYaxis()->SetTitle("Counts/(bin width)");
    phiPTSpectrum->SetMarkerStyle(3);
    phiPTSpectrum->SetMarkerSize(2.0);
    phiPTSpectrum->SetMarkerColor(4);
    phiPTSpectrum->Draw("C");

    printf("numphi: %i, numphiPTRange: %i\n\n\n",phiPTSpectrum->Integral("width"), phiPTSpectrum->Integral(4,10, "width")); 

    TPaveText *ptNumbers = new TPaveText(0.394, 0.268, 0.883, 0.380, "NDC NB");
    ptNumbers->SetFillColor(0);
    Int_t totPhi = (Int_t)phiPTSpectrum->Integral("width");
    Int_t ptRangePhi = (Int_t)phiPTSpectrum->Integral(9,11, "width");
    ptNumbers->AddText(Form("Total #phi(1020): %i", totPhi));
    ptNumbers->AddText(Form("In range (2 < p_{T}^{#phi} < 4) GeV/c^{2}: %i", ptRangePhi));
    ptNumbers->Draw();
}

void plotAllDistributions(THnSparse *kkLSDist, THnSparse *kkUSDist, THnSparse *trigDist){
    //first set the momentum range for the Assoc and Trigger particles, and use this range
    //to find the phi, eta (, and inv mass) distribution for trigger (assoc) particles
    kkLSDist->GetAxis(0)->SetRangeUser(2.0, 4.0);
    kkUSDist->GetAxis(0)->SetRangeUser(2.0, 4.0);
    trigDist->GetAxis(0)->SetRangeUser(4.0, 8.0);

    TH1D* hLSDist[3];
    TH1D* hUSDist[3];
    TH1D* hTrigDist[2];
    for(int i=0; i<3; i++){
        hLSDist[i] = kkLSDist->Projection(i+1);
        hUSDist[i] = kkUSDist->Projection(i+1);
        if(i<2){
            hTrigDist[i] = trigDist->Projection(i+1);
        }
    }

    //Next, reset the momentum ranges so that the full momentum distribution for trigger hadrons
    //and assoc KK (US and LS) can be plotted
    kkLSDist->GetAxis(0)->SetRangeUser(0.1, 10.1);
    kkLSDist->GetAxis(1)->SetRangeUser(0.98, 1.1);
    kkUSDist->GetAxis(0)->SetRangeUser(0.1, 10.1);
    kkUSDist->GetAxis(1)->SetRangeUser(0.98, 1.1);
    trigDist->GetAxis(0)->SetRangeUser(0.1, 10.1);

    TCanvas *cLSDist = new TCanvas("cLSDist", "cLSDist", 50, 50, 800, 800);
    cLSDist->Divide(2,2);
    TCanvas *cUSDist = new TCanvas("cUSDist", "cUSDist", 50, 50, 800, 800);
    cUSDist->Divide(2,2);
    TCanvas *cTrigDist = new TCanvas("cTrigDist", "cTrigDist", 50, 50, 800, 800);
    cTrigDist->Divide(2,2);

    for(int i=0; i<4; i++){
        if(i==0){
            cLSDist->cd(i+1);
            kkLSDist->Projection(i)->Draw("SAME");
            cUSDist->cd(i+1);
            kkUSDist->Projection(i)->Draw("SAME");
            cTrigDist->cd(i+1);
            trigDist->Projection(i)->Draw("SAME");
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
}

void plot2DCorrelations(TH3D* dEtadPhiDist, TH3D* dEtadPhiLSDist, TH3D* dEtadPhiMixed, TH3D* dEtadPhiLSMixed, TString suffix){

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
    
    TH2D *twoCorr = dEtadPhiDist->Project3D("xy")->Clone(Form("twoCorr_%s", suffix.Data()));
    twoCorr->GetXaxis()->SetTitle("#Delta#eta");
    twoCorr->GetYaxis()->SetTitle("#Delta#varphi");
    twoCorr->GetZaxis()->SetTitle("");
    twoCorr->SetTitle("Hadron-(K+K-) Correlations");

    twoCorr->GetYaxis()->SetRange(7,9);
    TH1D *dEta = twoCorr->ProjectionX();
    twoCorr->GetYaxis()->SetRange(0,0);

    TCanvas *twocanvas = new TCanvas(Form("2dcanvas_%s", suffix.Data()), Form("2dcanvas_%s", suffix.Data()), 50,50,800,800);
    twocanvas->cd();
    twoCorr->GetXaxis()->SetRangeUser(-1.5, 1.5);
    twoCorr->Draw("SURF1");

    /*
    TF1 *dEtaFit = new TF1(Form("dEtaFit_%s", suffix.Data()), "[0] - [1]*TMath::Abs(x)", -1.6, 1.6);
    dEtaFit->SetParameters(1000.0, 500.0);
   
    TCanvas *dEtaCanvas = new TCanvas(Form("cdeta_%s", suffix.Data()), Form("cdeta_%s", suffix.Data()), 50, 50, 800, 800);
    dEtaCanvas->cd();
    dEta->Fit("dEtaFit", "R");
    dEta->Draw();

    TH2D *correctedTwoCorr = twoCorr->Clone(Form("correctedTwoCorr_%s", suffix.Data()));
    for(int iphi = 1; iphi < 16+1; iphi++){
        for(int ieta = 6; ieta < 19; ieta++){
            Double_t value = correctedTwoCorr->GetBinContent(ieta, iphi);
            value = 2.0*value/(dEtaFit->Eval(-2.875+(0.25*(ieta-1))));
            correctedTwoCorr->SetBinContent(ieta, iphi, value);
        }
    }

    TCanvas *twoCorrcanvas = new TCanvas(Form("2dCorrcanvas_%s", suffix.Data()), Form("2dCorrcanvas_%s", suffix.Data()), 50,50,800,800);
    twoCorrcanvas->cd();
    correctedTwoCorr->GetXaxis()->SetRangeUser(-1.5, 1.5);
    correctedTwoCorr->Draw("SURF1");
*/
}

TH1D* plotPhiCorrelationsV1(THnSparse *dphiHPhi, THnSparse *dphiHKK){

    Int_t numBins = dphiHPhi->GetAxis(4)->GetNbins();
    Double_t invMin = dphiHPhi->GetAxis(4)->GetBinLowEdge(1);
    Double_t invMax = dphiHPhi->GetAxis(4)->GetBinUpEdge(numBins);
    Double_t binsPerMass = (Double_t)(numBins)/(invMax-invMin);

    TH1D *phiInvMassPerDPhi[16];
    TH1D *likesignInvMassPerDPhi[16]; 
    TH1D *corrPhiInvMassPerDPhi[16];
    TCanvas *cInvDPhi = new TCanvas("cInvDPhi", "cIncDPhi", 50,50,800,800);
    cInvDPhi->Divide(4,4);
    for(int i = 0; i < 16; i++){

        Double_t scaleFactor=0.0;
        cInvDPhi->cd(i+1);
        dphiHPhi->GetAxis(2)->SetRange((4*i)+1, 4*(i+1));
        dphiHKK->GetAxis(2)->SetRange((4*i)+1, 4*(i+1));
        
        phiInvMassPerDPhi[i] = dphiHPhi->Projection(4);
        phiInvMassPerDPhi[i]->Sumw2();
        likesignInvMassPerDPhi[i] = dphiHKK->Projection(4);
        likesignInvMassPerDPhi[i]->Sumw2();

        Double_t sideband = phiInvMassPerDPhi[i]->Integral((int)((1.04-invMin)*binsPerMass), (int)((1.06-invMin)*binsPerMass));
        Double_t likeSignSideBand = likesignInvMassPerDPhi[i]->Integral((int)((1.04-invMin)*binsPerMass), (int)((1.06-invMin)*binsPerMass));

        phiInvMassPerDPhi[i]->SetTitle(Form("Inv Mass dist. for %.2f < #Delta#varphi < %.2f", (-1.57 + i*(0.3925)), (-1.57+(i+1)*(0.3925))));
        phiInvMassPerDPhi[i]->Draw("HIST E");

        scaleFactor = sideband/likeSignSideBand;
        likesignInvMassPerDPhi[i]->Scale(scaleFactor);
        likesignInvMassPerDPhi[i]->SetLineColor(2);
        likesignInvMassPerDPhi[i]->Draw("SAME HIST E");
        
        //For corrected plots, need to rebin the distributions to get rid of some of the noise
        phiInvMassPerDPhi[i]->Rebin();
        likesignInvMassPerDPhi[i]->Rebin();
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
        corrPhiInvMassPerDPhi[i]->Fit(Form("dphif%i", i), "R");
        corrPhiInvMassPerDPhi[i]->Draw("SAME HIST E FUNC");
        numBins = (Double_t)corrPhiInvMassPerDPhi[i]->GetNbinsX();
        printf("correction: %f", dPhiFits[i]->GetParameter(3)*(int)((1.03-1.01)*numBins/(invMax-invMin)));
        numPhi[i] = corrPhiInvMassPerDPhi[i]->Integral((int)((1.01-invMin)*numBins/(invMax-invMin)), (int)((1.03-invMin)*numBins/(invMax-invMin))) - dPhiFits[i]->GetParameter(4)*(int)((1.03-1.01)*numBins/(invMax-invMin));
        printf("Corrected num Phis: %f\n", numPhi[i]);
        corrDPhi->SetBinContent(i+1, numPhi[i]);
    }

    TCanvas *cCorrDPhi = new TCanvas("cCorrDPhi", "cCorrDPhi", 30, 30, 600, 600);
    cCorrDPhi->cd();
    TF1 *corrFit = new TF1("corrFit", "gaus(0) + gaus(3) + pol0(6)", -1.4, 4.6);
    corrFit->SetParameter(6, corrDPhi->GetBinContent(9));
    corrFit->SetParameter(0, 1000);
    corrFit->SetParLimits(0, 1.0, 10000);
    corrFit->SetParameter(1, 0.0);
    corrFit->SetParameter(2, 1.0);
    corrFit->SetParameter(3, 500);
    corrFit->SetParLimits(3, 1.0, 10000);
    corrFit->SetParameter(4, 3.14);
    corrFit->SetParLimits(4, 3.0, 3.25);
    corrFit->SetParameter(5, 1.5);

    corrFit->SetLineColor(4);
    corrFit->SetLineWidth(3);
    corrFit->SetLineStyle(7);

    corrDPhi->Fit("corrFit", "R");
    corrDPhi->SetLineWidth(2);
    //corrDPhi->GetYaxis()->SetRangeUser(0, 1200);
    corrDPhi->GetXaxis()->SetTitle("#Delta#varphi");
    corrDPhi->Draw("H E");

    return corrDPhi;
}

TH1D* plotPhiCorrelationsV2(THnSparse *dphiHPhi, THnSparse *dphiHKK){

    Int_t numBins = dphiHPhi->GetAxis(4)->GetNbins();
    Double_t invMin = dphiHPhi->GetAxis(4)->GetBinLowEdge(1);
    Double_t invMax = dphiHPhi->GetAxis(4)->GetBinUpEdge(numBins);
    Double_t binsPerMass = (Double_t)(numBins)/(invMax-invMin);

// Setting & plotting Delta Phi distribution for Unlike sign and Likesign pairs in the
// sideband region (scaled by the ratio of the integral of the sideband region)
    dphiHPhi->GetAxis(4)->SetRangeUser(1.040,1.060);
    TH1D* dphiUSSideband = dphiHPhi->Projection(2);
    dphiUSSideband->SetName("dphiUSSideband");
    dphiUSSideband->Sumw2();
    dphiHKK->GetAxis(4)->SetRangeUser(1.040,1.060);
    TH1D* dphiLSSideband = dphiHKK->Projection(2);
    dphiLSSideband->SetName("dphiLSSideband");
    dphiLSSideband->Sumw2();
    dphiHKK->GetAxis(4)->SetRangeUser(0.98, 1.1);
    dphiHPhi->GetAxis(4)->SetRangeUser(0.98, 1.1);
    Double_t sideband = dphiHPhi->Projection(4)->Integral((int)((1.040-invMin)*binsPerMass), (int)((1.060-invMin)*binsPerMass));
    Double_t likeSignSideBand = dphiHKK->Projection(4)->Integral((int)((1.040-invMin)*binsPerMass), (int)((1.060-invMin)*binsPerMass));
    Double_t scaleFactor = sideband/likeSignSideBand;
    printf("sideband: %f, likeSignSideBand: %f, scaleFactor: %f \n", sideband, likeSignSideBand, scaleFactor);
    dphiLSSideband->Scale(scaleFactor);

    TCanvas *cDphiHPhi3 = new TCanvas("cDphiHPhi3", "cDphiHPhi3", 50, 50, 600, 600);
    cDphiHPhi3->cd();
    dphiUSSideband->SetLineWidth(2);
    dphiUSSideband->SetTitle("#Delta#varphi for Hadron-#Phi(1020) in Right Sideband Region");
    dphiUSSideband->GetXaxis()->SetTitle("#Delta#varphi");
    dphiUSSideband->Rebin(4);
    dphiLSSideband->Rebin(4);
    dphiUSSideband->Draw("SAME HIST E");
    dphiLSSideband->SetLineColor(2);
    dphiLSSideband->Draw("SAME HIST E");

    TPaveText *sidebandText = new TPaveText(0.537, 0.716, 0.852, 0.818, "NDC");
    sidebandText->AddText("1.040 < m_{kk} < 1.060 GeV/c^{2}"); 
    sidebandText->AddText(Form("LS scale factor: %f", scaleFactor));
    sidebandText->SetFillColor(0);
    sidebandText->SetBorderSize(1);
    sidebandText->Draw("SAME");

    dphiHPhi->GetAxis(4)->SetRangeUser(0.994,1.006);
    TH1D* dphiUSLSideband = dphiHPhi->Projection(2);
    dphiUSLSideband->SetName("dphiUSLSideband");
    dphiUSLSideband->Sumw2();
    dphiHKK->GetAxis(4)->SetRangeUser(0.994,1.006);
    TH1D* dphiLSLSideband = dphiHKK->Projection(2);
    dphiLSLSideband->SetName("dphiLSLSideband");
    dphiLSLSideband->Sumw2();
    dphiHKK->GetAxis(4)->SetRangeUser(0.98, 1.1);
    dphiHPhi->GetAxis(4)->SetRangeUser(0.98, 1.1);
    Double_t Lsideband = dphiHPhi->Projection(4)->Integral((int)((0.994-invMin)*binsPerMass), (int)((1.006-invMin)*binsPerMass));
    Double_t likeSignLSideBand = dphiHKK->Projection(4)->Integral((int)((0.994-invMin)*binsPerMass), (int)((1.006-invMin)*binsPerMass));
    Double_t LscaleFactor = Lsideband/likeSignLSideBand;
    printf("sideband: %f, likeSignSideBand: %f, scaleFactor: %f \n", Lsideband, likeSignLSideBand, LscaleFactor);
    dphiLSLSideband->Scale(LscaleFactor);
    
    TCanvas *cDphiHPhiL3 = new TCanvas("cDphiHPhiL3", "cDphiHPhiL3", 50, 50, 600, 600);
    cDphiHPhiL3->cd();
    dphiUSLSideband->SetLineWidth(2);
    dphiUSLSideband->SetTitle("#Delta#varphi for Hadron-#Phi(1020) in Left Sideband Region");
    dphiUSLSideband->GetXaxis()->SetTitle("#Delta#varphi");
    dphiUSLSideband->Rebin(4);
    dphiLSLSideband->Rebin(4);
    dphiUSLSideband->Draw("SAME HIST E");
    dphiLSLSideband->SetLineColor(2);
    dphiLSLSideband->Draw("SAME HIST E");

    TPaveText *LsidebandText = new TPaveText(0.537, 0.716, 0.852, 0.818, "NDC");
    LsidebandText->AddText("0.994 < m_{kk} < 1.006 GeV/c^{2}"); 
    LsidebandText->AddText(Form("LS scale factor: %f", LscaleFactor));
    LsidebandText->SetFillColor(0);
    LsidebandText->SetBorderSize(1);
    LsidebandText->Draw("SAME");

// Setting & plotting delta-phi distribution for Unlike and Like sign pairs in the phi mass peak
// region (scaled by the ratio of the integral of the sideband region).
    dphiHPhi->GetAxis(4)->SetRangeUser(1.010,1.030);
    TH1D* dphiUSPeak = dphiHPhi->Projection(2);
    dphiUSPeak->SetName("dphiUSPeak");
    dphiUSPeak->Sumw2();
    dphiHKK->GetAxis(4)->SetRangeUser(1.010,1.030);
    TH1D* dphiLSPeak = dphiHKK->Projection(2);
    dphiLSPeak->SetName("dphiLSPeak");
    dphiLSPeak->Sumw2();
    dphiLSPeak->Scale(scaleFactor);

    TCanvas *cDphiHPhi2 = new TCanvas("cDphiHPhi2", "cDphiHPhi2", 50, 50, 600, 600);
    cDphiHPhi2->cd();
    dphiUSPeak->SetLineWidth(2);
    dphiUSPeak->SetTitle("#Delta#varphi for Hadron-#Phi(1020) in Peak Region");
    dphiUSPeak->GetXaxis()->SetTitle("#Delta#varphi");
    dphiUSPeak->Rebin(4);
    dphiUSPeak->GetYaxis()->SetRangeUser(300,3000);
    dphiUSPeak->Draw("SAME HIST E");
    dphiLSPeak->SetLineColor(2);
    dphiLSPeak->Rebin(4);
    dphiLSPeak->Draw("SAME HIST E");

    TPaveText *peakText = new TPaveText(0.537, 0.716, 0.852, 0.818, "NDC");
    peakText->AddText("1.010 < m_{kk} < 1.030 GeV/c^{2}"); 
    peakText->AddText(Form("LS scale factor: %f", scaleFactor));
    peakText->SetFillColor(0);
    peakText->SetBorderSize(1);
    peakText->Draw("SAME");


//Calculating and plotting the difference between US and LS pairs for the phi mass peak region (second method)
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

    //corrDPhiPeakRegion->Rebin(17, "hnew", xbins);
    //corrDPhiPeakRegion->Rebin(4);
    corrDPhiPeakRegion->SetName("hnew");

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

    dphifit->SetLineColor(4);
    dphifit->SetLineStyle(7);
    dphifit->SetLineWidth(3);

    hnew->SetLineWidth(2);

    hnew->Fit("dphifit", "R");


    hnew->SetTitle("Corrected #Delta#varphi for Hadron-#Phi(1020) in Peak Region");
    hnew->GetXaxis()->SetTitle("#Delta#varphi");
    hnew->GetXaxis()->SetRangeUser(-1.3, 4.4);
    
    hnew->Draw("SAME HIST E");

    //plotting two sidebands on same plot:
    TH1D* rightSide = dphiUSSideband->Clone("rightSide");
    TH1D* leftSide = dphiUSLSideband->Clone("leftSide");
   
    rightSide->Scale(1.0/(rightSide->Integral()));
    leftSide->Scale(1.0/leftSide->Integral());

    TCanvas* ccompare = new TCanvas("ccompare", "ccompare", 50, 50, 600, 600);
    ccompare->cd();
    rightSide->SetLineColor(4);
    rightSide->SetLineWidth(2);
    leftSide->SetLineWidth(2);
    rightSide->Draw("SAME H E");
    leftSide->Draw("SAME H E");

    return hnew;
}

void fitZVtx(TH1D* zVtx){
    TF1 *fit = new TF1("zVtxFit", "gaus(0)",-20, 20);
    fit->SetParameters(1000000.0, 1.5, 10.0);

    TCanvas *cZFit = new TCanvas("cZFit", "cZFit", 50,50,600,600);
    cZFit->cd();
    zVtx->Fit("zVtxFit");
    zVtx->GetXaxis()->SetRangeUser(-20, 20);
    zVtx->Draw();

    Double_t x[11] = {-10.0, -6.15, -3.90, -2.13, -0.59, 0.86, 2.29, 3.77, 5.39, 7.30, 10.0};
    Int_t integral= 0;
    TLine *lines[11];
    for(int i=0; i<10; i++){
        lines[i] = new TLine(x[i], 0.0, x[i], zVtx->GetBinContent(zVtx->GetXaxis()->FindBin(x[i])));
        lines[i]->Draw("SAME");
        integral = zVtx->Integral(zVtx->GetXaxis()->FindBin(x[i]), zVtx->GetXaxis()->FindBin(x[i+1]));
        printf("zVtx bin%i: %i\n", i, integral);
    }
    lines[10] = new TLine(x[10], 0.0, x[10], zVtx->GetBinContent(zVtx->GetXaxis()->FindBin(x[i])));
    lines[10]->Draw("SAME");
}

/************************
 ***   MAIN FUNCTION  ***
 ************************/
void plot_phi_histo(string inputName){

    //Set global style stuff
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

    //open file, get necessary histograms
    TFile *histoFile = new TFile(inputName.c_str());
    histoFile->cd("PhiReconstruction");
    THnSparseF *fkkUSDist = (THnSparseF *)InvMass->FindObject("fkkUSDist");
    THnSparseF *fkkLSDist = (THnSparseF *)InvMass->FindObject("fkkLSDist");
    THnSparseF *fTrigDist = (THnSparseF *)InvMass->FindObject("fTrigDist");
    THnSparseF *dphiHPhi = (THnSparseF *)InvMass->FindObject("fDphiHPhi");
    THnSparseF *dphiHKK = (THnSparseF *)InvMass->FindObject("fDphiHKK");
    THnSparseF *dphiHPhiMixed = (THnSparseF *)InvMass->FindObject("fDphiHPhiMixed");
    THnSparseF *dphiHKKMixed = (THnSparseF *)InvMass->FindObject("fDphiHKKMixed");
    TH1D* zVtx = InvMass->FindObject("fVtxZ");


    //set titles for the various distribution historgram axes
    if(fkkUSDist && fkkLSDist && fTrigDist){
        fkkLSDist->GetAxis(0)->SetTitle("p_{T} (GeV/c)");
        fkkUSDist->GetAxis(0)->SetTitle("p_{T} (GeV/c)");
        fTrigDist->GetAxis(0)->SetTitle("p_{T} (GeV/c)");

        fkkLSDist->GetAxis(2)->SetTitle("#varphi");
        fkkUSDist->GetAxis(2)->SetTitle("#varphi");
        fTrigDist->GetAxis(1)->SetTitle("#varphi");

        fkkLSDist->GetAxis(3)->SetTitle("#eta");
        fkkUSDist->GetAxis(3)->SetTitle("#eta");
        fTrigDist->GetAxis(2)->SetTitle("#eta");
    }else{
      printf("couldn't open kk distributions!\n");
      return 0;
    }


    //plot the invariant mass distributions for different pT bins for US KK, LS KK, and Corrected
    //phi's, as well as a phi pT spectrum
    //plotRawAndCorrectedPT(fkkUSDist, fkkLSDist);

    // Plotting distributions for KK pairs (US and LS) and trigger particles 
//    plotAllDistributions(fkkLSDist, fkkUSDist, fTrigDist);
    
    TString suffix;
    //Set the Pt ranges for trigger and assoc particles
    dphiHPhi->GetAxis(0)->SetRangeUser(4.0, 8.0); 
    dphiHPhi->GetAxis(1)->SetRangeUser(2.0,4.0); 
    dphiHKK->GetAxis(0)->SetRangeUser(4.0, 8.0);
    dphiHKK->GetAxis(1)->SetRangeUser(2.0,4.0); 
    //Do 2D Correlation plot (and rough "corrected" 2D correaltion plot)
    suffix = "2_4";
    TH3D *dEtadPhiDist = dphiHPhi->Projection(2, 3, 4);
    dEtadPhiDist->Rebin3D(4,1,1);
    TH3D *dEtadPhiLSDist = dphiHKK->Projection(2, 3, 4);
    dEtadPhiLSDist->Rebin3D(4,1,1);

//    TH3D *dEtadPhiUS = mixedEventCorrection(dphiHPhi, dphiHPhiMixed);
//    TH3D *dEtadPhiLS = mixedEventCorrection(dphiHKK, dphiHKKMixed);
//d    plot2DCorrelations(dEtadPhiDist, dEtadPhiLSDist, suffix);
/*
    dphiHPhi->GetAxis(1)->SetRangeUser(1.0,2.0); 
    dphiHKK->GetAxis(1)->SetRangeUser(1.0,2.0); 
    //Do 2D Correlation plot (and rough "corrected" 2D correaltion plot)
    suffix = "1_2";
    TH3D *dEtadPhiDist2 = dphiHPhi->Projection(2, 3, 4);
    dEtadPhiDist2->Rebin3D(4,1,1);
    TH3D *dEtadPhiLSDist2 = dphiHKK->Projection(2, 3, 4);
    dEtadPhiLSDist2->Rebin3D(4,1,1);
    plot2DCorrelations(dEtadPhiDist2, dEtadPhiLSDist2, suffix);

    dphiHPhi->GetAxis(1)->SetRangeUser(2.0,3.0); 
    dphiHKK->GetAxis(1)->SetRangeUser(2.0,3.0); 
    //Do 2D Correlation plot (and rough "corrected" 2D correaltion plot)
    suffix = "2_3";
    TH3D *dEtadPhiDist3 = dphiHPhi->Projection(2, 3, 4);
    dEtadPhiDist3->Rebin3D(4,1,1);
    TH3D *dEtadPhiLSDist3 = dphiHKK->Projection(2, 3, 4);
    dEtadPhiLSDist3->Rebin3D(4,1,1);
    plot2DCorrelations(dEtadPhiDist3, dEtadPhiLSDist3, suffix);

    dphiHPhi->GetAxis(1)->SetRangeUser(3.0,4.0); 
    dphiHKK->GetAxis(1)->SetRangeUser(3.0,4.0); 
    //Do 2D Correlation plot (and rough "corrected" 2D correaltion plot)
    suffix = "3_4";
    TH3D *dEtadPhiDist4 = dphiHPhi->Projection(2, 3, 4);
    dEtadPhiDist4->Rebin3D(4,1,1);
    TH3D *dEtadPhiLSDist4 = dphiHKK->Projection(2, 3, 4);
    dEtadPhiLSDist4->Rebin3D(4,1,1);
    plot2DCorrelations(dEtadPhiDist4, dEtadPhiLSDist4, suffix);

    dphiHPhi->GetAxis(1)->SetRangeUser(4.0,5.0); 
    dphiHKK->GetAxis(1)->SetRangeUser(4.0,5.0); 
    //Do 2D Correlation plot (and rough "corrected" 2D correaltion plot)
    suffix = "4_5";
    TH3D *dEtadPhiDist5 = dphiHPhi->Projection(2, 3, 4);
    dEtadPhiDist5->Rebin3D(4,1,1);
    TH3D *dEtadPhiLSDist5 = dphiHKK->Projection(2, 3, 4);
    dEtadPhiLSDist5->Rebin3D(4,1,1);
    plot2DCorrelations(dEtadPhiDist5, dEtadPhiLSDist5, suffix);
*/

    TH1D* corrV1 = plotPhiCorrelationsV1(dphiHPhi, dphiHKK);
    //corrV1->Sumw2();

    //Reset the Delta-phi axis range after the dphi binned projections are done being created above.
    dphiHPhi->GetAxis(2)->SetRange(0,0);
    dphiHKK->GetAxis(2)->SetRange(0,0);

    TH1D* corrV2 = plotPhiCorrelationsV2(dphiHPhi, dphiHKK);
    //corrV2->Sumw2();

    TH1D* ratio = corrV1->Clone("ratio");
    ratio->Divide(corrV1, corrV2);

    TCanvas *cratio = new TCanvas("cratio", "cratio", 50, 50, 600, 600);
    cratio->cd();
    ratio->Draw("H E");
    TLine *line = new TLine(-1.57, 1.0, 4.7, 1.0);
    line->SetLineColor(2);
    line->SetLineStyle(7);
    line->Draw("SAME");

//    fitZVtx(zVtx);
}


