#include "TLegend.h"
#include <string>
#include <sstream>

/* This function takes as input the full path to a root file that contains the
 * THnSparse histograms (i.e. "/path/to/input.root"), opens up the file, and
 * loads the THnSparse histograms, and performs certain projections/cuts
 */

void sparseReader(string inputName){
    
    // Setting up some global style parameters (that are applied automatically when you draw a canvas)
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

    // Open up the file that contains the sparse histograms
    TFile *histoFile = new TFile(inputName.c_str());
   
    // Change to a specific "Directory" within the root file (you can see what it's called if you open up
    // the root file in a TBrowser). Once you've changed to this directory, you can use the variable
    // InvMass and the function InvMass->FindObject("histo_name") to load histograms
    histoFile->cd("PhiPIDCuts");

    //THnSparse that contains the following information for each K+K- pair:
    //  K+ TPCnSigma          (axis 0) from -5 to 5
    //  K- TPCnSigma          (axis 1) from -5 to 5
    //  K+ TOFnSigma          (axis 2) from -5 to 5
    //  K- TOFnSigma          (axis 3) from -5 to 5
    //  (K+K-) pT             (axis 4) from 0.1 to 10.0 GeV/c
    //  (K+K-) Invariant Mass (axis 5) from 0.98 to 1.1 GeV/c^2
    //  
    //  "US" is the unlike-sign kaons pairs (K+K-), while "LS" is the like-sign pairs (K+K+ and K-K-).
    //  "US" is where we are looking for the phi meson, "LS" can be used to estimate the background.
    THnSparseF *fkkUSDist = (THnSparseF *)InvMass->FindObject("fkkUSDist");
    THnSparseF *fkkLSDist = (THnSparseF *)InvMass->FindObject("fkkLSDist");

    /* Setting up the Titles for each of the axes of the thnsparse for US Kaon pairs*/
    if(fkkUSDist && fkkLSDist){
        fkkUSDist->GetAxis(0)->SetTitle("K+ n#sigma_{TPC}");
        fkkUSDist->GetAxis(1)->SetTitle("K- n#sigma_{TPC}");
        fkkUSDist->GetAxis(2)->SetTitle("K+ n#sigma_{TOF}");
        fkkUSDist->GetAxis(3)->SetTitle("K- n#sigma_{TOF}");
        fkkUSDist->GetAxis(4)->SetTitle("KK pT");
        fkkUSDist->GetAxis(5)->SetTitle("KK Invariant Mass");
    }else{
      printf("couldn't open KK distributions!\n");
      return 0;
    }

    // Setting the range of an axis using SetRange(firstBin, lastBin) uses the bin numbers
    // Setting the range of an axis using SetRangeUser(firstValue, lastValue) uses the number values from the histogram
    
    fkkUSDist->GetAxis(4)->SetRangeUser(1.0, 3.0); //sets the range of the pT axis to be between 1.0 and 3.0 GeV

    // Projecting the THnSparse onto the invariant mass axis (axis 5) and calling it a new 1D histogram called invariantMass
    TH1D *invariantMass = (TH1D*)fkkUSDist->Projection(5);

    TCanvas *canvas = new TCanvas("invmassCanvas", "invmassCanvas", 0,0,600,600);
    canvas->cd();
    invariantMass->Draw();
}


