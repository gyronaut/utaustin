#if !defined (__CINT__) || (defined(__MAKECINT__))
#include <iostream>
#include "TMath.h"
#include "TF1.h"

#endif
using namespace std;

Double_t BlastIntegrand(const Double_t *x,const Double_t *par){
  Double_t x0=x[0];
  Double_t m=par[0];
  Double_t t=fabs(par[1]);
  Double_t n=par[2];
  Double_t beta_max=par[3];
  Double_t pt=par[4];

  //Keep beta within reasonable limits.
  Double_t beta=beta_max*TMath::Power(x0,n);
  if(beta>0.9999999999999999) beta=0.9999999999999999;

  Double_t mt=TMath::Sqrt(m*m+pt*pt);
  Double_t rho0=TMath::ATanH(beta);
  Double_t a0=pt*TMath::SinH(rho0)/t;
  if(a0>700.) a0=700.;//avoid floating point exception
  Double_t a1=mt*TMath::CosH(rho0)/t;
  return x0*mt*TMath::BesselI0(a0)*TMath::BesselK1(a1);
}

Double_t Blast(Double_t *x,Double_t *par){
  static TF1* fint=0;
  if(!fint) fint=new TF1("fint",BlastIntegrand,0.,1.,5);

  fint->SetParameters(par[0],par[2],par[3],par[4],x[0]);
  return fint->Integral(0.,1.)*par[1];
}

Double_t  XBlast(Double_t *x,Double_t *par){
  return x[0]*Blast(x,par);
}


void BlastParameters(double& t,double& n,double & b){
  t=0.151;
  n=2.262;
  b=0.357;
  return;
}

Double_t GBC(TH1D* h, Int_t bin){
    return h->GetBinCenter(bin);
}

Double_t GBW(TH1D* h, Int_t bin){
    return h->GetXaxis()->GetBinWidth(bin);
}

void BlastWaveTests(){

    TH1D* h;//define your histogram

    double fmin,fmax;//define limits of fit
    fmin = 0.0;
    fmax = 10.0;
    TF1* bwPhi=new TF1("blast_phi",XBlast,fmin,fmax,5);
    TF1* bwPion = new TF1("blast_pion", XBlast, fmin, fmax, 5);
    TF1* bwProton = new TF1("blast_proton", XBlast, fmin, fmax, 5);

    //fill initial values for fit parameters
    double bt,bn,bb;
    BlastParameters(bt,bn,bb);

    //get integral of histogram, used to set initial scale
    double hint=0.;
    /*for(int pb=1;pb<=h->GetNbinsX();pb++) hint+=GBC(h,pb)*GBW(h,pb);*/
    hint = 100.0;

    //set initial values of parameters
    bwPhi->SetParError(0,0.);
    bwPhi->SetParameter(0,1.019);
    bwPhi->FixParameter(0,1.019);//fix the mass of the particle

    bwPhi->SetParError(1,0.1*hint);
    bwPhi->SetParameter(1,hint);
    bwPhi->FixParameter(1,hint);

    bwPhi->SetParError(2,0.1*bt);//temperature
    bwPhi->SetParameter(2,bt);
    bwPhi->FixParameter(2,bt);

    bwPhi->SetParError(3,0.1*bn);
    bwPhi->SetParameter(3,bn);
    bwPhi->FixParameter(3,bn);

    bwPhi->SetParError(4,0.1*bb);
    bwPhi->SetParameter(4,bb);
    bwPhi->FixParameter(4,bb);

    //values for pion BW
    bwPion->SetParError(0,0.);
    bwPion->SetParameter(0,0.1396);
    bwPion->FixParameter(0,0.1396);//fix the mass of the particle

    bwPion->SetParError(1,0.1*hint);
    bwPion->SetParameter(1,10.0*hint);
    bwPion->FixParameter(1,10.0*hint);

    bwPion->SetParError(2,0.1*bt);//temperature
    bwPion->SetParameter(2,bt);
    bwPion->FixParameter(2,bt);

    bwPion->SetParError(3,0.1*bn);
    bwPion->SetParameter(3,bn);
    bwPion->FixParameter(3,bn);

    bwPion->SetParError(4,0.1*bb);
    bwPion->SetParameter(4,bb);
    bwPion->FixParameter(4,bb);

    //values for proton BW

    bwProton->SetParError(0,0.);
    bwProton->SetParameter(0,0.938);
    bwProton->FixParameter(0,0.938);//fix the mass of the particle

    bwProton->SetParError(1,0.1*hint);
    bwProton->SetParameter(1,5.0*hint);
    bwProton->FixParameter(1,5.0*hint);

    bwProton->SetParError(2,0.1*bt);//temperature
    bwProton->SetParameter(2,bt);
    bwProton->FixParameter(2,bt);

    bwProton->SetParError(3,0.1*bn);
    bwProton->SetParameter(3,bn);
    bwProton->FixParameter(3,bn);

    bwProton->SetParError(4,0.1*bb);
    bwProton->SetParameter(4,bb);
    bwProton->FixParameter(4,bb);


    int j,status;
    TFitResultPtr fr;

    bwPhi->SetLineColor(kAzure+1);
    bwPion->SetLineColor(kRed+1);
    bwProton->SetLineColor(kGreen+2);

    TCanvas* c1 = new TCanvas("c1", "c1", 50, 50, 600, 600);
    c1->cd()->SetLogy();
    bwPion->Draw();
    bwPhi->Draw("SAME");
    bwProton->Draw("SAME");

    Double_t temps[] = {0.150, 0.140, 0.130, 0.120, 0.110, 0.100};

    Double_t pionInt = 0.0, protonInt = 0.0, phiInt = 0.0;
    for(int i = 0; i < 6; i++){
        printf("temp: %f\n", temps[i]);

        bwPion->SetParameter(2,temps[i]);
        bwPion->FixParameter(2,temps[i]);

        bwPhi->SetParameter(2,temps[i]);
        bwPhi->FixParameter(2,temps[i]);

        bwProton->SetParameter(2,temps[i]);
        bwProton->FixParameter(2,temps[i]);

        pionInt = bwPion->Integral(0.0, 10.0);
        protonInt = bwProton->Integral(0.0, 10.0);
        phiInt = bwPhi->Integral(0.0, 10.0);

        printf("total pion: %e, total phi: %e, ratio: %f\n", pionInt, phiInt, phiInt/pionInt);

        bwPion->SetParameter(1, 20.0/pionInt);
        bwProton->SetParameter(1, 10.0/protonInt);
        bwPhi->SetParameter(1, 1.0/phiInt);

        pionInt = bwPion->Integral(0.0, 10.0);
        protonInt = bwProton->Integral(0.0, 10.0);
        phiInt = bwPhi->Integral(0.0, 10.0);


        Double_t pion2_4Int = bwPion->Integral(2.0, 4.0);
        Double_t phi2_4Int = bwPhi->Integral(2.0, 4.0);

        printf("total pion: %e, total phi: %e, ratio: %f\n, 2-4 pion: %e, 2-4 phi: %e, 2-4 ratio: %f\n", pionInt, phiInt, phiInt/pionInt, pion2_4Int, phi2_4Int, phi2_4Int/pion2_4Int);
        printf("pion ratio: %2.2e%%, phi ratio: %2.2e%%, ratio ratio: %2.2e\n", 100.0*pion2_4Int/pionInt, 100.0*phi2_4Int/phiInt, phi2_4Int*pion2_4Int/(phiInt*pionInt));
        printf("===========================\n");
    }
    /*
    for(j=0;j<100;j++){
        fr=h->Fit(bwPhi,"RSQ0N");
        status=fr->Status();
        if(!status) break;
    }

    bwPhi->ReleaseParameter(2);
    bwPhi->SetParLimits(2,0.,10.);

    for(j=0;j<100;j++){
        fr=h->Fit(bwPhi,"RSQ0N");
        status=fr->Status();
        if(!status) break;
    }

    bwPhi->ReleaseParameter(3);
    bwPhi->SetParLimits(3,0.,100.);

    for(j=0;j<100;j++){
        fr=h->Fit(bwPhi,"RSQ0N");
        status=fr->Status();
        if(!status) break;
    }

    bwPhi->ReleaseParameter(4);
    bwPhi->SetParLimits(4,TMath::Min(0.6,0.8*bb),TMath::Min(TMath::Max(0.9,1.2*bb),0.99));

    for(j=0;j<100;j++){
        fr=h->Fit(bwPhi,"RSQ0N");
        status=fr->Status();
        if(!status) break;
    }

    for(j=0;j<100;j++){
        fr=h->Fit(bwPhi,"RSQ0NI");
        status=fr->Status();
        if(!status) break;
    }

    //These last two steps can be very time-consuminbwPhi, so you may want to skip them at first.
    for(j=0;j<100;j++){
        fr=h->Fit(bwPhi,"RSQ0NIE");
        status=fr->Status();
        if(!status) break;
    }

    for(j=0;j<2;j++){
        fr=h->Fit(bwPhi,"RSQ0NIEM");
        status=fr->Status();
        if(!status) break;
    }

    fr->Print();
    cerr<<"chi2/NDF = "<<bwPhi->GetChisquare()<<"/"<<bwPhi->GetNDF()<<" = "<<bwPhi->GetChisquare()/bwPhi->GetNDF()<<endl;
    if(status && status!=4000) cerr<<"BAD_FIT "<<h->GetName()<<endl;
    //The status should be 0 or 4000.
    */
}
