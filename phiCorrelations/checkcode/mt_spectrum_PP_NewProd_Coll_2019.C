#include <istream>
#include <fstream>
#include "math.h"
#include <string>
#include "TFile.h"
#include "TH1.h"
#include "TObjArray.h"
#include <algorithm>
#include <TPaveText.h>
//Root
#include "TNtuple.h"
#include "TStyle.h"
#include "TMatrix.h"
#include "TGraphErrors.h"
#include "TF1.h"
//#include "TH1D.h"
#include "TCanvas.h"
//#include "TH2F.h"
#include "TGaxis.h"

void fit2meanPT(TF1* boltzmann, Float_t meanpt){
    Float_t step = 0.1;
    Int_t count = 0;
    Int_t maxcount = 10000;
    Float_t calcmean = boltzmann->Moment(1, 0, 10);
    while(abs(calcmean-meanpt) > meanpt*0.005){
        boltzmann->SetParameter(1, boltzmann->GetParameter(1)*(1+((meanpt-calcmean)/abs(calcmean-meanpt)*step)));
        Float_t newmean = boltzmann->Moment(1, 0, 10);
        if((newmean-meanpt)*(calcmean-meanpt) > 0){
            step = step/2.0;
        }
        calcmean = newmean;
        count+=1;
        if(count>maxcount){
            printf("looped more than %d times for %s, diff = %f\n", maxcount, boltzmann->GetName(), newmean-meanpt);
            break;
        }
    }
}

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
//  if(t==0)printf("oh no, dividing by 0...");
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

Double_t  XXBlast(Double_t *x,Double_t *par){
  return x[0]*XBlast(x,par);
}


void BlastParameters(double& t,double& n,double & b){
  t=0.143;
  n=1.07;
  b=0.547;
  return;
}

Double_t GBC(TH1D* h, Int_t bin){
    return h->GetBinContent(bin);
}

Double_t GBW(TH1D* h, Int_t bin){
    return h->GetXaxis()->GetBinWidth(bin);
}


void mt_spectrum_PP_NewProd_Coll_2019()
{
    gROOT->Reset();
    gROOT->Time();
    Int_t i;
    Float_t nsig;
    Float_t nbac;
    Float_t fac;

    //gStyle->SetStatX(100);
    //gStyle->SetTitleX(100);

    gStyle->SetOptFit(0);
    gStyle->SetOptTitle(0);
    gStyle->SetOptStat(0);

    TCanvas *c = new TCanvas("c", "c",0,0,350,600);
    c->SetFillColor(10);
    c->SetBorderMode(0);
    c->Divide(1,3);
    //c->SetBorderSize(2);

    gStyle->SetTitleColor(10);
    //-------  stat box -------------------
    gStyle->SetStatColor(0);
    gStyle->SetStatH(0.16);
    gStyle->SetStatW(0.27);
    gStyle->SetStatX(0.9);
    gStyle->SetStatY(0.9);
    //gStyle->SetBorderSize(1);

    //gStyle->SetLabelOffset(1.2);
    gStyle->SetLabelFont(22);
    gStyle->SetFillColor(10);
    gStyle->SetCanvasColor(10);
    //------------pad------------------
    gStyle->SetPadLeftMargin(0.18);
    gStyle->SetPadBottomMargin(0.21);
    //--------- histogram---------------- 
    //gStyle->SetHistLineWidth(4);
    //gStyle->SetHistFillColor(1);
    //gStyle->SetLabelFont(22,"x");
    //gStyle->SetLabelFont(22,"y");
    //gStyle->SetLabelFont(22,"z");

    c->cd();
    c->Modified();
    c->cd();

    TH2F *h = new TH2F("pt","pt",100,0,2.5,100,0,100);
    TH2F *hraw = new TH2F("ptraw","ptraw",100,0,3,100,0,5000);
    TH2F *hraw2 = new TH2F("ptraw2","ptraw2",100,0,3,100,0,3000);

    TH2F *hacc = new TH2F("ptacc","ptacc",100,0,2,100,0,1);
    TH2F *hrawacc = new TH2F("ptrawacc","ptrawacc",100,0,2,100,0,1500);

    TH1F *histo = new TH1F("ptagain","ptagain",10,0,3);

    TH1F *phiyield = new TH1F("phiyield","phiyield",2,0,2);
    TH1F *pionyield = new TH1F("pionyield","pionyield",2,0,2);
    TH1F* yieldratio = new TH1F("ratioyield","ratioyield",2,0,2);

    c->cd(1);  
    c->SetFillColor(10);
    c->SetBorderSize(2);
    c->SetLeftMargin(0.19);
    c->SetBottomMargin(0.19); 

    h->SetLineColor(1);
    h->SetLineWidth(4);
    h->GetYaxis()->SetLabelSize(0.1);
    h->GetXaxis()->SetLabelSize(0.1);
    h->GetYaxis()->SetLabelOffset(0.01);
    h->GetXaxis()->SetLabelOffset(0.00);

    h->GetYaxis()->SetTitleSize(0.1);
    h->GetXaxis()->SetTitleSize(0.1);
    h->GetYaxis()->SetTitleOffset(1.2);
    h->GetXaxis()->SetTitleOffset(1.2);
    h->GetXaxis()->SetTitle("");
    h->GetYaxis()->SetTitle("");
    h->GetXaxis()->SetNdivisions(5);
    h->GetYaxis()->SetNdivisions(6);
    h->GetXaxis()->SetLabelFont(22);
    h->GetYaxis()->SetLabelFont(22);
    h->GetXaxis()->SetTitleFont(22);
    h->GetYaxis()->SetTitleFont(22);

    TGaxis *g = new TGaxis(0,0,6,0,0,6,6,""); 
    g->SetLabelOffset(h->GetXaxis()->GetLabelOffset());
    g->SetLabelSize(h->GetXaxis()->GetLabelSize());  
    g->SetTickSize(0.03);
    g->SetGridLength(0);
    g->SetTitleOffset(1);
    g->SetLabelOffset(0.005);
    g->SetLabelFont(20);
    g->SetName(" ");

    TGaxis *g1 = new TGaxis(0,0,0,0,0,3.5,7,""); 
    g1->SetLabelOffset(h->GetXaxis()->GetLabelOffset());
    g1->SetLabelSize(h->GetXaxis()->GetLabelSize());  
    g1->SetTickSize(0.03);
    g1->SetGridLength(0);
    g1->SetTitleOffset(1);
    g1->SetLabelOffset(0.005);
    g1->SetLabelFont(20);
    g1->SetName(" ");

    phiyield->SetLineColor(kGreen+2);
    phiyield->SetLineWidth(2);
    phiyield->SetMarkerColor(kGreen+2);
    phiyield->SetMarkerSize(2);
    phiyield->SetMarkerStyle(22);

    pionyield->SetLineColor(kBlue+1);
    pionyield->SetLineWidth(2);
    pionyield->SetMarkerColor(kBlue+1);
    pionyield->SetMarkerSize(2);
    pionyield->SetMarkerStyle(21);

    //////////////////// stat line /////////////

    Float_t  pt[5] = {0.4,0.8,1.2,1.6,2.0};

    ////////////////////////////////////////////////////
    // 0.8 cos theta 2cm dca 2.5 sigmadedx
    //lambda1520 + bar
    //Float_t  yield[5] ={244.0,947.0,1023.0,620.0,33.0};

    //lambda1520
    //Float_t  yield[5] ={159.0,603.0,569.0,320.0,27.0};
    //lambda1520bar
    //Float_t  yield[5] ={90.0,348.0,462.0,307.0,16.0};

    //lambda1520 + bar
    //Float_t  yielderr[5] = {57.0,125.0,142.0,96.0,35.0};
    //lambda1520
    //Float_t  yielderr[5] = {45.0,95.0,106.0,71.0,26.0};
    //lambda1520 bar
    //Float_t  yielderr[5] = {36.0,81.0,94.0,64.0,23.0};

    //correction 
    //Float_t  yield2[5] ={0.0217,0.0830,0.175,0.159,0.051};

    /////////////////////////////////////////////////////    
    ///////////all cos 0.9 dca 2 //////////////////////
    //Float_t  yield[5] ={217.0,946.0,1049.0,534.0,36.0};
    //Float_t  yielderr[5] = {54.0,103.0,114.0,91.0,35.0};
    //Float_t  yield2[5] ={0.0598377,0.191638,0.33207,0.253979,0.07389};
    ///////////////////////////////////////////////////////////
    ///////// all cos 0.9 dca 3 
    Float_t  yield[5] ={235.0,952.0,1156.0,620.0,36.0};
    //Float_t  yield[5] ={239.0,930.0,1061.0,619.0,65.0};

    Float_t  yielderr[5] = {63.0,112.0,123.0,107.0,37.0};

    Float_t  yield2[5] ={0.062373,0.193031,0.33249,0.253979,0.07389};
    ////////////////////////////////////////////////////



    Float_t  pterr[5] = {0.0,0.0,0.0,0.0,0.0};  
    Float_t  yielderr2[5] = {0.0,0.0,0.0,0.0,0.0};
    Float_t  pterr2[5] = {0.0,0.0,0.0,0.0,0.0};



    // y [-1,1]
    /*
       Float_t  yield2[5] ={0.0327,0.106,0.1847,0.1333,0.0390};
       Float_t  yield2[5] ={0.0583,0.1756,0.307,0.299,0.073};
       Float_t  yield2[5] ={0.0217,0.0830,0.175,0.159,0.051};
       Float_t  yield2[5] ={0.0387,0.139,0.288,0.269,0.0967};
       Float_t  yield2[5] ={0.0327,0.106,0.1847,0.1333,0.0390};
       Float_t  yield2[5] ={0.0217,0.0830,0.175,0.159,0.051};
       Float_t  yield2[5] ={0.0327,0.106,0.1847,0.1333,0.0390};
       */

    TH2F *histo6 = new TH2F("histo6","amplituds of fit",30,0,3,100,0,1500);

    histo6->SetLineColor(1);
    histo6->SetLineWidth(4);
    histo6->GetYaxis()->SetLabelSize(0.06);
    histo6->GetXaxis()->SetLabelSize(0.06);
    histo6->GetYaxis()->SetLabelOffset(0.01);
    histo6->GetXaxis()->SetLabelOffset(0.00);

    histo6->GetYaxis()->SetTitleSize(0.06);
    histo6->GetXaxis()->SetTitleSize(0.06);
    histo6->GetYaxis()->SetTitleOffset(1.5);
    histo6->GetXaxis()->SetTitleOffset(1.3);
    histo6->GetXaxis()->SetTitle("pt");
    histo6->GetYaxis()->SetTitle("Raw yield");
    histo6->GetXaxis()->SetNdivisions(6);
    histo6->GetYaxis()->SetNdivisions(6);
    histo6->GetXaxis()->SetLabelFont(22);
    histo6->GetYaxis()->SetLabelFont(22);
    histo6->GetXaxis()->SetTitleFont(22);
    histo6->GetYaxis()->SetTitleFont(22);

    c->cd(1);
    histo6->Draw();
    TGraph *graph6 = new TGraphErrors(5,pt,yield,pterr,yielderr);    
    graph6->SetMarkerStyle(20);     
    graph6->Draw("p");


    TH2F *histo7 = new TH2F("histo7","amplituds of fit",30,0,3,100,0,0.4);
    histo7->SetLineColor(1);
    histo7->SetLineWidth(4);
    histo7->GetYaxis()->SetLabelSize(0.06);
    histo7->GetXaxis()->SetLabelSize(0.06);
    histo7->GetYaxis()->SetLabelOffset(0.01);
    histo7->GetXaxis()->SetLabelOffset(0.00);

    histo7->GetYaxis()->SetTitleSize(0.06);
    histo7->GetXaxis()->SetTitleSize(0.06);
    histo7->GetYaxis()->SetTitleOffset(1.5);
    histo7->GetXaxis()->SetTitleOffset(1.3);
    histo7->GetXaxis()->SetTitle("pt");
    histo7->GetYaxis()->SetTitle("Acceptance+Efficiency");
    histo7->GetXaxis()->SetNdivisions(6);
    histo7->GetYaxis()->SetNdivisions(6);
    histo7->GetXaxis()->SetLabelFont(22);
    histo7->GetYaxis()->SetLabelFont(22);
    histo7->GetXaxis()->SetTitleFont(22);
    histo7->GetYaxis()->SetTitleFont(22);

    c->cd(2);
    histo7->Draw();
    TGraph *graph7 = new TGraphErrors(5,pt,yield2,pterr2,yielderr2);    
    graph7->SetMarkerStyle(20);     
    graph7->Draw("p");




    Float_t  yield3[5];
    Float_t  yielderr3[5];
    Float_t  pterr3[5] = {0.0,0.0,0.0,0.0,0.0};

    for(Int_t i=0 ;i<5;i++)
    {
        yield3[i] = yield[i]/ yield2[i];
        yielderr3[i] = yielderr[i]/  yield2[i];
    }


    TH2F *histo8 = new TH2F("histo8","amplituds of fit",30,0,3,100,0,0.01);

    histo8->SetLineColor(1);
    histo8->SetLineWidth(4);
    histo8->GetYaxis()->SetLabelSize(0.08);
    histo8->GetXaxis()->SetLabelSize(0.08);
    histo8->GetYaxis()->SetLabelOffset(0.01);
    histo8->GetXaxis()->SetLabelOffset(0.00);

    histo8->GetYaxis()->SetTitleSize(0.08);
    histo8->GetXaxis()->SetTitleSize(0.08);
    histo8->GetYaxis()->SetTitleOffset(1.2);
    histo8->GetXaxis()->SetTitleOffset(1.2);
    histo8->GetXaxis()->SetTitle("p_{T} [GeV]");
    histo8->GetYaxis()->SetTitle("dN/dy/dp_{T}");
    histo8->GetXaxis()->SetNdivisions(6);
    histo8->GetYaxis()->SetNdivisions(6);
    histo8->GetXaxis()->SetLabelFont(22);
    histo8->GetYaxis()->SetLabelFont(22);
    histo8->GetXaxis()->SetTitleFont(22);
    histo8->GetYaxis()->SetTitleFont(22);

    //Float_t scfac = (Float_t)((4.44/8416422)*(1.10)*1.10*2.5);

    Float_t scfac = (Float_t)((4.44/8416422)*(1.05)*1.01*2.5);


    for(Int_t i=0 ;i<5;i++)
    {
        yield3[i] = yield3[i] *scfac;
        yielderr3[i] = yielderr3[i]*scfac;
    }



    c->cd(3);
    histo8->Draw();
    TGraph *graph8 = new TGraphErrors(5,pt,yield3,pterr3,yielderr3);    
    graph8->SetMarkerStyle(20);     
    graph8->Draw("p");




    cout<< " TEST"<<endl; 

    Char_t  text6[256];
    sprintf(text6,"[0]*x * exp(- sqrt((x*x)+ (1.5195*1.5195)) * (1/[1])  )");

    // TF1 *ex6 = new TF1("ex6",text6,0,3);
    //ex6->SetParLimits(0,2,4);
    TF1 *ex6 = new TF1("ex6","[0]*x * exp(- sqrt((x*x)+ (1.5195*1.5195)) * (1/[1])  )",0,3);
    ex6->SetParLimits(1,0.300, 0.400);   
    ex6->SetLineColor(4);       
    //graph8->Fit("ex6","RMEI");

    //cout<< " integral "<< ex6->Integral(0,100)<<endl;
    //cout<< " integral "<< ex6->Integral(0,0.92)<<endl;

    // //cout<< " integral "<< ex6->Integral(0,100)*4.44/8416422<<endl;

    Char_t  fit[256];

    Double_t y0 = ex6->GetParameter(0);
    Double_t y1 = ex6->GetParameter(1);    

    cout<< y0<<endl;
    cout<< y1<<endl;



    // sprintf(fit,"%f*x * exp(- sqrt((x*x)+ (1.5195*1.5195)) * (1/%f)  )", y0,(y1+0.013));
    sprintf(fit,"%f*x * exp(- sqrt((x*x)+ (1.5195*1.5195)) * (1/%f)  )", y0,y1);

    TF1 *bw1 = new TF1("bw1",fit,0,3);
    //bw1->Integral(0,100);
    bw1->Draw("same");


    //////////////////////////////////////////   
    ////  new part
    ////////////////////////////////////////
    /// shift 
    ///////////////////////

    TCanvas *c2 = new TCanvas("c2", "c2",0,0,700,1000);
    c2->SetFillColor(10);
    c2->SetBorderMode(0);
    c2->Divide(2,3);
    c2->cd(1);

    Float_t mass8;
    mass8 = 1.0195;

    TH2F *histo_new = new TH2F("histo_new","amplituds of fit",60,0,6,100,0,1);
    histo_new->Draw();

    Char_t  fit_shift[256];
    Char_t  fit_shiftint[256];



    y0= 1;
    y1= 0.493;
    //sprintf(text9,"([0]/([1]*(1.5195+[1])))*x*exp(- (sqrt((x*x)+ (1.5195*1.5195))-1.5195) * (1/[1])  )");

    sprintf(fit_shift,"(%f/(%f*(1.0195+%f)))*x*exp(- (sqrt((x*x)+ (1.0195*1.0195))-1.0195) * (1/%f)  )",y0,y1,y1,y1);

    TF1 *bw1_shift = new TF1("bw1_shift",fit_shift,0,6);
    bw1_shift->SetLineColor(4);   
    //bw1->Integral(0,100);
    bw1_shift->Draw("same");

    cout<< " 1 integral all ----->"<< bw1_shift->Integral(0,100)<<endl;
    cout<< " 1 integral 2-4 ----->"<< bw1_shift->Integral(2,4)<<endl;



    // mean pt and 

    sprintf(fit_shiftint,"(%f/(%f*(1.0195+%f)))*x*x*exp(- (sqrt((x*x)+ (1.0195*1.0195))-1.0195) * (1/%f)  )",y0,y1,y1,y1);
    TF1 *bw1_shiftint = new TF1("bw1_shiftint",fit_shiftint,0,1000);
    cout<<"mean pt  --> "<<  bw1_shiftint->Integral(0,1000)/ bw1_shift->Integral(0,1000) <<endl;

    phiyield->SetBinContent(1, bw1_shift->Integral(2,4));


    /////////////////////////////////////////////////////
    cout<< "----------------------------------------->"<<endl;

    Char_t  fit_shift2[256];
    Char_t  fit_shift2int[256];

    y0= 1;
    y1= 0.595;
    sprintf(fit_shift2,"(%f/(%f*(1.0195+%f)))*x*exp(- (sqrt((x*x)+ (1.0195*1.0195))-1.0195) * (1/%f)  )",y0,y1,y1,y1);

    TF1 *bw1_shift2 = new TF1("bw1_shift2",fit_shift2,0,6);
    bw1_shift2->SetLineColor(3);   
    //bw1->Integral(0,100);
    bw1_shift2->Draw("same");

    cout<< " 2 integral all  ----->"<< bw1_shift2->Integral(0,100)<<endl;
    cout<< " 2 integral  2-4 ----->"<< bw1_shift2->Integral(2,4)<<endl;


    sprintf(fit_shift2int,"(%f/(%f*(1.0195+%f)))*x*x*exp(- (sqrt((x*x)+ (1.0195*1.0195))-1.0195) * (1/%f)  )",y0,y1,y1,y1);
    TF1 *bw1_shift2int = new TF1("bw1_shift2int",fit_shift2int,0,1000);
    cout<<"mean pt  --> "<<  bw1_shift2int->Integral(0,1000)/ bw1_shift2->Integral(0,1000) <<endl;
    cout<<"test 1st moment -->"<< bw1_shift2->Moment(1, 0, 10) <<endl;
    
    phiyield->SetBinContent(2, bw1_shift2->Integral(2,4));

    /////

    c2->cd(2)->SetLogy();
    TH2F *histo_new3 = new TH2F("histo_new3","amplituds of fit",60,0,6,100,0.001,1);

    histo_new3->Draw();
    bw1_shift->Draw("same");
    bw1_shift2->Draw("same");


    ////////////////////////////////////
    // TEST OF BOLTZMANN MEAN PT CALC //
    ////////////////////////////////////
    Float_t multbinlow[] = {100.0, 80.0, 60.0, 40.0, 20.0, 10.0, 5.0};
    Float_t multbinhigh[] = {80.0, 60.0, 40.0, 20.0, 10.0, 5.0, 0.0};
    Float_t multbincenters[] = {90, 70, 50, 30, 15, 7.5, 2.5};
    Float_t revbincenters[] = {10, 30.0, 50.0, 70.0, 85.0, 92.5, 97.5};
    Float_t revbinwidths[] = {10, 10, 10, 10, 5, 2.5, 2.5};
    Float_t cutrevbincenters[] = {30.0, 50.0, 70.0, 85.0, 92.5, 97.5};
    Float_t cutrevbinwidths[] = {10, 10, 10, 5, 2.5, 2.5};
    Float_t pionmeanpt[] = {0.4336, 0.4705, 0.4944, 0.5142, 0.5293, 0.5375, 0.5453};
    Float_t kaonmeanpt[] = {0.6809, 0.7722, 0.8317, 0.8764, 0.901, 0.9177, 0.9366};
    Float_t protonmeanpt[] = {0.8208, 0.9607, 1.053, 1.132, 1.186, 1.223, 1.248};
    Float_t phimeanpt[] = {1.055, 1.242, 1.310, 1.357, 1.421, 1.442, 1.437};

    TF1* phiboltz[7];
    TF1* pionboltz[7];
    TF1* kaonboltz[7];
    TF1* protonboltz[7];
    Float_t piratios[7];
    Float_t kratios[7];
    Float_t pratios[7];

    for(int i = 0; i < 7; i++){
        phiboltz[i] = new TF1(Form("phiboltz%d", i), "([0]/([1]*(1.0195 + [1])))*x*exp(-(sqrt((x*x) + (1.0195*1.0195))-1.0195)*(1/[1]))", 0, 6);
        phiboltz[i]->SetParameters(1.0, 0.44 + i*.01);
        phiboltz[i]->SetLineColor(kAzure+i);
        fit2meanPT(phiboltz[i], phimeanpt[i]);

        pionboltz[i] = new TF1(Form("pionboltz%d", i),  "([0]/([1]*(0.139 + [1])))*x*exp(-(sqrt((x*x) + (0.139*0.139))-0.139)*(1/[1]))", 0, 6);
        pionboltz[i]->SetParameters(1.0, 0.220 + i*.01);
        fit2meanPT(pionboltz[i], pionmeanpt[i]);
        pionboltz[i]->SetLineColor(kRed+1);

        kaonboltz[i] = new TF1(Form("kaonboltz%d", i),  "([0]/([1]*(0.494 + [1])))*x*exp(-(sqrt((x*x) + (0.494*0.494))-0.494)*(1/[1]))", 0, 6);
        kaonboltz[i]->SetParameters(1.0, 0.330 + i*.01);
        kaonboltz[i]->SetLineColor(kBlue+i);
        fit2meanPT(kaonboltz[i], kaonmeanpt[i]);

        protonboltz[i] = new TF1(Form("protonboltz%d", i),  "([0]/([1]*(1.007 + [1])))*x*exp(-(sqrt((x*x) + (1.007*1.007))-1.007)*(1/[1]))", 0, 6);
        protonboltz[i]->SetParameters(1.0, 0.34 + i*.015);
        fit2meanPT(protonboltz[i], protonmeanpt[i]);
        protonboltz[i]->SetLineColor(kGreen+i);

        Float_t phi= phiboltz[i]->Integral(2.0, 4.0);
        Float_t totphi = phiboltz[i]->Integral(0.0, 10.0);
        Float_t pi = pionboltz[i]->Integral(2.0, 4.0);
        Float_t k = kaonboltz[i]->Integral(2.0, 4.0);
        Float_t p = protonboltz[i]->Integral(2.0, 4.0);
        Float_t totp = protonboltz[i]->Integral(0.0, 10.0);

        piratios[i] = phi/pi;
        kratios[i] = phi/k;
        pratios[i] = phi/p;
    }
    Float_t normpiratios[7];
    Float_t normkratios[7];
    Float_t normpratios[7];
    for(int i = 0; i< 6; i++){
        normpiratios[i] = piratios[i+1]/piratios[1];
        normkratios[i] = kratios[i+1]/kratios[1];
        normpratios[i] = pratios[i+1]/pratios[1];
    }

    TGraphErrors *piratiograph = new TGraphErrors(6, cutrevbincenters, normpiratios, cutrevbinwidths);
    piratiograph->SetMarkerStyle(22);
    piratiograph->SetMarkerSize(2);
    piratiograph->SetMarkerColor(kBlue+1);
    piratiograph->SetLineColor(kBlue+1);
    piratiograph->GetYaxis()->SetTitle("#phi/#pi Ratio /  80-100%% Ratio");
    piratiograph->GetXaxis()->SetTitle("Multiplicity Pct. (V0A)");
    piratiograph->GetXaxis()->SetRangeUser(0.0, 100.0);

    TGraphErrors *kratiograph = new TGraphErrors(6, cutrevbincenters, normkratios, cutrevbinwidths);
    kratiograph->SetMarkerStyle(23);
    kratiograph->SetMarkerSize(2);
    kratiograph->SetMarkerColor(kRed+1);
    kratiograph->SetLineColor(kRed+1);
    kratiograph->GetYaxis()->SetTitle("#phi/K Ratio / 80-100%% Ratio");
    kratiograph->GetXaxis()->SetTitle("Multiplicity Pct. (V0A)");
    kratiograph->GetXaxis()->SetRangeUser(0.0, 100.0);

    TGraphErrors *pratiograph = new TGraphErrors(6, cutrevbincenters, normpratios, cutrevbinwidths);
    pratiograph->SetMarkerStyle(24);
    pratiograph->SetMarkerSize(2);
    pratiograph->SetMarkerColor(kGreen+1);
    pratiograph->SetLineColor(kGreen+1);
    pratiograph->GetYaxis()->SetTitle("#phi/p Ratio / 80-100%% Ratio");
    pratiograph->GetXaxis()->SetTitle("Multiplicity Pct. (V0A)");
    pratiograph->GetXaxis()->SetRangeUser(0.0, 100.0);


    TH1D* dummyh = new TH1D("dummyh", "dummyh", 10, 0, 100);
    dummyh->GetYaxis()->SetRangeUser(0.25, 1.25);
    dummyh->GetYaxis()->SetTitle("#phi/h Ratio / (60-80%%) Ratio");

    TLegend* ratleg = new TLegend(0.5, 0.7, 0.8, 0.9);
    ratleg->AddEntry(piratiograph, "#phi/#pi", "lp");
    ratleg->AddEntry(kratiograph, "#phi/K", "lp");
    ratleg->AddEntry(pratiograph, "#phi/p", "lp");
    TCanvas* cratios = new TCanvas("cratios", "cratios", 50, 50, 600, 600);
    cratios->cd();
    dummyh->Draw("AXIS");
    piratiograph->Draw("P SAME");
    kratiograph->Draw("P SAME");
    pratiograph->Draw("P SAME");
    ratleg->Draw();


    TF1* philowboltz = new TF1("philowboltz", "([0]/([1]*(1.0195 + [1])))*x*exp(-(sqrt((x*x) + (1.0195*1.0195))-1.0195)*(1/[1]))", 0, 6);
    philowboltz->SetParameters(1.0, 0.493);
    fit2meanPT(philowboltz, 1.242);
    printf("philowboltz 1st moment: %f\n", philowboltz->Moment(1, 0, 10));

    TF1* phihighboltz = new TF1("phihighboltz", "([0]/([1]*(1.0195 + [1])))*x*exp(-(sqrt((x*x) + (1.0195*1.0195))-1.0195)*(1/[1]))", 0, 6);
    phihighboltz->SetParameters(1.0, 0.595);
    fit2meanPT(phihighboltz, 1.437);
    printf("phihighboltz 1st moment: %f\n", phihighboltz->Moment(1, 0, 10));

    TF1* pionlowboltz = new TF1("pionlowboltz", "([0]/([1]*(0.139 + [1])))*x*exp(-(sqrt((x*x) + (0.139*0.139))-0.139)*(1/[1]))", 0, 6);
    pionlowboltz->SetParameters(1.0, 0.235);
    printf("before fit2meanpt, 1st moment = %f\n", pionlowboltz->Moment(1, 0, 10));
    fit2meanPT(pionlowboltz, 0.4705);
    printf("pionlowboltz 1st moment: %f\n", pionlowboltz->Moment(1, 0, 10));

    TF1* pionhighboltz = new TF1("pionhighboltz", "([0]/([1]*(0.139 + [1])))*x*exp(-(sqrt((x*x) + (0.139*0.139))-0.139)*(1/[1]))", 0, 6);
    pionhighboltz->SetParameters(1.0, 0.272);
    fit2meanPT(pionhighboltz, 0.5453);
    printf("pionhighboltz 1st moment: %f\n", pionhighboltz->Moment(1, 0, 10));

    Float_t pionT[7], kaonT[7], protonT[7], phiT[7];
    printf("\nboltzman params\n");
    printf("pion:\n");
    for(int i = 0; i<7; i++){
        pionT[i] = pionboltz[i]->GetParameter(1);
        printf("mult: %2.f-%2.f, <pT>=%f, T=%f\n",  multbinhigh[i], multbinlow[i], pionboltz[i]->Moment(1, 0, 10), pionboltz[i]->GetParameter(1));
    }
    printf("\nkaon:\n");
    for(int i = 0; i<7; i++){
        kaonT[i] = kaonboltz[i]->GetParameter(1);
        printf("mult: %2.f-%2.f, <pT>=%f, T=%f\n",  multbinhigh[i], multbinlow[i], kaonboltz[i]->Moment(1, 0, 10), kaonboltz[i]->GetParameter(1));
    }
    printf("\nproton:\n");
    for(int i = 0; i<7; i++){
        protonT[i] = protonboltz[i]->GetParameter(1);
        printf("mult: %2.f-%2.f, <pT>=%f, T=%f\n",  multbinhigh[i], multbinlow[i], protonboltz[i]->Moment(1, 0, 10), protonboltz[i]->GetParameter(1));
    }
    printf("\nphi:\n");
    for(int i = 0; i<7; i++){
        phiT[i] = phiboltz[i]->GetParameter(1);
        printf("mult: %2.f-%2.f, <pT>=%f, T=%f\n",  multbinhigh[i], multbinlow[i], phiboltz[i]->Moment(1, 0, 10), phiboltz[i]->GetParameter(1));
    }
    TGraphErrors *pionTgraph = new TGraphErrors(7, revbincenters, pionT, revbinwidths);
    pionTgraph->SetMarkerStyle(112);
    pionTgraph->SetMarkerSize(2);
    pionTgraph->SetMarkerColor(kBlack);
    pionTgraph->SetLineColor(kBlack);
    pionTgraph->GetYaxis()->SetTitle("Temperature (GeV)");
    pionTgraph->GetXaxis()->SetTitle("Multiplicity Pct. (V0A)");
    pionTgraph->GetXaxis()->SetRangeUser(0.0, 100.0);

    TGraphErrors *kaonTgraph = new TGraphErrors(7, revbincenters, kaonT, revbinwidths);
    kaonTgraph->SetMarkerStyle(34);
    kaonTgraph->SetMarkerSize(2);
    kaonTgraph->SetMarkerColor(kSpring-1);
    kaonTgraph->SetLineColor(kSpring-1);
    kaonTgraph->GetYaxis()->SetTitle("Temperature (GeV)");
    kaonTgraph->GetXaxis()->SetTitle("Multiplicity Pct. (V0A)");
    kaonTgraph->GetXaxis()->SetRangeUser(0.0, 100.0);
  
    TGraphErrors *protonTgraph = new TGraphErrors(7, revbincenters, protonT, revbinwidths);
    protonTgraph->SetMarkerStyle(33);
    protonTgraph->SetMarkerSize(2);
    protonTgraph->SetMarkerColor(kOrange-3);
    protonTgraph->SetLineColor(kOrange-3);
    protonTgraph->GetYaxis()->SetTitle("Temperature (GeV)");
    protonTgraph->GetXaxis()->SetTitle("Multiplicity Pct. (V0A)");
    protonTgraph->GetXaxis()->SetRangeUser(0.0, 100.0);

    TGraphErrors *phiTgraph = new TGraphErrors(7, revbincenters, phiT, revbinwidths);
    phiTgraph->SetMarkerStyle(20);
    phiTgraph->SetMarkerSize(2);
    phiTgraph->SetMarkerColor(kRed+1);
    phiTgraph->SetLineColor(kRed+1);
    phiTgraph->GetYaxis()->SetTitle("Temperature (GeV)");
    phiTgraph->GetXaxis()->SetTitle("Multiplicity Pct. (V0A)");
    phiTgraph->GetXaxis()->SetRangeUser(0.0, 100.0);


    philowboltz->SetParameter(0, 1.0/philowboltz->Integral(0, 100));
    phihighboltz->SetParameter(0, 1.0/phihighboltz->Integral(0, 100));
    pionlowboltz->SetParameter(0, 1.0/pionlowboltz->Integral(0, 100));
    pionhighboltz->SetParameter(0, 1.0/pionhighboltz->Integral(0, 100));

    TCanvas* ctemp = new TCanvas("ctemp", "ctemp", 50, 50, 600, 600);
    ctemp->cd();
    pionTgraph->GetYaxis()->SetRangeUser(0, 0.8);
    pionTgraph->GetXaxis()->SetRangeUser(0, 100);
    pionTgraph->Draw("AP");
    kaonTgraph->Draw("P SAME");
    protonTgraph->Draw("P SAME");
    phiTgraph->Draw("P SAME");

    TCanvas* cphiboltz = new TCanvas("cphiboltz", "cphiboltz", 50, 50, 600, 600);
    cphiboltz->cd();
    philowboltz->SetLineColor(kGreen+1);
    phihighboltz->SetLineColor(kRed+1);
    philowboltz->Draw();
    phihighboltz->Draw("SAME");
    for(int i=0; i<7; i++){
        phiboltz[i]->Draw("SAME");
    }

    TCanvas* cpionboltz = new TCanvas("cpionboltz", "cpionboltz", 50, 50, 600, 600);
    cpionboltz->cd();
    pionlowboltz->SetLineColor(kGreen+1);
    pionhighboltz->SetLineColor(kRed+1);
    pionlowboltz->Draw();
    pionhighboltz->Draw("SAME");
    for(int i = 0; i<7; i++){
        pionboltz[i]->Draw("SAME");
    }

    TCanvas* ckaonboltz = new TCanvas("ckboltz", "ckboltz", 50, 50, 600, 600);
    ckaonboltz->cd();
    kaonboltz[0]->Draw();
    for(int i=1; i< 7; i++){
        kaonboltz[i]->Draw("SAME");
    }


    TCanvas* cprotonboltz = new TCanvas("cpboltz", "cpboltz", 50, 50, 600, 600);
    cprotonboltz->cd();
    protonboltz[0]->Draw();
    for(int i=1; i< 7; i++){
        protonboltz[i]->Draw("SAME");
    }

    printf("\npercent yield in 2-4 region\n");
    printf("pion low: %f, pion high: %f\n", pionlowboltz->Integral(2.0, 4.0), pionhighboltz->Integral(2.0, 4.0));
    printf("phi low: %f, phi high: %f\n", philowboltz->Integral(2.0, 4.0), phihighboltz->Integral(2.0, 4.0));
    printf("\nphi/pion ratio in 2-4 region\n");
    printf("low: %f, high: %f, %% change: %2.f%%\n", philowboltz->Integral(2.0, 4.0)/pionlowboltz->Integral(2.0, 4.0), phihighboltz->Integral(2.0, 4.0)/pionhighboltz->Integral(2.0, 4.0), ((phihighboltz->Integral(2.0, 4.0)/pionhighboltz->Integral(2.0, 4.0))/(philowboltz->Integral(2.0, 4.0)/pionlowboltz->Integral(2.0, 4.0)) - 1.0)*100.0);
  

    
    ///////////////////////

    /// pion /////
    //////////////////////////

    cout<< "////////////////////////////////"<<endl;
    cout<< "/////////////PION///////////////////"<<endl;
    cout<< "////////////////////////////////"<<endl;

    c2->cd(3);
    TH2F *histo_new2 = new TH2F("histo_new2","amplituds of fit",60,0,6,100,0,1.6);
    histo_new2->Draw();

    Char_t  fit_shift_pion[256];
    Char_t  fit_shiftint_pion[256];



    y0= 1;
    y1= 0.235;
    //sprintf(text9,"([0]/([1]*(1.5195+[1])))*x*exp(- (sqrt((x*x)+ (1.5195*1.5195))-1.5195) * (1/[1])  )");

    sprintf(fit_shift_pion,"(%f/(%f*(0.139+%f)))*x*exp(- (sqrt((x*x)+ (0.139*0.139))-0.139) * (1/%f)  )",y0,y1,y1,y1);

    TF1 *bw1_shift_pion = new TF1("bw1_shift_pion",fit_shift_pion,0,6);
    bw1_shift_pion->SetLineColor(4);   
    //bw1->Integral(0,100);
    bw1_shift_pion->Draw("same");

    cout<< " 1 integral all ----->"<< bw1_shift_pion->Integral(0,100)<<endl;
    cout<< " 1 integral 2-4 ----->"<< bw1_shift_pion->Integral(2,4)<<endl;



    // mean pt and 

    sprintf(fit_shiftint_pion,"(%f/(%f*(0.139+%f)))*x*x*exp(- (sqrt((x*x)+ (0.139*0.139))-0.139) * (1/%f)  )",y0,y1,y1,y1);
    TF1 *bw1_shiftint_pion = new TF1("bw1_shiftint_pion",fit_shiftint_pion,0,1000);
    cout<<"mean pt  --> "<<  bw1_shiftint_pion->Integral(0,1000)/ bw1_shift_pion->Integral(0,1000) <<endl;

    pionyield->SetBinContent(1, bw1_shift_pion->Integral(2,4));
   
    /////////////////////////////////////////////////////
    cout<< "----------------------------------------->"<<endl;

    Char_t  fit_shift2_pion[256];
    Char_t  fit_shift2int_pion[256];

    y0= 1;
    y1= 0.272;
    sprintf(fit_shift2_pion,"(%f/(%f*(0.139+%f)))*x*exp(- (sqrt((x*x)+ (0.139*0.139))-0.139) * (1/%f)  )",y0,y1,y1,y1);

    TF1 *bw1_shift2_pion = new TF1("bw1_shift2_pion",fit_shift2_pion,0,6);
    bw1_shift2_pion->SetLineColor(3);   
    //bw1->Integral(0,100);
    bw1_shift2_pion->Draw("same");

    cout<< " 2 integral all  ----->"<< bw1_shift2_pion->Integral(0,100)<<endl;
    cout<< " 2 integral  2-4 ----->"<< bw1_shift2_pion->Integral(2,4)<<endl;


    sprintf(fit_shift2int_pion,"(%f/(%f*(0.139+%f)))*x*x*exp(- (sqrt((x*x)+ (0.139*0.139))-0.139) * (1/%f)  )",y0,y1,y1,y1);
    TF1 *bw1_shift2int_pion = new TF1("bw1_shift2int_pion",fit_shift2int_pion,0,1000);
    cout<<"mean pt  --> "<<  bw1_shift2int_pion->Integral(0,1000)/ bw1_shift2_pion->Integral(0,1000) <<endl;

    pionyield->SetBinContent(2, bw1_shift2_pion->Integral(2,4));
    c2->cd(4)->SetLogy();
    TH2F *histo_new4 = new TH2F("histo_new4","amplituds of fit",60,0,6,100,0.000001,5);

    histo_new4->Draw();


    bw1_shift_pion->Draw("same");
    bw1_shift2_pion->Draw("same");


    /////


    cout<< "////////////////////////////////"<<endl;
    cout<< "/////////////KAON///////////////////"<<endl;
    cout<< "////////////////////////////////"<<endl;

    c2->cd(5);
    TH2F *histo_new3again = new TH2F("histo_new3again","amplituds of fit",60,0,6,100,0,1.6);
    histo_new3again->Draw();

    Char_t  fit_shift_kaon[256];
    Char_t  fit_shiftint_kaon[256];



    y0= 1;
    y1= 0.343;
    //sprintf(text9,"([0]/([1]*(1.5195+[1])))*x*exp(- (sqrt((x*x)+ (1.5195*1.5195))-1.5195) * (1/[1])  )");

    sprintf(fit_shift_kaon,"(%f/(%f*(0.495+%f)))*x*exp(- (sqrt((x*x)+ (0.495*0.495))-0.495) * (1/%f)  )",y0,y1,y1,y1);

    TF1 *bw1_shift_kaon = new TF1("bw1_shift_kaon",fit_shift_kaon,0,6);
    bw1_shift_kaon->SetLineColor(4);   
    //bw1->Integral(0,100);
    bw1_shift_kaon->Draw("same");

    cout<< " 1 integral all ----->"<< bw1_shift_kaon->Integral(0,100)<<endl;
    cout<< " 1 integral 2-4 ----->"<< bw1_shift_kaon->Integral(2,4)<<endl;



    // mean pt and 

    sprintf(fit_shiftint_kaon,"(%f/(%f*(0.495+%f)))*x*x*exp(- (sqrt((x*x)+ (0.495*0.495))-0.495) * (1/%f)  )",y0,y1,y1,y1);
    TF1 *bw1_shiftint_kaon = new TF1("bw1_shiftint_kaon",fit_shiftint_kaon,0,1000);
    cout<<"mean pt  --> "<<  bw1_shiftint_kaon->Integral(0,1000)/ bw1_shift_kaon->Integral(0,1000) <<endl;

    /////////////////////////////////////////////////////
    cout<< "----------------------------------------->"<<endl;

    Char_t  fit_shift2_kaon[256];
    Char_t  fit_shift2int_kaon[256];

    y0= 1;
    y1= 0.417;
    sprintf(fit_shift2_kaon,"(%f/(%f*(0.495+%f)))*x*exp(- (sqrt((x*x)+ (0.495*0.495))-0.495) * (1/%f)  )",y0,y1,y1,y1);

    TF1 *bw1_shift2_kaon = new TF1("bw1_shift2_kaon",fit_shift2_kaon,0,6);
    bw1_shift2_kaon->SetLineColor(kGreen+1);   
    //bw1->Integral(0,100);
    bw1_shift2_kaon->Draw("same");

    cout<< " 2 integral all  ----->"<< bw1_shift2_kaon->Integral(0,100)<<endl;
    cout<< " 2 integral  2-4 ----->"<< bw1_shift2_kaon->Integral(2,4)<<endl;


    sprintf(fit_shift2int_kaon,"(%f/(%f*(0.495+%f)))*x*x*exp(- (sqrt((x*x)+ (0.495*0.495))-0.495) * (1/%f)  )",y0,y1,y1,y1);
    TF1 *bw1_shift2int_kaon = new TF1("bw1_shift2int_kaon",fit_shift2int_kaon,0,1000);
    cout<<"mean pt  --> "<<  bw1_shift2int_kaon->Integral(0,1000)/ bw1_shift2_kaon->Integral(0,1000) <<endl;


    c2->cd(6)->SetLogy();
    TH2F *histo_new4again = new TH2F("histo_new4again","amplituds of fit",60,0,6,100,0.000001,5);

    histo_new4again->Draw();
    bw1_shift_kaon->Draw("same");
    bw1_shift2_kaon->Draw("same");


    /////

    //cout<< " integral "<< ex6->Integral(0,100)<<endl;
    //cout<< " integral "<< ex6->Integral(0,0.92)<<endl;

    //cout<< " integral error "<< bw1->Integral(0,100)<<endl;
    //cout<< " integral error "<< bw1->Integral(0,0.89)<<endl;

    
    //Draw canvas with phi and pion yields for low mult (left) and high mult (right)
    TCanvas* cyields = new TCanvas("cyields", "cyields", 50, 50, 600, 600);
    cyields->cd();
    pionyield->SetTitle("BW Yields of Particles in 2-4 GeV/c p_{T}");
    pionyield->Draw("P");
    phiyield->Draw("SAME P");

    // =======================================================================
    // Loading phi spectra from HEPData files for fit comparison
    // =======================================================================

    TFile* philowmultFile = new TFile("~/alirepos/utaustin/phiCorrelations/results/paperresults/HEPData-ins1418181-v1-Table_14_mult_60_80.root");
    
    TDirectory* lowdir = (TDirectory*)philowmultFile->GetDirectory("Table 14");
    
    TH1D* philowmult = (TH1D*)(lowdir->Get("Hist1D_y1")->Clone("philowmult"));
    TH1D* philowerr = (TH1D*)(lowdir->Get("Hist1D_y1_e1")->Clone("philowerr"));
    for(int ibin = 1; ibin < philowmult->GetXaxis()->GetNbins(); ibin++){
        philowmult->SetBinError(ibin, philowerr->GetBinContent(ibin));
    }
    
    printf("everything loaded\n");
    philowmult->Scale(phiboltz[1]->Integral(1.0, 2.0)/philowmult->Integral(philowmult->GetXaxis()->FindBin(1.00), philowmult->GetXaxis()->FindBin(1.999), "binwidth"));

    TCanvas* ccompare = new TCanvas("ccomparelow", "ccomparelow", 50, 50, 600, 600);
    ccompare->cd();
    philowmult->GetYaxis()->SetRangeUser(0.0, 0.8);
    philowmult->GetXaxis()->SetRangeUser(0.0, 6.0);
    philowmult->GetXaxis()->SetTitle("p_{T} (GeV/c)");
    philowmult->GetYaxis()->SetTitle("Arb. Units (GeV/c)^{-1}");
    philowmult->Draw("P SAME");
    phiboltz[1]->Draw("SAME");
    //bw1_shift->Draw("SAME");

    TString phifit("([0]/([1]*(1.0195+[1])))*x*exp(-(sqrt((x*x)+(1.0195*1.0195))-1.0195)*(1/[1]))");
    TF1* bw1_lowfit = new TF1("bw1_lowfit",phifit.Data(),0.4,2.0);
    bw1_lowfit->SetParameters(1,0.4);

    TF1* bw1_lowfitwide = new TF1("bw1_lowfitwide",phifit.Data(),0.4,4.0);
    bw1_lowfitwide->SetParameters(1,0.4);


    //philowmult->Fit(bw1_lowfit, "REI");

    TF1* bw1_lowfitfull = new TF1("bw1_lowfitfull", phifit.Data(),0.0,10.0);
    bw1_lowfitfull->SetParameters(bw1_lowfit->GetParameter(0), bw1_lowfit->GetParameter(1));
    bw1_lowfitfull->SetParameter(0, bw1_lowfit->GetParameter(0)/bw1_lowfitfull->Integral(0.0, 100.0));

    //philowmult->Fit(bw1_lowfitwide, "REI");

    TF1* bw1_lowfitfullwide = new TF1("bw1_lowfitfullwide", phifit.Data(),0.0,10.0);
    bw1_lowfitfullwide->SetParameters(bw1_lowfitwide->GetParameter(0), bw1_lowfitwide->GetParameter(1));

    bw1_lowfit->SetLineColor(kRed+1);
    bw1_lowfitfull->SetLineColor(kRed+2);
    bw1_lowfitfullwide->SetLineColor(kViolet);
    
    //bw1_lowfitfull->Draw("SAME");
    //bw1_lowfitfullwide->Draw("SAME"); 

    //=============================================//
    // phi blastwave fits

    TF1* bwPhiHigh = new TF1("blast_phi_high", XBlast, 0.0, 10.0, 5);
    TF1* bwPhiLow = new TF1("blast_phi_low", XBlast, 0.0, 10.0, 5);

    //fill initial values for fit parameters
    double bt,bn,bb;
    BlastParameters(bt,bn,bb);
    bwPhiLow->SetParError(0,0.);
    bwPhiLow->SetParameter(0,1.0195);
    bwPhiLow->FixParameter(0,1.0195);//fix the mass of the particle
    bwPhiLow->SetParameter(1,1000.0);
    bwPhiLow->SetParameter(2,bt);
    bwPhiLow->SetParameter(3,bn);
    bwPhiLow->SetParameter(4,bb);
    bwPhiLow->SetLineColor(kMagenta-3);

    philowmult->Fit(bwPhiLow, "RI0");
    //bwPhiLow->Draw("SAME");

    bwPhiHigh->SetParError(0,0.);
    bwPhiHigh->SetParameter(0,1.0195);
    bwPhiHigh->FixParameter(0,1.0195);//fix the mass of the particle
    bwPhiHigh->SetParameter(1,1000.0);
    bwPhiHigh->SetParameter(2,bt);
    bwPhiHigh->SetParameter(3,bn);
    bwPhiHigh->SetParameter(4,bb);
    bwPhiHigh->SetLineColor(kAzure-2);  

    


    TFile* phihighmultFile = new TFile("~/alirepos/utaustin/phiCorrelations/results/paperresults/HEPData-ins1418181-v1-Table_9_mult_0_5.root");
    
    TDirectory* highdir = (TDirectory*)phihighmultFile->GetDirectory("Table 9");
    
    TH1D* phihighmult = (TH1D*)(highdir->Get("Hist1D_y1")->Clone("phihighmult"));
    TH1D* phihigherr = (TH1D*)(highdir->Get("Hist1D_y1_e1")->Clone("phihigherr"));
    for(int ibin = 1; ibin < phihighmult->GetXaxis()->GetNbins(); ibin++){
        phihighmult->SetBinError(ibin, phihigherr->GetBinContent(ibin));
    }

    printf("everything loaded\n");
    
    phihighmult->Scale(phiboltz[6]->Integral(1.0, 2.0)/phihighmult->Integral(phihighmult->GetXaxis()->FindBin(1.001), phihighmult->GetXaxis()->FindBin(1.9999), "binwidth"));

    
    TCanvas* ccomparehigh = new TCanvas("ccomparehigh", "ccomparehigh", 50, 50, 600, 600);
    ccomparehigh->cd();
    phihighmult->GetYaxis()->SetRangeUser(0.0, 0.8);
    phihighmult->GetXaxis()->SetRangeUser(0.0, 6.0);
    phihighmult->GetXaxis()->SetTitle("p_{T} (GeV/c)");
    phihighmult->GetYaxis()->SetTitle("Arb. Units (GeV/c)^{-1}");
    //phihighmult->Fit(bwPhiHigh, "RI0");
    phihighmult->Draw("P SAME");
    phiboltz[6]->Draw("SAME");
    //bw1_shift2->Draw("SAME");
    
    //copy phi blastwaves to normalize them to 1
    TF1* bwPhiHighScaled = new TF1("bwPhiHighScaled", XBlast, 0.0, 10.0, 5);
    bwPhiHighScaled->SetLineColor(kMagenta-3);
    TF1* bwPhiLowScaled = new TF1("bwPhiLowScaled", XBlast, 0.0, 10.0, 5);
    bwPhiLowScaled->SetLineColor(kAzure-2);
    for(int ipar = 0; ipar < 5; ipar++){
        bwPhiHighScaled->SetParameter(ipar, bwPhiHigh->GetParameter(ipar));
        bwPhiLowScaled->SetParameter(ipar, bwPhiLow->GetParameter(ipar));
    }
    Double_t bwphihighint = bwPhiHighScaled->Integral(0, 100.0);
    bwPhiHighScaled->SetParameter(1, bwPhiHighScaled->GetParameter(1)/bwphihighint);
    Double_t bwphilowint = bwPhiLowScaled->Integral(0, 100.0);
    bwPhiLowScaled->SetParameter(1, bwPhiLowScaled->GetParameter(1)/bwphilowint);
    
    // =======================================================================
    // Loading pi,K,p spectra from HEPData files for fit comparison
    // =======================================================================
    
    TFile* pionfile = new TFile("~/alirepos/utaustin/phiCorrelations/results/paperresults/HEPData-ins1415274-v1-Table_1.root");
    TDirectory* piondir = (TDirectory*)pionfile->GetDirectory("Table 1");
    TH1D* pionhighmult = (TH1D*)(piondir->Get("Hist1D_y1")->Clone("pionhighmult"));
    TH1D* pionhigherr = (TH1D*)(piondir->Get("Hist1D_y1_e1")->Clone("pionhighmult"));

    TH1D* pionlowmult = (TH1D*)(piondir->Get("Hist1D_y6")->Clone("pionlowmult"));
    TH1D* pionlowerr = (TH1D*)(piondir->Get("Hist1D_y6_e1")->Clone("pionlowerr"));
    
    for(int ibin = 1; ibin < pionhighmult->GetXaxis()->GetNbins(); ibin++){
        pionhighmult->SetBinError(ibin, pionhigherr->GetBinContent(ibin)*pionhighmult->GetBinCenter(ibin));
        pionhighmult->SetBinContent(ibin, pionhighmult->GetBinContent(ibin)*pionhighmult->GetBinCenter(ibin));
        //pionhighmult->SetBinError(ibin, pionhigherr->GetBinContent(ibin));

        pionlowmult->SetBinError(ibin, pionlowerr->GetBinContent(ibin)*pionlowmult->GetBinCenter(ibin));
        pionlowmult->SetBinContent(ibin, pionlowmult->GetBinContent(ibin)*pionlowmult->GetBinCenter(ibin));
        //pionlowmult->SetBinError(ibin, pionlowerr->GetBinContent(ibin));
    }
    
    pionhighmult->Scale(2.0*TMath::Pi());
    pionlowmult->Scale(2.0*TMath::Pi());

    //pionlowmult->Scale(bw1_shift_pion->Integral(0.1, 0.6)/pionlowmult->Integral(pionlowmult->GetXaxis()->FindBin(0.1001), pionlowmult->GetXaxis()->FindBin(0.5999), "binwidth"));

    pionhighmult->Scale(pionboltz[6]->Integral(0.2, 0.5)/pionhighmult->Integral(pionhighmult->GetXaxis()->FindBin(0.2001), pionhighmult->GetXaxis()->FindBin(0.4999), "binwidth"));


    TCanvas* clowpion = new TCanvas("clowpion", "clowpion", 50, 50, 600, 600);
    clowpion->cd();
    pionlowmult->GetXaxis()->SetRangeUser(0.0, 4.0);
    pionlowmult->Draw("HIST E");
    //bw1_shift_pion->Draw("SAME");

    printf("pion hist: %f, pion boltz: %f\n", pionhighmult->Integral(pionhighmult->GetXaxis()->FindBin(0.1), pionhighmult->GetXaxis()->FindBin(10.0), "binwidth"), pionboltz[6]->Integral(0.1, 10.0));
    TCanvas* chighpion = new TCanvas("chighpion", "chighpion", 50, 50, 600, 600);
    chighpion->cd();
    pionhighmult->GetXaxis()->SetRangeUser(0.0, 6.0);
    pionhighmult->Draw("P E");
    pionboltz[6]->Draw("SAME");
    //bw1_shift2_pion->Draw("SAME");

    TString pionfit("([0]/([1]*(0.139+[1])))*x*exp(-(sqrt((x*x)+(0.139*0.139))-0.139)*(1/[1]))");
    TF1* bw1_highfitpion = new TF1("bw1_highfit",pionfit.Data(),2.0,4.0);
    bw1_highfitpion->SetParameters(1,0.4);

    //pionhighmult->Fit(bw1_highfitpion, "R0");
    //bw1_highfitpion->SetParameter(1, 0.270);
    //bw1_highfitpion->Draw("SAME");

    TString pionfitint("([0]/([1]*(0.139+[1])))*x*x*exp(-(sqrt((x*x)+(0.139*0.139))-0.139)*(1/[1]))");
    TF1* bw1_pionfitint = new TF1("bw1_pionfitint", pionfitint.Data(), 2.0, 4.0);
    bw1_pionfitint->SetParameters(bw1_highfitpion->GetParameter(0), bw1_highfitpion->GetParameter(1));
    Float_t highentries = pionhighmult->Integral(0, pionhighmult->GetXaxis()->GetNbins());
    Float_t hightot = 0.0;
    for(int ibin = 1; ibin < pionhighmult->GetXaxis()->GetNbins(); ibin++){
        hightot += pionhighmult->GetBinContent(ibin)*pionhighmult->GetBinCenter(ibin);    
    }
    Float_t highmean = hightot/highentries;

    printf("pion high mult mean pT = %f\n", highmean);

    //printf("fit mean pT = %f\n", bw1_pionfitint->Integral(0.0, 100.0)/bw1_highfitpion->Integral(0.0, 100.0));

    // levy-tsallis fit (works great for high pT, but doesn't have a strong enough turnover at low pt
    TString levy_tsallis = "[0]*TMath::Power(1 + x/[1], -[2])";
    TF1* levyfit = new TF1("levyfit", levy_tsallis.Data(), 0.1, 10.0);
    levyfit->SetParLimits(1, 0.00001, 100000);
    levyfit->SetParLimits(2, 1.0, 10000000);
    levyfit->SetParameters(10.0, 5.0, 100.0);
    levyfit->SetLineColor(kOrange-1);
    printf("made it to after levy fit is defined\n");
    //pionhighmult->Fit(levyfit, "R0");
    //levyfit->Draw("SAME");

    // blast-wave fit

    double fmin,fmax;//define limits of fit
    fmin = 0.1;
    fmax = 3.0;
    TF1* bwPhi=new TF1("blast_phi",XBlast,fmin,fmax,5);
    TF1* bwPionHigh = new TF1("blast_pion_high", XBlast, fmin, fmax, 5);
    TF1* bwPionLow = new TF1("blast_pion_low", XBlast, fmin, fmax, 5);
    TF1* bwProton = new TF1("blast_proton", XBlast, fmin, fmax, 5);

    
    //get integral of histogram, used to set initial scale
    double hint=0.;
    for(int pb=1;pb<=pionhighmult->GetNbinsX();pb++) hint+=GBC(pionhighmult,pb)*GBW(pionhighmult,pb);

    bwPionHigh->SetParError(0,0.);
    bwPionHigh->SetParameter(0,0.1396);
    bwPionHigh->FixParameter(0,0.1396);//fix the mass of the particle

    printf("hint for high: %f\n", hint);
    bwPionHigh->SetParError(1,0.1*hint);
    bwPionHigh->SetParameter(1,1000.0*hint);
    //bwPionHigh->FixParameter(1,hint);

    bwPionHigh->SetParError(2,0.1*bt);//temperature
    bwPionHigh->SetParameter(2,bt);
    //bwPionHigh->FixParameter(2,bt);

    bwPionHigh->SetParError(3,0.1*bn);
    bwPionHigh->SetParameter(3,bn);
    //bwPionHigh->FixParameter(3,bn);

    bwPionHigh->SetParError(4,0.1*bb);
    bwPionHigh->SetParameter(4,bb);
    //bwPionHigh->FixParameter(4,bb);

    bwPionHigh->SetLineColor(kViolet-1);

    //pionhighmult->Fit(bwPionHigh, "R0");

    //bwPionHigh->Draw("SAME");

    TF1* bwPionHighfull = new TF1("blast_pion_highfull", XBlast, 0.0, 10.0, 5);
    for(int ipar = 0; ipar < 5; ipar++) bwPionHighfull->SetParameter(ipar, bwPionHigh->GetParameter(ipar));
    bwPionHighfull->SetLineColor(kViolet-1);
    //bwPionHighfull->Draw("SAME");

    clowpion->cd();
    hint = 0;
    for(int pb=1;pb<=pionlowmult->GetNbinsX();pb++) hint+=GBC(pionlowmult,pb)*GBW(pionlowmult,pb);

    bt = 0.169;
    bn = 2.4;
    bb = 0.98;
    bwPionLow->SetParError(0,0.);
    bwPionLow->SetParameter(0,0.1396);
    bwPionLow->FixParameter(0,0.1396);//fix the mass of the particle

    //bwPionLow->SetParError(1,0.1*hint);
    bwPionLow->SetParameter(1,70.0*hint);
    //bwPionLow->FixParameter(1,hint);

    //bwPionLow->SetParError(2,0.1*bt);//temperature
    bwPionLow->SetParameter(2,bt);
    //bwPionLow->FixParameter(2,bt);

    //bwPionLow->SetParError(3,0.1*bn);
    bwPionLow->SetParameter(3,bn);
    //bwPionLow->FixParameter(3,bn);

    //bwPionLow->SetParError(4,0.1*bb);
    bwPionLow->SetParameter(4,bb);
    //bwPionLow->FixParameter(4,bb);

    bwPionLow->SetLineColor(kGreen+3);

    //pionlowmult->Fit(bwPionLow, "R0");

    //bwPionLow->Draw("SAME");
    TF1* bwPionLowfull = new TF1("blast_pion_lowfull", XBlast, 0.0, 10.0, 5);
    for(int ipar = 0; ipar < 5; ipar++) bwPionLowfull->SetParameter(ipar, bwPionLow->GetParameter(ipar));
    bwPionLowfull->SetLineColor(kGreen+1);
    bwPionLowfull->Draw("SAME");

    //Double_t integral = bwPionLowfull->Integral(0.0, 100.0);
    //printf("hint: %f, function int: %f\n", hint, integral);

    //bwPionLowfull->SetParameter(1, bwPionLowfull->GetParameter(1)/integral);
    //Double_t newint = bwPionLowfull->Integral(0.0, 100.0);
    //printf("new func int: %f\n", newint);

    //Double_t highintegral = bwPionHighfull->Integral(0.0, 100.0);
    //bwPionHighfull->SetParameter(1, bwPionHighfull->GetParameter(1)/highintegral);

    TF1* bwPionLowScaled = new TF1("blast_pion_lowscaled", XBlast, 0.0, 10.0, 5);
    bwPionLowScaled->SetLineColor(kGreen+3);
    TF1* bwPionHighScaled = new TF1("blast_pion_highscaled", XBlast, 0.0, 10.0, 5);
    bwPionHighScaled->SetLineColor(kViolet-1);
    for(int ipar = 0; ipar < 5; ipar++){
        bwPionLowScaled->SetParameter(ipar, bwPionLowfull->GetParameter(ipar));
        bwPionHighScaled->SetParameter(ipar, bwPionHighfull->GetParameter(ipar));
    }

    //bwPionLowScaled->SetParameter(1, bwPionLowfull->GetParameter(1)/integral);
    //bwPionHighScaled->SetParameter(1, bwPionHighfull->GetParameter(1)/highintegral);

    //printf("high integral: %f\n", bwPionHighfull->Integral(0.0, 100.0));
    TCanvas* cscaled = new TCanvas("cscaled", "cscaled", 50, 50, 600, 600);
    cscaled->cd();
    bwPionLowScaled->SetRange(0.0, 4.0);
    bwPionLowScaled->Draw();
    bwPionHighScaled->Draw("SAME");

    bwPionLowScaled->GetHistogram()->GetXaxis()->SetTitle("p_{T} (GeV/c)");
    bwPionLowScaled->GetHistogram()->GetYaxis()->SetTitle("k*dN^2/dp_{T}dy (normalized)");
    cscaled->Modified();
    TLegend *pionlegend = new TLegend(0.525, 0.666, 0.868, 0.862);
    pionlegend->AddEntry(bwPionLowScaled, "60-80% mult.", "l");
    pionlegend->AddEntry(bwPionHighScaled, "0-5% mult.", "l");

    pionlegend->Draw();

    Double_t pionmidyieldhigh = bwPionHighScaled->Integral(2.0, 4.0);
    Double_t pionmidyieldlow = bwPionLowScaled->Integral(2.0, 4.0);
    Double_t phimidyieldhigh = bwPhiHighScaled->Integral(2.0, 4.0);
    Double_t phimidyieldlow = bwPhiLowScaled->Integral(2.0, 4.0);

    TH1D* phihighscaled = (TH1D*)phihighmult->Clone("phihighscaled");
    phihighscaled->Scale(bwPhiHighScaled->Integral(0, 100)/phihighscaled->Integral(0,phihighscaled->GetXaxis()->GetNbins(),"binwidth"));
    TH1D* philowscaled = (TH1D*)philowmult->Clone("philowscaled");
    philowscaled->Scale(bwPhiLowScaled->Integral(0, 100)/philowscaled->Integral(0,philowscaled->GetXaxis()->GetNbins(),"binwidth"));
    TH1D* pionhighscaled = (TH1D*)pionhighmult->Clone("pionhighscaled");
    pionhighscaled->Scale(bwPionHighScaled->Integral(0, 100)/pionhighscaled->Integral(0,pionhighscaled->GetXaxis()->GetNbins(),"binwidth"));
    TH1D* pionlowscaled = (TH1D*)pionlowmult->Clone("pionlowscaled");
    pionlowscaled->Scale(bwPionLowScaled->Integral(0, 100)/pionlowscaled->Integral(0,pionlowscaled->GetXaxis()->GetNbins(),"binwidth"));
    Double_t lowhistratio = philowscaled->Integral(philowscaled->GetXaxis()->FindBin(2.01), philowscaled->GetXaxis()->FindBin(3.99), "binwidth")/pionlowscaled->Integral(pionlowscaled->GetXaxis()->FindBin(2.01), pionlowscaled->GetXaxis()->FindBin(3.99), "binwidth");
    Double_t highhistratio = phihighscaled->Integral(phihighscaled->GetXaxis()->FindBin(2.01), phihighscaled->GetXaxis()->FindBin(3.99), "binwidth")/pionhighscaled->Integral(pionhighscaled->GetXaxis()->FindBin(2.01), pionhighscaled->GetXaxis()->FindBin(3.99), "binwidth");

    

    TCanvas* cscaledphi = new TCanvas("cscaledphi", "cscaledphi", 50, 50, 600, 600);
    cscaledphi->cd();
    bwPhiLowScaled->SetRange(0.0, 4.0);
    bwPhiLowScaled->Draw();
    bwPhiHighScaled->Draw("SAME");

    //take ratio of high to low spectrum for phi and pion to see where shift is occuring.
    bwPhiLowScaled->SetRange(0.0, 10.0);
    TH1* bwPhiLowHist = bwPhiLowScaled->CreateHistogram();
    bwPhiLowHist->GetXaxis()->SetRangeUser(0.0, 10.0);
    TH1* bwPhiHighHist = bwPhiHighScaled->CreateHistogram();
    bwPhiHighHist->GetXaxis()->SetRangeUser(0.0, 10.0);
    TH1* bwPhiRatio = (TH1*)bwPhiHighHist->Clone("bwPhiRatio");
    bwPhiRatio->Divide(bwPhiLowHist);

    TCanvas* cphiratio = new TCanvas("cphiratio", "cphiratio", 50, 50, 600, 600);
    cphiratio->cd();
    bwPhiRatio->Draw();

    bwPionLowScaled->SetRange(0.0, 10.0);
    TH1* bwPionLowHist = bwPionLowScaled->CreateHistogram();
    bwPionLowHist->GetXaxis()->SetRangeUser(0.0, 10.0);
    TH1* bwPionHighHist = bwPionHighScaled->CreateHistogram();
    bwPionHighHist->GetXaxis()->SetRangeUser(0.0, 10.0);
    TH1* bwPionRatio = (TH1*)bwPionHighHist->Clone("bwPionRatio");
    bwPionRatio->Divide(bwPionLowHist);

    TCanvas* cphihistratio = new TCanvas("cphihistratio", "cphihistratio", 50, 50, 600, 600);
    cphihistratio->cd();
    TH1D* phiratio = (TH1D*)phihighmult->Clone("phiratio");
    TH1D* philowmorebins = (TH1D*)phihighmult->Clone("phiratio");
    for(int i = 1; i < philowmult->GetXaxis()->GetNbins(); i++){
        philowmorebins->SetBinContent(i, philowmult->GetBinContent(i));
        philowmorebins->SetBinError(i, philowmult->GetBinError(i));
    }
    philowmorebins->SetBinContent(21, 1.0);   
    phiratio->Divide(philowmorebins);
    //phiratio->Scale(bwPhiLow->Integral(0, 100.0)/bwPhiHigh->Integral(0, 100.0));
    phiratio->Draw("HIST");

    TCanvas* cpionratio = new TCanvas("cpionratio", "cpionratio", 50, 50, 600, 600);
    cpionratio->cd();
    bwPionRatio->Draw();

    TCanvas* cpionhistratio = new TCanvas("cpionhistratio", "cpionhistratio", 50, 50, 600, 600);
    cpionhistratio->cd();
    pionhighmult->GetXaxis()->SetRange(0,0);
    pionlowmult->GetXaxis()->SetRange(0,0);
    TH1D* pionratio = (TH1D*)pionhighmult->Clone("pionratio");
    pionratio->Divide(pionlowmult);
    //pionratio->Scale(bwPionLowfull->Integral(0, 10.0)/bwPionHighfull->Integral(0, 10.0));
    pionratio->SetLineColor(kRed+1);
    pionratio->Draw("HIST");

    TCanvas* ccompareratio = new TCanvas("ccompareratio", "ccompareratio", 50, 50, 600, 600);
    ccompareratio->cd();
    phiratio->Draw("HIST");
    pionratio->Draw("SAME HIST");
    
    printf("testing phi ints: %f, %f\n", bwPhiLow->Integral(0, 10.0), bwPhiHigh->Integral(0, 10.0));

    TCanvas* cscaledphihist = new TCanvas("cscaledphihist", "cscaledphihist", 50, 50, 600, 600);
    cscaledphihist->cd();
    phihighscaled->Draw("HIST E");
    philowscaled->SetLineColor(kRed+1);
    philowscaled->Draw("HIST E SAME");

    TCanvas* cscaledpionhist = new TCanvas("cscaledpionhist", "cscaledpionhist", 50, 50, 600, 600);
    cscaledpionhist->cd();
    pionhighscaled->Draw("HIST E");
    pionlowscaled->SetLineColor(kRed+1);
    pionlowscaled->Draw("HIST E SAME");

    /*
    printf("test of all ints, highpion: %f, lowpion: %f, highpion0-7 - 1: %e, highphi: %f, lowphi: %f\n", bwPionHighScaled->Integral(0.0, 100.0), bwPionLowScaled->Integral(0.0, 100.0), bwPionHighScaled->Integral(0.0, 7.0) - 1.0, bwPhiHighScaled->Integral(0.0, 100.0), bwPhiLowScaled->Integral(0.0, 100.0));
    printf("blastwave ints, highpion: %f, lowpion: %f\n", pionmidyieldhigh, pionmidyieldlow);
    printf("boltzman ints, high phi: %f, low phi: %f\n", bw1_shift2->Integral(2,4), bw1_shift->Integral(2,4));
    printf("balswave phi ints, high phi: %f, low phi: %f\n", phimidyieldhigh, phimidyieldlow);
    printf("blastwave ratios, high: %f, low: %f, increase: %2.2f%%\n", bw1_shift2->Integral(2,4)/pionmidyieldhigh, bw1_shift->Integral(2,4)/pionmidyieldlow, ((bw1_shift2->Integral(2,4)/pionmidyieldhigh)/(bw1_shift->Integral(2,4)/pionmidyieldlow)-1.0)*100.0);
    printf("boltzman fit ratios, high: %f, low: %f, increase: %2.2f%%\n", bw1_shift2->Integral(2,4)/pionmidyieldhigh, bw1_lowfitfull->Integral(2,4)/pionmidyieldlow, ((bw1_shift2->Integral(2,4)/pionmidyieldhigh)/(bw1_lowfitfull->Integral(2,4)/pionmidyieldlow)-1.0)*100.0);
    printf("blastwave fit ratios, high: %f, low: %f, increase: %2.2f%%\n", phimidyieldhigh/pionmidyieldhigh, phimidyieldlow/pionmidyieldlow, ((phimidyieldhigh/pionmidyieldhigh)/(phimidyieldlow/pionmidyieldlow)-1.0)*100.0);
    printf("scaled hist ratios, high:%f, low:%f, increase %2.2f%%\n", highhistratio, lowhistratio, (highhistratio/lowhistratio - 1.0)*100.0);
    */

    //Compute and Draw ratio of phi over pion in this 2-4 GeV/c region for low (left) and high (right) mult
    yieldratio->Divide(phiyield, pionyield);
    yieldratio->SetBinContent(1, lowhistratio);
    yieldratio->SetBinContent(2, highhistratio);
    yieldratio->SetLineColor(kRed+1);
    yieldratio->SetLineWidth(2);
    yieldratio->SetMarkerColor(kRed+1);
    yieldratio->SetMarkerSize(2);
    yieldratio->SetMarkerStyle(23);
    yieldratio->GetXaxis()->SetBinLabel(1, "60-80%");
    yieldratio->GetXaxis()->SetBinLabel(2, "0-5%");
    yieldratio->GetYaxis()->SetTitle("#phi/#pi Ratio for p_{T} 2-4 GeV/c");

    yieldratio->SetTitle("BW Yield Ratio of #phi/#pi for p_{T} 2-4 GeV/c");
    TCanvas* cratio = new TCanvas("cratio", "cratio", 50, 50, 600, 600);
    cratio->cd();
    yieldratio->Draw("P");

    TH1D* normyieldratio = (TH1D*)yieldratio->Clone("normyieldratio");
    normyieldratio->SetBinContent(1, lowhistratio/lowhistratio);
    normyieldratio->SetBinContent(2, highhistratio/lowhistratio);
    normyieldratio->GetYaxis()->SetTitle("Normalized #phi/#pi Ratio for p_{T} 2-4 GeV/c");
    
    TCanvas* cnormratio = new TCanvas("cnormratio", "cnormratio", 50, 50, 600, 600);
    cnormratio->cd();
    normyieldratio->Draw("P");

    //set-up blastwaves with the paper fit values
    TF1* bwPionLowPaper = new TF1("bwPionLowPaper", XBlast, 0.0, 10.0, 5);
    TF1* bwKLowPaper = new TF1("bwKLowPaper", XBlast, 0.0, 10.0, 5);
    TF1* bwpLowPaper = new TF1("bwpLowPaper", XBlast, 0.0, 10.0, 5);
    TF1* bwPhiLowPaper = new TF1("bwPhiLowPaper", XBlast, 0.0, 10.0, 5);

    TF1* bwPionHighPaper = new TF1("bwPionHighPaper", XBlast, 0.0, 10.0, 5);
    TF1* bwKHighPaper = new TF1("bwKHighPaper", XBlast, 0.0, 10.0, 5);
    TF1* bwpHighPaper = new TF1("bwpHighPaper", XBlast, 0.0, 10.0, 5);
    TF1* bwPhiHighPaper = new TF1("bwPhiHighPaper", XBlast, 0.0, 10.0, 5);

    Double_t paperlowT = 0.169, paperlown = 2.4, paperlowB = (paperlown+2.0)*0.36/2.0;
    Double_t paperhighT = 0.143, paperhighn = 1.07, paperhighB = (paperhighn+2.0)*0.547/2.0;
    bwPionLowPaper->SetParameters(0.1396, 1.0, paperlowT, paperlown, paperlowB);
    bwPionHighPaper->SetParameters(0.1396, 1.0, paperhighT, paperhighn, paperhighB);
    bwPhiLowPaper->SetParameters(1.0195, 1.0, paperlowT, paperlown, paperlowB);
    bwPhiHighPaper->SetParameters(1.0195, 1.0, paperhighT, paperhighn, paperhighB);

    //scale fits to paper spectra in the (0.5 - 1.0) pT region same as paper fits
    //bwPionLowPaper->SetParameter(1, pionlowmult->Integral(pionlowmult->GetXaxis()->FindBin(0.5+0.0001), pionlowmult->GetXaxis()->FindBin(1.0-0.0001), "binwidth")/bwPionLowPaper->Integral(0.5, 1.0));
    //bwPionHighPaper->SetParameter(1, pionhighmult->Integral(pionhighmult->GetXaxis()->FindBin(0.5+0.0001), pionhighmult->GetXaxis()->FindBin(1.0-0.0001), "binwidth")/bwPionHighPaper->Integral(0.5, 1.0));
    //bwPhiLowPaper->SetParameter(1, philowmult->Integral(philowmult->GetXaxis()->FindBin(0.6+0.0001), philowmult->GetXaxis()->FindBin(3.0-0.0001), "binwidth")/bwPhiLowPaper->Integral(0.6, 3.0));
    //bwPhiHighPaper->SetParameter(1, phihighmult->Integral(phihighmult->GetXaxis()->FindBin(0.6+0.0001), phihighmult->GetXaxis()->FindBin(3.0-0.0001), "binwidth")/bwPhiHighPaper->Integral(0.6, 3.0));
    
    /*
    bwPionLowPaper->SetParameter(1, 1.0/bwPionLowPaper->Integral(0.0, 100.0));
    bwPionHighPaper->SetParameter(1, 1.0/bwPionHighPaper->Integral(0.0, 100.0));
    bwPhiLowPaper->SetParameter(1, 1.0/bwPhiLowPaper->Integral(0.0, 100.0));
    bwPhiHighPaper->SetParameter(1, 1.0/bwPhiHighPaper->Integral(0.0, 100.0));
    */

    TF1* bwPionLowXPaper = new TF1("bwPionLowXPaper", XXBlast, 0.0, 10.0, 5);
    TF1* bwPionHighXPaper = new TF1("bwPionHighXPaper", XXBlast, 0.0, 10.0, 5);
    TF1* bwPhiLowXPaper = new TF1("bwPhiLowXPaper", XXBlast, 0.0, 10.0, 5);
    TF1* bwPhiHighXPaper = new TF1("bwPhiHighXPaper", XXBlast, 0.0, 10.0, 5);
    for(int i = 0; i < 5; i++){
        bwPionLowXPaper->SetParameter(i, bwPionLowPaper->GetParameter(i));
        bwPionHighXPaper->SetParameter(i, bwPionHighPaper->GetParameter(i));
        bwPhiLowXPaper->SetParameter(i, bwPhiLowPaper->GetParameter(i));
        bwPhiHighXPaper->SetParameter(i, bwPhiHighPaper->GetParameter(i));
    }

    /*
    printf("paper stuff\n ==========================\n");
    printf("total pion ints - low: %f, high %f\n", bwPionLowPaper->Integral(0, 100.0), bwPionHighPaper->Integral(0, 100.0));
    printf("total phi ints - low: %f, high %f\n", bwPhiLowPaper->Integral(0, 100.0), bwPhiHighPaper->Integral(0, 100.0));
    printf("mid pion ints - low: %f, high %f\n", bwPionLowPaper->Integral(2.0, 4.0), bwPionHighPaper->Integral(2.0, 4.0));
    printf("mid phi ints - low: %f, high %f\n", bwPhiLowPaper->Integral(2.0, 4.0), bwPhiHighPaper->Integral(2.0, 4.0));
    printf("low ratio: %f, high ratio: %f, increase: %2.2f%%\n", bwPhiLowPaper->Integral(2.0, 4.0)/bwPionLowPaper->Integral(2.0, 4.0), bwPhiHighPaper->Integral(2.0, 4.0)/bwPionHighPaper->Integral(2.0, 4.0), ((bwPhiHighPaper->Integral(2.0, 4.0)/bwPionHighPaper->Integral(2.0, 4.0))/(bwPhiLowPaper->Integral(2.0, 4.0)/bwPionLowPaper->Integral(2.0, 4.0)) - 1.0)*100.0); 
    */ 
    bwPionLowPaper->SetLineColor(kBlue+1);
    bwPionHighPaper->SetLineColor(kGreen+1);
    bwPhiLowPaper->SetLineColor(kRed+1);
    bwPhiHighPaper->SetLineColor(kViolet+1);

    TCanvas* cpaperpionfits = new TCanvas("cpaperpionfits", "cpaperpionfits", 50, 50, 600, 600);
    cpaperpionfits->cd();
    bwPionLowPaper->SetRange(0.0, 6.0);
    bwPionLowPaper->Draw();
    pionlowmult->Draw("SAME HIST E");
    //bwPionHighPaper->Draw("SAME");
//    pionhighscaled->Draw("SAME HIST E");


    TCanvas* cpaperphifits = new TCanvas("cpaperphifits", "cpaperphifits", 50, 50, 600, 600);
    cpaperphifits->cd();
    bwPhiLowPaper->SetRange(0.0, 6.0);
    bwPhiLowPaper->Draw();
    bwPhiHighPaper->Draw("SAME");

    //compare paper values to paper spectrum (normalized to 1 vs norm to 0.5-1.0 pT?)
    TCanvas* ccomppaperpionlow = new TCanvas("ccomppaperpionlow", "ccomppaperpionlow", 50, 50, 600, 600);
    ccomppaperpionlow->cd();
    //pionlowmult->Draw("HIST E");
    bwPionLowPaper->Draw();
    bwPionHighPaper->Draw("SAME");
    //bwPionLowfull->Draw("SAME");

    TCanvas* ccomppaperpionhigh = new TCanvas("ccomppaperpionhigh", "ccomppaperpionhigh", 50, 50, 600, 600);
    ccomppaperpionhigh->cd();
    pionhighmult->Draw("HIST E");
    bwPionHighPaper->Draw("SAME");
    //bwPionHighfull->Draw("SAME");

    TCanvas* ccomppaperphilow = new TCanvas("ccomppaperphilow", "ccomppaperphilow", 50, 50, 600, 600);
    ccomppaperphilow->cd();
    //philowmult->Draw("HIST E");
    bwPhiLowPaper->Draw();
    bwPhiHighPaper->Draw("SAME");
    //bwPhiLow->Draw("SAME");

    TCanvas* ccomppaperphihigh = new TCanvas("ccomppaperphihigh", "ccomppaperphihigh", 50, 50, 600, 600);
    ccomppaperphihigh->cd();
    phihighmult->Draw("HIST E");
    bwPhiHighPaper->Draw("SAME");
    //bwPhiHigh->Draw("SAME");

    /*
    printf("ratio of paper blast wave fit to actual spectra:\n");
    printf("pion:\n");
    printf("low mult => 0.5-1.0: %f, 2-4: %f\n", bwPionLowPaper->Integral(0.5, 1.0)/pionlowmult->Integral(pionlowmult->GetXaxis()->FindBin(0.5+0.0001), pionlowmult->GetXaxis()->FindBin(1.0-0.0001), "binwidth"), bwPionLowPaper->Integral(2.0, 4.0)/pionlowmult->Integral(pionlowmult->GetXaxis()->FindBin(2.0+0.0001), pionlowmult->GetXaxis()->FindBin(4.0-0.0001), "binwidth"));
    printf("high mult => 0.5-1.0: %f, 2-4: %f\n", bwPionHighPaper->Integral(0.5, 1.0)/pionhighmult->Integral(pionhighmult->GetXaxis()->FindBin(0.5+0.0001), pionhighmult->GetXaxis()->FindBin(1.0-0.0001), "binwidth"), bwPionHighPaper->Integral(2.0, 4.0)/pionhighmult->Integral(pionhighmult->GetXaxis()->FindBin(2.0+0.0001), pionhighmult->GetXaxis()->FindBin(4.0-0.0001), "binwidth"));
    printf("phi:\n");
    printf("low mult => 0.6-3.0: %f, 2-4: %f\n", bwPhiLowPaper->Integral(0.6, 3.0)/philowmult->Integral(philowmult->GetXaxis()->FindBin(0.6+0.0001), philowmult->GetXaxis()->FindBin(3.0-0.0001), "binwidth"), bwPhiLowPaper->Integral(2.0, 4.0)/philowmult->Integral(philowmult->GetXaxis()->FindBin(2.0+0.0001), philowmult->GetXaxis()->FindBin(4.0-0.0001), "binwidth"));
    printf("high mult => 0.6-3.0: %f, 2-4: %f\n", bwPhiHighPaper->Integral(0.6, 3.0)/phihighmult->Integral(phihighmult->GetXaxis()->FindBin(0.6+0.0001), phihighmult->GetXaxis()->FindBin(3.0-0.0001), "binwidth"), bwPhiHighPaper->Integral(2.0, 4.0)/phihighmult->Integral(phihighmult->GetXaxis()->FindBin(2.0+0.0001), phihighmult->GetXaxis()->FindBin(4.0-0.0001), "binwidth"));

    printf("parameters of fits:\n");
    printf("phi low mult: T = %f, <B_T> = %f, n=%f\n", bwPhiLowPaper->GetParameter(2), bwPhiLowPaper->GetParameter(4)*2.0/(bwPhiLowPaper->GetParameter(3)+2.0), bwPhiLowPaper->GetParameter(3));
    printf("phi high mult: T = %f, <B_T> = %f, n=%f\n", bwPhiHighPaper->GetParameter(2), bwPhiHighPaper->GetParameter(4)*2.0/(bwPhiHighPaper->GetParameter(3)+2.0), bwPhiHighPaper->GetParameter(3));

    printf("pion low mult: T = %f, <B_T> = %f, n=%f\n", bwPionLowPaper->GetParameter(2), bwPionLowPaper->GetParameter(4)*2.0/(bwPionLowPaper->GetParameter(3)+2.0), bwPionLowPaper->GetParameter(3));
    printf("pion high mult: T = %f, <B_T> = %f, n=%f\n", bwPionHighPaper->GetParameter(2), bwPionHighPaper->GetParameter(4)*2.0/(bwPionHighPaper->GetParameter(3)+2.0), bwPionHighPaper->GetParameter(3));
    */

    //printf("mean pion pT 60-80%%: %f\n", bwPionLowXPaper->Integral(0, 10.0));
    //printf("mean pion pT 0-5%%: %f\n", bwPionHighXPaper->Integral(0, 10.0));
    //printf("mean phi pT 60-80%%: %f\n", bwPhiLowXPaper->Integral(0, 10.0));
    //printf("mean phi pT 0-5%%: %f\n", bwPhiHighXPaper->Integral(0, 10.0));
}
