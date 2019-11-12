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

    gStyle->SetOptFit(1);
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
  gStyle->SetLabelFont(22,"x");
  gStyle->SetLabelFont(22,"y");
  gStyle->SetLabelFont(22,"z");

   c->cd();
   c->Modified();
   c->cd();
     
   TH2F *h = new TH2F("pt","pt",100,0,2.5,100,0,100);
   TH2F *hraw = new TH2F("ptraw","ptraw",100,0,3,100,0,5000);
 TH2F *hraw2 = new TH2F("ptraw2","ptraw2",100,0,3,100,0,3000);

   TH2F *hacc = new TH2F("ptacc","ptacc",100,0,2,100,0,1);

   TH2F *hrawacc = new TH2F("ptrawacc","ptrawacc",100,0,2,100,0,1500);
   TH1F *histo = new TH1F("ptagain","ptagain",10,0,3);

   TH2F *hmtminusml = new TH2F("mt-m0","mt-m0",100,0,2,100,0,1500);
 
//hmtminusml->SetLogy(1);
  
   c->cd(1);  
//c_1->SetLogy(1);
//c_1->SetLogx(1);
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

 
   // h->Draw("");

   TGaxis *g = new TGaxis(0,0,6,0,0,6,6,""); 
   g->SetLabelOffset(h->GetXaxis()->GetLabelOffset());
   g->SetLabelSize(h->GetXaxis()->GetLabelSize());  
   g->SetTickSize(0.03);
   g->SetGridLength(0);
   g->SetTitleOffset(1);
   g->SetLabelOffset(0.005);
   g->SetLabelFont(20);
   g->SetName(" ");
//g->SetTitle(" ");

   //g->Draw();

   TGaxis *g1 = new TGaxis(0,0,0,0,0,3.5,7,""); 
   g1->SetLabelOffset(h->GetXaxis()->GetLabelOffset());
   g1->SetLabelSize(h->GetXaxis()->GetLabelSize());  
   g1->SetTickSize(0.03);
   g1->SetGridLength(0);
   g1->SetTitleOffset(1);
   g1->SetLabelOffset(0.005);
   g1->SetLabelFont(20);
   g1->SetName(" ");
//g->SetTitle(" ");

   //g->Draw();
 
   //c_1->GetFrame()->SetLineWidth(2); 
   //c_1->Modified();
  

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
   //Float_t  pt[5] = {0.4,0.8,1.2,1.6,2.0};
   //Float_t  yield[5] ={217.0,946.0,1049.0,534.0,36.0};
   //Float_t  yielderr[5] = {54.0,103.0,114.0,91.0,35.0};
   //Float_t  yield2[5] ={0.0598377,0.191638,0.33207,0.253979,0.07389};
///////////////////////////////////////////////////////////
///////// all cos 0.9 dca 3 
   //Float_t  pt[5] = {0.4,0.8,1.2,1.6,2.0};
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

       TF1 *ex6 = new TF1("ex6",text6,0,3);
        //ex6->SetParLimits(0,2,4);
       ex6->SetParLimits(1,0.300, 0.400);   
       ex6->SetLineColor(4);       
       graph8->Fit("ex6","RMEI");
      
       cout<< " integral "<< ex6->Integral(0,100)<<endl;
       cout<< " integral "<< ex6->Integral(0,0.92)<<endl;

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

/////

       c2->cd(2)->SetLogy();
      TH2F *histo_new3 = new TH2F("histo_new3","amplituds of fit",60,0,6,100,0.001,1);
     
       histo_new3->Draw();
       bw1_shift->Draw("same");
       bw1_shift2->Draw("same");


///////////////////////////////
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
       bw1_shift2_kaon->SetLineColor(3);   
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


      
}
