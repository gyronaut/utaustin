#include <stdio.h>

void KKInvMass(){
    gStyle->SetOptStat(0);
    gStyle->SetOptFit(0);

    TFile* file = new TFile("~/phiStudies/results_onlineEff/FAST/hphi_0_20_rapiditytest.root");
    TList* list = (TList*) file->Get("phiCorr_mult_0_20_");
    THnSparseF* kkUSDist= (THnSparseF*)list->FindObject("fkkUSDist");
    THnSparseF* kkLSDist= (THnSparseF*)list->FindObject("fkkLSDist");
    THnSparseF* trigDist = (THnSparseF*)list->FindObject("fTrigDist");

    TH1D* events = (TH1D*)list->FindObject("fNevents");
    
    Float_t epsilon = 0.0001;

    kkUSDist->GetAxis(0)->SetRangeUser(2.0 + epsilon, 4.0 - epsilon);
    kkLSDist->GetAxis(0)->SetRangeUser(2.0 + epsilon, 4.0 - epsilon);
    
    TH1D* USInvMass = kkUSDist->Projection(1);
    TH1D* LSInvMass = kkLSDist->Projection(1);
    TH1D* origLSInvMass = (TH1D*)LSInvMass->Clone("origLSInvMass");
    origLSInvMass->SetLineColor(kRed);
    origLSInvMass->SetLineWidth(2);
    
    Double_t sidebandUS = (Double_t)(USInvMass->Integral(USInvMass->GetXaxis()->FindBin(1.04), USInvMass->GetXaxis()->FindBin(1.06)) + USInvMass->Integral(USInvMass->GetXaxis()->FindBin(0.995), USInvMass->GetXaxis()->FindBin(1.005)));
    Double_t sidebandLS = (Double_t)(LSInvMass->Integral(LSInvMass->GetXaxis()->FindBin(1.04), LSInvMass->GetXaxis()->FindBin(1.06)) + LSInvMass->Integral(LSInvMass->GetXaxis()->FindBin(0.995), LSInvMass->GetXaxis()->FindBin(1.005)));
    Double_t scale = sidebandUS/sidebandLS;
    
    LSInvMass->Scale(scale);

    LSInvMass->SetLineColor(kRed);
    LSInvMass->SetLineWidth(2);
    USInvMass->SetLineWidth(3);
    
    USInvMass->SetTitle("");
    USInvMass->GetXaxis()->SetTitle("#it{m}_{KK} (GeV/#it{c}^{2})");
    USInvMass->GetXaxis()->SetTitleSize(0.05);
    USInvMass->SetLineColor(kBlack);
   
    TH1D* corrected = (TH1D*)USInvMass->Clone("corrected");
    corrected->Add(LSInvMass, -1.0);
    corrected->SetLineWidth(2);
    corrected->SetMarkerSize(2);
    corrected->SetMarkerStyle(34);
    corrected->GetXaxis()->SetTitle("#it{m}_{K^{+}K^{-}} (GeV/#it{c}^{2})");
    TF1* fit = new TF1("fit",  "[0]*TMath::Voigt(x - [1], [2], [3], 4) + pol2(4)",0.99, 1.07);
    fit->SetParameter(1, 1.020);
    fit->SetParameter(2, 0.0002);
    fit->SetParameter(0, 600);
    fit->FixParameter(3, 0.00426);
    fit->SetParLimits(1, 1.010, 1.030);
    fit->SetLineColor(kBlue);
    fit->SetLineStyle(7);
    fit->SetLineWidth(6);
    corrected->Fit(fit, "R");

    TF1* voigtFit = new TF1("voigtFit", "[0]*TMath::Voigt(x - [1], [2], [3], 4)", 0.99, 1.07);
    voigtFit->SetParameters(fit->GetParameter(0), fit->GetParameter(1), fit->GetParameter(2), fit->GetParameter(3));

    TF1* bgFit = new TF1("bgFit", "pol2(0)", 0.99, 1.07);
    bgFit->SetParameters(fit->GetParameter(4), fit->GetParameter(5), fit->GetParameter(6));

    Float_t fullInt = voigtFit->Integral(0.99, 1.07);
    Float_t onlyHistInt = corrected->Integral(corrected->GetXaxis()->FindBin(0.9901), corrected->GetXaxis()->FindBin(1.06999), "width");
    Float_t onlyBGInt = bgFit->Integral(0.99, 1.07);
    Float_t fullHistInt = corrected->Integral(corrected->GetXaxis()->FindBin(0.9901), corrected->GetXaxis()->FindBin(1.06999), "width") - bgFit->Integral(0.99, 1.07);
    Float_t widePct = voigtFit->Integral(1.01, 1.03);
    widePct = widePct/fullInt;
    Float_t peakInt = voigtFit->Integral(1.014, 1.026);
    Float_t bgpeakInt = bgFit->Integral(1.014, 1.026);
    Float_t bg2peak = bgpeakInt/peakInt;
    
    printf("=================\n\nfit integral: %e,    hist integral: %e,   only hist: %e,    only BG: %e,    ratio to residual: %4.2f%% \n=====================\n", fullInt, fullHistInt, onlyHistInt, onlyBGInt, bg2peak*100);

    Float_t diffPct[11];
    Float_t diffInt[11];
    Float_t sigBG[11];
    Float_t purity[11];
    TH1D* effHist = new TH1D("effHist", "Efficiency & Purity of Different Mass Ranges;Invariant Mass Range (MeV/c^{2});Percent", 11, 0., 11.);
    TH1D* purityHist = new TH1D("purityHist", "Purity;Invariant Mass Range;#frac{Signal}{Signal + BG}", 11, 0., 11.);

    for(int i = 0; i < 11; i++){
        diffPct[i] = Float_t(voigtFit->Integral(1.008 + 0.001*i, 1.032 - 0.001*i))/fullInt;
        diffInt[i] = Float_t(corrected->Integral(corrected->GetXaxis()->FindBin(1.0081 + 0.001*i), corrected->GetXaxis()->FindBin(1.032 - 0.0001 - 0.001*i), "width") - bgFit->Integral(1.008 + 0.001*i, 1.032 - 0.001*i))/fullHistInt;
        sigBG[i] = Float_t(corrected->Integral(corrected->GetXaxis()->FindBin(1.0081 + 0.001*i), corrected->GetXaxis()->FindBin(1.032 - 0.0001 - 0.001*i), "width") - bgFit->Integral(1.008 + 0.001*i, 1.032 - 0.001*i))/Float_t(LSInvMass->Integral(LSInvMass->GetXaxis()->FindBin(1.0081+0.001*i), LSInvMass->GetXaxis()->FindBin(1.0319 - 0.001*i),"width") + bgFit->Integral(1.008 + 0.001*i, 1.032 - 0.001*i));
        purity[i] = Float_t(corrected->Integral(corrected->GetXaxis()->FindBin(1.0081 + 0.001*i), corrected->GetXaxis()->FindBin(1.032 - 0.0001 - 0.001*i), "width") - bgFit->Integral(1.008 + 0.001*i, 1.032 - 0.001*i))/Float_t(USInvMass->Integral(USInvMass->GetXaxis()->FindBin(1.0081+0.001*i), USInvMass->GetXaxis()->FindBin(1.0319 - 0.001*i),"width") + bgFit->Integral(1.008 + 0.001*i, 1.032 - 0.001*i));
        effHist->SetBinContent(i+1, diffInt[i]*100.0);
        effHist->GetXaxis()->SetBinLabel(i+1, Form("%i", (int)(((1.032 - 0.001*i) - (1.008 + 0.001*i))*1000)));
        purityHist->SetBinContent(i+1, sigBG[i]);
    }


    TCanvas *ceff = new TCanvas("ceff", "ceff", 50, 50, 900, 600);
    ceff->cd();
    effHist->GetYaxis()->SetRangeUser(0.0, 100.0);
    effHist->SetLineWidth(2);
    effHist->Draw();
    ceff->Update();

    Float_t rightmax = 1.1*purityHist->GetMaximum();
    Float_t axisscale = gPad->GetUymax()/rightmax;
    purityHist->SetLineColor(kRed);
    purityHist->SetLineWidth(2);
    purityHist->Scale(axisscale);
    purityHist->Draw("HIST SAME");

    TGaxis *axis = new TGaxis(gPad->GetUxmax(), gPad->GetUymin(), gPad->GetUxmax(), gPad->GetUymax(),0,rightmax,510,"+L");
    axis->SetLineColor(kRed);
    axis->SetTextColor(kRed);
    axis->SetLabelColor(kRed);
    axis->SetLabelSize(0.035);
    axis->SetLabelFont(42);
    axis->SetTitle("#frac{Signal}{BG}");
    axis->Draw();
    ceff->Update();

    for(int j = 0; j < 11; j++){
        printf(" fit range %4.3f to %4.3f: %4.2f%%\n", 1.008+0.001*j, 1.032 -0.001*j, diffPct[j]*100.0);
        printf("hist range %4.3f to %4.3f: %4.2f%%, Signal/(Signal+ BG): %4.2f Signal/BG: %4.2f\n\n", 1.008+0.001*j, 1.032 -0.001*j, diffInt[j]*100.0, purity[j], sigBG[j]);
    }



    TLine* sbLine1 = new TLine(1.04, 0, 1.04, 44000);
    TLine* sbLine2 = new TLine(1.06, 0, 1.06, 44000);

    sbLine1->SetLineStyle(4);
    sbLine2->SetLineStyle(4);
    sbLine1->SetLineColor(kBlack);

    sbLine1->SetLineWidth(2);
    sbLine2->SetLineWidth(2);

    TLine* peak1 = new TLine(1.01, 0, 1.01, 44000);
    TLine* peak2 = new TLine(1.03, 0, 1.03, 44000);

    peak1->SetLineStyle(9);
    peak2->SetLineStyle(9);
    peak1->SetLineColor(kBlack);
    peak2->SetLineColor(kBlack);
    peak1->SetLineWidth(2);
    peak2->SetLineWidth(2);

    TLegend *lineLeg = new TLegend(0.4564, 0.5550, 0.8792, 0.6632);
    lineLeg->AddEntry(sbLine1, "Sideband Region", "l");
    //lineLeg->AddEntry(peak1, "#phi(1020) Peak Region", "l");

    TLegend *leg = new TLegend(0.5435, 0.507, 0.898, 0.6777);
    leg->AddEntry(USInvMass, "US Kaon Pairs", "le");
    leg->AddEntry(LSInvMass, "#splitline{Est. BG}{(using LS Kaon Pairs)}", "le");
    leg->SetLineWidth(0);

    TPaveText *text = new TPaveText(0.5953, 0.7622, 0.9147, 0.9103, "NDC");
    text->AddText("ALICE Preliminary");
    //text->AddText("Work In Progress");
    text->AddText("p-Pb #sqrt{s_{NN}} = 5.02 TeV");
    text->AddText("0-20% Multiplicity Class (V0A)");
    text->SetFillColor(kWhite);
    text->SetFillStyle(0);
    text->SetBorderSize(0);
    text->SetTextFont(42);

    TPaveText *pTText = new TPaveText(0.1589, 0.8031, 0.5134, 0.8902, "NDC");
    pTText->AddText("2.0 < #it{p}_{T}^{KK} < 4.0 GeV/#it{c}");
    pTText->SetFillColor(kWhite);
    pTText->SetBorderSize(0);
    pTText->SetFillStyle(0);
    pTText->SetTextFont(42);
   
    TCanvas *c = new TCanvas("c", "c", 50, 50, 600, 600);
    c->cd()->SetMargin(0.15, 0.05, 0.12, 0.06);
    //c->cd()->SetLeftMargin(0.15);
    //c->cd()->SetBottomMargin(0.12);
    USInvMass->GetYaxis()->SetTitle("Arb. Units");
    USInvMass->GetYaxis()->SetMaxDigits(2);
    USInvMass->GetYaxis()->SetRangeUser(0.0, 0.05E6);
    USInvMass->Draw("HIST E");
    LSInvMass->Draw("SAME E");
    //sbLine1->Draw("SAME");
    //sbLine2->Draw("SAME");
    //peak1->Draw("SAME");
    //peak2->Draw("SAME");
    //lineLeg->Draw();
    leg->Draw();
    text->Draw();
    pTText->Draw();


    TLegend *corrleg = new TLegend(0.46, 0.39, 0.88, 0.57);
    corrleg->AddEntry(corrected, "Corrected US Inv. Mass", "p");
    corrleg->AddEntry(fit, "Inv. Mass Fit", "l");
    corrleg->AddEntry(bgFit, "Residual BG Fit", "l");
    corrleg->SetLineWidth(0);

    TCanvas *c2 = new TCanvas("c2", "c2", 50, 50, 600, 600);
    c2->cd()->SetMargin(0.15, 0.05, 0.12, 0.06);
    //c2->cd()->SetLeftMargin(0.15);
    //c2->cd()->SetBottomMargin(0.12);
    corrected->GetYaxis()->SetTitle("Arb. Units");
    corrected->GetYaxis()->SetMaxDigits(2);
    corrected->GetYaxis()->SetRangeUser(-0.002E6, 0.04E6);
    corrected->Draw();
    bgFit->Draw("SAME");
    text->Draw();
    pTText->Draw();
    corrleg->Draw();

    TH1D* wideInvMass = (TH1D*)USInvMass->Clone("wideInvMass");
    wideInvMass->SetFillColor(36);
    TH1D* narrowInvMass=(TH1D*)USInvMass->Clone("narrowInvMass");
    narrowInvMass->SetFillColor(30);
    TH1D* LSBinvmass = (TH1D*)USInvMass->Clone("LSBinvmass");
    LSBinvmass->SetFillColor(kGray+2);
    TH1D* RSBinvmass = (TH1D*)USInvMass->Clone("RSBinvmass");
    RSBinvmass->SetFillColor(kGray+2);

    for(int i = 1; i<= wideInvMass->GetXaxis()->GetNbins(); i++){
        if(TMath::Abs(wideInvMass->GetXaxis()->GetBinCenter(i) - 1.020) > 0.010){
            wideInvMass->SetBinContent(i, 0.);
        }
        if(TMath::Abs(narrowInvMass->GetXaxis()->GetBinCenter(i) - 1.020) > 0.006){
            narrowInvMass->SetBinContent(i, 0.);
        }
        if(narrowInvMass->GetBinCenter(i) < 0.995 || narrowInvMass->GetBinCenter(i) > 1.005){
            LSBinvmass->SetBinContent(i, 0.);
        }
        if(narrowInvMass->GetBinCenter(i) < 1.040 || narrowInvMass->GetBinCenter(i) > 1.060){
            RSBinvmass->SetBinContent(i, 0.);
        }
    }

    TCanvas* c3 = new TCanvas("c3", "c3", 50, 50, 600, 600);
    c3->cd();
    USInvMass->Draw();
    wideInvMass->Draw("SAME");
    USInvMass->Draw("SAME");
    LSInvMass->Draw("SAME");
    leg->Draw();
    
    TCanvas* c4 = new TCanvas("c4", "c4", 50, 50, 600, 600);
    c4->cd();
    USInvMass->Draw();
    narrowInvMass->Draw("SAME");
    LSBinvmass->Draw("SAME");
    RSBinvmass->Draw("SAME");
    USInvMass->Draw("SAME");
    LSInvMass->Draw("SAME");
    leg->Draw();

    TLegend *leg2 = new TLegend(0.4581, 0.3927, 0.8809, 0.5637);
    leg2->AddEntry(USInvMass, "US Kaon Pairs", "l");
    
    TCanvas* c5 = new TCanvas("c5", "c5", 50, 50, 600, 600);
    c5->cd();
    USInvMass->Draw();
    narrowInvMass->Draw("SAME");
    LSBinvmass->Draw("SAME");
    RSBinvmass->Draw("SAME");
    USInvMass->Draw("SAME");
    leg2->Draw();

    TLegend *leg3 = new TLegend(0.4581, 0.3927, 0.8809, 0.5637);
    leg3->AddEntry(USInvMass, "US Kaon Pairs", "l");
    leg3->AddEntry(origLSInvMass, "LS Kaon Pairs", "l");

    TCanvas* c6 = new TCanvas("c6", "c6", 50, 50, 600, 600);
    c6->cd();
    USInvMass->Draw();
    narrowInvMass->Draw("SAME");
    LSBinvmass->Draw("SAME");
    RSBinvmass->Draw("SAME");
    USInvMass->Draw("SAME");
    origLSInvMass->Draw("SAME");
    leg3->Draw();


    //phi pT spectrum plotting
    //TFile* efffile = new TFile("~/alirepos/utaustin/efficiency/fits_17f2btriggerefficiency.root");
    TFile* efffile = new TFile("~/alirepos/utaustin/efficiency/fits_17f2b_newcodenewdenom.root");
    TF1* phiEff = (TF1*)efffile->Get("phiFit");
    TH1D* phiPT = new TH1D("phiPT", "phiPT", 8, 1.0, 5.0);
    TH1D* phiPTnoeff = new TH1D("phiPT", "phiPT", 8, 1.0, 5.0);

    //testing rapidity and pseudorapidity calculations...
    kkUSDist->GetAxis(3)->SetRangeUser(-0.8 + epsilon, 0.8 - epsilon);
    kkUSDist->GetAxis(0)->SetRangeUser(2.0 + epsilon, 4.0 - epsilon);
    kkUSDist->GetAxis(1)->SetRangeUser(1.014 + epsilon, 1.026 - epsilon);
    TH1D* kkrapidity = kkUSDist->Projection(4);
    kkrapidity->SetName("kkrapidity");

    kkUSDist->GetAxis(3)->SetRange(1, kkUSDist->GetAxis(3)->GetNbins());
    //kkUSDist->GetAxis(4)->SetRangeUser(-0.04 + epsilon, 0.46 - epsilon);
    TH1D* kkpseudorap = kkUSDist->Projection(3);
    kkpseudorap->SetName("kkpseudorap");
    trigDist->GetAxis(0)->SetRangeUser(2.0 + epsilon, 4.0 - epsilon);
    TH1D* trigEta = trigDist->Projection(2);
    trigEta->SetName("hybridEta");
    trigEta->SetMarkerColor(kRed);
    trigEta->SetLineColor(kRed);

    TCanvas* crap = new TCanvas("crap", "crap", 50, 50, 600, 600);
    crap->cd();
    kkrapidity->Draw();

    TCanvas* cpseudo  = new TCanvas("cpseudo", "cpseudo", 50, 50, 600, 600);
    cpseudo->cd();
    kkpseudorap->Draw();
    
    
    TCanvas* chybrideta  = new TCanvas("chybrideta", "chybrideta", 50, 50, 600, 600);
    chybrideta->cd();
    trigEta->Draw();

    
    TCanvas* cblank = new TCanvas("cblank", "cblank", 50, 50, 600, 600);
    cblank->cd();

    kkUSDist->GetAxis(4)->SetRange(0,0);
    kkUSDist->GetAxis(0)->SetRange(0,0);
    kkUSDist->GetAxis(1)->SetRange(0,0);
    


    //TCanvas* cmass[12];
    for(int i = 1; i <= 8; i++){
        kkUSDist->GetAxis(0)->SetRangeUser(1.0 + 0.5*(float(i-1.0)) +0.001, 1.0 + 0.5*(float(i)) - 0.001);
        kkLSDist->GetAxis(0)->SetRangeUser(1.0 + 0.5*(float(i-1.0)) +0.001, 1.0 + 0.5*(float(i)) - 0.001);
        //kkUSDist->GetAxis(3)->SetRangeUser(-0.8 + 0.0001, 0.8 - 0.0001);
        //kkLSDist->GetAxis(3)->SetRangeUser(-0.8 + 0.0001, 0.8 - 0.0001);
        kkUSDist->GetAxis(4)->SetRangeUser(0.0 +  0.0001, 0.6  - 0.0001);
        kkLSDist->GetAxis(4)->SetRangeUser(0.0 + 0.0001, 0.6 - 0.0001);
        TH1D* usinvmass = kkUSDist->Projection(1);
        usinvmass->SetName("usinvmass");
        TH1D* lsinvmass = kkLSDist->Projection(1);
        lsinvmass->SetName("lsinvmass");

        Float_t lsscale = usinvmass->Integral(usinvmass->GetXaxis()->FindBin(1.040+0.0001), usinvmass->GetXaxis()->FindBin(1.060-0.0001));
        lsscale = lsscale/lsinvmass->Integral(lsinvmass->GetXaxis()->FindBin(1.040+0.0001), lsinvmass->GetXaxis()->FindBin(1.060-0.0001));
        lsinvmass->Scale(lsscale);
        
        usinvmass->Add(lsinvmass, -1.0);
        for(int ibin = 1; ibin <= usinvmass->GetXaxis()->GetNbins(); ibin++){
            usinvmass->SetBinContent(ibin, usinvmass->GetBinContent(ibin)/usinvmass->GetXaxis()->GetBinWidth(ibin));
        }

        TF1* fitpt = new TF1("fitpt",  "[0]*TMath::Voigt(x - [1], [2], [3], 4) + pol2(4)",0.99, 1.07);
        fitpt->SetParameter(1, 1.020);
        fitpt->SetParameter(2, 0.0002);
        fitpt->SetParameter(0, 400);
        fitpt->FixParameter(3, 0.00426);
        fitpt->SetParLimits(1, 1.010, 1.030);
        fitpt->SetLineColor(kBlue);
        fitpt->SetLineStyle(7);
        fitpt->SetLineWidth(6);
        
        usinvmass->Fit(fitpt);


        //cmass[i-1] = new TCanvas(Form("cmass%d", i), Form("cmass%d", i), 50, 50, 600, 600);
        //cmass[i-1]->cd();
        //usinvmass->Draw();
        
        //Double_t numphis = usinvmass->Integral(usinvmass->GetXaxis()->FindBin(1.014+0.0001), usinvmass->GetXaxis()->FindBin(1.026-0.0001), "width"); //get numphis by histo integral
        Double_t error = 0.0;
        Double_t numphis = fitpt->Integral(1.014, 1.026);
        Double_t histint = usinvmass->IntegralAndError(usinvmass->GetXaxis()->FindBin(1.014 + 0.00001), usinvmass->GetXaxis()->FindBin(1.026 - 0.00001), error);
        
        numphis = numphis/0.802; //scale by signal pct.
        numphis = numphis/0.49; //scale by branching ratio
        numphis = numphis/phiEff->Eval(1.0 + 0.125 + 0.5*(float(i-1.0))); // scale by pt dependent efficiency
        numphis = numphis/0.5; //divide by pt bin width
        numphis = numphis/events->GetBinContent(3); //divide by number of events considered.
        phiPTnoeff->SetBinContent(i, numphis);
        numphis = numphis/0.6; //divide by rapidity range considered
        numphis = numphis/(2.0*TMath::Pi()); //scale for same factor 1/2pi(pT) as paper plot
        numphis = numphis/(1.0 + 0.25 + 0.5*(float(i-1.0)));
        phiPT->SetBinContent(i, numphis);
        phiPT->SetBinError(i, error/(0.802*0.49*0.5*0.6*2.0*TMath::Pi()*phiEff->Eval(1.0 + 0.125 + 0.5*(float(i-1.0)))*events->GetBinContent(3)*(1.0 + 0.25 + 0.5*(float(i-1.0)))));
                

    }

    //get published spectra
    TFile* pubfile = new TFile("~/phiStudies/pubRatioCheck/HEPData-ins1418181-v1-Table_11_mult_10_20.root");
    TDirectory* phidirmult = (TDirectory*) pubfile->GetDirectory("Table 11");
    TH1D* pubphihist = (TH1D*)(phidirmult->Get("Hist1D_y1")->Clone("pubphihist"));
    TH1D* pubphihistSyst = (TH1D*)(phidirmult->Get("Hist1D_y1_e2")->Clone("pubphihistSyst"));
    
    TFile* pubfile2 = new TFile("~/phiStudies/pubRatioCheck/HEPData-ins1418181-v1-Table_10_mult_5_10.root");
    TDirectory* phidirmult2 = (TDirectory*) pubfile2->GetDirectory("Table 10");
    TH1D* pubphihist2 = (TH1D*)(phidirmult2->Get("Hist1D_y1")->Clone("pubphihist2"));

    TFile* pubfile3 = new TFile("~/phiStudies/pubRatioCheck/HEPData-ins1418181-v1-Table_9_mult_0_5.root");
    TDirectory* phidirmult3 = (TDirectory*) pubfile3->GetDirectory("Table 9");
    TH1D* pubphihist3 = (TH1D*)(phidirmult3->Get("Hist1D_y1")->Clone("pubphihist3"));

    TFile* pubfile4 = new TFile("~/phiStudies/pubRatioCheck/HEPData-ins1418181-v1-Table_12_mult_20_40.root");
    TDirectory* phidirmult4 = (TDirectory*) pubfile4->GetDirectory("Table 12");
    TH1D* pubphihist4 = (TH1D*)(phidirmult4->Get("Hist1D_y1")->Clone("pubphihist4"));
    
    TFile* pubfile5 = new TFile("~/phiStudies/pubRatioCheck/HEPData-ins1418181-v1-Table_13_mult_40_60.root");
    TDirectory* phidirmult5 = (TDirectory*) pubfile5->GetDirectory("Table 13");
    TH1D* pubphihist5 = (TH1D*)(phidirmult5->Get("Hist1D_y1")->Clone("pubphihist5"));

    TFile* pubfile6 = new TFile("~/phiStudies/pubRatioCheck/HEPData-ins1418181-v1-Table_14_mult_60_80.root");
    TDirectory* phidirmult6 = (TDirectory*) pubfile6->GetDirectory("Table 14");
    TH1D* pubphihist6 = (TH1D*)(phidirmult6->Get("Hist1D_y1")->Clone("pubphihist6"));

    TFile* pubfile7 = new TFile("~/phiStudies/pubRatioCheck/HEPData-ins1418181-v1-Table_15_mult_80_100.root");
    TDirectory* phidirmult7 = (TDirectory*) pubfile7->GetDirectory("Table 15");
    TH1D* pubphihist7 = (TH1D*)(phidirmult7->Get("Hist1D_y1")->Clone("pubphihist7"));


    TH1D* pubphiMB = (TH1D*)pubphihist7->Clone("pubphiMB");
    pubphiMB->SetMarkerColor(kBlack);
    pubphiMB->SetMarkerStyle(23);
    TH1D* pubphi020 = (TH1D*)pubphihist7->Clone("pubphi020");

    //get relative syst errors for 10-20 mult bin
    pubphihistSyst->Divide(pubphihist);

    for(int i =1; i <= pubphiMB->GetXaxis()->GetNbins(); i++){
        Double_t binval = pubphihist3->GetBinContent(i)*0.05 + pubphihist2->GetBinContent(i)*0.05 + pubphihist->GetBinContent(i)*0.1 + pubphihist4->GetBinContent(i)*0.2 + pubphihist5->GetBinContent(i)*0.2 + pubphihist6->GetBinContent(i)*0.2 + pubphihist7->GetBinContent(i)*0.2;
        pubphiMB->SetBinContent(i, binval);
        Double_t binval020 = pubphihist3->GetBinContent(i)*0.25 + pubphihist2->GetBinContent(i)*0.25 + pubphihist->GetBinContent(i)*0.5;
        pubphi020->SetBinContent(i, binval020);
        //set systmatic errors as relative errors in 10-20 bin for the new 0-20 bin content
    }

    //adding together the 3 different mult bins to get 0-20%
    //pubphihist->Scale(0.5);
    //pubphihist2->Scale(0.25);
    //pubphihist3->Scale(0.25);
    
    //pubphihist->Add(pubphihist2);
    //pubphihist->Add(pubphihist3);

    for(int i = 1; i <= pubphihist->GetXaxis()->GetNbins(); i++){
        pubphihist->SetBinContent(i, pubphihist->GetBinContent(i)/(2.0*TMath::Pi()*pubphihist->GetXaxis()->GetBinCenter(i)));
        pubphihist2->SetBinContent(i, pubphihist2->GetBinContent(i)/(2.0*TMath::Pi()*pubphihist2->GetXaxis()->GetBinCenter(i)));
        pubphihist3->SetBinContent(i, pubphihist3->GetBinContent(i)/(2.0*TMath::Pi()*pubphihist3->GetXaxis()->GetBinCenter(i)));
        pubphihist4->SetBinContent(i, pubphihist4->GetBinContent(i)/(2.0*TMath::Pi()*pubphihist4->GetXaxis()->GetBinCenter(i)));
        pubphihist5->SetBinContent(i, pubphihist5->GetBinContent(i)/(2.0*TMath::Pi()*pubphihist5->GetXaxis()->GetBinCenter(i)));      
    }

    for(int i = 1; i <= pubphihist7->GetXaxis()->GetNbins(); i++){
        pubphihist6->SetBinContent(i, pubphihist6->GetBinContent(i)/(2.0*TMath::Pi()*pubphihist6->GetXaxis()->GetBinCenter(i)));
        pubphihist7->SetBinContent(i, pubphihist7->GetBinContent(i)/(2.0*TMath::Pi()*pubphihist7->GetXaxis()->GetBinCenter(i)));
        pubphiMB->SetBinContent(i, pubphiMB->GetBinContent(i)/(2.0*TMath::Pi()*pubphiMB->GetXaxis()->GetBinCenter(i)));
        pubphi020->SetBinContent(i, pubphi020->GetBinContent(i)/(2.0*TMath::Pi()*pubphi020->GetXaxis()->GetBinCenter(i)));
        pubphi020->SetBinError(i, pubphi020->GetBinContent(i)*pubphihistSyst->GetBinContent(i));

    }



    //scaling histograms to the paper plot scales
    pubphihist3->Scale(16.0);
    pubphihist2->Scale(8.0);
    pubphihist->Scale(4.0);
    pubphihist4->Scale(2.0);
    pubphihist5->Scale(0.5);
    pubphihist6->Scale(0.25);
    pubphihist7->Scale(0.125);

    //compute ratio of pub to mine
    TH1D* ptratio = new TH1D("ptratio", "ptratio", 6, 2.0, 5.0);     
    for(int i = 1; i <= 6; i++){
        ptratio->SetBinContent(i, phiPT->GetBinContent(phiPT->GetXaxis()->FindBin(2.0 + 0.25 + 0.5*float(i-1.0)))/pubphi020->GetBinContent(pubphihist->GetXaxis()->FindBin(2.0 + 0.25 + 0.5*float(i - 1.0))));
    }

    TCanvas* cpt = new TCanvas("cpt", "cpt", 50, 50, 600, 600);
    cpt->cd()->SetLogy();
    phiPT->SetMarkerStyle(23);
    phiPT->GetXaxis()->SetRangeUser(2.0, 4.5);
    phiPT->GetXaxis()->SetTitle("p_{T} (GeV/c)");
    phiPT->GetYaxis()->SetTitle("1/N_{evts} 1/2#pip_{T} d^{2}N/dp_{T}dy");
    phiPT->SetMarkerSize(2);
    phiPT->SetTitle("#phi p_{T} spectrum");
    //phiPT->Draw("P");
   /* pubphihist->SetMarkerStyle(22);
    pubphihist->SetMarkerColor(kYellow+1);
    pubphihist->GetYaxis()->SetRangeUser(2.0E-9, 8.0);
    pubphihist->Draw("P HIST");
    pubphihist2->SetMarkerStyle(22);
    pubphihist2->SetMarkerColor(kOrange);
    pubphihist2->Draw("P HIST SAME");
    pubphihist3->SetMarkerStyle(22);
    pubphihist3->SetMarkerColor(kRed);
    pubphihist3->Draw("P HIST SAME");
    pubphihist4->SetMarkerStyle(22);
    pubphihist4->SetMarkerColor(kGreen+1);
    pubphihist4->Draw("P HIST SAME");
    pubphihist5->SetMarkerStyle(22);
    pubphihist5->SetMarkerColor(kCyan+1);
    pubphihist5->Draw("P HIST SAME");
    pubphihist6->SetMarkerStyle(22);
    pubphihist6->SetMarkerColor(kViolet);
    pubphihist6->Draw("P HIST SAME");
    pubphihist7->SetMarkerStyle(22);
    pubphihist7->SetMarkerColor(kAzure);
    pubphihist7->Draw("P HIST SAME");
    */
    //pubphiMB->Draw("P HIST SAME");
   
    phiPT->Draw("P HIST E");
    pubphi020->SetMarkerStyle(22);
    pubphi020->SetMarkerSize(2);
    pubphi020->SetMarkerColor(kRed);
    pubphi020->SetLineColor(kRed);
    pubphi020->SetFillStyle(0);
    pubphi020->Draw("P E2 SAME");

    TLegend *ptleg = new TLegend(0.35, 0.70, 0.89, 0.88);
    ptleg->SetLineWidth(0);
    ptleg->AddEntry(phiPT, "My 0-20% #phi(1020) spectrum (rough)", "p");
    ptleg->AddEntry(pubphi020, "Published 0-20% #phi(1020) spectrum", "pl");
    ptleg->Draw();

    //pubphiMB->SetMarkerStyle(22);
    //pubphiMB->SetMarkerColor(kRed);
    //pubphiMB->Draw("P");
    
    //Double_t myYield = phiPT->Integral(phiPT->GetXaxis()->FindBin(2.0 + 0.0001), phiPT->GetXaxis()->FindBin(4.0 - 0.00001));
    //Double_t paperYield = pubphi020->Integral(pubphi020->GetXaxis()->FindBin(2.0 + 0.00001), pubphi020->GetXaxis()->FindBin(4.0 - 0.000001));
    Double_t myYield = 0.0;
    Double_t paperYield = 0.0;
    for(int i = 3; i<=6; i++){
        myYield += phiPT->GetBinContent(i);
        printf("myyield for %d: %f\n", i, phiPT->GetBinContent(i));
    }
    for(int i = 9; i<=12; i++){
        paperYield += pubphi020->GetBinContent(i);
        printf("paperyield for %d: %f\n", i, pubphi020->GetBinContent(i));
    }
    
    Double_t yieldratio = myYield/paperYield;

    printf("2 < pT < 4 yields\n");
    printf("my yield = %f\n", myYield);
    printf("paper yield = %f\n", paperYield);
    printf("mine/paper = %f\n", yieldratio);
    printf("my num yield = %f\n", phiPTnoeff->Integral(3, 6));
    
    TCanvas* cratio = new TCanvas("cratio", "cratio", 50, 50, 600, 600);
    cratio->cd();
    ptratio->Draw();

    //do trig dist comparison to measured p+K+pi spectra

    trigDist->GetAxis(0)->SetRange(0, 0);
    trigDist->GetAxis(2)->SetRangeUser(-0.8 + 0.0001, 0.8 - 0.00001);
    TH1D* mychargedpt = (TH1D*)trigDist->Projection(0);
    mychargedpt->Scale(1.0/events->GetBinContent(3));
    TH1D* myphipt = (TH1D*)phiPT->Clone("myphipt");
    for(int ibin = 1; ibin <= myphipt->GetXaxis()->GetNbins(); ibin++){
        myphipt->SetBinContent(ibin, myphipt->GetBinContent(ibin)*2.0*TMath::Pi()*myphipt->GetXaxis()->GetBinCenter(ibin)*myphipt->GetXaxis()->GetBinWidth(ibin));
    }

   /* 
    for(int ibin = 1; ibin < mychargept->GetXaxis()->GetNbins(); ibin++){
        mychargept->SetBinContent(ibin, mychargept->GetBinContent(ibin)/(2.0*TMath::Pi()*mychargept->GetXaxis()->GetBinWidth(ibin)8mychargept->GetXaxis()->GetBinCenter(ibin)));
    }
   */
    //Double_t myNumYield = myphipt->Integral(myphipt->GetXaxis()->FindBin(2.0 + 0.0001), myphipt->GetXaxis()->FindBin(4.0 - 0.00001));
    //Double_t mychargedyield = mychargedpt->Integral(mychargedpt->GetXaxis()->FindBin(2.1 + 0.00001), mychargedpt->GetXaxis()->FindBin(4.1  - 0.00001));
    Double_t myNumYield = 0.0;
    for(int i = 3; i<=6; i++){
        myNumYield += myphipt->GetBinContent(i)*0.5;
    }
    Double_t mychargedyield = 0.0;
    for(int i = 10; i<=19; i++){
        mychargedyield += mychargedpt->GetBinContent(i);
    }
    mychargedyield = mychargedyield/(1.6);

    TCanvas* cptcheck = new TCanvas("ptcheck", "ptcheck", 50, 50, 600, 600);
    cptcheck->cd();
    //myphipt->Draw("HIST SAME");
    //phiPTnoeff->Draw("HIST SAME");
    //mychargedpt->Draw("HIST SAME");
    myphipt->Divide(phiPTnoeff);
    myphipt->Draw("HIST");

    printf("my 0-20 phi yield: %f\n", float(myNumYield));
    printf("my 0-20 charge yield = %f\n", float(mychargedyield));
    printf("my 0-20 inclusive ratio: %f\n", float(myNumYield)/float(mychargedyield)); 


}
