void myPadSetUp(TPad *currentPad, float currentLeft, float currentTop, float currentRight, float currentBottom){
  currentPad->SetLeftMargin(currentLeft);
  currentPad->SetTopMargin(currentTop);
  currentPad->SetRightMargin(currentRight);
  currentPad->SetBottomMargin(currentBottom);
  return;
}

void SetStyle(Bool_t graypalette = kFALSE) {
  cout << "Setting style!" << endl;
  
  gStyle->Reset("Plain");
  gStyle->SetOptTitle(0);
  gStyle->SetOptStat(0);
  if(graypalette) gStyle->SetPalette(8,0);
  else gStyle->SetPalette(1);
  gStyle->SetCanvasColor(10);
  gStyle->SetCanvasBorderMode(0);
  gStyle->SetFrameLineWidth(1);
  gStyle->SetFrameFillColor(kWhite);
  gStyle->SetPadColor(10);
  gStyle->SetPadTickX(1);
  gStyle->SetPadTickY(1);
  gStyle->SetPadBottomMargin(0.15);
  gStyle->SetPadLeftMargin(0.15);
  gStyle->SetHistLineWidth(1);
  gStyle->SetHistLineColor(kRed);
  gStyle->SetFuncWidth(2);
  gStyle->SetFuncColor(kGreen);
  gStyle->SetLineWidth(2);
  gStyle->SetLabelSize(0.045,"xyz");
  gStyle->SetLabelOffset(0.01,"y");
  gStyle->SetLabelOffset(0.01,"x");
  gStyle->SetLabelColor(kBlack,"xyz");
  gStyle->SetTitleSize(0.05,"xyz");
  gStyle->SetTitleOffset(1.25,"y");
  gStyle->SetTitleOffset(1.2,"x");
  gStyle->SetTitleFillColor(kWhite);
  gStyle->SetTextSizePixels(26);
  gStyle->SetTextFont(42);
  //  gStyle->SetTickLength(0.04,"X");  gStyle->SetTickLength(0.04,"Y"); 

  gStyle->SetLegendBorderSize(0);
  gStyle->SetLegendFillColor(kWhite);
  //  gStyle->SetFillColor(kWhite);
  gStyle->SetLegendFont(42);


}

Double_t bin2val(Int_t bin){
    Double_t val = (-0.5*TMath::Pi())+((Double_t)bin)*(2.0*TMath::Pi()/16.0);
    return val;
}

Double_t flineStd(Double_t *x, Double_t *par){
    if((x[0] > bin2val(1) && x[0] < bin2val(7)) || (x[0] > bin2val(9) && x[0] < bin2val(15))){
        TF1::RejectPoint();
        return 0;
    }
    return par[0];
}

Double_t fline6bin(Double_t *x, Double_t *par){
    if((x[0] > bin2val(2) && x[0] < bin2val(6)) || (x[0] > bin2val(9) && x[0] < bin2val(15))){
        TF1::RejectPoint();
        return 0;
    }
    return par[0];
}

Double_t flineNoLast(Double_t *x, Double_t *par){
    if((x[0] > bin2val(1) && x[0] < bin2val(7)) || x[0] > bin2val(9)){
        TF1::RejectPoint();
        return 0;
    }
    return par[0];
}



TH1D* getHisto(TString filename, TString histotype, TString histoname, TString mult, Float_t etamin, Float_t etamax, Int_t color, Int_t markerstyle){
    TFile *hhFile = new TFile(filename.Data());
    TH2D* histo2D = (TH2D*)hhFile->Get(histoname.Data());
    TString newhistoname = histoname +"_"+ mult;
    histo2D->SetName(newhistoname);
    histo2D->Sumw2();
    TString histoname1D = histotype + "dphi_" + mult;
    TH1D* histo1D = (TH1D*)histo2D->ProjectionY(histoname1D.Data(), histo2D->GetXaxis()->FindBin(etamin+ 0.001), histo2D->GetXaxis()->FindBin(etamax - 0.001));
    if(histotype == "hh" && mult == "50_100"){
        //histo1D->Rebin();
    }
    histo1D->SetLineWidth(2);
    histo1D->SetLineColor(color);
    histo1D->SetMarkerColor(color);
    histo1D->SetMarkerStyle(markerstyle);
    histo1D->SetMarkerSize(2);
    histo1D->GetXaxis()->SetTitle("#Delta#varphi");
    histo1D->SetTitle("");
    histo1D->GetYaxis()->SetTitleOffset(1.60);
    histo1D->GetYaxis()->SetMaxDigits(2);
    histo1D->GetXaxis()->SetTitleSize(0.05);
    histo1D->GetXaxis()->SetTitleOffset(0.90);
    return histo1D;
}

TF1* setupFit(TString fitname, TH1D* hist, Int_t color, Int_t linestyle, Int_t bgMethod){
    //Do straight line fits for systematics
    TF1 *linefit = new TF1(Form("line%s", fitname.Data()), flineStd, -0.5*TMath::Pi(), 1.5*TMath::Pi(), 1); //use std. 4 points
    //TF1 *linefit = new TF1(Form("line%s", fitname.Data()), fline6bin, -0.5*TMath::Pi(), 1.5*TMath::Pi(), 1); //use 6 points
    //TF1 *linefit = new TF1(Form("line%s", fitname.Data()), flineNoLast, -0.5*TMath::Pi(), 1.5*TMath::Pi(), 1); //use 3 points (no last bin)

    hist->Fit(linefit);

    TF1 *basefit = new TF1(fitname, "gaus(0) + gaus(3) + ([3]/([5]))*exp(-(((x - [4] + 2.0*TMath::Pi())^2)/(2*[5]^2)))+ pol0(6)", -TMath::Pi()/2.0, 3.0*TMath::Pi()/2.0);
    
    switch(bgMethod){
        case 0: //4 bins for straight line
            basefit->FixParameter(6, 1.0/4.0*(hist->GetBinContent(8)+hist->GetBinContent(9)+hist->GetBinContent(16)+hist->GetBinContent(1)));
            break;
        case 1: //6 bins for straight line
            basefit->FixParameter(6, (1.0/6.0)*(hist->GetBinContent(8)+hist->GetBinContent(9)+hist->GetBinContent(16)+hist->GetBinContent(1)+hist->GetBinContent(2)+hist->GetBinContent(7)));
            break;
        case 2: //only around near peak (4)
            basefit->FixParameter(6, (1.0/4.0)*(hist->GetBinContent(8)+hist->GetBinContent(1)+hist->GetBinContent(2)+hist->GetBinContent(7)));
            break;
        case 3: //not last bin (3)
            basefit->FixParameter(6, (1.0/3.0)*(hist->GetBinContent(8)+hist->GetBinContent(9)+hist->GetBinContent(1)));
            break;
        case 4: //free fit
            basefit->SetParameter(0, hist->GetBinContent(hist->GetXaxis()->FindBin(0.0)) - basefit->GetParameter(6));
            basefit->SetParLimits(0, basefit->GetParameter(0)*0.1, basefit->GetParameter(0)*3.0);
            break;
        default: //fix straight line parameter to the line fit from above
            basefit->FixParameter(6, linefit->GetParameter(0));
            break;
    }
       
    /*if(fitname == "corrFit"){
        basefit->FixParameter(0, 1.68920E-01 - 3.31591E-04); 
    }else if(fitname == "corrFit2050"){
        basefit->FixParameter(0, 1.69976E-01 - 3.03244E-04);
    }else if(fitname == "corrFit50100"){
        basefit->FixParameter(0, 1.73948E-01 - 3.76013E-04);
    }*/
    basefit->FixParameter(1, 0.0);
    //basefit->SetParameter(1, 0.0);
    basefit->SetParLimits(1, -0.5, 0.5);
    basefit->SetParameter(2, 0.5);
    basefit->SetParLimits(2, 0.1, 1.0);
    basefit->SetParameter(3, hist->GetBinContent(hist->GetXaxis()->FindBin(3.14)) - basefit->GetParameter(6));
    basefit->SetParLimits(3, basefit->GetParameter(3)*0.1, basefit->GetParameter(3)*3.0);
    basefit->FixParameter(4, 3.14159);
    //basefit->SetParameter(4, 3.14);
    basefit->SetParLimits(4, 3.0, 3.25);
    basefit->SetParameter(5, 0.5);
    basefit->SetParLimits(5, 0.1, 1.0);

    basefit->SetParNames("peak1", "mean1", "sigma1", "peak2", "mean2", "sigma2", "BG");

    basefit->SetLineColor(color);
    basefit->SetLineWidth(6);
    basefit->SetLineStyle(linestyle);
    return basefit;
}

void intUSRatioPlot2mult(TString outputstring, TString input_0_20 = "", TString input_20_50 = "", TString input_50_80 = "", TString scaleSuffix = "avgscale", Int_t bgMethod=0){
    gStyle->SetOptStat(0);
    gStyle->SetOptFit(0);
    gStyle->SetErrorX(0);

    if(input_0_20.EqualTo("")) input_0_20 = "~/phiStudies/results_onlineEff/Combined/US_syst_trig_4_8_assoc_1.5_2.0_effcorr_hPhi_0_50_combined.root";
    if(input_20_50.EqualTo("")) input_20_50 = "~/phiStudies/results_onlineEff/Combined/US_syst_trig_4_8_assoc_1.5_2.0_effcorr_hPhi_50_80_combined.root";
    if(input_50_80.EqualTo("")) input_50_80 = "~/phiStudies/results_onlineEff/Combined/US_syst_trig_4_8_assoc_1.5_2.0_effcorr_hPhi_50_80_combined.root";
    
    TH1D* hhdphi_0_20 = getHisto("~/phiStudies/results_onlineEff/Combined/trig_4_8_assoc_1.5_2.0_effcorr_hh_0_50.root", "hh", "hh2D", "Eff_0_20", -1.2, 1.2, kBlue+2, 21);
    TH1D* hPhidphi_0_20 = getHisto(input_0_20, "hPhi", Form("AvgUSsubhPhi2Dpeak%s", scaleSuffix.Data()), "Eff_0_20", -1.2, 1.2, kRed+2, 22);


    //scale for inv. mass range
    //hPhidphi_0_20->Scale(1.0/0.897); //wide mass
    //hPhidphi_0_20->Scale(1.0/0.803); //narrow mass

    TF1 *corrFit = setupFit("corrFit", hhdphi_0_20, kBlue, 7, bgMethod);
    TF1 *corrFit2 = setupFit("corrFit2", hPhidphi_0_20, kRed, 7, bgMethod);

    hhdphi_0_20->Fit("corrFit", "R0");
    hPhidphi_0_20->Fit("corrFit2", "R0");

    TF1 *hhBG = new TF1("hhBG", "pol0(0)", -TMath::Pi()/2.0, 3.0*TMath::Pi()/2.0);
    hhBG->SetParLimits(0, 0.00001, 10000000.0);
    hhBG->SetParameter(0, 1.0*corrFit->GetParameter(6));
    hhBG->SetLineStyle(2);

    TF1 *hphiBG = new TF1("hphiBG", "pol0(0)", -TMath::Pi()/2.0, 3.0*TMath::Pi()/2.0);
    hphiBG->SetParLimits(0, 0.00001, 10000000.0);
    hphiBG->SetParameter(0, 1.0*corrFit2->GetParameter(6));
    hphiBG->SetLineStyle(2);

    Double_t near0_20hPhiError = 0;
    Double_t near0_20hhError = 0;
    Double_t away0_20hPhiError = 0;
    Double_t away0_20hhError = 0;
    //Double_t mid0_20hPhiError =(1.0/3.0)*TMath::Sqrt(TMath::Power(hPhidphi_0_20->GetBinError(8),2) + TMath::Power(hPhidphi_0_20->GetBinError(16),2) + TMath::Power(hPhidphi_0_20->GetBinError(1),2));
    Double_t mid0_20hPhiError = corrFit2->GetParError(6);
    //Double_t mid0_20hhError = (1.0/3.0)*TMath::Sqrt(TMath::Power(hhdphi_0_20->GetBinError(8),2) + TMath::Power(hhdphi_0_20->GetBinError(16),2) + TMath::Power(hhdphi_0_20->GetBinError(1),2));
    Double_t mid0_20hhError = corrFit->GetParError(6);
    Double_t total0_20hPhiError;
    Double_t total0_20hhError = 0;   

    Double_t near0_20hPhiYield = hPhidphi_0_20->IntegralAndError(1,8,near0_20hPhiError, "width") - hphiBG->Integral(hPhidphi_0_20->GetXaxis()->GetBinLowEdge(1), hPhidphi_0_20->GetXaxis()->GetBinUpEdge(8));
    //near0_20hPhiError = TMath::Sqrt(TMath::Power(near0_20hPhiError, 2) + TMath::Power(6.0*mid0_20hPhiError, 2));
    Double_t near0_20hhYield = hhdphi_0_20->IntegralAndError(1,8,near0_20hhError, "width") - hhBG->Integral(hPhidphi_0_20->GetXaxis()->GetBinLowEdge(1), hPhidphi_0_20->GetXaxis()->GetBinUpEdge(8));
    //near0_20hhError = TMath::Sqrt(TMath::Power(near0_20hhError, 2) + TMath::Power(6.0*mid0_20hhError, 2));
    Double_t away0_20hPhiYield = hPhidphi_0_20->IntegralAndError(9,16,away0_20hPhiError, "width") - hphiBG->Integral(hPhidphi_0_20->GetXaxis()->GetBinLowEdge(9), hPhidphi_0_20->GetXaxis()->GetBinUpEdge(16));
    //away0_20hPhiError = TMath::Sqrt(TMath::Power(away0_20hPhiError, 2) + TMath::Power(8.0*mid0_20hPhiError, 2));
    Double_t away0_20hhYield = hhdphi_0_20->IntegralAndError(9,16,away0_20hhError, "width")- hhBG->Integral(hPhidphi_0_20->GetXaxis()->GetBinLowEdge(9), hPhidphi_0_20->GetXaxis()->GetBinUpEdge(16));
    //away0_20hhError = TMath::Sqrt(TMath::Power(away0_20hhError, 2) + TMath::Power(8.0*mid0_20hhError, 2));
    Double_t total0_20hPhiYield = hPhidphi_0_20->IntegralAndError(1,16,total0_20hPhiError, "width");
    Double_t total0_20hhYield = hhdphi_0_20->IntegralAndError(1,16,total0_20hhError, "width");
    Double_t mid0_20hPhiYield = hphiBG->Integral(hPhidphi_0_20->GetXaxis()->GetBinLowEdge(1), hPhidphi_0_20->GetXaxis()->GetBinUpEdge(16));
    mid0_20hPhiError = mid0_20hPhiError*16.0;
    Double_t mid0_20hhYield = hhBG->Integral(hPhidphi_0_20->GetXaxis()->GetBinLowEdge(1), hPhidphi_0_20->GetXaxis()->GetBinUpEdge(16));
    mid0_20hhError = mid0_20hhError*16.0; 
    Double_t jet0_20hPhi = hPhidphi_0_20->Integral(1, 16) - mid0_20hPhiYield;
    Double_t jet0_20hh = hhdphi_0_20->Integral(1, 16) - mid0_20hhYield;

    Double_t near020 = near0_20hPhiYield/near0_20hhYield;
    Double_t near020Er = near020*TMath::Sqrt(TMath::Power(near0_20hPhiError/near0_20hPhiYield, 2) + TMath::Power(near0_20hhError/near0_20hhYield, 2));
    Double_t away020 = away0_20hPhiYield/away0_20hhYield;
    Double_t away020Er = away020*TMath::Sqrt(TMath::Power(away0_20hPhiError/away0_20hPhiYield, 2) + TMath::Power(away0_20hhError/away0_20hhYield, 2));
    Double_t mid020 = mid0_20hPhiYield/mid0_20hhYield;
    Double_t mid020Er = mid020*TMath::Sqrt(TMath::Power(mid0_20hPhiError/mid0_20hPhiYield, 2) + TMath::Power(mid0_20hhError/mid0_20hhYield, 2));
    Double_t total020 = total0_20hPhiYield/total0_20hhYield;
    Double_t total020Er = total020*TMath::Sqrt(TMath::Power(total0_20hPhiError/total0_20hPhiYield, 2) + TMath::Power(total0_20hhError/total0_20hhYield, 2));
    
    TH1D *ratios020 = new TH1D("ratios020", "(h-#phi / h-h) Ratios", 4, 0, 4);
    ratios020->GetXaxis()->SetBinLabel(1, "near-side");
    ratios020->SetBinContent(1, near020);
    ratios020->SetBinError(1, near020Er);
    ratios020->GetXaxis()->SetBinLabel(2, "mid");
    ratios020->SetBinContent(2, mid020);
    ratios020->SetBinError(2, mid020Er);
    ratios020->GetXaxis()->SetBinLabel(3, "away-side");
    ratios020->SetBinContent(3, away020);
    ratios020->SetBinError(3, away020Er);
    ratios020->GetXaxis()->SetBinLabel(4, "total");
    ratios020->SetBinContent(4, total020);
    ratios020->SetBinError(4, total020Er);
    ratios020->SetMarkerStyle(22);
    ratios020->SetMarkerColor(kRed+2);
    ratios020->SetLineColor(kRed+2);
    ratios020->SetMarkerSize(2);
    ratios020->SetLineWidth(2);

    TF1* paperfit020 = new TF2("paperfit020", "[0]*(1+2*0.12*0.1*cos(2.0*x))", -0.5*TMath::Pi(), 1.5*TMath::Pi());
    paperfit020->SetParameter(0, hphiBG->GetParameter(0));


   //20-50 section 
    TH1D* hhdphi_20_50 = getHisto("~/phiStudies/results_onlineEff/Combined/trig_4_8_assoc_1.5_2.0_effcorr_hh_50_80.root", "hh", "hh2D", "Eff_20_50", -1.2, 1.2, kBlue+2, 21);
    TH1D* hPhidphi_20_50 = getHisto(input_20_50, "hPhi", Form("AvgUSsubhPhi2Dpeak%s", scaleSuffix.Data()), "Eff_20_50", -1.2, 1.2, kRed+2, 22);


    //scale for inv. mass range
    //hPhidphi_20_50->Scale(1.0/0.897); //wide mass
    //hPhidphi_20_50->Scale(1.0/0.803); //narrow mass

    TF1 *corrFit2050 = setupFit("corrFit2050", hhdphi_20_50, kBlue, 7, bgMethod); 

    TF1 *corrFit2_2050 = setupFit("corrFit2_2050", hPhidphi_20_50, kRed, 7, bgMethod);

    hhdphi_20_50->Fit("corrFit2050", "R0");
    hPhidphi_20_50->Fit("corrFit2_2050", "R0");

    TF1 *hhBG_20_50 = new TF1("hhBG_20_50", "pol0(0)", -TMath::Pi()/2.0, 3.0*TMath::Pi()/2.0);
    hhBG_20_50->SetParLimits(0, 0.00001, 10000000.0);
    hhBG_20_50->SetParameter(0, 1.0*corrFit2050->GetParameter(6));
    hhBG_20_50->SetLineStyle(2);

    TF1 *hphiBG_20_50 = new TF1("hphiBG_20_50", "pol0(0)", -TMath::Pi()/2.0, 3.0*TMath::Pi()/2.0);
    hphiBG_20_50->SetParLimits(0, 0.00001, 10000000.0);
    hphiBG_20_50->SetParameter(0, 1.0*corrFit2_2050->GetParameter(6));
    hphiBG_20_50->SetLineStyle(2);

    Double_t near20_50hPhiError = 0;
    Double_t near20_50hhError = 0;
    Double_t away20_50hPhiError = 0;
    Double_t away20_50hhError = 0;
    //Double_t mid20_50hPhiError =(1.0/3.0)*TMath::Sqrt(TMath::Power(hPhidphi_20_50->GetBinError(8),2) + TMath::Power(hPhidphi_20_50->GetBinError(16),2) + TMath::Power(hPhidphi_20_50->GetBinError(1),2));
    Double_t mid20_50hPhiError = corrFit2_2050->GetParError(6);
    //Double_t mid20_50hhError = (1.0/3.0)*TMath::Sqrt(TMath::Power(hhdphi_20_50->GetBinError(8),2) + TMath::Power(hhdphi_20_50->GetBinError(16),2) + TMath::Power(hhdphi_20_50->GetBinError(1),2));
    Double_t mid20_50hhError = corrFit2050->GetParError(6);
    Double_t total20_50hPhiError = 0;
    Double_t total20_50hhError = 0;

    Double_t near20_50hPhiYield = hPhidphi_20_50->IntegralAndError(1,8,near20_50hPhiError, "width") - hphiBG_20_50->Integral(hPhidphi_20_50->GetXaxis()->GetBinLowEdge(1), hPhidphi_20_50->GetXaxis()->GetBinUpEdge(8));
    //near20_50hPhiError = TMath::Sqrt(TMath::Power(near20_50hPhiError, 2) + TMath::Power(6.0*mid20_50hPhiError, 2));
    Double_t near20_50hhYield = hhdphi_20_50->IntegralAndError(1,8,near20_50hhError, "width") - hhBG_20_50->Integral(hPhidphi_20_50->GetXaxis()->GetBinLowEdge(1), hPhidphi_20_50->GetXaxis()->GetBinUpEdge(8));
    //near20_50hhError = TMath::Sqrt(TMath::Power(near20_50hhError, 2) + TMath::Power(6.0*mid20_50hhError, 2));
    Double_t away20_50hPhiYield = hPhidphi_20_50->IntegralAndError(9,16,away20_50hPhiError, "width") - hphiBG_20_50->Integral(hPhidphi_20_50->GetXaxis()->GetBinLowEdge(9), hPhidphi_20_50->GetXaxis()->GetBinUpEdge(16));
    //away20_50hPhiError = TMath::Sqrt(TMath::Power(away20_50hPhiError, 2) + TMath::Power(8.0*mid20_50hPhiError, 2));
    Double_t away20_50hhYield = hhdphi_20_50->IntegralAndError(9,16,away20_50hhError, "width")- hhBG_20_50->Integral(hPhidphi_20_50->GetXaxis()->GetBinLowEdge(9), hPhidphi_20_50->GetXaxis()->GetBinUpEdge(16));
    //away20_50hhError = TMath::Sqrt(TMath::Power(away20_50hhError, 2) + TMath::Power(8.0*mid20_50hhError, 2));
    Double_t mid20_50hPhiYield = hphiBG_20_50->Integral(hPhidphi_20_50->GetXaxis()->GetBinLowEdge(1), hPhidphi_20_50->GetXaxis()->GetBinUpEdge(16));
    mid20_50hPhiError = mid20_50hPhiError*16.0;
    Double_t mid20_50hhYield = hhBG_20_50->Integral(hPhidphi_20_50->GetXaxis()->GetBinLowEdge(1), hPhidphi_20_50->GetXaxis()->GetBinUpEdge(16));
    mid20_50hhError = mid20_50hhError*16.0;
    Double_t total20_50hPhiYield = hPhidphi_20_50->IntegralAndError(1, 16,total20_50hPhiError, "width");
    Double_t total20_50hhYield = hhdphi_20_50->IntegralAndError(1, 16,total20_50hhError, "width");
    Double_t jet20_50hPhi = hPhidphi_20_50->Integral(1, 16) - mid20_50hPhiYield;
    Double_t jet20_50hh = hhdphi_20_50->Integral(1, 16) - mid20_50hhYield;


    Double_t near2050 = near20_50hPhiYield/near20_50hhYield;
    Double_t near2050Er = near2050*TMath::Sqrt(TMath::Power(near20_50hPhiError/near20_50hPhiYield, 2) + TMath::Power(near20_50hhError/near20_50hhYield, 2));
    Double_t away2050 = away20_50hPhiYield/away20_50hhYield;
    Double_t away2050Er = away2050*TMath::Sqrt(TMath::Power(away20_50hPhiError/away20_50hPhiYield, 2) + TMath::Power(away20_50hhError/away20_50hhYield, 2));
    Double_t mid2050 = mid20_50hPhiYield/mid20_50hhYield;
    Double_t mid2050Er = mid2050*TMath::Sqrt(TMath::Power(mid20_50hPhiError/mid20_50hPhiYield, 2) + TMath::Power(mid20_50hhError/mid20_50hhYield, 2));
    Double_t total2050 = total20_50hPhiYield/total20_50hhYield;
    Double_t total2050Er = total2050*TMath::Sqrt(TMath::Power(total20_50hPhiError/total20_50hPhiYield, 2) + TMath::Power(total20_50hhError/total20_50hhYield, 2));

    TH1D *ratios2050 = new TH1D("ratios2050", "(h-#phi / h-h) Ratios", 4, 0, 4);
    ratios2050->GetXaxis()->SetBinLabel(1, "near-side");
    ratios2050->SetBinContent(1, near2050);
    ratios2050->SetBinError(1, near2050Er);
    ratios2050->GetXaxis()->SetBinLabel(2, "mid");
    ratios2050->SetBinContent(2, mid2050);
    ratios2050->SetBinError(2, mid2050Er);
    ratios2050->GetXaxis()->SetBinLabel(3, "away-side");
    ratios2050->SetBinContent(3, away2050);
    ratios2050->SetBinError(3, away2050Er);
    ratios2050->GetXaxis()->SetBinLabel(4, "total");
    ratios2050->SetBinContent(4, total2050);
    ratios2050->SetBinError(4, total2050Er);
    ratios2050->SetMarkerStyle(21);
    ratios2050->SetMarkerColor(kBlue+2);
    ratios2050->SetLineColor(kBlue+2);
    ratios2050->SetMarkerSize(2);
    ratios2050->SetLineWidth(2);


    TF1* paperfit2050 = new TF2("paperfit2050", "[0]*(1+2*0.12*0.1*cos(2.0*x))", -0.5*TMath::Pi(), 1.5*TMath::Pi());
    paperfit2050->SetParameter(0, hphiBG_20_50->GetParameter(0));

    TLegend  *ratioslegend = new TLegend(0.3791, 0.1518, 0.8772, 0.2688);
    ratioslegend->SetMargin(0.35);
    ratioslegend->AddEntry(ratios020, "0-20%", "p");
    ratioslegend->AddEntry(ratios2050, "20-50%", "p");
    
    TLine *line = new TLine(3.0, 0.0, 3.0, 0.0040);
    line->SetLineStyle(7);
    line->SetLineWidth(2);
  
    TCanvas *testc = new TCanvas("test", "test",50, 50, 600, 600);
    testc->cd();
    ratios2050->GetYaxis()->SetRangeUser(0.0000, 0.0040);
    ratios2050->GetXaxis()->SetLabelSize(0.07);
    ratios2050->Draw("P SAME");
    ratios020->GetXaxis()->SetLimits(-0.05, 3.95);
    ratios020->Draw("P SAME");
    line->Draw("SAME");
    ratioslegend->Draw("SAME");

    //initialize all sytematic error values
    Double_t near0_20hPhiSystError = 0.068;
    Double_t near20_50hPhiSystError = 0.077;
    Double_t away0_20hPhiSystError = 0.097;
    Double_t away20_50hPhiSystError = 0.061;
    Double_t mid0_20hPhiSystError = 0.017;
    Double_t mid20_50hPhiSystError = 0.013;
    Double_t total0_20hPhiSystError = 0.014;
    Double_t total20_50hPhiSystError = 0.012;

    Double_t near0_20hhSystError = 0;
    Double_t near20_50hhSystError = 0;
    Double_t away0_20hhSystError = 0;
    Double_t away20_50hhSystError = 0;
    Double_t mid0_20hhSystError = 0;
    Double_t mid20_50hhSystError = 0;
    Double_t total0_20hhSystError = 0;
    Double_t total20_50hhSystError = 0;

    //total systematic errors from all sources for yields (plotSystErrors.cxx)
    Double_t nearsysthphi[3] = {0.082, 0.09, 0.14};
    Double_t awaysysthphi[3] = {0.108, 0.076, 0.12};
    Double_t uesysthphi[3] = {0.048, 0.047, 0.059};
    Double_t totalsysthphi[3] = {0.047, 0.047, 0.056};
    Double_t hhsyst[3] = {0.04, 0.04, 0.04};


    //Setup single yield arrays for different regions
    Double_t nearhPhiYieldArray[2] = {near20_50hPhiYield*300.0, near0_20hPhiYield*300.0};
    Double_t nearhPhiYieldArrayErr[2] = {near20_50hPhiError*300.0, near0_20hPhiError*300.0};
    Double_t nearhPhiYieldArraySystErr[2] = {near20_50hPhiYield*nearsysthphi[1]*300.0, near0_20hPhiYield*nearsysthphi[0]*300.0};   
    Double_t awayhPhiYieldArray[2] = {away20_50hPhiYield*300.0, away0_20hPhiYield*300.0};
    Double_t awayhPhiYieldArrayErr[2] = {away20_50hPhiError*300.0, away0_20hPhiError*300.0};
    Double_t awayhPhiYieldArraySystErr[2] = {away20_50hPhiYield*awaysysthphi[1]*300.0, away0_20hPhiYield*awaysysthphi[0]*300.0};
    Double_t bulkhPhiYieldArray[2] = {mid20_50hPhiYield*100.0, mid0_20hPhiYield*100.0};
    Double_t bulkhPhiYieldArrayErr[2] = {mid20_50hPhiError*100.0, mid0_20hPhiError*100.0};
    Double_t bulkhPhiYieldArraySystErr[2] = {mid20_50hPhiYield*uesysthphi[1]*100.0, mid0_20hPhiYield*uesysthphi[0]*100.0};
    Double_t totalhPhiYieldArray[2] = {total20_50hPhiYield*100.0, total0_20hPhiYield*100.0};
    Double_t totalhPhiYieldArrayErr[2] = {total20_50hPhiError*100.0, total0_20hPhiError*100.0};
    Double_t totalhPhiYieldArraySystErr[2] = {total20_50hPhiYield*totalsysthphi[1]*100.0, total0_20hPhiYield*totalsysthphi[0]*100.0};

    Double_t nearhhYieldArray[2] = {near20_50hhYield, near0_20hhYield};
    Double_t nearhhYieldArrayErr[2] = {near20_50hhError, near0_20hhError};
    Double_t nearhhYieldArraySystErr[2] = {near20_50hhYield*hhsyst[1], near0_20hhYield*hhsyst[0]};
    Double_t awayhhYieldArray[2] = {away20_50hhYield, away0_20hhYield};
    Double_t awayhhYieldArrayErr[2] = {away20_50hhError, away0_20hhError};
    Double_t awayhhYieldArraySystErr[2] = {away20_50hhYield*hhsyst[1], away0_20hhYield*hhsyst[0]};
    Double_t bulkhhYieldArray[2] = {mid20_50hhYield, mid0_20hhYield};
    Double_t bulkhhYieldArrayErr[2] = {mid20_50hhError, mid0_20hhError};
    Double_t bulkhhYieldArraySystErr[2] = {mid20_50hhYield*hhsyst[1], mid0_20hhYield*hhsyst[0]};
    Double_t totalhhYieldArray[2] = {total20_50hhYield, total0_20hhYield};
    Double_t totalhhYieldArrayErr[2] = {total20_50hhError, total0_20hhError};
    Double_t totalhhYieldArraySystErr[2] = {total20_50hhYield*hhsyst[1], total0_20hhYield*hhsyst[0]};


    //Plot ratio as a function of multiplicity for the different angular regions
    Double_t nearArray[2] = {ratios2050->GetBinContent(1), ratios020->GetBinContent(1)};
    Double_t nearArrayErr[2] = {ratios2050->GetBinError(1), ratios020->GetBinError(1)};
    Double_t awayArray[2] = {ratios2050->GetBinContent(3), ratios020->GetBinContent(3)};
    Double_t awayArrayErr[2] = {ratios2050->GetBinError(3), ratios020->GetBinError(3)};
    Double_t bulkArray[2] = {ratios2050->GetBinContent(2), ratios020->GetBinContent(2)};
    Double_t bulkArrayErr[2] = {ratios2050->GetBinError(2), ratios020->GetBinError(2)};
    Double_t totalArray[2] = {ratios2050->GetBinContent(4), ratios020->GetBinContent(4)};
    Double_t totalArrayErr[2] = {ratios2050->GetBinError(4), ratios020->GetBinError(4)};

    //ratio of Near Jet to Underlying Event vs. Multiplicity
    
    Double_t jet2UEhPhi[2] = {mid20_50hPhiYield/total20_50hPhiYield, mid0_20hPhiYield/total0_20hPhiYield};
    Double_t jet2UEhh[2] = {mid20_50hhYield/total20_50hhYield, mid0_20hhYield/total0_20hhYield}; 
   
    Double_t near2tothPhi[2] = {near20_50hPhiYield/total20_50hPhiYield, near0_20hPhiYield/total0_20hPhiYield};
    Double_t near2tothh[2] = {near20_50hhYield/total20_50hhYield, near0_20hhYield/total0_20hhYield}; 
    
    Double_t away2tothPhi[2] = {away20_50hPhiYield/total20_50hPhiYield, away0_20hPhiYield/total0_20hPhiYield};
    Double_t away2tothh[2] = {away20_50hhYield/total20_50hhYield, away0_20hhYield/total0_20hhYield}; 

    Double_t UE2tothPhi[2] = {mid20_50hPhiYield/total20_50hPhiYield, mid0_20hPhiYield/total0_20hPhiYield};
    Double_t UE2tothh[2] = {mid20_50hhYield/total20_50hhYield, mid0_20hhYield/total0_20hhYield}; 
   
    //systematic errors from the changing the fitting parameters
  /*  Double_t nearArraySystErr[3] = {ratios50100->GetBinContent(1)*0.068, ratios2050->GetBinContent(1)*0.077, ratios020->GetBinContent(1)*0.133};
    Double_t awayArraySystErr[3] = {ratios50100->GetBinContent(3)*0.097, ratios2050->GetBinContent(3)*0.061, ratios020->GetBinContent(3)*0.111};
    Double_t bulkArraySystErr[3] = {ratios50100->GetBinContent(2)*0.017, ratios2050->GetBinContent(2)*0.013, ratios020->GetBinContent(2)*0.037};
    Double_t totalArraySystErr[3] = {ratios50100->GetBinContent(4)*0.014, ratios2050->GetBinContent(4)*0.012, ratios020->GetBinContent(4)*0.032};
  */
    //Systematic errors for the ratio plot:
    Double_t nearsysttot[2];
    Double_t awaysysttot[2];
    Double_t uesysttot[2];
    Double_t totalsysttot[2];
    for(int i = 0; i < 2; i++){
        nearsysttot[i] = TMath::Sqrt(TMath::Power(nearsysthphi[i], 2.0) + TMath::Power(hhsyst[i], 2.0));
        awaysysttot[i] = TMath::Sqrt(TMath::Power(awaysysthphi[i], 2.0) + TMath::Power(hhsyst[i], 2.0));
        uesysttot[i] = TMath::Sqrt(TMath::Power(uesysthphi[i], 2.0) + TMath::Power(hhsyst[i], 2.0));
        totalsysttot[i] = TMath::Sqrt(TMath::Power(totalsysthphi[i], 2.0) + TMath::Power(hhsyst[i], 2.0));
    }

    Double_t nearArraySystErr[2] = {ratios2050->GetBinContent(1)*nearsysttot[1], ratios020->GetBinContent(1)*nearsysttot[0]};
    Double_t awayArraySystErr[2] = {ratios2050->GetBinContent(3)*awaysysttot[1], ratios020->GetBinContent(3)*awaysysttot[0]};
    Double_t bulkArraySystErr[2] = {ratios2050->GetBinContent(2)*uesysttot[1], ratios020->GetBinContent(2)*uesysttot[0]};
    Double_t totalArraySystErr[2] = {ratios2050->GetBinContent(4)*totalsysttot[1], ratios020->GetBinContent(4)*totalsysttot[0]};

   
    //getting the systematics from v2 uncertainty
    TFile* v2file = new TFile("~/phiStudies/results_onlineEff/Combined/v2syst.root");
    TGraphErrors* nearv2syst = (TGraphErrors*)v2file->Get("nearv2syst");
    nearv2syst->SetLineWidth(0);
    nearv2syst->SetLineColor(kWhite);
    nearv2syst->SetFillColor(kGray+2);
    TGraphErrors* awayv2syst = (TGraphErrors*)v2file->Get("awayv2syst");
    awayv2syst->SetFillColor(kGray+2);

    //binning for 0-20, 20-80 bins
    Double_t multArray[2] = {50.0, 90.0};
    Double_t multArrayErr[2] = {30.0, 10.0}; //normal errors that take up whole multiplicity range
    Double_t multArraySystErr[2] = {2.5, 2.5}; //box systematic errors that are small
  
    Double_t mult2Array[2] = {51.0, 91.0};
    Double_t mult2ArrayErr[2] = {30.0, 10.0};
    
    //sepearte binning for 0-50, 50-80
    /*Double_t multArray[2] = {35, 75.0};
    Double_t multArrayErr[2] = {15.0, 25.0}; //normal errors that take up whole multiplicity range
    Double_t multArraySystErr[2] = {2.5, 2.5}; //box systematic errors that are small
  
    Double_t mult2Array[3] = {36.0, 76.0};
    Double_t mult2ArrayErr[3] = {15.0, 25.0};
    */

    //jet vs UE ratios
    TGraph* jetratioshPhi = new TGraphErrors(2, multArray, jet2UEhPhi);
    jetratioshPhi->SetMarkerStyle(20);
    jetratioshPhi->SetMarkerSize(2);
    jetratioshPhi->SetMarkerColor(kCyan+1);
    jetratioshPhi->SetLineColor(kCyan+2);
    jetratioshPhi->SetLineWidth(2);
    jetratioshPhi->GetXaxis()->SetTitle("Multiplicity Percentile");
    jetratioshPhi->GetXaxis()->SetTitleSize(0.05);
    jetratioshPhi->GetXaxis()->SetLabelSize(0.04);
    jetratioshPhi->GetXaxis()->SetTitleOffset(0.9);
    jetratioshPhi->GetXaxis()->SetRangeUser(0.0, 100.0);
    jetratioshPhi->GetYaxis()->SetTitle("Yield Ratio #left(#frac{Near-Side}{Underlying Event}#right)");
    jetratioshPhi->GetYaxis()->SetTitleSize(0.04);
    jetratioshPhi->GetYaxis()->SetTitleOffset(1.5); 
    jetratioshPhi->GetYaxis()->SetRangeUser(0.0002, 0.0035);

    TGraph* jetratioshh = new TGraphErrors(2, multArray, jet2UEhh);
    jetratioshh->SetMarkerStyle(20);
    jetratioshh->SetMarkerSize(2);
    jetratioshh->SetMarkerColor(kOrange+1);
    jetratioshh->SetLineColor(kOrange+2);
    jetratioshh->SetLineWidth(2);
    jetratioshh->GetXaxis()->SetTitle("Multiplicity Percentile");
    jetratioshh->GetXaxis()->SetTitleSize(0.05);
    jetratioshh->GetXaxis()->SetLabelSize(0.04);
    jetratioshh->GetXaxis()->SetTitleOffset(0.9);
    jetratioshh->GetXaxis()->SetRangeUser(0.0, 100.0);
    jetratioshh->GetYaxis()->SetTitle("Yield Ratio #left(#frac{Near-Side}{Underlying Event}#right)");
    jetratioshh->GetYaxis()->SetTitleSize(0.04);
    jetratioshh->GetYaxis()->SetTitleOffset(1.5); 
    jetratioshh->GetYaxis()->SetRangeUser(0.0002, 0.0035);

    //trying instead with variable sized histograms:
    //Double_t binwidths[4] = {0.0, 20.0, 80.0, 100.0};
    Double_t binwidths[4] = {0.0, 20.0, 50.0, 100.0};
    TH1D* ratioNearHist = new TH1D("ratioNearHist", "", 3, binwidths);
    TH1D* ratioBulkHist = new TH1D("ratioBulkHist", "", 3, binwidths);
    TH1D* ratioJetHist = new TH1D("ratioJetHist", "", 3, binwidths);
    for(int i =0; i<2; i++){
        ratioNearHist->SetBinContent(i+2, nearArray[i]);
        ratioNearHist->SetBinError(i+2, nearArrayErr[i]);
        ratioBulkHist->SetBinContent(i+2, bulkArray[i]);
        ratioBulkHist->SetBinError(i+2, bulkArrayErr[i]);
        ratioJetHist->SetBinContent(i+2, jet2UEhPhi[i]);
    }
    ratioNearHist->SetMarkerStyle(20);
    ratioNearHist->SetMarkerSize(2);
    ratioNearHist->SetMarkerColor(kRed+1);
    ratioNearHist->SetLineColor(kRed+2);
    ratioNearHist->SetLineWidth(2);
    ratioNearHist->GetXaxis()->SetTitle("Multiplicity Percentile");
    ratioNearHist->GetXaxis()->SetTitleSize(0.05);
    ratioNearHist->GetXaxis()->SetLabelSize(0.04);
    ratioNearHist->GetXaxis()->SetTitleOffset(1.2);
    ratioNearHist->GetXaxis()->SetRangeUser(0.0, 100.0);
    ratioNearHist->GetYaxis()->SetTitle("1/#it{N}_{trig} d#it{N}_{assoc}/d#it{#Delta#varphi} Yield Ratios #left(#frac{h-#phi}{h-h}#right)");
    ratioNearHist->GetYaxis()->SetTitleSize(0.04);
    ratioNearHist->GetYaxis()->SetTitleOffset(1.5); 
    ratioNearHist->GetYaxis()->SetRangeUser(0.0, 0.025);

    ratioJetHist->SetMarkerStyle(20);
    ratioJetHist->SetMarkerSize(2);
    ratioJetHist->SetMarkerColor(kCyan+1);
    ratioJetHist->SetLineColor(kCyan+2);
    ratioJetHist->SetLineWidth(2);
    ratioJetHist->GetXaxis()->SetTitle("Multiplicity Percentile");
    ratioJetHist->GetXaxis()->SetTitleSize(0.05);
    ratioJetHist->GetXaxis()->SetLabelSize(0.04);
    ratioJetHist->GetXaxis()->SetTitleOffset(1.2);
    ratioJetHist->GetXaxis()->SetRangeUser(0.0, 100.0);
    ratioJetHist->GetYaxis()->SetTitle("Yield Ratio #left(#frac{Jet Yield}{Underlying Event}#right)");
    ratioJetHist->GetYaxis()->SetTitleSize(0.04);
    ratioJetHist->GetYaxis()->SetTitleOffset(1.5); 
    ratioJetHist->GetYaxis()->SetRangeUser(0.0, 0.6);


    ratioBulkHist->SetMarkerStyle(22);
    ratioBulkHist->SetMarkerSize(2);
    ratioBulkHist->SetMarkerColor(kGreen+2);
    ratioBulkHist->SetLineColor(kGreen+3);
    ratioBulkHist->SetLineWidth(2);
    ratioBulkHist->GetXaxis()->SetTitle("Multiplicity Percentile");
    ratioBulkHist->GetXaxis()->SetTitleSize(0.05);
    ratioBulkHist->GetXaxis()->SetLabelSize(0.04);
    ratioBulkHist->GetXaxis()->SetTitleOffset(1.2);
    ratioBulkHist->GetXaxis()->SetRangeUser(0.0, 100.0);
    ratioBulkHist->GetYaxis()->SetTitle("Yield Ratio #left(#frac{h-#phi}{h-h}#right)");
    ratioBulkHist->GetYaxis()->SetTitleSize(0.04);
    ratioBulkHist->GetYaxis()->SetTitleOffset(1.5); 
    ratioBulkHist->GetYaxis()->SetRangeUser(0.0002, 0.0035);



    TGraphErrors* ratiosNear = new TGraphErrors(2, multArray, nearArray, multArrayErr, nearArrayErr);
    ratiosNear->SetMarkerStyle(20);
    ratiosNear->SetMarkerSize(1.5);
    ratiosNear->SetMarkerColor(kRed+1);
    ratiosNear->SetLineColor(kRed+2);
    ratiosNear->SetLineWidth(2.0);
    ratiosNear->GetXaxis()->SetTitle("Multiplicity Percentile (V0A)");
    ratiosNear->GetXaxis()->SetTitleSize(0.05);
    ratiosNear->GetXaxis()->SetLabelSize(0.04);
    ratiosNear->GetXaxis()->SetTitleOffset(0.9);
    ratiosNear->GetXaxis()->SetRangeUser(0.0, 100.0);
    ratiosNear->GetYaxis()->SetTitle("Yield Ratio #left(#frac{h-#phi}{h-h}#right)");
    ratiosNear->GetYaxis()->SetTitleSize(0.04);
    ratiosNear->GetYaxis()->SetTitleOffset(1.5); 
    ratiosNear->GetYaxis()->SetRangeUser(0.0002, 0.0035);

    TGraphErrors* ratiosNearSyst = new TGraphErrors(2, multArray, nearArray, multArraySystErr, nearArraySystErr);
    ratiosNearSyst->SetMarkerStyle(20);
    ratiosNearSyst->SetMarkerSize(1);
    ratiosNearSyst->SetMarkerColor(kRed+1);
    ratiosNearSyst->SetLineColor(kRed+3);
    ratiosNearSyst->SetFillColorAlpha(kWhite, 0.0);
    ratiosNearSyst->SetLineWidth(2.0);
    ratiosNearSyst->GetXaxis()->SetTitle("Multiplicity Percentile");
    ratiosNearSyst->GetXaxis()->SetTitleSize(0.05);
    ratiosNearSyst->GetXaxis()->SetLabelSize(0.04);
    ratiosNearSyst->GetXaxis()->SetTitleOffset(0.9);
    ratiosNearSyst->GetXaxis()->SetRangeUser(0.0, 100.0);
    ratiosNearSyst->GetYaxis()->SetTitle("Yield Ratio #left(#frac{h-#phi}{h-h}#right)");
    ratiosNearSyst->GetYaxis()->SetTitleSize(0.04);
    ratiosNearSyst->GetYaxis()->SetTitleOffset(1.5); 
    ratiosNearSyst->GetYaxis()->SetRangeUser(0.0002, 0.0035);



    TGraphErrors* ratiosAway = new TGraphErrors(2, multArray, awayArray, multArrayErr, awayArrayErr);
    ratiosAway->SetMarkerStyle(21);
    ratiosAway->SetMarkerSize(1.5);
    ratiosAway->SetMarkerColor(kBlue+1);
    ratiosAway->SetLineColor(kBlue+2);
    ratiosAway->SetLineWidth(2.0);

    TGraphErrors* ratiosAwaySyst = new TGraphErrors(2, multArray, awayArray, multArraySystErr, awayArraySystErr);
    ratiosAwaySyst->SetMarkerStyle(21);
    ratiosAwaySyst->SetMarkerSize(1);
    ratiosAwaySyst->SetMarkerColor(kBlue+1);
    ratiosAwaySyst->SetLineColor(kBlue+3);
    ratiosAwaySyst->SetLineWidth(2.0);
    ratiosAwaySyst->SetFillColorAlpha(kWhite, 0.0);

    TGraphErrors* ratiosBulk = new TGraphErrors(2, multArray, bulkArray, multArrayErr, bulkArrayErr);
    ratiosBulk->SetMarkerStyle(kFullCross);
    ratiosBulk->SetMarkerSize(2.0);
    ratiosBulk->SetMarkerColor(kGreen+2);
    ratiosBulk->SetLineColor(kGreen+3);
    ratiosBulk->SetLineWidth(2.0);

    TGraphErrors* ratiosBulkSyst = new TGraphErrors(2, multArray, bulkArray, multArraySystErr, bulkArraySystErr);
    ratiosBulkSyst->SetMarkerStyle(kFullCross);
    ratiosBulkSyst->SetMarkerSize(1);
    ratiosBulkSyst->SetMarkerColor(kGreen+2);
    ratiosBulkSyst->SetLineColor(kGreen+3);
    ratiosBulkSyst->SetLineWidth(2);
    ratiosBulkSyst->SetFillColorAlpha(kWhite, .50);
    ratiosBulkSyst->SetFillStyle(0);
   
    TGraphErrors* ratiosTot = new TGraphErrors(2, multArray, totalArray, multArrayErr, totalArrayErr);
    ratiosTot->SetMarkerStyle(33);
    ratiosTot->SetMarkerSize(3);
    ratiosTot->SetMarkerColor(kMagenta+2);
    ratiosTot->SetLineColor(kMagenta+3);
    ratiosTot->SetLineWidth(2);
    ratiosTot->SetFillColor(kMagenta+1);
    ratiosTot->SetFillStyle(3144);

    TGraphErrors* ratiosTotSyst = new TGraphErrors(2, multArray, totalArray, multArraySystErr, totalArraySystErr);
    ratiosTotSyst->SetMarkerStyle(33);
    ratiosTotSyst->SetMarkerSize(3);
    ratiosTotSyst->SetMarkerColor(kMagenta+2);
    ratiosTotSyst->SetLineColor(kMagenta+3);
    ratiosTotSyst->SetLineWidth(2);
    ratiosTotSyst->SetFillColorAlpha(kWhite, 0.0);
    //ratiosTotSyst->SetFillStyle(3144);


    //setting up scaled TGraph's for comparing Efficiency and Non-efficiency corrected
   /* TH1D* ratioNearHistScaled = (TH1D*)ratioNearHist->Clone("ratioNearHistScaled");
    TGraphErrors* ratiosNearScaled = (TGraphErrors*)ratiosNear->Clone("ratiosNearScaled");
    TGraphErrors* ratiosAwayScaled = (TGraphErrors*)ratiosAway->Clone("ratiosAwayScaled");
    TGraphErrors* ratiosTotScaled = (TGraphErrors*)ratiosTot->Clone("ratiosTotScaled");
    TGraphErrors* ratiosBulkScaled = (TGraphErrors*)ratiosBulk->Clone("ratiosBulkScaled");

    ratioNearHistScaled->Scale(1.0/totalArray[2]);
    Double_t x,y;
    for(int i = 0; i<3; i++){
       ratiosNearScaled->GetPoint(i, x, y);
       ratiosNearScaled->SetPoint(i, x, y/totalArray[2]);
       ratiosNearScaled->SetPointError(i, ratiosNearScaled->GetErrorX(i), ratiosNearScaled->GetErrorY(i)/totalArray[2]);
       
       ratiosAwayScaled->GetPoint(i, x, y);
       ratiosAwayScaled->SetPoint(i, x, y/totalArray[2]);
       ratiosAwayScaled->SetPointError(i, ratiosAwayScaled->GetErrorX(i), ratiosAwayScaled->GetErrorY(i)/totalArray[2]);
       
       ratiosBulkScaled->GetPoint(i, x, y);
       ratiosBulkScaled->SetPoint(i, x, y/totalArray[2]);
       ratiosBulkScaled->SetPointError(i, ratiosBulkScaled->GetErrorX(i), ratiosBulkScaled->GetErrorY(i)/totalArray[2]);
       
       ratiosTotScaled->GetPoint(i, x, y);
       ratiosTotScaled->SetPoint(i, x, y/totalArray[2]);
       ratiosTotScaled->SetPointError(i, ratiosTotScaled->GetErrorX(i), ratiosTotScaled->GetErrorY(i)/totalArray[2]);
    }
*/
    TLegend  *ratiosMultlegend = new TLegend(0.1873, 0.696, 0.450, 0.939);
    ratiosMultlegend->SetMargin(0.35);
    ratiosMultlegend->AddEntry(ratiosBulk, "Underlying Event", "pl");
    ratiosMultlegend->AddEntry(ratiosAway, "Away-side (Jet)", "pl");
    ratiosMultlegend->AddEntry(ratiosNear, "Near-side (Jet)", "pl");
    //ratiosMultlegend->AddEntry(ratiosTot, "Total (Jet + UE)", "pl");
    ratiosMultlegend->SetLineWidth(0);

    TLegend *v2leg = new TLegend(0.209, 0.605, 0.392, 0.652);
    v2leg->AddEntry(nearv2syst, "syst unc. from v_{2}", "f");
    v2leg->SetLineWidth(0);

    TLegend *ratiosUEMultlegend = new TLegend(0.183, 0.686, 0.461, 0.928);
    ratiosUEMultlegend->SetMargin(0.35);
    ratiosUEMultlegend->AddEntry(ratiosBulk, "In U.E.", "pl");
    ratiosUEMultlegend->AddEntry(ratiosTot, "Total (Jet + UE)", "f");
    ratiosUEMultlegend->SetLineWidth(0);

    TLegend *ratiosJetMultlegend = new TLegend(0.183, 0.686, 0.461, 0.928);
    ratiosJetMultlegend->SetMargin(0.35);
    ratiosJetMultlegend->AddEntry(jetratioshh, "Jet/U.E. for (h-h)", "pl");
    ratiosJetMultlegend->AddEntry(jetratioshPhi, "Jet/U.E. for (h-#phi)", "pl");
    ratiosJetMultlegend->SetLineWidth(0);


    TPaveText *text2 = new TPaveText(0.5935, 0.6367, 0.9287, 0.7796, "NDC");
    text2->AddText("4.0 < #it{p}^{h}_{T,trig} < 8.0 GeV/#it{c}");
    text2->AddText("1.5 < #it{p}^{#phi}_{T,assoc} < 2.0 GeV/#it{c}");
    text2->SetTextSizePixels(18);
    text2->SetFillColor(kWhite);
    text2->SetBorderSize(0);
    text2->SetFillStyle(0);
    text2->SetTextFont(42);


    TPaveText *data = new TPaveText(0.6058, 0.8127, 0.9042, 0.9416, "NDC");
    data->AddText("ALICE Preliminary"); 
    data->AddText("p-Pb #sqrt{s_{NN}} = 5.02 TeV");
    data->GetLine(0)->SetTextSizePixels(32);
    data->GetLine(1)->SetTextSizePixels(24);
    data->SetBorderSize(0);
    data->SetFillColor(kWhite);
    data->SetFillStyle(0);
    data->SetTextFont(42);

    TPaveText *prelim = new TPaveText(0.5416, 0.7194, 0.9615, 0.8878, "NDC");
    prelim->AddText("ALICE Prelimenary"); 
    prelim->AddText("p-Pb #sqrt{s_{NN}} = 5.02 TeV");
    prelim->GetLine(0)->SetTextSizePixels(32);
    prelim->GetLine(1)->SetTextSizePixels(24);
    prelim->SetBorderSize(0);
    prelim->SetFillColor(kWhite);

    
    TCanvas* vsMultCanvas = new TCanvas("vsMultCanvas", "vsMultCanvas", 55, 55, 900, 600);
    vsMultCanvas->cd();
    vsMultCanvas->SetMargin(0.126, 0.05, 0.125, 0.05);
    //vsMultCanvas->SetMargin(0.4, 0.05, 0.125, 0.05); //testing margin stuff
    //TH1F* hist = ratiosNear->GetHistogram();
    gStyle->SetErrorX(0.5);
    ratioNearHist->GetYaxis()->SetRangeUser(0.0, 0.040);
    ratioNearHist->GetXaxis()->SetTitle("Multiplicity Percentile (V0A)");
    ratioNearHist->Draw("AXIS");

    ratioNearHist->GetXaxis()->SetLabelOffset(999);
    //ratioNearHist->GetXaxis()->SetTitleOffset(999);
    ratioNearHist->GetXaxis()->SetTickSize(0.0);

    //ratiosNear->Draw("P");
    gPad->Update();
    TGaxis* newaxis = new TGaxis(gPad->GetUxmax(),
            gPad->GetUymin(),
            gPad->GetUxmin(),
            gPad->GetUymin(),
            ratioNearHist->GetXaxis()->GetXmin(),
            ratioNearHist->GetXaxis()->GetXmax(),
            510,"-");
    newaxis->SetLabelOffset(-0.03);
    newaxis->SetLabelSize(0.04);
    newaxis->SetLabelFont(42);
    //newaxis->SetTitle("Multipliciy % (VOA)");
    //newaxis->SetTitleOffset(1.3);
    //newaxis->SetTextSize
    newaxis->Draw();
    gPad->Update();
    //ratiosNearSyst->Draw("5");
    //nearv2syst->Draw("2");   
    ratiosNear->Draw("P");
    //ratiosAwaySyst->Draw("5");
    //awayv2syst->Draw("2");
    ratiosAway->Draw("P"); 
   // ratiosTotSyst->Draw("5");
    //ratiosBulkSyst->Draw("5");
    ratiosBulk->Draw("P");
    ratiosTot->Draw("P");
    //ratiosTot->Draw("3");
    ratiosMultlegend->Draw();
    //v2leg->Draw();
    data->Draw();
    text2->Draw();
    //ratiosNear->Draw("PL");
    //newaxis->Draw();
    //gPad->Update();
   
   /* 
    printf("near jet ratio low: %f, near jet ratio high: %f, %%increase: %f\n", nearArray[0], nearArray[2], nearArray[2]/nearArray[0] - 1.0);
    printf("away jet ratio low: %f, away jet ratio high: %f, %%increase: %f\n", awayArray[0], awayArray[2], awayArray[2]/awayArray[0] - 1.0); 
    printf("bulk  ratio low: %f, bulk  ratio high: %f, %%increase: %f\n", bulkArray[0], bulkArray[2], bulkArray[2]/bulkArray[0] - 1.0);
    printf("total ratio low: %f, total ratio high: %f, %%increase: %f\n", totalArray[0], totalArray[2], totalArray[2]/totalArray[0] - 1.0);

    Double_t fakeTotHigh = near2tothh[2]*nearArray[0] + away2tothh[2]*awayArray[0] + UE2tothh[2]*bulkArray[2];
        
    printf("if jet was flat, high ratio would be: %f, total %% increase would be: %f\n", fakeTotHigh, fakeTotHigh/totalArray[0] - 1.0);
     */

    //Jet Ratio Canvas
    TCanvas* JetvsMultCanvas = new TCanvas("JetvsMultCanvas", "JetvsMultCanvas", 55, 55, 900, 600);
    JetvsMultCanvas->cd();
    JetvsMultCanvas->SetMargin(0.126, 0.05, 0.125, 0.05);
    gStyle->SetErrorX(0.5);
    ratioJetHist->Draw("P");

    ratioJetHist->GetXaxis()->SetLabelOffset(999);
    ratioJetHist->GetXaxis()->SetTickSize(0.0);
    gPad->Update();
    newaxis = new TGaxis(gPad->GetUxmax(),
            gPad->GetUymin(),
            gPad->GetUxmin(),
            gPad->GetUymin(),
            ratioNearHist->GetXaxis()->GetXmin(),
            ratioNearHist->GetXaxis()->GetXmax(),
            510,"-");
    newaxis->SetLabelOffset(-0.03);
    newaxis->SetLabelSize(0.04);
    newaxis->SetTextFont(42);
    newaxis->Draw();   
    jetratioshh->Draw("P L");
    jetratioshPhi->Draw("P L");
    data->Draw();
    ratiosJetMultlegend->Draw("SAME");

    //Double Ratio plots:
    Double_t doublenear020 = near020/mid020;
    Double_t doublenear020Er = doublenear020*TMath::Sqrt(TMath::Power((near020Er)/(near020), 2) + TMath::Power((mid020Er)/(mid020), 2));
    Double_t doubleaway020 = away020/mid020;
    Double_t doubleaway020Er = doubleaway020*TMath::Sqrt(TMath::Power((away020Er)/(away020), 2) + TMath::Power((mid020Er)/(mid020), 2));

    Double_t doublenear2050 = near2050/mid2050;
    Double_t doublenear2050Er = doublenear2050*TMath::Sqrt(TMath::Power((near2050Er)/(near2050), 2) + TMath::Power((mid2050Er)/(mid2050), 2));
    Double_t doubleaway2050 = away2050/mid2050;
    Double_t doubleaway2050Er = doubleaway2050*TMath::Sqrt(TMath::Power((away2050Er)/(away2050), 2) + TMath::Power((mid2050Er)/(mid2050), 2));

   
    TH1D *doubleratios020 = new TH1D("doubleratios020", "(h-#phi / h-h) Double Ratios", 2, 0, 2);
    doubleratios020->GetXaxis()->SetBinLabel(1, "#frac{near-side}{mid}");
    doubleratios020->SetBinContent(1, doublenear020);
    doubleratios020->SetBinError(1, doublenear020Er);
    doubleratios020->GetXaxis()->SetBinLabel(2, "#frac{away-side}{mid}");
    doubleratios020->SetBinContent(2, doubleaway020);
    doubleratios020->SetBinError(2, doubleaway020Er);
    doubleratios020->SetMarkerStyle(22);
    doubleratios020->SetMarkerColor(kRed+2);
    doubleratios020->SetLineColor(kRed+2);
    doubleratios020->SetLineWidth(2);
    doubleratios020->SetMarkerSize(2);

    TH1D *doubleratios2050 = new TH1D("doubleratios2050", "(h-#phi / h-h) Double Ratios", 2, 0, 2);
    doubleratios2050->GetXaxis()->SetBinLabel(1, "#frac{near-side}{mid}");
    doubleratios2050->SetBinContent(1, doublenear2050);
    doubleratios2050->SetBinError(1, doublenear2050Er);
    doubleratios2050->GetXaxis()->SetBinLabel(2, "#frac{away-side}{mid}");
    doubleratios2050->SetBinContent(2, doubleaway2050);
    doubleratios2050->SetBinError(2, doubleaway2050Er);
    doubleratios2050->SetMarkerStyle(21);
    doubleratios2050->SetMarkerColor(kBlue+2);
    doubleratios2050->SetLineColor(kBlue+2);
    doubleratios2050->SetLineWidth(2);
    doubleratios2050->SetMarkerSize(2);


    TCanvas *testDoublec = new TCanvas("testdouble", "testdouble",50, 50, 600, 600);
    testDoublec->cd();
    //doubleratios2050->GetYaxis()->SetRangeUser(0.0000, 0.0040);
    doubleratios2050->GetXaxis()->SetLabelSize(0.07);
    doubleratios2050->Draw("P SAME");
    doubleratios020->GetXaxis()->SetLimits(0.05, 2.05);
    doubleratios020->Draw("P SAME");
    //line->Draw("SAME");
    ratioslegend->Draw("SAME");

/*    printf("\n\n");
    printf(" h-Phi YIELDS  ||     0 - 20     ||    20 - 50     ||    50 - 100    ||     0 - 100    ||\n");
    printf("=========================================================================================\n");
    printf("      NEAR     ||  %E  ||  %E  ||  %E  ||  %E  ||\n", near0_20hPhiYield, near20_50hPhiYield, near50_100hPhiYield, near0_100hPhiYield);
    printf("      MID      ||  %E  ||  %E  ||  %E  ||  %E  ||\n", mid0_20hPhiYield, mid20_50hPhiYield, mid50_100hPhiYield, mid0_100hPhiYield);
    printf("      AWAY     ||  %E  ||  %E  ||  %E  ||  %E  ||\n", away0_20hPhiYield, away20_50hPhiYield, away50_100hPhiYield, away0_100hPhiYield);
    printf("      TOTAL    ||  %E  ||  %E  ||  %E  ||  %E  ||\n", total0_20hPhiYield, total20_50hPhiYield, total50_100hPhiYield, total0_100hPhiYield);

    printf("\n\n");
    printf(" h-Phi ERRORS  ||     0 - 20     ||    20 - 50     ||    50 - 100    ||     0 - 100    ||\n");
    printf("=========================================================================================\n");
    printf("      NEAR     ||  %E  ||  %E  ||  %E  ||  %E  ||\n", near0_20hPhiError, near20_50hPhiError, near50_100hPhiError, near0_100hPhiError);
    printf("      MID      ||  %E  ||  %E  ||  %E  ||  %E  ||\n", mid0_20hPhiError, mid20_50hPhiError, mid50_100hPhiError, mid0_100hPhiError);
    printf("      AWAY     ||  %E  ||  %E  ||  %E  ||  %E  ||\n", away0_20hPhiError, away20_50hPhiError, away50_100hPhiError, away0_100hPhiError);
    printf("      TOTAL    ||  %E  ||  %E  ||  %E  ||  %E  ||\n", total0_20hPhiError, total20_50hPhiError, total50_100hPhiError, total0_100hPhiError);

    printf("\n\n");
    printf(" h-Phi %%ERRORS ||     0 - 20     ||    20 - 50     ||    50 - 100    ||     0 - 100    ||\n");
    printf("=========================================================================================\n");
    printf("      NEAR     ||     %2.2f%%     ||     %2.2f%%     ||     %2.2f%%     ||     %2.2f%%     ||\n", 100*near0_20hPhiError/near0_20hPhiYield, 100*near20_50hPhiError/near20_50hPhiYield, 100*near50_100hPhiError/near50_100hPhiYield, 100*near0_100hPhiError/near0_100hPhiYield);
    printf("      MID      ||      %2.2f%%     ||      %2.2f%%     ||      %2.2f%%     ||      %2.2f%%     ||\n", 100*mid0_20hPhiError/mid0_20hPhiYield, 100*mid20_50hPhiError/mid20_50hPhiYield, 100*mid50_100hPhiError/mid50_100hPhiYield, 100*mid0_100hPhiError/mid0_100hPhiYield);
    printf("      AWAY     ||     %2.2f%%     ||     %2.2f%%     ||     %2.2f%%     ||     %2.2f%%     ||\n", 100*away0_20hPhiError/away0_20hPhiYield, 100*away20_50hPhiError/away20_50hPhiYield, 100*away50_100hPhiError/away50_100hPhiYield, 100*away0_100hPhiError/away0_100hPhiYield);
    printf("      TOTAL    ||      %2.2f%%     ||      %2.2f%%     ||      %2.2f%%     ||      %2.2f%%     ||\n", 100*total0_20hPhiError/total0_20hPhiYield, 100*total20_50hPhiError/total20_50hPhiYield, 100*total50_100hPhiError/total50_100hPhiYield, 100*total0_100hPhiError/total0_100hPhiYield);
   


    printf("\n\n");
    printf(" h-h   YIELDS  ||     0 - 20     ||    20 - 50     ||    50 - 100    ||     0 - 100    ||\n");
    printf("=========================================================================================\n");
    printf("      NEAR     ||  %E  ||  %E  ||  %E  ||  %E  ||\n", near0_20hhYield, near20_50hhYield, near50_100hhYield, near0_100hhYield);
    printf("      MID      ||  %E  ||  %E  ||  %E  ||  %E  ||\n", mid0_20hhYield, mid20_50hhYield, mid50_100hhYield, mid0_100hhYield);
    printf("      AWAY     ||  %E  ||  %E  ||  %E  ||  %E  ||\n", away0_20hhYield, away20_50hhYield, away50_100hhYield, away0_100hhYield);
    printf("      TOTAL    ||  %E  ||  %E  ||  %E  ||  %E  ||\n", total0_20hhYield, total20_50hhYield, total50_100hhYield, total0_100hhYield);

    printf("\n\n");
    printf(" h-h ERRORS  ||     0 - 20     ||    20 - 50     ||    50 - 100    ||     0 - 100    ||\n");
    printf("=========================================================================================\n");
    printf("      NEAR     ||  %E  ||  %E  ||  %E  ||  %E  ||\n", near0_20hhError, near20_50hhError, near50_100hhError, near0_100hhError);
    printf("      MID      ||  %E  ||  %E  ||  %E  ||  %E  ||\n", mid0_20hhError, mid20_50hhError, mid50_100hhError, mid0_100hhError);
    printf("      AWAY     ||  %E  ||  %E  ||  %E  ||  %E  ||\n", away0_20hhError, away20_50hhError, away50_100hhError, away0_100hhError);
    printf("      TOTAL    ||  %E  ||  %E  ||  %E  ||  %E  ||\n", total0_20hhError, total20_50hhError, total50_100hhError, total0_100hhError);
    printf("\n\n");
*/

    TH1D *yields020hPhi = new TH1D("yields020hPhiEff", "h-#phi Per Trigger Yields", 4, 0, 4);
    yields020hPhi->GetXaxis()->SetBinLabel(1, "near-side");
    yields020hPhi->SetBinContent(1, near0_20hPhiYield);
    yields020hPhi->SetBinError(1, near0_20hPhiError);
    yields020hPhi->GetXaxis()->SetBinLabel(2, "mid");
    yields020hPhi->SetBinContent(2, mid0_20hPhiYield);
    yields020hPhi->SetBinError(2, mid0_20hPhiError);
    yields020hPhi->GetXaxis()->SetBinLabel(3, "away-side");
    yields020hPhi->SetBinContent(3, away0_20hPhiYield);
    yields020hPhi->SetBinError(3, away0_20hPhiError);
    yields020hPhi->GetXaxis()->SetBinLabel(4, "total");
    yields020hPhi->SetBinContent(4, total0_20hPhiYield);
    yields020hPhi->SetBinError(4, total0_20hPhiError);
    yields020hPhi->SetMarkerStyle(22);
    yields020hPhi->SetMarkerColor(kRed+2);
    yields020hPhi->SetLineColor(kRed+2);
    yields020hPhi->SetMarkerSize(2);
    yields020hPhi->SetLineWidth(2);

    TH1D *yields2050hPhi = new TH1D("yields2050hPhiEff", "h-#phi Per Trigger Yields", 4, 0, 4);
    yields2050hPhi->GetXaxis()->SetBinLabel(1, "near-side");
    yields2050hPhi->SetBinContent(1, near20_50hPhiYield);
    yields2050hPhi->SetBinError(1, near20_50hPhiError);
    yields2050hPhi->GetXaxis()->SetBinLabel(2, "mid");
    yields2050hPhi->SetBinContent(2, mid20_50hPhiYield);
    yields2050hPhi->SetBinError(2, mid20_50hPhiError);
    yields2050hPhi->GetXaxis()->SetBinLabel(3, "away-side");
    yields2050hPhi->SetBinContent(3, away20_50hPhiYield);
    yields2050hPhi->SetBinError(3, away20_50hPhiError);
    yields2050hPhi->GetXaxis()->SetBinLabel(4, "total");
    yields2050hPhi->SetBinContent(4, total20_50hPhiYield);
    yields2050hPhi->SetBinError(4, total20_50hPhiError);
    yields2050hPhi->SetMarkerStyle(21);
    yields2050hPhi->SetMarkerColor(kBlue+2);
    yields2050hPhi->SetLineColor(kBlue+2);
    yields2050hPhi->SetMarkerSize(2);
    yields2050hPhi->SetLineWidth(2);

    

    TCanvas *yieldhPhic = new TCanvas("yieldhPhi", "yieldhPhi",50, 50, 600, 600);
    yieldhPhic->cd();
    yields2050hPhi->GetYaxis()->SetRangeUser(0.0000, 0.0040);
    yields2050hPhi->GetXaxis()->SetLabelSize(0.07);
    yields2050hPhi->Draw("P SAME");
    yields020hPhi->GetXaxis()->SetLimits(-0.05, 3.95);
    yields020hPhi->Draw("P SAME");
    line->Draw("SAME");
    ratioslegend->Draw("SAME");


    TH1D *yields020hh = new TH1D("yields020hhEff", "h-#phi Per Trigger Yields", 4, 0, 4);
    yields020hh->GetXaxis()->SetBinLabel(1, "near-side");
    yields020hh->SetBinContent(1, near0_20hhYield);
    yields020hh->SetBinError(1, near0_20hhError);
    yields020hh->GetXaxis()->SetBinLabel(2, "mid");
    yields020hh->SetBinContent(2, mid0_20hhYield);
    yields020hh->SetBinError(2, mid0_20hhError);
    yields020hh->GetXaxis()->SetBinLabel(3, "away-side");
    yields020hh->SetBinContent(3, away0_20hhYield);
    yields020hh->SetBinError(3, away0_20hhError);
    yields020hh->GetXaxis()->SetBinLabel(4, "total");
    yields020hh->SetBinContent(4, total0_20hhYield);
    yields020hh->SetBinError(4, total0_20hhError);
    yields020hh->SetMarkerStyle(22);
    yields020hh->SetMarkerColor(kRed+2);
    yields020hh->SetLineColor(kRed+2);
    yields020hh->SetMarkerSize(2);
    yields020hh->SetLineWidth(2);

    TH1D *yields2050hh = new TH1D("yields2050hhEff", "h-h Per Trigger Yields", 4, 0, 4);
    yields2050hh->GetXaxis()->SetBinLabel(1, "near-side");
    yields2050hh->SetBinContent(1, near20_50hhYield);
    yields2050hh->SetBinError(1, near20_50hhError);
    yields2050hh->GetXaxis()->SetBinLabel(2, "mid");
    yields2050hh->SetBinContent(2, mid20_50hhYield);
    yields2050hh->SetBinError(2, mid20_50hhError);
    yields2050hh->GetXaxis()->SetBinLabel(3, "away-side");
    yields2050hh->SetBinContent(3, away20_50hhYield);
    yields2050hh->SetBinError(3, away20_50hhError);
    yields2050hh->GetXaxis()->SetBinLabel(4, "total");
    yields2050hh->SetBinContent(4, total20_50hhYield);
    yields2050hh->SetBinError(4, total20_50hhError);
    yields2050hh->SetMarkerStyle(21);
    yields2050hh->SetMarkerColor(kBlue+2);
    yields2050hh->SetLineColor(kBlue+2);
    yields2050hh->SetMarkerSize(2);
    yields2050hh->SetLineWidth(2);

    TLine *linehh = new TLine(3.0, 0.0, 3.0, 10.000);
    linehh->SetLineStyle(7);
    linehh->SetLineWidth(2);

    TCanvas *yieldhhc = new TCanvas("yieldhh", "yieldhh",50, 50, 600, 600);
    yieldhhc->cd();
    yields2050hh->GetYaxis()->SetRangeUser(0.0000, 10.000);
    yields2050hh->GetXaxis()->SetLabelSize(0.07);
    yields2050hh->Draw("P SAME");
    yields020hh->GetXaxis()->SetLimits(-0.05, 3.95);
    yields020hh->Draw("P SAME");
    linehh->Draw("SAME");
    ratioslegend->Draw("SAME");



    
    TH1D *ratio = (TH1D*)hPhidphi_0_20->Clone("ratio");
    ratio->Divide(hhdphi_0_20);
    

    TLegend  *legend = new TLegend(0.3791, 0.1518, 0.8772, 0.2688);
    legend->SetMargin(0.15);
    legend->AddEntry(corrFit2, "Hadron-#phi(1020) Correlation", "l");
    legend->AddEntry(corrFit, "Hadron-hadron Correlations", "l");

    TPaveText *text = new TPaveText(0.4815, 0.7056, 0.8658, 0.8551, "NDC");
    text->AddText("ALICE");
    text->AddText("p-Pb #sqrt{s_{NN}} = 5 TeV");
    text->AddText("0%-20% Multiplicity");
    text->SetTextSizePixels(20);
    text->SetBorderSize(0);
    text->SetFillColor(kWhite);

    TPaveText *text2050 = new TPaveText(0.4815, 0.7056, 0.8658, 0.8551, "NDC");
    text2050->AddText("ALICE");
    text2050->AddText("p-Pb #sqrt{s_{NN}} = 5 TeV");
    text2050->AddText("20%-50% Multiplicity");
    text2050->SetTextSizePixels(20);
    text2050->SetBorderSize(0);
    text2050->SetFillColor(kWhite);

    TPaveText *text50100 = new TPaveText(0.4815, 0.7056, 0.8658, 0.8551, "NDC");
    text50100->AddText("ALICE");
    text50100->AddText("p-Pb #sqrt{s_{NN}} = 5 TeV");
    text50100->AddText("50%-80% Multiplicity");
    text50100->SetBorderSize(0);
    text50100->SetTextSizePixels(20);
    text50100->SetFillColor(kWhite); 
   
    Double_t dphi020syst = 0.049;
    Double_t dphi2050syst = 0.048;
    Double_t dphi5080syst = 0.068;

    Double_t hhdphi020syst = 0.035;
    Double_t hhdphi2050syst = 0.035;
    Double_t hhdphi5080syst = 0.035;


    TH1D* hPhidphi_0_20_syst = (TH1D*)hPhidphi_0_20->Clone("hPhidphi_0_20_syst");
    hPhidphi_0_20_syst->SetFillColor(kGray);
    TH1D* hPhidphi_20_50_syst = (TH1D*)hPhidphi_20_50->Clone("hPhidphi_20_50_syst");
    hPhidphi_20_50_syst->SetFillColor(kGray);

    TH1D* hhdphi_0_20_syst = (TH1D*)hhdphi_0_20->Clone("hhdphi_0_20_syst");
    hhdphi_0_20_syst->SetFillColor(kGray);
    TH1D* hhdphi_20_50_syst = (TH1D*)hhdphi_20_50->Clone("hhdphi_20_50_syst");
    hhdphi_20_50_syst->SetFillColor(kGray);

    for(int i = 1; i <= 16; i++){
        hPhidphi_0_20_syst->SetBinError(i, hPhidphi_0_20_syst->GetBinContent(i)*dphi020syst);
        hPhidphi_20_50_syst->SetBinError(i, hPhidphi_20_50_syst->GetBinContent(i)*dphi2050syst);
        hhdphi_0_20_syst->SetBinError(i, hhdphi_0_20_syst->GetBinContent(i)*hhdphi020syst);
        hhdphi_20_50_syst->SetBinError(i, hhdphi_20_50_syst->GetBinContent(i)*hhdphi2050syst);
    }

    TCanvas *c0_20 = new TCanvas("c0_20", "c0_20", 50, 50, 550, 600);
    c0_20->cd();
    c0_20->SetMargin(0.12, 0.05, 0.1, 0.05);
    //hhdphi_0_20->GetYaxis()->SetRangeUser(0.001, 0.25);
    //hhdphi_0_20->Draw("E0 X0");
    hPhidphi_0_20->GetYaxis()->SetTitle("Per Trigger (h-#phi) Pairs");
    hPhidphi_0_20_syst->Draw("E2 SAME");
    hPhidphi_0_20->Draw("E0 X0 HIST SAME");
    //corrFit->Draw("SAME");
    //corrFit2->Draw("SAME");
    //hhBG->Draw("SAME");
    hphiBG->Draw("SAME");
    text->Draw();
    text2->Draw();
    //legend->Draw();


    
    TCanvas *c0_20pp = new TCanvas("c0_20pp", "c0_20pp", 50, 50, 550, 600);
    c0_20pp->cd();
    c0_20pp->SetMargin(0.12, 0.05, 0.1, 0.05);
    //hhdphi_0_20->GetYaxis()->SetRangeUser(0.001, 0.25);
    hhdphi_0_20->GetYaxis()->SetTitle("Per Trigger (h-h) Pairs");
    hhdphi_0_20_syst->Draw("E2 SAME");
    hhdphi_0_20->Draw("E0 X0 HIST SAME");
    //hPhidphi_0_20->Draw("E0 X0 SAME");
    corrFit->Draw("SAME");
    //corrFit2->Draw("SAME");
    //hhBG->Draw("SAME");
    //hphiBG->Draw("SAME");
    text->Draw();
    text2->Draw();
    //legend->Draw();

    //"cartoon" plot showing different yield regions
    TH1D* regionHist = (TH1D*)hhdphi_0_20->Clone("regionHist");
    regionHist->SetLineWidth(3);
    regionHist->SetLineColor(kBlack);

    TH1D* nearHist = (TH1D*)regionHist->Clone("nearHist");
    nearHist->SetFillColor(kRed+2);
    TH1D* awayHist = (TH1D*)regionHist->Clone("awayHist");
    awayHist->SetFillColor(kBlue+1);
    TH1D* bulkHist = (TH1D*)regionHist->Clone("bulkHist");
    bulkHist->SetLineColor(kGreen+2);
    bulkHist->SetLineWidth(1);
    bulkHist->SetFillColor(kGreen+2);

    for(int ibin = 1; ibin <= nearHist->GetXaxis()->GetNbins(); ibin++){
        if(ibin <= nearHist->GetXaxis()->GetNbins()/2.0){
            awayHist->SetBinContent(ibin, 0.);
        }else{
            nearHist->SetBinContent(ibin, 0.);
        }
        bulkHist->SetBinContent(ibin, hhBG->GetParameter(0));
    }

    TLegend* regionLeg = new TLegend(0.5, 0.6, 0.8, 0.8);
    regionLeg->AddEntry(nearHist, "Near-side Jet", "f");
    regionLeg->AddEntry(awayHist, "Away-side Jet", "f");
    regionLeg->AddEntry(bulkHist, "Underlying Event", "f");
    TCanvas* cregions = new TCanvas("cregions", "cregions", 50, 50, 600, 600);
    cregions->cd();
    regionHist->GetYaxis()->SetRangeUser(0.0, 0.8);
    regionHist->Draw("HIST");
    nearHist->Draw("HIST SAME");
    awayHist->Draw("HIST SAME");
    bulkHist->Draw("HIST SAME");
    regionHist->Draw("HIST SAME");
    regionLeg->Draw();


    TCanvas *c20_50 = new TCanvas("c20_50", "c20_50", 50, 50, 550, 600);
    c20_50->cd();
    c20_50->SetMargin(0.12, 0.05, 0.1, 0.05);
    //hhdphi_0_20->GetYaxis()->SetRangeUser(0.03, 0.11);
    //hhdphi_20_50->Draw("E0 X0");
    hPhidphi_20_50_syst->Draw("E2 SAME");
    hPhidphi_20_50->Draw("E0 X0 HIST SAME");
    //corrFit2050->Draw("SAME");
    //corrFit2_2050->Draw("SAME");
    //hhBG_20_50->Draw("SAME");
    hphiBG_20_50->Draw("SAME");
    text2050->Draw();
    text2->Draw();
    //legend->Draw();

    TCanvas *c20_50pp = new TCanvas("c20_50pp", "c20_50pp", 50, 50, 550, 600);
    c20_50pp->cd();
    c20_50pp->SetMargin(0.12, 0.05, 0.1, 0.05);
    //hhdphi_0_20->GetYaxis()->SetRangeUser(0.03, 0.11);
    hhdphi_20_50_syst->Draw("E2 SAME");
    hhdphi_20_50->Draw("E0 X0 HIST SAME");
    //hPhidphi_20_50->Draw("E0 X0 SAME");
    //corrFit2050->Draw("SAME");
    //corrFit2_2050->Draw("SAME");
    hhBG_20_50->Draw("SAME");
    //hphiBG_20_50->Draw("SAME");
    text2050->Draw();
    text2->Draw();
    //legend->Draw();

    //Set-up and Draw the hh and hphi correlations on single canvases
    
    TH1D* hphi020 = (TH1D*)hPhidphi_0_20->Clone("hphi020");
    hphi020->SetMarkerColor(kRed);
    hphi020->SetMarkerSize(1.5);
    hphi020->SetMarkerStyle(kFullCross);
    hphi020->SetLineColor(kRed);
    hphi020->Add(hphiBG, -1.0);
    hphi020->Scale(1.0/(2.4*hphi020->GetXaxis()->GetBinWidth(1))); //scale by delta-eta and bin width

    TH1D* hphi020syst = (TH1D*)hPhidphi_0_20_syst->Clone("hphi020syst");
    hphi020syst->Add(hphiBG, -1.0);
    hphi020syst->Scale(1.0/(2.4*hphi020->GetXaxis()->GetBinWidth(1)));
    hphi020syst->SetFillColor(kRed);
    hphi020syst->SetFillStyle(3002);
    hphi020syst->SetMarkerSize(0);

    TH1D* hphi2050 = (TH1D*)hPhidphi_20_50->Clone("hphi2050");
    hphi2050->SetMarkerColor(kOrange+1);
    hphi2050->SetLineColor(kOrange+1);
    hphi2050->SetMarkerSize(1.5);
    hphi2050->SetMarkerStyle(kFullSquare);
    hphi2050->Add(hphiBG_20_50, -1.0);
    hphi2050->Scale(1.0/(2.4*hphi2050->GetXaxis()->GetBinWidth(1))); //scale by delta-eta and bin width

    TH1D* hphi2050syst = (TH1D*)hPhidphi_20_50_syst->Clone("hphi2050syst");
    hphi2050syst->Add(hphiBG_20_50, -1.0);
    hphi2050syst->Scale(1.0/(2.4*hphi2050->GetXaxis()->GetBinWidth(1)));
    hphi2050syst->SetFillColor(kOrange+1);
    hphi2050syst->SetFillStyle(3002);
    hphi2050syst->SetMarkerSize(0);


    TH1D* hh020 = (TH1D*)hhdphi_0_20->Clone("hh020");
    hh020->SetMarkerColor(kRed+1);
    hh020->SetMarkerSize(1.5);
    hh020->SetMarkerStyle(kOpenCross);
    hh020->SetLineColor(kRed+1);
    hh020->Add(hhBG, -1.0);
    hh020->Scale(1.0/(2.4*hh020->GetXaxis()->GetBinWidth(1))); //scale by delta-eta and bin width
    
    TH1D* hh020syst = (TH1D*)hhdphi_0_20_syst->Clone("hh020syst");
    hh020syst->Add(hhBG, -1.0);
    hh020syst->Scale(1.0/(2.4*hh020->GetXaxis()->GetBinWidth(1)));
    hh020syst->SetFillColor(kRed+1);
    hh020syst->SetFillStyle(3002);
    hh020syst->SetMarkerSize(0);

    TH1D* hh2050 = (TH1D*)hhdphi_20_50->Clone("hh2050");
    hh2050->SetMarkerColor(kOrange+2);
    hh2050->SetMarkerSize(1.5);
    hh2050->SetMarkerStyle(kOpenSquare);
    hh2050->SetLineColor(kOrange+2);
    hh2050->Add(hhBG_20_50, -1.0);
    hh2050->Scale(1.0/(2.4*hh2050->GetXaxis()->GetBinWidth(1))); //scale by delta-eta and bin width

    TH1D* hh2050syst = (TH1D*)hhdphi_20_50_syst->Clone("hh2050syst");
    hh2050syst->Add(hhBG_20_50, -1.0);
    hh2050syst->Scale(1.0/(2.4*hh2050->GetXaxis()->GetBinWidth(1)));
    hh2050syst->SetFillColor(kOrange+2);
    hh2050syst->SetFillStyle(3002);
    hh2050syst->SetMarkerSize(0);

    TF1* fline = new TF1("fline", "pol0", -0.5*TMath::Pi(), 1.5*TMath::Pi());
    fline->SetParameter(0, 0.0);
    fline->SetLineColor(kBlack);
    fline->SetLineWidth(3);
    fline->SetLineStyle(7);

    TLegend* hphileg = new TLegend(0.456, 0.641, 0.613, 0.937);
    hphileg->SetBorderSize(0);
    hphileg->AddEntry(hphi020, "V0A 0-20%", "lep");
    hphileg->AddEntry(hphi2050, "V0A 20-50%", "lep");

    TLegend* hhleg = new TLegend(0.456, 0.641, 0.613, 0.937);
    hhleg->SetBorderSize(0);
    hhleg->AddEntry(hh020, "V0A 0-20%", "lep");
    hhleg->AddEntry(hh2050, "V0A 20-50%", "lep");
    
    TPaveText* hphitext020 = new TPaveText(0.181, 0.827, 0.356, 0.925, "NDC");
    hphitext020->SetBorderSize(0);
    hphitext020->SetFillColor(kWhite);
    hphitext020->AddText("h-#phi");
    hphitext020->SetTextFont(42);
    TPaveText* hphitext2050 = new TPaveText(0.181, 0.827, 0.356, 0.925, "NDC");
    hphitext2050->SetBorderSize(0);
    hphitext2050->SetFillColor(kWhite);
    hphitext2050->AddText("h-#phi");
    hphitext2050->SetTextFont(42);
    TPaveText* hphitext5080 = new TPaveText(0.181, 0.827, 0.356, 0.925, "NDC");
    hphitext5080->SetBorderSize(0);
    hphitext5080->SetFillColor(kWhite);
    hphitext5080->AddText("h-#phi");

    TPaveText* hhtext020 = new TPaveText(0.2147, 0.7856, 0.3698, 0.8817, "NDC");
    hhtext020->SetBorderSize(0);
    hhtext020->SetFillColor(kWhite);
    hhtext020->AddText("h-h");
    //hhtext020->SetTextSize(36);
    hhtext020->SetTextFont(42);
    TPaveText* hhtext2050 = new TPaveText(0.2147, 0.7856, 0.3698, 0.8817, "NDC");
    hhtext2050->SetBorderSize(0);
    hhtext2050->SetFillColor(kWhite);
    hhtext2050->AddText("h-h");
    TPaveText* hhtext5080 = new TPaveText(0.2147, 0.7856, 0.3698, 0.8817, "NDC");
    hhtext5080->SetBorderSize(0);
    hhtext5080->SetFillColor(kWhite);
    hhtext5080->AddText("h-h");


    TCanvas* chphi = new TCanvas("chphi", "chphi", 50, 50, 550, 600);
    chphi->cd();
    chphi->SetMargin(0.12, 0.05, 0.1, 0.05);
    hphi020->GetYaxis()->SetTitle("1/N_{trig} dN/d#Delta#varphi per #Delta#eta - constant (rad^{-1})");
    hphi020->GetYaxis()->SetRangeUser(-0.5E-3, 2.0E-3);
    hphi020->Draw("E0 X0 P SAME");
    hphi2050->GetYaxis()->SetRangeUser(-0.5E-3, 2.0E-3);
    hphi2050->Draw("E0 X0 P SAME");
    fline->Draw("SAME");
    hphitext020->Draw();
    text2->Draw();
    hphileg->Draw();

    TCanvas* chh = new TCanvas("chh", "chh", 50, 50, 550, 600);
    chh->cd();
    chh->SetMargin(0.12, 0.05, 0.1, 0.05);
    hh020->Draw("E0 X0 P SAME");
    hh2050->Draw("E0 X0 P SAME");
    fline->Draw("SAME");
    hhtext020->Draw();
    text2->Draw();
    hhleg->Draw();


    TCanvas* chphi020 = new TCanvas("chphi020", "chphi020", 50, 50, 550, 600);
    chphi020->cd();
    chphi020->SetMargin(0.12, 0.05, 0.1, 0.05);
    hphi020syst->GetYaxis()->SetTitle("1/#it{N}_{trig} d#it{N}/d#Delta#varphi per #Delta#eta - constant (rad^{-1})");
    hphi020syst->GetYaxis()->SetRangeUser(-0.6E-3, 2.8E-3);
    hphi020syst->Draw("E2");
    hphi020->Draw("E0 X0 P SAME");
    fline->Draw("SAME");
    data->Draw();
    text2->Draw();
    hphitext020->Draw();
    

    TCanvas* chphi2050 = new TCanvas("chphi2050", "chphi2050", 50, 50, 550, 600);
    chphi2050->cd();
    chphi2050->SetMargin(0.12, 0.05, 0.1, 0.05);
    hphi2050syst->GetYaxis()->SetTitle("1/#it{N}_{trig} d#it{N}/d#Delta#varphi per #Delta#eta - constant (rad^{-1})");
    hphi2050syst->GetYaxis()->SetRangeUser(-0.6E-3, 2.8E-3);
    hphi2050syst->Draw("E2");
    hphi2050->Draw("E0 X0 P SAME");
    fline->Draw("SAME");
    data->Draw();
    hphitext2050->Draw();


    TCanvas* chh020 = new TCanvas("chh020", "chh020", 50, 50, 550, 600);
    chh020->cd();
    chh020->SetMargin(0.12, 0.05, 0.1, 0.05);
    hh020syst->GetYaxis()->SetTitle("1/#it{N}_{trig} d#it{N}/d#Delta#varphi per #Delta#eta - constant (rad^{-1})");
    hh020syst->GetYaxis()->SetRangeUser(-50E-3, 280E-3);
    hh020syst->Draw("E2");
    hh020->Draw("E0 X0 P SAME");
    fline->Draw("SAME");
    data->Draw();
    text2->Draw();
    hhtext020->Draw();
    

    TCanvas* chh2050 = new TCanvas("chh2050", "chh2050", 50, 50, 550, 600);
    chh2050->cd();
    chh2050->SetMargin(0.12, 0.05, 0.1, 0.05);
    hh2050syst->GetYaxis()->SetTitle("1/N_{trig} dN/d#Delta#varphi per #Delta#eta - constant (rad^{-1})");
    hh2050syst->GetYaxis()->SetRangeUser(-50E-3, 280E-3);
    hh2050syst->Draw("E2");
    hh2050->Draw("E0 X0 P SAME");
    fline->Draw("SAME");
    data->Draw();
    hhtext2050->Draw();



    TCanvas* chphiall = new TCanvas("chphiall", "chphiall", 50, 50, 1200, 525);
    //chphiall->Divide(3,1,0,0);
    //chphiall->cd(1);
    chphiall->cd()->SetMargin(0, 0, 0, 0);
    TPad* hphi020pad = new TPad("hphi020pad", "", 0, 0, 0.38, 1.0);
    hphi020pad->SetMargin(0.15, 0.0, 0.1, 0.1);
    hphi020pad->Draw();
    hphi020pad->cd();
    hphi020syst->GetYaxis()->SetTitleSize(0.05);
    hphi020syst->GetYaxis()->SetLabelSize(0.05);
    hphi020syst->GetYaxis()->SetTitleOffset(1.5);
    hphi020syst->GetXaxis()->SetLabelSize(0.05);
    hphi020syst->GetYaxis()->SetRangeUser(-0.5E-3, 2.0E-3);
    hphi020syst->Draw("E2");
    hphi020->Draw("E0 X0 P SAME");
    fline->Draw("SAME");
    data->Draw();
    //text2->Draw();
    hphitext020->Draw();
    chphiall->cd();
    TPad* hphi2050pad = new TPad("hphi2050pad", "", 0.38, 0, 0.68, 1.0);
    hphi2050pad->SetMargin(0.0, 0.0, 0.1, 0.1);
    hphi2050pad->Draw();
    hphi2050pad->cd();
    hphi2050syst->GetXaxis()->SetLabelSize(0.05);
    hphi2050syst->GetYaxis()->SetLabelSize(0.0);
    hphi2050syst->GetYaxis()->SetRangeUser(-0.5E-3, 2.0E-3);
    hphi2050syst->Draw("E2");
    hphi2050->Draw("E0 X0 P SAME");
    fline->Draw("SAME");
    text2->Draw();
    //data->Draw();
    //hphitext2050->Draw();
    chphiall->cd();

    TCanvas* chhall = new TCanvas("chhall", "chhall", 50, 50, 1200, 525);
    chhall->cd()->SetMargin(0, 0, 0.0, 0.0);
    //chhall->Divide(3,1,0.0,0.1);
    //chhall->cd(1)->SetMargin(0.15, 0.0, 0.1, 0.1);
    TPad* hh020pad = new TPad("hh020pad", "", 0, 0, 0.38, 1.0);
    hh020pad->SetMargin(0.15, 0.0, 0.1, 0.1);
    hh020pad->Draw();
    hh020pad->cd();
    hh020syst->GetYaxis()->SetTitleSize(0.05);
    hh020syst->GetYaxis()->SetLabelSize(0.05);
    hh020syst->GetYaxis()->SetTitleOffset(1.5);
    hh020syst->GetXaxis()->SetLabelSize(0.05);
    hh020syst->Draw("E2");
    hh020->Draw("E0 X0 P SAME");
    fline->Draw("SAME");
    data->Draw();
    //text2->Draw();
    hhtext020->Draw();
    //chhall->cd(2)->SetMargin(0.0,0.0, 0.1, 0.1);
    chhall->cd();
    TPad* hh2050pad = new TPad("hh2050pad", "", 0.38, 0, 0.68, 1.0);
    hh2050pad->SetMargin(0.0, 0.0, 0.1, 0.1);
    hh2050pad->Draw();
    hh2050pad->cd();
    //hh2050syst->GetXaxis()->SetLabelSize(0.06);
    //hh2050syst->GetXaxis()->SetLabelOffset(0.7);
    hh2050syst->GetXaxis()->SetTitleSize(0.065);
    //hh2050syst->GetXaxis()->SetTitleOffset(-0.003);
    hh2050syst->GetYaxis()->SetLabelSize(0.0);
    hh2050syst->Draw("E2");
    hh2050->Draw("E0 X0 P SAME");
    fline->Draw("SAME");
    //data->Draw();
    text2->Draw();
    //hhtext2050->Draw();
    //chhall->cd(3)->SetMargin(0.0, 0.1, 0.1, 0.1);
    chhall->cd();


 

    TH1F* hhyieldSyst_0_20 = new TH1F("hhyieldSyst_0_20", "hhyieldSyst_0_20", 500, 0.75, 1.25);
    TH1F* hhyieldSyst_20_50 = new TH1F("hhyieldSyst_20_50", "hhyieldSyst_20_50", 500, 0.75, 1.25);
    TH1F* hphiyieldSyst_0_20 = new TH1F("hphiyieldSyst_0_20", "hphiyieldSyst_0_20", 500, 0.75, 1.25);
    TH1F* hphiyieldSyst_20_50 = new TH1F("hphiyieldSyst_20_50", "hphiyieldSyst_20_50", 500, 0.75, 1.25);
    
    TH1F* hhnearyieldSyst_0_20 = new TH1F("hhnearyieldSyst_0_20", "hhnearyieldSyst_0_20", 500, 0.75, 1.25);
    TH1F* hhnearyieldSyst_20_50 = new TH1F("hhnearyieldSyst_20_50", "hhnearyieldSyst_20_50", 500, 0.75, 1.25);
    TH1F* hphinearyieldSyst_0_20 = new TH1F("hphinearyieldSyst_0_20", "hphinearyieldSyst_0_20", 500, 0.75, 1.25);
    TH1F* hphinearyieldSyst_20_50 = new TH1F("hphinearyieldSyst_20_50", "hphinearyieldSyst_20_50", 500, 0.75, 1.25);
   
    TH1F* hhawayyieldSyst_0_20 = new TH1F("hhawayyieldSyst_0_20", "hhawayyieldSyst_0_20", 500, 0.75, 1.25);
    TH1F* hhawayyieldSyst_20_50 = new TH1F("hhawayyieldSyst_20_50", "hhawayyieldSyst_20_50", 500, 0.75, 1.25);
    TH1F* hphiawayyieldSyst_0_20 = new TH1F("hphiawayyieldSyst_0_20", "hphiawayyieldSyst_0_20", 500, 0.75, 1.25);
    TH1F* hphiawayyieldSyst_20_50 = new TH1F("hphiawayyieldSyst_20_50", "hphiawayyieldSyst_20_50", 500, 0.75, 1.25);

    TH1F* hhbulkyieldSyst_0_20 = new TH1F("hhbulkyieldSyst_0_20", "hhbulkyieldSyst_0_20", 500, 0.75, 1.25);
    TH1F* hhbulkyieldSyst_20_50 = new TH1F("hhbulkyieldSyst_20_50", "hhbulkyieldSyst_20_50", 500, 0.75, 1.25);
    TH1F* hphibulkyieldSyst_0_20 = new TH1F("hphibulkyieldSyst_0_20", "hphibulkyieldSyst_0_20", 500, 0.75, 1.25);
    TH1F* hphibulkyieldSyst_20_50 = new TH1F("hphibulkyieldSyst_20_50", "hphibulkyieldSyst_20_50", 500, 0.75, 1.25);

    Float_t stdtotal0_20 = 5.850084E-02;
    Float_t stdtotal20_50 = 3.489754E-02;
    Float_t stdnear0_20 = 1.509740E-03;
    Float_t stdnear20_50 = 1.450373E-03;
    Float_t stdaway0_20 = 2.234562E-03;
    Float_t stdaway20_50 = 2.274811E-03;
    Float_t stdmid0_20 = 5.476692E-02;
    Float_t stdmid20_50 = 3.111197E-02;

    
    
    hhyieldSyst_0_20->Fill(total0_20hhYield);
    hhyieldSyst_20_50->Fill(total20_50hhYield);
    hhnearyieldSyst_0_20->Fill(near0_20hhYield);
    hhnearyieldSyst_20_50->Fill(near20_50hhYield);
    hhawayyieldSyst_0_20->Fill(away0_20hhYield);
    hhawayyieldSyst_20_50->Fill(away20_50hhYield);
    hhbulkyieldSyst_0_20->Fill(mid0_20hhYield);
    hhbulkyieldSyst_20_50->Fill(mid20_50hhYield);

    hphiyieldSyst_0_20->Fill(total0_20hPhiYield/stdtotal0_20);
    hphiyieldSyst_20_50->Fill(total20_50hPhiYield/stdtotal20_50);
    hphinearyieldSyst_0_20->Fill(near0_20hPhiYield/stdnear0_20);
    hphinearyieldSyst_20_50->Fill(near20_50hPhiYield/stdnear20_50);
    hphiawayyieldSyst_0_20->Fill(away0_20hPhiYield/stdaway0_20);
    hphiawayyieldSyst_20_50->Fill(away20_50hPhiYield/stdaway20_50);
    hphibulkyieldSyst_0_20->Fill(mid0_20hPhiYield/stdmid0_20);
    hphibulkyieldSyst_20_50->Fill(mid20_50hPhiYield/stdmid20_50);


    TFile* output = new TFile(Form("fitsyst_%s.root", outputstring.Data()), "RECREATE");

    hhyieldSyst_0_20->Write();
    hhyieldSyst_20_50->Write();
    hhnearyieldSyst_0_20->Write();
    hhnearyieldSyst_20_50->Write();
    hhawayyieldSyst_0_20->Write();
    hhawayyieldSyst_20_50->Write();
    hhbulkyieldSyst_0_20->Write();
    hhbulkyieldSyst_20_50->Write();

    hhdphi_0_20->Write();
    hPhidphi_0_20->Write();
    hhdphi_20_50->Write();
    hPhidphi_20_50->Write();
    

    hphiyieldSyst_0_20->Write();
    hphiyieldSyst_20_50->Write();
    hphinearyieldSyst_0_20->Write();
    hphinearyieldSyst_20_50->Write();
    hphiawayyieldSyst_0_20->Write();
    hphiawayyieldSyst_20_50->Write();
    hphibulkyieldSyst_0_20->Write();
    hphibulkyieldSyst_20_50->Write();

    ratioNearHist->Write();
    ratiosNear->Write("ratiosNear");
    ratiosAway->Write("ratiosAway");
    ratiosBulk->Write("ratiosBulk");
    ratiosTot->Write("ratiosTot");

  
}
