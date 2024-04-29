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
  gStyle->SetLabelSize(0.05,"xyz");
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
  gPad->SetTickx();
  gPad->SetTicky();

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
    histo1D->GetXaxis()->SetTitle("#Delta#varphi (rad)");
    histo1D->GetYaxis()->SetTitle("1/#it{N}_{trig} d#it{N}/d#Delta#varphi per #Delta#eta (rad^{-1})");
    histo1D->SetTitle("");
    histo1D->GetYaxis()->SetTitleOffset(1.60);
    histo1D->GetYaxis()->SetMaxDigits(2);
    histo1D->GetYaxis()->SetLabelSize(0.05);
    histo1D->GetYaxis()->SetLabelOffset(0.005);
    histo1D->GetYaxis()->SetDecimals(kTRUE);
    histo1D->GetXaxis()->SetTitleSize(0.06);
    histo1D->GetXaxis()->SetTitleOffset(0.80);
    histo1D->GetXaxis()->SetLabelSize(0.05);
    histo1D->GetXaxis()->SetLabelOffset(0.005);
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
            basefit->SetParameter(0, hist->GetBinContent(1));
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

Double_t calcJetFromBG(TH1D* h, Int_t firstbin, Int_t lastbin, Float_t BG, Double_t& error){
    Double_t sum = 0.0;
    error = 0.0;
    for(int ibin = firstbin; ibin <= lastbin; ibin++){
        Float_t pt = h->GetBinContent(ibin);
        if(pt >= BG){
            sum += pt - BG;
            error += TMath::Power(h->GetBinError(ibin), 2);
            //printf("ibin = %d, sum = %f\n", ibin, sum);
        }
    }
    error = TMath::Sqrt(error);
    return sum;
}

Double_t calcJetFromBG(TH1D* h, Int_t firstbin, Int_t lastbin, TF1* BG, Double_t& error){
    Double_t sum = 0.0;
    error = 0.0;
    for(int ibin = firstbin; ibin <= lastbin; ibin++){
        Float_t pt = h->GetBinContent(ibin);
        if(pt >= BG->Eval(h->GetBinCenter(ibin))){
            sum += pt - BG->Eval(h->GetBinCenter(ibin));
            error += TMath::Power(h->GetBinError(ibin), 2);
            //printf("ibin = %d, sum = %f\n", ibin, sum);
        }
    }
    error = TMath::Sqrt(error);
    return sum;

}

Double_t calcUE(TH1D* h, Int_t firstbin, Int_t lastbin, Float_t BG, Double_t& error){
    Double_t sum = 0.0;
    error = 0;
    for(int ibin = firstbin; ibin <= lastbin; ibin++){
        Float_t pt = h->GetBinContent(ibin);
        if(pt > BG){
            sum+=BG;
        }else{
            sum+=pt;
            error+=TMath::Power(h->GetBinError(ibin), 2);
        }
    }
    error = TMath::Sqrt(error);
    return sum;
}

Double_t calcUE(TH1D* h, Int_t firstbin, Int_t lastbin, TF1* BG, Double_t& error){
    Double_t sum = 0.0;
    error = 0;
    for(int ibin = firstbin; ibin <= lastbin; ibin++){
        Float_t pt = h->GetBinContent(ibin);
        if(pt > BG->Eval(h->GetBinCenter(ibin))){
            sum+=BG->Eval(h->GetBinCenter(ibin));
        }else{
            sum+=pt;
            error+=TMath::Power(h->GetBinError(ibin), 2);
        }
    }
    error = TMath::Sqrt(error);
    return sum;
}


Double_t ZYAM(TH1D* h){
    Double_t min = 0;
    for(int ibin = 1; ibin <= h->GetNbinsX(); ibin++){
        Double_t bc = h->GetBinContent(ibin);
        if(ibin == 1){
            min = bc;
        }else if(bc < min){
            min = bc;
        }
    }
    return min;
}


void intUSRatioPlotsystNewBase(TString outputstring, TString input_0_20 = "", TString input_20_50 = "", TString input_50_80 = "", TString scaleSuffix = "avgscale", Int_t bgMethod=1){

    bool IS_ZYAM = false;
    bool IS_AVERAGE_BG = false;

    gStyle->SetOptStat(0);
    gStyle->SetOptFit(0);
    gStyle->SetErrorX(0);
    gStyle->SetPadColor(10);
    gStyle->SetPadTickX(1);
    gStyle->SetPadTickY(1);


    if(input_0_20.EqualTo("")) input_0_20 = "~/phiStudies/results_onlineEff/Combined/US_syst_trig_4_8_assoc_2_4_effcorr_hPhi_0_20_combined.root";
    if(input_20_50.EqualTo("")) input_20_50 = "~/phiStudies/results_onlineEff/Combined/US_syst_trig_4_8_assoc_2_4_effcorr_hPhi_20_50_combined.root";
    if(input_50_80.EqualTo("")) input_50_80 = "~/phiStudies/results_onlineEff/Combined/US_syst_trig_4_8_assoc_2_4_effcorr_hPhi_50_80_combined.root";

    // 2019 approved results
    TString hhinput_0_20 = "~/phiStudies/results_onlineEff/Combined/trig_4_8_assoc_2_4_effcorr_hh_0_20.root";
    TString hhinput_20_50 = "~/phiStudies/results_onlineEff/Combined/trig_4_8_assoc_2_4_effcorr_hh_20_50.root";
    TString hhinput_50_80 = "~/phiStudies/results_onlineEff/Combined/trig_4_8_assoc_2_4_effcorr_hh_50_80.root";

    //2021 new code Combined
    
    input_0_20 = "~/phiStudies/results_onlineEff/Combined/2021code/US_syst_trig_4.0_8.0_assoc_2.0_4.0_effcorr_hPhi_0_20_2021code.root";
    input_20_50 = "~/phiStudies/results_onlineEff/Combined/2021code/US_syst_trig_4.0_8.0_assoc_2.0_4.0_effcorr_hPhi_20_50_2021code.root";
    input_50_80 = "~/phiStudies/results_onlineEff/Combined/2021code/US_syst_trig_4.0_8.0_assoc_2.0_4.0_effcorr_hPhi_50_80_moremix.root";
    hhinput_0_20 = "~/phiStudies/results_onlineEff/Combined/2021code/trig_4.0_8.0_assoc_2.0_4.0_effcorr_hh_0_20.root";
    hhinput_20_50 = "~/phiStudies/results_onlineEff/Combined/2021code/trig_4.0_8.0_assoc_2.0_4.0_effcorr_hh_20_50.root";
    hhinput_50_80 = "~/phiStudies/results_onlineEff/Combined/2021code/trig_4.0_8.0_assoc_2.0_4.0_effcorr_hh_50_80.root";

    //highest trigger test
//    input_0_20 = "~/phiStudies/results_onlineEff/FAST/highestTrig/US_syst_trig_4.0_8.0_assoc_2.0_4.0_effcorr_hPhi_0_20_hightrig.root";
//    input_20_50 = "~/phiStudies/results_onlineEff/FAST/highestTrig/US_syst_trig_4.0_8.0_assoc_2.0_4.0_effcorr_hPhi_20_50_hightrig.root";
//    input_50_80 = "~/phiStudies/results_onlineEff/FAST/highestTrig/US_syst_trig_4.0_8.0_assoc_2.0_4.0_effcorr_hPhi_50_80_hightrig.root";
//    hhinput_0_20 = "~/phiStudies/results_onlineEff/FAST/highestTrig/trig_4.0_8.0_assoc_2.0_4.0_effcorr_hh_0_20.root";
//    hhinput_20_50 = "~/phiStudies/results_onlineEff/FAST/highestTrig/trig_4.0_8.0_assoc_2.0_4.0_effcorr_hh_20_50.root";
//    hhinput_50_80 = "~/phiStudies/results_onlineEff/FAST/highestTrig/trig_4.0_8.0_assoc_2.0_4.0_effcorr_hh_50_80.root";


    double refoffset020 = 1.0;
    double refoffset2050 = 1.0;
    double refoffset5080 = 1.0;
    double refhhoffset020 = 1.0;
    double refhhoffset2050 = 1.0;
    double refhhoffset5080 = 1.0;

    double refoffset020avg = 1.0;
    double refoffset2050avg = 1.0;
    double refoffset5080avg = 1.0;
    double refhhoffset020avg = 1.0;
    double refhhoffset2050avg = 1.0;
    double refhhoffset5080avg = 1.0;

    double refoffset020std = 1.0;
    double refoffset2050std = 1.0;
    double refoffset5080std = 1.0;
    double refhhoffset020std = 1.0;
    double refhhoffset2050std = 1.0;
    double refhhoffset5080std = 1.0;


    //old online code with new offline code
/*    
    input_0_20 = "~/phiStudies/results_onlineEff/Combined/oldONnewOFF/US_syst_trig_4.0_8.0_assoc_2.0_4.0_effcorr_hPhi_0_20_oldONnewOFF.root";
    input_20_50 = "~/phiStudies/results_onlineEff/Combined/oldONnewOFF/US_syst_trig_4.0_8.0_assoc_2.0_4.0_effcorr_hPhi_20_50_oldONnewOFF.root";
    input_50_80 = "~/phiStudies/results_onlineEff/Combined/oldONnewOFF/US_syst_trig_4.0_8.0_assoc_2.0_4.0_effcorr_hPhi_50_80_oldONnewOFF.root";
    hhinput_0_20 = "~/phiStudies/results_onlineEff/Combined/oldONnewOFF/trig_4.0_8.0_assoc_2.0_4.0_effcorr_hh_0_20.root";
    hhinput_20_50 = "~/phiStudies/results_onlineEff/Combined/oldONnewOFF/trig_4.0_8.0_assoc_2.0_4.0_effcorr_hh_20_50.root";
    hhinput_50_80 = "~/phiStudies/results_onlineEff/Combined/oldONnewOFF/trig_4.0_8.0_assoc_2.0_4.0_effcorr_hh_50_80.root";
*/   


    TString assocpt = "2.0 < #it{p}^{#phi}_{T,assoc} < 4.0 GeV/#it{c}";
    TString hhassocpt = "2.0 < #it{p}^{h}_{T,assoc} < 4.0 GeV/#it{c}";
    TString genassocpt = "2.0 < #it{p}_{T,assoc} < 4.0 GeV/#it{c}";


    //pythia ratio from LF train 708 & 709 (LHC16_17_general_purpose)
    //values from 2-4 range {near, away, UE, total}
//    Float_t pythiaval[4] = {0.008781, 0.011406, 0.017875, 0.014989};
//    Float_t pythiazyam[4] = {0.008949, 0.011575, 0.017927, 0.014989};

    
    Float_t pythianch[1] = {12.3};
    Float_t pythianearnch[1] = { 0.008781};
    Float_t pythianearerrnch[1] = {0.000428};
    Float_t pythiaawaynch[1] = {0.011406};
    Float_t pythiaawayerrnch[1] = {0.000843};
    Float_t pythiatotnch[1] = {0.014989};
    Float_t pythiatoterrnch[1] = {0.000203};
    Float_t pythiaUEnch[1] = {0.017875};
    Float_t pythiaUEerrnch[1] = {0.000203};

    Float_t pythiahphinear[1] = {0.902};
    Float_t pythiahphinearerr[1] = {0.016};
    Float_t pythiahphiaway[1] = {0.644};
    Float_t pythiahphiawayerr[1] = {0.016};
    
    Float_t pythiahhnear[1] = {0.613};
    Float_t pythiahhnearerr[1] = {0.0004};
    Float_t pythiahhaway[1] = {0.340};
    Float_t pythiahhawayerr[1] = {0.0004};

/*    
Float_t pythiaSplitnch[2] = {10.1, 23.3};
 Float_t pythiaSplitnearnch[2] = {0.005891, 0.005850};
 Float_t pythiaSplitnearerrnch[2] = {0.00010, 0.000160};
 Float_t pythiaSplitawaynch[2] = {0.007562, 0.007897};
 Float_t pythiaSplitawayerrnch[2] = {0.000107, 0.000308};
 Float_t pythiaSplittotnch[2] = {0.010536, 0.010961};
 Float_t pythiaSplittoterrnch[2] = {0.000034, 0.000032};
 Float_t pythiaSplitUEnch[2] = {0.012796, 0.012170};
 Float_t pythiaSplitUEerrnch[2] = {0.00, 0.00};

 Float_t pythiaSplithphinear[2] = {0.902241, 0.886007};
 Float_t pythiaSplithphinearerr[2] = {0.015833, 0.024278};
 Float_t pythiaSplithphiaway[2] = {0.643532, 0.623499};
 Float_t pythiaSplithphiawayerr[2] = {0.015916, 0.024299};
 
 Float_t pythiaSplithhnear[2] = {0.612656, 0.605837};
 Float_t pythiaSplithhnearerr[2] = {0.000420, 0.000632};
 Float_t pythiaSplithhaway[2] = {0.340426, 0.315834};
 Float_t pythiaSplithhawayerr[2] = {0.000413, 0.000621};
*/

    
    //low pt
    //input_0_20 = "~/phiStudies/results_onlineEff/Combined/oldONnewOFF/US_syst_trig_4.0_8.0_assoc_1.5_2.5_effcorr_hPhi_0_20_oldONnewOFF.root";
    //input_20_50 = "~/phiStudies/results_onlineEff/Combined/oldONnewOFF/US_syst_trig_4.0_8.0_assoc_1.5_2.5_effcorr_hPhi_20_50_oldONnewOFF.root";
    //input_50_80 = "~/phiStudies/results_onlineEff/Combined/oldONnewOFF/US_syst_trig_4.0_8.0_assoc_1.5_2.5_effcorr_hPhi_50_80_oldONnewOFF.root";
    //hhinput_0_20 = "~/phiStudies/results_onlineEff/Combined/oldONnewOFF/trig_4.0_8.0_assoc_1.5_2.5_effcorr_hh_0_20.root";
    //hhinput_20_50 = "~/phiStudies/results_onlineEff/Combined/oldONnewOFF/trig_4.0_8.0_assoc_1.5_2.5_effcorr_hh_20_50.root";
    //hhinput_50_80 = "~/phiStudies/results_onlineEff/Combined/oldONnewOFF/trig_4.0_8.0_assoc_1.5_2.5_effcorr_hh_50_80.root";
 
    input_0_20 = "~/phiStudies/results_onlineEff/Combined/2021code/US_syst_trig_4.0_8.0_assoc_1.5_2.5_effcorr_hPhi_0_20_hphilow.root";
    input_20_50 = "~/phiStudies/results_onlineEff/Combined/2021code/US_syst_trig_4.0_8.0_assoc_1.5_2.5_effcorr_hPhi_20_50_hphilow.root";
    input_50_80 = "~/phiStudies/results_onlineEff/Combined/2021code/US_syst_trig_4.0_8.0_assoc_1.5_2.5_effcorr_hPhi_50_80_hphilow.root";
    hhinput_0_20 = "~/phiStudies/results_onlineEff/Combined/2021code/trig_4.0_8.0_assoc_1.5_2.5_effcorr_hh_0_20.root";
    hhinput_20_50 = "~/phiStudies/results_onlineEff/Combined/2021code/trig_4.0_8.0_assoc_1.5_2.5_effcorr_hh_20_50.root";
    hhinput_50_80 = "~/phiStudies/results_onlineEff/Combined/2021code/trig_4.0_8.0_assoc_1.5_2.5_effcorr_hh_50_80.root";
    
    assocpt = "1.5 < #it{p}^{#phi}_{T,assoc} < 2.5 GeV/#it{c}";
    hhassocpt = "1.5 < #it{p}^{h}_{T,assoc} < 2.5 GeV/#it{c}";
    genassocpt = "1.5 < #it{p}_{T,assoc} < 2.5 GeV/#it{c}";

    refoffset020 = 8.604E-3*2.4*TMath::Pi()/8.0;
    refoffset2050 = 4.99E-3*2.4*TMath::Pi()/8.0;
    refoffset5080 = 2.817E-3*2.4*TMath::Pi()/8.0;
    refhhoffset020 = 4.278E-1*2.4*TMath::Pi()/8.0;
    refhhoffset2050 = 2.549E-1*2.4*TMath::Pi()/8.0;
    refhhoffset5080 = 1.432E-1*2.4*TMath::Pi()/8.0;

    refoffset020avg = 0.995;
    refoffset2050avg = 0.972;
    refoffset5080avg = 0.965;
    refhhoffset020avg = 0.998;
    refhhoffset2050avg = 0.994;
    refhhoffset5080avg = 0.988;

    refoffset020std = 0.017;
    refoffset2050std = 0.033;
    refoffset5080std = 0.038;
    refhhoffset020std = 0.002;
    refhhoffset2050std = 0.005;
    refhhoffset5080std = 0.009;
    

    //pythia8 values for 1.5-2.5 range
//    Float_t pythiaval[4] = {0.005504, 0.007018, 0.012455, 0.010798};
//    Float_t pythiazyam[4] = {0.005742, 0.007486, 0.012440, 0.010798};

    pythianch[0] = {12.3};
    pythianearnch[0] = { 0.005504};
    pythianearerrnch[0] = {0.000428};
    pythiaawaynch[0] = {0.007018};
    pythiaawayerrnch[0] = {0.000843};
    pythiatotnch[0] = {0.010798};
    pythiatoterrnch[0] = {0.000203};
    pythiaUEnch[0] = {0.012455};
    pythiaUEerrnch[0] = {0.000203};

    //jet yields for pythia
    pythiahphinear[0] = {0.872};
    pythiahphinearerr[0] = {0.014};
    pythiahphiaway[0] = {0.611};
    pythiahphiawayerr[0] = {0.014};

    pythiahhnear[0] = {0.610};
    pythiahhnearerr[0] = {0.0004};
    pythiahhaway[0] = {0.330};
    pythiahhawayerr[0] = {0.0004};

 Float_t pythiaSplitnch[2] = {10.1, 23.3};
 Float_t pythiaSplitnearnch[2] = {0.005891, 0.005850};
 Float_t pythiaSplitnearerrnch[2] = {0.00010, 0.000160};
 Float_t pythiaSplitawaynch[2] = {0.007562, 0.007897};
 Float_t pythiaSplitawayerrnch[2] = {0.000107, 0.000308};
 Float_t pythiaSplittotnch[2] = {0.010536, 0.010961};
 Float_t pythiaSplittoterrnch[2] = {0.000034, 0.000032};
 Float_t pythiaSplitUEnch[2] = {0.012796, 0.012170};
 Float_t pythiaSplitUEerrnch[2] = {0.00, 0.00};

 Float_t pythiaSplithphinear[2] = {0.902241, 0.886007};
 Float_t pythiaSplithphinearerr[2] = {0.015833, 0.024278};
 Float_t pythiaSplithphiaway[2] = {0.643532, 0.623499};
 Float_t pythiaSplithphiawayerr[2] = {0.015916, 0.024299};
 
 Float_t pythiaSplithhnear[2] = {0.612656, 0.605837};
 Float_t pythiaSplithhnearerr[2] = {0.000420, 0.000632};
 Float_t pythiaSplithhaway[2] = {0.340426, 0.315834};
 Float_t pythiaSplithhawayerr[2] = {0.000413, 0.000621};

    TString ptprefix = "lowpt";
 
/*
    //high pt
   //input_0_20 = "~/phiStudies/results_onlineEff/Combined/oldONnewOFF/US_syst_trig_4.0_8.0_assoc_2.5_4.0_effcorr_hPhi_0_20_oldONnewOFF.root";
   //input_20_50 = "~/phiStudies/results_onlineEff/Combined/oldONnewOFF/US_syst_trig_4.0_8.0_assoc_2.5_4.0_effcorr_hPhi_20_50_oldONnewOFF.root";
   //input_50_80 = "~/phiStudies/results_onlineEff/Combined/oldONnewOFF/US_syst_trig_4.0_8.0_assoc_2.5_4.0_effcorr_hPhi_50_80_oldONnewOFF.root";
   //hhinput_0_20 = "~/phiStudies/results_onlineEff/Combined/oldONnewOFF/trig_4.0_8.0_assoc_2.5_4.0_effcorr_hh_0_20.root";
   //hhinput_20_50 = "~/phiStudies/results_onlineEff/Combined/oldONnewOFF/trig_4.0_8.0_assoc_2.5_4.0_effcorr_hh_20_50.root";
   //hhinput_50_80 = "~/phiStudies/results_onlineEff/Combined/oldONnewOFF/trig_4.0_8.0_assoc_2.5_4.0_effcorr_hh_50_80.root";
 
    input_0_20 = "~/phiStudies/results_onlineEff/Combined/2021code/US_syst_trig_4.0_8.0_assoc_2.5_4.0_effcorr_hPhi_0_20_hphihigh.root";
    input_20_50 = "~/phiStudies/results_onlineEff/Combined/2021code/US_syst_trig_4.0_8.0_assoc_2.5_4.0_effcorr_hPhi_20_50_hphihigh.root";
    input_50_80 = "~/phiStudies/results_onlineEff/Combined/2021code/US_syst_trig_4.0_8.0_assoc_2.5_4.0_effcorr_hPhi_50_80_hphihigh.root";
    hhinput_0_20 = "~/phiStudies/results_onlineEff/Combined/2021code/trig_4.0_8.0_assoc_2.5_4.0_effcorr_hh_0_20.root";
    hhinput_20_50 = "~/phiStudies/results_onlineEff/Combined/2021code/trig_4.0_8.0_assoc_2.5_4.0_effcorr_hh_20_50.root";
    hhinput_50_80 = "~/phiStudies/results_onlineEff/Combined/2021code/trig_4.0_8.0_assoc_2.5_4.0_effcorr_hh_50_80.root";

    assocpt = "2.5 < #it{p}^{#phi}_{T,assoc} < 4.0 GeV/#it{c}";
    hhassocpt = "2.5 < #it{p}^{h}_{T,assoc} < 4.0 GeV/#it{c}";
    genassocpt = "2.5 < #it{p}_{T,assoc} < 4.0 GeV/#it{c}";

    refoffset020 = 3.654E-3*2.4*TMath::Pi()/8.0;
    refoffset2050 = 2.109E-3*2.4*TMath::Pi()/8.0;
    refoffset5080 = 1.137E-3*2.4*TMath::Pi()/8.0;
    refhhoffset020 = 1.082E-1*2.4*TMath::Pi()/8.0;
    refhhoffset2050 = 6.595E-2*2.4*TMath::Pi()/8.0;
    refhhoffset5080 = 3.765E-2*2.4*TMath::Pi()/8.0;

    refoffset020avg = 0.988;
    refoffset2050avg = 1.0;
    refoffset5080avg = 0.944;
    refhhoffset020avg = 0.995;
    refhhoffset2050avg = 0.99;
    refhhoffset5080avg = 0.985;

    refoffset020std = 0.008;
    refoffset2050std = 0.014;
    refoffset5080std = 0.059;
    refhhoffset020std = 0.005;
    refhhoffset2050std = 0.008;
    refhhoffset5080std = 0.013;
    //pythia8 values for 2.5-4.0 range
//    Float_t pythiaval[4] = {0.010307, 0.013598, 0.020591, 0.016886};
//    Float_t pythiazyam[4] = {0.010356, 0.013457, 0.020784, 0.016886};

    pythianch[0] = {12.3};
    pythianearnch[0] = { 0.010307};
    pythianearerrnch[0] = {0.000428};
    pythiaawaynch[0] = {0.013598};
    pythiaawayerrnch[0] = {0.000843};
    pythiatotnch[0] = {0.016886};
    pythiatoterrnch[0] = {0.000203};
    pythiaUEnch[0] = {0.020591};
    pythiaUEerrnch[0] = {0.000203};

    //jet yields for pythia
    pythiahphinear[0] = {0.894};
    pythiahphinearerr[0] = {0.010};
    pythiahphiaway[0] = {0.549};
    pythiahphiawayerr[0] = {0.010};

    pythiahhnear[0] = {0.347};
    pythiahhnearerr[0] = {0.0002};
    pythiahhaway[0] = {0.164};
    pythiahhawayerr[0] = {0.0002};

    Float_t pythiaSplitnch[2] = {10.1, 23.3};
    Float_t pythiaSplitnearnch[2] = {0.009961, 0.010741};
    Float_t pythiaSplitnearerrnch[2] = {0.000137, 0.000203};
    Float_t pythiaSplitawaynch[2] = {0.012857, 0.014095};
    Float_t pythiaSplitawayerrnch[2] = {0.000278, 0.000446};
    Float_t pythiaSplittotnch[2] = {0.016032, 0.017451};
    Float_t pythiaSplittoterrnch[2] = {0.000067, 0.000068};
    Float_t pythiaSplitUEnch[2] = {0.021556, 0.020410};
    Float_t pythiaSplitUEerrnch[2] = {0.00, 0.00};

    Float_t pythiaSplithphinear[2] = {0.855750, 0.946524};
    Float_t pythiaSplithphinearerr[2] = {0.011723, 0.017834};
    Float_t pythiaSplithphiaway[2] = {0.542866, 0.558245};
    Float_t pythiaSplithphiawayerr[2] = {0.011724, 0.017622};
    
    Float_t pythiaSplithhnear[2] = {0.343656, 0.352479};
    Float_t pythiaSplithhnearerr[2] = {0.000252, 0.000369};
    Float_t pythiaSplithhaway[2] = {0.168893, 0.158429};
    Float_t pythiaSplithhawayerr[2] = {0.000242, 0.000355};

    TString ptprefix = "highpt";
  */ 
    //done with pt range  
    
    TH1D* hhdphi_0_20 = getHisto(hhinput_0_20, "hh", "hh2D", "Eff_0_20", -1.2, 1.2, kBlue+2, 21);
    TH1D* hPhidphi_0_20 = getHisto(input_0_20, "hPhi", Form("AvgUSsubhPhi2Dpeak%s", scaleSuffix.Data()), "Eff_0_20", -1.2, 1.2, kRed+2, 22);
    
    TH1D* hPhidphi_0_20_bg = getHisto(input_0_20, "hPhi", Form("AvgUSsubhPhi2Dpeak%s", scaleSuffix.Data()), "bg_0_20", -1.2, -0.8, kRed+2, 22);
    hPhidphi_0_20_bg->Add(getHisto(input_0_20, "hPhi", Form("AvgUSsubhPhi2Dpeak%s", scaleSuffix.Data()), "bg2_0_20", 0.8, 1.2, kRed+2, 22));
    hPhidphi_0_20_bg->Scale(1.0/(hPhidphi_0_20_bg->GetBinWidth(1)*0.8));

    TH1D* hhdphi_0_20_bg = getHisto(hhinput_0_20, "hh", "hh2D", "bg_0_20", -1.2, -0.8, kBlue+2, 21);
    hhdphi_0_20_bg->Add(getHisto(hhinput_0_20, "hh", "hh2D", "bg2_0_20", 0.8, 1.2, kBlue+2, 21));
    hhdphi_0_20_bg->Scale(1.0/(hhdphi_0_20_bg->GetBinWidth(1)*0.8));


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

    TF1 *hhBGv2 = new TF1("hhBGv2", "pol0(0) + [0]*2.0*(0.11)*(0.15)*(0.5 + cos(2.0*x))", -TMath::Pi()/2.0, 3.0*TMath::Pi()/2.0);
    hhBGv2->SetParLimits(0, 0.00001, 10000000.0);
    hhBGv2->SetParameter(0, 1.0*corrFit->GetParameter(6));
    hhBGv2->SetLineStyle(2);

    TF1 *hphiBGv2 = new TF1("hphiBGv2", "pol0(0) + [0]*2.0*(0.11)*(0.13)*(0.5 + cos(2.0*x))", -TMath::Pi()/2.0, 3.0*TMath::Pi()/2.0);
    hphiBGv2->SetParLimits(0, 0.00001, 10000000.0);
    hphiBGv2->SetParameter(0, 1.0*corrFit2->GetParameter(6));
    hphiBGv2->SetLineStyle(2);

    if(IS_AVERAGE_BG){
        hhBG->SetParameter(0, refhhoffset020*refhhoffset020avg);
        hphiBG->SetParameter(0, refoffset020*refoffset020avg);
        hhBGv2->SetParameter(0, refhhoffset020*refhhoffset020avg);
        hphiBGv2->SetParameter(0, refoffset020*refoffset020avg);
    }


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

    Double_t hphibg0_20 = hphiBG->GetParameter(0);
    Double_t hhbg0_20 = hhBG->GetParameter(0);
    
    //TF1* hphibg0_20 = hphiBGv2;
    //TF1* hhbg0_20 = hhBGv2;

    //ZYAM
    if(IS_ZYAM){ 
        hphibg0_20 = ZYAM(hPhidphi_0_20);
        hhbg0_20 = ZYAM(hhdphi_0_20);
        hphiBG->SetParameter(0, hphibg0_20);
        hhBG->SetParameter(0, hhbg0_20);
        hphiBGv2->SetParameter(0, hphibg0_20);
        hhBGv2->SetParameter(0, hhbg0_20);
    } 


    Double_t near0_20hPhiYield = calcJetFromBG(hPhidphi_0_20, 1, 8, hphibg0_20, near0_20hPhiError);
    Double_t near0_20hhYield = calcJetFromBG(hhdphi_0_20, 1, 8, hhbg0_20, near0_20hhError);
    Double_t away0_20hPhiYield = calcJetFromBG(hPhidphi_0_20, 9, 16, hphibg0_20, away0_20hPhiError);
    Double_t away0_20hhYield = calcJetFromBG(hhdphi_0_20, 9, 16, hhbg0_20, away0_20hhError); 
    Double_t total0_20hPhiYield = hPhidphi_0_20->IntegralAndError(1,16,total0_20hPhiError);
    Double_t total0_20hhYield = hhdphi_0_20->IntegralAndError(1,16,total0_20hhError);
    Double_t mid0_20hPhiYield = calcUE(hPhidphi_0_20, 1, 16, hphibg0_20, mid0_20hPhiError);
    Double_t mid0_20hhYield = calcUE(hhdphi_0_20, 1, 16, hhbg0_20, mid0_20hhError); 
    Double_t jet0_20hPhi = near0_20hPhiYield + away0_20hPhiYield;
    Double_t jet0_20hPhiErr = TMath::Sqrt(TMath::Power(near0_20hPhiError, 2.0) + TMath::Power(away0_20hPhiError, 2.0));
    Double_t jet0_20hh = near0_20hhYield + away0_20hhYield;
    Double_t jet0_20hhErr = TMath::Sqrt(TMath::Power(near0_20hhError, 2.0) + TMath::Power(away0_20hhError, 2.0));

    //calc yields with v2 assumption
    Double_t trash = 0;
    Double_t near0_20hPhiYieldv2 = calcJetFromBG(hPhidphi_0_20, 1, 8, hphiBGv2, trash);
    Double_t near0_20hhYieldv2 = calcJetFromBG(hhdphi_0_20, 1, 8, hhBGv2, trash);
    Double_t away0_20hPhiYieldv2 = calcJetFromBG(hPhidphi_0_20, 9, 16, hphiBGv2, trash);
    Double_t away0_20hhYieldv2 = calcJetFromBG(hhdphi_0_20, 9, 16, hhBGv2, trash); 
    Double_t mid0_20hPhiYieldv2 = calcUE(hPhidphi_0_20, 1, 16, hphiBGv2, trash);
    Double_t mid0_20hhYieldv2 = calcUE(hhdphi_0_20, 1, 16, hhBGv2, trash); 

    Double_t near020v2 = near0_20hPhiYieldv2/near0_20hhYieldv2;
    Double_t away020v2 = away0_20hPhiYieldv2/away0_20hhYieldv2;
    
    Double_t near020 = near0_20hPhiYield/near0_20hhYield;
    Double_t near020Er = near020*TMath::Sqrt(TMath::Power(near0_20hPhiError/near0_20hPhiYield, 2) + TMath::Power(near0_20hhError/near0_20hhYield, 2));
    Double_t away020 = away0_20hPhiYield/away0_20hhYield;
    Double_t away020Er = away020*TMath::Sqrt(TMath::Power(away0_20hPhiError/away0_20hPhiYield, 2) + TMath::Power(away0_20hhError/away0_20hhYield, 2));
    Double_t mid020 = mid0_20hPhiYield/mid0_20hhYield;
    Double_t mid020Er = mid020*TMath::Sqrt(TMath::Power(mid0_20hPhiError/mid0_20hPhiYield, 2) + TMath::Power(mid0_20hhError/mid0_20hhYield, 2));
    Double_t total020 = total0_20hPhiYield/total0_20hhYield;
    Double_t total020Er = total020*TMath::Sqrt(TMath::Power(total0_20hPhiError/total0_20hPhiYield, 2) + TMath::Power(total0_20hhError/total0_20hhYield, 2));
    
    TH1D *ratios020 = new TH1D("ratios020", "(h#font[122]{}-#phi / h#font[122]{-}h) Ratios", 4, 0, 4);
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
    TH1D* hhdphi_20_50 = getHisto(hhinput_20_50, "hh", "hh2D", "Eff_20_50", -1.2, 1.2, kBlue+2, 21);
    TH1D* hPhidphi_20_50 = getHisto(input_20_50, "hPhi", Form("AvgUSsubhPhi2Dpeak%s", scaleSuffix.Data()), "Eff_20_50", -1.2, 1.2, kRed+2, 22);
    
    TH1D* hPhidphi_20_50_bg = getHisto(input_20_50, "hPhi", Form("AvgUSsubhPhi2Dpeak%s", scaleSuffix.Data()), "bg_20_50", -1.2, -0.8, kRed+2, 22);
    hPhidphi_20_50_bg->Add(getHisto(input_20_50, "hPhi", Form("AvgUSsubhPhi2Dpeak%s", scaleSuffix.Data()), "bg2_20_50", 0.8, 1.2, kRed+2, 22));
    hPhidphi_20_50_bg->Scale(1.0/(hPhidphi_20_50_bg->GetBinWidth(1)*0.8));

    TH1D* hhdphi_20_50_bg = getHisto(hhinput_20_50, "hh", "hh2D", "bg_20_50", -1.2, -0.8, kBlue+2, 21);
    hhdphi_20_50_bg->Add(getHisto(hhinput_20_50, "hh", "hh2D", "bg2_20_50", 0.8, 1.2, kBlue+2, 21));
    hhdphi_20_50_bg->Scale(1.0/(hhdphi_20_50_bg->GetBinWidth(1)*0.8));



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

    TF1 *hhBGv2_20_50 = new TF1("hhBGv2_20_50", "pol0(0) + 0.8*0.8*[0]*2.0*(0.11)*(0.13)*(0.5 + cos(2.0*x))", -TMath::Pi()/2.0, 3.0*TMath::Pi()/2.0);
    hhBGv2_20_50->SetParLimits(0, 0.00001, 10000000.0);
    hhBGv2_20_50->SetParameter(0, 1.0*corrFit2050->GetParameter(6));
    hhBGv2_20_50->SetLineStyle(2);

    TF1 *hphiBGv2_20_50 = new TF1("hphiBGv2_20_50", "pol0(0) + 0.8*0.8*[0]*2.0*(0.11)*(0.13)*(0.5 + cos(2.0*x))", -TMath::Pi()/2.0, 3.0*TMath::Pi()/2.0);
    hphiBGv2_20_50->SetParLimits(0, 0.00001, 10000000.0);
    hphiBGv2_20_50->SetParameter(0, 1.0*corrFit2_2050->GetParameter(6));
    hphiBGv2_20_50->SetLineStyle(2);

    if(IS_AVERAGE_BG){
        hhBG_20_50->SetParameter(0, refhhoffset2050*refhhoffset2050avg);
        hphiBG_20_50->SetParameter(0, refoffset2050*refoffset2050avg);
        hhBGv2_20_50->SetParameter(0, refhhoffset2050*refhhoffset2050avg);
        hphiBGv2_20_50->SetParameter(0, refoffset2050*refoffset2050avg);
    }


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

    Double_t hphibg20_50 = hphiBG_20_50->GetParameter(0);
    Double_t hhbg20_50 = hhBG_20_50->GetParameter(0);

    //TF1* hphibg20_50 = hphiBGv2_20_50;
    //TF1* hhbg20_50 = hhBGv2_20_50;

    //ZYAM
    if(IS_ZYAM){  
        hphibg20_50 = ZYAM(hPhidphi_20_50);
        hhbg20_50 = ZYAM(hhdphi_20_50);
        hphiBG_20_50->SetParameter(0, hphibg20_50);
        hhBG_20_50->SetParameter(0, hhbg20_50);
        hphiBGv2_20_50->SetParameter(0, hphibg20_50);
        hhBGv2_20_50->SetParameter(0, hhbg20_50);
    } 

    Double_t near20_50hPhiYield = calcJetFromBG(hPhidphi_20_50, 1, 8, hphibg20_50, near20_50hPhiError);
    Double_t near20_50hhYield = calcJetFromBG(hhdphi_20_50, 1, 8, hhbg20_50, near20_50hhError);
    Double_t away20_50hPhiYield = calcJetFromBG(hPhidphi_20_50, 9, 16, hphibg20_50, away20_50hPhiError);
    Double_t away20_50hhYield = calcJetFromBG(hhdphi_20_50, 9, 16, hhbg20_50, away20_50hhError); 
    Double_t mid20_50hPhiYield = calcUE(hPhidphi_20_50, 1, 16, hphibg20_50, mid20_50hPhiError);
    Double_t mid20_50hhYield = calcUE(hhdphi_20_50, 1, 16, hhbg20_50, mid20_50hhError);
    Double_t total20_50hPhiYield = hPhidphi_20_50->IntegralAndError(1, 16,total20_50hPhiError);
    Double_t total20_50hhYield = hhdphi_20_50->IntegralAndError(1, 16,total20_50hhError);
    Double_t jet20_50hPhi = near20_50hPhiYield + away20_50hPhiYield;
    Double_t jet20_50hh = near20_50hhYield + away20_50hhYield;

    //calc yields with v2 assumption
    Double_t near20_50hPhiYieldv2 = calcJetFromBG(hPhidphi_20_50, 1, 8, hphiBGv2_20_50, trash);
    Double_t near20_50hhYieldv2 = calcJetFromBG(hhdphi_20_50, 1, 8, hhBGv2_20_50, trash);
    Double_t away20_50hPhiYieldv2 = calcJetFromBG(hPhidphi_20_50, 9, 16, hphiBGv2_20_50, trash);
    Double_t away20_50hhYieldv2 = calcJetFromBG(hhdphi_20_50, 9, 16, hhBGv2_20_50, trash); 
    Double_t mid20_50hPhiYieldv2 = calcUE(hPhidphi_20_50, 1, 16, hphiBGv2_20_50, trash);
    Double_t mid20_50hhYieldv2 = calcUE(hhdphi_20_50, 1, 16, hhBGv2_20_50, trash); 

    Double_t near2050v2 = near20_50hPhiYieldv2/near20_50hhYieldv2;
    Double_t away2050v2 = away20_50hPhiYieldv2/away20_50hhYieldv2;

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

    //50-100 section
    TH1D* hhdphi_50_100 = getHisto(hhinput_50_80, "hh", "hh2D", "Eff_50_80", -1.2, 1.2, kBlue+2, 21);
    TH1D* hPhidphi_50_100 = getHisto(input_50_80, "hPhi", Form("AvgUSsubhPhi2Dpeak%s", scaleSuffix.Data()), "Eff_50_80", -1.2, 1.2, kRed+2, 22);

    TH1D* hPhidphi_50_100_bg = getHisto(input_50_80, "hPhi", Form("AvgUSsubhPhi2Dpeak%s", scaleSuffix.Data()), "bg_50_100", -1.2, -0.8, kRed+2, 22);
    hPhidphi_50_100_bg->Add(getHisto(input_50_80, "hPhi", Form("AvgUSsubhPhi2Dpeak%s", scaleSuffix.Data()), "bg2_50_100", 0.8, 1.2, kRed+2, 22));
    hPhidphi_50_100_bg->Scale(1.0/(hPhidphi_50_100_bg->GetBinWidth(1)*0.8));

    TH1D* hhdphi_50_100_bg = getHisto(hhinput_50_80, "hh", "hh2D", "bg_50_100", -1.2, -0.8, kBlue+2, 21);
    hhdphi_50_100_bg->Add(getHisto(hhinput_50_80, "hh", "hh2D", "bg2_50_100", 0.8, 1.2, kBlue+2, 21));
    hhdphi_50_100_bg->Scale(1.0/(hhdphi_50_100_bg->GetBinWidth(1)*0.8));


    //scale for inv. mass range
    //hPhidphi_50_100->Scale(1.0/0.897); //wide mass
    //hPhidphi_50_100->Scale(1.0/0.803); //narrow mass

    TF1 *corrFit50100 = setupFit("corrFit50100", hhdphi_50_100, kBlue, 7, bgMethod);

    TF1 *corrFit2_50100 = setupFit("corrFit2_50100", hPhidphi_50_100, kRed, 7, bgMethod);

    hhdphi_50_100->Fit("corrFit50100", "R0");
    hPhidphi_50_100->Fit("corrFit2_50100", "R0");

    TF1 *hhBG_50_100 = new TF1("hhBG_50_100", "pol0(0)", -TMath::Pi()/2.0, 3.0*TMath::Pi()/2.0);
    hhBG_50_100->SetParLimits(0, 0.00001, 10000000.0);
    hhBG_50_100->SetParameter(0, 1.0*corrFit50100->GetParameter(6));
    hhBG_50_100->SetLineStyle(2);

    TF1 *hphiBG_50_100 = new TF1("hphiBG_50_100", "pol0(0)", -TMath::Pi()/2.0, 3.0*TMath::Pi()/2.0);
    hphiBG_50_100->SetParLimits(0, 0.00001, 10000000.0);
    hphiBG_50_100->SetParameter(0, 1.0*corrFit2_50100->GetParameter(6));
    hphiBG_50_100->SetLineStyle(2);

    TF1 *hhBGv2_50_100 = new TF1("hhBGv2_50_100", "pol0(0) + 0.5*0.5*[0]*2.0*(0.11)*(0.13)*(0.5 + cos(2.0*x))", -TMath::Pi()/2.0, 3.0*TMath::Pi()/2.0);
    hhBGv2_50_100->SetParLimits(0, 0.00001, 10000000.0);
    hhBGv2_50_100->SetParameter(0, 1.0*corrFit50100->GetParameter(6));
    hhBGv2_50_100->SetLineStyle(2);

    TF1 *hphiBGv2_50_100 = new TF1("hphiBGv2_50_100", "pol0(0) + 0.5*0.5*[0]*2.0*(0.11)*(0.13)*(0.5 + cos(2.0*x))", -TMath::Pi()/2.0, 3.0*TMath::Pi()/2.0);
    hphiBGv2_50_100->SetParLimits(0, 0.00001, 10000000.0);
    hphiBGv2_50_100->SetParameter(0, 1.0*corrFit2_50100->GetParameter(6));
    hphiBGv2_50_100->SetLineStyle(2);


    if(IS_AVERAGE_BG){
        hhBG_50_100->SetParameter(0, refhhoffset5080*refhhoffset5080avg);
        hphiBG_50_100->SetParameter(0, refoffset5080*refoffset5080avg);
        hhBGv2_50_100->SetParameter(0, refhhoffset5080*refhhoffset5080avg);
        hphiBGv2_50_100->SetParameter(0, refoffset5080*refoffset5080avg);
    }


    Double_t near50_100hPhiError = 0;
    Double_t near50_100hhError = 0;
    Double_t away50_100hPhiError = 0;
    Double_t away50_100hhError = 0;
    //Double_t mid50_100hPhiError =(1.0/3.0)*TMath::Sqrt(TMath::Power(hPhidphi_50_100->GetBinError(8),2) + TMath::Power(hPhidphi_50_100->GetBinError(16),2) + TMath::Power(hPhidphi_50_100->GetBinError(1),2));
    Double_t mid50_100hPhiError = corrFit2_50100->GetParError(6);
    //Double_t mid50_100hhError = (1.0/3.0)*TMath::Sqrt(TMath::Power(hhdphi_50_100->GetBinError(8),2) + TMath::Power(hhdphi_50_100->GetBinError(16),2) + TMath::Power(hhdphi_50_100->GetBinError(1),2));
    Double_t mid50_100hhError = corrFit50100->GetParError(6);
    Double_t total50_100hPhiError = 0;
    Double_t total50_100hhError = 0;

    Double_t hphibg50_100 = hphiBG_50_100->GetParameter(0);
    Double_t hhbg50_100 = hhBG_50_100->GetParameter(0);

    //TF1* hphibg50_100 = hphiBGv2_50_100;
    //TF1* hhbg50_100 = hhBGv2_50_100;


    //ZYAM
    if(IS_ZYAM){ 
        hphibg50_100 = ZYAM(hPhidphi_50_100);
        hhbg50_100 = ZYAM(hhdphi_50_100);
        hphiBG_50_100->SetParameter(0, hphibg50_100);
        hhBG_50_100->SetParameter(0, hhbg50_100);
        hphiBGv2_50_100->SetParameter(0, hphibg50_100);
        hhBGv2_50_100->SetParameter(0, hhbg50_100);
    }

    Double_t near50_100hPhiYield = calcJetFromBG(hPhidphi_50_100, 1, 8, hphibg50_100, near50_100hPhiError);
    Double_t near50_100hhYield = calcJetFromBG(hhdphi_50_100, 1, 8, hhbg50_100, near50_100hhError);
    Double_t away50_100hPhiYield = calcJetFromBG(hPhidphi_50_100, 9, 16, hphibg50_100, away50_100hPhiError);
    Double_t away50_100hhYield = calcJetFromBG(hhdphi_50_100, 9, 16, hhbg50_100, away50_100hhError); 
    Double_t mid50_100hPhiYield = calcUE(hPhidphi_50_100, 1, 16, hphibg50_100, mid50_100hPhiError);
    Double_t mid50_100hhYield = calcUE(hhdphi_50_100, 1, 16, hhbg50_100, mid50_100hhError);
    Double_t total50_100hPhiYield = hPhidphi_50_100->IntegralAndError(1, 16,total50_100hPhiError);
    Double_t total50_100hhYield = hhdphi_50_100->IntegralAndError(1, 16,total50_100hhError);
    Double_t jet50_100hPhi = near50_100hPhiYield + away50_100hPhiYield;
    Double_t jet50_100hPhiErr = TMath::Sqrt(TMath::Power(near50_100hPhiError, 2) +TMath::Power(away50_100hPhiError, 2) );
    Double_t jet50_100hh = near50_100hhYield + away50_100hhYield; 
    Double_t jet50_100hhErr = TMath::Sqrt(TMath::Power(near50_100hhError, 2) +TMath::Power(away50_100hhError, 2) );

    //calc yields with v2 assumption
    Double_t near50_100hPhiYieldv2 = calcJetFromBG(hPhidphi_50_100, 1, 8, hphiBGv2_50_100, trash);
    Double_t near50_100hhYieldv2 = calcJetFromBG(hhdphi_50_100, 1, 8, hhBGv2_50_100, trash);
    Double_t away50_100hPhiYieldv2 = calcJetFromBG(hPhidphi_50_100, 9, 16, hphiBGv2_50_100, trash);
    Double_t away50_100hhYieldv2 = calcJetFromBG(hhdphi_50_100, 9, 16, hhBGv2_50_100, trash); 
    Double_t mid50_100hPhiYieldv2 = calcUE(hPhidphi_50_100, 1, 16, hphiBGv2_50_100, trash);
    Double_t mid50_100hhYieldv2 = calcUE(hhdphi_50_100, 1, 16, hhBGv2_50_100, trash); 

    Double_t near50100v2 = near50_100hPhiYieldv2/near50_100hhYieldv2;
    Double_t away50100v2 = away50_100hPhiYieldv2/away50_100hhYieldv2;


    Double_t near50100 = near50_100hPhiYield/near50_100hhYield;
    Double_t near50100Er = near50100*TMath::Sqrt(TMath::Power(near50_100hPhiError/near50_100hPhiYield, 2) + TMath::Power(near50_100hhError/near50_100hhYield, 2));
    Double_t away50100 = away50_100hPhiYield/away50_100hhYield;
    Double_t away50100Er = away50100*TMath::Sqrt(TMath::Power(away50_100hPhiError/away50_100hPhiYield, 2) + TMath::Power(away50_100hhError/away50_100hhYield, 2));
    Double_t mid50100 = mid50_100hPhiYield/mid50_100hhYield;
    Double_t mid50100Er = mid50100*TMath::Sqrt(TMath::Power(mid50_100hPhiError/mid50_100hPhiYield, 2) + TMath::Power(mid50_100hhError/mid50_100hhYield, 2));
    Double_t total50100 = total50_100hPhiYield/total50_100hhYield;
    Double_t total50100Er = total50100*TMath::Sqrt(TMath::Power(total50_100hPhiError/total50_100hPhiYield, 2) + TMath::Power(total50_100hhError/total50_100hhYield, 2));

    TH1D *ratios50100 = new TH1D("ratios50100", "(h-#phi / h-h) Ratios", 4, 0, 4);
    ratios50100->GetXaxis()->SetBinLabel(1, "near-side");
    ratios50100->SetBinContent(1, near50100);
    ratios50100->SetBinError(1, near50100Er);
    ratios50100->GetXaxis()->SetBinLabel(2, "mid");
    ratios50100->SetBinContent(2, mid50100);
    ratios50100->SetBinError(2, mid50100Er);
    ratios50100->GetXaxis()->SetBinLabel(3, "away-side");
    ratios50100->SetBinContent(3, away50100);
    ratios50100->SetBinError(3, away50100Er);
    ratios50100->GetXaxis()->SetBinLabel(4, "total");
    ratios50100->SetBinContent(4, total50100);
    ratios50100->SetBinError(4, total50100Er);
    ratios50100->SetMarkerStyle(20);
    ratios50100->SetMarkerColor(kGreen+2);
    ratios50100->SetLineColor(kGreen+2);
    ratios50100->SetMarkerSize(2);
    ratios50100->SetLineWidth(2);

    TF1* paperfit50100 = new TF2("paperfit50100", "[0]*(1+2*0.12*0.1*cos(2.0*x))", -0.5*TMath::Pi(), 1.5*TMath::Pi());
    paperfit50100->SetParameter(0, hphiBG_50_100->GetParameter(0));
        
    TF1* paperhhfit50100 = new TF2("paperhhfit50100", "[0]*(1+2*0.12*0.1*cos(2.0*x))", -0.5*TMath::Pi(), 1.5*TMath::Pi());
    paperhhfit50100->SetParameter(0, hhBG_50_100->GetParameter(0));
 

    Double_t near0_100hPhiYield = near0_20hPhiYield + near20_50hPhiYield + near50_100hPhiYield;
    Double_t near0_100hPhiError = TMath::Sqrt((TMath::Power(near0_20hPhiError,2) + TMath::Power(near20_50hPhiError,2) + TMath::Power(near50_100hPhiError,2)));
    Double_t near0_100hhYield = near0_20hhYield + near20_50hhYield + near50_100hhYield;
    Double_t near0_100hhError = TMath::Sqrt((TMath::Power(near0_20hhError,2) + TMath::Power(near20_50hhError,2) + TMath::Power(near50_100hhError,2)));
    Double_t mid0_100hPhiYield = mid0_20hPhiYield + mid20_50hPhiYield + mid50_100hPhiYield;
    Double_t mid0_100hPhiError = TMath::Sqrt((TMath::Power(mid0_20hPhiError,2) + TMath::Power(mid20_50hPhiError,2) + TMath::Power(mid50_100hPhiError,2)));
    Double_t mid0_100hhYield = mid0_20hhYield + mid20_50hhYield + mid50_100hhYield;
    Double_t mid0_100hhError = TMath::Sqrt((TMath::Power(mid0_20hhError,2) + TMath::Power(mid20_50hhError,2) + TMath::Power(mid50_100hhError,2)));
    Double_t away0_100hPhiYield = away0_20hPhiYield + away20_50hPhiYield + away50_100hPhiYield;
    Double_t away0_100hPhiError = TMath::Sqrt((TMath::Power(away0_20hPhiError,2) + TMath::Power(away20_50hPhiError,2) + TMath::Power(away50_100hPhiError,2)));
    Double_t away0_100hhYield = away0_20hhYield + away20_50hhYield + away50_100hhYield;
    Double_t away0_100hhError = TMath::Sqrt((TMath::Power(away0_20hhError,2) + TMath::Power(away20_50hhError,2) + TMath::Power(away50_100hhError,2)));
    Double_t total0_100hPhiYield = total0_20hPhiYield + total20_50hPhiYield + total50_100hPhiYield;
    Double_t total0_100hPhiError = TMath::Sqrt((TMath::Power(total0_20hPhiError,2) + TMath::Power(total20_50hPhiError,2) + TMath::Power(total50_100hPhiError,2)));
    Double_t total0_100hhYield = total0_20hhYield + total20_50hhYield + total50_100hhYield;
    Double_t total0_100hhError = TMath::Sqrt((TMath::Power(total0_20hhError,2) + TMath::Power(total20_50hhError,2) + TMath::Power(total50_100hhError,2)));
   


    Double_t near0100 = (near0_20hPhiYield + near20_50hPhiYield + near50_100hPhiYield)/(near0_20hhYield + near20_50hhYield + near50_100hhYield);
    Double_t near0100Er = near0100*TMath::Sqrt((TMath::Power(near0_20hPhiError,2) + TMath::Power(near20_50hPhiError,2) + TMath::Power(near50_100hPhiError,2))/TMath::Power((near0_20hPhiYield + near20_50hPhiYield + near50_100hPhiYield), 2) + (TMath::Power(near0_20hhError,2) + TMath::Power(near20_50hhError,2) + TMath::Power(near50_100hhError,2))/TMath::Power((near0_20hhYield + near20_50hhYield + near50_100hhYield), 2));
    Double_t away0100 = (away0_20hPhiYield + away20_50hPhiYield + away50_100hPhiYield)/(away0_20hhYield + away20_50hhYield + away50_100hhYield);
    Double_t away0100Er = away0100*TMath::Sqrt((TMath::Power(away0_20hPhiError,2) + TMath::Power(away20_50hPhiError,2) + TMath::Power(away50_100hPhiError,2))/TMath::Power((away0_20hPhiYield + away20_50hPhiYield + away50_100hPhiYield), 2) + (TMath::Power(away0_20hhError,2) + TMath::Power(away20_50hhError,2) + TMath::Power(away50_100hhError,2))/TMath::Power((away0_20hhYield + away20_50hhYield + away50_100hhYield), 2));
    Double_t mid0100 = (mid0_20hPhiYield + mid20_50hPhiYield + mid50_100hPhiYield)/(mid0_20hhYield + mid20_50hhYield + mid50_100hhYield);
    Double_t mid0100Er = mid0100*TMath::Sqrt((TMath::Power(mid0_20hPhiError,2) + TMath::Power(mid20_50hPhiError,2) + TMath::Power(mid50_100hPhiError,2))/TMath::Power((mid0_20hPhiYield + mid20_50hPhiYield + mid50_100hPhiYield), 2) + (TMath::Power(mid0_20hhError,2) + TMath::Power(mid20_50hhError,2) + TMath::Power(mid50_100hhError,2))/TMath::Power((mid0_20hhYield + mid20_50hhYield + mid50_100hhYield), 2));
    Double_t total0100 = (total0_20hPhiYield + total20_50hPhiYield + total50_100hPhiYield)/(total0_20hhYield + total20_50hhYield + total50_100hhYield);
    Double_t total0100Er = total0100*TMath::Sqrt((TMath::Power(total0_20hPhiError,2) + TMath::Power(total20_50hPhiError,2) + TMath::Power(total50_100hPhiError,2))/TMath::Power((total0_20hPhiYield + total20_50hPhiYield + total50_100hPhiYield), 2) + (TMath::Power(total0_20hhError,2) + TMath::Power(total20_50hhError,2) + TMath::Power(total50_100hhError,2))/TMath::Power((total0_20hhYield + total20_50hhYield + total50_100hhYield), 2));
    
    printf("near020Er: %E\n", near020Er);
    printf("mid020Er: %E\n", mid020Er);
    printf("total0100Er: %E\n", total0100Er);
    printf("total020Er: %E\n", total020Er);
    printf("total0_20hPhiError: %E\n", total0_20hPhiError);
    TH1D *ratios0100 = new TH1D("ratios0100", "(h-#phi / h-h) Ratios", 4, 0, 4);
    ratios0100->GetXaxis()->SetBinLabel(1, "near-side");
    ratios0100->SetBinContent(1, near0100);
    ratios0100->SetBinError(1, near0100Er);
    ratios0100->GetXaxis()->SetBinLabel(2, "mid");
    ratios0100->SetBinContent(2, mid0100);
    ratios0100->SetBinError(2, mid0100Er);
    ratios0100->GetXaxis()->SetBinLabel(3, "away-side");
    ratios0100->SetBinContent(3, away0100);
    ratios0100->SetBinError(3, away0100Er);
    ratios0100->GetXaxis()->SetBinLabel(4, "total");
    ratios0100->SetBinContent(4, total0100);
    ratios0100->SetBinError(4, total0100Er);
    ratios0100->SetMarkerStyle(29);
    ratios0100->SetMarkerColor(kViolet-1);
    ratios0100->SetLineColor(kViolet-1);
    ratios0100->SetLineWidth(2);
    ratios0100->SetMarkerSize(2);

    TLegend  *ratioslegend = new TLegend(0.3791, 0.1518, 0.8772, 0.2688);
    ratioslegend->SetMargin(0.35);
    ratioslegend->AddEntry(ratios020, "0-20%", "p");
    ratioslegend->AddEntry(ratios2050, "20-50%", "p");
    ratioslegend->AddEntry(ratios50100, "50-100%", "p");
    ratioslegend->AddEntry(ratios0100, "0-100%", "p");
    
    TLine *line = new TLine(3.0, 0.0, 3.0, 0.0040);
    line->SetLineStyle(7);
    line->SetLineWidth(2);
  
    TCanvas *testc = new TCanvas("test", "test",50, 50, 600, 600);
    testc->cd();
    ratios2050->GetYaxis()->SetRangeUser(0.0000, 0.0040);
    ratios2050->GetXaxis()->SetLabelSize(0.07);
    ratios2050->Draw("P SAME");
    ratios50100->GetXaxis()->SetLimits(0.05, 4.05);
    ratios50100->Draw("P SAME");
    ratios020->GetXaxis()->SetLimits(-0.05, 3.95);
    ratios020->Draw("P SAME");
    ratios0100->GetXaxis()->SetLimits(-0.1, 3.9);
    ratios0100->Draw("P SAME");
    line->Draw("SAME");
    ratioslegend->Draw("SAME");

    //initialize all sytematic error values
    Double_t near0_20hPhiSystError = 0.068;
    Double_t near20_50hPhiSystError = 0.077;
    Double_t near50_100hPhiSystError = 0.133;
    Double_t away0_20hPhiSystError = 0.097;
    Double_t away20_50hPhiSystError = 0.061;
    Double_t away50_100hPhiSystError = 0.111;
    Double_t mid0_20hPhiSystError = 0.017;
    Double_t mid20_50hPhiSystError = 0.013;
    Double_t mid50_100hPhiSystError = 0.037;
    Double_t total0_20hPhiSystError = 0.014;
    Double_t total20_50hPhiSystError = 0.012;
    Double_t total50_100hPhiSystError = 0.032;

    Double_t near0_20hhSystError = 0;
    Double_t near20_50hhSystError = 0;
    Double_t near50_100hhSystError = 0;
    Double_t away0_20hhSystError = 0;
    Double_t away20_50hhSystError = 0;
    Double_t away50_100hhSystError = 0;
    Double_t mid0_20hhSystError = 0;
    Double_t mid20_50hhSystError = 0;
    Double_t mid50_100hhSystError = 0;
    Double_t total0_20hhSystError = 0;
    Double_t total20_50hhSystError = 0;
    Double_t total50_100hhSystError = 0;

    //total systematic errors from all sources for yields (plotSystErrors.cxx)
    Double_t nearsysthphi[3] = {0.082, 0.09, 0.14};
    Double_t awaysysthphi[3] = {0.108, 0.076, 0.12};
    Double_t uesysthphi[3] = {0.048, 0.047, 0.059};
    Double_t totalsysthphi[3] = {0.047, 0.047, 0.056};
    Double_t hhsyst[3] = {0.04, 0.04, 0.04};


    //Setup single yield arrays for different regions
    Double_t phijetscale = 250.0;
    Double_t phitotscale = 40.0;
    Double_t nearhPhiYieldArray[3] = {near50_100hPhiYield*phijetscale, near20_50hPhiYield*phijetscale, near0_20hPhiYield*phijetscale};
    Double_t nearhPhiYieldArrayErr[3] = {near50_100hPhiError*phijetscale, near20_50hPhiError*phijetscale, near0_20hPhiError*phijetscale};
    Double_t nearhPhiYieldArraySystErr[3] = {near50_100hPhiYield*nearsysthphi[2]*phijetscale, near20_50hPhiYield*nearsysthphi[1]*phijetscale, near0_20hPhiYield*nearsysthphi[0]*phijetscale};   
    Double_t awayhPhiYieldArray[3] = {away50_100hPhiYield*phijetscale, away20_50hPhiYield*phijetscale, away0_20hPhiYield*phijetscale};
    Double_t awayhPhiYieldArrayErr[3] = {away50_100hPhiError*phijetscale, away20_50hPhiError*phijetscale, away0_20hPhiError*phijetscale};
    Double_t awayhPhiYieldArraySystErr[3] = {away50_100hPhiYield*awaysysthphi[2]*phijetscale, away20_50hPhiYield*awaysysthphi[1]*phijetscale, away0_20hPhiYield*awaysysthphi[0]*phijetscale};
    Double_t bulkhPhiYieldArray[3] = {mid50_100hPhiYield*phitotscale, mid20_50hPhiYield*phitotscale, mid0_20hPhiYield*phitotscale};
    Double_t bulkhPhiYieldArrayErr[3] = {mid50_100hPhiError*phitotscale, mid20_50hPhiError*phitotscale, mid0_20hPhiError*phitotscale};
    Double_t bulkhPhiYieldArraySystErr[3] = {mid50_100hPhiYield*uesysthphi[2]*phitotscale, mid20_50hPhiYield*uesysthphi[1]*phitotscale, mid0_20hPhiYield*uesysthphi[0]*phitotscale};
    Double_t totalhPhiYieldArray[3] = {total50_100hPhiYield*phitotscale, total20_50hPhiYield*phitotscale, total0_20hPhiYield*phitotscale};
    Double_t totalhPhiYieldArrayErr[3] = {total50_100hPhiError*phitotscale, total20_50hPhiError*phitotscale, total0_20hPhiError*phitotscale};
    Double_t totalhPhiYieldArraySystErr[3] = {total50_100hPhiYield*totalsysthphi[2]*phitotscale, total20_50hPhiYield*totalsysthphi[1]*phitotscale, total0_20hPhiYield*totalsysthphi[0]*phitotscale};

    Double_t nearhPhiYieldArrayv2[3] = {near50_100hPhiYieldv2*phijetscale, near20_50hPhiYieldv2*phijetscale, near0_20hPhiYieldv2*phijetscale};
    Double_t awayhPhiYieldArrayv2[3] = {away50_100hPhiYieldv2*phijetscale, away20_50hPhiYieldv2*phijetscale, away0_20hPhiYieldv2*phijetscale};
    Double_t bulkhPhiYieldArrayv2[3] = {mid50_100hPhiYieldv2*phitotscale, mid20_50hPhiYieldv2*phitotscale, mid0_20hPhiYieldv2*phitotscale};

    Double_t nearhhYieldArray[3] = {near50_100hhYield, near20_50hhYield, near0_20hhYield};
    Double_t nearhhYieldArrayErr[3] = {near50_100hhError, near20_50hhError, near0_20hhError};
    Double_t nearhhYieldArraySystErr[3] = {near50_100hhYield*hhsyst[2], near20_50hhYield*hhsyst[1], near0_20hhYield*hhsyst[0]};
    Double_t awayhhYieldArray[3] = {away50_100hhYield, away20_50hhYield, away0_20hhYield};
    Double_t awayhhYieldArrayErr[3] = {away50_100hhError, away20_50hhError, away0_20hhError};
    Double_t awayhhYieldArraySystErr[3] = {away50_100hhYield*hhsyst[2], away20_50hhYield*hhsyst[1], away0_20hhYield*hhsyst[0]};
    Double_t bulkhhYieldArray[3] = {mid50_100hhYield, mid20_50hhYield, mid0_20hhYield};
    Double_t bulkhhYieldArrayErr[3] = {mid50_100hhError, mid20_50hhError, mid0_20hhError};
    Double_t bulkhhYieldArraySystErr[3] = {mid50_100hhYield*hhsyst[2], mid20_50hhYield*hhsyst[1], mid0_20hhYield*hhsyst[0]};
    Double_t totalhhYieldArray[3] = {total50_100hhYield, total20_50hhYield, total0_20hhYield};
    Double_t totalhhYieldArrayErr[3] = {total50_100hhError, total20_50hhError, total0_20hhError};
    Double_t totalhhYieldArraySystErr[3] = {total50_100hhYield*hhsyst[2], total20_50hhYield*hhsyst[1], total0_20hhYield*hhsyst[0]};

    Double_t nearhhYieldArrayv2[3] = {near50_100hhYieldv2, near20_50hhYieldv2, near0_20hhYieldv2};
    Double_t awayhhYieldArrayv2[3] = {away50_100hhYieldv2, away20_50hhYieldv2, away0_20hhYieldv2};
    Double_t bulkhhYieldArrayv2[3] = {mid50_100hhYieldv2, mid20_50hhYieldv2, mid0_20hhYieldv2};

    //print out the relative change due to the v2 assumption for the jet yields
    for(int i =0; i<3; i++){
        printf("h-phi:::%d  near change: %.2f, away change: %.2f\n", i, 100.*(nearhPhiYieldArray[i] - nearhPhiYieldArrayv2[i])/(nearhPhiYieldArray[i]), 100.*(awayhPhiYieldArray[i] - awayhPhiYieldArrayv2[i])/(awayhPhiYieldArray[i]));
        printf("h-h:::%d  near change: %.2f, away change: %.2f\n", i, 100.*(nearhhYieldArray[i] - nearhhYieldArrayv2[i])/(nearhhYieldArray[i]), 100.*(awayhhYieldArray[i] - awayhhYieldArrayv2[i])/(awayhhYieldArray[i]));
    }
    //Plot ratio as a function of multiplicity for the different angular regions
    Double_t nearArray[3] = {ratios50100->GetBinContent(1), ratios2050->GetBinContent(1), ratios020->GetBinContent(1)};
    Double_t nearArrayErr[3] = {ratios50100->GetBinError(1), ratios2050->GetBinError(1), ratios020->GetBinError(1)};
    Double_t awayArray[3] = {ratios50100->GetBinContent(3), ratios2050->GetBinContent(3), ratios020->GetBinContent(3)};
    Double_t awayArrayErr[3] = {ratios50100->GetBinError(3), ratios2050->GetBinError(3), ratios020->GetBinError(3)};
    Double_t bulkArray[3] = {ratios50100->GetBinContent(2), ratios2050->GetBinContent(2), ratios020->GetBinContent(2)};
    Double_t bulkArrayErr[3] = {ratios50100->GetBinError(2), ratios2050->GetBinError(2), ratios020->GetBinError(2)};
    Double_t totalArray[3] = {ratios50100->GetBinContent(4), ratios2050->GetBinContent(4), ratios020->GetBinContent(4)};
    Double_t totalArrayErr[3] = {ratios50100->GetBinError(4), ratios2050->GetBinError(4), ratios020->GetBinError(4)};

    //ratio of Near Jet to Underlying Event vs. Multiplicity
    //Double_t jet2UEhPhi[3] = {jet50_100hPhi/mid50_100hPhiYield, jet20_50hPhi/mid20_50hPhiYield, jet0_20hPhi/mid0_20hPhiYield};
    //Double_t jet2UEhh[3] = {jet50_100hh/mid50_100hhYield, jet20_50hh/mid20_50hhYield, jet0_20hh/mid0_20hhYield}; 
    Double_t jet2UEhPhi[3] = {mid50_100hPhiYield/total50_100hPhiYield, mid20_50hPhiYield/total20_50hPhiYield, mid0_20hPhiYield/total0_20hPhiYield};
    Double_t jet2UEhh[3] = {mid50_100hhYield/total50_100hhYield, mid20_50hhYield/total20_50hhYield, mid0_20hhYield/total0_20hhYield}; 
  
    Double_t jet2tothPhi[3] = {jet50_100hPhi/total50_100hPhiYield, jet20_50hPhi/total20_50hPhiYield, jet0_20hPhi/total0_20hPhiYield};
    Double_t jet2tothh[3] = {jet50_100hh/total50_100hhYield, jet20_50hh/total20_50hhYield, jet0_20hh/total0_20hhYield};

    Double_t near2tothPhi[3] = {near50_100hPhiYield/total50_100hPhiYield, near20_50hPhiYield/total20_50hPhiYield, near0_20hPhiYield/total0_20hPhiYield};
    Double_t near2tothh[3] = {near50_100hhYield/total50_100hhYield, near20_50hhYield/total20_50hhYield, near0_20hhYield/total0_20hhYield}; 
    
    Double_t away2tothPhi[3] = {away50_100hPhiYield/total50_100hPhiYield, away20_50hPhiYield/total20_50hPhiYield, away0_20hPhiYield/total0_20hPhiYield};
    Double_t away2tothh[3] = {away50_100hhYield/total50_100hhYield, away20_50hhYield/total20_50hhYield, away0_20hhYield/total0_20hhYield}; 

    Double_t UE2tothPhi[3] = {mid50_100hPhiYield/total50_100hPhiYield, mid20_50hPhiYield/total20_50hPhiYield, mid0_20hPhiYield/total0_20hPhiYield};
    Double_t UE2tothh[3] = {mid50_100hhYield/total50_100hhYield, mid20_50hhYield/total20_50hhYield, mid0_20hhYield/total0_20hhYield}; 
   
    //systematic errors from the changing the fitting parameters
  /*  Double_t nearArraySystErr[3] = {ratios50100->GetBinContent(1)*0.068, ratios2050->GetBinContent(1)*0.077, ratios020->GetBinContent(1)*0.133};
    Double_t awayArraySystErr[3] = {ratios50100->GetBinContent(3)*0.097, ratios2050->GetBinContent(3)*0.061, ratios020->GetBinContent(3)*0.111};
    Double_t bulkArraySystErr[3] = {ratios50100->GetBinContent(2)*0.017, ratios2050->GetBinContent(2)*0.013, ratios020->GetBinContent(2)*0.037};
    Double_t totalArraySystErr[3] = {ratios50100->GetBinContent(4)*0.014, ratios2050->GetBinContent(4)*0.012, ratios020->GetBinContent(4)*0.032};
  */
    //Systematic errors for the ratio plot:
    Double_t nearsysttot[3];
    Double_t awaysysttot[3];
    Double_t uesysttot[3];
    Double_t totalsysttot[3];
    for(int i = 0; i < 3; i++){
        nearsysttot[i] = TMath::Sqrt(TMath::Power(nearsysthphi[i], 2.0) + TMath::Power(hhsyst[i], 2.0));
        awaysysttot[i] = TMath::Sqrt(TMath::Power(awaysysthphi[i], 2.0) + TMath::Power(hhsyst[i], 2.0));
        uesysttot[i] = TMath::Sqrt(TMath::Power(uesysthphi[i], 2.0) + TMath::Power(hhsyst[i], 2.0));
        totalsysttot[i] = TMath::Sqrt(TMath::Power(totalsysthphi[i], 2.0) + TMath::Power(hhsyst[i], 2.0));
    }

    Double_t nearArraySystErr[3] = {ratios50100->GetBinContent(1)*nearsysttot[2], ratios2050->GetBinContent(1)*nearsysttot[1], ratios020->GetBinContent(1)*nearsysttot[0]};
    Double_t awayArraySystErr[3] = {ratios50100->GetBinContent(3)*awaysysttot[2], ratios2050->GetBinContent(3)*awaysysttot[1], ratios020->GetBinContent(3)*awaysysttot[0]};
    Double_t bulkArraySystErr[3] = {ratios50100->GetBinContent(2)*uesysttot[2], ratios2050->GetBinContent(2)*uesysttot[1], ratios020->GetBinContent(2)*uesysttot[0]};
    Double_t totalArraySystErr[3] = {ratios50100->GetBinContent(4)*totalsysttot[2], ratios2050->GetBinContent(4)*totalsysttot[1], ratios020->GetBinContent(4)*totalsysttot[0]};

    printf("near syst: %f, %f, %f\n away syst: %f, %f, %f\n", nearsysttot[0], nearsysttot[1], nearsysttot[2], awaysysttot[0], awaysysttot[1], awaysysttot[2]);
   
    //getting the systematics from v2 uncertainty
    //TFile* v2file = new TFile("~/phiStudies/results_onlineEff/Combined/v2syst.root");
    TFile* v2file = new TFile("~/alidock/alirepos/utaustin/phiCorrelations/plotcode/ratiocode/v2syst.root");
    TGraphErrors* nearv2syst = (TGraphErrors*)v2file->Get("nearv2syst");
    nearv2syst->SetLineWidth(2);
    nearv2syst->SetLineColor(kGray+2);
    nearv2syst->SetFillColor(kGray+2);
    nearv2syst->SetFillStyle(3144);
    TGraphErrors* awayv2syst = (TGraphErrors*)v2file->Get("awayv2syst");
    awayv2syst->SetFillColor(kGray+2);
    awayv2syst->SetFillStyle(3144);

    Double_t multArray[3] = {35.0, 65.0, 90.0};
    Double_t multArrayErr[3] = {15.0, 15.0, 10.0}; //normal errors that take up whole multiplicity range
    Double_t multArraySystErr[3] = {2.5, 2.5, 2.5}; //box systematic errors that are small

    Double_t nchArray[3] = {17.7, 27.6, 42.4};
    Double_t nchArrayErr[3] = {0.4, 0.5, 0.9};
    Double_t nchArraySystErr[3] = {0.4, 0.5, 0.9};

    //checking on the fly v2 assumption errors
    TGraphErrors* nearv2fly = (TGraphErrors*)nearv2syst->Clone("nearv2fly");
    TGraphErrors* awayv2fly = (TGraphErrors*)awayv2syst->Clone("awayv2fly");

    nearv2fly->SetPoint(2, nchArray[2], 0.5*(near020v2+near020));
    nearv2fly->SetPoint(1, nchArray[1], 0.5*(near2050v2+near2050));
    nearv2fly->SetPoint(0, nchArray[0], 0.5*(near50100v2+near50100));
    nearv2fly->SetPointError(2, nchArrayErr[2], -0.5*(near020v2-near020));
    nearv2fly->SetPointError(1, nchArrayErr[1], -0.5*(near2050v2-near2050));
    nearv2fly->SetPointError(0, nchArrayErr[0], -0.5*(near50100v2-near50100));
   
    awayv2fly->SetPoint(2, nchArray[2], 0.5*(away020+away020v2));
    awayv2fly->SetPoint(1, nchArray[1], 0.5*(away2050+away2050v2));
    awayv2fly->SetPoint(0, nchArray[0], 0.5*(away50100+away50100v2));
    awayv2fly->SetPointError(2, nchArrayErr[2], -0.5*(away020-away020v2));
    awayv2fly->SetPointError(1, nchArrayErr[1], -0.5*(away2050-away2050v2));
    awayv2fly->SetPointError(0, nchArrayErr[0], -0.5*(away50100-away50100v2));

    TGraphErrors* nearv2flyMult = (TGraphErrors*)nearv2syst->Clone("nearv2flyMult");
    TGraphErrors* awayv2flyMult = (TGraphErrors*)awayv2syst->Clone("awayv2flyMult");

    nearv2flyMult->SetPoint(2, multArray[2], 0.5*(near020v2+near020));
    nearv2flyMult->SetPoint(1, multArray[1], 0.5*(near2050v2+near2050));
    nearv2flyMult->SetPoint(0, multArray[0], 0.5*(near50100v2+near50100));
    nearv2flyMult->SetPointError(2, 2, -0.5*(near020v2-near020));
    nearv2flyMult->SetPointError(1, 2, -0.5*(near2050v2-near2050));
    nearv2flyMult->SetPointError(0, 2, -0.5*(near50100v2-near50100));
   
    awayv2flyMult->SetPoint(2, multArray[2], 0.5*(away020+away020v2));
    awayv2flyMult->SetPoint(1, multArray[1], 0.5*(away2050+away2050v2));
    awayv2flyMult->SetPoint(0, multArray[0], 0.5*(away50100+away50100v2));
    awayv2flyMult->SetPointError(2, 2, -0.5*(away020-away020v2));
    awayv2flyMult->SetPointError(1, 2, -0.5*(away2050-away2050v2));
    awayv2flyMult->SetPointError(0, 2, -0.5*(away50100-away50100v2));


    // making TGraph's for the ratios with v2 for linear fitting
    Double_t nearv2Array[3] = {near50100v2, near2050v2, near020v2};
    Double_t awayv2Array[3] = {away50100v2, away2050v2, away020v2};

    TGraphErrors* nearv2graph = new TGraphErrors(3, nchArray, nearv2Array, nchArrayErr, nearArrayErr);
    TGraphErrors* awayv2graph = new TGraphErrors(3, nchArray, awayv2Array, nchArrayErr, awayArrayErr);

    TF1* nearv2fit = new TF1("nearv2fit", "[0] + [1]*x", 0, 65);
    nearv2graph->Fit(nearv2fit);
    TF1* awayv2fit = new TF1("awayv2fit", "[0] + [1]*x", 0, 65);
    awayv2graph->Fit(awayv2fit);

    
    

    printf("ratio change w v2 for most central:near:  %.2f  away:  %.2f\n", 100.*(near020-near020v2)/(near020), 100.*(away020-away020v2)/(away020));

    //Double_t multArrayErr[3] = {0.0, 0.0, 0.0}; //normal errors, no x error length
    //Double_t multArrayErr2[3] = {15.0, 15.0, 15.0}; 
    //Double_t multArraySystErr[3] = {15.0, 15.0, 10.0};
    //Double_t multArraySystErr2[3] = {2.5, 2.5, 2.5};

    Double_t mult2Array[3] = {36.0, 66.0, 91.0};
    Double_t mult2ArrayErr[3] = {15.0, 15.0, 10.0};

    //jet vs UE ratios
    TGraph* jetratioshPhi = new TGraphErrors(3, multArray, jet2tothPhi);
    jetratioshPhi->SetMarkerStyle(20);
    jetratioshPhi->SetMarkerSize(2);
    jetratioshPhi->SetMarkerColor(kCyan+1);
    jetratioshPhi->SetLineColor(kCyan+2);
    jetratioshPhi->SetLineWidth(2);
    jetratioshPhi->GetXaxis()->SetTitle("Multiplicity Percentile");
    jetratioshPhi->GetXaxis()->SetTitleSize(0.05);
    jetratioshPhi->GetXaxis()->SetLabelSize(0.05);
    jetratioshPhi->GetXaxis()->SetTitleOffset(0.9);
    jetratioshPhi->GetXaxis()->SetRangeUser(0.0, 100.0);
    jetratioshPhi->GetYaxis()->SetTitle("Yield Ratio #left(#frac{Near + Away-side}{Total Pairs}#right)");
    jetratioshPhi->GetYaxis()->SetTitleSize(0.05);
    jetratioshPhi->GetYaxis()->SetTitleOffset(1.5); 
    jetratioshPhi->GetYaxis()->SetRangeUser(0.0002, 0.0035);

    TGraph* jetratioshh = new TGraph(3, multArray, jet2tothh);
    jetratioshh->SetMarkerStyle(20);
    jetratioshh->SetMarkerSize(2);
    jetratioshh->SetMarkerColor(kOrange+1);
    jetratioshh->SetLineColor(kOrange+2);
    jetratioshh->SetLineWidth(2);
    jetratioshh->GetXaxis()->SetTitle("Multiplicity Percentile");
    jetratioshh->GetXaxis()->SetTitleSize(0.05);
    jetratioshh->GetXaxis()->SetLabelSize(0.05);
    jetratioshh->GetXaxis()->SetTitleOffset(0.9);
    jetratioshh->GetXaxis()->SetRangeUser(0.0, 100.0);
    jetratioshh->GetYaxis()->SetTitle("Yield Ratio #left(#frac{Near + Away-side}{Total Pairs}#right)");
    jetratioshh->GetYaxis()->SetTitleSize(0.05);
    jetratioshh->GetYaxis()->SetTitleOffset(1.5); 
    jetratioshh->GetYaxis()->SetRangeUser(0.0002, 0.0035);

    TGraphErrors* jetratioshPhiNch = new TGraphErrors(3, nchArray, jet2tothPhi, nchArrayErr, bulkhPhiYieldArrayErr);
    jetratioshPhiNch->SetMarkerStyle(20);
    jetratioshPhiNch->SetMarkerSize(2);
    jetratioshPhiNch->SetMarkerColor(kCyan+1);
    jetratioshPhiNch->SetLineColor(kCyan+2);
    jetratioshPhiNch->SetLineWidth(2);
    jetratioshPhiNch->GetXaxis()->SetTitle("< #it{N}_{ch} >_{|#eta|<0.8}");
    jetratioshPhiNch->GetXaxis()->SetTitleSize(0.05);
    jetratioshPhiNch->GetXaxis()->SetLabelSize(0.05);
    jetratioshPhiNch->GetXaxis()->SetTitleOffset(0.9);
    jetratioshPhiNch->GetXaxis()->SetRangeUser(0.0, 55.0);
    jetratioshPhiNch->GetYaxis()->SetTitle("Yield Ratio #left(#frac{Near + Away-side}{Total Pairs}#right)");
    jetratioshPhiNch->GetYaxis()->SetTitleSize(0.05);
    jetratioshPhiNch->GetYaxis()->SetTitleOffset(1.5); 
    jetratioshPhiNch->GetYaxis()->SetRangeUser(0.0, 0.6);

    TGraphErrors* jetratioshhNch = new TGraphErrors(3, nchArray, jet2tothh, nchArrayErr, bulkhhYieldArrayErr);
    jetratioshhNch->SetMarkerStyle(20);
    jetratioshhNch->SetMarkerSize(2);
    jetratioshhNch->SetMarkerColor(kOrange+1);
    jetratioshhNch->SetLineColor(kOrange+2);
    jetratioshhNch->SetLineWidth(2);
    jetratioshhNch->SetTitle("");
    jetratioshhNch->GetXaxis()->SetTitle("< #it{N}_{ch} >_{|#eta|<0.8}");
    jetratioshhNch->GetXaxis()->SetTitleSize(0.05);
    jetratioshhNch->GetXaxis()->SetLabelSize(0.05);
    jetratioshhNch->GetXaxis()->SetTitleOffset(0.9);
    jetratioshhNch->GetXaxis()->SetRangeUser(0.0, 55.0);
    jetratioshhNch->GetYaxis()->SetTitle("Yield Ratio #left(#frac{Near + Away-side}{Total Pairs}#right)");
    jetratioshhNch->GetYaxis()->SetTitleSize(0.04);
    jetratioshhNch->GetYaxis()->SetTitleOffset(1.4); 
    jetratioshhNch->GetYaxis()->SetRangeUser(0.0, 0.6);


    TH1D* nchratioNearHist = new TH1D("nchratioNearHist", "", 4, 0.0, 55.0);
    nchratioNearHist->GetXaxis()->SetTitle("#LT#it{N}_{ch}#GT_{|#eta| < 0.8}");
    nchratioNearHist->GetXaxis()->SetTitleSize(0.05);
    nchratioNearHist->GetXaxis()->SetLabelSize(0.05);
    nchratioNearHist->GetXaxis()->SetTitleOffset(1.2);
    nchratioNearHist->GetXaxis()->SetRangeUser(0.0, 55.0);
    nchratioNearHist->GetYaxis()->SetTitle("1/#it{N}_{trig} d#it{N}_{assoc}/d#it{#Delta#varphi} yield ratios #left(#frac{h#font[122]{-}#phi}{h#font[122]{-}h}#right)");
    nchratioNearHist->GetYaxis()->SetTitleSize(0.05);
    nchratioNearHist->GetYaxis()->SetTitleOffset(1.5); 
    nchratioNearHist->GetYaxis()->SetRangeUser(0.0, 0.045);

    //trying instead with variable sized histograms:
    Double_t binwidths[5] = {0.0, 20.0, 50.0, 80.0, 100.0};
    TH1D* ratioNearHist = new TH1D("ratioNearHist", "", 4, binwidths);
    TH1D* ratioBulkHist = new TH1D("ratioBulkHist", "", 4, binwidths);
    TH1D* ratioJetHist = new TH1D("ratioJetHist", "", 4, binwidths);
    for(int i =0; i<3; i++){
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
    ratioJetHist->GetYaxis()->SetTitle("Yield Ratio #left(#frac{Jet Yield}{Total Pairs}#right)");
    ratioJetHist->GetYaxis()->SetTitleSize(0.04);
    ratioJetHist->GetYaxis()->SetTitleOffset(1.5); 
    ratioJetHist->GetYaxis()->SetRangeUser(0.0, 0.6);

    TH1D* ratioJetHistNch =(TH1D*)ratioJetHist->Clone("ratioJetHistNch");
    ratioJetHistNch->GetXaxis()->SetTitle("#LT #it{N}_{ch}#GT_{|#eta| < 0.8}");
    ratioJetHistNch->GetXaxis()->SetRangeUser(0.0, 55.0);

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



    TGraphErrors* ratiosNear = new TGraphErrors(3, multArray, nearArray, multArrayErr, nearArrayErr);
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

    TGraphErrors* ratiosNearSyst = new TGraphErrors(3, multArray, nearArray, multArraySystErr, nearArraySystErr);
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



    TGraphErrors* ratiosAway = new TGraphErrors(3, multArray, awayArray, multArrayErr, awayArrayErr);
    ratiosAway->SetMarkerStyle(21);
    ratiosAway->SetMarkerSize(1.5);
    ratiosAway->SetMarkerColor(kBlue+1);
    ratiosAway->SetLineColor(kBlue+2);
    ratiosAway->SetLineWidth(2.0);

    TGraphErrors* ratiosAwaySyst = new TGraphErrors(3, multArray, awayArray, multArraySystErr, awayArraySystErr);
    ratiosAwaySyst->SetMarkerStyle(21);
    ratiosAwaySyst->SetMarkerSize(1);
    ratiosAwaySyst->SetMarkerColor(kBlue+1);
    ratiosAwaySyst->SetLineColor(kBlue+3);
    ratiosAwaySyst->SetLineWidth(2.0);
    ratiosAwaySyst->SetFillColorAlpha(kWhite, 0.0);

    TGraphErrors* ratiosBulk = new TGraphErrors(3, multArray, bulkArray, multArrayErr, bulkArrayErr);
    ratiosBulk->SetMarkerStyle(kFullCross);
    ratiosBulk->SetMarkerSize(2.0);
    ratiosBulk->SetMarkerColor(kGreen+2);
    ratiosBulk->SetLineColor(kGreen+3);
    ratiosBulk->SetLineWidth(2.0);

    TGraphErrors* ratiosBulkSyst = new TGraphErrors(3, multArray, bulkArray, multArraySystErr, bulkArraySystErr);
    ratiosBulkSyst->SetMarkerStyle(kFullCross);
    ratiosBulkSyst->SetMarkerSize(1);
    ratiosBulkSyst->SetMarkerColor(kGreen+2);
    ratiosBulkSyst->SetLineColor(kGreen+3);
    ratiosBulkSyst->SetLineWidth(2);
    ratiosBulkSyst->SetFillColorAlpha(kWhite, .50);
    ratiosBulkSyst->SetFillStyle(0);
   
    TGraphErrors* ratiosTot = new TGraphErrors(3, multArray, totalArray, multArrayErr, totalArrayErr);
    ratiosTot->SetMarkerStyle(33);
    ratiosTot->SetMarkerSize(3);
    ratiosTot->SetMarkerColor(kMagenta+2);
    ratiosTot->SetLineColor(kMagenta+3);
    ratiosTot->SetLineWidth(2);
    ratiosTot->SetFillColor(kMagenta+1);
    ratiosTot->SetFillStyle(3144);

    TGraphErrors* ratiosTotSyst = new TGraphErrors(3, multArray, totalArray, multArraySystErr, totalArraySystErr);
    ratiosTotSyst->SetMarkerStyle(33);
    ratiosTotSyst->SetMarkerSize(3);
    ratiosTotSyst->SetMarkerColor(kMagenta+2);
    ratiosTotSyst->SetLineColor(kMagenta+3);
    ratiosTotSyst->SetLineWidth(2);
    ratiosTotSyst->SetFillColorAlpha(kWhite, 0.0);
    //ratiosTotSyst->SetFillStyle(3144);

    //plotting ratios as function of nch instead of mult percentile
    TGraphErrors* nchratiosNear = new TGraphErrors(3, nchArray, nearArray, nchArrayErr, nearArrayErr);
    nchratiosNear->SetMarkerStyle(20);
    nchratiosNear->SetMarkerSize(1.5);
    nchratiosNear->SetMarkerColor(kRed+1);
    nchratiosNear->SetLineColor(kRed+2);
    nchratiosNear->SetLineWidth(3);
    nchratiosNear->GetXaxis()->SetTitle("Multiplicity Percentile (V0A)");
    nchratiosNear->GetXaxis()->SetTitleSize(0.05);
    nchratiosNear->GetXaxis()->SetLabelSize(0.04);
    nchratiosNear->GetXaxis()->SetTitleOffset(0.9);
    nchratiosNear->GetXaxis()->SetRangeUser(0.0, 100.0);
    nchratiosNear->GetYaxis()->SetTitle("Yield Ratio #left(#frac{h-#phi}{h-h}#right)");
    nchratiosNear->GetYaxis()->SetTitleSize(0.04);
    nchratiosNear->GetYaxis()->SetTitleOffset(1.5); 
    nchratiosNear->GetYaxis()->SetRangeUser(0.0002, 0.0035);

    TGraphErrors* nchratiosNearSyst = new TGraphErrors(3, nchArray, nearArray, nchArraySystErr, nearArraySystErr);
    nchratiosNearSyst->SetMarkerStyle(20);
    nchratiosNearSyst->SetMarkerSize(1);
    nchratiosNearSyst->SetMarkerColor(kRed+1);
    nchratiosNearSyst->SetLineColor(kRed+3);
    nchratiosNearSyst->SetFillColorAlpha(kWhite, 0.0);
    nchratiosNearSyst->SetLineWidth(2);
    nchratiosNearSyst->GetXaxis()->SetTitle("Multiplicity Percentile");
    nchratiosNearSyst->GetXaxis()->SetTitleSize(0.05);
    nchratiosNearSyst->GetXaxis()->SetLabelSize(0.04);
    nchratiosNearSyst->GetXaxis()->SetTitleOffset(0.9);
    nchratiosNearSyst->GetXaxis()->SetRangeUser(0.0, 100.0);
    nchratiosNearSyst->GetYaxis()->SetTitle("Yield Ratio #left(#frac{h-#phi}{h-h}#right)");
    nchratiosNearSyst->GetYaxis()->SetTitleSize(0.04);
    nchratiosNearSyst->GetYaxis()->SetTitleOffset(1.5); 
    nchratiosNearSyst->GetYaxis()->SetRangeUser(0.0002, 0.0035);



    TGraphErrors* nchratiosAway = new TGraphErrors(3, nchArray, awayArray, nchArrayErr, awayArrayErr);
    nchratiosAway->SetMarkerStyle(21);
    nchratiosAway->SetMarkerSize(1.5);
    nchratiosAway->SetMarkerColor(kBlue+1);
    nchratiosAway->SetLineColor(kBlue+2);
    nchratiosAway->SetLineWidth(3);

    TGraphErrors* nchratiosAwaySyst = new TGraphErrors(3, nchArray, awayArray, nchArraySystErr, awayArraySystErr);
    nchratiosAwaySyst->SetMarkerStyle(21);
    nchratiosAwaySyst->SetMarkerSize(1);
    nchratiosAwaySyst->SetMarkerColor(kBlue+1);
    nchratiosAwaySyst->SetLineColor(kBlue+3);
    nchratiosAwaySyst->SetLineWidth(2);
    nchratiosAwaySyst->SetFillColorAlpha(kWhite, 0.0);

    TGraphErrors* nchratiosBulk = new TGraphErrors(3, nchArray, bulkArray, nchArrayErr, bulkArrayErr);
    nchratiosBulk->SetMarkerStyle(kFullCross);
    nchratiosBulk->SetMarkerSize(2.0);
    nchratiosBulk->SetMarkerColor(kGreen+2);
    nchratiosBulk->SetLineColor(kGreen+3);
    nchratiosBulk->SetLineWidth(3);

    TGraphErrors* nchratiosBulkSyst = new TGraphErrors(3, nchArray, bulkArray, nchArraySystErr, bulkArraySystErr);
    nchratiosBulkSyst->SetMarkerStyle(kFullCross);
    nchratiosBulkSyst->SetMarkerSize(1);
    nchratiosBulkSyst->SetMarkerColor(kGreen+2);
    nchratiosBulkSyst->SetLineColor(kGreen+3);
    nchratiosBulkSyst->SetLineWidth(2);
    nchratiosBulkSyst->SetFillColorAlpha(kWhite, .50);
    nchratiosBulkSyst->SetFillStyle(0);
   
    TGraphErrors* nchratiosTot = new TGraphErrors(3, nchArray, totalArray, nchArrayErr, totalArrayErr);
    nchratiosTot->SetMarkerStyle(33);
    nchratiosTot->SetMarkerSize(3);
    nchratiosTot->SetMarkerColor(kMagenta+2);
    nchratiosTot->SetLineColor(kMagenta+3);
    nchratiosTot->SetLineWidth(3);
    nchratiosTot->SetFillColor(kMagenta+1);
    nchratiosTot->SetFillStyle(3144);

    TGraphErrors* nchratiosTotSyst = new TGraphErrors(3, nchArray, totalArray, nchArraySystErr, totalArraySystErr);
    nchratiosTotSyst->SetMarkerStyle(33);
    nchratiosTotSyst->SetMarkerSize(3);
    nchratiosTotSyst->SetMarkerColor(kMagenta+2);
    nchratiosTotSyst->SetLineColor(kMagenta+3);
    nchratiosTotSyst->SetLineWidth(2);
    nchratiosTotSyst->SetFillColorAlpha(kWhite, 0.0);
    //nchratiosTotSyst->SetFillStyle(3144);


    //setting up scaled TGraph's for comparing Efficiency and Non-efficiency corrected
    TH1D* ratioNearHistScaled = (TH1D*)ratioNearHist->Clone("ratioNearHistScaled");
    TGraphErrors* ratiosNearScaled = (TGraphErrors*)ratiosNear->Clone("ratiosNearScaled");
    TGraphErrors* ratiosAwayScaled = (TGraphErrors*)ratiosAway->Clone("ratiosAwayScaled");
    TGraphErrors* ratiosTotScaled = (TGraphErrors*)ratiosTot->Clone("ratiosTotScaled");
    TGraphErrors* ratiosBulkScaled = (TGraphErrors*)ratiosBulk->Clone("ratiosBulkScaled");

    ratioNearHistScaled->Scale(1.0/totalArray[0]);
    Double_t x,y;
    for(int i = 0; i<3; i++){
       ratiosNearScaled->GetPoint(i, x, y);
       ratiosNearScaled->SetPoint(i, x, y/nearArray[0]);
       ratiosNearScaled->SetPointError(i, ratiosNearScaled->GetErrorX(i), ratiosNearScaled->GetErrorY(i)/nearArray[0]);
       
       ratiosAwayScaled->GetPoint(i, x, y);
       ratiosAwayScaled->SetPoint(i, x, y/awayArray[0]);
       ratiosAwayScaled->SetPointError(i, ratiosAwayScaled->GetErrorX(i), ratiosAwayScaled->GetErrorY(i)/awayArray[0]);
       
       ratiosBulkScaled->GetPoint(i, x, y);
       ratiosBulkScaled->SetPoint(i, x, y/bulkArray[0]);
       ratiosBulkScaled->SetPointError(i, ratiosBulkScaled->GetErrorX(i), ratiosBulkScaled->GetErrorY(i)/bulkArray[0]);
       
       ratiosTotScaled->GetPoint(i, x, y);
       ratiosTotScaled->SetPoint(i, x, y/totalArray[0]);
       ratiosTotScaled->SetPointError(i, ratiosTotScaled->GetErrorX(i), ratiosTotScaled->GetErrorY(i)/totalArray[0]);
    }

    printf("away percent increase: %f.2%% +- %f\n", (ratiosAwayScaled->GetPointY(2)-1)*100, ratiosAwayScaled->GetPointY(2)*TMath::Sqrt(TMath::Power(ratiosAway->GetErrorY(0)/ratiosAway->GetPointY(0), 2)+TMath::Power(ratiosAway->GetErrorY(2)/ratiosAway->GetPointY(2), 2))*100);

    printf("near percent increase: %f.2%% +- %f\n", (ratiosNearScaled->GetPointY(2)-1)*100, ratiosNearScaled->GetPointY(2)*TMath::Sqrt(TMath::Power(ratiosNear->GetErrorY(0)/ratiosNear->GetPointY(0), 2)+TMath::Power(ratiosNear->GetErrorY(2)/ratiosNear->GetPointY(2), 2))*100);
    //straight ratio from 50-80 data
    TF1* straightratio = new TF1("straightratio", "[0]", 20.0, 50.0);
    straightratio->SetLineColor(kGray+1);
    straightratio->SetLineWidth(2);
    straightratio->SetLineStyle(7);
    straightratio->SetParameter(0, 0.0252);

    TF1* trigratio = new TF1("trigratio", "[0]", 20.0, 50.0);
    trigratio->SetLineColor(kGray+3);
    trigratio->SetLineWidth(2);
    trigratio->SetLineStyle(7);
    trigratio->SetParameter(0, 0.0188);


    //ratio from pythia 18j2 (only taking phi/h directly, no correlation measurement)
    TF1* pythiaratio = new TF1("pythiaratio", "[0]", 0, 20.0);
    pythiaratio->SetLineColor(kBlack);
    pythiaratio->SetLineWidth(2);
    pythiaratio->SetLineStyle(7);
    pythiaratio->SetParameter(0, 0.0125);
/*   
    //ratio from pythia 18j2 with full correlations measurement
    TF1* pythianearratio = new TF1("pythianearratio", "[0]", 0, 20.0);
    pythianearratio->SetLineColor(kRed+1);
    pythianearratio->SetLineWidth(2);
    pythianearratio->SetLineStyle(9);
    pythianearratio->SetParameter(0, 0.008545);
    pythianearratio->SetParError(0,0.000428); 

    TF1* pythianearratiotop = new TF1("pythianearratiotop", "[0]", 0, 20.0);
    pythianearratiotop->SetLineColor(kRed+1);
    pythianearratiotop->SetLineWidth(2);
    pythianearratiotop->SetLineStyle(7);
    pythianearratiotop->SetParameter(0, 0.008545+0.000428);

    TF1* pythianearratiobot = new TF1("pythianearratiobot", "[0]", 0, 20.0);
    pythianearratiobot->SetLineColor(kRed+1);
    pythianearratiobot->SetLineWidth(2);
    pythianearratiobot->SetLineStyle(7);
    pythianearratiobot->SetParameter(0, 0.008545-0.000428);

    TF1* pythiaawayratio = new TF1("pythiaawayratio", "[0]", 0, 20.0);
    pythiaawayratio->SetLineColor(kBlue+1);
    pythiaawayratio->SetLineWidth(2);
    pythiaawayratio->SetLineStyle(9);
    pythiaawayratio->SetParameter(0, 0.012309);
    pythiaawayratio->SetParError(0, 0.000843);

    TF1* pythiaawayratiotop = new TF1("pythiaawayratiotop", "[0]", 0, 20.0);
    pythiaawayratiotop->SetLineColor(kBlue+1);
    pythiaawayratiotop->SetLineWidth(2);
    pythiaawayratiotop->SetLineStyle(7);
    pythiaawayratiotop->SetParameter(0, 0.012309+0.000843);

    TF1* pythiaawayratiobot = new TF1("pythiaawayratiobot", "[0]", 0, 20.0);
    pythiaawayratiobot->SetLineColor(kBlue+1);
    pythiaawayratiobot->SetLineWidth(2);
    pythiaawayratiobot->SetLineStyle(7);
    pythiaawayratiobot->SetParameter(0, 0.012309-0.000843);

    TF1* pythiatotalratio = new TF1("pythiatotalratio", "[0]", 0, 20.0);
    pythiatotalratio->SetLineColor(kMagenta+2);
    pythiatotalratio->SetLineWidth(2);
    pythiatotalratio->SetLineStyle(9);
    pythiatotalratio->SetParameter(0, 0.01497);
    pythiatotalratio->SetParError(0,0.000203); 

    TF1* pythiatotalratiotop = new TF1("pythiatotalratiotop", "[0]", 0, 20.0);
    pythiatotalratiotop->SetLineColor(kMagenta+2);
    pythiatotalratiotop->SetLineWidth(2);
    pythiatotalratiotop->SetLineStyle(7);
    pythiatotalratiotop->SetParameter(0, 0.01497 + 0.000203);

    TF1* pythiatotalratiobot = new TF1("pythiatotalratiobot", "[0]", 0, 20.0);
    pythiatotalratiobot->SetLineColor(kMagenta+2);
    pythiatotalratiobot->SetLineWidth(2);
    pythiatotalratiobot->SetLineStyle(7);
    pythiatotalratiobot->SetParameter(0, 0.01497 - 0.000203);

    TF1* pythiaUEratio = new TF1("pythiaUEratio", "[0]", 0, 20.0);
    pythiaUEratio->SetLineColor(kGreen+2);
    pythiaUEratio->SetLineWidth(2);
    pythiaUEratio->SetLineStyle(9);
    pythiaUEratio->SetParameter(0, 0.02011);
*/

    //ratio from pythia 16k5a (Pythia 8) with full correlations measurement
    
    //values from 2-4 range {near, away, UE, total}
//    Float_t pythiaval[4] = {0.009705, 0.012986, 0.018427, 0.01455};
//    Float_t pythiazyam[4] = {0.010483, 0.014339, 0.017562, 0.01455};
    //values for 1.5-2.5 range
//    Float_t pythiaval[4] = {0.005452, 0.007950, 0.013260, 0.010576};
//    Float_t pythiazyam[4] = {0.006123, 0.009269, 0.012812, 0.010576};
    //values for 2.5-4.0 range
//    Float_t pythiaval[4] = {0.011205, 0.013895, 0.021607, 0.016166};
//    Float_t pythiazyam[4] = {0.012993, 0.017533, 0.018414, 0.016166};

    /*

    TF1* pythianearratio = new TF1("pythianearratio", "[0]", 0, 20.0);
    pythianearratio->SetLineColor(kRed+1);
    pythianearratio->SetLineWidth(3);
    pythianearratio->SetLineStyle(9);
    pythianearratio->SetParameter(0, pythiaval[0]);
    
    TF1* pythianearratiotop = new TF1("pythianearratiotop", "[0]", 0, 20.0);
    pythianearratiotop->SetLineColor(kRed+1);
    pythianearratiotop->SetLineWidth(2);
    pythianearratiotop->SetLineStyle(7);
    pythianearratiotop->SetParameter(0, 0.00946+0.000428);

    TF1* pythianearratiobot = new TF1("pythianearratiobot", "[0]", 0, 20.0);
    pythianearratiobot->SetLineColor(kRed+1);
    pythianearratiobot->SetLineWidth(2);
    pythianearratiobot->SetLineStyle(7);
    pythianearratiobot->SetParameter(0, pythiazyam[0]);

    TF1* pythiaawayratio = new TF1("pythiaawayratio", "[0]", 0, 20.0);
    pythiaawayratio->SetLineColor(kBlue+1);
    pythiaawayratio->SetLineWidth(3);
    pythiaawayratio->SetLineStyle(9);
    pythiaawayratio->SetParameter(0, pythiaval[1]);

    TF1* pythiaawayratiotop = new TF1("pythiaawayratiotop", "[0]", 0, 20.0);
    pythiaawayratiotop->SetLineColor(kBlue+1);
    pythiaawayratiotop->SetLineWidth(2);
    pythiaawayratiotop->SetLineStyle(7);
    pythiaawayratiotop->SetParameter(0, 0.01275+0.000843);

    TF1* pythiaawayratiobot = new TF1("pythiaawayratiobot", "[0]", 0, 20.0);
    pythiaawayratiobot->SetLineColor(kBlue+1);
    pythiaawayratiobot->SetLineWidth(2);
    pythiaawayratiobot->SetLineStyle(7);
    pythiaawayratiobot->SetParameter(0, pythiazyam[1]);

    TF1* pythiatotalratio = new TF1("pythiatotalratio", "[0]", 0, 20.0);
    pythiatotalratio->SetLineColor(kMagenta+2);
    pythiatotalratio->SetLineWidth(3);
    pythiatotalratio->SetLineStyle(9);
    pythiatotalratio->SetParameter(0, pythiaval[3]);
    pythiatotalratio->SetParError(0,0.000203); 

    TF1* pythiatotalratiotop = new TF1("pythiatotalratiotop", "[0]", 0, 20.0);
    pythiatotalratiotop->SetLineColor(kMagenta+2);
    pythiatotalratiotop->SetLineWidth(2);
    pythiatotalratiotop->SetLineStyle(7);
    pythiatotalratiotop->SetParameter(0, 0.01455 + 0.000203);

    TF1* pythiatotalratiobot = new TF1("pythiatotalratiobot", "[0]", 0, 20.0);
    pythiatotalratiobot->SetLineColor(kMagenta+2);
    pythiatotalratiobot->SetLineWidth(2);
    pythiatotalratiobot->SetLineStyle(7);
    pythiatotalratiobot->SetParameter(0, 0.01455 - 0.000203);

    TF1* pythiaUEratio = new TF1("pythiaUEratio", "[0]", 0, 20.0);
    pythiaUEratio->SetLineColor(kGreen+2);
    pythiaUEratio->SetLineWidth(3);
    pythiaUEratio->SetLineStyle(9);
    pythiaUEratio->SetParameter(0, pythiaval[2]);

    TF1* pythiaUEratiobot = new TF1("pythiaUEratiobot", "[0]", 0, 20.0);
    pythiaUEratiobot->SetLineColor(kGreen+2);
    pythiaUEratiobot->SetLineWidth(2);
    pythiaUEratiobot->SetLineStyle(7);
    pythiaUEratiobot->SetParameter(0, pythiazyam[2]);
*/

    TGraphErrors* pyNchtot = new TGraphErrors(1, pythianch, pythiatotnch, {0}, pythiatoterrnch);
    pyNchtot->SetMarkerStyle(27);
    pyNchtot->SetMarkerSize(3);
    pyNchtot->SetMarkerColor(kMagenta+2);
    pyNchtot->SetLineColor(kMagenta+3);
    pyNchtot->SetLineWidth(2);
    

    TGraphErrors* pyNchUE = new TGraphErrors(1, pythianch, pythiaUEnch, {0}, pythiaUEerrnch);
    pyNchUE->SetMarkerStyle(kOpenCross);
    pyNchUE->SetMarkerSize(2.0);
    pyNchUE->SetMarkerColor(kGreen+2);
    pyNchUE->SetLineColor(kGreen+3);
    pyNchUE->SetLineWidth(2.0);

    TGraphErrors* pyNchaway = new TGraphErrors(1, pythianch, pythiaawaynch, {0}, pythiaawayerrnch);
    pyNchaway->SetMarkerStyle(kOpenSquare);
    pyNchaway->SetMarkerSize(1.5);
    pyNchaway->SetMarkerColor(kBlue+1);
    pyNchaway->SetLineColor(kBlue+2);
    pyNchaway->SetLineWidth(2.0);

    TGraphErrors* pyNchnear = new TGraphErrors(1, pythianch, pythianearnch, {0}, pythianearerrnch);
    pyNchnear->SetMarkerStyle(kOpenCircle);
    pyNchnear->SetMarkerSize(1.5);
    pyNchnear->SetMarkerColor(kRed+1);
    pyNchnear->SetLineColor(kRed+2);
    pyNchnear->SetLineWidth(2.0);

    // pythia 0-20 and 20-80 ratios
    TGraphErrors* pySplitNchtot = new TGraphErrors(2, pythiaSplitnch, pythiaSplittotnch, {0}, pythiaSplittoterrnch);
    pySplitNchtot->SetMarkerStyle(27);
    pySplitNchtot->SetMarkerSize(3);
    pySplitNchtot->SetMarkerColor(kMagenta+2);
    pySplitNchtot->SetLineColor(kMagenta+3);
    pySplitNchtot->SetLineWidth(3);
    pySplitNchtot->SetLineStyle(7);
    

    TGraphErrors* pySplitNchUE = new TGraphErrors(2, pythiaSplitnch, pythiaSplitUEnch, {0}, pythiaSplitUEerrnch);
    pySplitNchUE->SetMarkerStyle(kOpenCross);
    pySplitNchUE->SetMarkerSize(2.0);
    pySplitNchUE->SetMarkerColor(kGreen+2);
    pySplitNchUE->SetLineColor(kGreen+3);
    pySplitNchUE->SetLineWidth(3);
    pySplitNchUE->SetLineStyle(7);

    TGraphErrors* pySplitNchaway = new TGraphErrors(2, pythiaSplitnch, pythiaSplitawaynch, {0}, pythiaSplitawayerrnch);
    pySplitNchaway->SetMarkerStyle(kOpenSquare);
    pySplitNchaway->SetMarkerSize(1.5);
    pySplitNchaway->SetMarkerColor(kBlue+1);
    pySplitNchaway->SetLineColor(kBlue+2);
    pySplitNchaway->SetLineWidth(3);
    pySplitNchaway->SetLineStyle(7);

    TGraphErrors* pySplitNchnear = new TGraphErrors(2, pythiaSplitnch, pythiaSplitnearnch, {0}, pythiaSplitnearerrnch);
    pySplitNchnear->SetMarkerStyle(kOpenCircle);
    pySplitNchnear->SetMarkerSize(1.5);
    pySplitNchnear->SetMarkerColor(kRed+1);
    pySplitNchnear->SetLineColor(kRed+2);
    pySplitNchnear->SetLineWidth(3);
    pySplitNchnear->SetLineStyle(7);

// jet yields
    TGraphErrors* pyhphiaway = new TGraphErrors(1, pythianch, pythiahphiaway, {0}, pythiahphiawayerr);
    pyhphiaway->SetMarkerStyle(kFullStar);
    pyhphiaway->SetMarkerSize(1.5);
    pyhphiaway->SetMarkerColor(kBlue+1);
    pyhphiaway->SetLineColor(kBlue+2);
    pyhphiaway->SetLineWidth(2.0);

    TGraphErrors* pyhphinear = new TGraphErrors(1, pythianch, pythiahphinear, {0}, pythiahphinearerr);
    pyhphinear->SetMarkerStyle(kFullStar);
    pyhphinear->SetMarkerSize(1.5);
    pyhphinear->SetMarkerColor(kRed+1);
    pyhphinear->SetLineColor(kRed+2);
    pyhphinear->SetLineWidth(2.0);

    TGraphErrors* pyhhaway = new TGraphErrors(1, pythianch, pythiahhaway, {0}, pythiahhawayerr);
    pyhhaway->SetMarkerStyle(kOpenStar);
    pyhhaway->SetMarkerSize(1.5);
    pyhhaway->SetMarkerColor(kBlue+1);
    pyhhaway->SetLineColor(kBlue+2);
    pyhhaway->SetLineWidth(2.0);

    TGraphErrors* pyhhnear = new TGraphErrors(1, pythianch, pythiahhnear, {0}, pythiahhnearerr);
    pyhhnear->SetMarkerStyle(kOpenStar);
    pyhhnear->SetMarkerSize(1.5);
    pyhhnear->SetMarkerColor(kRed+1);
    pyhhnear->SetLineColor(kRed+2);
    pyhhnear->SetLineWidth(2.0);


/*
    //ratio from pythia 16k5b (Pythia 6) with full correlations measurement
    TF1* pythianearratio = new TF1("pythianearratio", "[0]", 0, 20.0);
    pythianearratio->SetLineColor(kRed+1);
    pythianearratio->SetLineWidth(3);
    pythianearratio->SetLineStyle(9);
    pythianearratio->SetParameter(0, 0.00723);
    
    TF1* pythianearratiotop = new TF1("pythianearratiotop", "[0]", 0, 20.0);
    pythianearratiotop->SetLineColor(kRed+1);
    pythianearratiotop->SetLineWidth(2);
    pythianearratiotop->SetLineStyle(7);
    pythianearratiotop->SetParameter(0, 0.00723+0.000428);

    TF1* pythianearratiobot = new TF1("pythianearratiobot", "[0]", 0, 20.0);
    pythianearratiobot->SetLineColor(kRed+1);
    pythianearratiobot->SetLineWidth(2);
    pythianearratiobot->SetLineStyle(7);
    pythianearratiobot->SetParameter(0, 0.00723-0.000428);

    TF1* pythiaawayratio = new TF1("pythiaawayratio", "[0]", 0, 20.0);
    pythiaawayratio->SetLineColor(kBlue+1);
    pythiaawayratio->SetLineWidth(3);
    pythiaawayratio->SetLineStyle(9);
    pythiaawayratio->SetParameter(0, 0.0109);
    pythiaawayratio->SetParError(0, 0.000843);

    TF1* pythiaawayratiotop = new TF1("pythiaawayratiotop", "[0]", 0, 20.0);
    pythiaawayratiotop->SetLineColor(kBlue+1);
    pythiaawayratiotop->SetLineWidth(2);
    pythiaawayratiotop->SetLineStyle(7);
    pythiaawayratiotop->SetParameter(0, 0.0109+0.000843);

    TF1* pythiaawayratiobot = new TF1("pythiaawayratiobot", "[0]", 0, 20.0);
    pythiaawayratiobot->SetLineColor(kBlue+1);
    pythiaawayratiobot->SetLineWidth(2);
    pythiaawayratiobot->SetLineStyle(7);
    pythiaawayratiobot->SetParameter(0, 0.0109-0.000843);

    TF1* pythiatotalratio = new TF1("pythiatotalratio", "[0]", 0, 20.0);
    pythiatotalratio->SetLineColor(kMagenta+2);
    pythiatotalratio->SetLineWidth(3);
    pythiatotalratio->SetLineStyle(9);
    pythiatotalratio->SetParameter(0, 0.01533);
    pythiatotalratio->SetParError(0,0.000203); 

    TF1* pythiatotalratiotop = new TF1("pythiatotalratiotop", "[0]", 0, 20.0);
    pythiatotalratiotop->SetLineColor(kMagenta+2);
    pythiatotalratiotop->SetLineWidth(2);
    pythiatotalratiotop->SetLineStyle(7);
    pythiatotalratiotop->SetParameter(0, 0.01533 + 0.000203);

    TF1* pythiatotalratiobot = new TF1("pythiatotalratiobot", "[0]", 0, 20.0);
    pythiatotalratiobot->SetLineColor(kMagenta+2);
    pythiatotalratiobot->SetLineWidth(2);
    pythiatotalratiobot->SetLineStyle(7);
    pythiatotalratiobot->SetParameter(0, 0.01533 - 0.000203);

    TF1* pythiaUEratio = new TF1("pythiaUEratio", "[0]", 0, 20.0);
    pythiaUEratio->SetLineColor(kGreen+2);
    pythiaUEratio->SetLineWidth(3);
    pythiaUEratio->SetLineStyle(9);
    pythiaUEratio->SetParameter(0, 0.02013);
*/

    TLegend  *ratiosMultlegend = new TLegend(0.1873, 0.696, 0.450, 0.939);
    ratiosMultlegend->SetMargin(0.2);
    ratiosMultlegend->SetHeader("ALICE p#font[122]{-}Pb #sqrt{s_{NN}} = 5 TeV","L");
    ratiosMultlegend->AddEntry(ratiosBulk, "Underlying event", "p");
    ratiosMultlegend->AddEntry(ratiosAway, "Away-side (jet)", "p");
    ratiosMultlegend->AddEntry(ratiosNear, "Near-side (jet)", "p");
    ratiosMultlegend->AddEntry(ratiosTot, "Total (jet + UE)", "p");
    //ratiosMultlegend->AddEntry(pythiaratio, "Simul. pp (Pythia8) h-#phi/h-h", "l");
    //ratiosMultlegend->AddEntry(straightratio, "#phi/h ratio", "l");
    //ratiosMultlegend->AddEntry(trigratio, "triggered #phi/h ratio", "l");
    ratiosMultlegend->SetLineWidth(0);

    TPaveText *pythiatext = new TPaveText(0.5935, 0.6367, 0.9287, 0.7796, "NDC");
    pythiatext->AddText("Pythia #phi/h");
    pythiatext->SetTextSizePixels(18);
    pythiatext->SetFillColor(kWhite);
    pythiatext->SetBorderSize(0);
    pythiatext->SetFillStyle(0);
    pythiatext->SetTextFont(42);


    TLegend *v2leg = new TLegend(0.209, 0.605, 0.392, 0.652);
    v2leg->AddEntry(nearv2syst, "syst unc. from #it{v}_{2}", "f");
    v2leg->SetLineWidth(0);

    TLegend *ratiosUEMultlegend = new TLegend(0.183, 0.686, 0.461, 0.928);
    ratiosUEMultlegend->SetMargin(0.35);
    ratiosUEMultlegend->AddEntry(ratiosBulk, "In U.E.", "p");
    ratiosUEMultlegend->AddEntry(ratiosTot, "Total (Jet + UE)", "f");
    ratiosUEMultlegend->SetLineWidth(0);

    TLegend *ratiosJetMultlegend = new TLegend(0.183, 0.686, 0.461, 0.928);
    ratiosJetMultlegend->SetMargin(0.35);
    ratiosJetMultlegend->AddEntry(jetratioshh, "Jet/Tot for (h-h)", "p");
    ratiosJetMultlegend->AddEntry(jetratioshPhi, "Jet/Tot for (h-#phi)", "p");
    ratiosJetMultlegend->SetLineWidth(0);


    TPaveText *text2 = new TPaveText(0.545, 0.671, 0.879, 0.815, "NDC");
    text2->AddText("4.0 < #it{p}^{h}_{T,trig} < 8.0 GeV/#it{c}");
    text2->AddText(assocpt.Data());
    text2->SetTextSizePixels(28);
    text2->SetFillColor(kWhite);
    text2->SetBorderSize(0);
    text2->SetFillStyle(0);
    text2->SetTextFont(42);

    TPaveText *text2hh = new TPaveText(0.545, 0.671, 0.879, 0.815, "NDC");
    text2hh->AddText("4.0 < #it{p}^{h}_{T,trig} < 8.0 GeV/#it{c}");
    text2hh->AddText(hhassocpt.Data());
    text2hh->SetTextSizePixels(28);
    text2hh->SetFillColor(kWhite);
    text2hh->SetBorderSize(0);
    text2hh->SetFillStyle(0);
    text2hh->SetTextFont(42);

    TPaveText *text2gen = new TPaveText(0.545, 0.671, 0.879, 0.815, "NDC");
    text2gen->AddText("4.0 < #it{p}^{h}_{T,trig} < 8.0 GeV/#it{c}");
    text2gen->AddText(genassocpt.Data());
    text2gen->SetTextSizePixels(28);
    text2gen->SetFillColor(kWhite);
    text2gen->SetBorderSize(0);
    text2gen->SetFillStyle(0);
    text2gen->SetTextFont(42);


//    TPaveText *data = new TPaveText(0.6058, 0.8127, 0.9042, 0.9416, "NDC"); // old position
//    TPaveText *data = new TPaveText(0.5657, 0.756, 0.8645, 0.883, "NDC");
    TPaveText *data = new TPaveText(0.467, 0.749, 0.764, 0.876, "NDC"); 
    //data->AddText("ALICE Preliminary"); 
    data->AddText("ALICE p#font[122]{-}Pb #sqrt{#it{s}_{NN}} = 5.02 TeV");
    data->GetLine(0)->SetTextSizePixels(32);
    //data->GetLine(1)->SetTextSizePixels(24);
    data->SetBorderSize(0);
    data->SetFillColor(kWhite);
    data->SetFillStyle(0);
    data->SetTextFont(42);

    TPaveText *highmulttext = new TPaveText(0.54, 0.84, 0.82, 0.94, "NDC");
    //highmulttext->AddText("ALICE Preliminary"); 
    highmulttext->AddText("0#font[122]{-}20% multiplicity class");
    highmulttext->GetLine(0)->SetTextSizePixels(32);
    //highmulttext->GetLine(1)->SetTextSizePixels(24);
    highmulttext->SetBorderSize(0);
    highmulttext->SetFillColor(kWhite);
    highmulttext->SetFillStyle(0);
    highmulttext->SetTextFont(42);

    TPaveText *midmulttext = new TPaveText(0.52, 0.84, 0.82, 0.939, "NDC");
    //midmulttext->AddText("ALICE Preliminary"); 
    midmulttext->AddText("20#font[122]{-}50% multiplicity class");
    midmulttext->GetLine(0)->SetTextSizePixels(32);
    //midmulttext->GetLine(1)->SetTextSizePixels(24);
    midmulttext->SetBorderSize(0);
    midmulttext->SetFillColor(kWhite);
    midmulttext->SetFillStyle(0);
    midmulttext->SetTextFont(42);

    TPaveText *lowmulttext = new TPaveText(0.52, 0.84, 0.82, 0.94, "NDC");
    //lowmulttext->AddText("ALICE Preliminary"); 
    lowmulttext->AddText("50#font[122]{-}80% multiplicity class");
    lowmulttext->GetLine(0)->SetTextSizePixels(32);
    //lowmulttext->GetLine(1)->SetTextSizePixels(24);
    lowmulttext->SetBorderSize(0);
    lowmulttext->SetFillColor(kWhite);
    lowmulttext->SetFillStyle(0);
    lowmulttext->SetTextFont(42);


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
    ratioNearHist->GetYaxis()->SetDecimals(kTRUE);
    ratioNearHist->GetXaxis()->SetTitle("Multiplicity Percentile (V0A)");
    ratioNearHist->GetYaxis()->SetTitleOffset(1.1);
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
    ratiosNearSyst->Draw("5");
    //nearv2syst->Draw("2");   
    nearv2flyMult->Draw("2");
    ratiosNear->Draw("P");
    ratiosAwaySyst->Draw("5");
    //awayv2syst->Draw("2");
    awayv2flyMult->Draw("2");
    ratiosAway->Draw("P"); 
    ratiosTotSyst->Draw("5");
    ratiosBulkSyst->Draw("5");
    ratiosBulk->Draw("P");
    ratiosTot->Draw("P");
    //ratiosTot->Draw("3");
/*    
    //pythianearratiotop->Draw("SAME");
    pythianearratiobot->Draw("SAME");
    pythianearratio->Draw("SAME");
    //pythiaawayratiotop->Draw("SAME");
    pythiaawayratiobot->Draw("SAME");
    pythiaawayratio->Draw("SAME");
    //pythiatotalratiotop->Draw("SAME");
    //pythiatotalratiobot->Draw("SAME");
    pythiatotalratio->Draw("SAME");
    pythiaUEratio->Draw("SAME");
    pythiaUEratiobot->Draw("SAME");
*/    
    //straightratio->Draw("SAME");
    //trigratio->Draw("SAME");
    ratiosMultlegend->Draw();
    v2leg->Draw();
    prelim->Draw();
    //data->Draw();
    text2->Draw();
    //pythiatext->Draw();


    TF1* nearlinearfit = new TF1("nearlinearfit", "pol1", 0, 65);
    nearlinearfit->SetLineColor(kRed+1);
    nearlinearfit->SetLineWidth(3);
    TF1* nearflatfit = new TF1("nearflatfit", "pol0", 0, 65);
    TF1* awaylinearfit = new TF1("awaylinearfit", "pol1", 0, 65);
    awaylinearfit->SetLineColor(kBlue+1);
    awaylinearfit->SetLineWidth(3);
    TF1* awayflatfit = new TF1("awayflatfit", "pol0", 0, 65);
    TF1* bulklinearfit = new TF1("bulklinearfit", "pol1", 0, 65);
    bulklinearfit->SetLineColor(kGreen+2);
    bulklinearfit->SetLineWidth(3);
    TF1* bulkflatfit = new TF1("bulkflatfit", "pol0", 0, 65);
    TF1* totallinearfit = new TF1("totallinearfit", "pol1", 0, 65);
    totallinearfit->SetLineColor(kMagenta+3);
    totallinearfit->SetLineWidth(3);
    TF1* totalflatfit = new TF1("totalflatfit", "pol0", 0, 65);

    TH1D* labelppb = new TH1D("labelppb", "labelppb", 1, 0., 1.);
    labelppb->SetMarkerColor(kBlack);
    labelppb->SetMarkerStyle(kFullDiamond);
    labelppb->SetLineColor(kBlack);
    labelppb->SetLineWidth(3);
    labelppb->SetMarkerSize(3);

    TH1D* labelpy = new TH1D("labelpy", "labelpy", 1, 0., 1.);
    labelpy->SetMarkerColor(kBlack);
    labelpy->SetMarkerStyle(kOpenDiamond);
    labelpy->SetLineColor(kBlack);
    labelpy->SetLineWidth(3);
    labelpy->SetLineStyle(7);
    labelpy->SetMarkerSize(3);

    TLegend* dataleg = new TLegend(0.547, 0.642, 0.777, 0.765);
    //dataleg->AddEntry(labelppb, "p#font[122]{-}Pb, #sqrt{#it{s}_{NN}} = 5.02 TeV", "pl");
    dataleg->SetMargin(0.12);
    dataleg->AddEntry(labelpy, "PYTHIA 8 Monash pp, #sqrt{#it{s}} = 13 TeV", "p");
    dataleg->SetLineWidth(0);


    //Drawing the Nch ratios instead of mult percentiles
    TCanvas* vsNchCanvas = new TCanvas("vsNchCanvas", "vsNchCanvas", 55, 55, 900, 600);
    vsNchCanvas->cd();
    vsNchCanvas->SetMargin(0.126, 0.05, 0.125, 0.05);
    //vsMultCanvas->SetMargin(0.4, 0.05, 0.125, 0.05); //testing margin stuff
    //TH1F* hist = ratiosNear->GetHistogram();
    gStyle->SetErrorX(0.5);
    //ratioNearHist->GetYaxis()->SetRangeUser(0.0, 0.040);
    //ratioNearHist->GetXaxis()->SetTitle("Multiplicity Percentile (V0A)");
    nchratioNearHist->GetYaxis()->SetDecimals(kTRUE);
    nchratioNearHist->GetYaxis()->SetMaxDigits(2);
    nchratioNearHist->Draw("AXIS");

    //ratioNearHist->GetXaxis()->SetLabelOffset(999);
    //ratioNearHist->GetXaxis()->SetTitleOffset(999);
    //ratioNearHist->GetXaxis()->SetTickSize(0.0);

    //nchratiosNear->Draw("P");
    //gPad->Update();
    /*TGaxis* newaxis = new TGaxis(gPad->GetUxmax(),
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
    */
    
    nchratiosNearSyst->Draw("5");
    nearv2fly->Draw("2");
    //printf("fits for near side:\n");
    //nchratiosNear->Fit("nearlinearfit");   
    nchratiosNear->Draw("PZ");
    
    nchratiosAwaySyst->Draw("5");
    awayv2fly->Draw("2");
    printf("fits for away side:\n");
    //nchratiosAway->Fit("awaylinearfit");
    nchratiosAway->Draw("PZ"); 
    //printf("fits for Total: \n");
   
    nchratiosTotSyst->Draw("5");
    //nchratiosTot->Fit("totallinearfit");
    nchratiosTot->Draw("PZ");
    
    nchratiosBulkSyst->Draw("5");
    printf("fits for UE: \n");
    //nchratiosBulk->Fit("bulklinearfit");
    nchratiosBulk->Draw("PZ");

   
   // single pythia 0-100 point 
//    pyNchnear->Draw("P");
//    pyNchaway->Draw("P");
//    pyNchtot->Draw("P");
//    pyNchUE->Draw("P");

    // split pythia 0-20, 20-80
    pySplitNchnear->Draw("PL");
    pySplitNchaway->Draw("PL");
    pySplitNchtot->Draw("PL");
    pySplitNchUE->Draw("PL");
    
    //ratiosTot->Draw("3");
    //pythianearratiotop->Draw("SAME");
    //pythianearratiobot->Draw("SAME");
    //pythianearratio->Draw("SAME");
    //pythiaawayratiotop->Draw("SAME");
    //pythiaawayratiobot->Draw("SAME");
    //pythiaawayratio->Draw("SAME");
    //pythiatotalratiotop->Draw("SAME");
    //pythiatotalratiobot->Draw("SAME");
    //pythiatotalratio->Draw("SAME");
    //pythiaUEratio->Draw("SAME");
    //straightratio->Draw("SAME");
    //trigratio->Draw("SAME");
    ratiosMultlegend->Draw();
    v2leg->Draw();
    //data->Draw();
    dataleg->Draw();
    //prelim->Draw();
    text2gen->Draw();
    //pythiatext->Draw();
    //separate canvas for scaled Total ratio

    //printout info on nch fits
    printf("nearside:\n");
    printf("slope: %e, error %e, chi2: %f, ndf: %d\n", nearlinearfit->GetParameter(1), nearlinearfit->GetParError(1), nearlinearfit->GetChisquare(), nearlinearfit->GetNDF());
    printf("awayside:\n");
    printf("slope: %e, error %e, chi2: %f, ndf: %d\n", awaylinearfit->GetParameter(1), awaylinearfit->GetParError(1), awaylinearfit->GetChisquare(), awaylinearfit->GetNDF());
    printf("UE:\n");
    printf("slope: %e, error %e, chi2: %f, ndf: %d\n", bulklinearfit->GetParameter(1), bulklinearfit->GetParError(1), bulklinearfit->GetChisquare(), bulklinearfit->GetNDF());
    printf("total:\n");
    printf("slope: %e, error %e, chi2: %f, ndf: %d\n", totallinearfit->GetParameter(1), totallinearfit->GetParError(1), totallinearfit->GetChisquare(), totallinearfit->GetNDF());

    printf("with v2\n");
    printf("nearside:\n");
    printf("slope: %e, error %e, chi2: %f, ndf: %d\n", nearv2fit->GetParameter(1), nearv2fit->GetParError(1), nearv2fit->GetChisquare(), nearv2fit->GetNDF());
    printf("awayside:\n");
    printf("slope: %e, error %e, chi2: %f, ndf: %d\n", awayv2fit->GetParameter(1), awayv2fit->GetParError(1), awayv2fit->GetChisquare(), awayv2fit->GetNDF());

    Double_t x0,y0;
    Double_t x1,y1;
    TF1* nearscalefit = new TF1("nearscalefit", "[2]*(x-[0])+[1]", 0, 100);
    ratiosNearScaled->GetPoint(0, x0, y0);
    ratiosNearScaled->GetPoint(2, x1, y1);
    nearscalefit->SetParameters(x0, y0, (y1-y0)/(x1-x0));
    nearscalefit->SetLineStyle(7);
    nearscalefit->SetLineWidth(2);
    nearscalefit->SetLineColor(kRed+1);

    TF1* awayscalefit = new TF1("awayscalefit", "[2]*(x-[0])+[1]", 0, 100);
    ratiosAwayScaled->GetPoint(0, x0, y0);
    ratiosAwayScaled->GetPoint(2, x1, y1);
    awayscalefit->SetParameters(x0, y0, (y1-y0)/(x1-x0));
    awayscalefit->SetLineStyle(7);
    awayscalefit->SetLineWidth(2);
    awayscalefit->SetLineColor(kBlue+1);

    TF1* bulkscalefit = new TF1("bulkscalefit", "[2]*(x-[0])+[1]", 0, 100);
    ratiosBulkScaled->GetPoint(0, x0, y0);
    ratiosBulkScaled->GetPoint(2, x1, y1);
    bulkscalefit->SetParameters(x0, y0, (y1-y0)/(x1-x0));
    bulkscalefit->SetLineStyle(7);
    bulkscalefit->SetLineWidth(2);
    bulkscalefit->SetLineColor(kGreen+2);

    TF1* totscalefit = new TF1("totscalefit", "[2]*(x-[0])+[1]", 0, 100);
    ratiosTotScaled->GetPoint(0, x0, y0);
    ratiosTotScaled->GetPoint(2, x1, y1);
    totscalefit->SetParameters(x0, y0, (y1-y0)/(x1-x0));
    totscalefit->SetLineStyle(7);
    totscalefit->SetLineWidth(2);
    totscalefit->SetLineColor(kViolet-1);

    printf("slopes\nnear: %f, away: %f, bulk: %f, tot: %f\n", nearscalefit->GetParameter(2), awayscalefit->GetParameter(2), bulkscalefit->GetParameter(2), totscalefit->GetParameter(2));

    TCanvas* vsMultCanvasScale = new TCanvas("vsMultCanvasScale", "vsMultCanvasScale", 55, 55, 900, 600);
    vsMultCanvasScale->cd();
    vsMultCanvasScale->SetMargin(0.126, 0.05, 0.125, 0.05);
    //vsMultCanvasScale->SetMargin(0.4, 0.05, 0.125, 0.05); //testing margin stuff
    //TH1F* hist = ratiosaear->GetHistogram();
    gStyle->SetErrorX(0.5);
    //ratioNearHist->GetYaxis()->SetRangeUser(0.5, 2.0);
    ratioNearHist->GetXaxis()->SetTitle("Multiplicity Percentile (V0A)");
//    ratioNearHist->GetYaxis()->SetTitle("(h-#phi)/(h-h) / 50-80% (h-#phi)/(h-h) Ratio");
//    ratioNearHist->GetYaxis()->SetRangeUser(0.5, 2.5);
    ratioNearHist->Draw("AXIS");

    ratioNearHist->GetXaxis()->SetLabelOffset(999);
    //ratioNearHist->GetXaxis()->SetTitleOffset(999);
    ratioNearHist->GetXaxis()->SetTickSize(0.0);

    //ratiosNear->Draw("P");
    gPad->Update();
    TGaxis* newaxis2 = new TGaxis(gPad->GetUxmax(),
            gPad->GetUymin(),
            gPad->GetUxmin(),
            gPad->GetUymin(),
            ratioNearHist->GetXaxis()->GetXmin(),
            ratioNearHist->GetXaxis()->GetXmax(),
            510,"-");
    newaxis2->SetLabelOffset(-0.03);
    newaxis2->SetLabelSize(0.04);
    newaxis2->SetLabelFont(42);
    //newaxis->SetTitle("Multipliciy % (VOA)");
    //newaxis->SetTitleOffset(1.3);
    //newaxis->SetTextSize
    newaxis2->Draw();
    gPad->Update();
    //ratiosNearSyst->Draw("5");
    //nearv2syst->Draw("2");   
    ratiosNearScaled->Draw("P");
    nearscalefit->Draw("SAME");
    //ratiosAwaySyst->Draw("5");
    //awayv2syst->Draw("2");
    ratiosAwayScaled->Draw("P");
    awayscalefit->Draw("SAME");
    //ratiosTotSyst->Draw("5");
    //ratiosBulkSyst->Draw("5");
    ratiosBulkScaled->Draw("P");
    bulkscalefit->Draw("SAME");
    ratiosTotScaled->Draw("P");
    totscalefit->Draw("SAME");
    //ratiosTotScaled->Draw("3");
    //pythiaratio->Draw("SAME");
    ratiosMultlegend->Draw();
    //v2leg->Draw();
    data->Draw();
    text2->Draw();
    //pythiatext->Draw();

    //ratiosNear->Draw("PL");
    //newaxis->Draw();
    //gPad->Update();
    
    printf("near jet ratio low: %f, mid: %f, near jet ratio high: %f, %%increase: %f\n", nearArray[0],nearArray[1], nearArray[2], nearArray[2]/nearArray[0] - 1.0);
    printf("away jet ratio low: %f, mid %f, away jet ratio high: %f, %%increase: %f\n", awayArray[0], awayArray[1], awayArray[2], awayArray[2]/awayArray[0] - 1.0); 
    printf("bulk ratio low: %f, mid: %f,  bulk ratio high: %f, %%increase: %f\n", bulkArray[0], bulkArray[1], bulkArray[2], bulkArray[2]/bulkArray[0] - 1.0);
    printf("total ratio low: %f, mid:%f, total ratio high: %f, %%increase: %f\n", totalArray[0], totalArray[1], totalArray[2], totalArray[2]/totalArray[0] - 1.0);

    Double_t fakeTotHigh = near2tothh[2]*nearArray[0] + away2tothh[2]*awayArray[0] + UE2tothh[2]*bulkArray[2];
        
    printf("if jet was flat, high ratio would be: %f, total %% increase would be: %f\n", fakeTotHigh, fakeTotHigh/totalArray[0] - 1.0);
      
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


    TCanvas* JetvsNchCanvas = new TCanvas("JetvsNchCanvas", "JetvsNchCanvas", 55, 55, 900, 600);
    JetvsNchCanvas->cd();
    JetvsNchCanvas->SetMargin(0.126, 0.05, 0.125, 0.05);
    jetratioshhNch->GetXaxis()->SetLimits(0, 65);
    ratioJetHistNch->GetXaxis()->SetLimits(0, 65);
//    ratioJetHistNch->Draw("AXIS");
    jetratioshhNch->Draw("P L A");
    jetratioshPhiNch->Draw("P L");
    jetratioshhNch->GetXaxis()->SetRangeUser(0, 65);
    jetratioshPhiNch->GetXaxis()->SetRangeUser(0, 65);
    data->Draw();
    text2hh->Draw("SAME");
    ratiosJetMultlegend->Draw("SAME");

    

 
    //scaled ratios
    TCanvas* vsMultCanvasScaled = new TCanvas("vsMultCanvasScaled", "vsMultCanvasScaled", 55, 55, 900, 600);
    vsMultCanvasScaled->cd();
    vsMultCanvasScaled->SetMargin(0.126, 0.05, 0.125, 0.05);
    //TH1F* hist = ratiosNear->GetHistogram();
    gStyle->SetErrorX(0.5);
    ratioNearHistScaled->Draw("PE");

    ratioNearHistScaled->GetXaxis()->SetLabelOffset(999);
    //ratioNearHist->GetXaxis()->SetTitleOffset(999);
    ratioNearHistScaled->GetXaxis()->SetTickSize(0.0);

    //ratiosNear->Draw("P");
    gPad->Update();
    TGaxis *newaxisScaled = new TGaxis(gPad->GetUxmax(),
            gPad->GetUymin(),
            gPad->GetUxmin(),
            gPad->GetUymin(),
            ratioNearHistScaled->GetXaxis()->GetXmin(),
            ratioNearHistScaled->GetXaxis()->GetXmax(),
            510,"-");
    newaxisScaled->SetLabelOffset(-0.03);
    //newaxis->SetTitle("Multipliciy % (VOA)");
    //newaxis->SetTitleOffset(1.3);
    newaxisScaled->Draw();   
    //ratiosNearSyst->Draw("[]");
    ratiosNearScaled->Draw("P");
    //ratiosAwaySyst->Draw("[]");
    ratiosAwayScaled->Draw("P");
    ratiosBulkScaled->Draw("P");
    //ratiosJet->Draw("P");
    ratiosTotScaled->Draw("2");
    //ratiosTot->Draw("3");
    ratiosMultlegend->Draw();
    data->Draw();
    //ratiosNear->Draw("PL");
    //newaxis->Draw();
    //gPad->Update();
  
    //Just Underlying Event 
    TCanvas* vsUEMultCanvas = new TCanvas("vsUEMultCanvas", "vsUEMultCanvas", 55, 55, 900, 600);
    vsUEMultCanvas->cd();
    vsUEMultCanvas->SetMargin(0.126, 0.05, 0.125, 0.05);
    //TH1F* hist = ratiosNear->GetHistogram();
    gStyle->SetErrorX(0.5);
    ratioBulkHist->Draw("PE");

    ratioBulkHist->GetXaxis()->SetLabelOffset(999);
    ratioBulkHist->GetXaxis()->SetTitleOffset(999);
    ratioBulkHist->GetXaxis()->SetTickSize(0.0);

    gPad->Update(); 
    newaxis->Draw();   
    //ratiosBulk->Draw("P");
    ratiosTot->Draw("2");
    ratiosUEMultlegend->Draw();
    data->Draw();


    //yield graphs
    TGraphErrors* yieldsNear = new TGraphErrors(3, nchArray, nearhPhiYieldArray, nchArrayErr, nearhPhiYieldArrayErr);
    yieldsNear->SetMarkerStyle(20);
    yieldsNear->SetMarkerSize(2);
    yieldsNear->SetMarkerColor(kRed+1);
    yieldsNear->SetLineColor(kRed+2);
    yieldsNear->SetLineWidth(2);
    yieldsNear->GetXaxis()->SetTitle("Multiplicity Percentile");
    yieldsNear->GetXaxis()->SetTitleSize(0.05);
    yieldsNear->GetXaxis()->SetLabelSize(0.04);
    yieldsNear->GetXaxis()->SetTitleOffset(0.9);
    yieldsNear->GetXaxis()->SetRangeUser(0.0, 100.0);
    yieldsNear->GetYaxis()->SetTitle("Yield Ratio #left(#frac{h-#phi}{h-h}#right)");
    yieldsNear->GetYaxis()->SetTitleSize(0.04);
    yieldsNear->GetYaxis()->SetTitleOffset(1.5); 
    yieldsNear->GetYaxis()->SetRangeUser(0.0002, 0.0035);

    TGraphErrors* yieldsNearSyst = new TGraphErrors(3, nchArray, nearhPhiYieldArray, nchArraySystErr, nearhPhiYieldArraySystErr);
    yieldsNearSyst->SetMarkerStyle(20);
    yieldsNearSyst->SetMarkerSize(1);
    yieldsNearSyst->SetMarkerColor(kRed+1);
    yieldsNearSyst->SetLineColor(kRed+3);
    //yieldsNearSyst->SetFillColor(kWhite);
    yieldsNearSyst->SetLineWidth(2);
    yieldsNearSyst->GetXaxis()->SetTitle("Multiplicity Percentile");
    yieldsNearSyst->GetXaxis()->SetTitleSize(0.05);
    yieldsNearSyst->GetXaxis()->SetLabelSize(0.04);
    yieldsNearSyst->GetXaxis()->SetTitleOffset(0.9);
    yieldsNearSyst->GetXaxis()->SetRangeUser(0.0, 100.0);
    yieldsNearSyst->GetYaxis()->SetTitle("Yield Ratio #left(#frac{h-#phi}{h-h}#right)");
    yieldsNearSyst->GetYaxis()->SetTitleSize(0.04);
    yieldsNearSyst->GetYaxis()->SetTitleOffset(1.5); 
    yieldsNearSyst->GetYaxis()->SetRangeUser(0.0002, 0.0035);

    TGraphErrors* yieldsNearv2 = new TGraphErrors(3, nchArray, nearhPhiYieldArrayv2, nchArrayErr, nearhPhiYieldArrayErr);
    yieldsNearv2->SetMarkerStyle(20);
    yieldsNearv2->SetMarkerSize(2);
    yieldsNearv2->SetMarkerColor(kRed+1);
    yieldsNearv2->SetLineColor(kRed+2);
    yieldsNearv2->SetLineWidth(2);
    yieldsNearv2->GetXaxis()->SetTitle("Multiplicity Percentile");
    yieldsNearv2->GetXaxis()->SetTitleSize(0.05);
    yieldsNearv2->GetXaxis()->SetLabelSize(0.04);
    yieldsNearv2->GetXaxis()->SetTitleOffset(0.9);
    yieldsNearv2->GetXaxis()->SetRangeUser(0.0, 100.0);
    yieldsNearv2->GetYaxis()->SetTitle("Yield Ratio #left(#frac{h-#phi}{h-h}#right)");
    yieldsNearv2->GetYaxis()->SetTitleSize(0.04);
    yieldsNearv2->GetYaxis()->SetTitleOffset(1.5); 
    yieldsNearv2->GetYaxis()->SetRangeUser(0.0002, 0.0035);

    TGraphErrors* yieldsNearSystv2 = (TGraphErrors*)yieldsNearSyst->Clone("yieldsNearSystv2");
    yieldsNearSystv2->SetPoint(0, nchArray[0], 0.5*(nearhPhiYieldArray[0]+nearhPhiYieldArrayv2[0]));
    yieldsNearSystv2->SetPoint(1, nchArray[1], 0.5*(nearhPhiYieldArray[1]+nearhPhiYieldArrayv2[1]));
    yieldsNearSystv2->SetPoint(2, nchArray[2], 0.5*(nearhPhiYieldArray[2]+nearhPhiYieldArrayv2[2]));
    yieldsNearSystv2->SetPointError(0, 2.0, -0.5*(nearhPhiYieldArray[0]-nearhPhiYieldArrayv2[0]));
    yieldsNearSystv2->SetPointError(1, 2.0, -0.5*(nearhPhiYieldArray[1]-nearhPhiYieldArrayv2[1]));
    yieldsNearSystv2->SetPointError(2, 2.0, -0.5*(nearhPhiYieldArray[2]-nearhPhiYieldArrayv2[2]));
    yieldsNearSystv2->SetFillColor(kRed);
    yieldsNearSystv2->SetFillStyle(3245);
//    yieldsNearSystv2->SetLineWidth(0);

    TGraphErrors* yieldsAway = new TGraphErrors(3, nchArray, awayhPhiYieldArray, nchArrayErr, awayhPhiYieldArrayErr);
    yieldsAway->SetMarkerStyle(21);
    yieldsAway->SetMarkerSize(2);
    yieldsAway->SetMarkerColor(kBlue+1);
    yieldsAway->SetLineColor(kBlue+2);
    yieldsAway->SetLineWidth(2);

    TGraphErrors* yieldsAwaySyst = new TGraphErrors(3, nchArray, awayhPhiYieldArray, nchArraySystErr, awayhPhiYieldArraySystErr);
    yieldsAwaySyst->SetMarkerStyle(21);
    yieldsAwaySyst->SetMarkerSize(1);
    yieldsAwaySyst->SetMarkerColor(kBlue+1);
    yieldsAwaySyst->SetLineColor(kBlue+3);
    yieldsAwaySyst->SetLineWidth(2);
    //yieldsAwaySyst->SetFillColor(kWhite);
    
    TGraphErrors* yieldsAwayv2 = new TGraphErrors(3, nchArray, awayhPhiYieldArrayv2, nchArrayErr, awayhPhiYieldArrayErr);
    yieldsAwayv2->SetMarkerStyle(21);
    yieldsAwayv2->SetMarkerSize(2);
    yieldsAwayv2->SetMarkerColor(kBlue+1);
    yieldsAwayv2->SetLineColor(kBlue+2);
    yieldsAwayv2->SetLineWidth(2);


    TGraphErrors* yieldsAwaySystv2 = (TGraphErrors*)yieldsAwaySyst->Clone("yieldsAwaySystv2");
    yieldsAwaySystv2->SetPoint(0, nchArray[0], 0.5*(awayhPhiYieldArray[0]+awayhPhiYieldArrayv2[0]));
    yieldsAwaySystv2->SetPoint(1, nchArray[1], 0.5*(awayhPhiYieldArray[1]+awayhPhiYieldArrayv2[1]));
    yieldsAwaySystv2->SetPoint(2, nchArray[2], 0.5*(awayhPhiYieldArray[2]+awayhPhiYieldArrayv2[2]));
    yieldsAwaySystv2->SetPointError(0, 2.0, -0.5*(awayhPhiYieldArray[0]-awayhPhiYieldArrayv2[0]));
    yieldsAwaySystv2->SetPointError(1, 2.0, -0.5*(awayhPhiYieldArray[1]-awayhPhiYieldArrayv2[1]));
    yieldsAwaySystv2->SetPointError(2, 2.0, -0.5*(awayhPhiYieldArray[2]-awayhPhiYieldArrayv2[2]));
    yieldsAwaySystv2->SetFillColor(kBlue);
    yieldsAwaySystv2->SetFillStyle(3254);
//    yieldsAwaySystv2->SetLineWidth(0);


    TGraphErrors* yieldsBulk = new TGraphErrors(3, nchArray, bulkhPhiYieldArray, nchArrayErr, bulkhPhiYieldArrayErr);
    yieldsBulk->SetMarkerStyle(22);
    yieldsBulk->SetMarkerSize(3);
    yieldsBulk->SetMarkerColor(kGreen+2);
    yieldsBulk->SetLineColor(kGreen+3);
    yieldsBulk->SetLineWidth(4);

    TGraphErrors* yieldsBulkSyst = new TGraphErrors(3, nchArray, bulkhPhiYieldArray, nchArraySystErr, bulkhPhiYieldArraySystErr);
    yieldsBulkSyst->SetMarkerStyle(22);
    yieldsBulkSyst->SetMarkerSize(1);
    yieldsBulkSyst->SetMarkerColor(kGreen+2);
    yieldsBulkSyst->SetLineColor(kGreen+3);
    yieldsBulkSyst->SetLineWidth(2);
    yieldsBulkSyst->SetFillColor(kWhite);

   
    TGraphErrors* yieldsTot = new TGraphErrors(3, nchArray, totalhPhiYieldArray, nchArrayErr, totalhPhiYieldArrayErr);
    yieldsTot->SetMarkerStyle(29);
    yieldsTot->SetMarkerSize(3);
    yieldsTot->SetMarkerColor(kMagenta+2);
    yieldsTot->SetLineColor(kMagenta+3);
    yieldsTot->SetLineWidth(2);
    yieldsTot->SetFillColor(kMagenta+1);
    yieldsTot->SetFillStyle(3144);

    TGraphErrors* yieldsTotSyst = new TGraphErrors(3, nchArray, totalhPhiYieldArray, nchArrayErr, totalhPhiYieldArraySystErr);
    yieldsTotSyst->SetMarkerStyle(29);
    yieldsTotSyst->SetMarkerSize(1);
    yieldsTotSyst->SetMarkerColor(kMagenta+2);
    yieldsTotSyst->SetLineColor(kMagenta+3);
    yieldsTotSyst->SetLineWidth(2);
    //yieldsTot->SetFillColor(kMagenta+1);
    //yieldsTot->SetFillStyle(3144);


    TGraphErrors* yieldshhNear = new TGraphErrors(3, nchArray, nearhhYieldArray, nchArrayErr, nearhhYieldArrayErr);
    yieldshhNear->SetMarkerStyle(24);
    yieldshhNear->SetMarkerSize(2);
    yieldshhNear->SetMarkerColor(kRed+2);
    yieldshhNear->SetLineColor(kRed+3);
    yieldshhNear->SetLineWidth(2);
    yieldshhNear->GetXaxis()->SetTitle("Multiplicity Percentile");
    yieldshhNear->GetXaxis()->SetTitleSize(0.05);
    yieldshhNear->GetXaxis()->SetLabelSize(0.04);
    yieldshhNear->GetXaxis()->SetTitleOffset(0.9);
    yieldshhNear->GetXaxis()->SetRangeUser(0.0, 100.0);
    yieldshhNear->GetYaxis()->SetTitle("Yield Ratio #left(#frac{h-#phi}{h-h}#right)");
    yieldshhNear->GetYaxis()->SetTitleSize(0.04);
    yieldshhNear->GetYaxis()->SetTitleOffset(1.5); 
    yieldshhNear->GetYaxis()->SetRangeUser(0.0002, 0.0035);

    TGraphErrors* yieldshhNearSyst = new TGraphErrors(3, nchArray, nearhhYieldArray, nchArraySystErr, nearhhYieldArraySystErr);
    yieldshhNearSyst->SetMarkerStyle(24);
    yieldshhNearSyst->SetMarkerSize(2);
    yieldshhNearSyst->SetMarkerColor(kRed+2);
    yieldshhNearSyst->SetLineColor(kRed+3);
//    yieldshhNearSyst->SetFillColor(kWhite);
    yieldshhNearSyst->SetLineWidth(2);
    yieldshhNearSyst->GetXaxis()->SetTitle("Multiplicity Percentile");
    yieldshhNearSyst->GetXaxis()->SetTitleSize(0.05);
    yieldshhNearSyst->GetXaxis()->SetLabelSize(0.04);
    yieldshhNearSyst->GetXaxis()->SetTitleOffset(0.9);
    yieldshhNearSyst->GetXaxis()->SetRangeUser(0.0, 100.0);
    yieldshhNearSyst->GetYaxis()->SetTitle("Yield Ratio #left(#frac{h-#phi}{h-h}#right)");
    yieldshhNearSyst->GetYaxis()->SetTitleSize(0.04);
    yieldshhNearSyst->GetYaxis()->SetTitleOffset(1.5); 
    yieldshhNearSyst->GetYaxis()->SetRangeUser(0.0002, 0.0035);

    TGraphErrors* yieldshhNearv2 = new TGraphErrors(3, nchArray, nearhhYieldArrayv2, nchArrayErr, nearhhYieldArrayErr);
    yieldshhNearv2->SetMarkerStyle(25);
    yieldshhNearv2->SetMarkerSize(2);
    yieldshhNearv2->SetMarkerColor(kBlue+2);
    yieldshhNearv2->SetLineColor(kBlue+3);
    yieldshhNearv2->SetLineWidth(2);


    TGraphErrors* yieldshhNearSystv2 = (TGraphErrors*)yieldshhNearSyst->Clone("yieldshhNearSystv2");
    yieldshhNearSystv2->SetPoint(0, nchArray[0], 0.5*(nearhhYieldArray[0]+nearhhYieldArrayv2[0]));
    yieldshhNearSystv2->SetPoint(1, nchArray[1], 0.5*(nearhhYieldArray[1]+nearhhYieldArrayv2[1]));
    yieldshhNearSystv2->SetPoint(2, nchArray[2], 0.5*(nearhhYieldArray[2]+nearhhYieldArrayv2[2]));
    yieldshhNearSystv2->SetPointError(0, 2.0, -0.5*(nearhhYieldArray[0]-nearhhYieldArrayv2[0]));
    yieldshhNearSystv2->SetPointError(1, 2.0, -0.5*(nearhhYieldArray[1]-nearhhYieldArrayv2[1]));
    yieldshhNearSystv2->SetPointError(2, 2.0, -0.5*(nearhhYieldArray[2]-nearhhYieldArrayv2[2]));
    yieldshhNearSystv2->SetFillColor(kRed+2);
    yieldshhNearSystv2->SetFillStyle(3254);
//    yieldshhNearSystv2->SetLineWidth(0);



    TGraphErrors* yieldshhAway = new TGraphErrors(3, nchArray, awayhhYieldArray, nchArrayErr, awayhhYieldArrayErr);
    yieldshhAway->SetMarkerStyle(25);
    yieldshhAway->SetMarkerSize(2);
    yieldshhAway->SetMarkerColor(kBlue+2);
    yieldshhAway->SetLineColor(kBlue+3);
    yieldshhAway->SetLineWidth(2);

    TGraphErrors* yieldshhAwaySyst = new TGraphErrors(3, nchArray, awayhhYieldArray, nchArraySystErr, awayhhYieldArraySystErr);
    yieldshhAwaySyst->SetMarkerStyle(25);
    yieldshhAwaySyst->SetMarkerSize(2);
    yieldshhAwaySyst->SetMarkerColor(kBlue+2);
    yieldshhAwaySyst->SetLineColor(kBlue+3);
    yieldshhAwaySyst->SetLineWidth(2);
//    yieldshhAwaySyst->SetFillColor(kWhite);
//
    TGraphErrors* yieldshhAwayv2 = new TGraphErrors(3, nchArray, awayhhYieldArrayv2, nchArrayErr, awayhhYieldArrayErr);
    yieldshhAwayv2->SetMarkerStyle(25);
    yieldshhAwayv2->SetMarkerSize(2);
    yieldshhAwayv2->SetMarkerColor(kBlue+2);
    yieldshhAwayv2->SetLineColor(kBlue+3);
    yieldshhAwayv2->SetLineWidth(2);


    TGraphErrors* yieldshhAwaySystv2 = (TGraphErrors*)yieldshhAwaySyst->Clone("yieldshhAwaySystv2");
    yieldshhAwaySystv2->SetPoint(0, nchArray[0], 0.5*(awayhhYieldArray[0]+awayhhYieldArrayv2[0]));
    yieldshhAwaySystv2->SetPoint(1, nchArray[1], 0.5*(awayhhYieldArray[1]+awayhhYieldArrayv2[1]));
    yieldshhAwaySystv2->SetPoint(2, nchArray[2], 0.5*(awayhhYieldArray[2]+awayhhYieldArrayv2[2]));
    yieldshhAwaySystv2->SetPointError(0, 2.0, -0.5*(awayhhYieldArray[0]-awayhhYieldArrayv2[0]));
    yieldshhAwaySystv2->SetPointError(1, 2.0, -0.5*(awayhhYieldArray[1]-awayhhYieldArrayv2[1]));
    yieldshhAwaySystv2->SetPointError(2, 2.0, -0.5*(awayhhYieldArray[2]-awayhhYieldArrayv2[2]));
    yieldshhAwaySystv2->SetFillColor(kBlue+2);
    yieldshhAwaySystv2->SetFillStyle(3245);
//    yieldshhAwaySystv2->SetLineWidth(0);


    TGraphErrors* yieldshhBulk = new TGraphErrors(3, nchArray, bulkhhYieldArray, nchArrayErr, bulkhhYieldArrayErr);
    yieldshhBulk->SetMarkerStyle(26);
    yieldshhBulk->SetMarkerSize(3);
    yieldshhBulk->SetMarkerColor(kGreen+3);
    yieldshhBulk->SetLineColor(kGreen+4);
    yieldshhBulk->SetLineWidth(4);

    TGraphErrors* yieldshhBulkSyst = new TGraphErrors(3, nchArray, bulkhhYieldArray, nchArraySystErr, bulkhhYieldArraySystErr);
    yieldshhBulkSyst->SetMarkerStyle(26);
    yieldshhBulkSyst->SetMarkerSize(2);
    yieldshhBulkSyst->SetMarkerColor(kGreen+3);
    yieldshhBulkSyst->SetLineColor(kGreen+4);
    yieldshhBulkSyst->SetLineWidth(2);
//    yieldshhBulkSyst->SetFillColor(kWhite);

   
    TGraphErrors* yieldshhTot = new TGraphErrors(3, nchArray, totalhhYieldArray, nchArrayErr, totalhhYieldArrayErr);
    yieldshhTot->SetMarkerStyle(30);
    yieldshhTot->SetMarkerSize(3);
    yieldshhTot->SetMarkerColor(kMagenta+3);
    yieldshhTot->SetLineColor(kMagenta+4);
    yieldshhTot->SetLineWidth(2);
    yieldshhTot->SetFillColor(kMagenta+1);
    yieldshhTot->SetFillStyle(3144);

    TGraphErrors* yieldshhTotSyst = new TGraphErrors(3, nchArray, totalhhYieldArray, nchArraySystErr, totalhhYieldArraySystErr);
    yieldshhTotSyst->SetMarkerStyle(30);
    yieldshhTotSyst->SetMarkerSize(3);
    yieldshhTotSyst->SetMarkerColor(kMagenta+3);
    yieldshhTotSyst->SetLineColor(kMagenta+4);
    yieldshhTotSyst->SetLineWidth(2);
    //yieldshhTot->SetFillColor(kMagenta+1);
    //yieldshhTot->SetFillStyle(3144);

    // Do Fitting of jet yields for v2 and flat UE
//    TF1* yfit = new TF1("yfit", "[0]+[1]*x", 0, 65);
    TF1* yfit = new TF1("yfit", "[0]", 0, 65);
    TFitResultPtr yn1 = yieldsNear->Fit(yfit, "0S");
    printf("near fit %e %e %f, prob: %f\n", yfit->GetParameter(1)/250., yfit->GetParError(1)/250., yfit->GetChisquare(), yn1->Prob());
    printf("%.2e\% $\\pm$ %2.e (%2.e\% $\\pm$ %2.e)\n", 100.*yfit->GetParameter(1)/yieldsNear->GetPointY(0), 100.*yfit->GetParError(1)/yieldsNear->GetPointY(0));
    TFitResultPtr yn2 = yieldsNearv2->Fit(yfit, "0S");
    printf("nearv2 fit %e %e %f, prob: %f\n", yfit->GetParameter(1)/250., yfit->GetParError(1)/250., yfit->GetChisquare(), yn2->Prob());
    printf("%.2e\% $\\pm$ %2.e (%2.e\% $\\pm$ %2.e)\n", 100.*yfit->GetParameter(1)/yieldsNearv2->GetPointY(0), 100.*yfit->GetParError(1)/yieldsNearv2->GetPointY(0));
    TFitResultPtr ya1 = yieldsAway->Fit(yfit, "0S");
    printf("away fit %e %e %f\n, prob: %f\n", yfit->GetParameter(1)/250., yfit->GetParError(1)/250., yfit->GetChisquare(), ya1->Prob());
    printf("%.2e\% $\\pm$ %2.e (%2.e\% $\\pm$ %2.e)\n", 100.*yfit->GetParameter(1)/yieldsAway->GetPointY(0), 100.*yfit->GetParError(1)/yieldsAway->GetPointY(0));
    TFitResultPtr ya2 = yieldsAwayv2->Fit(yfit, "0S");
    printf("awayv2 fit %e %e %f, prob: %f\n", yfit->GetParameter(1)/250., yfit->GetParError(1)/250., yfit->GetChisquare(), ya2->Prob());
    printf("%.2f\% $\\pm$ %2.e (%2.e\% $\\pm$ %2.e)\n", 100.*yfit->GetParameter(1)/yieldsAwayv2->GetPointY(0), 100.*yfit->GetParError(1)/yieldsAwayv2->GetPointY(0));

    TF1* yfit3 = new TF1("yfit3", "[0]+[1]*x", 0, 65);
    yfit3->SetParameter(0, 0.1);
    yieldshhNear->Fit(yfit3, "0");
    printf("hh near fit %e %e %f\n", yfit3->GetParameter(1), yfit3->GetParError(1), yfit3->GetChisquare());
    printf("%.2f%% $\\pm$ %3.f (%2.f%% $\\pm$ %3.f)\n", 100.*yfit3->GetParameter(1)/yieldshhNear->GetPointY(0), 100.*yfit3->GetParError(1)/yieldshhNear->GetPointY(0));
    yieldshhNearv2->Fit(yfit3, "0");
    printf("hh nearv2 fit %e %e %f\n", yfit3->GetParameter(1), yfit3->GetParError(1), yfit3->GetChisquare());
    printf("%.2f%% $\\pm$ %3.f (%2.f%% $\\pm$ %3.f)\n", 100.*yfit3->GetParameter(1)/yieldshhNearv2->GetPointY(0), 100.*yfit3->GetParError(1)/yieldshhNearv2->GetPointY(0));

    TF1* yfit2 = new TF1("yfit2", "[0]+[1]*x", 0, 65);
    yieldshhAway->Fit(yfit2, "0");
    printf("hh away fit %e %e %f\n", yfit2->GetParameter(1), yfit2->GetParError(1), yfit2->GetChisquare());
    printf("%.2f%% $\\pm$ %2.f (%2.f%% $\\pm$ %2.f)\n", 100.*yfit2->GetParameter(1)/yieldshhAway->GetPointY(0), 100.*yfit2->GetParError(1)/yieldshhAway->GetPointY(0));
    yfit->SetParameter(0, 0.1);
    yieldshhAwayv2->Fit(yfit2, "0");
    printf("hh awayv2 fit %e %e %f\n", yfit2->GetParameter(1), yfit2->GetParError(1), yfit2->GetChisquare());
    printf("%.2f%% $\\pm$ %3.f (%2.f%% $\\pm$ %3.f)\n", 100.*yfit2->GetParameter(1)/yieldshhAwayv2->GetPointY(0), 100.*yfit2->GetParError(1)/yieldshhAwayv2->GetPointY(0));
 
  //print out jet yield percent change and errors (and for v2 assumption)  
    double nearhphipercent = 100.*(nearhPhiYieldArray[2]/nearhPhiYieldArray[0]-1);
    double nearhphipercentErr = 100.*nearhPhiYieldArray[2]/nearhPhiYieldArray[0]*(TMath::Sqrt(TMath::Power(nearhPhiYieldArrayErr[2]/nearhPhiYieldArray[2], 2) + TMath::Power(nearhPhiYieldArrayErr[0]/nearhPhiYieldArray[0], 2)));

    double nearhphiv2percent = 100.*(nearhPhiYieldArrayv2[2]/nearhPhiYieldArrayv2[0]-1);
    double nearhphiv2percentErr = 100.*nearhPhiYieldArrayv2[2]/nearhPhiYieldArrayv2[0]*(TMath::Sqrt(TMath::Power(nearhPhiYieldArrayErr[2]/nearhPhiYieldArrayv2[2], 2) + TMath::Power(nearhPhiYieldArrayErr[0]/nearhPhiYieldArrayv2[0], 2)));

    double awayhphipercent = 100.*(awayhPhiYieldArray[2]/awayhPhiYieldArray[0]-1);
    double awayhphipercentErr = 100.*awayhPhiYieldArray[2]/awayhPhiYieldArray[0]*(TMath::Sqrt(TMath::Power(awayhPhiYieldArrayErr[2]/awayhPhiYieldArray[2], 2) + TMath::Power(awayhPhiYieldArrayErr[0]/awayhPhiYieldArray[0], 2)));

    double awayhphiv2percent = 100.*(awayhPhiYieldArrayv2[2]/awayhPhiYieldArrayv2[0]-1);
    double awayhphiv2percentErr = 100.*awayhPhiYieldArrayv2[2]/awayhPhiYieldArrayv2[0]*(TMath::Sqrt(TMath::Power(awayhPhiYieldArrayErr[2]/awayhPhiYieldArrayv2[2], 2) + TMath::Power(awayhPhiYieldArrayErr[0]/awayhPhiYieldArrayv2[0], 2)));


    double nearhhpercent = 100.*(nearhhYieldArray[2]/nearhhYieldArray[0]-1);
    double nearhhpercentErr = 100.*nearhhYieldArray[2]/nearhhYieldArray[0]*(TMath::Sqrt(TMath::Power(nearhhYieldArrayErr[2]/nearhhYieldArray[2], 2) + TMath::Power((nearhhYieldArrayErr[0] + 0.05*nearhhYieldArray[0])/nearhhYieldArray[0], 2)));

    double nearhhv2percent = 100.*(nearhhYieldArrayv2[2]/nearhhYieldArrayv2[0]-1);
    double nearhhv2percentErr = 100.*nearhhYieldArrayv2[2]/nearhhYieldArrayv2[0]*(TMath::Sqrt(TMath::Power(nearhhYieldArrayErr[2]/nearhhYieldArrayv2[2], 2) + TMath::Power((nearhhYieldArrayErr[0] + 0.05*nearhhYieldArrayv2[0])/nearhhYieldArrayv2[0], 2)));

    double awayhhpercent = 100.*(awayhhYieldArray[2]/awayhhYieldArray[0]-1);
    double awayhhpercentErr = 100.*awayhhYieldArray[2]/awayhhYieldArray[0]*(TMath::Sqrt(TMath::Power(awayhhYieldArrayErr[2]/awayhhYieldArray[2], 2) + TMath::Power((awayhhYieldArrayErr[0] + 0.05*awayhhYieldArray[0])/awayhhYieldArray[0], 2)));

    double awayhhv2percent = 100.*(awayhhYieldArrayv2[2]/awayhhYieldArrayv2[0]-1);
    double awayhhv2percentErr = 100.*awayhhYieldArrayv2[2]/awayhhYieldArrayv2[0]*(TMath::Sqrt(TMath::Power((awayhhYieldArrayErr[2])/awayhhYieldArrayv2[2], 2) + TMath::Power((awayhhYieldArrayErr[0]+ 0.05*awayhhYieldArrayv2[0])/awayhhYieldArrayv2[0], 2)));

    printf("%.1f\\%% $\\pm$ %.2f\\%% (%.1f\\%% $\\pm$ %.2f\\%%)\n", nearhphipercent, nearhphipercentErr, nearhphiv2percent, nearhphiv2percentErr);

    printf("%.1f\\%% $\\pm$ %.2f\\%% (%.1f\\%% $\\pm$ %.2f\\%%)\n", awayhphipercent, awayhphipercentErr, awayhphiv2percent, awayhphiv2percentErr);

    printf("%.1f\\%% $\\pm$ %.2f\\%% (%.1f\\%% $\\pm$ %.2f\\%%)\n", nearhhpercent, nearhhpercentErr, nearhhv2percent, nearhhv2percentErr);

    printf("%.1f\\%% $\\pm$ %.2f\\%% (%.1f\\%% $\\pm$ %.2f\\%%)\n", awayhhpercent, awayhhpercentErr, awayhhv2percent, awayhhv2percentErr);

    //set-up jet/total ratio histograms for h-phi and h-h

    TH1D* jet2totalhPhi = (TH1D*)ratioNearHist->Clone("jet2totalhPhi");
    TH1D* jet2totalhh = (TH1D*)ratioNearHist->Clone("jet2totalhh");

    for(int i = 0; i < 3; i++){
        jet2totalhPhi->SetBinContent(i+2, (nearhPhiYieldArray[i] + awayhPhiYieldArray[i])/totalhPhiYieldArray[i]);
        jet2totalhh->SetBinContent(i+2, (nearhhYieldArray[i] + awayhhYieldArray[i])/totalhhYieldArray[i]);
    }

    TH1D* jet2totalhPhiSyst = (TH1D*)jet2totalhPhi->Clone("jet2totalhPhiSyst");
    TH1D* jet2totalhhSyst = (TH1D*)jet2totalhh->Clone("jet2totalhhSyst");

    for(int i = 0; i < 3; i++){
        jet2totalhPhiSyst->SetBinError(i+2, TMath::Sqrt(TMath::Power(nearhPhiYieldArraySystErr[i]/nearhPhiYieldArray[i], 2.0) + TMath::Power(totalhPhiYieldArraySystErr[i]/totalhPhiYieldArray[i], 2.0)));
        jet2totalhhSyst->SetBinError(i+2, TMath::Sqrt(TMath::Power(nearhhYieldArraySystErr[i]/nearhhYieldArray[i], 2.0) + TMath::Power(totalhhYieldArraySystErr[i]/totalhhYieldArray[i], 2.0)));
    }

    jet2totalhPhi->SetLineColor(kMagenta+2);
    jet2totalhPhiSyst->SetLineColor(kMagenta+2);
    jet2totalhPhiSyst->SetFillStyle(0);
    jet2totalhPhi->SetMarkerColor(kMagenta+2);
    jet2totalhPhi->GetYaxis()->SetTitle("(Jet/Total) pair-yield Ratio");

    jet2totalhh->SetLineColor(kCyan+2);
    jet2totalhhSyst->SetLineColor(kCyan+2);
    jet2totalhhSyst->SetFillStyle(0);
    jet2totalhh->SetMarkerColor(kCyan+2);

    TLegend *jetlegend = new TLegend(0.186, 0.759, 0.529, 0.868);
    jetlegend->AddEntry(jet2totalhPhi, "Jet/Total (h-#phi)", "p");
    jetlegend->AddEntry(jet2totalhh, "Jet/Total (h-h)", "p");
    jetlegend->SetLineWidth(0);

    TCanvas* jet2totalcanvas = new TCanvas("jet2totalcanvas", "jet2totalcanvas", 55, 55, 900, 600);
    jet2totalcanvas->cd();
    jet2totalcanvas->SetMargin(0.126, 0.05, 0.125, 0.05);
    gStyle->SetErrorX(0.5);
    jet2totalhPhi->Draw("AXIS");
    jet2totalhPhi->GetXaxis()->SetLabelOffset(999);
    jet2totalhPhi->GetXaxis()->SetTickSize(0.0);
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
    //jet2totalhPhiSyst->Draw("E2 SAME");
    jet2totalhPhi->Draw("P E SAME");
    //jet2totalhhSyst->Draw("E2 SAME");
    jet2totalhh->Draw("P E SAME");
    jetlegend->Draw();
    //data->Draw();
    text2->Draw();

    //draw near-side yields
    TLegend *nearYieldMultlegend = new TLegend(0.164, 0.527, 0.498, 0.908);
    nearYieldMultlegend->SetMargin(0.12);
    nearYieldMultlegend->AddEntry(yieldshhNear, "(h#font[122]{-}h) in near-side jet", "p");
    nearYieldMultlegend->AddEntry(yieldsNear, Form("(h#font[122]{-}#phi) in near-side jet (#times %d)", int(phijetscale)), "p");
    //nearYieldMultlegend->AddEntry(pyhhnear, "PYTHIA8 (h-h) in Near Jet", "pl");
    //nearYieldMultlegend->AddEntry(pyhphinear, Form("PYTHIA8 (h-#phi) in Near Jet (#times %d)", int(phijetscale)), "pl");
    nearYieldMultlegend->AddEntry(yieldsNearSystv2, "yields with #it{v}_{2} assumption", "f");
    nearYieldMultlegend->AddEntry(yieldshhAway, "(h#font[122]{-}h) in away-side jet", "p");
    nearYieldMultlegend->AddEntry(yieldsAway, Form("(h#font[122]{-}#phi) in away-side jet (#times %d)", int(phijetscale)), "p");
    nearYieldMultlegend->AddEntry(yieldsAwaySystv2, "yields with #it{v}_{2} assumption", "f");

     
    nearYieldMultlegend->SetLineWidth(0);

    TH1D* nchJetYieldHist = (TH1D*) nchratioNearHist->Clone("nchJetYieldHist");
    nchJetYieldHist->GetYaxis()->SetTitle("Per-trigger pair yields");
    nchJetYieldHist->GetYaxis()->SetRangeUser(0.0, 2.0);

    TCanvas* nearYieldvsMultCanvas = new TCanvas("nearYieldvsMultCanvas", "nearYieldvsMultCanvas", 55, 55, 900, 600);
    nearYieldvsMultCanvas->cd();
    nearYieldvsMultCanvas->SetMargin(0.126, 0.05, 0.125, 0.05);
    nchJetYieldHist->GetYaxis()->SetTitleOffset(1.0);
    nchJetYieldHist->GetYaxis()->SetDecimals(kTRUE);
    //gStyle->SetErrorX(0.5);
    //ratioJetHist->Draw("AXIS");
    nchJetYieldHist->Draw("AXIS");
    //ratioJetHist->GetYaxis()->SetTitle("Per-trigger Pair Yields");
/*    ratioJetHist->GetXaxis()->SetLabelOffset(999);
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
*/
    //yieldsNearSyst->Draw("5");   
    yieldsNearSystv2->Draw("3");
    yieldsNear->Draw("PZL");
    yieldsNearv2->Draw("LX");
    //yieldshhNearSyst->Draw("5");
    yieldshhNear->Draw("PZL");
    yieldshhNearSystv2->Draw("3");
    yieldshhNearv2->Draw("LX");
    //yieldsAwaySyst->Draw("5");   
    yieldsAway->Draw("PZL");
    yieldsAwaySystv2->Draw("3");
    yieldsAwayv2->Draw("LX");
    //yieldshhAwaySyst->Draw("5");
    yieldshhAway->Draw("PZL");
    yieldshhAwaySystv2->Draw("3");
    yieldshhAwayv2->Draw("LX");
    nearYieldMultlegend->Draw("SAME");
   data->Draw();
    //pyhphinear->Draw("P");
    //pyhphiaway->Draw("P");
    //pyhhnear->Draw("P");
    //pyhhaway->Draw("P");
    text2gen->Draw();

    //draw away-side yields
    TLegend *awayYieldMultlegend = new TLegend(0.1737, 0.781, 0.533, 0.919);
    awayYieldMultlegend->SetMargin(0.20);
    awayYieldMultlegend->AddEntry(yieldshhAway, "(h#font[122]{-}h) in Away Jet", "p");
    awayYieldMultlegend->AddEntry(yieldsAway, "(h#font[122]{-}#phi) in Away Jet (#times 300)", "p");
    awayYieldMultlegend->SetLineWidth(0);


    TCanvas* awayYieldvsMultCanvas = new TCanvas("awayYieldvsMultCanvas", "awayYieldvsMultCanvas", 55, 55, 900, 600);
    awayYieldvsMultCanvas->cd();
    awayYieldvsMultCanvas->SetMargin(0.126, 0.05, 0.125, 0.05);
    gStyle->SetErrorX(0.5);
    ratioJetHist->Draw("AXIS");
    ratioJetHist->GetYaxis()->SetTitle("Per-trigger Pair Yields");
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
    yieldsAwaySyst->Draw("5");   
    yieldsAway->Draw("P");
    yieldshhAwaySyst->Draw("5");
    yieldshhAway->Draw("P");
    data->Draw();
    text2->Draw();
    awayYieldMultlegend->Draw("SAME");


    //draw bulk yields
    TLegend *bulkYieldMultlegend = new TLegend(0.1737, 0.781, 0.533, 0.919);
    bulkYieldMultlegend->SetMargin(0.30);
    bulkYieldMultlegend->AddEntry(yieldshhBulk, "(h#font[122]{-}h) in UE", "p");
    bulkYieldMultlegend->AddEntry(yieldsBulk, Form("(h#font[122]{-}#phi) in UE (#times %d)", int(phitotscale)), "p");
    bulkYieldMultlegend->SetLineWidth(0);


    TLegend *totalYieldMultlegend = new TLegend(0.1737, 0.781, 0.533, 0.919);
    totalYieldMultlegend->SetMargin(0.30);
    totalYieldMultlegend->AddEntry(yieldshhTot, "Total (h#font[122]{-}h) pairs", "p");
    totalYieldMultlegend->AddEntry(yieldsTot, Form("Total (h#font[122]{-}#phi) pairs (#times %d)", int(phitotscale)), "p");
    totalYieldMultlegend->SetLineWidth(0);


    TCanvas* bulkYieldvsMultCanvas = new TCanvas("bulkYieldvsMultCanvas", "bulkYieldvsMultCanvas", 55, 55, 900, 600);
    bulkYieldvsMultCanvas->cd();
    bulkYieldvsMultCanvas->SetMargin(0.126, 0.05, 0.125, 0.05);
    nchJetYieldHist->Draw("AXIS");
    yieldsBulkSyst->Draw("5");   
    yieldsBulk->Draw("PC");
    yieldshhBulkSyst->Draw("5");
    yieldshhBulk->Draw("PC");
    //yieldsTotSyst->Draw("5");
    //yieldsTot->Draw("P");
    //yieldshhTotSyst->Draw("5");
    //yieldshhTot->Draw("P");
    data->Draw();
    text2->Draw();
    bulkYieldMultlegend->Draw("SAME");
    //totalYieldMultlegend->Draw("SAME");

    //draw total yields
    
    TCanvas* totalYieldvsMultCanvas = new TCanvas("totalYieldvsMultCanvas", "totalYieldvsMultCanvas", 55, 55, 900, 600);
    totalYieldvsMultCanvas->cd();
    totalYieldvsMultCanvas->SetMargin(0.126, 0.05, 0.125, 0.05);
    gStyle->SetErrorX(0.5);
    ratioJetHist->Draw("AXIS");
    ratioJetHist->GetYaxis()->SetTitle("Per-trigger Pair Yields");
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
    yieldsTotSyst->Draw("5");
    yieldsTot->Draw("P");
    yieldshhTotSyst->Draw("5");
    yieldshhTot->Draw("P");
    //data->Draw();
    //text2->Draw();
    //totalYieldMultlegend->Draw("SAME");


    //Double Ratio plots:
    Double_t doublenear020 = near020/mid020;
    Double_t doublenear020Er = doublenear020*TMath::Sqrt(TMath::Power((near020Er)/(near020), 2) + TMath::Power((mid020Er)/(mid020), 2));
    Double_t doubleaway020 = away020/mid020;
    Double_t doubleaway020Er = doubleaway020*TMath::Sqrt(TMath::Power((away020Er)/(away020), 2) + TMath::Power((mid020Er)/(mid020), 2));

    Double_t doublenear2050 = near2050/mid2050;
    Double_t doublenear2050Er = doublenear2050*TMath::Sqrt(TMath::Power((near2050Er)/(near2050), 2) + TMath::Power((mid2050Er)/(mid2050), 2));
    Double_t doubleaway2050 = away2050/mid2050;
    Double_t doubleaway2050Er = doubleaway2050*TMath::Sqrt(TMath::Power((away2050Er)/(away2050), 2) + TMath::Power((mid2050Er)/(mid2050), 2));

    Double_t doublenear50100 = near50100/mid50100;
    Double_t doublenear50100Er = doublenear50100*TMath::Sqrt(TMath::Power((near50100Er)/(near50100), 2) + TMath::Power((mid50100Er)/(mid50100), 2));
    Double_t doubleaway50100 = away50100/mid50100;
    Double_t doubleaway50100Er = doubleaway50100*TMath::Sqrt(TMath::Power((away50100Er)/(away50100), 2) + TMath::Power((mid50100Er)/(mid50100), 2));

    Double_t doublenear0100 = near0100/mid0100;
    Double_t doublenear0100Er = doublenear0100*TMath::Sqrt(TMath::Power((near0100Er)/(near0100), 2) + TMath::Power((mid0100Er)/(mid0100), 2));
    Double_t doubleaway0100 = away0100/mid0100;
    Double_t doubleaway0100Er = doubleaway0100*TMath::Sqrt(TMath::Power((away0100Er)/(away0100), 2) + TMath::Power((mid0100Er)/(mid0100), 2));

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

    TH1D *doubleratios50100 = new TH1D("doubleratios50100", "(h-#phi / h-h) Double Ratios", 2, 0, 2);
    doubleratios50100->GetXaxis()->SetBinLabel(1, "#frac{near-side}{mid}");
    doubleratios50100->SetBinContent(1, doublenear50100);
    doubleratios50100->SetBinError(1, doublenear50100Er);
    doubleratios50100->GetXaxis()->SetBinLabel(2, "#frac{away-side}{mid}");
    doubleratios50100->SetBinContent(2, doubleaway50100);
    doubleratios50100->SetBinError(2, doubleaway50100Er);
    doubleratios50100->SetMarkerStyle(20);
    doubleratios50100->SetMarkerColor(kGreen+2);
    doubleratios50100->SetLineColor(kGreen+2);
    doubleratios50100->SetLineWidth(2);
    doubleratios50100->SetMarkerSize(2);

    TH1D *doubleratios0100 = new TH1D("doubleratios0100", "(h-#phi / h-h) Double Ratios", 2, 0, 2);
    doubleratios0100->GetXaxis()->SetBinLabel(1, "#frac{near-side}{mid}");
    doubleratios0100->SetBinContent(1, doublenear0100);
    doubleratios0100->SetBinError(1, doublenear0100Er);
    doubleratios0100->GetXaxis()->SetBinLabel(2, "#frac{away-side}{mid}");
    doubleratios0100->SetBinContent(2, doubleaway0100);
    doubleratios0100->SetBinError(1, doubleaway0100Er);
    doubleratios0100->SetMarkerStyle(29);
    doubleratios0100->SetMarkerColor(kViolet-1);
    doubleratios0100->SetLineColor(kViolet-1);
    doubleratios0100->SetLineWidth(2);
    doubleratios0100->SetMarkerSize(2);


    TCanvas *testDoublec = new TCanvas("testdouble", "testdouble",50, 50, 600, 600);
    testDoublec->cd();
    //doubleratios2050->GetYaxis()->SetRangeUser(0.0000, 0.0040);
    doubleratios2050->GetXaxis()->SetLabelSize(0.07);
    doubleratios2050->Draw("P SAME");
    doubleratios50100->GetXaxis()->SetLimits(-0.05, 1.95);
    doubleratios50100->Draw("P SAME");
    doubleratios020->GetXaxis()->SetLimits(0.05, 2.05);
    doubleratios020->Draw("P SAME");
    doubleratios0100->GetXaxis()->SetLimits(-0.1,1.9);
    doubleratios0100->Draw("P SAME");
    //line->Draw("SAME");
    ratioslegend->Draw("SAME");

    printf("\n\n");
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

    TH1D *yields50100hPhi = new TH1D("yields50100hPhiEff", "h-#phi Per Trigger Yields", 4, 0, 4);
    yields50100hPhi->GetXaxis()->SetBinLabel(1, "near-side");
    yields50100hPhi->SetBinContent(1, near50_100hPhiYield);
    yields50100hPhi->SetBinError(1, near50_100hPhiError);
    yields50100hPhi->GetXaxis()->SetBinLabel(2, "mid");
    yields50100hPhi->SetBinContent(2, mid50_100hPhiYield);
    yields50100hPhi->SetBinError(2, mid50_100hPhiError);
    yields50100hPhi->GetXaxis()->SetBinLabel(3, "away-side");
    yields50100hPhi->SetBinContent(3, away50_100hPhiYield);
    yields50100hPhi->SetBinError(3, away50_100hPhiError);
    yields50100hPhi->GetXaxis()->SetBinLabel(4, "total");
    yields50100hPhi->SetBinContent(4, total50_100hPhiYield);
    yields50100hPhi->SetBinError(4, total50_100hPhiError);
    yields50100hPhi->SetMarkerStyle(20);
    yields50100hPhi->SetMarkerColor(kGreen+2);
    yields50100hPhi->SetLineColor(kGreen+2);
    yields50100hPhi->SetMarkerSize(2);
    yields50100hPhi->SetLineWidth(2);

    TH1D *yields0100hPhi = new TH1D("yields0100hPhi", "h-#phi Per Trigger Yields", 4, 0, 4);
    yields0100hPhi->GetXaxis()->SetBinLabel(1, "near-side");
    yields0100hPhi->SetBinContent(1, near0_100hPhiYield);
    yields0100hPhi->SetBinError(1, near0_100hPhiError);
    yields0100hPhi->GetXaxis()->SetBinLabel(2, "mid");
    yields0100hPhi->SetBinContent(2, mid0_100hPhiYield);
    yields0100hPhi->SetBinError(2, mid0_100hPhiError);
    yields0100hPhi->GetXaxis()->SetBinLabel(3, "away-side");
    yields0100hPhi->SetBinContent(3, away0_100hPhiYield);
    yields0100hPhi->SetBinError(3, away0_100hPhiError);
    yields0100hPhi->GetXaxis()->SetBinLabel(4, "total");
    yields0100hPhi->SetBinContent(4, total0_100hPhiYield);
    yields0100hPhi->SetBinError(4, total0_100hPhiError);
    yields0100hPhi->SetMarkerStyle(29);
    yields0100hPhi->SetMarkerColor(kViolet-1);
    yields0100hPhi->SetLineColor(kViolet-1);
    yields0100hPhi->SetMarkerSize(2);
    yields0100hPhi->SetLineWidth(2);


    TCanvas *yieldhPhic = new TCanvas("yieldhPhi", "yieldhPhi",50, 50, 600, 600);
    yieldhPhic->cd();
    yields2050hPhi->GetYaxis()->SetRangeUser(0.0000, 0.0040);
    yields2050hPhi->GetXaxis()->SetLabelSize(0.07);
    yields2050hPhi->Draw("P SAME");
    yields50100hPhi->GetXaxis()->SetLimits(0.05, 4.05);
    yields50100hPhi->Draw("P SAME");
    yields020hPhi->GetXaxis()->SetLimits(-0.05, 3.95);
    yields020hPhi->Draw("P SAME");
    yields0100hPhi->GetXaxis()->SetLimits(-0.1, 3.9);
    yields0100hPhi->Draw("P SAME");
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

    TH1D *yields50100hh = new TH1D("yields50100hhEff", "h-#phi Per Trigger Yields", 4, 0, 4);
    yields50100hh->GetXaxis()->SetBinLabel(1, "near-side");
    yields50100hh->SetBinContent(1, near50_100hhYield);
    yields50100hh->SetBinError(1, near50_100hhError);
    yields50100hh->GetXaxis()->SetBinLabel(2, "mid");
    yields50100hh->SetBinContent(2, mid50_100hhYield);
    yields50100hh->SetBinError(2, mid50_100hhError);
    yields50100hh->GetXaxis()->SetBinLabel(3, "away-side");
    yields50100hh->SetBinContent(3, away50_100hhYield);
    yields50100hh->SetBinError(3, away50_100hhError);
    yields50100hh->GetXaxis()->SetBinLabel(4, "total");
    yields50100hh->SetBinContent(4, total50_100hhYield);
    yields50100hh->SetBinError(4, total50_100hhError);
    yields50100hh->SetMarkerStyle(20);
    yields50100hh->SetMarkerColor(kGreen+2);
    yields50100hh->SetLineColor(kGreen+2);
    yields50100hh->SetMarkerSize(2);
    yields50100hh->SetLineWidth(2);

    TH1D *yields0100hh = new TH1D("yields0100hh", "h-#phi Per Trigger Yields", 4, 0, 4);
    yields0100hh->GetXaxis()->SetBinLabel(1, "near-side");
    yields0100hh->SetBinContent(1, near0_100hhYield);
    yields0100hh->SetBinError(1, near0_100hhError);
    yields0100hh->GetXaxis()->SetBinLabel(2, "mid");
    yields0100hh->SetBinContent(2, mid0_100hhYield);
    yields0100hh->SetBinError(2, mid0_100hhError);
    yields0100hh->GetXaxis()->SetBinLabel(3, "away-side");
    yields0100hh->SetBinContent(3, away0_100hhYield);
    yields0100hh->SetBinError(3, away0_100hhError);
    yields0100hh->GetXaxis()->SetBinLabel(4, "total");
    yields0100hh->SetBinContent(4, total0_100hhYield);
    yields0100hh->SetBinError(4, total0_100hhError);
    yields0100hh->SetMarkerStyle(29);
    yields0100hh->SetMarkerColor(kViolet-1);
    yields0100hh->SetLineColor(kViolet-1);
    yields0100hh->SetMarkerSize(2);
    yields0100hh->SetLineWidth(2);

    TLine *linehh = new TLine(3.0, 0.0, 3.0, 10.000);
    linehh->SetLineStyle(7);
    linehh->SetLineWidth(2);

    TCanvas *yieldhhc = new TCanvas("yieldhh", "yieldhh",50, 50, 600, 600);
    yieldhhc->cd();
    yields2050hh->GetYaxis()->SetRangeUser(0.0000, 10.000);
    yields2050hh->GetXaxis()->SetLabelSize(0.07);
    yields2050hh->Draw("P SAME");
    yields50100hh->GetXaxis()->SetLimits(0.05, 4.05);
    yields50100hh->Draw("P SAME");
    yields020hh->GetXaxis()->SetLimits(-0.05, 3.95);
    yields020hh->Draw("P SAME");
    yields0100hh->GetXaxis()->SetLimits(-0.1, 3.9);
    yields0100hh->Draw("P SAME");
    linehh->Draw("SAME");
    ratioslegend->Draw("SAME");



    
    TH1D *ratio = (TH1D*)hPhidphi_0_20->Clone("ratio");
    ratio->Divide(hhdphi_0_20);
    

    TLegend  *legend = new TLegend(0.3791, 0.1518, 0.8772, 0.2688);
    legend->SetMargin(0.15);
    legend->AddEntry(corrFit2, "Hadron-#phi(1020) Correlation", "l");
    legend->AddEntry(corrFit, "Hadron-hadron Correlations", "l");

    TPaveText *text = new TPaveText(0.4815, 0.7056, 0.8658, 0.8551, "NDC");
//    text->AddText("ALICE");
    text->AddText("ALICE p#font[122]{-}Pb #sqrt{s_{NN}} = 5 TeV");
    text->AddText("0%#font[122]{-}20% Multiplicity");
    text->SetTextSizePixels(20);
    text->SetBorderSize(0);
    text->SetFillColor(kWhite);

    TPaveText *text2050 = new TPaveText(0.4815, 0.7056, 0.8658, 0.8551, "NDC");
//    text2050->AddText("ALICE");
    text2050->AddText("ALICE p#font[122]{-}Pb #sqrt{s_{NN}} = 5 TeV");
    text2050->AddText("20%#font[122]{-}50% Multiplicity");
    text2050->SetTextSizePixels(20);
    text2050->SetBorderSize(0);
    text2050->SetFillColor(kWhite);

    TPaveText *text50100 = new TPaveText(0.4815, 0.7056, 0.8658, 0.8551, "NDC");
//    text50100->AddText("ALICE");
    text50100->AddText("ALICE p#font[122]{-}Pb #sqrt{s_{NN}} = 5 TeV");
    text50100->AddText("50%#font[122]{-}80% Multiplicity");
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
    TH1D* hPhidphi_50_100_syst = (TH1D*)hPhidphi_50_100->Clone("hPhidphi_50_100_syst");
    hPhidphi_50_100_syst->SetFillColor(kGray);

    TH1D* hhdphi_0_20_syst = (TH1D*)hhdphi_0_20->Clone("hhdphi_0_20_syst");
    hhdphi_0_20_syst->SetFillColor(kGray);
    TH1D* hhdphi_20_50_syst = (TH1D*)hhdphi_20_50->Clone("hhdphi_20_50_syst");
    hhdphi_20_50_syst->SetFillColor(kGray);
    TH1D* hhdphi_50_100_syst = (TH1D*)hhdphi_50_100->Clone("hhdphi_50_100_syst");
    hhdphi_50_100_syst->SetFillColor(kGray);

    for(int i = 1; i <= 16; i++){
        hPhidphi_0_20_syst->SetBinError(i, hPhidphi_0_20_syst->GetBinContent(i)*dphi020syst);
        hPhidphi_20_50_syst->SetBinError(i, hPhidphi_20_50_syst->GetBinContent(i)*dphi2050syst);
        hPhidphi_50_100_syst->SetBinError(i, hPhidphi_50_100_syst->GetBinContent(i)*dphi5080syst);
        hhdphi_0_20_syst->SetBinError(i, hhdphi_0_20_syst->GetBinContent(i)*hhdphi020syst);
        hhdphi_20_50_syst->SetBinError(i, hhdphi_20_50_syst->GetBinContent(i)*hhdphi2050syst);
        hhdphi_50_100_syst->SetBinError(i, hhdphi_50_100_syst->GetBinContent(i)*hhdphi5080syst);

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
    hphiBGv2->Draw("SAME");
    text->Draw();
    text2->Draw();
    //legend->Draw();


    
    TCanvas *c0_20pp = new TCanvas("c0_20pp", "c0_20pp", 50, 50, 550, 600);
    c0_20pp->cd();
    c0_20pp->SetMargin(0.12, 0.05, 0.1, 0.05);
    //hhdphi_0_20->GetYaxis()->SetRangeUser(0.001, 0.25);
    //hhdphi_0_20->GetYaxis()->SetTitle("Per Trigger (h-h) Pairs");
    hhdphi_0_20_syst->Draw("E2 SAME");
    hhdphi_0_20->Draw("E0 X0 HIST SAME");
    //hPhidphi_0_20->Draw("E0 X0 SAME");
    corrFit->Draw("SAME");
    //corrFit2->Draw("SAME");
    hhBG->Draw("SAME");
    hhBGv2->Draw("SAME");
    //hphiBG->Draw("SAME");
    text->Draw();
    text2hh->Draw();
    //legend->Draw();

    //"cartoon" plot showing different yield regions
    TH1D* regionHist = (TH1D*)hhdphi_0_20->Clone("regionHist");
    regionHist->SetLineWidth(3);
    regionHist->SetLineColor(kBlack);
    regionHist->SetMarkerColor(kBlack);
    regionHist->SetMarkerStyle(kFullCircle);
    regionHist->SetMarkerSize(1);
    regionHist->SetFillColor(kViolet-9);

    TH1D* nearHist = (TH1D*)regionHist->Clone("nearHist");
    nearHist->SetFillColor(kRed+1);
    TH1D* awayHist = (TH1D*)regionHist->Clone("awayHist");
    awayHist->SetFillColor(kBlue);
    TH1D* bulkHist = (TH1D*)regionHist->Clone("bulkHist");
    bulkHist->SetLineColor(kGreen+1);
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
    //regionLeg->AddEntry(nearHist, "Near-side Jet", "f");
    //regionLeg->AddEntry(awayHist, "Away-side Jet", "f");
    //regionLeg->AddEntry(bulkHist, "Underlying Event", "f");
    regionLeg->AddEntry(regionHist, "Total Yield", "f");
    TCanvas* cregions = new TCanvas("cregions", "cregions", 50, 50, 600, 600);
    cregions->cd();
    regionHist->GetYaxis()->SetRangeUser(0.0, 0.8);
    regionHist->Draw("HIST");
    //nearHist->Draw("HIST SAME");
    //awayHist->Draw("HIST SAME");
    //bulkHist->Draw("HIST SAME");
    regionHist->Draw("HIST P SAME");
    //regionHist->Draw("HIST SAME");
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
    hhBGv2_20_50->Draw("SAME");
    //hphiBG_20_50->Draw("SAME");
    text2050->Draw();
    text2hh->Draw();
    //legend->Draw();

    TCanvas *c50_100 = new TCanvas("c50_100", "c50_100", 50, 50, 550, 600);
    c50_100->cd();
    c50_100->SetMargin(0.12, 0.05, 0.1, 0.05);
    //hhdphi_50_100->Draw("E0 X0");
    hPhidphi_50_100_syst->Draw("E2 SAME");
    hPhidphi_50_100->Draw("E0 X0 HIST SAME");
    //corrFit50100->Draw("SAME");
    //corrFit2_50100->Draw("SAME");
    //hhBG_50_100->Draw("SAME");
    hphiBG_50_100->Draw("SAME");
    text50100->Draw();
    text2->Draw();
    //legend->Draw();

    TCanvas *c50_100pp = new TCanvas("c50_100pp", "c50_100pp", 50, 50, 550, 600);
    c50_100pp->cd();
    c50_100pp->SetMargin(0.12, 0.05, 0.1, 0.05);
    hhdphi_50_100_syst->Draw("E2 SAME");
    hhdphi_50_100->Draw("E0 X0 HIST SAME");
    //hPhidphi_50_100->Draw("E0 X0 SAME");
    //corrFit50100->Draw("SAME");
    //corrFit2_50100->Draw("SAME");
    hhBG_50_100->Draw("SAME");
    hhBGv2_50_100->Draw("SAME");
    //hphiBG_50_100->Draw("SAME");
    text50100->Draw();
    text2hh->Draw();
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

    


    TH1D* hphi5080 = (TH1D*)hPhidphi_50_100->Clone("hphi5080");
    hphi5080->SetMarkerColor(kAzure+1);
    hphi5080->SetMarkerSize(1.5);
    hphi5080->SetMarkerStyle(kFullCircle);
    hphi5080->SetLineColor(kAzure+1);
    hphi5080->Add(hphiBG_50_100, -1.0);
    hphi5080->Scale(1.0/(2.4*hphi5080->GetXaxis()->GetBinWidth(1))); //scale by delta-eta and bin width

    TH1D* hphi5080syst = (TH1D*)hPhidphi_50_100_syst->Clone("hphi5080syst");
    hphi5080syst->Add(hphiBG_50_100, -1.0);
    hphi5080syst->Scale(1.0/(2.4*hphi5080->GetXaxis()->GetBinWidth(1)));
    hphi5080syst->SetFillColor(kAzure+1);
    hphi5080syst->SetFillStyle(3002);
    hphi5080syst->SetMarkerSize(0);


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


    TH1D* hh5080 = (TH1D*)hhdphi_50_100->Clone("hh5080");
    hh5080->SetMarkerColor(kAzure+2);
    hh5080->SetMarkerSize(1.5);
    hh5080->SetMarkerStyle(kOpenCircle);
    hh5080->SetLineColor(kAzure+2);
    hh5080->Add(hhBG_50_100, -1.0);
    hh5080->Scale(1.0/(2.4*hh5080->GetXaxis()->GetBinWidth(1))); //scale by delta-eta and bin width

    TH1D* hh5080syst = (TH1D*)hhdphi_50_100_syst->Clone("hh5080syst");
    hh5080syst->Add(hhBG_50_100, -1.0);
    hh5080syst->Scale(1.0/(2.4*hh5080->GetXaxis()->GetBinWidth(1)));
    hh5080syst->SetFillColor(kAzure+2);
    hh5080syst->SetFillStyle(3002);
    hh5080syst->SetMarkerSize(0);



    // alternative plots, not BG subtracted
    TH1D* althphi020 = (TH1D*)hPhidphi_0_20->Clone("althphi020");
    althphi020->SetMarkerColor(kRed);
    althphi020->SetLineColor(kRed);
    althphi020->SetMarkerSize(1.5);
    althphi020->SetMarkerStyle(kFullCross);
    althphi020->Scale(1.0/(2.4*althphi020->GetXaxis()->GetBinWidth(1))); //scale by delta-eta and bin width

    TH1D* althphi020syst = (TH1D*)hPhidphi_0_20_syst->Clone("althphi020syst");
    althphi020syst->SetLineColor(kRed);
    althphi020syst->Scale(1.0/(2.4*althphi020->GetXaxis()->GetBinWidth(1)));
    althphi020syst->SetFillColor(kRed);
    althphi020syst->SetFillStyle(0);
    althphi020syst->SetMarkerSize(0);

    TH1D* althphi2050 = (TH1D*)hPhidphi_20_50->Clone("althphi2050");
    althphi2050->SetMarkerColor(kOrange+1);
    althphi2050->SetLineColor(kOrange+1);
    althphi2050->SetMarkerSize(1.5);
    althphi2050->SetMarkerStyle(kFullSquare);
    althphi2050->Scale(1.0/(2.4*althphi2050->GetXaxis()->GetBinWidth(1))); //scale by delta-eta and bin width

    TH1D* althphi2050syst = (TH1D*)hPhidphi_20_50_syst->Clone("althphi2050syst");
    althphi2050syst->Scale(1.0/(2.4*althphi2050->GetXaxis()->GetBinWidth(1)));
    althphi2050syst->SetLineColor(kOrange+1);
    althphi2050syst->SetFillColor(kOrange+1);
    althphi2050syst->SetFillStyle(0);
    althphi2050syst->SetMarkerSize(0);

    TH1D* althphi5080 = (TH1D*)hPhidphi_50_100->Clone("althphi5080");
    althphi5080->SetMarkerColor(kAzure+1);
    althphi5080->SetLineColor(kAzure+1);
    althphi5080->SetMarkerSize(1.5);
    althphi5080->SetMarkerStyle(kFullCircle);
    althphi5080->Scale(1.0/(2.4*althphi5080->GetXaxis()->GetBinWidth(1))); //scale by delta-eta and bin width

    TH1D* althphi5080syst = (TH1D*)hPhidphi_50_100_syst->Clone("althphi5080syst");
    althphi5080syst->Scale(1.0/(2.4*althphi5080->GetXaxis()->GetBinWidth(1)));
    althphi5080syst->SetLineColor(kAzure+1);
    althphi5080syst->SetFillColor(kAzure+1);
    althphi5080syst->SetFillStyle(0);
    althphi5080syst->SetMarkerSize(0);

    TH1D* althh020 = (TH1D*)hhdphi_0_20->Clone("althh020");
    althh020->SetMarkerColor(kRed);
    althh020->SetMarkerSize(1.5);
    althh020->SetMarkerStyle(kFullCross);
    althh020->SetLineColor(kRed);
    althh020->Scale(1.0/(2.4*althh020->GetXaxis()->GetBinWidth(1))); //scale by delta-eta and bin width

    TH1D* althh020syst = (TH1D*)hhdphi_0_20_syst->Clone("althh020syst");
    althh020syst->Scale(1.0/(2.4*althh020->GetXaxis()->GetBinWidth(1)));
    althh020syst->SetLineColor(kRed);
    althh020syst->SetFillColor(kRed);
    althh020syst->SetFillStyle(0);
    althh020syst->SetMarkerSize(0);

    TH1D* althh2050 = (TH1D*)hhdphi_20_50->Clone("althh2050");
    althh2050->SetMarkerColor(kOrange+2);
    althh2050->SetMarkerSize(1.5);
    althh2050->SetMarkerStyle(kOpenSquare);
    althh2050->SetLineColor(kOrange+2);
    althh2050->Scale(1.0/(2.4*althh2050->GetXaxis()->GetBinWidth(1))); //scale by delta-eta and bin width

    TH1D* althh2050syst = (TH1D*)hhdphi_20_50_syst->Clone("althh2050syst");
    althh2050syst->Scale(1.0/(2.4*althh2050->GetXaxis()->GetBinWidth(1)));
    althh2050syst->SetLineColor(kOrange+2);
    althh2050syst->SetFillColor(kOrange+2);
    althh2050syst->SetFillStyle(0);
    althh2050syst->SetMarkerSize(0);

    TH1D* althh5080 = (TH1D*)hhdphi_50_100->Clone("althh5080");
    althh5080->SetMarkerColor(kAzure+1);
    althh5080->SetMarkerSize(1.5);
    althh5080->SetMarkerStyle(kOpenCircle);
    althh5080->SetLineColor(kAzure+1);
    althh5080->Scale(1.0/(2.4*althh5080->GetXaxis()->GetBinWidth(1))); //scale by delta-eta and bin width

    TH1D* althh5080syst = (TH1D*)hhdphi_50_100_syst->Clone("althh5080syst");
    althh5080syst->Scale(1.0/(2.4*althh5080->GetXaxis()->GetBinWidth(1)));
    althh5080syst->SetLineColor(kAzure+1);
    althh5080syst->SetFillColor(kAzure+1);
    althh5080syst->SetFillStyle(0);
    althh5080syst->SetMarkerSize(0);

    // Set-up histograms for the BG shaded systematic errors
    TH1D* hphiBG020hist = new TH1D("hphiBG020hist", "hphiBG020hist", 16, -TMath::Pi()/2.0, 3.0*TMath::Pi()/2.0);
    hphiBG020hist->SetFillColor(kGray+1);
    hphiBG020hist->SetFillStyle(3244);
    
    TH1D* hphiBG2050hist = (TH1D*) hphiBG020hist->Clone("hphiBG2050hist");
    TH1D* hphiBG5080hist = (TH1D*) hphiBG020hist->Clone("hphiBG5080hist");
    TH1D* hhBG020hist = (TH1D*) hphiBG020hist->Clone("hhBG020hist");
    TH1D* hhBG2050hist = (TH1D*) hphiBG020hist->Clone("hhBG2050hist");
    TH1D* hhBG5080hist = (TH1D*) hphiBG020hist->Clone("hhBG5080hist");
    for(int i = 1; i<=18; i++){
        if(IS_AVERAGE_BG){
            hphiBG020hist->SetBinContent(i, hphiBG->GetParameter(0)/(2.4*hphiBG020hist->GetBinWidth(1)));
            hphiBG2050hist->SetBinContent(i, hphiBG_20_50->GetParameter(0)/(2.4*hphiBG2050hist->GetBinWidth(1)));
            hphiBG5080hist->SetBinContent(i, hphiBG_50_100->GetParameter(0)/(2.4*hphiBG5080hist->GetBinWidth(1)));
            hhBG020hist->SetBinContent(i, hhBG->GetParameter(0)/(2.4*hhBG020hist->GetBinWidth(1)));
            hhBG2050hist->SetBinContent(i, hhBG_20_50->GetParameter(0)/(2.4*hhBG2050hist->GetBinWidth(1)));
            hhBG5080hist->SetBinContent(i, hhBG_50_100->GetParameter(0)/(2.4*hhBG5080hist->GetBinWidth(1)));
        }else{
            hphiBG020hist->SetBinContent(i, refoffset020avg*refoffset020/(2.4*hphiBG020hist->GetBinWidth(1)));
            hphiBG2050hist->SetBinContent(i, refoffset2050avg*refoffset2050/(2.4*hphiBG020hist->GetBinWidth(1)));
            hphiBG5080hist->SetBinContent(i, refoffset5080avg*refoffset5080/(2.4*hphiBG020hist->GetBinWidth(1)));
            hhBG020hist->SetBinContent(i, refhhoffset020avg*refhhoffset020/(2.4*hphiBG020hist->GetBinWidth(1)));
            hhBG2050hist->SetBinContent(i, refhhoffset2050avg*refhhoffset2050/(2.4*hphiBG020hist->GetBinWidth(1)));
            hhBG5080hist->SetBinContent(i, refhhoffset5080avg*refhhoffset5080/(2.4*hphiBG020hist->GetBinWidth(1)));
        }
       

        
        
        hhBG5080hist->SetBinError(i, hhBG5080hist->GetBinContent(i)*refhhoffset5080std);
        hhBG020hist->SetBinError(i, hhBG020hist->GetBinContent(i)*refhhoffset020std);
        hhBG2050hist->SetBinError(i, hhBG2050hist->GetBinContent(i)*refhhoffset2050std);
        hphiBG5080hist->SetBinError(i, hphiBG5080hist->GetBinContent(i)*refoffset5080std);
        hphiBG020hist->SetBinError(i, hphiBG020hist->GetBinContent(i)*refoffset020std);
        hphiBG2050hist->SetBinError(i, hphiBG2050hist->GetBinContent(i)*refoffset2050std);

}
    


    TF1* fline = new TF1("fline", "pol0", -0.5*TMath::Pi(), 1.5*TMath::Pi());
    fline->SetParameter(0, 0.0);
    fline->SetLineColor(kBlack);
    fline->SetLineWidth(3);
    fline->SetLineStyle(7);

//    TLegend* hphileg = new TLegend(0.456, 0.641, 0.613, 0.937);
    TLegend* hphileg = new TLegend(0.5219, 0.6243, 0.9653, 0.9148);
    hphileg->SetBorderSize(0);
    hphileg->AddEntry(hphi020, "V0A 0#font[122]{-}20%", "lep");
    hphileg->AddEntry(hphi2050, "V0A 20#font[122]{-}50%", "lep");
    hphileg->AddEntry(hphi5080, "V0A 50#font[122]{-}80%", "lep");

//    TLegend* hhleg = new TLegend(0.456, 0.641, 0.613, 0.937);
    TLegend* hhleg = new TLegend(0.5219, 0.6243, 0.9653, 0.9148);
    hhleg->SetBorderSize(0);
    hhleg->AddEntry(hh020, "V0A 0#font[122]{-}20%", "lep");
    hhleg->AddEntry(hh2050, "V0A 20#font[122]{-}50%", "lep");
    hhleg->AddEntry(hh5080, "V0A 50#font[122]{-}80%", "lep");
    
    TPaveText* hphitext020 = new TPaveText(0.208, 0.786, 0.383, 0.883, "NDC");
    hphitext020->SetBorderSize(0);
    hphitext020->SetFillColor(kWhite);
    hphitext020->AddText("h#font[122]{-}#phi");
    hphitext020->SetTextFont(42);
    TPaveText* hphitext2050 = new TPaveText(0.181, 0.827, 0.356, 0.925, "NDC");
    hphitext2050->SetBorderSize(0);
    hphitext2050->SetFillColor(kWhite);
    hphitext2050->AddText("h#font[122]{-}#phi");
    hphitext2050->SetTextFont(42);
    TPaveText* hphitext5080 = new TPaveText(0.181, 0.80, 0.356, 0.90, "NDC");
    hphitext5080->SetBorderSize(0);
    hphitext5080->SetFillColor(kWhite);
    hphitext5080->AddText("h#font[122]{-}#phi");
    hphitext5080->SetTextFont(42);

    TPaveText* hhtext020 = new TPaveText(0.2147, 0.7856, 0.3698, 0.8817, "NDC");
    hhtext020->SetBorderSize(0);
    hhtext020->SetFillColor(kWhite);
    hhtext020->AddText("h#font[122]{-}h");
    //hhtext020->SetTextSize(36);
    hhtext020->SetTextFont(42);
    TPaveText* hhtext2050 = new TPaveText(0.2147, 0.7856, 0.3698, 0.8817, "NDC");
    hhtext2050->SetBorderSize(0);
    hhtext2050->SetFillColor(kWhite);
    hhtext2050->AddText("h#font[122]{-}h");
    hhtext2050->SetTextFont(42);
    TPaveText* hhtext5080 = new TPaveText(0.2147, 0.7856, 0.3698, 0.8817, "NDC");
    hhtext5080->SetBorderSize(0);
    hhtext5080->SetFillColor(kWhite);
    hhtext5080->AddText("h#font[122]{-}h");
    hhtext5080->SetTextFont(42);


    TCanvas* chphi = new TCanvas("chphi", "chphi", 50, 50, 550, 600);
    chphi->cd();
    chphi->SetMargin(0.12, 0.05, 0.1, 0.05);
    hphi020->GetYaxis()->SetTitle("1/N_{trig} dN/d#Delta#varphi per #Delta#eta - constant (rad^{-1})");
    hphi020->Draw("E0 X0 P SAME");
    hphi2050->Draw("E0 X0 P SAME");
    hphi5080->Draw("E0 X0 P SAME");
    fline->Draw("SAME");
    hphitext020->Draw();
    text2->Draw();
    hphileg->Draw();

    TCanvas* chh = new TCanvas("chh", "chh", 50, 50, 550, 600);
    chh->cd();
    chh->SetMargin(0.12, 0.05, 0.1, 0.05);
    hh020->Draw("E0 X0 P SAME");
    hh2050->Draw("E0 X0 P SAME");
    hh5080->Draw("E0 X0 P SAME");
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
   
    TLegend* bgleg = new TLegend(0.555, 0.683, 0.96, 0.817);
    bgleg->SetMargin(0.4);
    bgleg->AddEntry(hphiBG, "Flat U.E.", "l");
    bgleg->AddEntry(hphiBGv2, "#it{v}_{2} Est.", "l");
    bgleg->SetLineWidth(0);

    TCanvas* calthphi020 = new TCanvas("calthphi020", "calthphi020", 50, 50, 550, 600);
    calthphi020->cd();
    calthphi020->SetMargin(0.166, 0.02, 0.108, 0.054);
    //althphi020syst->GetYaxis()->SetTitle("1/#it{N}_{trig} d#it{N}/d#Delta#varphi per #Delta#eta (rad^{-1})");
    double offset = hphiBG->GetParameter(0)/(2.4*althphi020->GetBinWidth(1));
    double highlim = offset + (althphi020->GetBinContent(althphi020->GetMaximumBin()) - offset)*2.;
    double lowlim = offset - (althphi020->GetBinContent(althphi020->GetMaximumBin()) - offset)*.75;
    althphi020syst->GetYaxis()->SetRangeUser(lowlim, highlim);
    althphi020syst->GetYaxis()->SetTitleSize(0.05);
    //althphi020syst->SetFillColor(kGray+2);
    althphi020syst->Draw("AXIS");
    hphiBG020hist->Draw("E2 SAME");
    althphi020syst->Draw("E2 SAME");
    //althphi020->SetMarkerColor(kBlack);
    //althphi020->SetLineColor(kBlack);
    althphi020->Draw("EX0 P SAME");
    //hPhidphi_0_20_bg->Draw("E0 X0 P HIST SAME"); 
    hphiBG->SetParameter(0, hphiBG->GetParameter(0)/(2.4*althphi020->GetBinWidth(1)));
    hphiBG->SetLineColor(kBlack);
    hphiBG->SetLineWidth(3);
    hphiBG->SetLineStyle(9);
    hphiBG->Draw("SAME");
    //hphiBG020hist->Draw("E2 SAME");
    hphiBGv2->SetParameter(0, hphiBGv2->GetParameter(0)/(2.4*althphi020->GetBinWidth(1)));
    hphiBGv2->SetLineColor(kBlack);
    hphiBGv2->SetLineWidth(3);
    hphiBGv2->SetLineStyle(7);
    hphiBGv2->Draw("SAME");
    data->Draw();
    highmulttext->Draw();
    //hphitext020->Draw();
    calthphi020->SaveAs(Form("%s_openerror_hphi_020_CR1.pdf", ptprefix.Data()));


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

    TCanvas* calthphi2050 = new TCanvas("calthphi2050", "calthphi2050", 50, 50, 550, 600);
    calthphi2050->cd();
    calthphi2050->SetMargin(0.166, 0.02, 0.108, 0.054);
    //althphi2050syst->GetYaxis()->SetTitle("1/#it{N}_{trig} d#it{N}/d#Delta#varphi per #Delta#eta (rad^{-1})");
    double offset2050 = hphiBG_20_50->GetParameter(0)/(2.4*althphi2050->GetBinWidth(1));
    double highlim2050 = offset2050 + (althphi2050->GetBinContent(althphi2050->GetMaximumBin()) - offset2050)*2.;
    double lowlim2050 = offset2050 - (althphi2050->GetBinContent(althphi2050->GetMaximumBin()) - offset2050)*1.;

    //althphi2050syst->GetYaxis()->SetRangeUser(lowlim2050, highlim2050);
    althphi2050syst->GetYaxis()->SetRangeUser(offset2050-(offset-lowlim), offset2050+(highlim-offset));
    althphi2050syst->GetYaxis()->SetTitleSize(0.05);
    //althphi2050syst->SetFillColor(kGray+2);
    althphi2050syst->Draw("AXIS");
    hphiBG2050hist->Draw("E2 SAME");
    althphi2050syst->Draw("E2 SAME");
    //althphi2050->SetMarkerColor(kBlack);
    //althphi2050->SetLineColor(kBlack);
    althphi2050->Draw("EX0 P SAME");
    //hPhidphi_20_50_bg->Draw("E0 X0 P HIST SAME"); 
    hphiBG_20_50->SetParameter(0, hphiBG_20_50->GetParameter(0)/(2.4*althphi2050->GetBinWidth(1)));
    hphiBG_20_50->SetLineColor(kBlack);
    hphiBG_20_50->SetLineWidth(3);
    hphiBG_20_50->SetLineStyle(9);
    hphiBG_20_50->Draw("SAME");
    //hphiBG2050hist->Draw("E2 SAME");
    hphiBGv2_20_50->SetParameter(0, hphiBGv2_20_50->GetParameter(0)/(2.4*althphi2050->GetBinWidth(1)));
    hphiBGv2_20_50->SetLineColor(kBlack);
    hphiBGv2_20_50->SetLineWidth(3);
    hphiBGv2_20_50->SetLineStyle(7);
    hphiBGv2_20_50->Draw("SAME");
    //data->Draw();
    //hphitext2050->Draw();
    //text2->Draw();
    bgleg->Draw("SAME");
    midmulttext->Draw();
    calthphi2050->SaveAs(Form("%s_openerror_hphi_2050_CR1.pdf", ptprefix.Data()));


    TCanvas* chphi5080 = new TCanvas("chphi5080", "chphi5080", 50, 50, 550, 600);
    chphi5080->cd();
    chphi5080->SetMargin(0.12, 0.05, 0.1, 0.05);
    hphi5080syst->GetYaxis()->SetTitle("1/#it{N}_{trig} d#it{N}/d#Delta#varphi per #Delta#eta - constant (rad^{-1})");
    hphi5080syst->GetYaxis()->SetRangeUser(-0.6E-3, 2.8E-3);
    hphi5080syst->Draw("E2");
    hphi5080->Draw("E0 X0 P SAME");
    hphileg->Draw();
    fline->Draw("SAME");
    data->Draw();
    hphitext5080->Draw();

    TCanvas* calthphi5080 = new TCanvas("calthphi5080", "calthphi5080", 50, 50, 550, 600);
    calthphi5080->cd();
    calthphi5080->SetMargin(0.166, 0.02, 0.108, 0.054);
    //althphi5080syst->GetYaxis()->SetTitle("1/#it{N}_{trig} d#it{N}/d#Delta#varphi per #Delta#eta (rad^{-1})");

    double offset5080 = hphiBG_50_100->GetParameter(0)/(2.4*althphi5080->GetBinWidth(1));
    double highlim5080 = offset5080 + (althphi5080->GetBinContent(althphi5080->GetMaximumBin()) - offset5080)*2.;
    double lowlim5080 = offset5080 - (althphi5080->GetBinContent(althphi5080->GetMaximumBin()) - offset5080)*1.;
    //althphi5080syst->GetYaxis()->SetRangeUser(lowlim5080, highlim5080);
    althphi5080syst->GetYaxis()->SetRangeUser(offset5080-(offset-lowlim), offset5080+(highlim-offset));
    althphi5080syst->GetYaxis()->SetTitleSize(0.05);
    //althphi5080syst->SetFillColor(kGray+2);
    althphi5080syst->Draw("AXIS");
    hphiBG5080hist->Draw("E2 SAME");
    althphi5080syst->Draw("E2 SAME");
    //althphi5080->SetMarkerColor(kBlack);
    //althphi5080->SetLineColor(kBlack);
    althphi5080->Draw("EX0 P SAME");
    //hPhidphi_50_100_bg->Draw("E0 X0 P HIST SAME"); 
    hphiBG_50_100->SetParameter(0, hphiBG_50_100->GetParameter(0)/(2.4*althphi5080->GetBinWidth(1)));
    hphiBG_50_100->SetLineColor(kBlack);
    hphiBG_50_100->SetLineWidth(3);
    hphiBG_50_100->SetLineStyle(9);
    hphiBG_50_100->Draw("SAME");
    //hphiBG5080hist->Draw("E2 SAME");
    hphiBGv2_50_100->SetParameter(0, hphiBGv2_50_100->GetParameter(0)/(2.4*althphi5080->GetBinWidth(1)));
    hphiBGv2_50_100->SetLineColor(kBlack);
    hphiBGv2_50_100->SetLineWidth(3);
    hphiBGv2_50_100->SetLineStyle(7);
    hphiBGv2_50_100->Draw("SAME");
    //data->Draw();
    hphitext5080->Draw();
    //hphileg->Draw();
    //bgleg->Draw();
    text2->Draw();
    lowmulttext->Draw();
    calthphi5080->SaveAs(Form("%s_openerror_hphi_5080_CR1.pdf", ptprefix.Data()));

    TCanvas* calthphiall = new TCanvas("calthphiall", "calthphiall", 50, 50, 550, 600);
    calthphiall->cd();
    calthphiall->SetMargin(0.12, 0.05, 0.1, 0.05);
    althphi020syst->Draw("E2");
    althphi020->Draw("E0 X0 P SAME");
    althphi2050syst->Draw("E2 SAME");
    althphi2050->Draw("E0 X0 P SAME");
    althphi5080syst->Draw("E2 SAME");
    althphi5080->Draw("E0 X0 P SAME");
    hphiBG->Draw("SAME");
    hphiBG_20_50->Draw("SAME");
    hphiBG_50_100->Draw("SAME");
    data->Draw();
    hphitext5080->Draw();



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
    
    TCanvas* calthh020 = new TCanvas("calthh020", "calthh020", 50, 50, 550, 600);
    calthh020->cd();
    calthh020->SetMargin(0.166, 0.02, 0.108, 0.054);
    //althh020syst->GetYaxis()->SetTitle("1/#it{N}_{trig} d#it{N}/d#Delta#varphi per #Delta#eta (rad^{-1})");
    //althh020syst->SetFillColor(kGray+3);
    althh020syst->GetYaxis()->SetTitleSize(0.05);
    althh020syst->Draw("E2");
    //althh020->SetMarkerColor(kBlack);
    //althh020->SetLineColor(kBlack);
    //althh020->SetMarkerStyle(22);
    althh020->Draw("EX0 P SAME");
    //hhdphi_0_20_bg->Draw("E0 X0 P HIST SAME");
    hhBG->SetParameter(0, hhBG->GetParameter(0)/(2.4*althh020->GetBinWidth(1)));
    double hhoffset = hhBG->GetParameter(0);
    double hhhighlim = hhoffset + (althh020->GetBinContent(althh020->GetMaximumBin()) - hhoffset)*2.;
    double hhlowlim = hhoffset - (althh020->GetBinContent(althh020->GetMaximumBin()) - hhoffset)*.75;
    althh020syst->GetYaxis()->SetRangeUser(hhlowlim, hhhighlim);
    hhBG->SetLineColor(kBlack);
    hhBG->SetLineWidth(4);
    hhBG->SetLineStyle(9);
    hhBG->Draw("SAME");
    //hhBG020hist->Draw("E2 SAME");
    hhBGv2->SetParameter(0, hhBGv2->GetParameter(0)/(2.4*althh020->GetBinWidth(1)));
    hhBGv2->SetLineColor(kBlack);
    hhBGv2->SetLineWidth(3);
    hhBGv2->SetLineStyle(7);
    hhBGv2->Draw("SAME");
    data->Draw();
    highmulttext->Draw();
    //hhtext020->Draw();
    calthh020->SaveAs(Form("%s_openerror_hh_020_CR1.pdf", ptprefix.Data()));

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

    TCanvas* calthh2050 = new TCanvas("calthh2050", "calthh2050", 50, 50, 550, 600);
    calthh2050->cd();
    calthh2050->SetMargin(0.166, 0.02, 0.108, 0.054);
    //althh2050syst->GetYaxis()->SetTitle("1/#it{N}_{trig} d#it{N}/d#Delta#varphi per #Delta#eta (rad^{-1})");
    //althh2050syst->SetFillColor(kGray+3);
    althh2050syst->GetYaxis()->SetTitleSize(0.05);
    althh2050syst->Draw("E2");
    //althh2050->SetMarkerColor(kBlack);
    //althh2050->SetLineColor(kBlack);
    //althh2050->SetMarkerStyle(22);
    althh2050->Draw("EX0 P SAME");
    //hhdphi_20_50_bg->Draw("E0 X0 P HIST SAME");
    hhBG_20_50->SetParameter(0, hhBG_20_50->GetParameter(0)/(2.4*althphi2050->GetBinWidth(1)));
    double hhoffset2050 = hhBG_20_50->GetParameter(0);
    althh2050syst->GetYaxis()->SetRangeUser(hhoffset2050-(hhoffset - hhlowlim), hhoffset2050 + (hhhighlim - hhoffset));
    hhBG_20_50->SetLineColor(kBlack);
    hhBG_20_50->SetLineWidth(4);
    hhBG_20_50->SetLineStyle(9);
    hhBG_20_50->Draw("SAME");
    //hhBG2050hist->Draw("E2 SAME");
    hhBGv2_20_50->SetParameter(0, hhBGv2_20_50->GetParameter(0)/(2.4*althh2050->GetBinWidth(1)));
    hhBGv2_20_50->SetLineColor(kBlack);
    hhBGv2_20_50->SetLineWidth(3);
    hhBGv2_20_50->SetLineStyle(7);
    hhBGv2_20_50->Draw("SAME");
    //data->Draw();
    //hhtext2050->Draw();
    //text2hh->Draw();
    bgleg->Draw("SAME");
    midmulttext->Draw();
    calthh2050->SaveAs(Form("%s_openerror_hh_2050_CR1.pdf", ptprefix.Data()));
    
    TCanvas* chh5080 = new TCanvas("chh5080", "chh5080", 50, 50, 550, 600);
    chh5080->cd();
    chh5080->SetMargin(0.12, 0.05, 0.1, 0.05);
    hh5080syst->GetYaxis()->SetTitle("1/N_{trig} dN/d#Delta#varphi per #Delta#eta - constant (rad^{-1})");
    hh5080syst->GetYaxis()->SetRangeUser(-50E-3, 280E-3);
    hh5080syst->Draw("E2");
    hh5080->Draw("E0 X0 P SAME");
    hhleg->Draw();
    fline->Draw("SAME");
    data->Draw();
    hhtext5080->Draw();

    TCanvas* calthh5080 = new TCanvas("calthh5080", "calthh5080", 50, 50, 550, 600);
    calthh5080->cd();
    calthh5080->SetMargin(0.166, 0.02, 0.108, 0.054);
    //althh5080syst->GetYaxis()->SetTitle("1/#it{N}_{trig} d#it{N}/d#Delta#varphi per #Delta#eta (rad^{-1})");
    //althh5080syst->SetFillColor(kGray+3);
    althh5080syst->GetYaxis()->SetTitleSize(0.05);
    althh5080syst->Draw("E2");
    //althh5080->SetMarkerColor(kBlack);
    //althh5080->SetLineColor(kBlack);
    //althh5080->SetMarkerStyle(22);
    althh5080->Draw("EX0 P SAME");
    //hhdphi_50_100_bg->Draw("E0 X0 P HIST SAME");
    hhBG_50_100->SetParameter(0, hhBG_50_100->GetParameter(0)/(2.4*althh5080->GetBinWidth(1)));
    double hhoffset5080 = hhBG_50_100->GetParameter(0);
    althh5080syst->GetYaxis()->SetRangeUser(hhoffset5080 - (hhoffset - hhlowlim), hhoffset5080 + (hhhighlim - hhoffset));
    hhBG_50_100->SetLineColor(kBlack);
    hhBG_50_100->SetLineWidth(4);
    hhBG_50_100->SetLineStyle(9);
    hhBG_50_100->Draw("SAME");
    //hhBG5080hist->Draw("E2 SAME");
    hhBGv2_50_100->SetParameter(0, hhBGv2_50_100->GetParameter(0)/(2.4*althh5080->GetBinWidth(1)));
    hhBGv2_50_100->SetLineColor(kBlack);
    hhBGv2_50_100->SetLineWidth(3);
    hhBGv2_50_100->SetLineStyle(7);
    hhBGv2_50_100->Draw("SAME");
    //data->Draw();
    hhtext5080->Draw();
    //hhleg->Draw();
    //bgleg->Draw();
    text2hh->Draw();
    lowmulttext->Draw();
    calthh5080->SaveAs(Form("%s_openerror_hh_5080_CR1.pdf", ptprefix.Data()));


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
    hphi2050syst->Draw("E2");
    hphi2050->Draw("E0 X0 P SAME");
    fline->Draw("SAME");
    text2->Draw();
    //data->Draw();
    //hphitext2050->Draw();
    chphiall->cd();
    TPad* hphi5080pad = new TPad("hphi5080pad", "", 0.68, 0, 1.0, 1.0);
    hphi5080pad->SetMargin(0.0, 0.05, 0.1, 0.1);
    hphi5080pad->Draw();
    hphi5080pad->cd();
    hphi5080syst->GetXaxis()->SetLabelSize(0.05);
    hphi5080syst->GetYaxis()->SetLabelSize(0.0);
    hphi5080syst->Draw("E2");
    hphi5080->Draw("E0 X0 P SAME");
    hphileg->Draw();
    fline->Draw("SAME");
    //data->Draw();
    //hphitext5080->Draw();

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
    text2hh->Draw();
    //hhtext2050->Draw();
    //chhall->cd(3)->SetMargin(0.0, 0.1, 0.1, 0.1);
    chhall->cd();
    TPad* hh5080pad = new TPad("hh5080pad", "", 0.68, 0, 1.0, 1.0);
    hh5080pad->SetMargin(0.0, 0.05, 0.1, 0.1);
    hh5080pad->Draw();
    hh5080pad->cd();
    hh5080syst->GetXaxis()->SetLabelSize(0.05);
    hh5080syst->GetYaxis()->SetLabelSize(0.0);
    hh5080syst->Draw("E2");
    hh5080->Draw("E0 X0 P SAME");
    hhleg->Draw();
    fline->Draw("SAME");
    //data->Draw();
    //hhtext5080->Draw();


    // Histograms to calc the sytematic error on the UE estimation method
    TH1F* hoff020 = new TH1F("hoff020", "hoff020", 1000, 0.5, 1.5);
    TH1F* hoff2050 = new TH1F("hoff2050", "hoff2050", 1000, 0.5, 1.5);
    TH1F* hoff5080 = new TH1F("hoff5080", "hoff5080", 1000, 0.5, 1.5);
    TH1F* hhhoff020 = new TH1F("hhhoff020", "hhhoff020", 1000, 0.5, 1.5);
    TH1F* hhhoff2050 = new TH1F("hhhoff2050", "hhhoff2050", 1000, 0.5, 1.5);
    TH1F* hhhoff5080 = new TH1F("hhhoff5080", "hhhoff5080", 1000, 0.5, 1.5);

    hoff020->Fill(offset/refoffset020);
    hoff2050->Fill(offset2050/refoffset2050);
    hoff5080->Fill(offset5080/refoffset5080);
    hhhoff020->Fill(hhoffset/refhhoffset020);
    hhhoff2050->Fill(hhoffset2050/refhhoffset2050);
    hhhoff5080->Fill(hhoffset5080/refhhoffset5080);




    /*
    TCanvas *cratio = new TCanvas("cratio", "cratio", 50, 50, 550, 600);
    cratio->cd();
    cratio->SetMargin(0.12, 0.05, 0.1, 0.05);
    ratio->Draw("E0 X0");
*/

    //Do plots of the yields for the systematics:
/*    TH1F* hhyieldSyst_0_20 = new TH1F("hhyieldSyst_0_20", "hhyieldSyst_0_20", 1000, 0.0, 10.0);
    TH1F* hhyieldSyst_20_50 = new TH1F("hhyieldSyst_20_50", "hhyieldSyst_20_50", 1000, 0.0, 10.0);
    TH1F* hhyieldSyst_50_100 = new TH1F("hhyieldSyst_50_100", "hhyieldSyst_50_100", 1000, 0.0, 10.0);
    TH1F* hphiyieldSyst_0_20 = new TH1F("hphiyieldSyst_0_20", "hphiyieldSyst_0_20", 700, 3.0E-02, 10.0E-02);
    TH1F* hphiyieldSyst_20_50 = new TH1F("hphiyieldSyst_20_50", "hphiyieldSyst_20_50", 600, 1.0E-02, 8.0E-02);
    TH1F* hphiyieldSyst_50_100 = new TH1F("hphiyieldSyst_50_100", "hphiyieldSyst_50_100", 700, 0.5E-02, 7.5E-02);
    
    TH1F* hhnearyieldSyst_0_20 = new TH1F("hhnearyieldSyst_0_20", "hhnearyieldSyst_0_20", 1000, 0.0, 2.0);
    TH1F* hhnearyieldSyst_20_50 = new TH1F("hhnearyieldSyst_20_50", "hhnearyieldSyst_20_50", 1000, 0.0, 2.0);
    TH1F* hhnearyieldSyst_50_100 = new TH1F("hhnearyieldSyst_50_100", "hhnearyieldSyst_50_100", 1000, 0.0, 2.0);
    TH1F* hphinearyieldSyst_0_20 = new TH1F("hphinearyieldSyst_0_20", "hphinearyieldSyst_0_20", 800, 0.5E-03, 4.5E-03);
    TH1F* hphinearyieldSyst_20_50 = new TH1F("hphinearyieldSyst_20_50", "hphinearyieldSyst_20_50", 800, 0.5E-03, 4.5E-03);
    TH1F* hphinearyieldSyst_50_100 = new TH1F("hphinearyieldSyst_50_100", "hphinearyieldSyst_50_100", 800, 0.5E-03, 4.5E-03);
   
    TH1F* hhawayyieldSyst_0_20 = new TH1F("hhawayyieldSyst_0_20", "hhawayyieldSyst_0_20", 1000, 0.0, 2.0);
    TH1F* hhawayyieldSyst_20_50 = new TH1F("hhawayyieldSyst_20_50", "hhawayyieldSyst_20_50", 1000, 0.0, 2.0);
    TH1F* hhawayyieldSyst_50_100 = new TH1F("hhawayyieldSyst_50_100", "hhawayyieldSyst_50_100", 1000, 0.0, 2.0);
    TH1F* hphiawayyieldSyst_0_20 = new TH1F("hphiawayyieldSyst_0_20", "hphiawayyieldSyst_0_20", 800, 0.5E-03, 4.5E-03);
    TH1F* hphiawayyieldSyst_20_50 = new TH1F("hphiawayyieldSyst_20_50", "hphiawayyieldSyst_20_50", 800, 0.5E-03, 4.5E-03);
    TH1F* hphiawayyieldSyst_50_100 = new TH1F("hphiawayyieldSyst_50_100", "hphiawayyieldSyst_50_100", 800, 0.5E-03, 4.5E-03);

    TH1F* hhbulkyieldSyst_0_20 = new TH1F("hhbulkyieldSyst_0_20", "hhbulkyieldSyst_0_20", 1000, 0.0, 10.0);
    TH1F* hhbulkyieldSyst_20_50 = new TH1F("hhbulkyieldSyst_20_50", "hhbulkyieldSyst_20_50", 1000, 0.0, 10.0);
    TH1F* hhbulkyieldSyst_50_100 = new TH1F("hhbulkyieldSyst_50_100", "hhbulkyieldSyst_50_100", 1000, 0.0, 10.0);
    TH1F* hphibulkyieldSyst_0_20 = new TH1F("hphibulkyieldSyst_0_20", "hphibulkyieldSyst_0_20", 700, 3.0E-02, 10.0E-02);
    TH1F* hphibulkyieldSyst_20_50 = new TH1F("hphibulkyieldSyst_20_50", "hphibulkyieldSyst_20_50", 600, 1.0E-02, 8.0E-02);
    TH1F* hphibulkyieldSyst_50_100 = new TH1F("hphibulkyieldSyst_50_100", "hphibulkyieldSyst_50_100", 700, 0.5E-02, 7.5E-02);
*/

    TH1F* hhyieldSyst_0_20 = new TH1F("hhyieldSyst_0_20", "hhyieldSyst_0_20", 500, 0.75, 1.25);
    TH1F* hhyieldSyst_20_50 = new TH1F("hhyieldSyst_20_50", "hhyieldSyst_20_50", 500, 0.75, 1.25);
    TH1F* hhyieldSyst_50_100 = new TH1F("hhyieldSyst_50_100", "hhyieldSyst_50_100", 500, 0.75, 1.25);
    TH1F* hphiyieldSyst_0_20 = new TH1F("hphiyieldSyst_0_20", "hphiyieldSyst_0_20", 500, 0.75, 1.25);
    TH1F* hphiyieldSyst_20_50 = new TH1F("hphiyieldSyst_20_50", "hphiyieldSyst_20_50", 500, 0.75, 1.25);
    TH1F* hphiyieldSyst_50_100 = new TH1F("hphiyieldSyst_50_100", "hphiyieldSyst_50_100", 500, 0.75, 1.25);
    
    TH1F* hhnearyieldSyst_0_20 = new TH1F("hhnearyieldSyst_0_20", "hhnearyieldSyst_0_20", 500, 0.75, 1.25);
    TH1F* hhnearyieldSyst_20_50 = new TH1F("hhnearyieldSyst_20_50", "hhnearyieldSyst_20_50", 500, 0.75, 1.25);
    TH1F* hhnearyieldSyst_50_100 = new TH1F("hhnearyieldSyst_50_100", "hhnearyieldSyst_50_100", 500, 0.75, 1.25);
    TH1F* hphinearyieldSyst_0_20 = new TH1F("hphinearyieldSyst_0_20", "hphinearyieldSyst_0_20", 500, 0.75, 1.25);
    TH1F* hphinearyieldSyst_20_50 = new TH1F("hphinearyieldSyst_20_50", "hphinearyieldSyst_20_50", 500, 0.75, 1.25);
    TH1F* hphinearyieldSyst_50_100 = new TH1F("hphinearyieldSyst_50_100", "hphinearyieldSyst_50_100", 500, 0.75, 1.25);
   
    TH1F* hhawayyieldSyst_0_20 = new TH1F("hhawayyieldSyst_0_20", "hhawayyieldSyst_0_20", 500, 0.75, 1.25);
    TH1F* hhawayyieldSyst_20_50 = new TH1F("hhawayyieldSyst_20_50", "hhawayyieldSyst_20_50", 500, 0.75, 1.25);
    TH1F* hhawayyieldSyst_50_100 = new TH1F("hhawayyieldSyst_50_100", "hhawayyieldSyst_50_100", 500, 0.75, 1.25);
    TH1F* hphiawayyieldSyst_0_20 = new TH1F("hphiawayyieldSyst_0_20", "hphiawayyieldSyst_0_20", 500, 0.75, 1.25);
    TH1F* hphiawayyieldSyst_20_50 = new TH1F("hphiawayyieldSyst_20_50", "hphiawayyieldSyst_20_50", 500, 0.75, 1.25);
    TH1F* hphiawayyieldSyst_50_100 = new TH1F("hphiawayyieldSyst_50_100", "hphiawayyieldSyst_50_100", 500, 0.75, 1.25);

    TH1F* hhbulkyieldSyst_0_20 = new TH1F("hhbulkyieldSyst_0_20", "hhbulkyieldSyst_0_20", 500, 0.75, 1.25);
    TH1F* hhbulkyieldSyst_20_50 = new TH1F("hhbulkyieldSyst_20_50", "hhbulkyieldSyst_20_50", 500, 0.75, 1.25);
    TH1F* hhbulkyieldSyst_50_100 = new TH1F("hhbulkyieldSyst_50_100", "hhbulkyieldSyst_50_100", 500, 0.75, 1.25);
    TH1F* hphibulkyieldSyst_0_20 = new TH1F("hphibulkyieldSyst_0_20", "hphibulkyieldSyst_0_20", 500, 0.75, 1.25);
    TH1F* hphibulkyieldSyst_20_50 = new TH1F("hphibulkyieldSyst_20_50", "hphibulkyieldSyst_20_50", 500, 0.75, 1.25);
    TH1F* hphibulkyieldSyst_50_100 = new TH1F("hphibulkyieldSyst_50_100", "hphibulkyieldSyst_50_100", 500, 0.75, 1.25);

    Float_t stdtotal0_20 = 5.850084E-02;
    Float_t stdtotal20_50 = 3.489754E-02;
    Float_t stdtotal50_100 = 1.881542E-02;
    Float_t stdnear0_20 = 1.509740E-03;
    Float_t stdnear20_50 = 1.450373E-03;
    Float_t stdnear50_100 = 8.015203E-04;
    Float_t stdaway0_20 = 2.234562E-03;
    Float_t stdaway20_50 = 2.274811E-03;
    Float_t stdaway50_100 = 1.119789E-03;
    Float_t stdmid0_20 = 5.476692E-02;
    Float_t stdmid20_50 = 3.111197E-02;
    Float_t stdmid50_100 = 1.686431E-02;

    
    
    hhyieldSyst_0_20->Fill(total0_20hhYield);
    hhyieldSyst_20_50->Fill(total20_50hhYield);
    hhyieldSyst_50_100->Fill(total50_100hhYield);
    hhnearyieldSyst_0_20->Fill(near0_20hhYield);
    hhnearyieldSyst_20_50->Fill(near20_50hhYield);
    hhnearyieldSyst_50_100->Fill(near50_100hhYield);
    hhawayyieldSyst_0_20->Fill(away0_20hhYield);
    hhawayyieldSyst_20_50->Fill(away20_50hhYield);
    hhawayyieldSyst_50_100->Fill(away50_100hhYield);
    hhbulkyieldSyst_0_20->Fill(mid0_20hhYield);
    hhbulkyieldSyst_20_50->Fill(mid20_50hhYield);
    hhbulkyieldSyst_50_100->Fill(mid50_100hhYield);

    hphiyieldSyst_0_20->Fill(total0_20hPhiYield/stdtotal0_20);
    hphiyieldSyst_20_50->Fill(total20_50hPhiYield/stdtotal20_50);
    hphiyieldSyst_50_100->Fill(total50_100hPhiYield/stdtotal50_100);
    hphinearyieldSyst_0_20->Fill(near0_20hPhiYield/stdnear0_20);
    hphinearyieldSyst_20_50->Fill(near20_50hPhiYield/stdnear20_50);
    hphinearyieldSyst_50_100->Fill(near50_100hPhiYield/stdnear50_100);
    hphiawayyieldSyst_0_20->Fill(away0_20hPhiYield/stdaway0_20);
    hphiawayyieldSyst_20_50->Fill(away20_50hPhiYield/stdaway20_50);
    hphiawayyieldSyst_50_100->Fill(away50_100hPhiYield/stdaway50_100);
    hphibulkyieldSyst_0_20->Fill(mid0_20hPhiYield/stdmid0_20);
    hphibulkyieldSyst_20_50->Fill(mid20_50hPhiYield/stdmid20_50);
    hphibulkyieldSyst_50_100->Fill(mid50_100hPhiYield/stdmid50_100);



/*
    TCanvas *csyst = new TCanvas("csyst", "csyst", 50, 50, 1000, 500);
    csyst->Divide(1,3);
    csyst->cd(1);
    hhyieldSyst_0_20->Draw();
    csyst->cd(2);
    hhyieldSyst_20_50->Draw();
    csyst->cd(3);
    hhyieldSyst_50_100->Draw();
*/
    TFile* output = new TFile(Form("fitsyst_%s.root", outputstring.Data()), "RECREATE");

    hhyieldSyst_0_20->Write();
    hhyieldSyst_20_50->Write();
    hhyieldSyst_50_100->Write();
    hhnearyieldSyst_0_20->Write();
    hhnearyieldSyst_20_50->Write();
    hhnearyieldSyst_50_100->Write();
    hhawayyieldSyst_0_20->Write();
    hhawayyieldSyst_20_50->Write();
    hhawayyieldSyst_50_100->Write();
    hhbulkyieldSyst_0_20->Write();
    hhbulkyieldSyst_20_50->Write();
    hhbulkyieldSyst_50_100->Write();

    hhdphi_0_20->Write();
    hhdphi_0_20_syst->Write();
    hPhidphi_0_20->Write();
    hPhidphi_0_20_syst->Write();
    hhdphi_20_50->Write();
    hhdphi_20_50_syst->Write();
    hPhidphi_20_50->Write();
    hPhidphi_20_50_syst->Write();
    hhdphi_50_100->Write();
    hhdphi_50_100_syst->Write();
    hPhidphi_50_100->Write();
    hPhidphi_50_100_syst->Write();
    

    hphiyieldSyst_0_20->Write();
    hphiyieldSyst_20_50->Write();
    hphiyieldSyst_50_100->Write();
    hphinearyieldSyst_0_20->Write();
    hphinearyieldSyst_20_50->Write();
    hphinearyieldSyst_50_100->Write();
    hphiawayyieldSyst_0_20->Write();
    hphiawayyieldSyst_20_50->Write();
    hphiawayyieldSyst_50_100->Write();
    hphibulkyieldSyst_0_20->Write();
    hphibulkyieldSyst_20_50->Write();
    hphibulkyieldSyst_50_100->Write();

    ratioNearHist->Write();
    ratiosNear->Write("ratiosNear");
    ratiosNearSyst->Write("ratiosNearSyst");
    ratiosAway->Write("ratiosAway");
    ratiosAwaySyst->Write("ratiosAwaySyst");
    ratiosBulk->Write("ratiosBulk");
    ratiosBulkSyst->Write("ratiosBulkSyst");
    ratiosTot->Write("ratiosTot");
    ratiosTotSyst->Write("ratiosTotSyst");
    jetratioshh->Write("jetratioshh");
    jetratioshPhi->Write("jetratioshPhi");

    hoff020->Write();
    hoff2050->Write();
    hoff5080->Write();
    hhhoff020->Write();
    hhhoff2050->Write();
    hhhoff5080->Write();

    yieldshhNear->SetName("yieldshhNear");
    yieldshhNear->Write();
    yieldshhAway->SetName("yieldshhAway");
    yieldshhAway->Write();

    vsMultCanvas->Write("ratiocanvas");

    printf("offsets: \n       h-phi     h-h\n 0-20   %e   %e\n 20-50   %e   %e\n 50-80   %e   %e\n", offset, hhoffset, offset2050, hhoffset2050, offset5080, hhoffset5080);

    printf("high: %.2e $\\pm$ %.2e ($\\pm$ %.2e)   %.2e $\\pm$ %.2e ($\\pm$ %.2e)  %.2e $\\pm$ %.2e ($\\pm$ %.2e)  %.2e $\\pm$ %.2e ($\\pm$ %.2e)\nmid: %.2e $\\pm$ %.2e ($\\pm$ %.2e)   %.2e $\\pm$ %.2e ($\\pm$ %.2e)   %.2e $\\pm$ %.2e ($\\pm$ %.2e)   %.2e $\\pm$ %.2e ($\\pm$ %.2e)\nlow: %.2e $\\pm$ %.2e ($\\pm$ %.2e)   %.2e $\\pm$ %.2e ($\\pm$ %.2e)   %.2e $\\pm$ %.2e ($\\pm$ %.2e)   %.2e $\\pm$ %.2e ($\\pm$ %.2e)\n", nearArray[2], nearArrayErr[2], nearArraySystErr[2], awayArray[2], awayArrayErr[2], awayArraySystErr[2], bulkArray[2], bulkArrayErr[2], bulkArraySystErr[2], totalArray[2], totalArrayErr[2], totalArraySystErr[2], nearArray[1], nearArrayErr[1], nearArraySystErr[1], awayArray[1], awayArrayErr[1], awayArraySystErr[1], bulkArray[1], bulkArrayErr[1], bulkArraySystErr[1], totalArray[1], totalArrayErr[1], totalArraySystErr[1], nearArray[0], nearArrayErr[0], nearArraySystErr[0], awayArray[0], awayArrayErr[0], awayArraySystErr[0], bulkArray[0], bulkArrayErr[0], bulkArraySystErr[0], totalArray[0], totalArrayErr[0], totalArraySystErr[0]);

}
