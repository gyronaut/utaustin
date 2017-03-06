#include <string>

void plotCalibration(string filename){
    
    /*setting up some global style variables */ 
    gStyle->SetLabelSize(0.05, "xyz");
    gStyle->SetLabelOffset(0.015, "xyz");
    gStyle->SetTitleSize(0.05, "xyz");
    gStyle->SetTitleSize(0.055, "h");
    gStyle->SetTitleOffset(1.1, "x");
    gStyle->SetTitleOffset(1.2, "y");
    gStyle->SetPadTopMargin(0.08);
    gStyle->SetPadBottomMargin(0.12);
    gStyle->SetPadLeftMargin(0.14);
    gStyle->SetPadRightMargin(0.1);
    gStyle->SetOptStat(0);
    gROOT->ForceStyle();


    TFile *file = new TFile(filename.c_str());
    
    if(!file){
        printf("file not loaded!\n");
        return;
    }

    TH1F *pass3[4];
    TH1F *calib[4];
    TH2F *rawTime[4];
    TH1D *rawProjection[4];
    TH2F *corrTime[4];
    TH1D *corrProjection[4];

    TH1F *compare[4];


    //Initialize OADB container
    AliOADBContainer *container = new AliOADBContainer("");
    container->InitFromFile("$ALICE_PHYSICS/OADB/EMCAL/EMCALTimeCalib.root", "AliEMCALTimeCalib");
    if(!container){
        fprintf(stderr, "No OADB container!");
        return;
    }

    //Find container in OADB
    TObject *oadbCalib = (TObject*) container->GetObject(0, "TimeCalib10d");
    if(!oadbCalib){
        fprintf(stderr, "No calibration container found in OADB!");
        return;
    }

    //Get Array for pass1 (13g) and pass4 (13b and 13c) and set up histograms
    TObject *oadbCalibPass3 = (TObject*) oadbCalib->FindObject("pass3");
    if(!oadbCalibPass3){
        fprintf(stderr, "No pass3!");
        return;
    }

    //Set-up all histograms
    for(int ibc=0; ibc<4; ibc++){
        pass3[ibc] = (TH1F*)oadbCalibPass3->FindObject(Form("hAllTimeAvBC%d", ibc));
        if(!pass3[ibc]){
            printf("no pass3 histogram!\n");
            return;
        }
        pass3[ibc]->GetXaxis()->SetTitle("Cell AbsID");
        pass3[ibc]->GetXaxis()->SetRangeUser(0, 4607);
        pass3[ibc]->GetYaxis()->SetTitle("Calibration Offset (ns)");
        pass3[ibc]->GetYaxis()->SetRangeUser(575, 675);

        calib[ibc] = (TH1F*)file->Get(Form("hAllTimeAvLGBC%d", ibc));
        calib[ibc]->GetXaxis()->SetTitle("Cell AbsID");
        calib[ibc]->GetXaxis()->SetRangeUser(0, 4607);
        calib[ibc]->GetYaxis()->SetTitle("Calibration Offset (ns)");
        calib[ibc]->GetYaxis()->SetRangeUser(575, 675);
       
        rawTime[ibc] = (TH2F*)file->Get(Form("RawTimeVsIdLGBC%d", ibc));
        rawTime[ibc]->GetXaxis()->SetTitle("Cell AbsID");
        rawTime[ibc]->GetXaxis()->SetRangeUser(0, 4607);
        rawTime[ibc]->GetYaxis()->SetTitle("Time (ns)");
        rawTime[ibc]->GetYaxis()->SetRangeUser(450, 800);
        rawTime[ibc]->GetZaxis()->SetRangeUser(1, 150);
        rawTime[ibc]->SetTitle(Form("Raw Time vs. Cell ID for LHC10d pass4, BC%d", ibc));

        rawProjection[ibc] = rawTime[ibc]->ProjectionY(Form("rawprojBC%d", ibc), 0, 4607);

        corrTime[ibc] = (TH2F*)file->Get(Form("CorrectedTimeBC%d", ibc));
        corrTime[ibc]->GetXaxis()->SetTitle("Cell AbsID");
        corrTime[ibc]->GetXaxis()->SetRangeUser(0, 4607);
        corrTime[ibc]->GetYaxis()->SetTitle("Corrected Time (ns)");
        corrTime[ibc]->GetYaxis()->SetRangeUser(-150, 150);
        corrTime[ibc]->GetZaxis()->SetRangeUser(1, 150);
        corrTime[ibc]->SetTitle(Form("Corrected Time for LHC10d pass4, BC%d", ibc));

        corrProjection[ibc] = corrTime[ibc]->ProjectionY(Form("corrprojBC%d", ibc), 0, 4607);

        compare[ibc] =(TH1F*)calib[ibc]->Clone(Form("compareOADBBC%d", ibc));        
        for(int ibin=0; ibin < compare[ibc]->GetNbinsX(); ibin++){
            compare[ibc]->SetBinContent(ibin, pass3[ibc]->GetBinContent(ibin) - compare[ibc]->GetBinContent(ibin));
        }
        compare[ibc]->GetXaxis()->SetTitle("Cell AbsID");
        compare[ibc]->GetXaxis()->SetRangeUser(0, 4607);
        compare[ibc]->GetYaxis()->SetTitle("Difference in Calibration (ns)");
        compare[ibc]->GetYaxis()->SetRangeUser(-3, 3);
        compare[ibc]->SetTitle(Form("Difference (Old OADB Pass3 - New Calibration Pass4) for LHC10d, BC%d", ibc));

    }


    int bcNum=0;
    
    //raw time
    TCanvas* cRaw = new TCanvas("rawtime", "rawtime", 1000,600);
    cRaw->cd();
    gPad->SetRightMargin(0.15);
    gPad->SetLogz();
    rawTime[bcNum]->Draw("SAME COLZ");

    //raw time projection
    TCanvas* cRawProj =  new TCanvas("rawproj", "rawproj", 1000, 300);
    cRawProj->cd();
    gPad->SetLogy();
    gPad->SetLeftMargin(0.1);
    gPad->SetRightMargin(0.05);
    gPad->SetBottomMargin(0.2);
    gPad->SetTopMargin(0.1);
    rawProjection[bcNum]->GetYaxis()->SetLabelSize(0.1);
    rawProjection[bcNum]->GetXaxis()->SetLabelSize(0.1);
    rawProjection[bcNum]->GetXaxis()->SetTitle("Raw Time (ns)");
    rawProjection[bcNum]->GetXaxis()->SetTitleSize(0.09);
    rawProjection[bcNum]->SetTitle(Form("Raw Time Projection for LHC10d pass4, BC%d", bcNum));
//    rawProjection[bcNum]->SetTitleSize(0.2);
//    rawProjection[bcNum]->SetTitleOffset(0.2, "h");
    rawProjection[bcNum]->SetLineColor(1);
    rawProjection[bcNum]->SetFillColor(38);
    rawProjection[bcNum]->Draw("SAME");

    //corrected time
    TCanvas* cCorr = new TCanvas("corrtime", "corrtime", 1000,600);
    cCorr->cd();
    gPad->SetRightMargin(0.15);
    gPad->SetLogz();
    corrTime[bcNum]->Draw("SAME COLZ");

    //corrected time projection
    TCanvas* cCorrProj = new TCanvas("corrproj", "corrproj", 1000,300);
    cCorrProj->cd();
    gPad->SetLogy();
    gPad->SetLeftMargin(0.1);
    gPad->SetRightMargin(0.05);
    gPad->SetBottomMargin(0.2);
    gPad->SetTopMargin(0.1);
    corrProjection[bcNum]->GetYaxis()->SetLabelSize(0.1);
    corrProjection[bcNum]->GetXaxis()->SetLabelSize(0.1);
    corrProjection[bcNum]->GetXaxis()->SetTitle("Corrected Time (ns)");
    corrProjection[bcNum]->GetXaxis()->SetTitleSize(0.09);
    corrProjection[bcNum]->SetTitle(Form("Corrected Time Projection for LHC10d pass4, BC%d", bcNum));
//    corrProjection[bcNum]->SetTitleSize(.2);
//    corrProjection[bcNum]->SetTitleOffset(0.2, "h");
    corrProjection[bcNum]->SetLineColor(1);
    corrProjection[bcNum]->SetFillColor(38);
    corrProjection[bcNum]->Draw("SAME");

    //OADB offset vs. calibration
    TCanvas* cCalib = new TCanvas("calib","calib", 1000,600);
    cCalib->cd();
    pass3[bcNum]->SetMarkerColor(1);
    pass3[bcNum]->SetMarkerStyle(20);
    pass3[bcNum]->SetMarkerSize(1);
    pass3[bcNum]->SetTitle(Form("Calibration Offsets for LHC10d, BC%d", bcNum));
    pass3[bcNum]->Draw("P SAME");
    calib[bcNum]->SetMarkerColorAlpha(2, 0.5);
    calib[bcNum]->SetMarkerStyle(3);
    calib[bcNum]->SetMarkerSize(0.6);
    calib[bcNum]->Draw("P SAME");

    //just calibration
    TCanvas* cNewCalib = new TCanvas("newCalib", "newCalib", 1000,600);
    cNewCalib->cd();
    calib[bcNum]->SetMarkerStyle(20);
    calib[bcNum]->SetMarkerSize(1);
    calib[bcNum]->SetMarkerColorAlpha(1,1);
    calib[bcNum]->GetXaxis()->SetRangeUser(0,1150);
    calib[bcNum]->GetYaxis()->SetRangeUser(550,675);
    calib[bcNum]->SetTitle(Form("Step 2: Calibration Offsets for LHC10d, BC%d SM0", bcNum));
    calib[bcNum]->Draw("P SAME");

    //comparison betwen oadb and calibration
    TCanvas* cCompare = new TCanvas("compare", "compare", 1000,600);
    cCompare->cd();
    compare[bcNum]->SetLineColor(2);
    compare[bcNum]->Draw("L SAME");

    //Raw Time for single SM
    TCanvas* cSM = new TCanvas("SM","SM", 1000,600);
    cSM->cd();
    gPad->SetLogz();
    rawTime[bcNum]->GetXaxis()->SetRangeUser(0, 1150);
    rawTime[bcNum]->GetYaxis()->SetRangeUser(550,675);
    rawTime[bcNum]->SetTitle(Form("Raw Time vs. Cell ID for LHC10d pass4, BC%d SM0", ibc));
    rawTime[bcNum]->Draw("COLZ SAME");

    //Raw Time for All BC
    TCanvas* cRawAll = new TCanvas("rawAll","rawAll",1000,600);
    cRawAll->Divide(2,2);
    for(ibc = 0; ibc<4; ibc++){
        cRawAll->cd(ibc+1);
        gPad->SetRightMargin(0.15);
        gPad->SetLogz();
        rawTime[ibc]->GetXaxis()->SetRangeUser(0, 4607);
        rawTime[ibc]->GetYaxis()->SetRangeUser(450,800);
        rawTime[ibc]->Draw("SAME COLZ");
    }
       












}
