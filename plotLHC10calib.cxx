#include <string>

void plotLHC10calib(){
    
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

    /*  OUTLINE:
     *
     *  pointers for all histograms needed:
     *d      Corrected time for each period, BC (2D)
     *d      Projection of corrected for each SM, each period, BC (1D)
     *d      calibration offset for each period, BC
     *
     *  For each period:
     *d      Open file
     *d      Get Calibration offsets in histogram
     *d      Get Corrected Time in histogram
     *      For Each SM:
     *d          project corrected time onto 1D histogram (log plot)
     *
     *
     *  Plot Calibration constants on same plot
     *
     *  Plot each SM projection on same plot (BC merged as well?)
     *
     */

    //Cell ID limits for first 4 SM's
    Int_t lowerlimit[] = { 0, 1152, 2304, 3456};
    Int_t upperlimit[] = { 1151, 2303, 3455, 4607};

    TH2F* rawTime[5][4]; //indices: period, bc
    TH1F* calib[5][4]; //indices: Period, BC
    TH2F* corrTime[5]; //indices: Period
    TH1F* corrProjection[5][4]; //indices: Period, SM


    TFile* file[5];
    TString period;
    TString periodName[5];

    //Set-up all histograms
    for(int iperiod=0; iperiod<5; iperiod++){
        switch(iperiod){
            case 0: {
                        period = "calibration_LHC10b.root";
                        periodName[iperiod]="LHC10b";
                        break;
                    }
            case 1: {
                        period = "calibration_LHC10c.root";
                        periodName[iperiod]="LHC10c";
                        break;
                    }
            case 2: {
                        period = "calibration_LHC10d_pass4_corr.root";
                        periodName[iperiod]="LHC10d";
                        break;
                    }
            case 3: {
                        period = "calibration_LHC10e_pass4_corr.root";
                        periodName[iperiod]="LHC10e";
                        break;
                    }
            case 4: {
                        period = "calibration_LHC10f_pass4_corr.root";
                        periodName[iperiod]="LHC10f";
                        break;
                    }
            default: period = "default";break;
        }
        //printf("name: %s\n\n", period.Data());
        file[iperiod] = new TFile(period.Data());
        //if(file[iperiod])printf("loaded file!\n");
        for(int ibc = 0; ibc < 4; ibc++){
            rawTime[iperiod][ibc] = (TH2F*)file[iperiod]->Get(Form("RawTimeVsIdLGBC%d", ibc)); 
            calib[iperiod][ibc] = (TH1F*)file[iperiod]->Get(Form("hAllTimeAvLGBC%d", ibc));          
            if(calib[iperiod][ibc]){
                //printf("calib hist loaded\n");
            }else{
                printf("not loaded!\n");
            }
            calib[iperiod][ibc]->GetXaxis()->SetTitle("Cell AbsID");
            calib[iperiod][ibc]->GetXaxis()->SetRangeUser(0, 4607);
            calib[iperiod][ibc]->GetYaxis()->SetTitle("Calibration Offset (ns)");
            calib[iperiod][ibc]->GetYaxis()->SetRangeUser(600, 800);
            //printf("got past calib!\n");
            if(ibc == 0){
                corrTime[iperiod] = (TH2F*)file[iperiod]->Get(Form("CorrectedTimeBC%d", ibc));
                corrTime[iperiod]->GetXaxis()->SetTitle("Cell AbsID");
                corrTime[iperiod]->GetXaxis()->SetRangeUser(0, 4607);
                corrTime[iperiod]->GetYaxis()->SetTitle("Corrected Time (ns)");
                corrTime[iperiod]->GetYaxis()->SetRangeUser(-150, 150);
                corrTime[iperiod]->GetZaxis()->SetRangeUser(1, 150);
                corrTime[iperiod]->SetTitle(Form("Corrected Time vs Cell ID, BC%d", ibc));
            }else{
                corrTime[iperiod]->Add((TH2F*)file[iperiod]->Get(Form("CorrectedTimeBC%d", ibc)));
            }
            for(int iSM = 0; iSM<4; iSM++){
                corrProjection[iperiod][iSM] = (TH1F*)corrTime[iperiod]->ProjectionY(Form("corrprojPeriod%dSM%d",iperiod,iSM), lowerlimit[iSM], upperlimit[iSM]);
            }
        }
    }

    int bcNum = 0;

    //raw time
    TCanvas* cRaw[5];
    for(int iperiod = 0; iperiod<5; iperiod++){ 
        cRaw[iperiod] = new TCanvas(Form("rawtime%d", iperiod), Form("rawtime%d", iperiod), 1000,600);
        cRaw[iperiod]->Divide(2,2);
        for(int ibc = 0; ibc < 4; ibc++){
            cRaw[iperiod]->cd(ibc+1);
            gPad->SetRightMargin(0.15);
            gPad->SetLogz();
            rawTime[iperiod][ibc]->SetTitle(Form("Raw Time vs Cell ID for %s, BC%d", periodName[iperiod].Data(), ibc));
            rawTime[iperiod][ibc]->GetXaxis()->SetRangeUser(0, 4607);
            rawTime[iperiod][ibc]->GetYaxis()->SetRangeUser(400, 800);
            rawTime[iperiod][ibc]->Draw("SAME COLZ");
        }
    }

    TLegend *corrLegend = new TLegend(0.75,0.6, 0.9, 0.9);
    corrLegend->SetHeader("Period:");
    corrLegend->SetTextAlign(22);
    
    //corrected time projection (all periods together)
    TCanvas* cCorrProj = new TCanvas("corrproj", "corrproj", 800,600);
    cCorrProj->Divide(2,2);
    for(int i =0; i<4; i++){
        cCorrProj->cd(i+1);
        gPad->SetLogy();
        gPad->SetLeftMargin(0.1);
        gPad->SetRightMargin(0.05);
        gPad->SetBottomMargin(0.2);
        gPad->SetTopMargin(0.1);
        for(int jperiod=0; jperiod<5; jperiod++){
            corrProjection[jperiod][i]->GetYaxis()->SetLabelSize(0.1);
            //corrProjection[jperiod][i]->GetYaxis()->SetRangeUser(0,100000);
            corrProjection[jperiod][i]->GetXaxis()->SetLabelSize(0.1);
            corrProjection[jperiod][i]->GetXaxis()->SetTitle("Corrected Time (ns)");
            corrProjection[jperiod][i]->GetXaxis()->SetTitleSize(0.09);
            corrProjection[jperiod][i]->SetTitle(Form("Corrected Time Projection for SM%d", i));
            //    corrProjection[bcNum]->SetTitleSize(.2);
            //    corrProjection[bcNum]->SetTitleOffset(0.2, "h");
            if (jperiod==0){
                corrProjection[jperiod][i]->SetLineColor(1);
            }else{
                corrProjection[jperiod][i]->SetLineColor((jperiod)*2);
            }
            corrProjection[jperiod][i]->Draw("SAME");
            if(i==2){ //only draw legend for 1 of the BC
                corrLegend->AddEntry(corrProjection[jperiod][i], periodName->Data(), "l");
            }
        }
        if(i==2){
            corrLegend->Draw("SAME HIST");
        }
    }


    //draw corrected time projection for each period separtely, fit with gaussian curves
    TF1* fits[5][4];
    TCanvas *cCorrProjPeriod[5];
    for(int iperiod = 0; iperiod<5; iperiod++){
        cCorrProjPeriod[iperiod] = new TCanvas(Form("corrproj%d", iperiod), Form("corrproj%d", iperiod), 800, 600);
        cCorrProjPeriod[iperiod]->Divide(2,2);
        for(int i =0; i<4; i++){
            cCorrProjPeriod[iperiod]->cd(i+1);
            gPad->SetLeftMargin(0.1);
            gPad->SetRightMargin(0.05);
            gPad->SetBottomMargin(0.2);
            gPad->SetTopMargin(0.1);
            fits[iperiod][i] = new TF1(Form("fCorr%sSM%d", periodName[iperiod].Data(), i), "gaus(0)", -75, 75);
            fits[iperiod][i]->SetParameter(1, 0.0);
            fits[iperiod][i]->SetParLimits(1,-25,25);
            fits[iperiod][i]->SetParName(1, "mean");
            fits[iperiod][i]->SetParameter(2, 5.0);
            fits[iperiod][i]->SetParLimits(2, 2, 20);
            fits[iperiod][i]->SetParName(2, "std dev");
            fits[iperiod][i]->SetParameter(0, 50000);
            fits[iperiod][i]->SetParName(0, "max");

            //fits[iperiod][i]->SetParLimits(0,10000,200000);
            corrProjection[iperiod][i]->GetYaxis()->SetLabelSize(0.1);
            corrProjection[iperiod][i]->GetXaxis()->SetLabelSize(0.1);
            corrProjection[iperiod][i]->GetXaxis()->SetTitle("Corrected Time (ns)");
            corrProjection[iperiod][i]->GetXaxis()->SetTitleSize(0.09);
            corrProjection[iperiod][i]->SetTitle(Form("Corrected Time Projection for %s SM%d",periodName[iperiod].Data(), i));
            if (iperiod==0){
                corrProjection[iperiod][i]->SetLineColor(1);
            }else{
                corrProjection[iperiod][i]->SetLineColor((iperiod)*2);
            }
            printf("Period: %s, SM: %d:\n", periodName[iperiod].Data(), i);
            corrProjection[iperiod][i]->Fit(Form("fCorr%sSM%d", periodName[iperiod].Data(), i), "R");
            printf("\n");
            gPad->SetLogy();
            corrProjection[iperiod][i]->Draw("SAME HIST FUNC");
        }
    }

    TLegend *calibLegend = new TLegend(0.7,0.7, 0.9, 0.9);
    calibLegend->SetHeader("Period:");
    calibLegend->SetTextAlign(22);
    
    //just plot calibration (all periods together, separated by BC)
    TCanvas* cNewCalib = new TCanvas("newCalib", "newCalib", 1000,600);
    cNewCalib->Divide(2,2);
    int color=0;
    for(int j = 0; j < 4; j++){
        cNewCalib->cd(j+1);
        for(int k =0; k<5; k++){
            calib[k][j]->SetMarkerStyle(20);
            calib[k][j]->SetMarkerSize(1.);
            if(k==0){
                color=1;
            }else{
                color=2*k;
            }
            calib[k][j]->SetMarkerColorAlpha(color,0.7);
            calib[k][j]->GetXaxis()->SetRangeUser(0,4607);
            calib[k][j]->GetYaxis()->SetRangeUser(500,800);
            calib[k][j]->SetTitle(Form("Calibration Constansts for BC%d", j));
            calib[k][j]->Draw("P SAME");
            if(j==2){ //only need legend entries for one BC
                calibLegend->AddEntry(calib[k][j], periodName[k].Data(), "p");
            }
        }
        if(j==2){
            calibLegend->Draw("SAME");
        }
    }
}
       

