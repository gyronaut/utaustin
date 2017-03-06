#include <string>

void energyHistoCompare(string dataset){

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


    //Initialize OADB container
    AliOADBContainer *container = new AliOADBContainer("");
    container->InitFromFile("$ALICE_PHYSICS/OADB/EMCAL/EMCALTimeCalib.root", "AliEMCALTimeCalib");
    if(!container){
        fprintf(stderr, "No OADB container!");
        return;
    }
    
    //Find container in OADB
    TObject *calibration = (TObject*) container->GetObject(0, "TimeCalib13");
    if(!calibration){
        fprintf(stderr, "No calibration container found in OADB!");
        return;
    }

    //Get Array for pass1 (13g) and pass4 (13b and 13c) and set up histograms
    TObject *calibPass1 = (TObject*) calibration->FindObject("pass1");
    if(!calibPass1){
        fprintf(stderr, "No pass1!");
        return;
    }
    TH1F *pass1[4];

    TObject* calibPass4 = (TObject*) calibration->FindObject("pass4");
    if(!calibPass4){
        fprintf(stderr, "No pass4!");
        return;
    }
    TH1F *pass4[4];

    //Set up the file to check against OADB
    TFile *checkFile[6];
    checkFile[0] = new TFile(Form("Calibration_%s_1000_400.root", dataset.c_str()));
    checkFile[1] = new TFile(Form("Calibration_%s_1000_200.root", dataset.c_str()));
    checkFile[2] = new TFile(Form("Calibration_%s_800_400.root", dataset.c_str()));
    checkFile[3] = new TFile(Form("Calibration_%s_800_200.root", dataset.c_str()));
    checkFile[4] = new TFile(Form("Calibration_%s_400_400.root", dataset.c_str()));
    checkFile[5] = new TFile(Form("Calibration_%s_400_200.root", dataset.c_str()));

    for(int i=0; i<6; i++){
        if(!checkFile[i]){
            fprintf(stderr, "Check file %i wasn't loaded!", i);
            return;
        }
    }
//  printf("dataset: %s, strcomp(\"13g\", %s) = %i\n", dataset.c_str(), dataset.c_str(), strcmp("13g", dataset.c_str()));
    TH1F **check[6];
    TH1F **errors[6];
    TH1F **errorsRebin[6];
    TH1F **compare[6];
    TH1F **compareRebin[6];
    TH1F **compare1000[6];
    TH1F **compare1000Rebin[6];
    for(int nfile = 0; nfile<6; nfile++){
        check[nfile] = new TH1F*[4];
        compare[nfile] = new TH1F*[4];
        compareRebin[nfile] = new TH1F*[4];
        errors[nfile] = new TH1F*[4];
        errorsRebin[nfile] = new TH1F*[4];
        compare1000[nfile] = new TH1F*[4];
        compare1000Rebin[nfile] = new TH1F*[4];
    }
    TH1F *compareSingle[4];
    TH1F *errors1000400upper[4];
    TH1F *errors1000400lower[4];

    //Get histograms from files (not LG for dataset 13g)
    for(int ibc =0; ibc<4; ibc++){
       for(int nfile = 0; nfile<6; nfile++){ 
            if(strcmp("13g", dataset.c_str())==0){
                check[nfile][ibc] = (TH1F*)checkFile[nfile]->Get(Form("hAllTimeAvBC%d", ibc));
                errors[nfile][ibc] = (TH1F*)checkFile[nfile]->Get(Form("hAllTimeRMSBC%d", ibc));
            }else{
                check[nfile][ibc] = (TH1F*)checkFile[nfile]->Get(Form("hAllTimeAvLGBC%d", ibc));
                errors[nfile][ibc] = (TH1F*)checkFile[nfile]->Get(Form("hAllTimeRMSLGBC%d", ibc));
            }
            switch(nfile){
                case 0 : check[nfile][ibc]->SetName(Form("1000mev400mevBC%d", ibc));
                    break;
                case 1 : check[nfile][ibc]->SetName(Form("1000mev200mevBC%d", ibc));
                        break;
                case 2 : check[nfile][ibc]->SetName(Form("800mev400mevBC%d", ibc));
                       break;
                case 3 : check[nfile][ibc]->SetName(Form("800mev200mevBC%d", ibc));
                         break;
                case 4 : check[nfile][ibc]->SetName(Form("400mev400mevBC%d", ibc));
                         break;
                case 5 : check[nfile][ibc]->SetName(Form("400mev200mevBC%d", ibc));
                         break;
                default: break;
            }
            compare[nfile][ibc] = check[nfile][ibc]->Clone(Form("compare%s", check[nfile][ibc]->GetName()));
            compare1000[nfile][ibc] = check[0][ibc]->Clone(Form("compare1000_%s", check[nfile][ibc]->GetName()));
      }
        pass1[ibc] = (TH1F*)calibPass1->FindObject(Form("hAllTimeAvBC%d", ibc));
        pass4[ibc] = (TH1F*)calibPass4->FindObject(Form("hAllTimeAvBC%d", ibc));
    }
   
    //Generate the comparison histograms
    for(int nfile=0; nfile<6; nfile++){
        for(int jbc=0; jbc<4; jbc++){
            for(int bin=0; bin < compare[nfile][jbc]->GetNbinsX(); bin++){
                if(strcmp("13g", dataset.c_str())==0){
                    compare[nfile][jbc]->SetBinContent(bin, compare[nfile][jbc]->GetBinContent(bin) - pass1[jbc]->GetBinContent(bin));
                    if(compare[nfile][jbc]->GetBinContent(bin)<-200){
                        compare[nfile][jbc]->SetBinContent(bin, 0.);
                    } 
                }else{
                    compare[nfile][jbc]->SetBinContent(bin, compare[nfile][jbc]->GetBinContent(bin) - pass4[jbc]->GetBinContent(bin));
                    if(compare[nfile][jbc]->GetBinContent(bin)<-200){
                        compare[nfile][jbc]->SetBinContent(bin, 0.);
                    }
                }
                compare1000[nfile][jbc]->SetBinContent(bin, check[nfile][jbc]->GetBinContent(bin) - compare1000[nfile][jbc]->GetBinContent(bin));
            }
            if(nfile==0){
                compareSingle[jbc] = (TH1F*)compare[nfile][jbc]->Clone("singlecomp");
            }
            compareRebin[nfile][jbc] = compare[nfile][jbc]->Rebin(100, Form("compareR%s", check[nfile][jbc]->GetName()));
            compareRebin[nfile][jbc]->Scale(0.01);
            compare1000Rebin[nfile][jbc]= compare1000[nfile][jbc]->Rebin(100, Form("compare1000R%s", check[nfile][jbc]->GetName()));
            compare1000Rebin[nfile][jbc]->Scale(0.01);
        }
    }

    Double_t variance = 0.0;
    //Generate the Errors histograms using standard deviation squared
    for(int nfile = 0; nfile<6; nfile++){
        for(int jbc=0; jbc<4; jbc++){
            for(int bin=0; bin < check[nfile][jbc]->GetNbinsX(); bin++){
                variance = TMath::Power(errors[nfile][jbc]->GetBinContent(bin), 2) - TMath::Power(check[nfile][jbc]->GetBinContent(bin), 2);
                if(variance >= 0){
                    errors[nfile][jbc]->SetBinContent(bin, (TMath::Sqrt(variance)));
                }else{
                    errors[nfile][jbc]->SetBinContent(bin, 0.0);
                }
            }
            //rebin
            errorsRebin[nfile][jbc] = errors[nfile][jbc]->Rebin(100, Form("%sR", errors[nfile][jbc]->GetName()));
            errorsRebin[nfile][jbc]->Scale(0.01);
            if(nfile==0){
                errors1000400upper[jbc] = (TH1F*)errors[nfile][jbc]->Clone("errorUpper");
                errors1000400lower[jbc] = (TH1F*)errors[nfile][jbc]->Clone("errorLower");
                errors1000400lower[jbc]->Scale(-1.0);
            }
        }
    }

   errors[0][0]->Print("base"); 
    /*    
//Setup Output file
    TFile *output = new TFile("comparison.root", "RECREATE");
    output->cd();
    check_1000_200[0]->Write();
    output->Close();
*/

    //draw all calibration on same plot
    TCanvas* c0 = new TCanvas("Calibration Offsets BC0", "calibrationoffsetsbc0", 1000,600);
    c0->cd();

    TLegend *legend = new TLegend(0.1,0.7,0.3,0.9);
    legend->SetHeader("Calibration for Cluster/Cell energy cuts");
   
    // bc number to plot 
    int bcNum = 0;

    if(strcmp("13g",dataset.c_str())==0){
        pass1[bcNum]->SetMarkerStyle(20);
        pass1[bcNum]->SetMarkerColor(1);
        pass1[bcNum]->SetMarkerSize(1);
        pass1[bcNum]->SetTitle(Form("Calibration for Different Cluster/Cell Energy Cuts for LHC%s BC%d", dataset.c_str(), bcNum));
        pass1[bcNum]->GetXaxis()->SetTitle("Cell AbsID");
        pass1[bcNum]->GetXaxis()->SetRangeUser(0,11650);
        pass1[bcNum]->GetYaxis()->SetTitle("Calibrated Offset (ns)");
        pass1[bcNum]->GetYaxis()->SetRangeUser(550,650);
        pass1[bcNum]->SetStats(kFALSE);
        pass1[bcNum]->Draw("P SAME");
        legend->AddEntry(pass1[bcNum], "Current OADB Calibration", "p");
    }else{
        pass4[bcNum]->SetMarkerStyle(20);
        pass4[bcNum]->SetMarkerColor(1);
        pass4[bcNum]->SetMarkerSize(1);
        pass4[bcNum]->SetTitle(Form("Calibration for Different Cluster/Cell Energy Cuts for LHC%s BC%d", dataset.c_str(), bcNum));
        pass4[bcNum]->GetXaxis()->SetTitle("Cell AbsID");
        pass4[bcNum]->GetXaxis()->SetRangeUser(0,11650);
        pass4[bcNum]->GetYaxis()->SetTitle("Calibrated Offset (ns)");
        pass4[bcNum]->GetYaxis()->SetRangeUser(550,650);
        pass4[bcNum]->SetStats(kFALSE);
        pass4[bcNum]->Draw("P SAME");
        legend->AddEntry(pass4[bcNum], "Current OADB Calibration", "p");
    }

    for(int nfile=0; nfile<6; nfile++){
        check[nfile][bcNum]->SetMarkerStyle(3);
        check[nfile][bcNum]->SetMarkerSize(0.6);
        int color = nfile+2;
        check[nfile][bcNum]->SetMarkerColorAlpha(color, 0.3);

        switch(nfile){
            case 0 : legend->AddEntry((TH1F*)check[nfile][bcNum], "1000MeV / 400MeV", "p");
                     break;
            case 1 : legend->AddEntry((TH1F*)check[nfile][bcNum], "1000MeV / 200MeV", "p");
                     break;
            case 2 : legend->AddEntry((TH1F*)check[nfile][bcNum], "800MeV / 400MeV", "p");
                     break;
            case 3 : legend->AddEntry((TH1F*)check[nfile][bcNum], "800MeV / 200MeV", "p");
                     break;
            case 4 : legend->AddEntry((TH1F*)check[nfile][bcNum], "400MeV / 400MeV", "p");
                     break;
            case 5 : legend->AddEntry((TH1F*)check[nfile][bcNum], "400MeV / 200MeV", "p");
                     break;
            default: break;
        }
        check[nfile][bcNum]->SetTitle(Form("Calibration for Different Cluster/Cell Energy Cuts for LHC%s BC%d", dataset.c_str(), bcNum));
        check[nfile][bcNum]->GetXaxis()->SetTitle("Cell AbsID");
        check[nfile][bcNum]->GetXaxis()->SetRangeUser(0,11650);
        check[nfile][bcNum]->GetYaxis()->SetTitle("Calibrated Offset (ns)");
        check[nfile][bcNum]->GetYaxis()->SetRangeUser(550,650);
        check[nfile][bcNum]->SetStats(kFALSE);
        check[nfile][bcNum]->Draw("P SAME");
    }
    legend->SetTextSizePixels(10);
    legend->Draw();

    //draw just OADB 1000/400 MeV cut on same plot

    TCanvas* c400 = new TCanvas("Calibration Offsets 1000MeV/400MeV BC0", "calibrationoffsets400bc0", 1000,600);
    c400->cd();

    TLegend *legend400 = new TLegend(0.1,0.7,0.3,0.9);
    legend400->SetHeader("Calibration for Cluster/Cell energy cuts");



    nfile=0;    
    int color = nfile+2;
    if(strcmp("13g", dataset.c_str())==0){
        pass1[bcNum]->SetTitle(Form("Time Calibration for LHC%s BC%d", dataset.c_str(), bcNum));
        pass1[bcNum]->Draw("P SAME");
        legend400->AddEntry((TH1F*)pass1[bcNum], "OADB Calibration", "p");
    }else{
        pass4[bcNum]->SetTitle(Form("Time Calibration for LHC%s BC%d", dataset.c_str(), bcNum));
        pass4[bcNum]->Draw("P SAME");
        legend400->AddEntry((TH1F*)pass4[bcNum], "OADB Calibration", "p");
    }

    legend400->AddEntry((TH1F*)check[nfile][bcNum], "1000MeV / 400MeV", "p");
    check[nfile][bcNum]->SetTitle(Form("Time Calibration for LHC%s BC%d", dataset.c_str(), bcNum));
    check[nfile][bcNum]->GetYaxis()->SetRangeUser(550,650);
    check[nfile][bcNum]->GetXaxis()->SetRangeUser(0,11650);
    check[nfile][bcNum]->Draw("P SAME");

    legend400->SetTextSizePixels(12);
    legend400->Draw();


    //draw all standard deviation on same plot
    TCanvas* cError = new TCanvas("Calibration Errors BC0", "calibrationerrorsbc0", 1000,600);
    cError->cd();

    TLegend *legendError = new TLegend(0.1,0.7,0.4,0.9);
    legendError->SetHeader("Calibration for Cluster/Cell Energy cuts");

    for(int nfile=0; nfile<6; nfile++){
        errors[nfile][bcNum]->SetMarkerStyle(3);
        errors[nfile][bcNum]->SetMarkerSize(0.6);
        int color = nfile+2;
        errors[nfile][bcNum]->SetMarkerColorAlpha(color, 0.3);
        errors[nfile][bcNum]->SetLineWidth(2);
        errors[nfile][bcNum]->SetLineColor(color);
        switch(nfile){
            case 0 : legendError->AddEntry((TH1F*)errors[nfile][bcNum], "1000MeV / 400MeV", "l");
                     break;
            case 1 : legendError->AddEntry((TH1F*)errors[nfile][bcNum], "1000MeV / 200MeV", "l");
                     break;
            case 2 : legendError->AddEntry((TH1F*)errors[nfile][bcNum], "800MeV / 400MeV", "l");
                     break;
            case 3 : legendError->AddEntry((TH1F*)errors[nfile][bcNum], "800MeV / 200MeV", "l");
                     break;
            case 4 : legendError->AddEntry((TH1F*)errors[nfile][bcNum], "400MeV / 400MeV", "l");
                     break;
            case 5 : legendError->AddEntry((TH1F*)errors[nfile][bcNum], "400MeV / 200MeV", "l");
                     break;
            default: break;
        }
        errors[nfile][bcNum]->SetTitle(Form("Standard Deviation of Calibration for LHC%s BC%d", dataset.c_str(), bcNum));
        errors[nfile][bcNum]->GetXaxis()->SetRangeUser(0, 11650);
        errors[nfile][bcNum]->GetXaxis()->SetTitle("Cell AbsID");
        errors[nfile][bcNum]->GetYaxis()->SetRangeUser(-2,8);
        errors[nfile][bcNum]->GetYaxis()->SetTitle("Standard Deviation of Calibration (ns)");
        errors[nfile][bcNum]->SetStats(kFALSE);
        errors[nfile][bcNum]->Draw("L SAME");
    }
    legendError->SetTextSizePixels(10);
//    legendError->Draw();
 
//Draw comparison with OADB plot (all cuts)
    TCanvas* c1 = new TCanvas("(OADB - Calibration)", "comparisonbc0", 1000,600);
    c1->cd();

    TLegend *compLegend = new TLegend(0.1,0.7,0.4,0.9);
    compLegend->SetHeader("Difference between OADB and Calibration for Cluster/Cell energy cuts");

    for(int nfile=0; nfile<6; nfile++){
        compare[nfile][bcNum]->SetMarkerStyle(3);
        compare[nfile][bcNum]->SetMarkerSize(0.6);
        int color = nfile+2;
        compare[nfile][bcNum]->SetMarkerColorAlpha(color, 0.3);
        compare[nfile][bcNum]->SetLineWidth(2);
        compare[nfile][bcNum]->SetLineColor(color);
        switch(nfile){
            case 0 : compLegend->AddEntry((TH1F*)compare[nfile][bcNum], "(OADB - 1000MeV/400MeV)", "l");
                     break;
            case 1 : compLegend->AddEntry((TH1F*)compare[nfile][bcNum], "(OADB - 1000MeV/200MeV)", "l");
                     break;
            case 2 : compLegend->AddEntry((TH1F*)compare[nfile][bcNum], "(OADB - 800MeV/400MeV)", "l");
                     break;
            case 3 : compLegend->AddEntry((TH1F*)compare[nfile][bcNum], "(OADB - 800MeV/200MeV)", "l");
                     break;
            case 4 : compLegend->AddEntry((TH1F*)compare[nfile][bcNum], "(OADB - 400MeV/400MeV)", "l");
                     break;
            case 5 : compLegend->AddEntry((TH1F*)compare[nfile][bcNum], "(OADB - 400MeV/200MeV)", "l");
                     break;
            default: break;
        }
        compare[nfile][bcNum]->SetTitle(Form("(OADB - Calibration) for LHC%s BC%d", dataset.c_str(), bcNum));
        compare[nfile][bcNum]->GetXaxis()->SetRangeUser(0,11650);
        compare[nfile][bcNum]->GetXaxis()->SetTitle("Cell AbsID");
        compare[nfile][bcNum]->GetYaxis()->SetRangeUser(-3,10);
        compare[nfile][bcNum]->GetYaxis()->SetTitle("Difference in Calibration (ns)"); 
        compare[nfile][bcNum]->SetStats(kFALSE);
        compare[nfile][bcNum]->Draw("L SAME");
    }
    compLegend->SetTextSizePixels(10);
//    compLegend->Draw();

//Draw comparison with OADB (400 cell cut, all cluster cuts)
    TCanvas* c2 = new TCanvas("(OADB - Calibration_400)", "comparison_400_bc0", 1000,600);
    c2->cd();

    TLegend *compLegend2 = new TLegend(0.1,0.7,0.4,0.9);
    compLegend2->SetHeader("Difference between OADB and Calibration for Cluster/Cell energy cuts");

    for(int nfile=0; nfile<6; nfile+=2){
        compare[nfile][bcNum]->SetMarkerStyle(3);
        compare[nfile][bcNum]->SetMarkerSize(0.8);
        int color = nfile+2;
        compare[nfile][bcNum]->SetMarkerColorAlpha(color, 0.3);
        compare[nfile][bcNum]->SetLineColor(color);
        switch(nfile){
            case 0 : compLegend2->AddEntry((TH1F*)compare[nfile][bcNum], "(OADB - 1000MeV/400MeV)", "l");
                     break;
            case 1 : compLegend2->AddEntry((TH1F*)compare[nfile][bcNum], "(OADB - 1000MeV/200MeV)", "l");
                     break;
            case 2 : compLegend2->AddEntry((TH1F*)compare[nfile][bcNum], "(OADB - 800MeV/400MeV)", "l");
                     break;
            case 3 : compLegend2->AddEntry((TH1F*)compare[nfile][bcNum], "(OADB - 800MeV/200MeV)", "l");
                     break;
            case 4 : compLegend2->AddEntry((TH1F*)compare[nfile][bcNum], "(OADB - 400MeV/400MeV)", "l");
                     break;
            case 5 : compLegend2->AddEntry((TH1F*)compare[nfile][bcNum], "(OADB - 400MeV/200MeV)", "l");
                     break;
            default: break;
        }
        compare[nfile][bcNum]->Draw("L SAME");
    }
    compLegend2->SetTextSizePixels(10);
//    compLegend2->Draw();

    //Draw comparison with OADB (200 cell cut, cluster cuts)
    TCanvas* c3 = new TCanvas("(OADB - Calibration_200)", "comparison_200_bc0", 1000,600);
    c3->cd();

    TLegend *compLegend3 = new TLegend(0.1,0.7,0.4,0.9);
    compLegend3->SetHeader("Difference between OADB and Calibration for Cluster/Cell energy cuts");

    for(int nfile=1; nfile<6; nfile+=2){
        compare[nfile][bcNum]->SetMarkerStyle(3);
        compare[nfile][bcNum]->SetMarkerSize(0.8);
        int color = nfile+2;
        compare[nfile][bcNum]->SetMarkerColorAlpha(color, 0.3);
        compare[nfile][bcNum]->SetLineColor(color);
        switch(nfile){
            case 0 : compLegend3->AddEntry((TH1F*)compare[nfile][bcNum], "(OADB - 1000MeV/400MeV)", "l");
                     break;
            case 1 : compLegend3->AddEntry((TH1F*)compare[nfile][bcNum], "(OADB - 1000MeV/200MeV)", "l");
                     break;
            case 2 : compLegend3->AddEntry((TH1F*)compare[nfile][bcNum], "(OADB - 800MeV/400MeV)", "l");
                     break;
            case 3 : compLegend3->AddEntry((TH1F*)compare[nfile][bcNum], "(OADB - 800MeV/200MeV)", "l");
                     break;
            case 4 : compLegend3->AddEntry((TH1F*)compare[nfile][bcNum], "(OADB - 400MeV/400MeV)", "l");
                     break;
            case 5 : compLegend3->AddEntry((TH1F*)compare[nfile][bcNum], "(OADB - 400MeV/200MeV)", "l");
                     break;
            default: break;
        }
        compare[nfile][bcNum]->Draw("L SAME");
    }
    compLegend3->SetTextSizePixels(10);
//    compLegend3->Draw();

//Draw comparison with OADB and 1000/400 cut
    TCanvas* c400comp = new TCanvas("(OADB - Calibration_1000/400)", "comparison_1000400_bc0", 1000,600);
    c400comp->cd();

    nfile=0;
    compareSingle[bcNum]->SetMarkerStyle(3);
    compareSingle[bcNum]->SetMarkerSize(0.8);
    int color = nfile+2;
    compareSingle[bcNum]->SetMarkerColorAlpha(color, 0.3);
    compareSingle[bcNum]->SetLineColor(color);
    compareSingle[bcNum]->SetTitle(Form("(OADB - Calibration) for LHC%s BC%d", dataset.c_str(), bcNum));
    compareSingle[bcNum]->GetXaxis()->SetRangeUser(0,11650);
    compareSingle[bcNum]->GetXaxis()->SetTitle("Cell AbsID");
    compareSingle[bcNum]->GetYaxis()->SetRangeUser(-3,10);
    compareSingle[bcNum]->GetYaxis()->SetTitle("Difference in Calibration (ns)"); 
    compareSingle[bcNum]->SetStats(kFALSE); 
    compareSingle[bcNum]->Draw("L SAME");

//Draw comparison with 1000_400 plot (all cuts)
    TCanvas* c1000 = new TCanvas("(Calibration - 1000MeV/400MeV)", "comparison1000bc0", 1000,600);
    c1000->cd();

    TLegend *comp1000Legend = new TLegend(0.1,0.7,0.4,0.9);
    comp1000Legend->SetHeader("Difference between 1000MeV/400MeV cut and other Cluster/Cell energy cuts");

    for(int nfile=0; nfile<6; nfile++){
        compare1000[nfile][bcNum]->SetMarkerStyle(3);
        compare1000[nfile][bcNum]->SetMarkerSize(0.6);
        int color = nfile+2;
        compare1000[nfile][bcNum]->SetMarkerColorAlpha(color, 0.3);
        compare1000[nfile][bcNum]->SetLineWidth(2);
        compare1000[nfile][bcNum]->SetLineColor(color);
        switch(nfile){
            case 0 : comp1000Legend->AddEntry((TH1F*)compare1000[nfile][bcNum], "(1000MeV/400MeV - 1000MeV/400MeV)", "l");
                     break;
            case 1 : comp1000Legend->AddEntry((TH1F*)compare1000[nfile][bcNum], "(1000MeV/200MeV - 1000MeV/400MeV)", "l");
                     break;
            case 2 : comp1000Legend->AddEntry((TH1F*)compare1000[nfile][bcNum], "(800MeV/400MeV - 1000MeV/400MeV)", "l");
                     break;
            case 3 : comp1000Legend->AddEntry((TH1F*)compare1000[nfile][bcNum], "(800MeV/200MeV - 1000MeV/400MeV)", "l");
                     break;
            case 4 : comp1000Legend->AddEntry((TH1F*)compare1000[nfile][bcNum], "(400MeV/400MeV - 1000MeV/400MeV)", "l");
                     break;
            case 5 : comp1000Legend->AddEntry((TH1F*)compare1000[nfile][bcNum], "(400MeV/200MeV - 1000MeV/400MeV)", "l");
                     break;
            default: break;
        }
        compare1000[nfile][bcNum]->SetStats(kFALSE);
        compare1000[nfile][bcNum]->GetXaxis()->SetRangeUser(0,11650);
        compare1000[nfile][bcNum]->GetXaxis()->SetTitle("Cell AbsID");
        compare1000[nfile][bcNum]->GetYaxis()->SetRangeUser(-5, 15);
        compare1000[nfile][bcNum]->GetYaxis()->SetTitle("Difference in Calibration (ns)");
        compare1000[nfile][bcNum]->SetTitle(Form("(Calibration - 1000MeV/400MeV) for LHC%s BC%d", dataset.c_str(), bcNum));
        compare1000[nfile][bcNum]->SetStats(kFALSE);
        compare1000[nfile][bcNum]->Draw("L SAME");
    }
    comp1000Legend->SetTextSizePixels(10);
//    comp1000Legend->Draw();



//**REBIN draw all standard deviation on same plot
    TCanvas* cErrorR = new TCanvas("Calibration Errors BC0 Rebin", "calibrationerrorsbc0R", 1000,600);
    cErrorR->cd();

    TLegend *legendErrorR = new TLegend(0.1,0.7,0.4,0.9);
    legendError->SetHeader("Calibration for Cluster/Cell Energy cuts");

    for(int nfile=0; nfile<6; nfile++){
        errorsRebin[nfile][bcNum]->SetMarkerStyle(3);
        errorsRebin[nfile][bcNum]->SetMarkerSize(0.6);
        int color = nfile+2;
        errorsRebin[nfile][bcNum]->SetMarkerColorAlpha(color, 0.3);
        errorsRebin[nfile][bcNum]->SetLineWidth(2);
        errorsRebin[nfile][bcNum]->SetLineColor(color);
        switch(nfile){
            case 0 : legendError->AddEntry((TH1F*)errorsRebin[nfile][bcNum], "1000MeV / 400MeV", "l");
                     break;
            case 1 : legendError->AddEntry((TH1F*)errorsRebin[nfile][bcNum], "1000MeV / 200MeV", "l");
                     break;
            case 2 : legendError->AddEntry((TH1F*)errorsRebin[nfile][bcNum], "800MeV / 400MeV", "l");
                     break;
            case 3 : legendError->AddEntry((TH1F*)errorsRebin[nfile][bcNum], "800MeV / 200MeV", "l");
                     break;
            case 4 : legendError->AddEntry((TH1F*)errorsRebin[nfile][bcNum], "400MeV / 400MeV", "l");
                     break;
            case 5 : legendError->AddEntry((TH1F*)errorsRebin[nfile][bcNum], "400MeV / 200MeV", "l");
                     break;
            default: break;
        }
        errorsRebin[nfile][bcNum]->SetTitle(Form("Standard Deviation of Calibration for LHC%s BC%d", dataset.c_str(), bcNum));
        errorsRebin[nfile][bcNum]->GetXaxis()->SetRangeUser(0, 11650);
        errorsRebin[nfile][bcNum]->GetXaxis()->SetTitle("Cell AbsID (Rebinned by 100)");
        errorsRebin[nfile][bcNum]->GetYaxis()->SetRangeUser(-2,8);
        errorsRebin[nfile][bcNum]->GetYaxis()->SetTitle("Standard Deviation of Calibration (ns)");
        errorsRebin[nfile][bcNum]->SetStats(kFALSE);
        errorsRebin[nfile][bcNum]->Draw("L SAME");
    }
    legendError->SetTextSizePixels(10);
//    legendError->Draw();
 

//**REBIN Draw comparison with OADB plot (all cuts)
    TCanvas* c1R = new TCanvas("(OADB - Calibration) Rebin", "comparisonbc0R", 1000,600);
    c1R->cd();

    TLegend *compLegendR = new TLegend(0.1,0.7,0.4,0.9);
    compLegend->SetHeader("Difference between OADB and Calibration for Cluster/Cell energy cuts");

    for(int nfile=0; nfile<6; nfile++){
        compareRebin[nfile][bcNum]->SetMarkerStyle(3);
        compareRebin[nfile][bcNum]->SetMarkerSize(0.6);
        int color = nfile+2;
        compareRebin[nfile][bcNum]->SetMarkerColorAlpha(color, 0.3);
        compareRebin[nfile][bcNum]->SetLineWidth(2);
        compareRebin[nfile][bcNum]->SetLineColor(color);
        switch(nfile){
            case 0 : compLegend->AddEntry((TH1F*)compareRebin[nfile][bcNum], "(OADB - 1000MeV/400MeV)", "l");
                     break;
            case 1 : compLegend->AddEntry((TH1F*)compareRebin[nfile][bcNum], "(OADB - 1000MeV/200MeV)", "l");
                     break;
            case 2 : compLegend->AddEntry((TH1F*)compareRebin[nfile][bcNum], "(OADB - 800MeV/400MeV)", "l");
                     break;
            case 3 : compLegend->AddEntry((TH1F*)compareRebin[nfile][bcNum], "(OADB - 800MeV/200MeV)", "l");
                     break;
            case 4 : compLegend->AddEntry((TH1F*)compareRebin[nfile][bcNum], "(OADB - 400MeV/400MeV)", "l");
                     break;
            case 5 : compLegend->AddEntry((TH1F*)compareRebin[nfile][bcNum], "(OADB - 400MeV/200MeV)", "l");
                     break;
            default: break;
        }
        compareRebin[nfile][bcNum]->SetTitle(Form("(OADB - Calibration) for LHC%s BC%d", dataset.c_str(), bcNum));
        compareRebin[nfile][bcNum]->GetXaxis()->SetRangeUser(0,11650);
        compareRebin[nfile][bcNum]->GetXaxis()->SetTitle("Cell AbsID (Rebinned by 100)");
        compareRebin[nfile][bcNum]->GetYaxis()->SetRangeUser(-3,10);
        compareRebin[nfile][bcNum]->GetYaxis()->SetTitle("Difference in Calibration (ns)"); 
        compareRebin[nfile][bcNum]->SetStats(kFALSE);
        compareRebin[nfile][bcNum]->Draw("L SAME");
    }
    compLegend->SetTextSizePixels(10);
//    compLegend->Draw();


//**REBIN Draw comparison with OADB (400 cell cut, all cluster cuts)
    TCanvas* c2R = new TCanvas("(OADB - Calibration_400) Rebin", "comparison_400_bc0R", 1000,600);
    c2R->cd();

    TLegend *compLegend2R = new TLegend(0.1,0.7,0.4,0.9);
    compLegend2->SetHeader("Difference between OADB and Calibration for Cluster/Cell energy cuts");

    for(int nfile=0; nfile<6; nfile+=2){
        compareRebin[nfile][bcNum]->SetMarkerStyle(3);
        compareRebin[nfile][bcNum]->SetMarkerSize(0.8);
        int color = nfile+2;
        compareRebin[nfile][bcNum]->SetMarkerColorAlpha(color, 0.3);
        compareRebin[nfile][bcNum]->SetLineColor(color);
        switch(nfile){
            case 0 : compLegend2->AddEntry((TH1F*)compareRebin[nfile][bcNum], "(OADB - 1000MeV/400MeV)", "l");
                     break;
            case 1 : compLegend2->AddEntry((TH1F*)compareRebin[nfile][bcNum], "(OADB - 1000MeV/200MeV)", "l");
                     break;
            case 2 : compLegend2->AddEntry((TH1F*)compareRebin[nfile][bcNum], "(OADB - 800MeV/400MeV)", "l");
                     break;
            case 3 : compLegend2->AddEntry((TH1F*)compareRebin[nfile][bcNum], "(OADB - 800MeV/200MeV)", "l");
                     break;
            case 4 : compLegend2->AddEntry((TH1F*)compareRebin[nfile][bcNum], "(OADB - 400MeV/400MeV)", "l");
                     break;
            case 5 : compLegend2->AddEntry((TH1F*)compareRebin[nfile][bcNum], "(OADB - 400MeV/200MeV)", "l");
                     break;
            default: break;
        }
        compareRebin[nfile][bcNum]->Draw("L SAME");
    }
    compLegend2->SetTextSizePixels(10);
//    compLegend2->Draw();

//**REBIN Draw comparison with OADB (200 cell cut, cluster cuts)
    TCanvas* c3R = new TCanvas("(OADB - Calibration_200) Rebin", "comparison_200_bc0R", 1000,600);
    c3R->cd();

    TLegend *compLegend3R = new TLegend(0.1,0.7,0.4,0.9);
    compLegend3->SetHeader("Difference between OADB and Calibration for Cluster/Cell energy cuts");

    for(int nfile=1; nfile<6; nfile+=2){
        compareRebin[nfile][bcNum]->SetMarkerStyle(3);
        compareRebin[nfile][bcNum]->SetMarkerSize(0.8);
        int color = nfile+2;
        compareRebin[nfile][bcNum]->SetMarkerColorAlpha(color, 0.3);
        compareRebin[nfile][bcNum]->SetLineColor(color);
        switch(nfile){
            case 0 : compLegend3->AddEntry((TH1F*)compareRebin[nfile][bcNum], "(OADB - 1000MeV/400MeV)", "l");
                     break;
            case 1 : compLegend3->AddEntry((TH1F*)compareRebin[nfile][bcNum], "(OADB - 1000MeV/200MeV)", "l");
                     break;
            case 2 : compLegend3->AddEntry((TH1F*)compareRebin[nfile][bcNum], "(OADB - 800MeV/400MeV)", "l");
                     break;
            case 3 : compLegend3->AddEntry((TH1F*)compareRebin[nfile][bcNum], "(OADB - 800MeV/200MeV)", "l");
                     break;
            case 4 : compLegend3->AddEntry((TH1F*)compareRebin[nfile][bcNum], "(OADB - 400MeV/400MeV)", "l");
                     break;
            case 5 : compLegend3->AddEntry((TH1F*)compareRebin[nfile][bcNum], "(OADB - 400MeV/200MeV)", "l");
                     break;
            default: break;
        }
        compareRebin[nfile][bcNum]->Draw("L SAME");
    }
    compLegend3->SetTextSizePixels(10);
//    compLegend3->Draw();


//**REBIN Draw comparison with 1000_400 plot (all cuts)
    TCanvas* c1000R = new TCanvas("(Calibration - 1000MeV/400MeV) Rebin", "comparison1000bc0R", 1000,600);
    c1000R->cd();

    TLegend *comp1000LegendR = new TLegend(0.1,0.7,0.4,0.9);
    comp1000LegendR->SetHeader("Difference between 1000MeV/400MeV cut and other Cluster/Cell energy cuts");

    for(int nfile=0; nfile<6; nfile++){
        compare1000Rebin[nfile][bcNum]->SetMarkerStyle(3);
        compare1000Rebin[nfile][bcNum]->SetMarkerSize(0.6);
        int color = nfile+2;
        compare1000Rebin[nfile][bcNum]->SetMarkerColorAlpha(color, 0.3);
        compare1000Rebin[nfile][bcNum]->SetLineWidth(2);
        compare1000Rebin[nfile][bcNum]->SetLineColor(color);
        switch(nfile){
            case 0 : comp1000Legend->AddEntry((TH1F*)compare1000Rebin[nfile][bcNum], "(1000MeV/400MeV - 1000MeV/400MeV)", "l");
                     break;
            case 1 : comp1000Legend->AddEntry((TH1F*)compare1000Rebin[nfile][bcNum], "(1000MeV/200MeV - 1000MeV/400MeV)", "l");
                     break;
            case 2 : comp1000Legend->AddEntry((TH1F*)compare1000Rebin[nfile][bcNum], "(800MeV/400MeV - 1000MeV/400MeV)", "l");
                     break;
            case 3 : comp1000Legend->AddEntry((TH1F*)compare1000Rebin[nfile][bcNum], "(800MeV/200MeV - 1000MeV/400MeV)", "l");
                     break;
            case 4 : comp1000Legend->AddEntry((TH1F*)compare1000Rebin[nfile][bcNum], "(400MeV/400MeV - 1000MeV/400MeV)", "l");
                     break;
            case 5 : comp1000Legend->AddEntry((TH1F*)compare1000Rebin[nfile][bcNum], "(400MeV/200MeV - 1000MeV/400MeV)", "l");
                     break;
            default: break;
        }
        compare1000Rebin[nfile][bcNum]->SetStats(kFALSE);
        compare1000Rebin[nfile][bcNum]->GetXaxis()->SetRangeUser(0,11650);
        compare1000Rebin[nfile][bcNum]->GetXaxis()->SetTitle("Cell AbsID (Rebinned by 100)");
        compare1000Rebin[nfile][bcNum]->GetYaxis()->SetRangeUser(-5, 15);
        compare1000Rebin[nfile][bcNum]->GetYaxis()->SetTitle("Difference in Calibration (ns)");
        compare1000Rebin[nfile][bcNum]->SetTitle(Form("(Calibration - 1000MeV/400MeV) for LHC%s BC%d", dataset.c_str(), bcNum));
        compare1000Rebin[nfile][bcNum]->SetStats(kFALSE);
        compare1000Rebin[nfile][bcNum]->Draw("L SAME");
    }
    comp1000Legend->SetTextSizePixels(10);
//    comp1000Legend->Draw();


}

