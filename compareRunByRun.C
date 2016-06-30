#include <string>

void compareRunByRun(){
    TString filebase = "/Users/jtblair/Downloads/AnalysisResults_LHC10f_";
    TString runlist[] = {"133006", "133007", "133010", "133282", "133327", "133329", "133330", "133414", "133419", "133563", "133762", "133800", "133920", "133924", "133969", "133982", "133985", "134198", "134204", "134297", "134299", "134301", "134304"}; //Runlist for 10f_pass4
    int nruns = sizeof(runlist)/sizeof(*runlist);
    TString filename = "";

    //Set up 4 TCanvas (one for each BC)
    TCanvas *c0 = new TCanvas("c0","c0", 1000, 600);
    c0->Divide(6, 4);
    c0->SetLogz();

    TCanvas *c1 = new TCanvas("c1","c1", 1000, 600);
    c1->Divide(6, 4);
    c0->SetLogz();

    TCanvas *c2 = new TCanvas("c2","c2", 1000, 600);
    c2->Divide(6, 4);
    c2->SetLogz();

    TCanvas *c3 = new TCanvas("c3","c3", 1000, 600);
    c3->Divide(6, 4);
    c3->SetLogz();

    TH2F *rawtime[4];
    //Loop over all the files to get the histograms
    for(int i=0; i<nruns; i++){
        //Set up the file to check
        filename = filebase + runlist[i];
        filename += ".root";

        TFile *checkFile = new TFile(filename);
        if(!checkFile){
            fprintf(stderr, "Check file wasn't loaded!");
            return;
        }
        TList *list = (TList*)checkFile->Get("chistolist");
        if(!list){
            printf("Didn't find list! /n");
        }

        //Get histograms from files
        for(int j =0; j<4; j++){ 
            rawtime[j] = (TH2F*)list->FindObject(Form("RawTimeVsIdLGBC%d", j));
            if(!rawtime[j]){
                printf("No histogram for RawTimeVsIdLGBC%d in run %d \n", j, i);
            }
            rawtime[j]->GetXaxis()->SetRangeUser(0, 4607);
            rawtime[j]->SetOption("COL");
            rawtime[j]->SetStats(kFALSE);
            rawtime[j]->GetZaxis()->SetRangeUser(1,10);
            rawtime[j]->SetTitle("");
        }
        c0->cd(i+1);
        gPad->SetLogz();
        rawtime[0]->Draw();
        c1->cd(i+1);
        rawtime[1]->Draw();
        c2->cd(i+1);
        rawtime[2]->Draw();
        c3->cd(i+1);
        rawtime[3]->Draw();

        printf("Made it through run number %d!\n", i);
        checkFile->Close();
    }
}
