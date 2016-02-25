#include <string>

void histoCompare(string checkFileName){
    //Initialize OADB container
    AliOADBContainer *container = new AliOADBContainer("");
    container->InitFromFile("$ALICE_PHYSICS/OADB/EMCAL/EMCALTimeCalib.root", "AliEMCALTimeCalib");
    if(!container){
        fprintf(stderr, "No OADB container!");
        return;
    }
    
    //Find 10d list of objects
    TObject *calibration = (TObject*) container->GetObject(0, "TimeCalib10d");
    if(!calibration){
        fprintf(stderr, "No 10d calibration list!");
        return;
    }

    //Get Array for each pass and set up histograms
    TObject *calibPass1 = (TObject*) calibration->FindObject("pass1");
    if(!calibPass1){
        fprintf(stderr, "No pass1!");
        return;
    }
    TH1F *pass1[4];

    TObject *calibPass2 = (TObject*) calibration->FindObject("pass2");
     if(!calibPass2){
        fprintf(stderr, "No pass2!");
        return;
    }  
 
    TH1F *pass2[4];

    TObject* calibPass3 = (TObject*) calibration->FindObject("pass3");
    if(!calibPass3){
        fprintf(stderr, "No pass3!");
        return;
    }

    TH1F *pass3[4];

    //Set up the file to check
    TFile *checkFile = new TFile(checkFileName.c_str());
    if(!checkFile){
        fprintf(stderr, "Check file wasn't loaded!");
        return;
    }
    TH1F *check[4];
    //Get histograms from files
    for(int i =0; i<4; i++){ 
        check[i] = (TH1F*)checkFile->Get(Form("hAllTimeAvLGBC%d", i));
        pass1[i] = (TH1F*)calibPass1->FindObject(Form("hAllTimeAvBC%d", i));
        pass2[i] = (TH1F*)calibPass2->FindObject(Form("hAllTimeAvBC%d", i));
        pass3[i] = (TH1F*)calibPass3->FindObject(Form("hAllTimeAvBC%d", i));
    }

    //Generate the comparison histograms
    TH1F *compare[12];
    for(int j=0; j<4; j++){
        compare[j] = (TH1F*)check[j%4]->Clone(Form("comparepass4pass1BC%d", j%4));
        for(int bin=0; bin < compare[j]->GetNbinsX(); bin++){
            compare[j]->SetBinContent(bin, compare[j]->GetBinContent(bin) - pass1[j%4]->GetBinContent(bin));
        }
        compare[j]->SetTitle(Form("Difference between pass4 and pass1, BC%d", j%4));
    }
    for(int j=4; j<8; j++){
        compare[j] = (TH1F*)check[j%4]->Clone(Form("comparepass4pass2BC%d", j%4));
        for(int bin=0; bin < compare[j]->GetNbinsX(); bin++){
            compare[j]->SetBinContent(bin, compare[j]->GetBinContent(bin) - pass2[j%4]->GetBinContent(bin));
        }
        compare[j]->SetTitle(Form("Difference between pass4 and pass2, BC%d", j%4));
    }
    for(int j=8; j<12; j++){
        compare[j] = (TH1F*)check[j%4]->Clone(Form("comparepass4pass3BC%d", j%4));
        for(int bin=0; bin < compare[j]->GetNbinsX(); bin++){
            compare[j]->SetBinContent(bin, compare[j]->GetBinContent(bin) - pass3[j%4]->GetBinContent(bin));
        }
        compare[j]->SetTitle(Form("Difference between pass4 and pass3, BC%d", j%4));
    }
    
    //Setup Output file
    TFile *output = new TFile("comparison.root", "RECREATE");
    output->cd();
    for(int j=0; j<12; j++){
        compare[j]->Write();
    }
    output->Close();
}
