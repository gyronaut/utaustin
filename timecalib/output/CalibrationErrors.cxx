
//________________________________________________________________________
/// Calculate calibration constants with errors
/// input - root file with constants in histgrams
/// output - root file with constants and errors
void CalibrationErrors(TString inputFile,TString outputFile)
{
  TFile *file =new TFile(inputFile.Data());
  if(file==0x0) {
    //AliWarning("Input file does not exist!");
    return;
  }

  //high gain
  TH1F *hAllTimeAvBC[4];
  TH1F *hAllTimeRMSBC[4];

  //low gain
  TH1F *hAllTimeAvLGBC[4];
  TH1F *hAllTimeRMSLGBC[4];

  //output histos with Errors
  TH1F *hCalibrationWithErrorsBC[4];
  TH1F *hCalibrationWithErrorsLGBC[4];

  for(Int_t i=0;i<4;i++){
      hAllTimeAvBC[i] = (TH1F *)file->Get(Form("hAllTimeAvBC%d", i));
      hAllTimeRMSBC[i] = (TH1F *)file->Get(Form("hAllTimeRMSBC%d", i));

      hAllTimeAvLGBC[i] = (TH1F *)file->Get(Form("hAllTimeAvLGBC%d", i));
      hAllTimeRMSLGBC[i] = (TH1F *)file->Get(Form("hAllTimeRMSLGBC%d", i));
  
      if(!hAllTimeAvBC[i]){
        printf("No average histo found!!\n");
      }
      if(!hAllTimeRMSBC[i]){
          printf("no RMS histo found!!!\n");
      }
  }
  //AliWarning("Input histograms read.");

  for(Int_t i=0;i<4;i++){
    hCalibrationWithErrorsBC[i] = (TH1F*) hAllTimeAvBC[i]->Clone();
    hCalibrationWithErrorsLGBC[i] = (TH1F*) hAllTimeAvLGBC[i]->Clone();
  }
  
  //AliWarning("New histograms booked.");

  //important remark: we use 'underflow bin' for absid=0 in OADB  . That's why there is j-1 below.
  for(Int_t i=0;i<4;i++){
    for(Int_t j=1;j<=hCalibrationWithErrorsBC[i]->GetNbinsX();j++){
      //high gain
        hCalibrationWithErrorsBC[i]->SetBinError(j, TMath::Sqrt((hAllTimeRMSBC[i]->GetBinContent(j)*(hAllTimeRMSBC[i]->GetBinContent(j)) - (hAllTimeAvBC[i]->GetBinContent(j)*(hAllTimeAvBC[i]->GetBinContent(j))))));    
      //low gain
        hCalibrationWithErrorsLGBC[i]->SetBinError(j, TMath::Sqrt((hAllTimeRMSLGBC[i]->GetBinContent(j)*(hAllTimeRMSLGBC[i]->GetBinContent(j)) - (hAllTimeAvLGBC[i]->GetBinContent(j)*(hAllTimeAvLGBC[i]->GetBinContent(j))))));    

    }
  }

  //AliWarning("Average and rms calculated.");
  TFile *fileNew=new TFile(outputFile.Data(),"recreate");
  for(Int_t i=0;i<4;i++){
    hCalibrationWithErrorsBC[i]->Write();
    hCalibrationWithErrorsLGBC[i]->Write();   
  }

  //AliWarning(Form("Histograms saved in %s file.",outputFile.Data()));

  fileNew->Close();
  delete fileNew;

  for(Int_t i=0;i<4;i++){
    delete hAllTimeAvBC[i];
    delete hAllTimeRMSBC[i];
    delete hAllTimeAvLGBC[i];
    delete hAllTimeRMSLGBC[i];

    delete hCalibrationWithErrorsBC[i];
    delete hCalibrationWithErrorsLGBC[i];
  }
  
  file->Close();
  delete file;

  //AliWarning("Pointers deleted. Memory cleaned.");
}


