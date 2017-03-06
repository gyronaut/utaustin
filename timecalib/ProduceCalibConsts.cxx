
//________________________________________________________________________
/// Calculate calibration constants
/// input - root file with histograms 
/// output - root file with constants in historams
/// isFinal - flag: kFALSE-first iteration, kTRUE-final iteration
void ProduceCalibConsts(TString inputFile,TString outputFile,Bool_t isFinal, TString containerName = "chistolist")
{
  TFile *file =new TFile(inputFile.Data());
  if(file==0x0) {
    //AliWarning("Input file does not exist!");
    return;
  }

  TList *list=(TList*)file->Get(containerName.Data());
  if(list==0x0) 
  {
    //AliWarning("List chistolist does not exist in file!");
    return;
  }

  //high gain
  TH1F *h1[4];
  TH1F *h2[4];
  TH1F *h3[4];
  TH1F *hAllTimeAvBC[4];
  TH1F *hAllTimeRMSBC[4];

  //low gain
  TH1F *h4[4];
  TH1F *h5[4];
  TH1F *h6[4];
  TH1F *hAllTimeAvLGBC[4];
  TH1F *hAllTimeRMSLGBC[4];

  //raw and corrected histos
  TH2F *raw[4];
  TH2F *corrected[4];

  if(isFinal==kFALSE){//first itereation
    for(Int_t i=0;i<4;i++){
      h1[i]=(TH1F *)list->FindObject(Form("RawTimeSumBC%d",i));
      h2[i]=(TH1F *)list->FindObject(Form("RawTimeEntriesBC%d",i));
      h3[i]=(TH1F *)list->FindObject(Form("RawTimeSumSqBC%d",i));
      
      h4[i]=(TH1F *)list->FindObject(Form("RawTimeSumLGBC%d",i));
      h5[i]=(TH1F *)list->FindObject(Form("RawTimeEntriesLGBC%d",i));
      h6[i]=(TH1F *)list->FindObject(Form("RawTimeSumSqLGBC%d",i));
    
//      raw[i]=(TH2F *)list->FindObject(Form("RawTimeVsIdLGBC%d", i));
    }
  } else {//final iteration
    for(Int_t i=0;i<4;i++){
      h1[i]=(TH1F *)list->FindObject(Form("hTimeSum%d",i));
      h2[i]=(TH1F *)list->FindObject(Form("hTimeEnt%d",i));
      h3[i]=(TH1F *)list->FindObject(Form("hTimeSumSq%d",i));
      
      h4[i]=(TH1F *)list->FindObject(Form("hTimeLGSum%d",i));
      h5[i]=(TH1F *)list->FindObject(Form("hTimeLGEnt%d",i));
      h6[i]=(TH1F *)list->FindObject(Form("hTimeLGSumSq%d",i));

      raw[i]=(TH2F *)list->FindObject(Form("RawTimeVsIdLGBC%d", i)); 
    }
  }
  //AliWarning("Input histograms read.");

  for(Int_t i=0;i<4;i++){
    hAllTimeAvBC[i]=new TH1F(Form("hAllTimeAvBC%d",i),Form("hAllTimeAvBC%d",i),h1[i]->GetNbinsX(),h1[i]->GetXaxis()->GetXmin(),h1[i]->GetXaxis()->GetXmax());
    hAllTimeRMSBC[i]=new TH1F(Form("hAllTimeRMSBC%d",i),Form("hAllTimeRMSBC%d",i),h3[i]->GetNbinsX(),h3[i]->GetXaxis()->GetXmin(),h3[i]->GetXaxis()->GetXmax());

    hAllTimeAvLGBC[i]=new TH1F(Form("hAllTimeAvLGBC%d",i),Form("hAllTimeAvLGBC%d",i),h4[i]->GetNbinsX(),h4[i]->GetXaxis()->GetXmin(),h4[i]->GetXaxis()->GetXmax());
    hAllTimeRMSLGBC[i]=new TH1F(Form("hAllTimeRMSLGBC%d",i),Form("hAllTimeRMSLGBC%d",i),h6[i]->GetNbinsX(),h6[i]->GetXaxis()->GetXmin(),h6[i]->GetXaxis()->GetXmax());
  
    if(isFinal==kTRUE){  
        corrected[i]=new TH2F(Form("CorrectedTimeBC%d",i),Form("Corrected Time for BC%d",i),raw[i]->GetNbinsX(),raw[i]->GetXaxis()->GetXmin(),raw[i]->GetXaxis()->GetXmax(),100,-100,100);
    }
  }
  
  //AliWarning("New histograms booked.");

  //important remark: we use 'underflow bin' for absid=0 in OADB  . That's why there is j-1 below.
  for(Int_t i=0;i<4;i++){
    for(Int_t j=1;j<=h1[i]->GetNbinsX();j++){
      //high gain
      if(h2[i]->GetBinContent(j)!=0){
	hAllTimeAvBC[i]->SetBinContent(j-1,h1[i]->GetBinContent(j)/h2[i]->GetBinContent(j));
	hAllTimeRMSBC[i]->SetBinContent(j-1,TMath::Sqrt(h3[i]->GetBinContent(j)/h2[i]->GetBinContent(j)) );
      } else {
	hAllTimeAvBC[i]->SetBinContent(j-1,0.);
	hAllTimeRMSBC[i]->SetBinContent(j-1,0.);
      }
      //low gain
      if(h5[i]->GetBinContent(j)!=0){
	hAllTimeAvLGBC[i]->SetBinContent(j-1,h4[i]->GetBinContent(j)/h5[i]->GetBinContent(j));
	hAllTimeRMSLGBC[i]->SetBinContent(j-1,TMath::Sqrt(h6[i]->GetBinContent(j)/h5[i]->GetBinContent(j)) );
      } else {
	hAllTimeAvLGBC[i]->SetBinContent(j-1,0.);
	hAllTimeRMSLGBC[i]->SetBinContent(j-1,0.);
      }
      if(isFinal==kTRUE){
          double min_raw = (double)(raw[i]->GetYaxis()->GetXmin());
          double min_corr = (double)(corrected[i]->GetYaxis()->GetXmin());
          double scale_raw = (double)((raw[i]->GetYaxis()->GetXmax() - raw[i]->GetYaxis()->GetXmin())/((double)raw[i]->GetNbinsY()));
          double scale_corr = (double)(corrected[i]->GetYaxis()->GetXmax() - corrected[i]->GetYaxis()->GetXmin())/((double)(corrected[i]->GetNbinsY()));
          for(Int_t iy=0; iy<100; iy++){
              calc_ybin = (int)((scale_corr*iy-min_raw+min_corr+hAllTimeAvLGBC[i]->GetBinContent(j-1))/scale_raw);
              if(calc_ybin > 0 && calc_ybin < 200){
                  //printf("calc_ybin=%f\n", calc_ybin);
                  corrected[i]->SetBinContent(j-1,iy,raw[i]->GetBinContent(j-1,calc_ybin));
              }else{
                  corrected[i]->SetBinContent(j-1,iy,0.);
              }
          }
      }

    }
  }

  //AliWarning("Average and rms calculated.");
  TFile *fileNew=new TFile(outputFile.Data(),"recreate");
  for(Int_t i=0;i<4;i++){
    hAllTimeAvBC[i]->Write();
    hAllTimeRMSBC[i]->Write();
    hAllTimeAvLGBC[i]->Write();
    hAllTimeRMSLGBC[i]->Write();
    
    if(isFinal==kTRUE){
        raw[i]->Write();
        corrected[i]->Write();
    }
  }

  //AliWarning(Form("Histograms saved in %s file.",outputFile.Data()));

  fileNew->Close();
  delete fileNew;

  for(Int_t i=0;i<4;i++){
    delete hAllTimeAvBC[i];
    delete hAllTimeRMSBC[i];
    delete hAllTimeAvLGBC[i];
    delete hAllTimeRMSLGBC[i];

    delete h1[i];
    delete h2[i];
    delete h3[i];
    delete h4[i];
    delete h5[i];
    delete h6[i];

    if(isFinal==kTRUE){
        delete raw[i];
        delete corrected[i];
    }
  }
  
  file->Close();
  delete file;

  //AliWarning("Pointers deleted. Memory cleaned.");
}


