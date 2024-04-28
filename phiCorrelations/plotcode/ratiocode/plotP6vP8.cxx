Double_t calcJetFromBG(TH1D* h, Int_t firstbin, Int_t lastbin, Float_t BG, Double_t& error){
    Double_t sum = 0.0;
    error = 0.0;
    for(int ibin = firstbin; ibin <= lastbin; ibin++){
        Float_t pt = h->GetBinContent(ibin);
        if(pt >= BG){
            sum += pt - BG;
            error += TMath::Power(h->GetBinError(ibin), 2);
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

Double_t ZYAM(TH1D* h, Float_t &bgerr){
    Double_t min = 0;
    for(int ibin = 1; ibin <= h->GetNbinsX(); ibin++){
        Double_t bc = h->GetBinContent(ibin);
        if(ibin == 1){
            min = bc;
            bgerr = h->GetBinError(ibin);
        }else if(bc < min){
            min = bc;
            bgerr = h->GetBinError(ibin);
        }
    }
    return min;
}


void plotP6vP8(){
    //primary hadron distribution p6
    //TFile* rawfile = new TFile("~/phiStudies/LHC18j2_FAST/LHC18j2_hphi_alltrig_0_100.root");
    TFile* rawfile = new TFile("~/phiStudies/LHC16k5_MC/LHC16k5a_hphi_0_100.root");
    TList* list = (TList*)rawfile->Get("phiCorr_mult_0_100_");
    THnSparseF* hdist = (THnSparseF*)list->FindObject("fTruePrimHDist");
    hdist->GetAxis(3)->SetRange(hdist->GetAxis(3)->FindBin(-0.49999), hdist->GetAxis(3)->FindBin(0.4999));
    hdist->GetAxis(4)->SetRange(1,1);
    TH1D* pipt = (TH1D*)hdist->Projection(0)->Clone("pipt");
    hdist->GetAxis(4)->SetRange(2,2);
    TH1D* kpt = (TH1D*)hdist->Projection(0)->Clone("kpt");
    hdist->GetAxis(4)->SetRange(3,3);
    TH1D* ppt = (TH1D*)hdist->Projection(0)->Clone("ppt");
    hdist->GetAxis(3)->SetRange(0, -1); //do pseudorapiity cuts instead
    hdist->GetAxis(2)->SetRange(hdist->GetAxis(2)->FindBin(-0.7999), hdist->GetAxis(2)->FindBin(0.7999));
    hdist->GetAxis(4)->SetRange(1,1);
    TH1D* pipt_etacut = (TH1D*)(hdist->Projection(0))->Clone("pipt_etacut");
    hdist->GetAxis(4)->SetRange(2,2);
    TH1D* kpt_etacut = (TH1D*)(hdist->Projection(0))->Clone("kpt_etacut");
    hdist->GetAxis(4)->SetRange(3,3);
    TH1D* ppt_etacut = (TH1D*)(hdist->Projection(0))->Clone("ppt_etacut");

    //primary hadron in triggered events
    THnSparseF* trighdist = (THnSparseF*)list->FindObject("fTriggeredTruePrimHDist");
    trighdist->GetAxis(3)->SetRange(trighdist->GetAxis(3)->FindBin(-0.49999), trighdist->GetAxis(3)->FindBin(0.4999));
    trighdist->GetAxis(4)->SetRange(1,1);
    TH1D* trigpipt = (TH1D*)trighdist->Projection(0)->Clone("trigpipt");
    trighdist->GetAxis(4)->SetRange(2,2);
    TH1D* trigkpt = (TH1D*)trighdist->Projection(0)->Clone("trigkpt");
    trighdist->GetAxis(4)->SetRange(3,3);
    TH1D* trigppt = (TH1D*)trighdist->Projection(0)->Clone("trigppt");
    trighdist->GetAxis(3)->SetRange(0, -1); //do pseudorapiity cuts instead
    trighdist->GetAxis(2)->SetRange(trighdist->GetAxis(2)->FindBin(-0.7999), trighdist->GetAxis(2)->FindBin(0.7999));
    trighdist->GetAxis(4)->SetRange(1,1);
    TH1D* trigpipt_etacut = (TH1D*)(trighdist->Projection(0))->Clone("trigpipt_etacut");
    trighdist->GetAxis(4)->SetRange(2,2);
    TH1D* trigkpt_etacut = (TH1D*)(trighdist->Projection(0))->Clone("trigkpt_etacut");
    trighdist->GetAxis(4)->SetRange(3,3);
    TH1D* trigppt_etacut = (TH1D*)(trighdist->Projection(0))->Clone("trigppt_etacut");

    trighdist->GetAxis(4)->SetRange(1,4);
    trighdist->GetAxis(2)->SetRange(trighdist->GetAxis(2)->FindBin(-0.79999), trighdist->GetAxis(2)->FindBin(0.7999));
    trighdist->GetAxis(3)->SetRange(0,-1);
    TH1D* trighpt = (TH1D*)trighdist->Projection(0)->Clone("trighpt");

    hdist->GetAxis(4)->SetRange(1,4);
    hdist->GetAxis(2)->SetRange(hdist->GetAxis(2)->FindBin(-0.79999), hdist->GetAxis(2)->FindBin(0.7999));
    hdist->GetAxis(3)->SetRange(0,-1);
    TH1D* hpt = (TH1D*)hdist->Projection(0)->Clone("hpt");

    
    THnSparseF* phidist = (THnSparseF*)list->FindObject("fTruePhiDist");
    phidist->GetAxis(4)->SetRange(phidist->GetAxis(4)->FindBin(-0.49999), phidist->GetAxis(4)->FindBin(0.49999));
    TH1D* phipt = (TH1D*)phidist->Projection(0)->Clone("phipt");
    phipt->Scale(1.0/0.49); 
    phidist->GetAxis(4)->SetRange(0, -1);
    phidist->GetAxis(3)->SetRange(phidist->GetAxis(3)->FindBin(-0.79999), phidist->GetAxis(3)->FindBin(0.79999));
    TH1D* phipt_etacut = (TH1D*)(phidist->Projection(0))->Clone("phipt_etacut");
    phipt_etacut->Scale(1.0/0.49);
    THnSparseF* trigphidist = (THnSparseF*)list->FindObject("fTriggeredTruePhiDist");
    trigphidist->GetAxis(4)->SetRange(trigphidist->GetAxis(4)->FindBin(-0.49999), trigphidist->GetAxis(4)->FindBin(0.49999));
    TH1D* trigphipt = (TH1D*)trigphidist->Projection(0)->Clone("trigphipt");
    trigphipt->Scale(1.0/0.49);
    trigphidist->GetAxis(4)->SetRange(0, -1);
    trigphidist->GetAxis(3)->SetRange(trigphidist->GetAxis(3)->FindBin(-0.79999), trigphidist->GetAxis(3)->FindBin(0.79999));
    TH1D* trigphipt_etacut = (TH1D*)(trigphidist->Projection(0))->Clone("trigphipt_etacut");
    trigphipt_etacut->Scale(1.0/0.49);
    

    TH1D* phipiratio = (TH1D*)phipt->Clone("phipiratio");
    phipiratio->Divide(pipt);
    phipiratio->GetYaxis()->SetRangeUser(0, 0.975);
    phipiratio->SetMarkerColor(kSpring+2);
    phipiratio->SetMarkerStyle(23);
    phipiratio->SetMarkerSize(2);

    TH1D* phikratio = (TH1D*)phipt->Clone("phikratio");
    phikratio->Divide(kpt);
    phikratio->GetYaxis()->SetRangeUser(0, 0.975);
    phikratio->SetMarkerColor(kAzure+2);
    phikratio->SetMarkerStyle(22);
    phikratio->SetMarkerSize(2);

    TH1D* pphiratio = (TH1D*)ppt->Clone("pphiratio");
    pphiratio->Divide(phipt);
    pphiratio->Scale(0.5);
    pphiratio->GetYaxis()->SetRangeUser(0, 5.5);
    pphiratio->SetMarkerColor(kBlack);
    pphiratio->SetMarkerStyle(21);
    pphiratio->SetMarkerSize(2);

    TH1D* phipiratio_etacut = (TH1D*)phipt_etacut->Clone("phipiratio_etacut");
    phipiratio_etacut->Divide(pipt_etacut);
    phipiratio_etacut->GetYaxis()->SetRangeUser(0, 0.975);
    phipiratio_etacut->SetMarkerColor(kSpring-7);
    phipiratio_etacut->SetMarkerStyle(23);
    phipiratio_etacut->SetMarkerSize(2);

    TH1D* phikratio_etacut = (TH1D*)phipt_etacut->Clone("phikratio_etacut");
    phikratio_etacut->Divide(kpt_etacut);
    phikratio_etacut->GetYaxis()->SetRangeUser(0, 0.975);
    phikratio_etacut->SetMarkerColor(kAzure+3);
    phikratio_etacut->SetMarkerStyle(22);
    phikratio_etacut->SetMarkerSize(2);

    TCanvas* cphipt = new TCanvas("cphipt", "cphipt", 50, 50, 600, 600);
    cphipt->cd();
    phipt->Draw("HIST");
    phipt_etacut->SetLineColor(kRed);
    phipt_etacut->Draw("SAME HIST");

    TCanvas* cpipt = new TCanvas("cpipt", "cpipt", 50, 50, 600, 600);
    cpipt->cd();
    pipt->Draw("HIST");
    pipt_etacut->SetLineColor(kRed);
    pipt_etacut->Draw("SAME HIST");


    TCanvas* cptratio = new TCanvas("cptratio", "cptratio", 50, 50, 600, 600);
    cptratio->cd();
    phipiratio->Draw("P");
    //phikratio->Draw("SAME P");
    phipiratio_etacut->Draw("SAME P");
    //phikratio->Draw("SAME P");

    TCanvas* cpphiratio = new TCanvas("cpphiratio", "cpphiratio", 50, 50, 600, 600);
    cpphiratio->cd();
    pphiratio->Draw("P");

    //test of straight phi/h ratio (h = pi+k+p)
    Float_t numphi_etacut = phipt_etacut->Integral(phipt_etacut->FindBin(2.001), phipt_etacut->FindBin(3.999));
    Float_t numh_etacut =  pipt_etacut->Integral(pipt_etacut->FindBin(2.001), pipt_etacut->FindBin(3.999));
    //numh_etacut +=  kpt_etacut->Integral(kpt_etacut->FindBin(2.001), kpt_etacut->FindBin(3.999));
    //numh_etacut +=  ppt_etacut->Integral(ppt_etacut->FindBin(2.001), ppt_etacut->FindBin(3.999));
    Float_t ratio_etacut = numphi_etacut/numh_etacut;
    printf("integrated ratio phi/pi with eta cut: %f\n", ratio_etacut);

    //test of straight phi/h ratio (h = pi+k+p) in rapidity range -0.5 to 0.5
    Float_t numphi = phipt->Integral(phipt->FindBin(2.001), phipt->FindBin(3.999));
    Float_t numh =  pipt->Integral(pipt->FindBin(2.001), pipt->FindBin(3.999));
    numh +=  kpt->Integral(kpt->FindBin(2.001), kpt->FindBin(3.999));
    numh +=  ppt->Integral(ppt->FindBin(2.001), ppt->FindBin(3.999));
    Float_t ratio = numphi/numh;
    printf("integrated ratio w/ rapidity cut: %f\n", ratio);

    //test of straight phi/h ratio (h = pi+k+p) in triggered events
    Float_t trignumphi_etacut = trigphipt_etacut->Integral(trigphipt_etacut->FindBin(2.001), trigphipt_etacut->FindBin(3.999));
    Float_t trignumh_etacut =  trigpipt_etacut->Integral(trigpipt_etacut->FindBin(2.001), trigpipt_etacut->FindBin(3.999));
    //trignumh_etacut +=  trigkpt_etacut->Integral(trigkpt_etacut->FindBin(2.001), trigkpt_etacut->FindBin(3.999));
    //trignumh_etacut +=  trigppt_etacut->Integral(trigppt_etacut->FindBin(2.001), trigppt_etacut->FindBin(3.999));
    Float_t trigratio_etacut = trignumphi_etacut/trignumh_etacut;
    printf("triggered integrated phi/pi ratio w/ eta cut: %f\n", trigratio_etacut);

    //test of straight phi/h ratio (h = pi+k+p) in triggered events with rapidity range -0.5 to 0.5
    Float_t trignumphi = trigphipt->Integral(trigphipt->FindBin(2.001), trigphipt->FindBin(3.999));
    Float_t trignumh =  trigpipt->Integral(trigpipt->FindBin(2.001), trigpipt->FindBin(3.999));
    trignumh +=  trigkpt->Integral(trigkpt->FindBin(2.001), trigkpt->FindBin(3.999));
    trignumh +=  trigppt->Integral(trigppt->FindBin(2.001), trigppt->FindBin(3.999));
    Float_t trigratio = trignumphi/trignumh;
    printf("triggered integrated ratio w/ rapidity cut: %f\n", trigratio);


    printf("triggered phi/h ratio in low pT: %f, high pT: %f\n", trigphipt_etacut->Integral(trigphipt_etacut->FindBin(1.5001), trigphipt_etacut->FindBin(2.499))/trighpt->Integral(trighpt->FindBin(1.50001), trighpt->FindBin(2.499)), trigphipt_etacut->Integral(trigphipt_etacut->FindBin(2.5001), trigphipt_etacut->FindBin(3.999))/trighpt->Integral(trighpt->FindBin(2.50001), trighpt->FindBin(3.999)));

    printf("phi/h ratio in low pT: %f, high pT: %f\n", phipt_etacut->Integral(phipt_etacut->FindBin(1.5001), phipt_etacut->FindBin(2.499))/hpt->Integral(hpt->FindBin(1.50001), hpt->FindBin(2.499)), phipt_etacut->Integral(phipt_etacut->FindBin(2.5001), phipt_etacut->FindBin(3.999))/hpt->Integral(hpt->FindBin(2.50001), hpt->FindBin(3.999)));


    //h-pion correlations
    TFile* hpifile = new TFile("~/phiStudies/LHC18j2_FAST/trig_4.0_8.0_assoc_2.0_4.0_MC_hpi_0_100.root");
    TH2D* hpi2D = (TH2D*)hpifile->Get("hh2D");
    hpi2D->SetName("hpi2D");
    TH1D* hpidphi = (TH1D*)hpi2D->ProjectionY("hpidphi", hpi2D->GetXaxis()->FindBin(-1.199), hpi2D->GetXaxis()->FindBin(1.199));

    Float_t hpiUELine = 0.25*(hpidphi->GetBinContent(1) + hpidphi->GetBinContent(8) +hpidphi->GetBinContent(9) + hpidphi->GetBinContent(16));
    Float_t hpiNear = hpidphi->Integral(1, 8) - hpiUELine*8.0;
    Float_t hpiAway = hpidphi->Integral(9, 16) - hpiUELine*8.0;
    Float_t hpiTot = hpidphi->Integral(1, 16);
    Float_t hpiUE = hpiUELine*16.0;
   
    //h-p correlations
    TFile* hpfile = new TFile("~/phiStudies/LHC18j2_FAST/trig_4.0_8.0_assoc_2.0_4.0_hp_0_100.root");
    TH2D* hp2D = (TH2D*)hpfile->Get("hh2D");
    hp2D->SetName("hp2D");
    TH1D* hpdphi = (TH1D*)hp2D->ProjectionY("hpdphi", hp2D->GetXaxis()->FindBin(-1.199), hp2D->GetXaxis()->FindBin(1.199));

    Float_t hpUELine = 0.25*(hpdphi->GetBinContent(1) + hpdphi->GetBinContent(8) +hpdphi->GetBinContent(9) + hpdphi->GetBinContent(16));
    Float_t hpNear = hpdphi->Integral(1, 8) - hpUELine*8.0;
    Float_t hpAway = hpdphi->Integral(9, 16) - hpUELine*8.0;
    Float_t hpTot = hpdphi->Integral(1, 16);
    Float_t hpUE = hpUELine*16.0;

    //h-h correlations p8
    //TFile* hhfile = new TFile("~/phiStudies/LHC16k5_MC/trig_4_8_assoc_2_4_hh_16k5a_0_100.root");
    //2-4 pt range
//    TFile* hhfile = new TFile("~/phiStudies/LHC16k5_MC/trig_4_8_assoc_2_4_hh_09_29__0_100.root");
    //1.5-2.5 range
//    TFile* hhfile = new TFile("~/phiStudies/LHC16k5_MC/trig_4_8_assoc_1_2_hh_pythia8_0_100.root");
    //2.5-4.0 range
//    TFile* hhfile = new TFile("~/phiStudies/LHC16k5_MC/trig_4_8_assoc_2_4_hh_pythia8_0_100.root");

    //LF train runs 708 & 709
    //2-4 pt range
//    TFile* hhfile = new TFile("~/phiStudies/LFtrain/708_709/trig_4_8_assoc_2.0_4.0_hh_pythia8_0_100.root");
    //1.5-2.5 range
//    TFile* hhfile = new TFile("~/phiStudies/LFtrain/708_709/trig_4_8_assoc_1_2_hh_pythia820_80.root");
    TFile* hhfile = new TFile("~/phiStudies/LFtrain/708_709/trig_4_8_assoc_1_2_hh_pythia80_100.root");
//    TFile* hhfile = new TFile("~/phiStudies/LFtrain/708_709/trig_4_8_assoc_1_2_hh_pythia80_20.root");
    //2.5-4.0 range
//    TFile* hhfile = new TFile("~/phiStudies/LFtrain/708_709/trig_4_8_assoc_2.5_4.0_hh_pythia80_100.root");
//    TFile* hhfile = new TFile("~/phiStudies/LFtrain/708_709/trig_4_8_assoc_2.5_4.0_hh_pythia820_80.root");
//    TFile* hhfile = new TFile("~/phiStudies/LFtrain/708_709/trig_4_8_assoc_2.5_4.0_hh_pythia80_20.root");

    TH2D* hh2D = (TH2D*)hhfile->Get("hh2D");
    TH1D* hhdphi = (TH1D*)hh2D->ProjectionY("hhdphi", hh2D->GetXaxis()->FindBin(-1.199), hh2D->GetXaxis()->FindBin(1.199))->Clone("hhdphi");
    hhdphi->SetMarkerStyle(28);
    hhdphi->SetMarkerColor(kBlue+1);
    hhdphi->SetMarkerSize(2);
    hhdphi->SetLineWidth(2);
    hhdphi->SetLineColor(kBlue+1);
    hhdphi->GetYaxis()->SetTitle("Per Trigger (h-h) pairs");
    hhdphi->GetYaxis()->SetMaxDigits(2);
    hhdphi->SetStats(kFALSE);
    hhdphi->SetTitle("h-h #Delta#varphi Correlations");
    hhdphi->GetXaxis()->SetTitle("#Delta#varphi");

    Float_t hhbgerr = 0.0;
    Float_t hhbg = ZYAM(hhdphi, hhbgerr);

    Double_t hhNearError = 0, hhAwayError = 0, hhTotError = 0, hhJetError = 0, hhUEerror = 0;

    Double_t hhUELine = 0.25*(hhdphi->GetBinContent(1) + hhdphi->GetBinContent(8) +hhdphi->GetBinContent(9) + hhdphi->GetBinContent(16));
    //hhbg=hhUELine;
    Double_t hhNear = calcJetFromBG(hhdphi, 1, 8, hhbg, hhNearError);
    Double_t hhAway = calcJetFromBG(hhdphi, 9, 16, hhbg, hhAwayError);
    Double_t hhTot = hhdphi->IntegralAndError(1, 16, hhTotError);
    Double_t hhUE = calcUE(hhdphi, 1, 16, hhbg, hhUEerror);
    Double_t hhJet = hhNear + hhAway;
    hhJetError = TMath::Sqrt(TMath::Power(hhNearError, 2) + TMath::Power(hhAwayError, 2));

    TF1* hhdphiUE = new TF1("hhdphiUE", "pol0(0)", -TMath::Pi()/2.0, 1.5*TMath::Pi());
    hhdphiUE->SetParameter(0, hhbg/(2.0*2.4*TMath::Pi()/16.0));
    hhdphiUE->SetLineStyle(9);
    hhdphiUE->SetLineWidth(4);
    hhdphiUE->SetLineColor(kBlue);

    TCanvas* chh = new TCanvas("chh", "chh", 50, 50, 600, 600);
    chh->cd();
    hhdphi->Scale(1.0/(2.4*2.0*TMath::Pi()/16.0));
    hhdphi->Draw();
    hhdphiUE->Draw("SAME");

    //h-phi correlations p8
    //TFile* hphifile = new TFile("~/phiStudies/LHC16k5_MC/trig_4_8_assoc_2_4_hphi_16k5a_0_100.root");
    //2-4 pt range
//    TFile* hphifile = new TFile("~/phiStudies/LHC16k5_MC/trig_4_8_assoc_2_4_MC_hphi_09_29__0_100.root");
    //1.5-2.5 pt range
//    TFile* hphifile = new TFile("~/phiStudies/LHC16k5_MC/trig_4_8_assoc_1_2_MC_hphi_pythia8_0_100.root");
    //2.5-4.0 pt range
//    TFile* hphifile = new TFile("~/phiStudies/LHC16k5_MC/trig_4_8_assoc_2_4_MC_hphi_pythia8_0_100.root");

    //LF train runs 708 & 709
    //2-4 pt range
//    TFile* hphifile = new TFile("~/phiStudies/LFtrain/708_709/trig_4_8_assoc_2.0_4.0_MC_hphi_pythia8_0_100.root");
    //1.5-2.5 pt range
//    TFile* hphifile = new TFile("~/phiStudies/LFtrain/708_709/trig_4_8_assoc_1_2_MC_hphi_pythia820_80.root");
//    TFile* hphifile = new TFile("~/phiStudies/LFtrain/708_709/trig_4_8_assoc_1_2_MC_hphi_pythia80_20.root");
    TFile* hphifile = new TFile("~/phiStudies/LFtrain/708_709/trig_4_8_assoc_1_2_MC_hphi_pythia80_100.root");
    //2.5-4.0 pt range
//    TFile* hphifile = new TFile("~/phiStudies/LFtrain/708_709/trig_4_8_assoc_2.5_4.0_MC_hphi_pythia80_100.root");
//    TFile* hphifile = new TFile("~/phiStudies/LFtrain/708_709/trig_4_8_assoc_2.5_4.0_MC_hphi_pythia820_80.root");
//    TFile* hphifile = new TFile("~/phiStudies/LFtrain/708_709/trig_4_8_assoc_2.5_4.0_MC_hphi_pythia80_20.root");

    TH2D* hphi2D = (TH2D*)hphifile->Get("MChPhi2D");
    TH1D* hphidphi = (TH1D*)hphi2D->ProjectionY("hphidphi", hphi2D->GetXaxis()->FindBin(-1.199), hphi2D->GetXaxis()->FindBin(1.199));
    hphidphi->SetMarkerStyle(34);
    hphidphi->SetMarkerColor(kBlue+1);
    hphidphi->SetMarkerSize(2);
    hphidphi->SetLineColor(kBlue+1);
    hphidphi->SetLineWidth(2);
    hphidphi->GetYaxis()->SetTitle("Per Trigger (h-#phi) pairs");
    hphidphi->GetYaxis()->SetMaxDigits(2);
    hphidphi->SetStats(kFALSE);
    hphidphi->GetXaxis()->SetTitle("#Delta#varphi");
    hphidphi->SetTitle("h-#phi #Delta#varphi Correlations");

    TF1 *hphifit = new TF1("hphifit", "[0] + [1]*TMath::Exp(TMath::Power(x - [2], 2)/[3]) + [4]*TMath::Exp(TMath::Power(x-[5], 2)/[6])", -TMath::Pi()/2.0, 1.5*TMath::Pi());
    hphifit->SetParameters(0.008, 0.005, 0.0, 0.0005, 0.0005, 3.14, 0.0005);
    hphifit->SetParLimits(0, 0, 0.01);
    hphifit->SetParLimits(1, 0.0001, 0.02);
    //hphidphi->Fit(hphifit);

    Float_t hphibgerr = 0.0;
    Float_t hphibg = ZYAM(hphidphi, hphibgerr);

    Double_t hphiNearError = 0, hphiAwayError = 0, hphiTotError = 0, hphiJetError = 0, hphiUEerror=0;
    Double_t hphiUELine = 0.25*(hphidphi->GetBinContent(1) + hphidphi->GetBinContent(8) +hphidphi->GetBinContent(9) + hphidphi->GetBinContent(16));
    //hphibg=hphiUELine;
    Double_t hphiNear = calcJetFromBG(hphidphi, 1, 8, hphibg, hphiNearError);
    Double_t hphiAway = calcJetFromBG(hphidphi, 9, 16, hphibg, hphiAwayError);
    Double_t hphiTot = hphidphi->IntegralAndError(1, 16, hphiTotError);
    Double_t hphiUE = calcUE(hphidphi, 1, 16, hphibg, hphiUEerror);
    Double_t hphiJet = hphiNear+hphiAway;
    hphiJetError = TMath::Sqrt(TMath::Power(hphiNearError, 2) + TMath::Power(hphiAwayError, 2));
    
    //Calculate ratio errors for h-phi/h-h:
    Double_t nearError = (hphiNear/hhNear)*TMath::Sqrt(TMath::Power(hphiNearError/hphiNear, 2) + TMath::Power(hhNearError/hhNear, 2));
    Double_t awayError = (hphiAway/hhAway)*TMath::Sqrt(TMath::Power(hphiAwayError/hphiAway, 2) + TMath::Power(hhAwayError/hhAway, 2));
    Double_t totError = (hphiTot/hhTot)*TMath::Sqrt(TMath::Power(hphiTotError/hphiTot, 2) + TMath::Power(hhTotError/hhTot, 2));
    Double_t jetError = (hphiJet/hhJet)*TMath::Sqrt(TMath::Power(hphiJetError/hphiJet, 2) + TMath::Power(hhJetError/hhJet, 2));

    TF1* hphidphiUE = new TF1("hphidphiUE", "pol0(0)", -TMath::Pi()/2.0, 1.5*TMath::Pi());
    hphidphiUE->SetParameter(0, hphiUELine/(2.4*2.0*TMath::Pi()/16.0));
    hphidphiUE->SetLineStyle(9);
    hphidphiUE->SetLineWidth(4);
    hphidphiUE->SetLineColor(kBlue);

    TF1* hphidphiZYAM = new TF1("hphidphiZYAM", "pol0(0)", -TMath::Pi()/2.0, 1.5*TMath::Pi());
    hphidphiZYAM->SetParameter(0, hphibg/(2.4*2.0*TMath::Pi()/16.0));
    hphidphiZYAM->SetLineStyle(9);
    hphidphiZYAM->SetLineWidth(4);
    hphidphiZYAM->SetLineColor(kBlue);

    TCanvas* c1 = new TCanvas("c1", "c1", 50, 50, 600, 600);
    c1->cd();
    hphidphi->Scale(1.0/(2.4*2.0*TMath::Pi()/16.0));
    hphidphi->Draw();
    hphidphiUE->Draw("SAME");
    hphidphiZYAM->Draw("SAME");
   
    printf("pythia8\n");
    printf("h-phi/h-h  near: %f, away: %f, UE: %f, Total: %f\n", hphiNear/hhNear, hphiAway/hhAway, hphiUE/hhUE, hphiTot/hhTot);
    printf("errors near: %f, away: %f, UE %f, Total: %f\n", nearError, awayError, 0.0, totError);
//    printf("h-phi/h-p  near: %f, away: %f, UE: %f, Total: %f\n", hphiNear/hpNear, hphiAway/hpAway, hphiUE/hpUE, hphiTot/hpTot);
    printf("h-phi/h-pi  near: %f, away: %f, UE: %f, Total: %f\n", hphiNear/hpiNear, hphiAway/hpiAway, hphiUE/hpiUE, hphiTot/hpiTot);
    printf("per trigger yields\n");
    printf("h-phi (x250) near: %f, away: %f, UE: %f, Total: %f\n", hphiNear*250., hphiAway*250., hphiUE*250., hphiTot*250.);
    printf("h-phi (x250)errors  near: %f, away: %f, UE: %f, Total: %f\n", hphiNearError*250., hphiAwayError*250., hphiUEerror*250., hphiTotError*250.);
    printf("h-h near: %f, away: %f, UE: %f, Total: %f\n", hhNear, hhAway, hhUE, hhTot);
    printf("h-h errors near: %f, away: %f, UE: %f, Total: %f\n", hhNearError, hhAwayError, hhUEerror, hhTotError);
//    printf("h-p/h-pi  near: %f, away: %f, UE: %f, Total: %f\n", hpNear/hpiNear, hpAway/hpiAway, hpUE/hpiUE, hpTot/hpiTot);
 
    printf("h-h Jet/Tot = %f\nh-phi Jet/Tot = %f\n", hhJet/hhTot, hphiJet/hphiTot);

    //h-h correlations p6
    TFile* hh6file = new TFile("~/phiStudies/LHC16k5_MC/trig_4_8_assoc_2_4_hh_16k5b_0_100.root");
    TH2D* hh62D = (TH2D*)hh6file->Get("hh2D");
    TH1D* hh6dphi = (TH1D*)hh62D->ProjectionY("hh6dphi", hh62D->GetXaxis()->FindBin(-1.199), hh62D->GetXaxis()->FindBin(1.199))->Clone("hh6dphi");
    hh6dphi->SetMarkerStyle(28);
    hh6dphi->SetMarkerColor(kGreen+1);
    hh6dphi->SetMarkerSize(2);
    hh6dphi->SetLineWidth(2);
    hh6dphi->SetLineColor(kBlue+1);
    hh6dphi->GetYaxis()->SetTitle("Per Trigger (h-h) pairs");
    hh6dphi->GetYaxis()->SetMaxDigits(2);
    hh6dphi->SetStats(kFALSE);
    hh6dphi->SetTitle("h-h #Delta#varphi Correlations");
    hh6dphi->GetXaxis()->SetTitle("#Delta#vaprhi");

    
    Double_t hh6NearError = 0, hh6AwayError = 0, hh6TotError = 0, hh6JetError = 0;

    Double_t hh6UELine = 0.25*(hh6dphi->GetBinContent(1) + hh6dphi->GetBinContent(8) +hh6dphi->GetBinContent(9) + hh6dphi->GetBinContent(16));
    Double_t hh6Near = hh6dphi->IntegralAndError(1, 8, hh6NearError) - hh6UELine*8.0;
    Double_t hh6Away = hh6dphi->IntegralAndError(9, 16, hh6AwayError) - hh6UELine*8.0;
    Double_t hh6Tot = hh6dphi->IntegralAndError(1, 16, hh6TotError);
    Double_t hh6UE = hh6UELine*16.0;
    Double_t hh6Jet = hh6Near + hh6Away;
    hh6JetError = TMath::Sqrt(TMath::Power(hh6NearError, 2) + TMath::Power(hh6AwayError, 2));

    TF1* hh6dphiUE = new TF1("hh6dphiUE", "pol0(0)", -TMath::Pi()/2.0, 1.5*TMath::Pi());
    hh6dphiUE->SetParameter(0, hh6UELine);
    hh6dphiUE->SetLineStyle(9);
    hh6dphiUE->SetLineWidth(4);
    hh6dphiUE->SetLineColor(kGreen);

    chh->cd();
    //hh6dphi->Draw("SAME");
    //hh6dphiUE->Draw("SAME");

    //h-phi correlations p6
    TFile* hphi6file = new TFile("~/phiStudies/LHC16k5_MC/trig_4_8_assoc_2_4_hphi_16k5b_0_100.root");
    TH2D* hphi62D = (TH2D*)hphi6file->Get("MChPhi2D");
    TH1D* hphi6dphi = (TH1D*)hphi62D->ProjectionY("hphi6dphi", hphi62D->GetXaxis()->FindBin(-1.199), hphi62D->GetXaxis()->FindBin(1.199));
    hphi6dphi->SetMarkerStyle(34);
    hphi6dphi->SetMarkerColor(kGreen+1);
    hphi6dphi->SetMarkerSize(2);
    hphi6dphi->SetLineColor(kBlue+1);
    hphi6dphi->SetLineWidth(2);
    hphi6dphi->GetYaxis()->SetTitle("Per Trigger (h-#phi) pairs");
    hphi6dphi->GetYaxis()->SetMaxDigits(2);
    hphi6dphi->SetStats(kFALSE);
    hphi6dphi->GetXaxis()->SetTitle("#Delta#varphi");
    hphi6dphi->SetTitle("h-#phi #Delta#varphi Correlations");

    TF1 *hphi6fit = new TF1("hphi6fit", "[0] + [1]*TMath::Exp(TMath::Power(x - [2], 2)/[3]) + [4]*TMath::Exp(TMath::Power(x-[5], 2)/[6])", -TMath::Pi()/2.0, 1.5*TMath::Pi());
    hphi6fit->SetParameters(0.008, 0.005, 0.0, 0.0005, 0.0005, 3.14, 0.0005);
    hphi6fit->SetParLimits(0, 0, 0.01);
    hphi6fit->SetParLimits(1, 0.0001, 0.02);
    //hphi6dphi->Fit(hphi6fit);

    Double_t hphi6NearError = 0, hphi6AwayError = 0, hphi6TotError = 0, hphi6JetError = 0;
    Double_t hphi6UELine = 0.25*(hphi6dphi->GetBinContent(1) + hphi6dphi->GetBinContent(8) +hphi6dphi->GetBinContent(9) + hphi6dphi->GetBinContent(16));
    Double_t hphi6Near = hphi6dphi->IntegralAndError(1, 8, hphi6NearError) - hphi6UELine*8.0;
    Double_t hphi6Away = hphi6dphi->IntegralAndError(9, 16, hphi6AwayError) - hphi6UELine*8.0;
    Double_t hphi6Tot = hphi6dphi->IntegralAndError(1, 16, hphi6TotError);
    Double_t hphi6UE = hphi6UELine*16.0;
    Double_t hphi6Jet = hphi6Near+hphi6Away;
    hphi6JetError = TMath::Sqrt(TMath::Power(hphi6NearError, 2) + TMath::Power(hphi6AwayError, 2));

    TF1* hphi6dphiUE = new TF1("hphi6dphiUE", "pol0(0)", -TMath::Pi()/2.0, 1.5*TMath::Pi());
    hphi6dphiUE->SetParameter(0, hphi6UELine);
    hphi6dphiUE->SetLineStyle(9);
    hphi6dphiUE->SetLineWidth(4);
    hphi6dphiUE->SetLineColor(kGreen);

    TLegend *hphileg = new TLegend(0.5, 0.8, 0.6, 0.9);
    //hphileg->AddEntry(hphi6dphi, "Pythia 6 Perguia2011", "lp");
    hphileg->AddEntry(hphidphi, "Pythia 8 Monash", "lp");

    c1->cd();
    //hphi6dphi->Draw("SAME");
    //hphi6dphiUE->Draw("SAME");
    hphileg->Draw("SAME");

    printf("pythia6\n");
    printf("h-phi/h-h  near: %f, away: %f, UE: %f, Total: %f\n", hphi6Near/hh6Near, hphi6Away/hh6Away, hphi6UE/hh6UE, hphi6Tot/hh6Tot);


}
