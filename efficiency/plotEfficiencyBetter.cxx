void plotMultEff(THnSparse* reco, THnSparse* real, TH1D** recoVar, TH1D** realVar, TH1D** eff, TH1D** ratio, Float_t* mult, Int_t numMultBins, Int_t axis, TCanvas* ceff, TCanvas* cratio, Bool_t isSingle = kFALSE){

    TString var = "";
    TString label = "";
    Int_t rebin = 1;
    Int_t numbins = 10;
    Double_t binedge[11] = {0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 5.0, 7.0, 10.0};

    Int_t colors[] = {kGreen+4,kGreen-2,kSpring+5,kRed+1};
    
    Int_t multAxis = 6;
    if(isSingle) multAxis = 5;

    printf("num mult bins: %i\n", numMultBins);
    switch(axis){
        case 0: var = "PT";
                label = "p_{T}";
                //rebin = 5;
                break;
        case 1: var = "Phi";
                label = "#varphi";
                rebin = 4;
                break;
        case 2: var = "Eta";
                label = "#eta";
                rebin = 2;
                break;
        case 3: var = "Y";
                label = "y";
                rebin = 2;
                break;
        case 4: var = "Z";
                label = "Z Vertex (cm)";
                break;
        default: var = "Var";
                 label = "var";
                 break;
    }

    for(int imult = 0; imult < numMultBins; imult++){
        reco->GetAxis(multAxis)->SetRangeUser(mult[imult], mult[imult+1]);
        real->GetAxis(multAxis)->SetRangeUser(mult[imult], mult[imult+1]);
        recoVar[imult] = reco->Projection(axis);
        recoVar[imult]->SetStats(kFALSE);
        recoVar[imult]->SetName(Form("%sreco%s_mult%i%i", reco->GetTitle(), var.Data(), int(mult[imult]), int(mult[imult+1])));
        realVar[imult] = real->Projection(axis);
        realVar[imult]->SetName(Form("%sreco%s_mult%i%i", real->GetTitle(), var.Data(), int(mult[imult]), int(mult[imult+1])));
        if(axis!=0){
            recoVar[imult]->Rebin(rebin);
            realVar[imult]->Rebin(rebin);
        }else{
            recoVar[imult]= (TH1D*)recoVar[imult]->Rebin(numbins, Form("%sreco_rebin%s_mult%i%i", reco->GetTitle(), var.Data(), int(mult[imult]), int(mult[imult+1])), binedge);
            realVar[imult]= (TH1D*)realVar[imult]->Rebin(numbins, Form("%sreal_rebin%s_mult%i%i", real->GetTitle(), var.Data(), int(mult[imult]), int(mult[imult+1])), binedge);
        } 
        eff[imult] = (TH1D*)recoVar[imult]->Clone(Form("%s_%s_eff%s_mult%i%i", reco->GetTitle(), real->GetTitle(), var.Data(), int(mult[imult]), int(mult[imult+1])));
        eff[imult]->Divide(eff[imult], realVar[imult], 1., 1., "B");
        eff[imult]->SetTitle(Form("(%s)/(%s) Efficiency vs. %s for Mult. Bins;%s;Efficiency", reco->GetTitle(), real->GetTitle(), label.Data(), label.Data()));
        eff[imult]->SetLineColor(colors[imult]);
        eff[imult]->SetMarkerColor(colors[imult]);
        eff[imult]->SetMarkerSize(1);
        eff[imult]->SetMarkerStyle(22);
        eff[imult]->GetYaxis()->SetRangeUser(0.0, eff[imult]->GetBinContent(eff[imult]->GetMaximumBin())*1.25);
    }
 
    reco->GetAxis(multAxis)->SetRangeUser(0.0, 90.0);
    real->GetAxis(multAxis)->SetRangeUser(0.0, 90.0);
    recoVar[numMultBins] = reco->Projection(axis);
    recoVar[numMultBins]->SetName(Form("%sreco%s_mult090", reco->GetTitle(), var.Data()));
    realVar[numMultBins] = real->Projection(axis);
    realVar[numMultBins]->SetName(Form("%sreal%s_mult090", real->GetTitle(), var.Data()));
    if(axis!=0){
        recoVar[numMultBins]->Rebin(rebin);
        realVar[numMultBins]->Rebin(rebin);
    }else{
        recoVar[numMultBins]= (TH1D*)recoVar[numMultBins]->Rebin(numbins, Form("%sreco_rebin%s_mult090", reco->GetTitle(), var.Data()), binedge);
        realVar[numMultBins]= (TH1D*)realVar[numMultBins]->Rebin(numbins, Form("%sreal_rebin%s_mult090", real->GetTitle(), var.Data()), binedge);

    }
    eff[numMultBins] = (TH1D*)recoVar[numMultBins]->Clone(Form("%s_%s_eff%s_mult090",reco->GetTitle(), real->GetTitle(), var.Data()));
    eff[numMultBins]->Divide(eff[numMultBins], realVar[numMultBins], 1., 1., "B");
    eff[numMultBins]->SetTitle(Form("(%s)/(%s) Efficiency vs. %s for Mult. Bins", reco->GetTitle(), real->GetTitle(), label.Data()));
    eff[numMultBins]->SetLineColor(kRed+1);
    eff[numMultBins]->SetMarkerColor(kRed+1);
    eff[numMultBins]->SetMarkerSize(1);
    eff[numMultBins]->SetMarkerStyle(22);

    cratio->SetLeftMargin(0.13);
    cratio->cd();
    for(int i=0; i<numMultBins; i++){
        ratio[i] = (TH1D*)eff[i]->Clone(Form("%s_%s_ratio%s_mult%i", reco->GetTitle(), real->GetTitle(), var.Data(), i));
        ratio[i]->SetTitle(Form("(%s)/(%s) Efficiency vs %s Ratio to 0-100%% Efficiency;%s;Ratio", reco->GetTitle(), real->GetTitle(), label.Data(), label.Data()));
        ratio[i]->Divide(ratio[i], eff[3], 1., 1., "B");
        ratio[i]->GetYaxis()->SetRangeUser(0.75, 1.25);
        ratio[i]->Draw("P SAME");
    } 
    //cratio->Print(Form("plots/binom_%s_%s_%sratio.pdf", var.Data(), reco->GetName(), real->GetName()));

    ceff->SetLeftMargin(0.13);
    ceff->cd();
    for(int j=0; j<=numMultBins; j++){
        eff[j]->Draw("P SAME");
    }
    //ceff->Print(Form("plots/binom_%s_%s_%seff.pdf", var.Data(), reco->GetName(), real->GetName()));
    TLegend* effleg = new TLegend(0.2, 0.6, 0.6, 0.8);
    effleg->AddEntry(eff[3], "0-100% Mult. Percentile", "lpe");
    effleg->AddEntry(eff[0], "0-20% Mult. Percentile", "lpe");
    effleg->AddEntry(eff[1], "20-50% Mult. Percentile", "lpe");
    effleg->AddEntry(eff[2], "50-80% Mult. Percentile", "lpe");
    effleg->Draw();

}; 

void plotEfficiencyBetter(TString input){

    Float_t mult[4] = {0.0, 20.0, 50.0, 80.0};
    
    TFile* infile = new TFile(input.Data());
    TList* list = (TList*)infile->Get("phiEff_mult_0_100");
    
    //single track histos
    THnSparseF* realCharged = (THnSparseF*)list->FindObject("fRealChargedDist");
    THnSparseF* recoCharged = (THnSparseF*)list->FindObject("fRecoChargedDist");
    THnSparseF* recop = (THnSparseF*)list->FindObject("fRecopDist");
    THnSparseF* realp = (THnSparseF*)list->FindObject("fRealpDist");
    THnSparseF* recoPi = (THnSparseF*)list->FindObject("fRecoPiDist");
    THnSparseF* realPi = (THnSparseF*)list->FindObject("fRealPiDist");
    THnSparseF* recoK = (THnSparseF*)list->FindObject("fRecoKDist");
    THnSparseF* realK = (THnSparseF*)list->FindObject("fRealKDist");
    THnSparseF* recoe = (THnSparseF*)list->FindObject("fRecoeDist");
    THnSparseF* reale = (THnSparseF*)list->FindObject("fRealeDist");
    THnSparseF* recoMuon = (THnSparseF*)list->FindObject("fRecoMuonDist");
    THnSparseF* realMuon = (THnSparseF*)list->FindObject("fRealMuonDist");

    realCharged->GetAxis(0)->SetRangeUser(0.5, 10.0);
    recoCharged->GetAxis(0)->SetRangeUser(0.5, 10.0);
    recop->GetAxis(0)->SetRangeUser(0.5, 10.0);
    realp->GetAxis(0)->SetRangeUser(0.5, 10.0);
    recoPi->GetAxis(0)->SetRangeUser(0.5, 10.0);
    realPi->GetAxis(0)->SetRangeUser(0.5, 10.0);
    recoK->GetAxis(0)->SetRangeUser(0.5, 10.0);
    realK->GetAxis(0)->SetRangeUser(0.5, 10.0);
    recoe->GetAxis(0)->SetRangeUser(0.5, 10.0);
    reale->GetAxis(0)->SetRangeUser(0.5, 10.0);
    recoMuon->GetAxis(0)->SetRangeUser(0.5, 10.0);
    realMuon->GetAxis(0)->SetRangeUser(0.5, 10.0);

    recoCharged->GetAxis(2)->SetRangeUser(-0.8, 0.8);
    realCharged->GetAxis(2)->SetRangeUser(-0.8, 0.8);
    recop->GetAxis(2)->SetRangeUser(-0.8, 0.8);
    realp->GetAxis(2)->SetRangeUser(-0.8, 0.8);
    recoPi->GetAxis(2)->SetRangeUser(-0.8, 0.8);
    realPi->GetAxis(2)->SetRangeUser(-0.8, 0.8);
    recoK->GetAxis(2)->SetRangeUser(-0.8, 0.8);
    realK->GetAxis(2)->SetRangeUser(-0.8, 0.8);
    recoe->GetAxis(2)->SetRangeUser(-0.8, 0.8);
    reale->GetAxis(2)->SetRangeUser(-0.8, 0.8);
    recoMuon->GetAxis(2)->SetRangeUser(-0.8, 0.8);
    realMuon->GetAxis(2)->SetRangeUser(-0.8, 0.8);

    recoCharged->GetAxis(4)->SetRangeUser(-10.0, 10.0);
    realCharged->GetAxis(4)->SetRangeUser(-10.0, 10.0);
    recop->GetAxis(4)->SetRangeUser(-10.0, 10.0);
    realp->GetAxis(4)->SetRangeUser(-10.0, 10.0);
    recoPi->GetAxis(4)->SetRangeUser(-10.0, 10.0);
    realPi->GetAxis(4)->SetRangeUser(-10.0, 10.0);
    recoK->GetAxis(4)->SetRangeUser(-10.0, 10.0);
    realK->GetAxis(4)->SetRangeUser(-10.0, 10.0);
    recoe->GetAxis(4)->SetRangeUser(-10.0, 10.0);
    reale->GetAxis(4)->SetRangeUser(-10.0, 10.0);
    recoMuon->GetAxis(4)->SetRangeUser(-10.0, 10.0);
    realMuon->GetAxis(4)->SetRangeUser(-10.0, 10.0);

    recoCharged->Sumw2();
    realCharged->Sumw2();
    recop->Sumw2();
    realp->Sumw2();
    recoPi->Sumw2();
    realPi->Sumw2();
    recoK->Sumw2();
    realK->Sumw2();
    recoe->Sumw2();
    reale->Sumw2();
    recoMuon->Sumw2();
    realMuon->Sumw2();
    
    TH1D* recoCharged_PT_mult[4];
    TH1D* realCharged_PT_mult[4];
    TH1D* recop_PT_mult[4];
    TH1D* realp_PT_mult[4];
    TH1D* recoPi_PT_mult[4];
    TH1D* realPi_PT_mult[4];
    TH1D* recoK_PT_mult[4];
    TH1D* realK_PT_mult[4];
    TH1D* recoe_PT_mult[4];
    TH1D* reale_PT_mult[4];
    TH1D* recoMuon_PT_mult[4];
    TH1D* realMuon_PT_mult[4];

    TH1D* effCharged_PT_mult[4];
    TH1D* effp_PT_mult[4];
    TH1D* effPi_PT_mult[4];
    TH1D* effK_PT_mult[4];
    TH1D* effe_PT_mult[4];
    TH1D* effMuon_PT_mult[4];

    TH1D* ratioCharged_PT_mult[3];
    TH1D* ratiop_PT_mult[3];
    TH1D* ratioPi_PT_mult[3];
    TH1D* ratioK_PT_mult[3];
    TH1D* ratioe_PT_mult[3];
    TH1D* ratioMuon_PT_mult[3];

    TCanvas* ceffCharged_PT = new TCanvas("ceffCharged_PT", "effCharged_PT", 50, 50, 600, 600);
    TCanvas* ceffp_PT = new TCanvas("ceffp_PT", "effp_PT", 50, 50, 600, 600);
    TCanvas* ceffPi_PT = new TCanvas("ceffPi_PT", "effPi_PT", 50, 50, 600, 600);
    TCanvas* ceffK_PT = new TCanvas("ceffK_PT", "effK_PT", 50, 50, 600, 600);
    TCanvas* ceffe_PT = new TCanvas("ceffe_PT", "effe_PT", 50, 50, 600, 600);
    TCanvas* ceffMuon_PT = new TCanvas("ceffMuon_PT", "effMuon_PT", 50, 50, 600, 600);

    TCanvas* cratioCharged_PT = new TCanvas("cratioCharged_PT", "ratioCharged_PT", 50, 50, 600, 600);
    TCanvas* cratiop_PT = new TCanvas("cratiop_PT", "ratiop_PT", 50, 50, 600, 600);
    TCanvas* cratioPi_PT = new TCanvas("cratioPi_PT", "ratioPi_PT", 50, 50, 600, 600);
    TCanvas* cratioK_PT = new TCanvas("cratioK_PT", "ratioK_PT", 50, 50, 600, 600);
    TCanvas* cratioe_PT = new TCanvas("cratioe_PT", "ratioe_PT", 50, 50, 600, 600);
    TCanvas* cratioMuon_PT = new TCanvas("cratioMuon_PT", "ratioMuon_PT", 50, 50, 600, 600);

    plotMultEff(recoCharged, realCharged, recoCharged_PT_mult, realCharged_PT_mult, effCharged_PT_mult, ratioCharged_PT_mult, mult, 3, 0, ceffCharged_PT, cratioCharged_PT, kTRUE);
    plotMultEff(recop, realp, recop_PT_mult, realp_PT_mult, effp_PT_mult, ratiop_PT_mult, mult, 3, 0, ceffp_PT, cratiop_PT, kTRUE);
    plotMultEff(recoPi, realPi, recoPi_PT_mult, realPi_PT_mult, effPi_PT_mult, ratioPi_PT_mult, mult, 3, 0, ceffPi_PT, cratioPi_PT, kTRUE);
    plotMultEff(recoK, realK, recoK_PT_mult, realK_PT_mult, effK_PT_mult, ratioK_PT_mult, mult, 3, 0, ceffK_PT, cratioK_PT, kTRUE);
    plotMultEff(recoe, reale, recoe_PT_mult, reale_PT_mult, effe_PT_mult, ratioe_PT_mult, mult, 3, 0, ceffe_PT, cratioe_PT, kTRUE);
    plotMultEff(recoMuon, realMuon, recoMuon_PT_mult, realMuon_PT_mult, effMuon_PT_mult, ratioMuon_PT_mult, mult, 3, 0, ceffMuon_PT, cratioMuon_PT, kTRUE);

    //trigger distribution efficiency histos 

    //phi(1020) histos
    THnSparseF* recoPhi = (THnSparseF*)list->FindObject("fRecoPhiDist");
    THnSparseF* realPhi = (THnSparseF*)list->FindObject("fRealPhiDist");
    THnSparseF* trackPhi = (THnSparseF*)list->FindObject("fTrackRecoPhiDist");
    THnSparseF* TOFPhi = (THnSparseF*)list->FindObject("fTOFRecoPhiDist");
    THnSparseF* TPCTrackPhi = (THnSparseF*)list->FindObject("fTPCPIDTrackRecoPhiDist");
    THnSparseF* TPCTOFPhi = (THnSparseF*)list->FindObject("fTPCPIDRecoPhiDist");
    THnSparseF* PIDPhi = (THnSparseF*)list->FindObject("fPIDRecoPhiDist");

    realPhi->SetTitle("real #phi");
    realPhi->SetName("real");
    recoPhi->SetTitle("reco");
    recoPhi->SetName("reco");
    trackPhi->SetTitle("track cuts");
    trackPhi->SetName("track");
    TOFPhi->SetTitle("track cuts + TOF hit");
    TOFPhi->SetName("TOF");
    TPCTrackPhi->SetTitle("track cuts + TPC PID");
    TPCTrackPhi->SetName("TPCtrack");
    TPCTOFPhi->SetTitle("track cuts + TOF hit + TPC PID");
    TPCTOFPhi->SetName("TPCTOF");
    PIDPhi->SetTitle("track cuts + TOF PID + TPC PID");
    PIDPhi->SetName("PID");
    
    //set standard "reconstructed" histo to be TOF hit + track cuts
    recoPhi = (THnSparseF*)TOFPhi->Clone("reco");

    
    recoPhi->GetAxis(0)->SetRangeUser(0.5, 10.0);
    realPhi->GetAxis(0)->SetRangeUser(0.5, 10.0);
    trackPhi->GetAxis(0)->SetRangeUser(0.5, 10.0);
    TOFPhi->GetAxis(0)->SetRangeUser(0.5, 10.0);
    TPCTrackPhi->GetAxis(0)->SetRangeUser(0.5, 10.0);
    TPCTOFPhi->GetAxis(0)->SetRangeUser(0.5, 10.0);
    PIDPhi->GetAxis(0)->SetRangeUser(0.5, 10.0);
       
    recoPhi->GetAxis(2)->SetRangeUser(-0.8, 0.8);
    realPhi->GetAxis(2)->SetRangeUser(-0.8, 0.8);
    trackPhi->GetAxis(2)->SetRangeUser(-0.8, 0.8);
    TOFPhi->GetAxis(2)->SetRangeUser(-0.8, 0.8);
    TPCTrackPhi->GetAxis(2)->SetRangeUser(-0.8, 0.8);
    TPCTOFPhi->GetAxis(2)->SetRangeUser(-0.8, 0.8);
    PIDPhi->GetAxis(2)->SetRangeUser(-0.8, 0.8);

    recoPhi->GetAxis(5)->SetRangeUser(1.01, 1.03);
    realPhi->GetAxis(5)->SetRangeUser(1.01, 1.03);
    trackPhi->GetAxis(5)->SetRangeUser(1.01, 1.03);
    TOFPhi->GetAxis(5)->SetRangeUser(1.01, 1.03);
    TPCTrackPhi->GetAxis(5)->SetRangeUser(1.01, 1.03);
    TPCTOFPhi->GetAxis(5)->SetRangeUser(1.01, 1.03);
    PIDPhi->GetAxis(5)->SetRangeUser(1.01, 1.03);

    recoPhi->GetAxis(4)->SetRangeUser(-10.0, 10.0);
    realPhi->GetAxis(4)->SetRangeUser(-10.0, 10.0);
    trackPhi->GetAxis(4)->SetRangeUser(-10.0, 10.0);
    TOFPhi->GetAxis(4)->SetRangeUser(-10.0, 10.0);
    TPCTrackPhi->GetAxis(4)->SetRangeUser(-10.0, 10.0);
    TPCTOFPhi->GetAxis(4)->SetRangeUser(-10.0, 10.0);
    PIDPhi->GetAxis(4)->SetRangeUser(-10.0, 10.0);

    recoPhi->Sumw2();
    realPhi->Sumw2();
    trackPhi->Sumw2();
    TOFPhi->Sumw2();
    TPCTrackPhi->Sumw2();
    TPCTOFPhi->Sumw2();
    PIDPhi->Sumw2();

    // pT vs mult
    TH1D* recoPhi_PT_mult[4];
    TH1D* realPhi_PT_mult[4];
    TH1D* trackPhi_PT_mult[4];
    TH1D* TOFPhi_PT_mult[4];
    TH1D* TPCTrackPhi_PT_mult[4];
    TH1D* TPCTOFPhi_PT_mult[4];
    TH1D* PIDPhi_PT_mult[4];

    TH1D* effPT_mult[4];
    TH1D* trackeffPT_mult[4];
    TH1D* TOF_TrackeffPT_mult[4];
    TH1D* TPC_TOFtrackeffPT_mult[4];
    TH1D* PID_TOFtrackeffPT_mult[4];
    TH1D* TPC_TrackeffPT_mult[4];

    TH1D* ratioPT_mult[3];
    TH1D* trackratioPT_mult[3];
    TH1D* TOF_TrackratioPT_mult[3];
    TH1D* TPC_TOFtrackratioPT_mult[3];
    TH1D* PID_TOFtrackratioPT_mult[3];
    TH1D* TPC_TrackratioPT_mult[3];
     
    TCanvas* cratioPT = new TCanvas("cratioPT", "cratioPT", 55, 55, 600, 600);
    TCanvas* ctrackratioPT = new TCanvas("ctrackratioPT", "ctrackratioPT", 55, 55, 600, 600);
    TCanvas* cTOF_TrackratioPT = new TCanvas("cTOF_TrackratioPT", "cTOF_TrackratioPT", 55, 55, 600, 600);
    TCanvas* cTPC_TOFtrackratioPT = new TCanvas("cTPC_TOFtrackratioPT", "cTPC_TOFtrackratioPT", 55, 55, 600, 600);
    TCanvas* cPID_TOFtrackratioPT = new TCanvas("cPID_TOFtrackratioPT", "cPID_TOFtrackratioPT", 55, 55, 600, 600);
    TCanvas* cTPC_TrackratioPT = new TCanvas("cTPC_TrackratioPT", "cTPC_TrackratioPT", 55, 55, 600, 600);

    TCanvas* ceffPT = new TCanvas("ceffPT", "ceffPT", 55, 55, 600, 600);
    TCanvas* ctrackeffPT = new TCanvas("ctrackeffPT", "ctrackeffPT", 55, 55, 600, 600);
    TCanvas* cTOF_TrackeffPT = new TCanvas("cTOF_TrackeffPT", "cTOF_TrackeffPT", 55, 55, 600, 600);
    TCanvas* cTPC_TOFtrackeffPT = new TCanvas("cTPC_TOFtrackeffPT", "cTPC_TOFtrackeffPT", 55, 55, 600, 600);
    TCanvas* cPID_TOFtrackeffPT = new TCanvas("cPID_TOFtrackeffPT", "cPID_TOFtrackeffPT", 55, 55, 600, 600);
    TCanvas* cTPC_TrackeffPT = new TCanvas("cTPC_TrackeffPT", "cTPC_TrackeffPT", 55, 55, 600, 600);

    plotMultEff(recoPhi, realPhi, recoPhi_PT_mult, realPhi_PT_mult, effPT_mult, ratioPT_mult, mult, 3, 0, ceffPT, cratioPT, kFALSE);
    plotMultEff(trackPhi, realPhi, trackPhi_PT_mult, realPhi_PT_mult, trackeffPT_mult, trackratioPT_mult, mult, 3, 0, ctrackeffPT, ctrackratioPT, kFALSE);
    plotMultEff(TOFPhi, trackPhi, TOFPhi_PT_mult, trackPhi_PT_mult, TOF_TrackeffPT_mult, TOF_TrackratioPT_mult, mult, 3, 0, cTOF_TrackeffPT, cTOF_TrackratioPT, kFALSE);
    plotMultEff(TPCTrackPhi, trackPhi, TPCTrackPhi_PT_mult, trackPhi_PT_mult, TPC_TrackeffPT_mult, TPC_TrackratioPT_mult, mult, 3, 0, cTPC_TrackeffPT, cTPC_TrackratioPT, kFALSE);
    plotMultEff(TPCTOFPhi, TOFPhi, TPCTOFPhi_PT_mult, TOFPhi_PT_mult, TPC_TOFtrackeffPT_mult, TPC_TOFtrackratioPT_mult, mult, 3, 0, cTPC_TOFtrackeffPT, cTPC_TOFtrackratioPT, kFALSE); 
    plotMultEff(PIDPhi, TPCTOFPhi, PIDPhi_PT_mult, TPCTOFPhi_PT_mult, PID_TOFtrackeffPT_mult, PID_TOFtrackratioPT_mult, mult, 3, 0, cPID_TOFtrackeffPT, cPID_TOFtrackratioPT, kFALSE);

    TH1D* totalPTeff0100 = (TH1D*)trackeffPT_mult[3]->Clone("totalPTeff0100");
    totalPTeff0100->Multiply(TOF_TrackeffPT_mult[3]);
    totalPTeff0100->Multiply(TPC_TOFtrackeffPT_mult[3]);
    totalPTeff0100->Multiply(TPC_TOFtrackeffPT_mult[3]); //estimate real TOF PID eff with TPC PID eff since it's flat at ~99.7%

    TCanvas* ctotal = new TCanvas("ctotal", "ctotal", 100, 50, 600, 600);
    ctotal->cd();
    totalPTeff0100->Draw();

    //for 0.5 bins
    printf("    pT    |    2.0-2.5    |    2.5-3.0    |    3.0-3.5    |    3.5-4.0    |\n");
    printf("total eff |    %6f    |    %6f    |    %6f    |    %6f    |\n", totalPTeff0100->GetBinContent(4), totalPTeff0100->GetBinContent(5), totalPTeff0100->GetBinContent(6), totalPTeff0100->GetBinContent(7));

    //for wide bins
    //printf("    pT    |    2.0-3.0    |    3.0-4.0    |\n");
    //printf("total eff |    %6f    |    %6f    |\n", totalPTeff0100->GetBinContent(2), totalPTeff0100->GetBinContent(3));
    

    TFile* output = new TFile("17f2befficiency.root", "RECREATE");
    effPT_mult[3]->SetName("phiPTEff");
    effPT_mult[3]->Write();
    trackeffPT_mult[3]->SetName("hadronPTEff");
    trackeffPT_mult[3]->Write();


    recoPhi->GetAxis(0)->SetRangeUser(1.0, 10.0);
    realPhi->GetAxis(0)->SetRangeUser(1.0, 10.0);
    //phi vs. Mult
    TH1D* recoPhi_Phi_mult[4];
    TH1D* trackPhi_Phi_mult[4];
    TH1D* realPhi_Phi_mult[4];
    TH1D* effPhi_mult[4];
    TH1D* trackeffPhi_mult[4];
    TH1D* TOFeffPhi_mult[4];
    TH1D* ratioPhi_mult[3];
    TH1D* trackratioPhi_mult[3];
    TH1D* TOFratioPhi_mult[3];
    TCanvas* ceffphi = new TCanvas("ceffphi", "ceffphimult", 55, 55, 600, 600);
    TCanvas* cratiophi = new TCanvas("cratioPhi", "cratioPhi", 55, 55, 600, 600);
    TCanvas* ctrackeffphi = new TCanvas("ctrackeffphi", "ctrackeffphi", 55, 55, 600, 600);
    TCanvas* ctrackratiophi = new TCanvas("ctrackratiophi", "ctrackratiophi", 60, 60, 600, 600);
    TCanvas* cTOFeffphi = new TCanvas("cTOFeffphi", "cTOFeffphi", 55, 55, 600, 600);
    TCanvas* cTOFratiophi = new TCanvas("cTOFratiophi", "cTOFratiophi", 60, 60, 600, 600);

    //recoPhi->GetAxis(2)->SetRangeUser(-0.3, 0.3);
    //realPhi->GetAxis(2)->SetRangeUser(-0.3, 0.3);
    //recoPhi->GetAxis(4)->SetRangeUser(-2., 2.);
    plotMultEff(recoPhi, realPhi, recoPhi_Phi_mult, realPhi_Phi_mult, effPhi_mult, ratioPhi_mult, mult, 3, 1, ceffphi, cratiophi);
    plotMultEff(trackPhi, realPhi, trackPhi_Phi_mult, realPhi_Phi_mult, trackeffPhi_mult, trackratioPhi_mult, mult, 3, 1, ctrackeffphi, ctrackratiophi);
    plotMultEff(recoPhi, trackPhi, recoPhi_Phi_mult, trackPhi_Phi_mult, TOFeffPhi_mult, TOFratioPhi_mult, mult, 3, 1, cTOFeffphi, cTOFratiophi);
    
    //recoPhi->GetAxis(2)->SetRangeUser(-0.8, 0.8);
    //realPhi->GetAxis(2)->SetRangeUser(-0.8, 0.8);
    //recoPhi->GetAxis(4)->SetRangeUser(-10., 10.);
    
    //eta vs. Mult
    TH1D* recoPhi_Eta_mult[4];
    TH1D* trackPhi_Eta_mult[4];
    TH1D* TOFPhi_Eta_mult[4];
    TH1D* realPhi_Eta_mult[4];
    TH1D* effEta_mult[4];
    TH1D* trackeffEta_mult[4];
    TH1D* TOFeffEta_mult[4];
    TH1D* ratioEta_mult[3];
    TH1D* trackratioEta_mult[3];
    TH1D* TOFratioEta_mult[3];
    TCanvas* ceffeta = new TCanvas("ceffeta", "ceffetamult", 55, 55, 600, 600);
    TCanvas* ctrackeffeta = new TCanvas("ctrackeffeta", "ctrackeffetamult", 55, 55, 600, 600);
    TCanvas* cTOFeffeta = new TCanvas("cTOFeffeta", "cTOFeffetamult", 55, 55, 600, 600);
    TCanvas* cratioEta = new TCanvas("cratioEta", "cratioEta", 55, 55, 600, 600);
    TCanvas* ctrackratioEta = new TCanvas("ctrackratioEta", "ctrackratioEta", 55, 55, 600, 600);
    TCanvas* cTOFratioEta = new TCanvas("cTOFratioEta", "cTOFratioEta", 55, 55, 600, 600);

    plotMultEff(recoPhi, realPhi, recoPhi_Eta_mult, realPhi_Eta_mult, effEta_mult, ratioEta_mult, mult, 3, 2, ceffeta, cratioEta);
    plotMultEff(trackPhi, realPhi, trackPhi_Eta_mult, realPhi_Eta_mult, trackeffEta_mult, trackratioEta_mult, mult, 3, 2, ctrackeffeta, ctrackratioEta);
    plotMultEff(recoPhi, trackPhi, TOFPhi_Eta_mult, trackPhi_Eta_mult, TOFeffEta_mult, TOFratioEta_mult, mult, 3, 2, cTOFeffeta, cTOFratioEta);


    //Z vs. Mult
    TH1D* recoPhi_Z_mult[4];
    TH1D* trackPhi_Z_mult[4];
    TH1D* TOFPhi_Z_mult[4];
    TH1D* realPhi_Z_mult[4];
    TH1D* effZ_mult[4];
    TH1D* trackeffZ_mult[4];
    TH1D* TOFeffZ_mult[4];
    TH1D* ratioZ_mult[4];
    TH1D* trackratioZ_mult[4];
    TH1D* TOFratioZ_mult[4];
    TCanvas* ceffZ = new TCanvas("ceffZ", "ceffZmult", 55, 55, 600, 600);
    TCanvas* ctrackeffZ = new TCanvas("ctrackeffZ", "ctrackeffZmult", 55, 55, 600, 600);
    TCanvas* cTOFeffZ = new TCanvas("cTOFeffZ", "cTOFeffZmult", 55, 55, 600, 600);
    TCanvas* cratioZ = new TCanvas("cratioZ", "cratioZ", 55, 55, 600, 600);
    TCanvas* ctrackratioZ = new TCanvas("ctrackratioZ", "ctrackratioZ", 55, 55, 600, 600);
    TCanvas* cTOFratioZ = new TCanvas("cTOFratioZ", "cTOFratioZ", 55, 55, 600, 600);
 
    plotMultEff(recoPhi, realPhi, recoPhi_Z_mult, realPhi_Z_mult, effZ_mult, ratioZ_mult, mult, 3, 4, ceffZ, cratioZ);
    plotMultEff(trackPhi, realPhi, trackPhi_Z_mult, realPhi_Z_mult, trackeffZ_mult, trackratioZ_mult, mult, 3, 4, ctrackeffZ, ctrackratioZ);
    plotMultEff(recoPhi, trackPhi, TOFPhi_Z_mult, trackPhi_Z_mult, TOFeffZ_mult, TOFratioZ_mult, mult, 3, 4, cTOFeffZ, cTOFratioZ);
 
    // 0-100% mult with Z vertex dependency 
    realPhi->GetAxis(6)->SetRangeUser(0.0, 100.0);
    recoPhi->GetAxis(6)->SetRangeUser(0.0, 100.0);

    TH1D* recoPhi_PT = recoPhi->Projection(0);
    recoPhi_PT->SetName("recoPT");
    recoPhi_PT->Rebin(5);
    TH1D* realPhi_PT = realPhi->Projection(0);
    realPhi_PT->SetName("realPT");
    realPhi_PT->Rebin(5);

    TH1D* effPT = (TH1D*)recoPhi_PT->Clone("Efficiency vs. p_{T};p_{T}");
    effPT->Divide(realPhi_PT);
    effPT->SetTitle("Efficiency vs. p_{T};p_{T} (GeV/c)");
   
    TCanvas* c1 = new TCanvas("c1", "recon pT", 50, 50, 600, 600);
    c1->cd();
    recoPhi_PT->Draw();

    TCanvas* c2 = new TCanvas("c2", "real pT", 50, 50, 600, 600);
    c2->cd();
    realPhi_PT->Draw();

    TCanvas* c3 = new TCanvas ("c3", "eff PT", 50, 50, 600, 600);
    c3->cd();
    effPT->Draw();

    //recoPhi->GetAxis(0)->SetRangeUser(2.0, 4.0);
    //realPhi->GetAxis(0)->SetRangeUser(2.0, 4.0);

    TH1D* recoPhi_Phi = recoPhi->Projection(1);
    recoPhi_Phi->SetName("recoPhiPhi");

    TH1D* realPhi_Phi = realPhi->Projection(1);
    realPhi_Phi->SetName("realPhiPhi");


    TH1D* effPhi = (TH1D*)recoPhi_Phi->Clone("Efficiency vs. #varphi for 2.0 < p^{#phi}_{T} < 4.0; #phi");
    effPhi->Divide(realPhi_Phi);
    effPhi->SetTitle("Efficiency vs. #varphi for 2.0 < p^{#phi}_{T} < 4.0 GeV/c; #phi");

    TCanvas* c4 = new TCanvas("c4" , "recon phi", 50, 50, 600, 600);
    c4->cd();
    recoPhi_Phi->Draw();

    TCanvas* c5 = new TCanvas("c5", "real phi", 50, 50, 600, 600);
    c5->cd();
    realPhi_Phi->Draw();

    TCanvas* c6 = new TCanvas("c6", "eff phi", 55, 55, 600, 600);
    c6->cd();
    effPhi->Draw();

    TH1D* recoPhi_Eta = recoPhi->Projection(2);
    recoPhi_Eta->SetName("recoPhiEta");
    TH1D* realPhi_Eta = realPhi->Projection(2);
    realPhi_Eta->SetName("realPhiEta");

    TH1D* effEta = (TH1D*)recoPhi_Eta->Clone("Efficiency vs. #eta for 2.0 < p^{#phi}_{T} < 4.0 GeV/c; #eta");
    effEta->Divide(realPhi_Eta);
    effEta->SetTitle("Efficiency vs. #eta for 2.0 < p^{#phi}_{T} < 4.0; #eta");

    TCanvas* c7 = new TCanvas("c7" , "recon eta", 50, 50, 600, 600);
    c7->cd();
    recoPhi_Eta->Draw();

    TCanvas* c8 = new TCanvas("c8", "real eta", 50, 50, 600, 600);
    c8->cd();
    realPhi_Eta->Draw();

    TCanvas* c9 = new TCanvas("c9", "eff eta", 55, 55, 600, 600);
    c9->cd();
    effEta->Draw();

    
    TH1D* recoPhi_Z = recoPhi->Projection(4);
    recoPhi_Z->SetName("recoPhiZ");
    TH1D* realPhi_Z = realPhi->Projection(4);
    realPhi_Z->SetName("realPhiZ");

    TH1D* effZ = (TH1D*)recoPhi_Z->Clone("Efficiency vs. ZVertex for 2.0 < p^{#phi}_{T} < 4.0 GeV/c; Z Vertex (cm)");
    effZ->Divide(realPhi_Z);
    effZ->SetTitle("Efficiency vs. Zvertex for 2.0 < p^{#phi}_{T} < 4.0; Z Vertex (cm)");

    TCanvas* c10 = new TCanvas("c10" , "recon Z", 50, 50, 600, 600);
    c10->cd();
    recoPhi_Z->Draw();

    TCanvas* c11 = new TCanvas("c11", "real Z", 50, 50, 600, 600);
    c11->cd();
    realPhi_Z->Draw();

    TCanvas* c12 = new TCanvas("c12", "eff Z", 55, 55, 600, 600);
    c12->cd();
    effZ->Draw();

    


    recoPhi->GetAxis(0)->SetRange(0,0);
    realPhi->GetAxis(0)->SetRange(0,0);
    TH1D* recoPhi_PT_Z[10];
    TH1D* realPhi_PT_Z[10];
    TH1D* effPT_Z[10];

    //Int_t colors[10] = {kRed+1, kOrange+2, kOrange-2, kSpring+9, kGreen+2, kTeal+2, kCyan+2, kBlue+1, kViolet+7, kViolet-3};
    Int_t colors[10] = {kRed-10, kRed-9, kRed-7, kRed-3, kRed-4, kRed, kRed+1, kRed+2, kRed+3, kRed+4};


    TCanvas* cptz = new TCanvas("cptz", "eff pt z", 60, 60, 600, 600);
    cptz->cd();
    for(int iz = 0; iz < 10; iz++){
        recoPhi->GetAxis(4)->SetRange(iz+1, iz+1);
        realPhi->GetAxis(4)->SetRange(iz+1, iz+1);

        recoPhi_PT_Z[iz] = recoPhi->Projection(0);
        recoPhi_PT_Z[iz]->SetName(Form("reco_PT_Z_%i", iz));
        recoPhi_PT_Z[iz]->Rebin(5);

        realPhi_PT_Z[iz] = realPhi->Projection(0);
        realPhi_PT_Z[iz]->SetName(Form("real_PT_Z_%i", iz));
        realPhi_PT_Z[iz]->Rebin(5);

        effPT_Z[iz] = (TH1D*)recoPhi_PT_Z[iz]->Clone(Form("effPT_Z_%i", iz));
        effPT_Z[iz]->Divide(realPhi_PT_Z[iz]);
        effPT_Z[iz]->SetTitle("Efficiency vs. p_{T} for ZVertex bins:p_{T} (GeV/c)");
        effPT_Z[iz]->SetLineColor(colors[iz]);
        effPT_Z[iz]->SetMarkerColor(colors[iz]);
        effPT_Z[iz]->SetMarkerSize(2);
        effPT_Z[iz]->SetMarkerStyle(22);
        //effPT_Z[iz]->SetLineWidth(2);
        effPT_Z[iz]->Draw("P SAME");
    }

    //recoPhi->GetAxis(0)->SetRangeUser(2.0,4.0);
    //realPhi->GetAxis(0)->SetRange(2.0,4.0);
    TH1D* recoPhi_Phi_Z[10];
    TH1D* realPhi_Phi_Z[10];
    TH1D* effPhi_Z[10];

    TCanvas* cphiz = new TCanvas("cphiz", "eff phi z", 60, 60, 600, 600);
    cphiz->cd();
    for(int iz = 0; iz < 10; iz++){
        recoPhi->GetAxis(4)->SetRange(iz+1, iz+1);
        realPhi->GetAxis(4)->SetRange(iz+1, iz+1);

        recoPhi_Phi_Z[iz] = recoPhi->Projection(1);
        recoPhi_Phi_Z[iz]->SetName(Form("reco_Phi_Z_%i", iz));
        recoPhi_Phi_Z[iz]->Rebin(4);

        realPhi_Phi_Z[iz] = realPhi->Projection(1);
        realPhi_Phi_Z[iz]->SetName(Form("real_Phi_Z_%i", iz));
        realPhi_Phi_Z[iz]->Rebin(4);

        effPhi_Z[iz] = (TH1D*)recoPhi_Phi_Z[iz]->Clone(Form("effPhi_Z_%i", iz));
        effPhi_Z[iz]->Divide(realPhi_Phi_Z[iz]);
        effPhi_Z[iz]->SetTitle("Efficiency vs. #varphi for ZVertex bins;#varphi");
        effPhi_Z[iz]->SetLineColor(colors[iz]);
        effPhi_Z[iz]->SetMarkerColor(colors[iz]);
        effPhi_Z[iz]->SetMarkerStyle(23);
        effPhi_Z[iz]->SetMarkerSize(2);
        //effPhi_Z[iz]->SetLineWidth(2);
        effPhi_Z[iz]->Draw("P SAME");
    }

    recoPhi->GetAxis(2)->SetRange(0,0);
    realPhi->GetAxis(2)->SetRange(0,0);
    TH1D* recoPhi_Eta_Z[10];
    TH1D* realPhi_Eta_Z[10];
    TH1D* effEta_Z[10];

    TCanvas* cetaz = new TCanvas("cetaz", "eff eta z", 60, 60, 600, 600);
    cetaz->cd();
    for(int iz = 0; iz < 10; iz++){
        recoPhi->GetAxis(4)->SetRange(iz+1, iz+1);
        realPhi->GetAxis(4)->SetRange(iz+1, iz+1);

        recoPhi_Eta_Z[iz] = recoPhi->Projection(2);
        recoPhi_Eta_Z[iz]->SetName(Form("reco_Eta_Z_%i", iz));
        recoPhi_Eta_Z[iz]->Rebin(4);

        realPhi_Eta_Z[iz] = realPhi->Projection(2);
        realPhi_Eta_Z[iz]->SetName(Form("real_Eta_Z_%i", iz));
        realPhi_Eta_Z[iz]->Rebin(4);

        effEta_Z[iz] = (TH1D*)recoPhi_Eta_Z[iz]->Clone(Form("effEta_Z_%i", iz));
        effEta_Z[iz]->Divide(realPhi_Eta_Z[iz]);
        effEta_Z[iz]->SetTitle("Efficiency vs. #eta for ZVertex bins;#eta");
        effEta_Z[iz]->SetLineColor(colors[iz]);
        effEta_Z[iz]->SetMarkerColor(colors[iz]);
        effEta_Z[iz]->SetMarkerStyle(23);
        effEta_Z[iz]->SetMarkerSize(2);
        //effEta_Z[iz]->SetLineWidth(2);
        effEta_Z[iz]->Draw("P SAME");
    }

    //2D eta phi efficiency map for 0-100%
    recoPhi->GetAxis(0)->SetRangeUser(0.5, 10.0);
    recoPhi->GetAxis(4)->SetRangeUser(-10.0, 10.0);
    recoPhi->GetAxis(5)->SetRangeUser(0.0, 100.0);
    realPhi->GetAxis(0)->SetRangeUser(0.5, 10.0);
    realPhi->GetAxis(4)->SetRangeUser(-10.0, 10.0);
    realPhi->GetAxis(5)->SetRangeUser(0.0, 100.0);

    TH2D* recoEtaPhi = (TH2D*)recoPhi->Projection(2, 1);
    recoEtaPhi->SetTitle("Reconstructed Phi vs. #eta and #varphi; #varphi; #eta");
    recoEtaPhi->SetName("reco_eta_phi");

    TH2D* realEtaPhi = (TH2D*)realPhi->Projection(2, 1);
    realEtaPhi->SetTitle("Real Phi vs. #eta and #varphi; #varphi; #eta");
    realEtaPhi->SetName("real_eta_phi");

    TH2D* effEtaPhi = (TH2D*)recoEtaPhi->Clone("effEtaPhi");
    effEtaPhi->Divide(realEtaPhi);

    TCanvas* cetaphi = new TCanvas("cetaphi", "cetaphi", 50, 50, 600, 600);
    cetaphi->cd();
    effEtaPhi->Draw("colz");


    TCanvas* clegend = new TCanvas("clegend", "clegend", 75, 75, 300, 100);
    clegend->cd();
    TLegend *legend = new TLegend(0.0, 0.0, 1.0, 1.0);
    legend->AddEntry(effPhi_mult[3], "0-100% Efficiency", "ep");
    legend->AddEntry(effPhi_mult[0], "0-20% Efficiency", "ep");
    legend->AddEntry(effPhi_mult[1], "20-50% Efficiency", "ep");
    legend->AddEntry(effPhi_mult[2], "50-80% Efficiency", "ep");
    legend->Draw();

}
