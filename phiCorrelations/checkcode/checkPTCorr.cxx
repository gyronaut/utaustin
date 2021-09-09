void checkPTCorr(){

    TFile *histofile = new TFile("~/phistudies/results_onlineEff/Combined/TOTAL_hphi_0_20_50_80.root");
    TList* list = (TList*)histofile->Get("phiCorr_mult_0_20");

    TH1D* vtxZmixbins = (TH1D*)list->FindObject("fVtxZmixbins");
    Int_t numbinsZvtx = vtxZmixbins->GetXaxis()->GetNbins();

    THnSparseF *dphiHPhi[numbinsZvtx];
    THnSparseF *dphiHKK[numbinsZvtx];
    TH2D *trigphiPTcorr;

    for(int izvtx = 0; izvtx < numbinsZvtx; izvtx++){
        dphiHPhi[izvtx] = (THnSparseF *)list->FindObject(Form("fDphiHPhiz%i", izvtx));
        dphiHKK[izvtx] = (THnSparseF *)list->FindObject(Form("fDphiHKKz%i", izvtx));

        dphiHPhi[izvtx]->GetAxis(4)->SetRangeUser(1.014, 1.026);
        if(izvtx == 0){
            trigphiPTcorr = (TH2D*)dphiHPhi[izvtx]->Projection(0,1);
        }else{
            trigphiPTcorr->Add((TH2D*)dphiHPhi[izvtx]->Projection(0,1));
        }
    }

    TH1D *trigPTHighPhiPT, *trigPTMidPhiPT, *trigPTLowPhiPT;
    trigPTHighPhiPT = (TH1D*)trigphiPTcorr->ProjectionY("trigPTHighPhiPT", trigphiPTcorr->GetXaxis()->FindBin(4.001), trigphiPTcorr->GetXaxis()->FindBin(4.999));
    trigPTMidPhiPT = (TH1D*)trigphiPTcorr->ProjectionY("trigPTMidPhiPT", trigphiPTcorr->GetXaxis()->FindBin(3.001), trigphiPTcorr->GetXaxis()->FindBin(3.999));
    trigPTLowPhiPT = (TH1D*)trigphiPTcorr->ProjectionY("trigPTLowPhiPT", trigphiPTcorr->GetXaxis()->FindBin(2.001), trigphiPTcorr->GetXaxis()->FindBin(2.999));

    TCanvas* c2D = new TCanvas("c2D", "c2D", 50, 50, 600, 600);
    c2D->cd();
    trigphiPTcorr->Draw("COLZ");

    TCanvas* clow = new TCanvas("clow", "clow", 50, 50, 600, 600);
    clow->cd();
    trigPTLowPhiPT->Draw();

    TCanvas* cmid = new TCanvas("cmid", "cmid", 50, 50, 600, 600);
    cmid->cd();
    trigPTMidPhiPT->Draw();

    TCanvas* chigh = new TCanvas("chigh", "chigh", 50, 50, 600, 600);
    chigh->cd();
    trigPTHighPhiPT->Draw();

    TH1D *highptnorm, *midptnorm, *lowptnorm;
    highptnorm = (TH1D*)trigPTHighPhiPT->Clone("highptnorm");
    highptnorm->Scale(1.0/highptnorm->Integral());
    highptnorm->SetLineColor(kBlue+1);
    highptnorm->SetLineWidth(2);
    highptnorm->SetMarkerColor(kBlue+1);
    highptnorm->SetMarkerSize(2);
    highptnorm->SetMarkerStyle(22);
    midptnorm = (TH1D*)trigPTMidPhiPT->Clone("midptnorm");
    midptnorm->Scale(1.0/midptnorm->Integral());
    midptnorm->SetLineColor(kRed+1);
    midptnorm->SetLineWidth(2);
    midptnorm->SetMarkerColor(kRed+1);
    midptnorm->SetMarkerSize(2);
    midptnorm->SetMarkerStyle(24);
    lowptnorm = (TH1D*)trigPTLowPhiPT->Clone("lowptnorm");
    lowptnorm->Scale(1.0/lowptnorm->Integral());
    lowptnorm->SetLineColor(kGreen+1);
    lowptnorm->SetLineWidth(2);
    lowptnorm->SetMarkerColor(kGreen+1);
    lowptnorm->SetMarkerSize(2);
    lowptnorm->SetMarkerStyle(23);

    TCanvas* cptcomp = new TCanvas("cptcomp", "cptcomp", 50, 50, 600, 600);
    cptcomp->cd();
//    highptnorm->GetXaxis()->SetRangeUser(0.0, 10.0);
    lowptnorm->Draw("HIST PL");
    highptnorm->Draw("SAME HIST PL");
    midptnorm->Draw("SAME HIST PL");




}

