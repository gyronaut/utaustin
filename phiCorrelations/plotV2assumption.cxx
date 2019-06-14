void plotV2assumption(){
    TFile* filev2 = new TFile("~/phiStudies/results_onlineEff/Combined/fitsyst_v2.root");
    TGraphErrors* nearv2 = (TGraphErrors*)filev2->Get("ratiosNear")->Clone("nearv2");
    TGraphErrors* awayv2 = (TGraphErrors*)filev2->Get("ratiosAway")->Clone("awayv2");
    TGraphErrors* uev2 = (TGraphErrors*)filev2->Get("ratiosBulk")->Clone("uev2");

    TFile* fileno = new TFile("~/phiStudies/results_onlineEff/Combined/fitsyst_nov2.root");
    TGraphErrors* nearno = (TGraphErrors*)fileno->Get("ratiosNear")->Clone("nearno");
    TGraphErrors* awayno = (TGraphErrors*)fileno->Get("ratiosAway")->Clone("awayno");
    TGraphErrors* ueno = (TGraphErrors*)fileno->Get("ratiosBulk")->Clone("ueno");

    TGraphErrors* nearavg = (TGraphErrors*)nearv2->Clone("nearavg");
    nearavg->SetMarkerSize(0);
    nearavg->SetFillColor(kGray+1);
    nearavg->SetFillStyle(3001);
    TGraphErrors* awayavg = (TGraphErrors*)awayv2->Clone("awayavg");
    awayavg->SetMarkerSize(0);
    awayavg->SetFillColor(kGray+1);
    awayavg->SetFillStyle(3001);
    TGraphErrors* ueavg = (TGraphErrors*)uev2->Clone("ueavg");
    ueavg->SetMarkerSize(0);
    ueavg->SetFillColor(kGray+1);
    ueavg->SetFillStyle(3001);
  
    Double_t xpoint = 0.0, ypoint = 0.0;
    Double_t xpointno = 0.0, ypointno = 0.0;
    
    for(int i =0; i<3; i++){
        nearv2->GetPoint(i, xpoint, ypoint);
        nearno->GetPoint(i, xpointno, ypointno);
        nearavg->SetPoint(i, xpoint, 0.5*(ypointno-ypoint) + ypoint);
        nearavg->SetPointError(i, 2.5, TMath::Abs(0.5*(ypointno-ypoint)));
        awayv2->GetPoint(i, xpoint, ypoint);
        awayno->GetPoint(i, xpointno, ypointno);
        awayavg->SetPoint(i, xpoint, 0.5*(ypointno-ypoint) + ypoint);
        awayavg->SetPointError(i, 2.5, TMath::Abs(0.5*(ypointno-ypoint)));
        uev2->GetPoint(i, xpoint, ypoint);
        ueno->GetPoint(i, xpointno, ypointno);
        ueavg->SetPoint(i, xpoint, 0.5*(ypointno-ypoint) + ypoint);
        ueavg->SetPointError(i, 2.5, TMath::Abs(0.5*(ypointno-ypoint)));
    }   

    TCanvas* c1 = new TCanvas("c1", "c1", 50, 50, 600, 600);
    c1->cd();
    nearv2->Draw("ALP");
    nearno->Draw("LP");
    nearavg->Draw("2");

    TCanvas* c2 = new TCanvas("c2", "c2", 50, 50, 600, 600);
    c2->cd();
    awayv2->Draw("ALP");
    awayno->Draw("LP");
    awayavg->Draw("2");

    TCanvas* c3 = new TCanvas("c3", "c3", 50, 50, 600, 600);
    c3->cd();
    uev2->Draw("ALP");
    ueno->Draw("LP");
    ueavg->Draw("2");

    TFile* outfile = new TFile("v2syst.root", "RECREATE");
    nearavg->SetName("nearv2syst");
    nearavg->Write();
    awayavg->SetName("awayv2syst");
    awayavg->Write();
}
