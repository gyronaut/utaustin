void plotJetRatio(){

    TFile* lowptfile = new TFile("/Users/jtblair/alidock/alirepos/utaustin/phiCorrelations/plotcode/ratiocode/fitsyst_testlowpt.root");

    TGraphErrors* lowptjethh = (TGraphErrors*)lowptfile->Get("jetratioshh");
    lowptjethh->SetName("lowptjethh");
    TGraphErrors* lowptjethphi = (TGraphErrors*)lowptfile->Get("jetratioshPhi");
    lowptjethphi->SetName("lowptjethphi");

    TGraphErrors* lowptjethhnorm = (TGraphErrors*)lowptjethh->Clone("lowptjethhnorm");
    lowptjethhnorm->SetPoint(0, lowptjethhnorm->GetPointX(0), lowptjethhnorm->GetPointY(0)/lowptjethh->GetPointY(0));
    lowptjethhnorm->SetPoint(1, lowptjethhnorm->GetPointX(1), lowptjethhnorm->GetPointY(1)/lowptjethh->GetPointY(0));
    lowptjethhnorm->SetPoint(2, lowptjethhnorm->GetPointX(2), lowptjethhnorm->GetPointY(2)/lowptjethh->GetPointY(0));
    lowptjethhnorm->GetYaxis()->SetTitle("Yield Ratio (#frac{Jet Yield}{Total})");
    lowptjethhnorm->SetTitle("Low p_{T} Jet Ratio");

    TGraphErrors* lowptjethphinorm = (TGraphErrors*)lowptjethphi->Clone("lowptjethphinorm");
    lowptjethphinorm->SetPoint(0, lowptjethphinorm->GetPointX(0), lowptjethphinorm->GetPointY(0)/lowptjethphi->GetPointY(0));
    lowptjethphinorm->SetPoint(1, lowptjethphinorm->GetPointX(1), lowptjethphinorm->GetPointY(1)/lowptjethphi->GetPointY(0));
    lowptjethphinorm->SetPoint(2, lowptjethphinorm->GetPointX(2), lowptjethphinorm->GetPointY(2)/lowptjethphi->GetPointY(0));


    
    TFile* highptfile = new TFile("/Users/jtblair/alidock/alirepos/utaustin/phiCorrelations/plotcode/ratiocode/fitsyst_test.root");

    TGraphErrors* highptjethh = (TGraphErrors*)highptfile->Get("jetratioshh");
    highptjethh->SetName("highptjethh");
    TGraphErrors* highptjethphi = (TGraphErrors*)highptfile->Get("jetratioshPhi");
    highptjethphi->SetName("highptjethphi");

    TGraphErrors* highptjethhnorm = (TGraphErrors*)highptjethh->Clone("highptjethhnorm");
    highptjethhnorm->SetPoint(0, highptjethhnorm->GetPointX(0), highptjethhnorm->GetPointY(0)/highptjethh->GetPointY(0));
    highptjethhnorm->SetPoint(1, highptjethhnorm->GetPointX(1), highptjethhnorm->GetPointY(1)/highptjethh->GetPointY(0));
    highptjethhnorm->SetPoint(2, highptjethhnorm->GetPointX(2), highptjethhnorm->GetPointY(2)/highptjethh->GetPointY(0));
    highptjethhnorm->GetYaxis()->SetTitle("Yield Ratio (#frac{JetYield}{Total})");
    highptjethhnorm->SetTitle("High p_{T} Jet Ratio");

    TGraphErrors* highptjethphinorm = (TGraphErrors*)highptjethphi->Clone("highptjethphinorm");
    highptjethphinorm->SetPoint(0, highptjethphinorm->GetPointX(0), highptjethphinorm->GetPointY(0)/highptjethphi->GetPointY(0));
    highptjethphinorm->SetPoint(1, highptjethphinorm->GetPointX(1), highptjethphinorm->GetPointY(1)/highptjethphi->GetPointY(0));
    highptjethphinorm->SetPoint(2, highptjethphinorm->GetPointX(2), highptjethphinorm->GetPointY(2)/highptjethphi->GetPointY(0));


    TCanvas *clow = new TCanvas("clow", "clow", 50, 50, 600, 600);
    clow->cd();
    lowptjethhnorm->Draw();
    lowptjethphinorm->Draw("SAME PL");

    TCanvas *chigh = new TCanvas("chigh", "chigh", 50, 50, 600, 600);
    chigh->cd();
    highptjethhnorm->Draw();
    highptjethphinorm->Draw("SAME PL");

}
