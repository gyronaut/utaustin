void plotUSCorrected(string inputname){
    gStyle->SetOptStat(0);
    gStyle->SetOptFit(0);

    TFile* eta20File = new TFile(inputname.c_str());

    TH2D* eta20peak = (TH2D*)eta20File->Get("AvgUSsubhPhi2Dpeakavgscale");
    TH2D* eta20RSB = (TH2D*)eta20File->Get("RLSsubhPhi2DRside");
    TH2D* eta20LSB = (TH2D*)eta20File->Get("LLSsubhPhi2DLside");
    
    TH1D* scales = (TH1D*)eta20File->Get("scales");
    TH2D* uncorr2D = (TH2D*)eta20File->Get("uncorrectedhPhi2Dpeak");
    TH2D* BG2D = (TH2D*)eta20File->Get("hPhiBGPeakregion");
    BG2D->Scale(scales->GetBinContent(6)); 

    TH2D* residual = (TH2D*)eta20peak->Clone("residual");
    TH2D* peak = (TH2D*)eta20peak->Clone("peak");
    TH2D* hfits[4];
    TH2D* hotherfits[4];
    for(int i =0; i<4; i++){
        hfits[i] = (TH2D*)eta20peak->Clone(Form("hfit%d", i));
        hotherfits[i] = (TH2D*)eta20peak->Clone(Form("hotherfit%d", i));

    }

    Float_t epsilon = 0.001;

    eta20peak->GetXaxis()->SetTitle("#Delta#eta");
    eta20peak->GetXaxis()->SetTitleSize(0.05);
    eta20peak->GetXaxis()->SetTitleOffset(1.3);
    eta20peak->GetYaxis()->SetTitle("#Delta#varphi");
    eta20peak->GetYaxis()->SetTitleSize(0.05);
    eta20peak->GetYaxis()->SetTitleOffset(1.3);
    eta20peak->SetTitle("");
    //eta20peak->SetStats(kFALSE);
    //eta20peak->Scale(1.0/(eta20peak->Integral(eta20peak->GetXaxis()->FindBin(-1.2 + epsilon), eta20peak->GetXaxis()->FindBin(1.2 - epsilon), 1, eta20peak->GetYaxis()->GetNbins())));

    eta20RSB->GetXaxis()->SetTitle("#Delta#eta");
    eta20RSB->GetXaxis()->SetTitleSize(0.05);
    eta20RSB->GetXaxis()->SetTitleOffset(1.3);
    eta20RSB->GetYaxis()->SetTitle("#Delta#varphi");
    eta20RSB->GetYaxis()->SetTitleSize(0.05);
    eta20RSB->GetYaxis()->SetTitleOffset(1.3);
    eta20RSB->SetTitle("");
    eta20RSB->SetStats(kFALSE);

    eta20LSB->GetXaxis()->SetTitle("#Delta#eta");
    eta20LSB->GetXaxis()->SetTitleSize(0.05);
    eta20LSB->GetXaxis()->SetTitleOffset(1.3);
    eta20LSB->GetYaxis()->SetTitle("#Delta#varphi");
    eta20LSB->GetYaxis()->SetTitleSize(0.05);
    eta20LSB->GetYaxis()->SetTitleOffset(1.3);
    eta20LSB->SetTitle("");
    eta20LSB->SetStats(kFALSE);

    TH1D* eta20peakEta = eta20peak->ProjectionX("eta20peakEta", 1, eta20peak->GetYaxis()->GetNbins());
    eta20peakEta->GetXaxis()->SetTitleOffset(1.0);
    eta20peakEta->GetXaxis()->SetRangeUser(-1.2, 1.2);
    eta20peakEta->SetStats(kFALSE);
    TH1D* eta20peakPhi = eta20peak->ProjectionY("eta20peakPhi", eta20peak->GetXaxis()->FindBin(-1.5 + epsilon), eta20peak->GetXaxis()->FindBin(1.5 - epsilon));
    eta20peakPhi->GetXaxis()->SetTitleOffset(1.0);
    eta20peakPhi->SetStats(kFALSE);
    eta20peakPhi->SetLineColor(kBlue);
    TH1D* eta20peakPhiNarrow = eta20peak->ProjectionY("eta20peakPhiNarrow", eta20peak->GetXaxis()->FindBin(-1.2 + epsilon), eta20peak->GetXaxis()->FindBin(1.2 - epsilon));
    eta20peakPhiNarrow->SetLineColor(kViolet);
    eta20peakPhiNarrow->SetStats(kFALSE);
    eta20peakPhiNarrow->GetXaxis()->SetTitleOffset(1.0);
    TH1D* eta20peakPhiNarrowest = eta20peak->ProjectionY("eta20peakPhiNarrowest", eta20peak->GetXaxis()->FindBin(-1.0 + epsilon), eta20peak->GetXaxis()->FindBin(1.0 - epsilon));
    eta20peakPhiNarrowest->SetLineColor(kRed);

    TH1D* eta20RSBEta = eta20RSB->ProjectionX("eta20RSBEta", 1, eta20RSB->GetYaxis()->GetNbins());
    eta20RSBEta->GetXaxis()->SetTitleOffset(1.0);
    eta20RSBEta->GetXaxis()->SetRangeUser(-1.2, 1.2);
    eta20RSBEta->SetStats(kFALSE);
    TH1D* eta20RSBPhi = eta20RSB->ProjectionY("eta20RSBPhi", eta20RSB->GetXaxis()->FindBin(-1.5 + epsilon), eta20RSB->GetXaxis()->FindBin(1.5 - epsilon));
    eta20RSBPhi->GetXaxis()->SetTitleOffset(1.0);
    eta20RSBPhi->SetLineColor(kBlue);
    eta20RSBPhi->SetStats(kFALSE);
    TH1D* eta20RSBPhiNarrow = eta20RSB->ProjectionY("eta20RSBPhiNarrow", eta20RSB->GetXaxis()->FindBin(-1.2 + epsilon), eta20RSB->GetXaxis()->FindBin(1.2 - epsilon));
    eta20RSBPhiNarrow->SetLineColor(kViolet);
    TH1D* eta20RSBPhiNarrowest = eta20RSB->ProjectionY("eta20RSBPhiNarrowest", eta20RSB->GetXaxis()->FindBin(-1.0 + epsilon), eta20RSB->GetXaxis()->FindBin(1.0 - epsilon));
    eta20RSBPhiNarrowest->SetLineColor(kRed);

    TH1D* eta20LSBEta = eta20LSB->ProjectionX("eta20LSBEta", 1, eta20LSB->GetYaxis()->GetNbins());
    eta20LSBEta->GetXaxis()->SetTitleOffset(1.0);
    eta20LSBEta->GetXaxis()->SetRangeUser(-1.2, 1.2);
    eta20LSBEta->SetStats(kFALSE);
    TH1D* eta20LSBPhi = eta20LSB->ProjectionY("eta20LSBPhi", eta20LSB->GetXaxis()->FindBin(-1.5 + epsilon), eta20LSB->GetXaxis()->FindBin(1.5 - epsilon));
    eta20LSBPhi->GetXaxis()->SetTitleOffset(1.0);
    eta20LSBPhi->SetLineColor(kBlue);
    eta20LSBPhi->SetStats(kFALSE);
    TH1D* eta20LSBPhiNarrow = eta20LSB->ProjectionY("eta20LSBPhiNarrow", eta20LSB->GetXaxis()->FindBin(-1.2 + epsilon), eta20LSB->GetXaxis()->FindBin(1.2 - epsilon));
    eta20LSBPhiNarrow->SetLineColor(kViolet);
    TH1D* eta20LSBPhiNarrowest = eta20LSB->ProjectionY("eta20LSBPhiNarrowest", eta20LSB->GetXaxis()->FindBin(-1.0 + epsilon), eta20LSB->GetXaxis()->FindBin(1.0 - epsilon));
    eta20LSBPhiNarrowest->SetLineColor(kRed);


    //reset eta range to narrow view for 2D plotting and rebin
//    eta20peak->Rebin2D(2,2);
//    eta20RSB->Rebin2D(2,2);
//    eta20LSB->Rebin2D(2,2);
    eta20peak->GetXaxis()->SetRangeUser(-1.2, 1.2);
    eta20RSB->GetXaxis()->SetRangeUser(-1.2, 1.2);
    eta20LSB->GetXaxis()->SetRangeUser(-1.2, 1.2);

    //set up legend for delta-phi projection
    TLegend* legend = new TLegend(0.5612, 0.6716, 0.89, 0.8712);
    legend->SetMargin(0.15);
//    legend->AddEntry(eta20peakPhi, "-1.5 < #Delta#eta < 1.5", "le");
    legend->AddEntry(eta20peakPhiNarrow, "-1.2 < #Delta#eta < 1.2", "le");
    legend->AddEntry(eta20peakPhiNarrowest, "-1.0 < #Delta#eta < 1.0", "le");
    legend->SetBorderSize(0);

    TCanvas* ceta20peak = new TCanvas("ceta20peak", "ceta20peak", 50, 50, 800, 800);
    ceta20peak->Divide(2,2);
    ceta20peak->cd(1)->SetTheta(23);
    ceta20peak->cd(1)->SetPhi(65);
    eta20peak->Draw("SURF1");
    ceta20peak->cd(2);
    eta20peakEta->Draw("H");
    ceta20peak->cd(3);
    eta20peakPhiNarrow->GetYaxis()->SetRangeUser(0.5*eta20peakPhiNarrowest->GetMinimum(), 1.2*eta20peakPhiNarrow->GetMaximum());
    eta20peakPhi->Draw("H");
    eta20peakPhiNarrow->Draw("H");
    eta20peakPhiNarrowest->Draw("H SAME");
    legend->Draw();

    TCanvas* ceta20RSB = new TCanvas("ceta20RSB", "ceta20RSB", 60, 60, 800, 800);
    ceta20RSB->Divide(2,2);
    ceta20RSB->cd(1)->SetTheta(23);
    ceta20RSB->cd(1)->SetPhi(65);
    eta20RSB->Draw("SURF1");
    ceta20RSB->cd(2);
    eta20RSBEta->Draw("H");
    ceta20RSB->cd(3);
    //eta20RSBPhi->GetYaxis()->SetRangeUser(1.5*eta20RSBPhiNarrowest->GetMinimum(), 1.5*eta20RSBPhi->GetMaximum());
    eta20RSBPhi->Draw("H");
    eta20RSBPhiNarrow->Draw("H SAME");
    eta20RSBPhiNarrowest->Draw("H SAME");
    legend->Draw();

    TCanvas* ceta20LSB = new TCanvas("ceta20LSB", "ceta20LSB", 70, 70, 800, 800);
    ceta20LSB->Divide(2,2);
    ceta20LSB->cd(1)->SetTheta(23);
    ceta20LSB->cd(1)->SetPhi(65);
    eta20LSB->Draw("SURF1");
    ceta20LSB->cd(2);
    eta20LSBEta->Draw("H");
    ceta20LSB->cd(3);
    //eta20LSBPhi->GetYaxis()->SetRangeUser(1.5*eta20LSBPhiNarrowest->GetMinimum(), 1.5* eta20LSBPhi->GetMaximum());
    eta20LSBPhi->Draw("H");
    eta20LSBPhiNarrow->Draw("H SAME");
    eta20LSBPhiNarrowest->Draw("H SAME");
    legend->Draw();

    //fit with dipole for away
    TF2 *fit2D = new TF2("fit2D", "[7] + [7]*2.0*[0]*cos(2.0*y) + [1]*cos(y) + ([2]/([5]*[6]))*exp(-(((x - [3])^2)/(2*[5]^2) + ((y - [4])^2)/(2*[6]^2)))", -1.2, 1.2, -0.5*TMath::Pi(), 1.5*TMath::Pi());
    fit2D->SetParameters(0.0005, -0.0002, 0.000001, 0.0, 0.0, 1.0, 1.0, 0.00005);
    fit2D->SetParLimits(0, 0.00001, 0.01);
    //fit2D->FixParameter(0, 0.0);
    fit2D->SetParLimits(2, 0.0000005, 0.05);
    fit2D->SetParLimits(3, -0.1, 0.1);
    fit2D->SetParLimits(4, -0.2, 0.2);
    fit2D->SetParLimits(5, 0.01, 2.0);
    fit2D->SetParLimits(6, 0.01, 2.0);
    fit2D->SetParNames("Quadrupole", "Dipole", "jet peak Amp.", "jet peak mean #Delta#eta", "jet peak mean #Delta#varphi", "peak sigma deta", "peak sigma dphi", "flat BG");

    //fit with 1D gaus for away (include nearside also at 2pi and awayside also at -pi)

    TF2 *fitother2D = new TF2("fitother2D", "[9] + [9]*2.0*[0]*cos(2.0*y) + ([1]/([3]))*exp(-(((y - [2])^2)/(2*[3]^2))) + ([1]/([3]))*exp(-(((y - [2] + 2.0*TMath::Pi())^2)/(2*[3]^2)))+ ([4]/([7]*[8]))*exp(-(((x - [5])^2)/(2*[7]^2) + ((y - [6])^2)/(2*[8]^2))) + ([4]/([7]*[8]))*exp(-(((x - [5])^2)/(2*[7]^2) + ((y - [6] - 2.0*TMath::Pi())^2)/(2*[8]^2)))", -1.2, 1.2, -0.5*TMath::Pi(), 1.5*TMath::Pi());
    fitother2D->SetParameters(0.0005, 0.1, 3.14, 1.0, 0.0, 0.0, 1.0, 1.0, 1.0, 0.0005);
    fitother2D->SetParLimits(0, 0.00001, 0.01);
    //fitother2D->FixParameter(0, 0.0);
    fitother2D->SetParLimits(1, 0.0000005, 0.05);
    fitother2D->SetParLimits(2, 3.1, 3.2);
    fitother2D->SetParLimits(3, 0.01, 1.0);
    fitother2D->SetParLimits(4, 0.0000005, 0.05);
    fitother2D->SetParLimits(5, -0.1, 0.1);
    fitother2D->SetParLimits(6, -0.1, 0.1);
    fitother2D->SetParLimits(7, 0.01, 1.0);
    fitother2D->SetParLimits(8, 0.01, 1.0);
    fitother2D->SetParNames("Quadrupole", "awayside Amp", "awayside Mean #Delta#varphi", "awayside sigma dphi", "jet peak Amp.", "jet peak mean #Delta#eta", "jet peak mean #Delta#varphi", "peak sigma deta", "peak sigma dphi", "flat BG");


    
    TCanvas *fitCanvas = new TCanvas("fitcanvas", "fitcanvas", 80, 80, 800, 800);
    fitCanvas->cd();
    //peak->Fit(fit2D, "R0");
    peak->Fit(fitother2D, "R0");
    fit2D->SetTitle("Mult. 0-20\% 2D Correletion Fit");
    //fit2D->Draw("SURF1");
    fitother2D->SetMinimum(0.15E-3);
    fitother2D->SetMaximum(0.41E-3);
    fitother2D->Draw("SURF1");

    residual->Add(fitother2D, -1.0);
    residual->Divide(eta20peak);

    TF2 *fitParts[4];
    fitParts[0] = new TF2("fitpart0", "[0]+0.0*x", -1.2, 1.2, -0.5*TMath::Pi(), 1.5*TMath::Pi());
    fitParts[0]->SetParameter(0, (Double_t)fit2D->GetParameter(7));
    fitParts[0]->SetTitle("Constant");
    fitParts[0]->SetLineWidth(0);
    fitParts[1] = new TF2("fitpart1", "[0]*cos(2.0*y)", -1.2, 1.2, -0.5*TMath::Pi(), 1.5*TMath::Pi());
    fitParts[1]->SetParameter(0, (Double_t)fit2D->GetParameter(7)*2.0*fit2D->GetParameter(0));
    fitParts[1]->SetTitle("Quadrupole");
    fitParts[1]->SetLineWidth(0);
    fitParts[2] = new TF2("fitpart2", "[0]*cos(y)", -1.2, 1.2, -0.5*TMath::Pi(), 1.5*TMath::Pi());
    fitParts[2]->SetParameter(0, (Double_t)fit2D->GetParameter(1));
    fitParts[2]->SetTitle("Dipole (away-side jet)");
    fitParts[2]->SetLineWidth(0);
    fitParts[3] = new TF2("fitpart3", "([0]/([3]*[4]))*exp(-(((x - [1])^2)/(2*[3]^2) + ((y - [2])^2)/(2*[4]^2)))", -1.2, 1.2, -0.5*TMath::Pi(), 1.5*TMath::Pi());
    fitParts[3]->SetParameters(fit2D->GetParameter(2), fit2D->GetParameter(3), fit2D->GetParameter(4), fit2D->GetParameter(5), fit2D->GetParameter(6));
    fitParts[3]->SetTitle("2D Gaussian (near-side jet)");
    fitParts[3]->SetLineWidth(0);

    TF2 *otherfitParts[4];
    otherfitParts[0] = new TF2("fitpart0", "[0]+0.0*x", -1.2, 1.2, -0.5*TMath::Pi(), 1.5*TMath::Pi());
    otherfitParts[0]->SetParameter(0, (Double_t)fitother2D->GetParameter(9));
    otherfitParts[0]->SetTitle("Constant");
    otherfitParts[0]->SetLineWidth(0);
    otherfitParts[1] = new TF2("fitpart1", "[0]*cos(2.0*y)", -1.2, 1.2, -0.5*TMath::Pi(), 1.5*TMath::Pi());
    otherfitParts[1]->SetParameter(0, (Double_t)fitother2D->GetParameter(9)*2.0*fitother2D->GetParameter(0));
    otherfitParts[1]->SetTitle("Quadrupole");
    otherfitParts[1]->SetLineWidth(0);
    otherfitParts[2] = new TF2("fitpart2", "([0]/[2])*exp(-(y-[1])^2/(2*[2]^2))", -1.2, 1.2, -0.5*TMath::Pi(), 1.5*TMath::Pi());
    otherfitParts[2]->SetParameter(0, (Double_t)fitother2D->GetParameter(1));
    otherfitParts[2]->SetParameter(1, (Double_t)fitother2D->GetParameter(2));
    otherfitParts[2]->SetParameter(2, (Double_t)fitother2D->GetParameter(3));
    otherfitParts[2]->SetTitle("Away-side jet)");
    otherfitParts[2]->SetLineWidth(0);
    otherfitParts[3] = new TF2("fitpart3", "([0]/([3]*[4]))*exp(-(((x - [1])^2)/(2*[3]^2) + ((y - [2])^2)/(2*[4]^2)))", -1.2, 1.2, -0.5*TMath::Pi(), 1.5*TMath::Pi());
    otherfitParts[3]->SetParameters(fitother2D->GetParameter(4), fitother2D->GetParameter(5), fitother2D->GetParameter(6), fitother2D->GetParameter(7), fitother2D->GetParameter(8));
    otherfitParts[3]->SetTitle("2D Gaussian (near-side jet)");
    otherfitParts[3]->SetLineWidth(0);


    TCanvas *resCanvas = new TCanvas("resCanvas", "resCanvas", 85, 85, 800, 800);
    resCanvas->cd();
    residual->GetXaxis()->SetRangeUser(-1.2, 1.2);
    residual->Draw("SURF1");

    TCanvas *fitPartCanvas = new TCanvas("fitpartcanvas", "fitpartcanvas", 85, 85, 1200, 400);
    fitPartCanvas->Divide(4,1);
    for(int i = 0; i<4; i++){
        fitPartCanvas->cd(i+1);
        hfits[i]->Eval(fitParts[i]);
        hotherfits[i]->Eval(otherfitParts[i]);
        if(i==0){
            hfits[i]->GetZaxis()->SetLabelSize(0.07);
        }else{
            float max = 0.000030;
            float min = -0.000006;
            hfits[i]->GetZaxis()->SetRangeUser(min, max);
            hotherfits[i]->GetZaxis()->SetRangeUser(min, max);
        }
        //hfits[i]->SetRange(-1.2, -0.5*TMath::Pi(), 1.2, 1.5*TMath::Pi());
        //hfits[i]->GetZaxis()->SetRangeUser(-0.00005, 0.00005);
        hfits[i]->Draw("SURF1");
        //hotherfits[i]->Draw("SURF1");
    }

    TCanvas *otherfitPartCanvas = new TCanvas("otherfitpartcanvas", "otherfitpartcanvas", 85, 85, 1200, 400);
    otherfitPartCanvas->Divide(4,1);
    for(int i = 0; i<4; i++){
        otherfitPartCanvas->cd(i+1);
        hotherfits[i]->Draw("SURF1");
    }

    TH2D* hfitnojet = (TH2D*)hotherfits[0]->Clone("nojet");
    hfitnojet->Add(hotherfits[1]);
    //hfitnojet->Add(hfits[2]);
    
    TF2* paperfit = new TF2("paperfit", "[0]*(1+2*0.12*0.1*cos(2.0*y))", -1.2, 1.2, -0.5*TMath::Pi(), 1.5*TMath::Pi());
    paperfit->SetParameter(0, fitParts[0]->GetParameter(0));
    TH2D* hpaperfit = (TH2D*)eta20peak->Clone("hpaperfit");
    hpaperfit->Eval(paperfit);

    //TH2D* hfitnonear = (TH2D*)hfits[0]->Clone("nonear");
    //hfitnonear->Add(hfits[1]);
    //hfitnonear->Add(hfits[2]);
    TH2D* hfitnonear = (TH2D*)hotherfits[0]->Clone("nonear");
    hfitnonear->Add(hotherfits[1]);
    hfitnonear->Add(hotherfits[2]);

    TH2D* hfitnoaway = (TH2D*)hotherfits[0]->Clone("noaway");
    hfitnoaway->Add(hotherfits[1]);
    hfitnoaway->Add(hotherfits[3]);


    TH2D* hfullfit = (TH2D*)eta20peak->Clone("hfullfit");
    //hfullfit->Eval(fit2D);
    hfullfit->Eval(fitother2D);

    TPaveText *data = new TPaveText(0.5415, 0.8234, 0.9914, 0.9911, "NDC");
    data->AddText("ALICE Preliminary");
    data->AddText("p-Pb #sqrt{s_{NN}} = 5.02 TeV");
    data->AddText("0-20% Multiplicity Class (V0A)");
    data->SetFillStyle(0);
    data->SetBorderSize(0);
    data->SetTextFont(42);

    TPaveText *othertext = new TPaveText(0.0057, 0.8264, 0.4054, 0.9941, "NDC");
    //othertext->AddText("h-#phi Correlation");
    othertext->AddText("4.0 < #it{p}^{h}_{T,trig} < 8.0 GeV/#it{c}");
    othertext->AddText("2.0 < #it{p}^{#phi}_{T,assoc} < 4.0 GeV/#it{c}");
    othertext->SetFillStyle(0);
    othertext->SetBorderSize(0);
    othertext->SetTextFont(42);


    TCanvas *peakCanvas = new TCanvas("peakCanvas", "peakCanvas", 90, 90, 700, 700);
    //peakCanvas->cd()->SetTheta(23.0);
    //peakCanvas->cd()->SetPhi(65.0);
    peakCanvas->cd()->SetTheta(16.2);
    peakCanvas->cd()->SetPhi(-65.0);
    peakCanvas->cd()->SetMargin(0.15, 0.05, 0.12, 0.18);
    peak->GetXaxis()->SetRangeUser(-1.2, 1.2);
    peak->Scale(1.0/(peak->GetXaxis()->GetBinWidth(1)*peak->GetYaxis()->GetBinWidth(1)));
    hfullfit->GetXaxis()->SetRangeUser(-1.0, 1.0);
    hfullfit->SetLineColor(kRed);
    hfullfit->SetLineWidth(2);
    hfitnojet->GetXaxis()->SetRangeUser(-1.0, 1.0);
    hfitnojet->SetLineColor(kRed);
    hfitnonear->GetXaxis()->SetRangeUser(-1.0, 1.0);
    hfitnonear->SetLineColor(kMagenta);
    //peak->GetZaxis()->SetRangeUser(0.000040, 0.000090);
    hfullfit->GetZaxis()->SetRangeUser(0.000040, 0.000090);
    hfitnojet->GetZaxis()->SetRangeUser(0.000040, 0.000090);
    hpaperfit->GetZaxis()->SetRangeUser(0.000040, 0.000090);
    hfitnonear->GetZaxis()->SetRangeUser(0.000040, 0.000090);
    peak->SetTitle("");
    peak->GetXaxis()->SetTitle("#Delta#eta");
    peak->GetXaxis()->SetTitleOffset(1.5);
    peak->GetXaxis()->SetTitleSize(0.05);
    peak->GetXaxis()->CenterTitle();
    peak->GetYaxis()->SetTitle("#Delta#varphi");
    peak->GetYaxis()->SetTitleOffset(1.5);
    peak->GetYaxis()->SetTitleSize(0.05);
    peak->GetYaxis()->CenterTitle();
    peak->GetZaxis()->SetTitle("1/N_{trig} d^{2}#it{N}_{assoc}/d#it{#Delta#varphi}d#it{#Delta#eta}");
    peak->GetZaxis()->SetTitleOffset(1.7);
    peak->GetZaxis()->SetNdivisions(409);
    peak->GetZaxis()->SetMaxDigits(2);
    //TGaxis::SetExponentOffset(0.25, 0.05, "xy");
    peak->Draw("SAME SURF1");
    othertext->Draw();
    data->Draw();
    //hfullfit->Draw("SAME SURF");
    //hpaperfit->Draw("SAME SURF");
    //hfitnojet->Draw("SAME SURF");
    //hfitnonear->Draw("SAME SURF");
    
    //Project fitparts and correlation onto just dphi
    TH1D* hnojet1D = hfitnojet->ProjectionY("hnojet1D", hfitnojet->GetXaxis()->FindBin(-1.2 + epsilon), hfitnojet->GetXaxis()->FindBin(1.2 - epsilon));
    hnojet1D->SetLineWidth(4);
    hnojet1D->SetLineStyle(7);
    hnojet1D->SetLineColor(kBlue);

    TH1D* hnonear1D = hfitnonear->ProjectionY("hnonear1D", hfitnonear->GetXaxis()->FindBin(-1.2 + epsilon), hfitnonear->GetXaxis()->FindBin(1.2 - epsilon));
    hnonear1D->SetLineWidth(4);
    hnonear1D->SetLineStyle(6);
    hnonear1D->SetLineColor(kGreen+2);

    TH1D* hfullfit1D = hfullfit->ProjectionY("hfullfit1D", hfullfit->GetXaxis()->FindBin(-1.2 + epsilon), hfullfit->GetXaxis()->FindBin(1.2 - epsilon));
    hfullfit1D->SetLineWidth(4);
    hfullfit1D->SetLineStyle(6);
    hfullfit1D->SetLineColor(kGreen+2);

    TH1D* hpeak1D = peak->ProjectionY("hpeak1D", peak->GetXaxis()->FindBin(-1.2 + epsilon), peak->GetXaxis()->FindBin(1.2 - epsilon));
    hpeak1D->SetLineWidth(3);
    TF1* line = new TF1("line", "[0]", -1.57, 4.71);
    line->SetParameter(0, 0.25*(hpeak1D->GetBinContent(8)+hpeak1D->GetBinContent(9)+hpeak1D->GetBinContent(16)+hpeak1D->GetBinContent(1)));
    line->SetLineStyle(4);
    line->SetLineWidth(4);
    line->SetLineColor(kRed);

    TCanvas *c1D = new TCanvas("c1D", "c1D", 65, 65, 800, 600);
    c1D->cd();
    hpeak1D->GetYaxis()->SetRangeUser(0.0004, 0.0007);
    hpeak1D->Draw("P E");
    hnojet1D->Draw("HIST C SAME");
    //hnonear1D->Draw("HIST C SAME");
    hfullfit1D->Draw("HIST C SAME");
    line->Draw("SAME");

    Float_t aboveline = hpeak1D->Integral(1, 8) - line->GetParameter(0)*8.0;
    Float_t abovefit = hpeak1D->Integral(1, 8) - hnojet1D->Integral(1, 8);

    Float_t awayaboveline = hpeak1D->Integral(9, 16) - line->GetParameter(0)*8.0;
    Float_t awayabovefit = hpeak1D->Integral(9, 16) - hnojet1D->Integral(9, 16);

    TString intline(Form("Near Peak Integral above Flat BG: %E", aboveline));
    TString intfit(Form("Near Peak Integral above 2D Fit BG: %E", abovefit));
    TString awayintline(Form("Away Peak Integral above Flat BG: %E", awayaboveline));
    TString awayintfit(Form("Away Peak Integral above 2D Fit BG: %E", awayabovefit));


    TPaveText* text =  new TPaveText(0.345, 0.731, 0.891, 0.89, "NDC");
    text->AddText(intline.Data());
    text->AddText(intfit.Data());
    text->AddText("");
    text->AddText(awayintline.Data());
    text->AddText(awayintfit.Data());
    text->SetTextAlign(32);
    text->SetShadowColor(kWhite);
    text->SetBorderSize(3);
    text->Draw();

    TH1D* dphiuncorr = uncorr2D->ProjectionY("dphiuncorr", uncorr2D->GetXaxis()->FindBin(-1.2 + epsilon), uncorr2D->GetXaxis()->FindBin(1.2 - epsilon));
    TH1D* dphiBG = BG2D->ProjectionY("dphiBG", BG2D->GetXaxis()->FindBin(-1.2 + epsilon), BG2D->GetXaxis()->FindBin(1.2 - epsilon));



    TCanvas* cBGcompare = new TCanvas("cBGcompare", "cBGcompare", 50, 50, 600, 600);
    cBGcompare->cd()->SetMargin(0.12, 0.03, 0.1, 0.04);
    dphiuncorr->SetLineColor(kBlue+1);
    dphiuncorr->SetLineWidth(2);
    dphiuncorr->SetTitle("");
    dphiuncorr->GetXaxis()->SetTitle("#Delta#varphi");
    dphiuncorr->GetYaxis()->SetMaxDigits(3);
    dphiuncorr->GetYaxis()->SetTitle("1/#it{N}_{trig} d#it{N}_{assoc}/d#it{#Delta#varphi}");
    dphiuncorr->Scale(1.0/(0.49*0.82*dphiuncorr->GetXaxis()->GetBinWidth(1)));
    dphiuncorr->Draw("HIST E");
    dphiBG->SetLineColor(kRed+1);
    dphiBG->SetLineWidth(2);
    dphiBG->Scale(1.0/(0.49*0.82*dphiBG->GetXaxis()->GetBinWidth(1)));
    dphiBG->Draw("SAME");
    TLegend* bglegend = new TLegend(0.4, 0.6, 0.5, 0.7);
    bglegend->AddEntry(dphiBG, "BG Estimate", "le");
    bglegend->AddEntry(dphiuncorr, "Signal + BG", "le");
    bglegend->SetLineWidth(0);
    bglegend->Draw();
    data->Draw();
    othertext->Draw();

   }
