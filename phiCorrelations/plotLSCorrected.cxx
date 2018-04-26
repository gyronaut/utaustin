void plotLSCorrected(string inputname){
    gStyle->SetOptStat(0);
    gStyle->SetOptFit(1);

    TFile* eta20File = new TFile(inputname.c_str());

    TH2D* eta20peak = RLSsubhPhi2Dpeak->Clone("eta20peak");
    TH2D* eta20RSB = RLSsubhPhi2DRside->Clone("eta20RSB");
    TH2D* eta20LSB = LLSsubhPhi2DLside->Clone("eta20LSB");

    TH2D* residual = (TH2D*)eta20peak->Clone("residual");
    TH2D* peak = (TH2D*)eta20peak->Clone("peak");
    TH2D* hfits[4];
    for(int i =0; i<4; i++){
        hfits[i] = (TH2D*)eta20peak->Clone(Form("hfit%d", i));
    }

    eta20peak->GetXaxis()->SetTitle("#Delta#eta");
    eta20peak->GetXaxis()->SetTitleSize(0.05);
    eta20peak->GetXaxis()->SetTitleOffset(1.3);
    eta20peak->GetYaxis()->SetTitle("#Delta#varphi");
    eta20peak->GetYaxis()->SetTitleSize(0.05);
    eta20peak->GetYaxis()->SetTitleOffset(1.3);
    eta20peak->SetTitle("");
    //eta20peak->SetStats(kFALSE);
    //eta20peak->Scale(1.0/(eta20peak->Integral(eta20peak->GetXaxis()->FindBin(-1.2), eta20peak->GetXaxis()->FindBin(1.2), 1, eta20peak->GetYaxis()->GetNbins())));

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
    TH1D* eta20peakPhi = eta20peak->ProjectionY("eta20peakPhi", eta20peak->GetXaxis()->FindBin(-1.5), eta20peak->GetXaxis()->FindBin(1.5));
    eta20peakPhi->GetXaxis()->SetTitleOffset(1.0);
    eta20peakPhi->SetStats(kFALSE);
    eta20peakPhi->SetLineColor(kBlue);
    TH1D* eta20peakPhiNarrow = eta20peak->ProjectionY("eta20peakPhiNarrow", eta20peak->GetXaxis()->FindBin(-1.2), eta20peak->GetXaxis()->FindBin(1.2));
    eta20peakPhiNarrow->SetLineColor(kViolet);
    eta20peakPhiNarrow->SetStats(kFALSE);
    eta20peakPhiNarrow->GetXaxis()->SetTitleOffset(1.0);
    TH1D* eta20peakPhiNarrowest = eta20peak->ProjectionY("eta20peakPhiNarrowest", eta20peak->GetXaxis()->FindBin(-1.0), eta20peak->GetXaxis()->FindBin(1.0));
    eta20peakPhiNarrowest->SetLineColor(kRed);

    TH1D* eta20RSBEta = eta20RSB->ProjectionX("eta20RSBEta", 1, eta20RSB->GetYaxis()->GetNbins());
    eta20RSBEta->GetXaxis()->SetTitleOffset(1.0);
    eta20RSBEta->GetXaxis()->SetRangeUser(-1.2, 1.2);
    eta20RSBEta->SetStats(kFALSE);
    TH1D* eta20RSBPhi = eta20RSB->ProjectionY("eta20RSBPhi", eta20RSB->GetXaxis()->FindBin(-1.5), eta20RSB->GetXaxis()->FindBin(1.5));
    eta20RSBPhi->GetXaxis()->SetTitleOffset(1.0);
    eta20RSBPhi->SetLineColor(kBlue);
    eta20RSBPhi->SetStats(kFALSE);
    TH1D* eta20RSBPhiNarrow = eta20RSB->ProjectionY("eta20RSBPhiNarrow", eta20RSB->GetXaxis()->FindBin(-1.2), eta20RSB->GetXaxis()->FindBin(1.2));
    eta20RSBPhiNarrow->SetLineColor(kViolet);
    TH1D* eta20RSBPhiNarrowest = eta20RSB->ProjectionY("eta20RSBPhiNarrowest", eta20RSB->GetXaxis()->FindBin(-1.0), eta20RSB->GetXaxis()->FindBin(1.0));
    eta20RSBPhiNarrowest->SetLineColor(kRed);

    TH1D* eta20LSBEta = eta20LSB->ProjectionX("eta20LSBEta", 1, eta20LSB->GetYaxis()->GetNbins());
    eta20LSBEta->GetXaxis()->SetTitleOffset(1.0);
    eta20LSBEta->GetXaxis()->SetRangeUser(-1.2, 1.2);
    eta20LSBEta->SetStats(kFALSE);
    TH1D* eta20LSBPhi = eta20LSB->ProjectionY("eta20LSBPhi", eta20LSB->GetXaxis()->FindBin(-1.5), eta20LSB->GetXaxis()->FindBin(1.5));
    eta20LSBPhi->GetXaxis()->SetTitleOffset(1.0);
    eta20LSBPhi->SetLineColor(kBlue);
    eta20LSBPhi->SetStats(kFALSE);
    TH1D* eta20LSBPhiNarrow = eta20LSB->ProjectionY("eta20LSBPhiNarrow", eta20LSB->GetXaxis()->FindBin(-1.2), eta20LSB->GetXaxis()->FindBin(1.2));
    eta20LSBPhiNarrow->SetLineColor(kViolet);
    TH1D* eta20LSBPhiNarrowest = eta20LSB->ProjectionY("eta20LSBPhiNarrowest", eta20LSB->GetXaxis()->FindBin(-1.0), eta20LSB->GetXaxis()->FindBin(1.0));
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
    ceta20peak->cd(1)->SetTheta(50);
    ceta20peak->cd(1)->SetPhi(50);
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
    ceta20RSB->cd(1)->SetTheta(50);
    ceta20RSB->cd(1)->SetPhi(50);
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
    ceta20LSB->cd(1)->SetTheta(50);
    ceta20LSB->cd(1)->SetPhi(50);
    eta20LSB->Draw("SURF1");
    ceta20LSB->cd(2);
    eta20LSBEta->Draw("H");
    ceta20LSB->cd(3);
    //eta20LSBPhi->GetYaxis()->SetRangeUser(1.5*eta20LSBPhiNarrowest->GetMinimum(), 1.5* eta20LSBPhi->GetMaximum());
    eta20LSBPhi->Draw("H");
    eta20LSBPhiNarrow->Draw("H SAME");
    eta20LSBPhiNarrowest->Draw("H SAME");
    legend->Draw();

    TF2 *fit2D = new TF2("fit2D", "[7] + [7]*2.0*[0]*cos(2.0*y) + [1]*cos(y) + ([2]/([5]*[6]))*exp(-(((x - [3])^2)/(2*[5]^2) + ((y - [4])^2)/(2*[6]^2)))", -1.2, 1.2, -0.5*TMath::Pi(), 1.5*TMath::Pi());
    fit2D->SetParameters(0.0005, -0.0002, 0.000001, 0.0, 0.0, 1.0, 1.0, 0.00005);
    fit2D->SetParLimits(0, 0.00001, 0.01);
    //fit2D->FixParameter(0, 0.0);
    fit2D->SetParLimits(2, 0.0000005, 0.05);
    fit2D->SetParLimits(3, -0.4, 0.4);
    fit2D->SetParLimits(4, -0.2, 0.2);
    fit2D->SetParLimits(5, 0.01, 2.0);
    fit2D->SetParLimits(6, 0.01, 2.0);
    fit2D->SetParNames("Quadrupole", "Dipole", "jet peak Amp.", "jet peak mean #Delta#eta", "jet peak mean #Delta#varphi", "peak sigma deta", "peak sigma dphi", "flat BG");

    
    TCanvas *fitCanvas = new TCanvas("fitcanvas", "fitcanvas", 80, 80, 800, 800);
    fitCanvas->cd();
    peak->Fit(fit2D, "R0");
    fit2D->SetTitle("Mult. 0-20\% 2D Correletion Fit");
    fit2D->Draw("SURF1");

    residual->Add(fit2D, -1.0);
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

    TCanvas *resCanvas = new TCanvas("resCanvas", "resCanvas", 85, 85, 800, 800);
    resCanvas->cd();
    residual->GetXaxis()->SetRangeUser(-1.2, 1.2);
    residual->Draw("SURF1");

    TCanvas *fitPartCanvas = new TCanvas("fitpartcanvas", "fitpartcanvas", 85, 85, 1200, 400);
    fitPartCanvas->Divide(4,1);
    for(int i = 0; i<4; i++){
        fitPartCanvas->cd(i+1);
        hfits[i]->Eval(fitParts[i]);
        if(i==0){
            hfits[i]->GetZaxis()->SetLabelSize(0.07);
        }else{
            float max = 0.00003;
            float min = -0.000006;
            hfits[i]->GetZaxis()->SetRangeUser(min, max);
        }
        //hfits[i]->SetRange(-1.2, -0.5*TMath::Pi(), 1.2, 1.5*TMath::Pi());
        //hfits[i]->GetZaxis()->SetRangeUser(-0.00005, 0.00005);
        hfits[i]->Draw("SURF1");
    }

    hfitnojet = (TH2D*)hfits[0]->Clone("nojet");
    hfitnojet->Add(hfits[1]);
    //hfitnojet->Add(hfits[2]);
    
    TF2* paperfit = new TF2("paperfit", "[0]*(1+2*0.12*0.1*cos(2.0*y))", -1.2, 1.2, -0.5*TMath::Pi(), 1.5*TMath::Pi());
    paperfit->SetParameter(0, fitParts[0]->GetParameter(0));
    hpaperfit = (TH2D*)eta20peak->Clone("hpaperfit");
    hpaperfit->Eval(paperfit);

    TH2D* hfitnonear = (TH2D*)hfits[0]->Clone("nonear");
    hfitnonear->Add(hfits[1]);
    hfitnonear->Add(hfits[2]);


    TH2D* hfullfit = (TH2D*)eta20peak->Clone("hfullfit");
    hfullfit->Eval(fit2D);

    TCanvas *peakCanvas = new TCanvas("peakCanvas", "peakCanvas", 90, 90, 800, 800);
    peakCanvas->cd()->SetTheta(50.0);
    peakCanvas->cd()->SetPhi(50.0);
    peak->GetXaxis()->SetRangeUser(-1.0, 1.0);
    hfullfit->GetXaxis()->SetRangeUser(-1.0, 1.0);
    hfullfit->SetLineColor(kRed);
    hfullfit->SetLineWidth(2);
    hfitnojet->GetXaxis()->SetRangeUser(-1.0, 1.0);
    hfitnojet->SetLineColor(kRed);
    hfitnonear->GetXaxis()->SetRangeUser(-1.0, 1.0);
    hfitnonear->SetLineColor(kMagenta);
    peak->GetZaxis()->SetRangeUser(0.000035, 0.00008);
    hfullfit->GetZaxis()->SetRangeUser(0.000035, 0.00008);
    hfitnojet->GetZaxis()->SetRangeUser(0.000035, 0.00008);
    hpaperfit->GetZaxis()->SetRangeUser(0.000035, 0.00008);
    hfitnonear->GetZaxis()->SetRangeUser(0.000035, 0.00008);
    peak->SetTitle("");
    peak->GetXaxis()->SetTitle("#Delta#eta");
    peak->GetXaxis()->SetTitleOffset(1.4);
    peak->GetYaxis()->SetTitle("#Delta#varphi");
    peak->GetYaxis()->SetTitleOffset(1.3);
    peak->Draw("SAME SURF1");
    //hfullfit->Draw("SAME SURF");
    //hpaperfit->Draw("SAME SURF");
    //hfitnojet->Draw("SAME SURF");
    //hfitnonear->Draw("SAME SURF"); 
   }
