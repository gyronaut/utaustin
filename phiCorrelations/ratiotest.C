{

    auto canvas = new TCanvas("c1", "c1", 50, 50 , 900, 600);
    canvas->cd();

    gROOT->SetStyle("Plain");
    double p9065_d30x1y1_xval[] = { 2.5, 7.5, 15.0, 30.0, 50.0, 70.0, 90.0 };
    double p9065_d30x1y1_xerrminus[] = { 2.5, 2.5, 5.0, 10.0, 10.0, 10.0, 10.0 };
    double p9065_d30x1y1_xerrplus[] = { 2.5, 2.5, 5.0, 10.0, 10.0, 10.0, 10.0 };
    double p9065_d30x1y1_yval[] = { 0.0185, 0.01743, 0.01741, 0.01704, 0.01612, 0.01474, 0.01426 };
    double p9065_d30x1y1_yerrminus[] = { 0.0016835082417380675, 0.0014335968749965941, 0.0013756816492197603, 0.0013773525329413671, 0.0013773525329413671, 0.0012300812981262661, 0.0016333401360402555 };
    double p9065_d30x1y1_yerrplus[] = { 0.0016835082417380675, 0.0014335968749965941, 0.0013756816492197603, 0.0013773525329413671, 0.0013773525329413671, 0.0012300812981262661, 0.0016333401360402555 };
    double p9065_d30x1y1_ystatminus[] = { 1.8E-4, 1.8E-4, 1.4E-4, 1.1E-4, 1.1E-4, 1.5E-4, 2.1E-4 };
    double p9065_d30x1y1_ystatplus[] = { 1.8E-4, 1.8E-4, 1.4E-4, 1.1E-4, 1.1E-4, 1.5E-4, 2.1E-4 };
    int p9065_d30x1y1_numpoints = 7;
    
    double p9065_d30x1y2_xval[] = { 2.5, 7.5, 15.0, 30.0, 50.0, 70.0, 90.0 };
    double p9065_d30x1y2_xerrminus[] = { 2.5, 2.5, 5.0, 10.0, 10.0, 10.0, 10.0 };
    double p9065_d30x1y2_xerrplus[] = { 2.5, 2.5, 5.0, 10.0, 10.0, 10.0, 10.0 };
    double p9065_d30x1y2_yval[] = { 0.129, 0.1241, 0.1254, 0.125, 0.12132, 0.1143, 0.116 };
    double p9065_d30x1y2_yerrminus[] = { 0.014771932845772079, 0.012634080892569906, 0.012251122397560151, 0.011967455869983393, 0.011520598942763352, 0.010078690391117289, 0.013363008643265931 };
    double p9065_d30x1y2_yerrplus[] = { 0.014771932845772079, 0.012634080892569906, 0.012251122397560151, 0.011967455869983393, 0.011520598942763352, 0.010078690391117289, 0.013363008643265931 };
    double p9065_d30x1y2_ystatminus[] = { 0.0013, 0.0013, 0.001, 8.0E-4, 8.5E-4, 0.0011, 0.0018 };
    double p9065_d30x1y2_ystatplus[] = { 0.0013, 0.0013, 0.001, 8.0E-4, 8.5E-4, 0.0011, 0.0018 };
    int p9065_d30x1y2_numpoints = 7;
    
    double p9065_d30x1y3_xval[] = { 2.5, 7.5, 15.0, 30.0, 50.0, 70.0, 90.0 };
    double p9065_d30x1y3_xerrminus[] = { 2.5, 2.5, 5.0, 10.0, 10.0, 10.0, 10.0 };
    double p9065_d30x1y3_xerrplus[] = { 2.5, 2.5, 5.0, 10.0, 10.0, 10.0, 10.0 };
    double p9065_d30x1y3_yval[] = { 0.3311, 0.311, 0.3096, 0.3032, 0.286, 0.2609, 0.2666 };
    double p9065_d30x1y3_yerrminus[] = { 0.03469884724310017, 0.030316002374983414, 0.029192122225011323, 0.028818917398125837, 0.027386127875258306, 0.024718818742003025, 0.03304330491945381 };
    double p9065_d30x1y3_yerrplus[] = { 0.03469884724310017, 0.030316002374983414, 0.029192122225011323, 0.028818917398125837, 0.027386127875258306, 0.024718818742003025, 0.03304330491945381 };
    double p9065_d30x1y3_ystatminus[] = { 0.0033, 0.0033, 0.0025, 0.0019, 0.002, 0.0026, 0.0041 };
    double p9065_d30x1y3_ystatplus[] = { 0.0033, 0.0033, 0.0025, 0.0019, 0.002, 0.0026, 0.0041 };
    int p9065_d30x1y3_numpoints = 7;
    auto p9065_d30x1y3 = TGraphAsymmErrors(p9065_d30x1y3_numpoints, p9065_d30x1y3_xval, p9065_d30x1y3_yval, p9065_d30x1y3_xerrminus, p9065_d30x1y3_xerrplus, p9065_d30x1y3_yerrminus, p9065_d30x1y3_yerrplus);
    p9065_d30x1y3.SetName("/HepData/9065/d30x1y3");
    p9065_d30x1y3.SetTitle("/HepData/9065/d30x1y3");

 
    double p9065_d30x1y1_xvalrev[7];
    double p9065_d30x1y1_xerrminusrev[7];
    double p9065_d30x1y1_xerrplusrev[7];
    double p9065_d30x1y1_yvalrev[7];
    double p9065_d30x1y1_yerrminusrev[7];
    double p9065_d30x1y1_yerrplusrev[7];
    double p9065_d30x1y1_ystatminusrev[7];
    double p9065_d30x1y1_ystatplusrev[7];

    double p9065_d30x1y2_xvalrev[7];
    double p9065_d30x1y2_xerrminusrev[7];
    double p9065_d30x1y2_xerrplusrev[7];
    double p9065_d30x1y2_yvalrev[7];
    double p9065_d30x1y2_yerrminusrev[7];
    double p9065_d30x1y2_yerrplusrev[7];
    double p9065_d30x1y2_ystatminusrev[7];
    double p9065_d30x1y2_ystatplusrev[7];

    double p9065_d30x1y3_xvalrev[7];
    double p9065_d30x1y3_xerrminusrev[7];
    double p9065_d30x1y3_xerrplusrev[7];
    double p9065_d30x1y3_yvalrev[7];
    double p9065_d30x1y3_yerrminusrev[7];
    double p9065_d30x1y3_yerrplusrev[7];
    double p9065_d30x1y3_ystatminusrev[7];
    double p9065_d30x1y3_ystatplusrev[7];

    for(int i = 0; i < 7; i++){
        p9065_d30x1y1_yvalrev[i] = p9065_d30x1y1_yval[i]*12.0;
        p9065_d30x1y1_yerrminusrev[i] = p9065_d30x1y1_yerrminus[i]*12.0;
        p9065_d30x1y1_yerrplusrev[i] = p9065_d30x1y1_yerrplus[i]*12.0;
        p9065_d30x1y2_yvalrev[i] = p9065_d30x1y2_yval[i]*2.0;
        p9065_d30x1y2_yerrminusrev[i] = p9065_d30x1y2_yerrminus[i]*2.0;
        p9065_d30x1y2_yerrplusrev[i] = p9065_d30x1y2_yerrplus[i]*2.0;
        p9065_d30x1y3_yvalrev[i] = p9065_d30x1y3_yval[i];
        p9065_d30x1y3_yerrminusrev[i] = p9065_d30x1y3_yerrminus[i];
        p9065_d30x1y3_yerrplusrev[i] = p9065_d30x1y3_yerrplus[i];
        
        p9065_d30x1y1_xvalrev[i] = 100.0 - p9065_d30x1y1_xval[i];
        p9065_d30x1y1_xvalrev[i] = 100.0 - p9065_d30x1y1_xval[i];
        p9065_d30x1y1_xvalrev[i] = 100.0 - p9065_d30x1y1_xval[i];

        p9065_d30x1y2_xvalrev[i] = 100.0 - p9065_d30x1y2_xval[i];
        p9065_d30x1y2_xvalrev[i] = 100.0 - p9065_d30x1y2_xval[i];
        p9065_d30x1y2_xvalrev[i] = 100.0 - p9065_d30x1y2_xval[i];

        p9065_d30x1y3_xvalrev[i] = 100.0 - p9065_d30x1y3_xval[i];
        p9065_d30x1y3_xvalrev[i] = 100.0 - p9065_d30x1y3_xval[i];
        p9065_d30x1y3_xvalrev[i] = 100.0 - p9065_d30x1y3_xval[i];

    }


    TH1F* hist = new TH1F("hist", "#phi(1020) Ratios to Stable Particles in p-Pb Collisions", 20, 0, 100);
    hist->GetXaxis()->SetTitle("Multiplicity Pctl.");
    hist->GetXaxis()->SetTitleOffset(1.3);
    hist->GetXaxis()->SetTickLength(0);
    hist->GetXaxis()->SetLabelOffset(999);
    hist->GetYaxis()->SetTitle("Particle Ratio");
    hist->GetYaxis()->SetRangeUser(0.15, 0.38);
    hist->SetNdivisions(405);
    hist->SetStats(kFALSE);
    hist->Draw();

    gPad->Update();
    TGaxis *newaxis = new TGaxis(gPad->GetUxmax(),
            gPad->GetUymin(),
            gPad->GetUxmin(),
            gPad->GetUymin(),
            hist->GetXaxis()->GetXmin(),
            hist->GetXaxis()->GetXmax(),
            405,"-");
    newaxis->SetLabelOffset(-0.03);
    //newaxis->SetTitle("Multipliciy % (VOA)");
    //newaxis->SetTitleOffset(1.3);
    newaxis->Draw();

    auto p9065_d30x1y1 = TGraphAsymmErrors(p9065_d30x1y1_numpoints, p9065_d30x1y1_xvalrev, p9065_d30x1y1_yvalrev, p9065_d30x1y1_xerrminus, p9065_d30x1y1_xerrplus, p9065_d30x1y1_yerrminusrev, p9065_d30x1y1_yerrplusrev);
    p9065_d30x1y1.SetName("/HepData/9065/d30x1y1");
    p9065_d30x1y1.SetTitle("#phi(1020) Ratios to Stable Particles in p-Pb Collisions");
    p9065_d30x1y1.GetXaxis()->SetTitle("Multiplicity Pctl.");
    p9065_d30x1y1.GetYaxis()->SetTitle("Particle Ratio");
    p9065_d30x1y1.GetYaxis()->SetRangeUser(0.15, 0.38);
    p9065_d30x1y1.SetMarkerStyle(20);
    p9065_d30x1y1.SetMarkerSize(2);
    p9065_d30x1y1.SetMarkerColor(kRed+1);
    p9065_d30x1y1.SetLineColor(kRed+1);
    p9065_d30x1y1.Draw("P");

    auto p9065_d30x1y2 = TGraphAsymmErrors(p9065_d30x1y2_numpoints, p9065_d30x1y2_xvalrev, p9065_d30x1y2_yvalrev, p9065_d30x1y2_xerrminus, p9065_d30x1y2_xerrplus, p9065_d30x1y2_yerrminusrev, p9065_d30x1y2_yerrplusrev);
    p9065_d30x1y2.SetName("/HepData/9065/d30x1y2");
    p9065_d30x1y2.SetTitle("/HepData/9065/d30x1y2");
    p9065_d30x1y2.SetMarkerStyle(21);
    p9065_d30x1y2.SetMarkerSize(2);
    p9065_d30x1y2.SetMarkerColor(kBlue+2);
    p9065_d30x1y2.SetLineColor(kBlue+2);
    p9065_d30x1y2.Draw("P");

    auto p9065_d30x1y3 = TGraphAsymmErrors(p9065_d30x1y3_numpoints, p9065_d30x1y3_xvalrev, p9065_d30x1y3_yvalrev, p9065_d30x1y3_xerrminus, p9065_d30x1y3_xerrplus, p9065_d30x1y3_yerrminusrev, p9065_d30x1y3_yerrplusrev);
    p9065_d30x1y3.SetName("/HepData/9065/d30x1y3");
    p9065_d30x1y3.SetTitle("/HepData/9065/d30x1y3");
    p9065_d30x1y3.SetMarkerStyle(23);
    p9065_d30x1y3.SetMarkerSize(2);
    p9065_d30x1y3.SetMarkerColor(kGreen+2);
    p9065_d30x1y3.SetLineColor(kGreen+2);
    p9065_d30x1y3.Draw("P");


    TLegend* legend = new TLegend(0.134, 0.667, 0.424, 0.8866);
    legend->AddEntry(&p9065_d30x1y1, "2#phi/(#pi^{+} + #pi^{-}) (x 12)", "p");
    legend->AddEntry(&p9065_d30x1y2, "2#phi/(K^{+} + K^{-}) (x 2)", "p");
    legend->AddEntry(&p9065_d30x1y3, "2#phi/(p + #bar{p})", "p");
    legend->Draw("SAME");
 
    value1 = p9065_d30x1y1_yval[0];
    value2 = p9065_d30x1y2_yval[0];
    value3 = p9065_d30x1y3_yval[0];

    for(int i=0; i<7; i++){
        p9065_d30x1y1_yval[i] *= 1.0/value1;
        p9065_d30x1y2_yval[i] *= 1.0/value2;
        p9065_d30x1y3_yval[i] *= 1.0/value3;
    }

    auto ratio1 = TGraph(7, p9065_d30x1y1_xval, p9065_d30x1y1_yval);
    auto ratio2 = TGraph(7, p9065_d30x1y2_xval, p9065_d30x1y2_yval);
    auto ratio3 = TGraph(7, p9065_d30x1y3_xval, p9065_d30x1y3_yval);

    auto canvas2 = TCanvas("c2", "c2", 50, 50, 800, 600);
    canvas2.cd();
   
    ratio1.SetTitle("#phi(1020) Ratios to Stable Particles in p-Pb Collisions");
    ratio1.GetXaxis()->SetTitle("Multiplicity");
    ratio1.GetYaxis()->SetTitle("(Particle Ratio)/(0-5% Ratio)");
    ratio1.GetYaxis()->SetRangeUser(0.5, 1.5);
    ratio1.SetMarkerStyle(20);
    ratio1.SetMarkerSize(2);
    ratio1.SetMarkerColor(kRed+1);
    ratio1.SetLineColor(kRed+1);
    ratio1.Draw("AP");


    ratio2.SetMarkerStyle(21);
    ratio2.SetMarkerSize(2);
    ratio2.SetMarkerColor(kBlue+2);
    ratio2.SetLineColor(kBlue+2);
    ratio2.Draw("P");


    ratio3.SetMarkerStyle(23);
    ratio3.SetMarkerSize(2);
    ratio3.SetMarkerColor(kGreen+2);
    ratio3.SetLineColor(kGreen+2);
    ratio3.Draw("P");

    legend->Draw("SAME");


}

    
    auto p9065_d30x1y1 = TGraphAsymmErrors(p9065_d30x1y1_numpoints, p9065_d30x1y1_xval, p9065_d30x1y1_yval, p9065_d30x1y1_xerrminus, p9065_d30x1y1_xerrplus, p9065_d30x1y1_yerrminus, p9065_d30x1y1_yerrplus);
    p9065_d30x1y1.SetName("/HepData/9065/d30x1y1");
    p9065_d30x1y1.SetTitle("#phi(1020) Ratios to Stable Particles in p-Pb Collisions");
    p9065_d30x1y1.GetXaxis()->SetTitle("Multiplicity");
    p9065_d30x1y1.GetYaxis()->SetTitle("Particle Ratio");
    p9065_d30x1y1.GetYaxis()->SetRangeUser(0.15, 0.38);
    p9065_d30x1y1.SetMarkerStyle(20);
    p9065_d30x1y1.SetMarkerSize(2);
    p9065_d30x1y1.SetMarkerColor(kRed+1);
    p9065_d30x1y1.SetLineColor(kRed+1);
    p9065_d30x1y1.Draw("AP");

    auto p9065_d30x1y2 = TGraphAsymmErrors(p9065_d30x1y2_numpoints, p9065_d30x1y2_xval, p9065_d30x1y2_yval, p9065_d30x1y2_xerrminus, p9065_d30x1y2_xerrplus, p9065_d30x1y2_yerrminus, p9065_d30x1y2_yerrplus);
    p9065_d30x1y2.SetName("/HepData/9065/d30x1y2");
    p9065_d30x1y2.SetTitle("/HepData/9065/d30x1y2");
    p9065_d30x1y2.SetMarkerStyle(21);
    p9065_d30x1y2.SetMarkerSize(2);
    p9065_d30x1y2.SetMarkerColor(kBlue+2);
    p9065_d30x1y2.SetLineColor(kBlue+2);
    p9065_d30x1y2.Draw("P");

    p9065_d30x1y3.SetMarkerStyle(23);
    p9065_d30x1y3.SetMarkerSize(2);
    p9065_d30x1y3.SetMarkerColor(kGreen+2);
    p9065_d30x1y3.SetLineColor(kGreen+2);
    p9065_d30x1y3.Draw("P");


    TLegend* legend = new TLegend(0.2, 0.5, 0.5, 0.9);
    legend->AddEntry(&p9065_d30x1y1, "2#phi/(#pi^{+} + #pi^{-}) (x 12)", "p");
    legend->AddEntry(&p9065_d30x1y2, "2#phi/(K^{+} + K^{-}) (x 2)", "p");
    legend->AddEntry(&p9065_d30x1y3, "2#phi/(p + #bar{p})", "p");
    legend->Draw("SAME");
 
    value1 = p9065_d30x1y1_yval[0];
    value2 = p9065_d30x1y2_yval[0];
    value3 = p9065_d30x1y3_yval[0];

    for(int i=0; i<7; i++){
        p9065_d30x1y1_yval[i] *= 1.0/value1;
        p9065_d30x1y2_yval[i] *= 1.0/value2;
        p9065_d30x1y3_yval[i] *= 1.0/value3;
    }

    auto ratio1 = TGraph(7, p9065_d30x1y1_xval, p9065_d30x1y1_yval);
    auto ratio2 = TGraph(7, p9065_d30x1y2_xval, p9065_d30x1y2_yval);
    auto ratio3 = TGraph(7, p9065_d30x1y3_xval, p9065_d30x1y3_yval);

    auto canvas2 = TCanvas("c2", "c2", 50, 50, 800, 600);
    canvas2.cd();
   
    ratio1.SetTitle("#phi(1020) Ratios to Stable Particles in p-Pb Collisions");
    ratio1.GetXaxis()->SetTitle("Multiplicity");
    ratio1.GetYaxis()->SetTitle("(Particle Ratio)/(0-5% Ratio)");
    ratio1.GetYaxis()->SetRangeUser(0.5, 1.5);
    ratio1.SetMarkerStyle(20);
    ratio1.SetMarkerSize(2);
    ratio1.SetMarkerColor(kRed+1);
    ratio1.SetLineColor(kRed+1);
    ratio1.Draw("AP");


    ratio2.SetMarkerStyle(21);
    ratio2.SetMarkerSize(2);
    ratio2.SetMarkerColor(kBlue+2);
    ratio2.SetLineColor(kBlue+2);
    ratio2.Draw("P");


    ratio3.SetMarkerStyle(23);
    ratio3.SetMarkerSize(2);
    ratio3.SetMarkerColor(kGreen+2);
    ratio3.SetLineColor(kGreen+2);
    ratio3.Draw("P");

    legend->Draw("SAME");


}
