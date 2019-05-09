void plotFitSyst(){
    Float_t stdtotal0_20 = 5.850084E-02;
    Float_t stdtotal20_50 = 3.489754E-02;
    Float_t stdtotal50_100 = 1.881542E-02;
    Float_t stdnear0_20 = 1.509740E-03;
    Float_t stdnear20_50 = 1.450373E-03;
    Float_t stdnear50_100 = 8.015203E-04;
    Float_t stdaway0_20 = 2.234562E-03;
    Float_t stdaway20_50 = 2.274811E-03;
    Float_t stdaway50_100 = 1.119789E-03;
    Float_t stdbulk0_20 = 5.476692E-02;
    Float_t stdbulk20_50 = 3.111197E-02;
    Float_t stdbulk50_100 = 1.686431E-02;


    TFile* peakfile = new TFile("masstest/fitsyst_allpeak.root");
    TH1F* peaktotal0_20 = (TH1F*)(peakfile->Get("hphiyieldSyst_0_20")->Clone("peaktotal0_20"));
    TH1F* peaktotal20_50 = (TH1F*)(peakfile->Get("hphiyieldSyst_20_50")->Clone("peaktotal20_50"));
    TH1F* peaktotal50_100 = (TH1F*)(peakfile->Get("hphiyieldSyst_50_100")->Clone("peaktotal50_100"));

    TH1F* peaknear0_20 = (TH1F*)(peakfile->Get("hphinearyieldSyst_0_20")->Clone("peaknear0_20"));
    TH1F* peaknear20_50 = (TH1F*)(peakfile->Get("hphinearyieldSyst_20_50")->Clone("peaknear20_50"));
    TH1F* peaknear50_100 = (TH1F*)(peakfile->Get("hphinearyieldSyst_50_100")->Clone("peaknear50_100"));
 
    TH1F* peakaway0_20 = (TH1F*)(peakfile->Get("hphiawayyieldSyst_0_20")->Clone("peakaway0_20"));
    TH1F* peakaway20_50 = (TH1F*)(peakfile->Get("hphiawayyieldSyst_20_50")->Clone("peakaway20_50"));
    TH1F* peakaway50_100 = (TH1F*)(peakfile->Get("hphiawayyieldSyst_50_100")->Clone("peakaway50_100"));

    TH1F* peakbulk0_20 = (TH1F*)(peakfile->Get("hphibulkyieldSyst_0_20")->Clone("peakbulk0_20"));
    TH1F* peakbulk20_50 = (TH1F*)(peakfile->Get("hphibulkyieldSyst_20_50")->Clone("peakbulk20_50"));
    TH1F* peakbulk50_100 = (TH1F*)(peakfile->Get("hphibulkyieldSyst_50_100")->Clone("peakbulk50_100"));

    //LSB
    TFile* LSBfile = new TFile("masstest/fitsyst_allLSB.root");
    TH1F* LSBtotal0_20 = (TH1F*)(LSBfile->Get("hphiyieldSyst_0_20")->Clone("LSBtotal0_20"));
    TH1F* LSBtotal20_50 = (TH1F*)(LSBfile->Get("hphiyieldSyst_20_50")->Clone("LSBtotal20_50"));
    TH1F* LSBtotal50_100 = (TH1F*)(LSBfile->Get("hphiyieldSyst_50_100")->Clone("LSBtotal50_100"));

    TH1F* LSBnear0_20 = (TH1F*)(LSBfile->Get("hphinearyieldSyst_0_20")->Clone("LSBnear0_20"));
    TH1F* LSBnear20_50 = (TH1F*)(LSBfile->Get("hphinearyieldSyst_20_50")->Clone("LSBnear20_50"));
    TH1F* LSBnear50_100 = (TH1F*)(LSBfile->Get("hphinearyieldSyst_50_100")->Clone("LSBnear50_100"));
 
    TH1F* LSBaway0_20 = (TH1F*)(LSBfile->Get("hphiawayyieldSyst_0_20")->Clone("LSBaway0_20"));
    TH1F* LSBaway20_50 = (TH1F*)(LSBfile->Get("hphiawayyieldSyst_20_50")->Clone("LSBaway20_50"));
    TH1F* LSBaway50_100 = (TH1F*)(LSBfile->Get("hphiawayyieldSyst_50_100")->Clone("LSBaway50_100"));

    TH1F* LSBbulk0_20 = (TH1F*)(LSBfile->Get("hphibulkyieldSyst_0_20")->Clone("LSBbulk0_20"));
    TH1F* LSBbulk20_50 = (TH1F*)(LSBfile->Get("hphibulkyieldSyst_20_50")->Clone("LSBbulk20_50"));
    TH1F* LSBbulk50_100 = (TH1F*)(LSBfile->Get("hphibulkyieldSyst_50_100")->Clone("LSBbulk50_100"));

    //RSB
    TFile* RSBfile = new TFile("masstest/fitsyst_allRSB.root");
    TH1F* RSBtotal0_20 = (TH1F*)(RSBfile->Get("hphiyieldSyst_0_20")->Clone("RSBtotal0_20"));
    TH1F* RSBtotal20_50 = (TH1F*)(RSBfile->Get("hphiyieldSyst_20_50")->Clone("RSBtotal20_50"));
    TH1F* RSBtotal50_100 = (TH1F*)(RSBfile->Get("hphiyieldSyst_50_100")->Clone("RSBtotal50_100"));

    TH1F* RSBnear0_20 = (TH1F*)(RSBfile->Get("hphinearyieldSyst_0_20")->Clone("RSBnear0_20"));
    TH1F* RSBnear20_50 = (TH1F*)(RSBfile->Get("hphinearyieldSyst_20_50")->Clone("RSBnear20_50"));
    TH1F* RSBnear50_100 = (TH1F*)(RSBfile->Get("hphinearyieldSyst_50_100")->Clone("RSBnear50_100"));
 
    TH1F* RSBaway0_20 = (TH1F*)(RSBfile->Get("hphiawayyieldSyst_0_20")->Clone("RSBaway0_20"));
    TH1F* RSBaway20_50 = (TH1F*)(RSBfile->Get("hphiawayyieldSyst_20_50")->Clone("RSBaway20_50"));
    TH1F* RSBaway50_100 = (TH1F*)(RSBfile->Get("hphiawayyieldSyst_50_100")->Clone("RSBaway50_100"));

    TH1F* RSBbulk0_20 = (TH1F*)(RSBfile->Get("hphibulkyieldSyst_0_20")->Clone("RSBbulk0_20"));
    TH1F* RSBbulk20_50 = (TH1F*)(RSBfile->Get("hphibulkyieldSyst_20_50")->Clone("RSBbulk20_50"));
    TH1F* RSBbulk50_100 = (TH1F*)(RSBfile->Get("hphibulkyieldSyst_50_100")->Clone("RSBbulk50_100"));

    //scale
    TFile* scalefile = new TFile("masstest/fitsyst_allscale.root");
    TH1F* scaletotal0_20 = (TH1F*)(scalefile->Get("hphiyieldSyst_0_20")->Clone("scaletotal0_20"));
    TH1F* scaletotal20_50 = (TH1F*)(scalefile->Get("hphiyieldSyst_20_50")->Clone("scaletotal20_50"));
    TH1F* scaletotal50_100 = (TH1F*)(scalefile->Get("hphiyieldSyst_50_100")->Clone("scaletotal50_100"));

    TH1F* scalenear0_20 = (TH1F*)(scalefile->Get("hphinearyieldSyst_0_20")->Clone("scalenear0_20"));
    TH1F* scalenear20_50 = (TH1F*)(scalefile->Get("hphinearyieldSyst_20_50")->Clone("scalenear20_50"));
    TH1F* scalenear50_100 = (TH1F*)(scalefile->Get("hphinearyieldSyst_50_100")->Clone("scalenear50_100"));
 
    TH1F* scaleaway0_20 = (TH1F*)(scalefile->Get("hphiawayyieldSyst_0_20")->Clone("scaleaway0_20"));
    TH1F* scaleaway20_50 = (TH1F*)(scalefile->Get("hphiawayyieldSyst_20_50")->Clone("scaleaway20_50"));
    TH1F* scaleaway50_100 = (TH1F*)(scalefile->Get("hphiawayyieldSyst_50_100")->Clone("scaleaway50_100"));

    TH1F* scalebulk0_20 = (TH1F*)(scalefile->Get("hphibulkyieldSyst_0_20")->Clone("scalebulk0_20"));
    TH1F* scalebulk20_50 = (TH1F*)(scalefile->Get("hphibulkyieldSyst_20_50")->Clone("scalebulk20_50"));
    TH1F* scalebulk50_100 = (TH1F*)(scalefile->Get("hphibulkyieldSyst_50_100")->Clone("scalebulk50_100"));
    
    //fits
    TFile* fitsfile = new TFile("fitsyst_allfits.root");
    TH1F* fitstotal0_20 = (TH1F*)(fitsfile->Get("hphiyieldSyst_0_20")->Clone("fitstotal0_20"));
    TH1F* fitstotal20_50 = (TH1F*)(fitsfile->Get("hphiyieldSyst_20_50")->Clone("fitstotal20_50"));
    TH1F* fitstotal50_100 = (TH1F*)(fitsfile->Get("hphiyieldSyst_50_100")->Clone("fitstotal50_100"));

    TH1F* fitsnear0_20 = (TH1F*)(fitsfile->Get("hphinearyieldSyst_0_20")->Clone("fitsnear0_20"));
    TH1F* fitsnear20_50 = (TH1F*)(fitsfile->Get("hphinearyieldSyst_20_50")->Clone("fitsnear20_50"));
    TH1F* fitsnear50_100 = (TH1F*)(fitsfile->Get("hphinearyieldSyst_50_100")->Clone("fitsnear50_100"));
 
    TH1F* fitsaway0_20 = (TH1F*)(fitsfile->Get("hphiawayyieldSyst_0_20")->Clone("fitsaway0_20"));
    TH1F* fitsaway20_50 = (TH1F*)(fitsfile->Get("hphiawayyieldSyst_20_50")->Clone("fitsaway20_50"));
    TH1F* fitsaway50_100 = (TH1F*)(fitsfile->Get("hphiawayyieldSyst_50_100")->Clone("fitsaway50_100"));

    TH1F* fitsbulk0_20 = (TH1F*)(fitsfile->Get("hphibulkyieldSyst_0_20")->Clone("fitsbulk0_20"));
    TH1F* fitsbulk20_50 = (TH1F*)(fitsfile->Get("hphibulkyieldSyst_20_50")->Clone("fitsbulk20_50"));
    TH1F* fitsbulk50_100 = (TH1F*)(fitsfile->Get("hphibulkyieldSyst_50_100")->Clone("fitsbulk50_100")); 


    printf("          0-20%% | 20-50%% | 50-80%% |\n");
    printf("peaknear: %4.2f %4.2f %4.2f \n", peaknear0_20->GetRMS()*100.0, peaknear20_50->GetRMS()*100.0, peaknear50_100->GetRMS()*100.0);
    printf("peakaway: %4.2f %4.2f %4.2f \n", peakaway0_20->GetRMS()*100.0, peakaway20_50->GetRMS()*100.0, peakaway50_100->GetRMS()*100.0);
    printf("peakbulk: %4.2f %4.2f %4.2f \n", peakbulk0_20->GetRMS()*100.0, peakbulk20_50->GetRMS()*100.0, peakbulk50_100->GetRMS()*100.0);
    printf("peaktotal: %4.2f %4.2f %4.2f \n\n", peaktotal0_20->GetRMS()*100.0, peaktotal20_50->GetRMS()*100.0, peaktotal50_100->GetRMS()*100.0);

    printf("LSBnear: %4.2f %4.2f %4.2f \n", LSBnear0_20->GetRMS()*100.0, LSBnear20_50->GetRMS()*100.0, LSBnear50_100->GetRMS()*100.0);
    printf("LSBaway: %4.2f %4.2f %4.2f \n", LSBaway0_20->GetRMS()*100.0, LSBaway20_50->GetRMS()*100.0, LSBaway50_100->GetRMS()*100.0);
    printf("LSBbulk: %4.2f %4.2f %4.2f \n", LSBbulk0_20->GetRMS()*100.0, LSBbulk20_50->GetRMS()*100.0, LSBbulk50_100->GetRMS()*100.0);
    printf("LSBtotal: %4.2f %4.2f %4.2f \n\n", LSBtotal0_20->GetRMS()*100.0, LSBtotal20_50->GetRMS()*100.0, LSBtotal50_100->GetRMS()*100.0);

    printf("RSBnear: %4.2f %4.2f %4.2f \n", RSBnear0_20->GetRMS()*100.0, RSBnear20_50->GetRMS()*100.0, RSBnear50_100->GetRMS()*100.0);
    printf("RSBaway: %4.2f %4.2f %4.2f \n", RSBaway0_20->GetRMS()*100.0, RSBaway20_50->GetRMS()*100.0, RSBaway50_100->GetRMS()*100.0);
    printf("RSBbulk: %4.2f %4.2f %4.2f \n", RSBbulk0_20->GetRMS()*100.0, RSBbulk20_50->GetRMS()*100.0, RSBbulk50_100->GetRMS()*100.0);
    printf("RSBtotal: %4.2f %4.2f %4.2f \n\n", RSBtotal0_20->GetRMS()*100.0, RSBtotal20_50->GetRMS()*100.0, RSBtotal50_100->GetRMS()*100.0);

    printf("scalenear: %4.2f %4.2f %4.2f \n", scalenear0_20->GetRMS()*100.0, scalenear20_50->GetRMS()*100.0, scalenear50_100->GetRMS()*100.0);
    printf("scaleaway: %4.2f %4.2f %4.2f \n", scaleaway0_20->GetRMS()*100.0, scaleaway20_50->GetRMS()*100.0, scaleaway50_100->GetRMS()*100.0);
    printf("scalebulk: %4.2f %4.2f %4.2f \n", scalebulk0_20->GetRMS()*100.0, scalebulk20_50->GetRMS()*100.0, scalebulk50_100->GetRMS()*100.0);
    printf("scaletotal: %4.2f %4.2f %4.2f \n\n", scaletotal0_20->GetRMS()*100.0, scaletotal20_50->GetRMS()*100.0, scaletotal50_100->GetRMS()*100.0);

    printf("fitsnear: %4.2f %4.2f %4.2f \n", fitsnear0_20->GetRMS()*100.0, fitsnear20_50->GetRMS()*100.0, fitsnear50_100->GetRMS()*100.0);
    printf("fitsaway: %4.2f %4.2f %4.2f \n", fitsaway0_20->GetRMS()*100.0, fitsaway20_50->GetRMS()*100.0, fitsaway50_100->GetRMS()*100.0);
    printf("fitsbulk: %4.2f %4.2f %4.2f \n", fitsbulk0_20->GetRMS()*100.0, fitsbulk20_50->GetRMS()*100.0, fitsbulk50_100->GetRMS()*100.0);
    printf("fitstotal: %4.2f %4.2f %4.2f \n\n", fitstotal0_20->GetRMS()*100.0, fitstotal20_50->GetRMS()*100.0, fitstotal50_100->GetRMS()*100.0);


}

