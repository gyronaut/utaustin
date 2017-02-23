#include <string>
#include <TH1F> 

int makeInvMassHistos(){
    int NUM_MASS_BINS = 100;
    double MASS_LOW = 0.0;
    double MASS_HIGH = 2.0;
    string particles [8];
    particles[0] = "K*^{+} + K*^{0}";
    particles[1] = "K*^{-} + #bar{K}*^{0}";
    particles[2] = "K*^{+}";
    particles[3] = "K*^{0}";
    particles[4] = "K*^{-}";
    particles[5] = "#bar{K}*^{0}";
    particles[6] = "K*^{0} + #bar{K}*^{0}";
    particles[7] = "K*^{+} + K*^{-}";

    string folder = "/Users/jtblair/Downloads/invm/pt02/";
    string files[20];
    files[0] = "invm_[0.0,0.2].dat";
    files[1] = "invm_[0.2,0.4].dat";
    files[2] = "invm_[0.4,0.6].dat";
    files[3] = "invm_[0.6,0.8].dat";
    files[4] = "invm_[0.8,1.0].dat";
    files[5] = "invm_[1.0,1.2].dat";
    files[6] = "invm_[1.2,1.4].dat";
    files[7] = "invm_[1.4,1.6].dat";   
    files[8] = "invm_[1.6,1.8].dat";
    files[9] = "invm_[1.8,2.0].dat";
    files[10] = "invm_[2.0,2.2].dat";
    files[11] = "invm_[2.2,2.4].dat";
    files[12] = "invm_[2.4,2.6].dat";
    files[13] = "invm_[2.6,2.8].dat";
    files[14] = "invm_[2.8,3.0].dat";
    files[15] = "invm_[3.0,3.2].dat";
    files[16] = "invm_[3.2,3.4].dat";
    files[17] = "invm_[3.4,3.6].dat";
    files[18] = "invm_[3.6,3.8].dat";
    files[19] = "invm_[3.8,4.0].dat";



    TFile *output = new TFile("output_invm_bggaus_23_02_2017.root", "RECREATE");

    TH1D *kstar0mass = new TH1D("kstar0mass", "K*^{0} Mass vs. p_{T}", 20, 0.0, 4.0);
    TH1D *kstar0massBG = new TH1D("kstar0massBG", "K*^{0} Mass vs. p_{T} with BG", 20, 0.0, 4.0);
    TH1D *kstar0width = new TH1D("kstar0width", "K*^{0} Width vs. p_{T}", 20, 0.0, 4.0);
    kstar0mass->GetXaxis()->SetTitle("p_{T} (GeV/c)");
    kstar0mass->GetYaxis()->SetTitle("Mass (GeV/c^{2})");
    kstar0width->GetXaxis()->SetTitle("p_{T} (GeV/c)");
    kstar0width->GetYaxis()->SetTitle("Width (GeV/c^2)");

    double mass = 0.0, width = 0.0, massBG=0.0;
    double massError = 0.0, widthError = 0.0, massBGError=0.0;

    TCanvas *canvas[9];

    for(int nfile = 0; nfile < 20; nfile++){
        string filename = folder+files[nfile];
        string ptLower = filename.substr(filename.find("[")+1, 3);
        string ptHigher = filename.substr(filename.find(",")+1, 3);   
        TH1D* histos[8];
        TH1D* bg[8];
        for(int i=0; i<8; i++){
            if(nfile<5){
                histos[i] = new TH1D(Form("ptbin0%dparticle%d",nfile*2+1, i), Form("Invariant Mass for (%s), %s < p_{T} < %s",particles[i].c_str(), ptLower.c_str(), ptHigher.c_str()), 100, 0.0, 2.0);
            }else{
                histos[i] = new TH1D(Form("ptbin%dparticle%d",nfile*2+1, i), Form("Invariant Mass for (%s), %s < p_{T} < %s",particles[i].c_str(), ptLower.c_str(), ptHigher.c_str()), 100, 0.0, 2.0);
            }
            histos[i]->GetXaxis()->SetTitle("Invariant Mass (GeV/c^{2})");
            histos[i]->GetYaxis()->SetTitle("Counts");
        }

        ifstream input;
        input.open(filename.c_str());
        string line = "";
        if(input.good()){
            getline(input, line);
        }

        double massBin=0.0;
        double invMass[8];
        for(int i=0; i<8; i++){
            invMass[i] = 0.0;
        }
        int lineNumber = 1;
        while(1){
            input >> massBin >> invMass[0] >> invMass[1] >> invMass[2] >> invMass[3] >> invMass[4] >> invMass[5] >> invMass[6] >> invMass[7];
            if(!input.good())break;
            for(int i =0; i<8; i++){
                histos[i]->SetBinContent(lineNumber, invMass[i]);
            }
            lineNumber++;
        }
        
        for(int i =0; i<8; i++){
            if(nfile==0){
                canvas[i] = new TCanvas(Form("c%i", i),Form("c%i", i), 0,0,900,900);
                canvas[i]->Divide(5,4);
            }
            canvas[i]->cd(nfile+1);
            //canvas[i]->SetLogy();
            histos[i]->Sumw2();
            histos[i]->SetLineColor(1);
            histos[i]->SetLineWidth(1);
            histos[i]->GetXaxis()->SetRangeUser(0.6, 1.2);
            //histos[i]->GetYaxis()->SetRangeUser(0, 100);
            histos[i]->Draw("H SAME");
        }
 
        printf("****** Fits for file: %s ******\n", filename.c_str());
        for(int i=0; i<8; i++){
            if(nfile<5){
                TF1 *fit = new TF1(Form("fitPTbin0%dparticle%d", nfile*2+1, i), "[0]*TMath::BreitWigner(x, [1], [2])", 0.6, 1.2);
                //TF1 *secondFit = new TF1(Form("secondFitPTbin0%dparticle%d", nfile*2+1, i), "[0]*TMath::BreitWigner(x, [1], [2]) + ([3] + [4]*x + [5]*x*x)*0.5*(1 + TMath::Sign(1, [3] + [4]*x + [5]*x*x))", 0.6, 1.2);
                //TF1 *bgFit = new TF1(Form("bgFitPTbin0%dparticle%d", nfile*2+1, i), "[0] + [1]*x + [2]*x*x", 0.6, 0.8);
                TF1 *secondFit = new TF1(Form("secondFitPTbin0%dparticle%d", nfile*2+1, i), "[0]*TMath::BreitWigner(x, [1], [2]) + gaus(3)", 0.6, 1.2);
                TF1 *bgFit = new TF1(Form("bgFitPTbin0%dparticle%d", nfile*2+1, i), "gaus(0)", 0.6, 1.2);
                bg[i] = (TH1D*)histos[i]->Clone(Form("bgPTbin0%dparticle%d", nfile*2+1, i));              
           }else{
                TF1 *fit = new TF1(Form("fitPTbin%dparticle%d", nfile*2+1, i), "[0]*TMath::BreitWigner(x, [1], [2])", 0.6, 1.2);
                //TF1 *secondFit = new TF1(Form("secondFitPTbin%dparticle%d", nfile*2+1, i), "[0]*TMath::BreitWigner(x, [1], [2]) + ([3] + [4]*x + [5]*x*x)*0.5*(1 + TMath::Sign(1, [3] + [4]*x + [5]*x*x))", 0.6, 1.2);
                //TF1 *bgFit = new TF1(Form("bgFitPTbin%dparticle%d", nfile*2+1, i), "[0] + [1]*x + [2]*x*x", 0.6, 0.8);
                TF1 *secondFit = new TF1(Form("secondFitPTbin%dparticle%d", nfile*2+1, i), "[0]*TMath::BreitWigner(x, [1], [2]) + gaus(3)", 0.6, 1.2);
                TF1 *bgFit = new TF1(Form("bgFitPTbin%dparticle%d", nfile*2+1, i), "gaus(0)", 0.6, 1.2);
                bg[i] = (TH1D*)histos[i]->Clone(Form("bgPTbin%dparticle%d", nfile*2+1, i));              
          }
            fit->SetParameter(0, 10.0);
            fit->SetParameter(1, 0.9);
            fit->SetParameter(2, 0.0474);
            //fit->FixParameter(2, 0.0474);
            fit->SetParLimits(0, 0.1, 1000.0);
            fit->SetParLimits(1, 0.6, 1.2);
            fit->SetParLimits(2, 0.001, 0.2);
            fit->SetParNames("scale", "M", "Gamma");
            fit->SetLineColor(2);

            bgFit->SetParameter(0, 10);
            bgFit->SetParLimits(0, 0.01, 30);
            bgFit->SetParameter(1, 0.8);
            bgFit->SetParLimits(1, 0.75, 1.2);
            bgFit->SetParameter(2,1.0);
            bgFit->SetParLimits(2, 0.05, 4.0);


            printf("%s\n", fit->GetName());
            if(nfile<5){
                histos[i]->Fit(Form("fitPTbin0%dparticle%d", nfile*2+1, i), "R", "SAME");
                //histos[i]->Fit(Form("bgFitPTbin0%dpatriclt%d", nfile*2+1, i), "R");
            }else{
                histos[i]->Fit(Form("fitPTbin%dparticle%d", nfile*2+1, i), "R", "SAME");
                //histos[i]->Fit(Form("bgFitPTbin%dparticle%d", nfile*2+1, i), "R");
            }
           
            for(int j=0; j<100; j++){
                bg[i]->SetBinContent(j, bg[i]->GetBinContent(j) - fit->Eval(bg[i]->GetBinCenter(j)));
            } 

            if(nfile<5){
                bg[i]->Fit(Form("bgFitPTbin0%dparticle%d", nfile*2+1, i), "R", "SAME");
            }else{
                bg[i]->Fit(Form("bgFitPTbin%dparticle%d", nfile*2+1, i), "R", "SAME");
            }
 
            secondFit->SetParameter(0, fit->GetParameter(0));
            secondFit->SetParameter(1, fit->GetParameter(0));
            secondFit->SetParameter(2, fit->GetParameter(0));
            secondFit->SetParameter(3, bgFit->GetParameter(0));
            secondFit->SetParLimits(3, bgFit->GetParameter(0)*0.8, bgFit->GetParameter(0)*1.2);
            secondFit->SetParameter(4, bgFit->GetParameter(1));
            secondFit->SetParLimits(4, bgFit->GetParameter(1)*0.8, bgFit->GetParameter(1)*1.2);
            secondFit->SetParameter(5, bgFit->GetParameter(2));
            secondFit->SetParLimits(5, bgFit->GetParameter(2)*0.8, bgFit->GetParameter(2)*1.2);
            secondFit->SetParLimits(0, 0.1, 1000.0);
            secondFit->SetParLimits(1, 0.6, 1.2);
            secondFit->SetParLimits(2, 0.001, 0.2);
            secondFit->SetParNames("scale", "M", "Gamma", "c", "b", "a");
            secondFit->SetLineColor(2);
            secondFit->SetLineStyle(2);

            if(nfile<5){
                histos[i]->Fit(Form("secondFitPTbin0%dparticle%d", nfile*2+1, i), "R", "SAME");
            }else{
                histos[i]->Fit(Form("secondFitPTbin%dparticle%d", nfile*2+1, i), "R", "SAME");
            }
        

            printf("\n");
            histos[i]->Write();
            bg[i]->Write();
            fit->Write();
            secondFit->Write();
            bgFit->Write();
            //Do mass and width vs. pT plots just for K*0
            if(i==3){
                mass = fit->GetParameter(1);
                massBG = secondFit->GetParameter(1);
                massError = fit->GetParError(1);
                massBGError = secondFit->GetParError(1);
                width = fit->GetParameter(2);
                widthError = fit->GetParError(2);

                kstar0mass->SetBinContent(nfile+1, mass);
                kstar0mass->SetBinError(nfile+1, massError);

                kstar0massBG->SetBinContent(nfile+1, massBG);
                kstar0massBG->SetBinError(nfile+1, massBGError);

                kstar0width->SetBinContent(nfile+1, width);
                kstar0width->SetBinError(nfile+1, widthError);
            }
        }
        printf("************************************************************\n");
         
    }

    for(int i=0; i<8; i++){
        canvas[i]->Write();
    }
    kstar0mass->Write();
    kstar0massBG->Write();
    kstar0width->Write();
    //output->Close();
}
