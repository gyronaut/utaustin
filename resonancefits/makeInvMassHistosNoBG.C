#include <string>
#include <TH1F> 

// Relativisitic Breit-Wigner from Subash (and accompanying functions)
// To call, use the TF1 with FitFunRelBW
Double_t s(Double_t x0, Double_t x1, Double_t x2)
{
    return pow(x0*x0-x1*x1-x2*x2,2.0)-4.*x1*x1*x2*x2;
}

Double_t PS(Double_t m, Double_t pT, Double_t T)
{
    Double_t mT = sqrt(m*m+pT*pT);
    return m/mT*exp(-mT/T);
}

Double_t bw(Double_t m, Double_t m0, Double_t Gamma)
{
    return m*m*Gamma/(pow(m*m-m0*m0,2.0)+m*m*Gamma*Gamma);
    //return m*m*Gamma/(pow(m*m-m0*m0,2.0)+m*m*m*m*Gamma*Gamma/(m0*m0));
}

Double_t bw1(Double_t *x, Double_t *par)
{
    const Double_t MassK = 0.49368;
    const Double_t MassPi = 0.13957;
    const Double_t MassKstar = 0.892;
    const Double_t Gamma0 = 0.042;
    Double_t Gamma = Gamma0*pow(MassKstar/x[0],5.0);
    Gamma *= pow(s(x[0],MassK,MassPi)/s(MassKstar,MassK,MassPi),1.5);
    Double_t Gamma += par[2];
    //Double_t Gamma = par[2];
    return bw(x[0],par[1],Gamma);
}

Double_t KstarFun(Double_t *x, Double_t *par)
{
    const Double_t Temp = 0.160;
    //return par[0]*bw1(x, par)*PS(x[0], par[3], Temp);
    return par[0]*bw1(x, par);
    //return par[0]*bw1(x, par)*PS(x[0], par[3], Temp)*1.e6;
    
}

Double_t FitFunRelBWGaus(Double_t *x, Double_t *par)
{
    return KstarFun(x,&par[3])+par[0]*TMath::Gaus(*x,par[1],par[2]);
}

Double_t FitFunRelBW(Double_t *x, Double_t *par)
{
    return KstarFun(x,par);
}

//My attempts at RBW, doesn't include pT dependent term...
Double_t relativisticBW(Double_t *x, Double_t *par){
    Double_t xx = x[0];
    Double_t S = par[0];
    Double_t M = par[1];
    Double_t G = par[2];
    Double_t g = TMath::Sqrt(M*M*(M*M + G*G));
    Double_t k = 2.0*TMath::Sqrt(2)*M*G*g/(TMath::Pi()*TMath::Sqrt(M**2 + g));
    Double_t rbw = S*k/((xx**2 - M**2)**2 + (M**2)*(G**2));
    return rbw;
}

Double_t rbwWithBG(Double_t *x, Double_t *par){
    Double_t xx = x[0];
    Double_t S = par[0];
    Double_t M = par[1];
    Double_t G = par[2];
    Double_t g = TMath::Sqrt(M*M*(M*M + G*G));
    Double_t k = 2.0*TMath::Sqrt(2)*M*G*g/(TMath::Pi()*TMath::Sqrt(M**2 + g));
    Double_t rbw = S*k/((xx**2 - M**2)**2 + (M**2)*(G**2));
    
    Double_t height = par[3];
    Double_t mean = par[4];
    Double_t sigma = par[5];

    Double_t bg = height*TMath::Gaus(xx, mean, sigma, kFALSE);
    return rbw+bg;
}


int makeInvMassHistosNoBG(){
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

    string folder = "/Users/jtblair/Downloads/invm_decayed/pt02/";
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



    TFile *output = new TFile("output_invm_norescatter_relBW_noPS_simplewidthvar_2017_03_03.root", "RECREATE");

    TH1D *kstar0mass = new TH1D("kstar0mass", "K*^{0} Mass vs. p_{T}", 20, 0.0, 4.0);
    TH1D *kstar0massBG = new TH1D("kstar0massBG", "K*^{0} Mass vs. p_{T} with BG", 20, 0.0, 4.0);
    TH1D *kstar0width = new TH1D("kstar0width", "K*^{0} Width vs. p_{T}", 20, 0.0, 4.0);
    TH1D *kstar0widthBG = new TH1D("kstar0widthBG", "K*^{0} Width vs. p_{T}", 20, 0.0, 4.0);
    kstar0mass->GetXaxis()->SetTitle("p_{T} (GeV/c)");
    kstar0mass->GetYaxis()->SetTitle("Mass (GeV/c^{2})");
    kstar0width->GetXaxis()->SetTitle("p_{T} (GeV/c)");
    kstar0width->GetYaxis()->SetTitle("Width (GeV/c^2)");

    double mass = 0.0, width = 0.0, massBG=0.0;
    double massError = 0.0, widthError = 0.0, massBGError=0.0;

    TCanvas *canvas[9];

    for(int nfile = 0; nfile < 20; nfile++){
        double meanPT = (double)(nfile*2+1)/10.0;
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
        

        printf("****** Fits for file: %s ******\n", filename.c_str());
        for(int i=3; i<4; i++){
           
            //histos[i]->Sumw2(); 
            if(nfile==0){
                canvas[i] = new TCanvas(Form("c%i", i),Form("c%i", i), 0,0,900,900);
                canvas[i]->Divide(5,4);
            }
            canvas[i]->cd(nfile+1);
            //canvas[i]->SetLogy();
            histos[i]->SetLineColor(1);
            histos[i]->SetLineWidth(1);
            histos[i]->GetXaxis()->SetRangeUser(0.7, 1.2);
            histos[i]->GetYaxis()->SetRangeUser(0, 1.5*histos[i]->GetBinContent(histos[i]->GetMaximumBin()));
            histos[i]->SetStats(kFALSE);
            //histos[i]->Draw("HIST");
            //secondFit->Draw("SAME");
 
/*
            if(nfile<5){
                //TF1 *fit = new TF1(Form("fitPTbin0%dparticle%d", nfile*2+1, i), "[0]*TMath::BreitWigner(x, [1], [2])", 0.6, 1.2);
                TF1 *fit = new TF1(Form("fitPTbin0%dparticle%d", nfile*2+1, i), relativisticBW, 0.6, 1.2, 3);
                //TF1 *secondFit = new TF1(Form("secondFitPTbin0%dparticle%d", nfile*2+1, i), "[0]*TMath::BreitWigner(x, [1], [2]) + ([3] + [4]*x + [5]*x*x)*0.5*(1 + TMath::Sign(1, [3] + [4]*x + [5]*x*x))", 0.6, 1.2);
                //TF1 *secondFit = new TF1(Form("secondFitPTbin0%dparticle%d", nfile*2+1, i), "[0]*TMath::BreitWigner(x, [1], [2]) + gaus(3)", 0.6, 1.2);
                TF1 *secondFit = new TF1(Form("secondFitPTbin0%dparticle%d", nfile*2+1, i), rbwWithBG, 0.6, 1.2, 6); 
                //TF1 *bgFit = new TF1(Form("bgFitPTbin0%dparticle%d", nfile*2+1, i), "[0] + [1]*x + [2]*x*x", 0.6, 0.8);
                TF1 *bgFit = new TF1(Form("bgFitPTbin0%dparticle%d", nfile*2+1, i), "gaus(0)", 0.6, 1.2);
                bg[i] = (TH1D*)histos[i]->Clone(Form("bgPTbin0%dparticle%d", nfile*2+1, i));              
            }else{
                //TF1 *fit = new TF1(Form("fitPTbin%dparticle%d", nfile*2+1, i), "[0]*TMath::BreitWigner(x, [1], [2])", 0.6, 1.2);
                TF1 *fit = new TF1(Form("fitPTbin%dparticle%d", nfile*2+1, i), relativisticBW, 0.6, 1.2, 3);
                //TF1 *secondFit = new TF1(Form("secondFitPTbin%dparticle%d", nfile*2+1, i), "[0]*TMath::BreitWigner(x, [1], [2]) + ([3] + [4]*x + [5]*x*x)*0.5*(1 + TMath::Sign(1, [3] + [4]*x + [5]*x*x))", 0.6, 1.2);
                //TF1 *secondFit = new TF1(Form("secondFitPTbin%dparticle%d", nfile*2+1, i), "[0]*TMath::BreitWigner(x, [1], [2]) + gaus(3)", 0.6, 1.2);
                TF1 *secondFit = new TF1(Form("secondFitPTbin%dparticle%d", nfile*2+1, i), rbwWithBG, 0.6, 1.2, 6);               
                //TF1 *bgFit = new TF1(Form("bgFitPTbin%dparticle%d", nfile*2+1, i), "[0] + [1]*x + [2]*x*x", 0.6, 0.8);
                TF1 *bgFit = new TF1(Form("bgFitPTbin%dparticle%d", nfile*2+1, i), "gaus(0)", 0.6, 1.2);
                bg[i] = (TH1D*)histos[i]->Clone(Form("bgPTbin%dparticle%d", nfile*2+1, i));              
            }
*/
            //bg[i] = (TH1D*)histos[i]->Clone(Form("bgPTbin%d00particle%d", nfile*2+1, i));

            printf("mean PT: %f", meanPT);

            TF1 *fit = new TF1(Form("fitPTbin%d00particle%d", nfile*2+1, i), FitFunRelBW, 0.7, 1.1, 4);
            //TF1 *bgFit = new TF1(Form("bgFitPTbin%d00particle%d", nfile*2+1, i), "gaus(0)", 0.6, 0.95);
            //TF1 *secondFit = new TF1(Form("secondFitPTbin%d00particle%d", nfile*2+1, i), FitFunRelBWGaus, 0.75, 1.1, 7);

            fit->SetParNames("BW Area", "Mass", "Width", "PT");
            fit->SetParameters(1.0, 0.89, 0.0474, 0.5);
            fit->SetParLimits(0, .00001, 1.e3);
            fit->SetParLimits(1, 0.80, 1.0);
            fit->SetParLimits(2, 0.02, 0.1);
            fit->FixParameter(3, meanPT);
            fit->SetLineColor(2);
/*
            bgFit->SetParameter(0, 10);
            bgFit->SetParLimits(0, 0.00001, 100);
            bgFit->SetParameter(1, 0.8);
            bgFit->SetParLimits(1, 0.7, 1.1);
            bgFit->SetParameter(2, 1.0);
            bgFit->SetParLimits(2, 0.08, 1000.0);
*/

            printf("%s\n", fit->GetName());
            histos[i]->Fit(Form("fitPTbin%d00particle%d", nfile*2+1, i), "R", "SAME");
/*
            for(int j=0; j<100; j++){
                bg[i]->SetBinContent(j, bg[i]->GetBinContent(j) - fit->Eval(bg[i]->GetBinCenter(j)));
            } 

            bg[i]->Fit(Form("bgFitPTbin%d00particle%d", nfile*2+1, i), "R N", "SAME");

            secondFit->SetParNames("Gaus Scale", "Gaus Mean", "Gaus Width", "BW Area", "Mass", "Width", "PT");
            secondFit->SetParameters(bgFit->GetParameter(0), bgFit->GetParameter(1), bgFit->GetParameter(2), fit->GetParameter(0), fit->GetParameter(1), fit->GetParameter(2), fit->GetParameter(3));
            secondFit->SetParLimits(3, 0.00001, 1.e3);
            secondFit->SetParLimits(4, 0.87, 0.92);
            //secondFit->SetParLimits(5, 0.035, 0.065);
            secondFit->FixParameter(5, 0.0474);
            secondFit->FixParameter(6, fit->GetParameter(3));
            secondFit->SetParLimits(0, 0.0001, bgFit->GetParameter(0)*1.2);
            secondFit->SetParLimits(1, 0.7, 1.1);
            secondFit->SetParLimits(2, bgFit->GetParameter(2)*0.8, bgFit->GetParameter(2)*1.2);
            secondFit->SetLineColor(2);
            secondFit->SetLineStyle(2);

            printf("set parameters for second fit!\n\n");

            histos[i]->Fit(Form("secondFitPTbin%d00particle%d", nfile*2+1, i), "R", "SAME");
            */
            histos[i]->Draw("HIST SAME");
            fit->Draw("SAME");
            //secondFit->Draw("SAME");

            printf("\n");
            histos[i]->Write();
            //bg[i]->Write();
            fit->Write();
            //secondFit->Write();
            //bgFit->Write();
            //Do mass and width vs. pT plots just for K*0
            if(i==3){
                mass = fit->GetParameter(1);
                //massBG = secondFit->GetParameter(4);
                massError = fit->GetParError(1);
                //massBGError = secondFit->GetParError(4);
                width = fit->GetParameter(2);
                widthError = fit->GetParError(2);
                //widthBG = secondFit->GetParameter(5);
                //widthBGError = secondFit->GetParError(5);

                kstar0mass->SetBinContent(nfile+1, mass);
                kstar0mass->SetBinError(nfile+1, massError);

                //kstar0massBG->SetBinContent(nfile+1, massBG);
                //kstar0massBG->SetBinError(nfile+1, massBGError);

                kstar0width->SetBinContent(nfile+1, width);
                kstar0width->SetBinError(nfile+1, widthError);

                //kstar0widthBG->SetBinContent(nfile+1, widthBG);
                //kstar0widthBG->SetBinError(nfile+1, widthBGError);
                
                if(nfile==4){
                    TCanvas *singlecanvas = new TCanvas("singlecanvas", "singlecanvas", 0,0,600,600);
                    singlecanvas->cd();
                    printf("Got here! \n");
                    histos[i]->Draw("HIST SAME");
                    /*
                    bgFit->SetParameter(0, secondFit->GetParameter(0));
                    bgFit->SetParameter(1, secondFit->GetParameter(1));
                    bgFit->SetParameter(2, secondFit->GetParameter(2));
                    bgFit->SetLineColor(9);
                    bgFit->SetLineStyle(4);

                    fit->SetParameter(0, secondFit->GetParameter(3));
                    fit->SetParameter(1, secondFit->GetParameter(4));
                    fit->SetParameter(2, secondFit->GetParameter(5));
                    */
                    fit->SetLineColor(8);
                    fit->SetLineStyle(1);
                    //bgFit->Draw("SAME");
                    
                    fit->Draw("SAME");
                    //secondFit->Draw("SAME");
                }
            }
        }
        printf("************************************************************\n");
         
    }

    for(int i=3; i<4; i++){
        canvas[i]->Write();
    }
    kstar0mass->Write();
    //kstar0massBG->Write();
    kstar0width->Write();
    //kstar0widthBG->Write();
    //output->Close();
}
