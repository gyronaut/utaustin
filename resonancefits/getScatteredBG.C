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
    return m*m0*Gamma/(pow(m*m-m0*m0,2.0)+m0*m0*Gamma*Gamma);
    //return m*m*Gamma/(pow(m*m-m0*m0,2.0)+m*m*m*m*Gamma*Gamma/(m0*m0));
}

Double_t bw1(Double_t *x, Double_t *par)
{
    const Double_t MassK = 0.49368;
    const Double_t MassPi = 0.13957;
    //Double_t Gamma = par[2]*pow(par[1]/x[0],4.0);
    //Gamma *= pow(s(x[0],MassK,MassPi)/s(par[1],MassK,MassPi),1.5);
    Double_t Gamma = par[2];
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


int getScatteredBG(){
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

    string rescatterFolder = "/Users/jtblair/Downloads/invm/pt02/";
    string folder= "/Users/jtblair/Downloads/invm_decayed/pt02/";
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



    TFile *output = new TFile("onlyScatteredBG_20170426.root", "RECREATE");

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
        string rescatterFilename = rescatterFolder+files[nfile];
        string filename = folder+files[nfile];
        string ptLower = filename.substr(filename.find("[")+1, 3);
        string ptHigher = filename.substr(filename.find(",")+1, 3);   
        TH1D* histos[8];
        TH1D* decayHistos[8];
        TH1D* bgHistos[8];
        TH1D* bg[8];
        for(int i=0; i<8; i++){
            if(nfile<5){
                histos[i] = new TH1D(Form("ptbin0%dparticle%d",nfile*2+1, i), Form("Invariant Mass for (%s), %s < p_{T} < %s",particles[i].c_str(), ptLower.c_str(), ptHigher.c_str()), 100, 0.0, 2.0);
                decayHistos[i] = new TH1D(Form("DECAYptbin0%dparticle%d",nfile*2+1, i), Form("Invariant Mass for (%s), %s < p_{T} < %s",particles[i].c_str(), ptLower.c_str(), ptHigher.c_str()), 100, 0.0, 2.0);
                bgHistos[i] = new TH1D(Form("BGptbin0%dparticle%d",nfile*2+1, i), Form("Invariant Mass for (%s), %s < p_{T} < %s",particles[i].c_str(), ptLower.c_str(), ptHigher.c_str()), 100, 0.0, 2.0);

            }else{
                histos[i] = new TH1D(Form("ptbin%dparticle%d",nfile*2+1, i), Form("Invariant Mass for (%s), %s < p_{T} < %s",particles[i].c_str(), ptLower.c_str(), ptHigher.c_str()), 100, 0.0, 2.0);
                decayHistos[i] = new TH1D(Form("DECAYptbin%dparticle%d",nfile*2+1, i), Form("Invariant Mass for (%s), %s < p_{T} < %s",particles[i].c_str(), ptLower.c_str(), ptHigher.c_str()), 100, 0.0, 2.0);
                bgHistos[i] = new TH1D(Form("BGptbin%dparticle%d",nfile*2+1, i), Form("Invariant Mass for (%s), %s < p_{T} < %s",particles[i].c_str(), ptLower.c_str(), ptHigher.c_str()), 100, 0.0, 2.0);
            }
            histos[i]->GetXaxis()->SetTitle("Invariant Mass (GeV/c^{2})");
            histos[i]->GetYaxis()->SetTitle("Counts");
        }

        ifstream rescatterInput;
        rescatterInput.open(rescatterFilename.c_str());
        string line = "";
        if(rescatterInput.good()){
            getline(rescatterInput, line);
        }

        double massBin=0.0;
        double invMass[8];
        for(int i=0; i<8; i++){
            invMass[i] = 0.0;
        }
        int lineNumber = 1;
        while(1){
            rescatterInput >> massBin >> invMass[0] >> invMass[1] >> invMass[2] >> invMass[3] >> invMass[4] >> invMass[5] >> invMass[6] >> invMass[7];
            if(!rescatterInput.good())break;
            for(int i =0; i<8; i++){
                histos[i]->SetBinContent(lineNumber, invMass[i]);
            }
            lineNumber++;
        }
        
        ifstream input;
        input.open(filename.c_str());
        if(input.good()){
            getline(input, line);
        }
        lineNumber=1;
        while(1){
            input >> massBin >> invMass[0] >> invMass[1] >> invMass[2] >> invMass[3] >> invMass[4] >> invMass[5] >> invMass[6] >> invMass[7];
            if(!input.good())break;
            for(int i =0; i<8; i++){
                decayHistos[i]->SetBinContent(lineNumber, invMass[i]);
                bgHistos[i]->SetBinContent(lineNumber, histos[i]->GetBinContent(lineNumber) - invMass[i]);
            }
            lineNumber++;
        }


        printf("****** Fits for file: %s ******\n", filename.c_str());
        for(int i=3; i<4; i++){
           
            //histos[i]->Sumw2(); 
            //canvas[i]->SetLogy();
            histos[i]->SetLineColor(1);
            histos[i]->SetLineWidth(1);
            histos[i]->GetXaxis()->SetRangeUser(0.7, 1.2);
            histos[i]->GetYaxis()->SetRangeUser(0, 1.5*histos[i]->GetBinContent(histos[i]->GetMaximumBin()));
            histos[i]->SetStats(kFALSE);

            bgHistos[i]->SetLineColor(1);
            bgHistos[i]->SetLineWidth(1);
            bgHistos[i]->GetXaxis()->SetRangeUser(0.7, 1.2);
            //bgHistos[i]->GetYaxis()->SetRangeUser(0, 1.5*bgHistos[i]->GetBinContent(bgHistos[i]->GetMaximumBin()));
            bgHistos[i]->SetStats(kFALSE);
 
            decayHistos[i]->SetLineColor(1);
            decayHistos[i]->SetLineWidth(1);
            decayHistos[i]->GetXaxis()->SetRangeUser(0.7, 1.2);
            decayHistos[i]->GetYaxis()->SetRangeUser(0, 1.5*decayHistos[i]->GetBinContent(decayHistos[i]->GetMaximumBin()));
            decayHistos[i]->SetStats(kFALSE);

            printf("\n");
            histos[i]->Write();
            decayHistos[i]->Write();
            bgHistos[i]->Write();
        }
    }
 //output->Close();
}
