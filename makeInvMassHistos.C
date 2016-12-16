#include <string>
#include <TH1F> 

int makeInvMassHistos(string filename){
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
    
    string ptLower = filename.substr(filename.find("[")+1, filename.find(",")-filename.find("[")-1);
    string ptHigher = filename.substr(filename.find(",")+1, filename.find("]")-filename.find(",")-1);

    TH1D* histos[8];
    for(int i=0; i<8; i++){
        histos[i] = new TH1D(Form("histo%d", i), Form("%s < p_{T} < %s, %s", ptLower.c_str(), ptHigher.c_str(), particles[i].c_str()), 100, 0.0, 2.0);
    }

    ifstream input;
    input.open(filename.c_str());
    string line = "";
    if(input.good()){
        getline(input, line);
    }

    printf("line: %s\n", line.c_str());

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
    printf("linenumber: %d\n", lineNumber);

    TFile *output = new TFile("output.root", "RECREATE");
    for(int i=0; i<8; i++){
        histos[i]->Write();
    }
    output->Close();
}
