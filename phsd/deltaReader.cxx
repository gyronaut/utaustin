#include <string>

struct particle{
    int pdg;
    int process;
    int parent;
    TLorentzVector pvec;
};

double calcMass(double px, double py, double pz, double E){
    return (double)TMath::Sqrt(E*E - (px*px + py*py + pz*pz));
}
double pT(double px, double py){
    return (double)TMath::Sqrt(px*px + py*py);
}
int nextline(ifstream* input); 

void deltaReader(string dataFileName){

    TH2D* massuuuDelta = new TH2D("massuuuDelta", "massuuuDelta", 800, 1.00, 1.4, 5, 0.5, 3.0);
    TH2D* massuudDelta = new TH2D("massuudDelta", "massuudDelta", 800, 1.00, 1.4, 5, 0.5, 3.0);
    TH2D* massuddDelta = new TH2D("massuddDelta", "massuddDelta", 800, 1.00, 1.4, 5, 0.5, 3.0);
    TH2D* massdddDelta = new TH2D("massdddDelta", "massdddDelta", 800, 1.00, 1.4, 5, 0.5, 3.0);
    TH2D* massantiuuuDelta = new TH2D("massantiuuuDelta", "massantiuuuDelta", 800, 1.00, 1.4, 5, 0.5, 3.0);
    TH2D* massantiuudDelta = new TH2D("massantiuudDelta", "massantiuudDelta", 800, 1.00, 1.4, 5, 0.5, 3.0);
    TH2D* massantiuddDelta = new TH2D("massantiuddDelta", "massantiuddDelta", 800, 1.00, 1.4, 5, 0.5, 3.0);
    TH2D* massantidddDelta = new TH2D("massantidddDelta", "massantidddDelta", 800, 1.00, 1.4, 5, 0.5, 3.0);
    TH2D* massallDelta = new TH2D("massallDelta", "massallDelta", 800, 1.00, 1.4, 5, 0.5, 3.0);
    TH2D* massallantiDelta = new TH2D("massallantiDelta", "massallantiDelta", 800, 1.00, 1.4, 5, 0.5, 3.0);
    TH2D* massTotalDelta = new TH2D("massTotalDelta", "massTotalDelta", 800, 1.00, 1.4, 5, 0.5, 3.0);

    TH1D* massEarlyDelta = new TH1D("massEarlyDelta", "massEarlyDelta", 800, 1.00, 1.4);
    TH1D* massLateDelta = new TH1D("massLateDelta", "massLateDelta", 800, 1.00, 1.4); 
 
    TH1D* mass1Delta = new TH1D("mass1Delta", "mass1Delta", 800, 1.00, 1.4);    
    TH1D* mass17Delta = new TH1D("mass17Delta", "mass17Delta", 800, 1.00, 1.4); 
    TH1D* mass19Delta = new TH1D("mass19Delta", "mass19Delta", 800, 1.00, 1.4); 
    TH1D* mass92Delta = new TH1D("mass92Delta", "mass92Delta", 800, 1.00, 1.4); 
    TH1D* mass93Delta = new TH1D("mass93Delta", "mass93Delta", 800, 1.00, 1.4); 
    TH1D* mass98Delta = new TH1D("mass98Delta", "mass98Delta", 800, 1.00, 1.4); 
    TH1D* mass200Delta = new TH1D("mass200Delta", "mass200Delta", 800, 1.00, 1.4); 
    TH1D* mass201Delta = new TH1D("mass201Delta", "mass201Delta", 800, 1.00, 1.4);
    TH1D* massOtherDelta = new TH1D("massOtherDelta", "massOtherDelta", 800, 1.00, 1.4);
    TH1D* massNot92Delta = new TH1D("massNot92Delta", "massNot92Delta", 800, 1.00, 1.4); 


    TH1D* prodHist = new TH1D("prodHist", "prodHist", 500, 0, 500);

    std::vector<particle> uuuDelta;
    std::vector<particle> uudDelta;
    std::vector<particle> uddDelta;
    std::vector<particle> dddDelta;
    std::vector<particle> antiuuuDelta;
    std::vector<particle> antiuudDelta;
    std::vector<particle> antiuddDelta;
    std::vector<particle> antidddDelta;


    //next read in the output file
    ifstream input;
    input.open(dataFileName.c_str());
    if(!input.good()){
        printf("Invalid data file! Check name: %s\n", dataFileName.c_str());
    }
    double px = 0.0, py = 0.0, pz = 0.0, E = 0.0, mass = 0.0, calcpt = 0.0, prodtime=0.0, decaytime=0.0, impact=0.0;
    int id = 0, charge = 0, isub = 0, irun = 0, prodchannel = 0;
    while(input.good()){
        input >> id >> charge >> isub >> irun >> px >> py >> pz >> E >> mass >> prodtime >> decaytime >> impact >> prodchannel;
        if(nextline(&input)!=0) break;
        if(input.good()){
            calcpt = pT(px, py);
            if(id > 0){
                switch(charge){
                    case 0: massuddDelta->Fill(mass, calcpt); break;
                    case 1: massuudDelta->Fill(mass, calcpt); break;
                    case 2: massuuuDelta->Fill(mass, calcpt); break;
                    case -1: massdddDelta->Fill(mass, calcpt); break;
                }
                massallDelta->Fill(mass, calcpt);
            }else{
                switch(charge){
                    case 0: massantiuddDelta->Fill(mass, calcpt); break;
                    case -1: massantiuudDelta->Fill(mass, calcpt); break;
                    case -2: massantiuuuDelta->Fill(mass, calcpt); break;
                    case 1: massantidddDelta->Fill(mass, calcpt); break;
                }
                massallantiDelta->Fill(mass, calcpt);
            }
            massTotalDelta->Fill(mass, calcpt);
            if(decaytime > 1000.0){
                massLateDelta->Fill(mass);
            }else{
                massEarlyDelta->Fill(mass);
            }
            switch(prodchannel){
                case 1: mass1Delta->Fill(mass); break;
                case 17: mass17Delta->Fill(mass); break;
                case 19: mass19Delta->Fill(mass); break;
                case 92: mass92Delta->Fill(mass); break;
                case 93: mass93Delta->Fill(mass); break;
                case 98: mass98Delta->Fill(mass); break;
                case 200: mass200Delta->Fill(mass); break;
                case 201: mass201Delta->Fill(mass); break;
                default: massOtherDelta->Fill(mass); break;
            }
            if(prodchannel !=92) massNot92Delta->Fill(mass);
            prodHist->Fill(prodchannel);
        }
    }
    
    input.close();

    TFile* output = new TFile("testdeltaoutput.root", "RECREATE");
    massuuuDelta->Write();
    massuudDelta->Write();
    massuddDelta->Write();
    massdddDelta->Write();
    massantiuuuDelta->Write();
    massantiuudDelta->Write();
    massantiuddDelta->Write();
    massantidddDelta->Write();
    massallDelta->Write();
    massallantiDelta->Write();
    massTotalDelta->Write();

    massEarlyDelta->Write();
    massLateDelta->Write();
    mass1Delta->Write();
    mass17Delta->Write();
    mass19Delta->Write();
    mass92Delta->Write();
    mass93Delta->Write();
    mass98Delta->Write();
    mass200Delta->Write();
    mass201Delta->Write();
    massNot92Delta->Write();
    massOtherDelta->Write();

    prodHist->Write();

}

int nextline(ifstream* input){
    input->ignore(1000, '\n');
    if(input->eof()){
        return 1;
    }
    if(input->bad()){
        return -1;
    }
    return 0;
}

