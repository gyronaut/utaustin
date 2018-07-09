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
    double px = 0.0, py = 0.0, pz = 0.0, E = 0.0, mass = 0.0, calcpt = 0.0;
    int id = 0, charge = 0, isub = 0, irun = 0;
    while(input.good()){
        input >> id >> charge >> isub >> irun >> px >> py >> pz >> E >> mass;
        if(nextline(&input)!=0) break;
        if(input.good()){
            calcpt = pT(px, py);
            if(id > 0){
                switch(charge){
                    case 0: massuddDelta->Fill(mass, calcpt);
                    case 1: massuudDelta->Fill(mass, calcpt);
                    case 2: massuuuDelta->Fill(mass, calcpt);
                    case -1: massdddDelta->Fill(mass, calcpt);
                }
                massallDelta->Fill(mass, calcpt);
            }else{
                switch(charge){
                    case 0: massantiuddDelta->Fill(mass, calcpt);
                    case -1: massantiuudDelta->Fill(mass, calcpt);
                    case -2: massantiuuuDelta->Fill(mass, calcpt);
                    case 1: massantidddDelta->Fill(mass, calcpt);
                }
                massallantiDelta->Fill(mass, calcpt);
            }
            massTotalDelta->Fill(mass, calcpt);
        }
    }
    
    input.close();

    TFile* output = new TFile("deltaoutput.root", "RECREATE");
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

