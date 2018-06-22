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

int nextline(ifstream* input); 

void phsdReader(string parFileName, string dataFileName){

    TH1D* massKst = new TH1D("massKst", "massKst", 100 ,0.0, 2.0);
    TH1D* massRealKst = new TH1D("massRealKst", "massRealKst", 100, 0.0, 2.0);

    std::vector<particle> chargedPions;
    std::vector<particle> chargedKaons;
    std::vector<particle> neutralPions;
    std::vector<particle> neutralKaons;
    
    //first read in the input parameter file
    ifstream input;
    input.open(parFileName.c_str());
    if(!input.good()){
        printf("Invalid parameter file! Check name: %s\n", parFileName.c_str());
        return;
    }
    //skip first few parameters
    for(int i=0; i<8; i++){
        if(nextline(&input)!=0){
            printf("problem skipping lines in par file!\n");
            return;
        }
    }

    int numPar=0, numRuns=0;
    input >> numPar;
    nextline(&input);
    input >> numRuns;

    input.close();

    printf("NUM: %d, ISUBS: %d\n", numPar, numRuns);

    //next read in the output file
    input.open(dataFileName.c_str());
    if(!input.good()){
        printf("Invalid data file! Check name: %s\n", dataFileName.c_str());
    }
    //for each event (total events = numPar*numRuns), loop over all particle in the event (nPart)
    int nPart=0, pdg=0, charge=0, process=0, parent=0, nkst = 0;
    double px=0., py=0., pz=0., E=0., kpx=0., kpy=0., kpz=0., kE=0.;
    for(int ievt=1; ievt<numPar*numRuns+1; ievt++){
        input >> nPart;
        nextline(&input);
        nextline(&input);
        printf("Event %d, nPart = %d\n", ievt, nPart);
        chargedPions.clear();
        chargedKaons.clear();
        neutralPions.clear();
        neutralKaons.clear();
        int npi = 0;
        for(int ipart = 0; ipart < nPart; ipart++){
            input >> pdg >> charge >> px >> py >> pz >> E >> process >> parent;
            nextline(&input);
            particle part;
            part.process = process;
            part.parent = parent;
            part.pdg = pdg;
            part.pvec.SetPx(px);
            part.pvec.SetPy(py);
            part.pvec.SetPz(pz);
            part.pvec.SetE(E);
            if(TMath::Abs(pdg) == 211){
                chargedPions.push_back(part);
            }
            if(pdg == 111){
                neutralPions.push_back(part);
            }
            if(TMath::Abs(pdg) == 321){
                chargedKaons.push_back(part);
            }
            if(TMath::Abs(pdg) == 311){
                neutralKaons.push_back(part);
            }
        }
        printf("Num Kaons: %d,     Num Pions: %d\n", (int)(chargedKaons.size() + neutralKaons.size()), (int)(chargedPions.size() + neutralPions.size()));
        if(chargedPions.size() == 0 && (neutralPions.size() == 0 || chargedKaons.size() == 0)){
           continue;
        }

        nkst = 0;
        //loop over the kaons and pions to look for pairs that came from the same K*
        for(int ikch = 0; ikch < chargedKaons.size(); ikch++){
            if(chargedKaons[ikch].pdg!=321) continue;
            for(int ipi0 = 0; ipi0 < neutralPions.size(); ipi0++){

                kpx = chargedKaons[ikch].pvec.Px() + neutralPions[ipi0].pvec.Px();
                kpy = chargedKaons[ikch].pvec.Py() + neutralPions[ipi0].pvec.Py();
                kpz = chargedKaons[ikch].pvec.Pz() + neutralPions[ipi0].pvec.Pz();
                kE = chargedKaons[ikch].pvec.E() + neutralPions[ipi0].pvec.E();
                double mass = calcMass(kpx, kpy, kpz, kE);
                float y = 0.5*TMath::Log((kE - kpz)/(kE + kpz));      
                if(TMath::Abs(y) <= 0.5){
                    massKst->Fill(mass, 1.0/double (0.02*numPar*numRuns));
                    if(chargedKaons[ikch].parent == neutralPions[ipi0].parent){
                        nkst+=1;
                        massRealKst->Fill(mass, 1.0/double(0.02*numPar*numRuns));
                    }
                }
            }
        }

        for(int ik0 = 0; ik0 < neutralKaons.size(); ik0++){
            if(neutralKaons[ik0].pdg != 311) continue;
            for(int ipich = 0; ipich < chargedPions.size(); ipich++){
                if(chargedPions[ipich].pdg != 211) continue;
                kpx = neutralKaons[ik0].pvec.Px() + chargedPions[ipich].pvec.Px();
                kpy = neutralKaons[ik0].pvec.Py() + chargedPions[ipich].pvec.Py();
                kpz = neutralKaons[ik0].pvec.Pz() + chargedPions[ipich].pvec.Pz();
                kE = neutralKaons[ik0].pvec.E() + chargedPions[ipich].pvec.E();
                double mass = calcMass(kpx, kpy, kpz, kE);
                float y = 0.5*TMath::Log((kE - kpz)/(kE + kpz));

                if(TMath::Abs(y) <= 0.5){
                    massKst->Fill(mass, 1.0/double(0.02*numPar*numRuns));
                    if(neutralKaons[ik0].parent == chargedPions[ipich].parent){
                        nkst+=1;
                        massRealKst->Fill(mass, 1.0/double(0.02*numPar*numRuns));
                    }
                }
            }
        }
        printf("Kstars: %d\n", nkst);

    } 
    
    input.close();

    TFile* output = new TFile("testoutput.root", "RECREATE");
    massKst->Write();
    massRealKst->Write();
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

