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

void deltaReader(string parFileName, string dataFileName){

    TH1D* massPhi = new TH1D("massPhi", "massPhi", 200, 0.9, 1.1);
    TH1D* massRealPhi = new TH1D("massRealPhi", "massRealPhi", 200, 0.9, 1.1);
    TH1D* massKst = new TH1D("massKst", "massKst", 100 ,0.0, 2.0);
    TH1D* massRealKst = new TH1D("massRealKst", "massRealKst", 100, 0.0, 2.0);
    TH1D* pTpicharged = new TH1D("ptpicharged", "ptpicharged", 100, 0.0, 10.0);

    std::vector<particle> uuuDelta;
    std::vector<particle> uudDelta;
    std::vector<particle> uddDelta;
    std::vector<particle> dddDelta;
    std::vector<particle> antiuuuDelta;
    std::vector<particle> antiuudDelta;
    std::vector<particle> antiuddDelta;
    std::vector<particle> antidddDelta;


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
    int pid=0, charge=0, process=0, isub=0, irun=0; // process = ID5 = ID(j,5); parent = ID3 = ID(j,3)
    double px=0., py=0., pz=0., E=0., m=0., prodTime=0., decayTime=0., impact=0.;
    for(int ievt=1; ievt<numPar*numRuns+1; ievt++){
        posDelta.clear();
        negDelta.clear();
        neutralDelta.clear();
        int npi = 0;
        for(int ipart = 0; ipart < nPart; ipart++){
            input >> pid >> charge >> isub >> irun >> px >> py >> pz >> E >> m >> prodTime >> decayTime >> impact >> process;
            particle part;
            part.process = process;
            part.parent = parent;
            part.pdg = pdg;
            part.pvec.SetPx(px);
            part.pvec.SetPy(py);
            part.pvec.SetPz(pz);
            part.pvec.SetE(E);
            if(pdg == 211){
                posPions.push_back(part);
            }else if(pdg == -211){
                negPions.push_back(part);
            }else if(pdg == 111){
                neutralPions.push_back(part);
            }else if(pdg == 321){
                posKaons.push_back(part);
            }else if(pdg == -321){
                negKaons.push_back(part);
            }else if(TMath::Abs(pdg) == 311){
                neutralKaons.push_back(part);
            }
        }
        printf("Num Kaons: %d,     Num Pions: %d\n", (int)(posKaons.size() + negKaons.size() + neutralKaons.size()), (int)(posPions.size() + negPions.size() + neutralPions.size()));
        if(posPions.size() == 0 && (neutralPions.size() == 0 || posKaons.size() == 0)){
           continue;
        }

        nkst = 0;
        //loop over the kaons and pions to look for pairs that came from the same K*
        for(int ikch = 0; ikch < posKaons.size(); ikch++){
            for(int ipi0 = 0; ipi0 < neutralPions.size(); ipi0++){

                kpx = posKaons[ikch].pvec.Px() + neutralPions[ipi0].pvec.Px();
                kpy = posKaons[ikch].pvec.Py() + neutralPions[ipi0].pvec.Py();
                kpz = posKaons[ikch].pvec.Pz() + neutralPions[ipi0].pvec.Pz();
                kE = posKaons[ikch].pvec.E() + neutralPions[ipi0].pvec.E();
                double mass = calcMass(kpx, kpy, kpz, kE);
                float y = 0.5*TMath::Log((kE - kpz)/(kE + kpz));      
                if(TMath::Abs(y) <= 0.5){
                    massKst->Fill(mass, 1.0/double (0.02*numPar*numRuns));
                    if(posKaons[ikch].parent == neutralPions[ipi0].parent){
                        nkst+=1;
                        massRealKst->Fill(mass, 1.0/double(0.02*numPar*numRuns));
                    }
                }
            }
            for(int ikneg = 0; ikneg < negKaons.size(); ikneg++){
                kpx = posKaons[ikch].pvec.Px() + negKaons[ikneg].pvec.Px();
                kpy = posKaons[ikch].pvec.Py() + negKaons[ikneg].pvec.Py();
                kpz = posKaons[ikch].pvec.Pz() + negKaons[ikneg].pvec.Pz();
                kE = posKaons[ikch].pvec.E() + negKaons[ikneg].pvec.E();
                double mass = calcMass(kpx, kpy, kpz, kE);
                float y = 0.5*TMath::Log((kE - kpz)/(kE + kpz));
                if(TMath::Abs(y) <=0.5){
                    massPhi->Fill(mass, 1.0/(0.02*numPar*numRuns));
                    if(posKaons[ikch].parent == negKaons[ikneg].parent && (posKaons[ikch].process == 5 && negKaons[ikch].process == 5)){
                        massRealPhi->Fill(mass, 1.0/(0.02*numPar*numRuns));
                    }
                }

            }
        }

        for(int ik0 = 0; ik0 < neutralKaons.size(); ik0++){
            if(neutralKaons[ik0].pdg != 311) continue;
            for(int ipich = 0; ipich < posPions.size(); ipich++){
                if(posPions[ipich].pdg != 211) continue;
                kpx = neutralKaons[ik0].pvec.Px() + posPions[ipich].pvec.Px();
                kpy = neutralKaons[ik0].pvec.Py() + posPions[ipich].pvec.Py();
                kpz = neutralKaons[ik0].pvec.Pz() + posPions[ipich].pvec.Pz();
                kE = neutralKaons[ik0].pvec.E() + posPions[ipich].pvec.E();
                double mass = calcMass(kpx, kpy, kpz, kE);
                float y = 0.5*TMath::Log((kE - kpz)/(kE + kpz));

                if(TMath::Abs(y) <= 0.5){
                    massKst->Fill(mass, 1.0/double(0.02*numPar*numRuns));
                    if(neutralKaons[ik0].parent == posPions[ipich].parent){
                        nkst+=1;
                        massRealKst->Fill(mass, 1.0/double(0.02*numPar*numRuns));
                    }
                }
            }
        }
        printf("Kstars: %d\n", nkst);
//charged pion pT
        double piPx=0.0, piPy=0.0, piPt = 0.0;
        for(int ipich = 0; ipich < posPions.size(); ipich++){
            piPx = posPions[ipich].pvec.Px();
            piPy = posPions[ipich].pvec.Py();
            piPt = TMath::Sqrt(piPx*piPx + piPy*piPy);
            pTpicharged->Fill(piPt, 1.0/double(0.1*numPar*numRuns));
        }
        for(int ipich = 0; ipich < negPions.size(); ipich++){
            piPx = negPions[ipich].pvec.Px();
            piPy = negPions[ipich].pvec.Py();
            piPt = TMath::Sqrt(piPx*piPx + piPy*piPy);
            pTpicharged->Fill(piPt, 1.0/double(0.1*numPar*numRuns));
        }

            
    } 
    
    input.close();

    TFile* output = new TFile("testoutput.root", "RECREATE");
    massKst->Write();
    massRealKst->Write();
    massPhi->Write();
    massRealPhi->Write();
    pTpicharged->Write();
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

