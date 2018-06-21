#include <string>

int nextline(ifstream* input); 

void phsdReader(string parFileName, string dataFileName){

    std::vector<TLorentzVector> chargedPions;
    std::vector<TLorentzVector> chargedKaons;
    std::vector<TLorentzVector> neutralPions;
    std::vector<TLorentzVector> neutralKaons;

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
    int nPart=0, pdg=0, charge=0, parent=0;
    float px=0., py=0., pz=0., E=0.;
    for(int ievt=+1; ievt<numPar*numRuns+1; ievt++){
        input >> nPart;
        nextline(&input);
        nextline(&input);
        printf("Event %d, nPart = %d\n", ievt, nPart);
        for(int ipart = 0; ipart < nPart; ipart++){
            input >> pdg >> charge >> px >> py >> pz >> E >> parent;
            nextline(&input);
        }
    } 
    
    input.close();

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
