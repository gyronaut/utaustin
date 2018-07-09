#include <string>

struct particle{
    int pdg;
    int process;
    int parent;
    int tracknum;
    TLorentzVector pvec;
};

struct resonance{
    int daughter1;
    int daughter2;
    TLorentzVector pvec;
};

double calcMass(double px, double py, double pz, double E){
    return (double)TMath::Sqrt(E*E - (px*px + py*py + pz*pz));
};
double pT(double px, double py){
    return (double)TMath::Sqrt(px*px+py*py);
};
double eta(double px, double py, double pz){
    double p = TMath::Sqrt(px*px + py*py + pz*pz);
    return (double)0.5*TMath::Log((p + pz)/(p - pz));
};
double calcphi(double px, double py){
    return TMath::ATan2(py, px);
};

int nextline(ifstream* input); 

void phsdReader(string parFileName, string dataFileName){

    TH1D* massPhi = new TH1D("massPhi", "massPhi", 200, 0.9, 1.1);
    TH1D* massRealPhi = new TH1D("massRealPhi", "massRealPhi", 300, 0.9, 1.2);
    TH1D* massDecayPhi = new TH1D("massDecayPhi", "massDecayPhi", 300, 0.9, 1.2);
    TH1D* massKst = new TH1D("massKst", "massKst", 100 ,0.0, 2.0);
    TH1D* massRealKst = new TH1D("massRealKst", "massRealKst", 100, 0.0, 2.0);
    TH1D* massDecayKst = new TH1D("massDecayKst", "massDecayKst", 100, 0.0, 2.0);
    TH1D* pTpicharged = new TH1D("ptpicharged", "ptpicharged", 100, 0.0, 10.0);
    TH1D* pTkcharged = new TH1D("ptkcharged", "ptkcharged", 100, 0.0, 10.0);
    TH1D* pTproton = new TH1D("ptproton", "ptproton", 100, 0.0, 10.0);
    TH1D* pTcharged = new TH1D("ptcharged", "ptcharged", 100, 0.0, 10.0);
    TH1D* pTphi = new TH1D("ptphi", "ptphi", 100, 0.0, 10.0);
    TH1D* pTdecayphi = new TH1D ("ptdecahphi", "ptdecayphi", 100, 0.0, 10.0);
    TH1D* etaphi = new TH1D("etaphi", "etaphi", 100, -3.0, 3.0);
    TH1D* etah = new TH1D("etah", "etah", 100, -3.0, 3.0);
    TH1D* phiphi = new TH1D("phiphi", "phiphi", 100, -3.1416, 3.1416);
    TH1D* eventMult = new TH1D("eventMult", "eventMult", 2000, 0, 2000);
    TH1D* hphi = new TH1D("hphi", "hphi", 32, -3.14159/2.0, 3.0*3.1416/2.0);
    TH1D* hh = new TH1D("hh", "hh", 32, -3.14159/2.0, 3.0*3.1416/2.0);
    
    TH3D* chargedDist = new TH3D("chargedDist", "chargedDist; p_{T}; eta; phi",100, 0.0, 10.0, 100, -3.0, 3.0, 64, -3.1416, 3.1416); 
    TH3D* phiDist = new TH3D("phiDist", "phiDist; p_{T}; eta; phi",100, 0.0, 10.0, 100, -3.0, 3.0, 64, -3.1416, 3.1416); 
  
    Int_t dphi_bins[4]=    {20,     20,    64,   64};
    Double_t dphi_min[4] = {0.0,   0.0, -1.57, -2.0};
    Double_t dphi_max[4] = {10.0, 10.0,  4.71,  2.0};
    THnSparseF* hphiSparse = new THnSparseF("hphiSparse", "hphiSparse", 4, dphi_bins, dphi_min, dphi_max);
    THnSparseF* hhSparse = new THnSparseF("hhSparse", "hhSparse", 4, dphi_bins, dphi_min, dphi_max);

    //different particles are stored in separate vectors to be looped over after initial read-in
    std::vector<particle> posPions;
    std::vector<particle> negPions;
    std::vector<particle> posKaons;
    std::vector<particle> negKaons;
    std::vector<particle> neutralPions;
    std::vector<particle> neutralKaons;
    std::vector<particle> charged;
    std::vector<resonance> realPhi;
    
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
    //get numPar, numRuns from the inputPHSD File to know how many events to loop over
    int numPar=0, numRuns=0, numCharged=0;
    input >> numPar;
    nextline(&input);
    input >> numRuns;

    input.close();

    printf("NUM: %d, ISUBS: %d\n", numPar, numRuns);

    //next read in the output file phsd.dat
    input.open(dataFileName.c_str());
    if(!input.good()){
        printf("Invalid data file! Check name: %s\n", dataFileName.c_str());
    }
    //for each event (total events = numPar*numRuns), loop over all particle in the event (nPart)
    int nPart=0, pdg=0, charge=0, process=0, parent=0, nkst = 0; // process = ID5 = ID(j,5); parent = ID3 = ID(j,3)
    double px=0., py=0., pz=0., E=0., kpx=0., kpy=0., kpz=0., kE=0.;
    for(int ievt=1; ievt<numPar*numRuns+1; ievt++){ //loop over each event, where the number of events is numPar*numRuns
        input >> nPart;
        nextline(&input);
        nextline(&input);
        //printf("Event %d, nPart = %d\n", ievt, nPart);
        posPions.clear();
        negPions.clear();
        posKaons.clear();
        negKaons.clear();
        neutralPions.clear();
        neutralKaons.clear();
        charged.clear();
        realPhi.clear();
        numCharged=0;
        int npi = 0;

        for(int ipart = 0; ipart < nPart; ipart++){ //for each event, loop over each particle in the event and fill relevant histos
            input >> pdg >> charge >> px >> py >> pz >> E >> process >> parent;
            nextline(&input);
            particle part;
            part.process = process;
            part.parent = parent;
            part.pdg = pdg;
            part.tracknum = ipart;
            part.pvec.SetPx(px);
            part.pvec.SetPy(py);
            part.pvec.SetPz(pz);
            part.pvec.SetE(E);
            //fill the different particle vectors with the different species
            if(pdg == 211){
                if(pT(px, py) > 0.15 && TMath::Abs(eta(px, py, pz)) < 2.0){
                    posPions.push_back(part);
                    pTpicharged->Fill(pT(px, py));
                    pTcharged->Fill(pT(px, py));
                    etah->Fill(eta(px, py, pz));
                    charged.push_back(part);
                }
            }else if(pdg == -211){

                if(pT(px, py) > 0.15 && TMath::Abs(eta(px, py, pz)) < 2.0){
                    charged.push_back(part);
                    pTpicharged->Fill(pT(px, py));
                    etah->Fill(eta(px, py, pz));
                    pTcharged->Fill(pT(px, py));
                    negPions.push_back(part);
                }
            }else if(pdg == 111){
                neutralPions.push_back(part);
            }else if(pdg == 321){

                if(pT(px, py) > 0.15 && TMath::Abs(eta(px, py, pz)) < 2.0){
                    charged.push_back(part);
                    posKaons.push_back(part);
                    pTkcharged->Fill(pT(px, py));
                    etah->Fill(eta(px, py, pz));
                    pTcharged->Fill(pT(px, py));
                }
            }else if(pdg == -321){
                if(pT(px, py) > 0.15 && TMath::Abs(eta(px, py, pz)) < 2.0){
                    charged.push_back(part);
                    negKaons.push_back(part);
                    pTpicharged->Fill(pT(px, py));
                    etah->Fill(eta(px, py, pz));
                    pTcharged->Fill(pT(px, py));
                }
            }else if(TMath::Abs(pdg) == 311){
                neutralKaons.push_back(part);
            }else if(TMath::Abs(pdg) == 2212){
                if(pT(px, py) > 0.15 && TMath::Abs(eta(px, py, pz)) < 2.0){
                    charged.push_back(part);
                    pTproton->Fill(pT(px, py));
                    etah->Fill(eta(px, py, pz));
                    pTcharged->Fill(pT(px, py));
                }
            }
        }
        eventMult->Fill(charged.size());
        //printf("Num Kaons: %d,     Num Pions: %d\n", (int)(posKaons.size() + negKaons.size() + neutralKaons.size()), (int)(posPions.size() + negPions.size() + neutralPions.size()));

        nkst = 0;
        // after vectors are filled and all particles read in, loop over the kaons and pions to look
        // for pairs that came from the same K*, and loop of unlike sign kaons to reconstruct all phi mesons
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
                        if(posKaons[ikch].process == 7 && neutralPions[ipi0].process == 7){
                            massDecayKst->Fill(mass, 1.0/double(0.02*numPar*numRuns));
                        }
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
                    massPhi->Fill(mass);
                    if(posKaons[ikch].parent == negKaons[ikneg].parent){ 
                        massRealPhi->Fill(mass);
                        pTphi->Fill(pT(kpx, kpy));
                        etaphi->Fill(eta(kpx, kpy, kpz));
                        phiphi->Fill(calcphi(kpx, kpy));
                        if(posKaons[ikch].process == 5 && negKaons[ikneg].process == 5){
                            massDecayPhi->Fill(mass);
                            pTdecayphi->Fill(pT(kpx, kpy));
                            resonance phi;
                            phi.pvec.SetPx(kpx);
                            phi.pvec.SetPy(kpy);
                            phi.pvec.SetPz(kpz);
                            phi.pvec.SetE(kE);
                            phi.daughter1 = posKaons[ikch].tracknum;
                            phi.daughter2 = negKaons[ikneg].tracknum;
                            realPhi.push_back(phi); //add phi to phi vector
                            phiDist->Fill(pT(kpx, kpy), eta(kpx, kpy, kpz), calcphi(kpx, kpy));
                        }
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
                        if(posKaons[ik0].process == 7 && neutralPions[ipich].process == 7){
                            massDecayKst->Fill(mass, 1.0/double(0.02*numPar*numRuns));
                        }
                    }
                }
            }
        }
        //printf("Kstars: %d\n", nkst);
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

        //do simple correlations for h-phi and hh using charged hadron trigger
        for(int itrig = 0; itrig < charged.size(); itrig++){
            double trigpt = pT(charged[itrig].pvec.Px(), charged[itrig].pvec.Py());
            double trigeta = eta(charged[itrig].pvec.Px(), charged[itrig].pvec.Py(), charged[itrig].pvec.Pz());
            double trigphi = calcphi(charged[itrig].pvec.Px(), charged[itrig].pvec.Py());
            chargedDist->Fill(trigpt, trigeta, trigphi);
            //do h-phi correlations
            if(realPhi.size() !=0){
                for(int iphi = 0; iphi < realPhi.size(); iphi++){
                    if(charged[itrig].tracknum == realPhi[iphi].daughter1 || charged[itrig].tracknum == realPhi[iphi].daughter2) continue;
                    double assocpt = pT(realPhi[iphi].pvec.Px(), realPhi[iphi].pvec.Py());
                    double assoceta = eta(realPhi[iphi].pvec.Px(), realPhi[iphi].pvec.Py(), realPhi[iphi].pvec.Pz());
                    double assocphi = calcphi(realPhi[iphi].pvec.Px(), realPhi[iphi].pvec.Py());
                    double deltaphi = trigphi - assocphi; 
                    if(deltaphi > TMath::Pi()*3.0/2.0){
                        deltaphi += -2.0*TMath::Pi();
                    }else if(deltaphi < -0.5*TMath::Pi()){
                        deltaphi += 2.0*TMath::Pi();
                    }
                    double deltaeta = trigeta - assoceta;
                    Double_t point[4] = {trigpt, assocpt, deltaphi, deltaeta};
                    hphi->Fill(deltaphi);
                    hphiSparse->Fill(point);
                }
            }
            //do h-h correlations
            for(int iassoc = 0; iassoc < charged.size(); iassoc++){
                if(itrig == iassoc) continue;
                double assocPT = pT(charged[iassoc].pvec.Px(), charged[iassoc].pvec.Py());
                double assoceta = eta(charged[iassoc].pvec.Px(), charged[iassoc].pvec.Py(), charged[iassoc].pvec.Pz());
                double assocphi = calcphi(charged[iassoc].pvec.Px(), charged[iassoc].pvec.Py());
                double deltaeta = trigeta - assoceta;
                double deltaphi = trigphi - assocphi;
                if(deltaphi > TMath::Pi()*3.0/2.0){
                    deltaphi += -2.0*TMath::Pi();
                }else if(deltaphi < -0.5*TMath::Pi()){
                    deltaphi += 2.0*TMath::Pi();
                }

                if( assocPT < 3.0 && assocPT > 1.0){
                    hh->Fill(deltaphi);
                }
                Double_t point[4] = {trigpt, assocPT, deltaphi, deltaeta};
                hhSparse->Fill(point);
            }
        }

    }//end event loop 

        
    input.close();

    //open an output file to write out histograms
    TFile* output = new TFile("testoutput.root", "RECREATE");
    massKst->Write();
    massRealKst->Write();
    massDecayKst->Write();
    massPhi->Write();
    massRealPhi->Write();
    massDecayPhi->Write();
    pTpicharged->Write();
    pTkcharged->Write();
    pTproton->Write();
    pTcharged->Write();
    pTphi->Write();
    pTdecayphi->Write();
    etaphi->Write();
    phiphi->Write();
    eventMult->Write();
    hphi->Write();
    hh->Write();
    hphiSparse->Write();
    hhSparse->Write();

    //do some quick projections of the sparses and write out the results to file for ease of access
    hhSparse->GetAxis(0)->SetRangeUser(4.0, 10.0);
    hhSparse->GetAxis(1)->SetRangeUser(1.5, 4.0);
    TH1D* hhtrig4assoc1o5_4 = hhSparse->Projection(2);

    hphiSparse->GetAxis(0)->SetRangeUser(4.0, 10.0);
    hphiSparse->GetAxis(1)->SetRangeUser(1.5, 10.0);
    TH1D* hphitrig4assoc1o5_10 = hphiSparse->Projection(2);

    hhtrig4assoc1o5_4->Write();
    hphitrig4assoc1o5_10->Write();
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

