using namespace std;

#include <string>
#include <TROOT.h>
#include "TH1F.h"
#include "TTree.h"
#include "TSystem.h"
#include "TClassTable.h"
#include "TObjArray.h"

void histo_gen(string inputDir, string inputFile, string outputDir, string outputFile){

    if (gClassTable->GetID("AliRun") < 0) {
        gROOT->LoadMacro("loadlibs.C");
        loadlibs();
    }

    gSystem->SetIncludePath("-I$ROOTSYS/include -I$ALICE_ROOT/include -I$ALICE_ROOT");
    gSystem->Load("libhijing.so");
    gSystem->Load("libTHijing.so");     // AliGenHijing is defined here
    gSystem->Load("liblhapdf.so");      // Parton density functions
    gSystem->Load("libEGPythia6.so");   // TGenerator interface
    gSystem->Load("libpythia6.so");     // Pythia
    gSystem->Load("libAliPythia6.so");  // ALICE specific implementations

    TH1F *hadronPt = new TH1F("HadronPt", "hadron pt distribution", 1000, 0, 50);
    TH1F *phiPt = new TH1F("PhiPt", "phi pt distribution", 1000, 0, 50);
    TH1F *phiPhiDist = new TH1F("PhiPhiDist", "phi meson angular phi distribution", 100, -0.1, 6.29);
    TH1F *hPhiDist = new TH1F("HPhiDist", "hadron angular phi distribution", 100, -0.1, 6.29);
    TH3F *DphiHPhi = new TH3F("DphiHPhi", "#Delta#phi correlation for Hadron-Phi", 1000, 0, 50, 1000, 0, 50, 64, -1.57, 4.71);
    string inputFullPath = inputDir+"/"+inputFile;
    AliRunLoader* rl = AliRunLoader::Open(inputFullPath.c_str());

    if(!rl){
        fprintf(stderr, "Run Loader not initialized, stopping program.");
        return;
    }

    Int_t numEvents = rl->GetNumberOfEvents();
    TDatabasePDG* DataBase = new TDatabasePDG();

    rl->LoadKinematics();
    rl->LoadHeader();
    Int_t count = 0;
    //Loop over all events
    for(Int_t nev=0; nev<numEvents; nev++){
        //if(nev%100 == 0)printf("Event: %d\n", nev);
        rl->GetEvent(nev);
        AliStack* stack = rl->Stack();
        Int_t npart = stack->GetNprimary();
        Double_t pt=-1, E=-1, theta=-1, phi=-1, eta=-99, y=-99, eta_calc=-99, p_tot = 0, p_z = 0;
        Double_t assoPt=-1, assoPhi=-1, assoEta=-9, Dphi=-99;
        Int_t parPdg = 0, assoPdg=0;
        TParticle *particle, *hAsso;
        //loop over all particles in stack
        for(Int_t part=0; part<npart; part++){
            particle = stack->Particle(part);
            pt = particle->Pt();
            E = particle->Energy();
            theta = particle->Theta();
            phi = particle->Phi();
            eta = particle->Eta();
            p_tot = particle->P();
            p_z = particle->Pz();
            parPdg = particle->GetPdgCode();
            eta_calc = 0.5*TMath::Log((p_tot + p_z)/(p_tot - p_z));
            //select only hadrons within the eta range -0.9 < eta <0.9 and pt > 150 MeV
            //if(nev%1000==0 && TMath::Abs(eta_calc)<1) printf("  Particle eta: %f, Calculated eta: %f, Pt: %f, PDG: %d\n", eta, eta_calc, pt, parPdg);
            if((TMath::Abs(parPdg)==2212 || TMath::Abs(parPdg)==211 || TMath::Abs(parPdg)==321 || TMath::Abs(parPdg)==11 || TMath::Abs(parPdg)==333)&&(TMath::Abs(eta_calc)<0.9)&&(pt > 0.150)){
                hadronPt->Fill(pt);
                count++;
                if(TMath::Abs(parPdg)==333){                  
                    phiPt->Fill(pt);
                    phiPhiDist->Fill(phi);
                }else{
                    hPhiDist->Fill(phi);
                    //loop over all particles that aren't this hadron to compute delta-phi
                    for(Int_t apart = 0; apart<npart; apart++){
                        if(apart == part) continue;
                        hAsso = stack->Particle(apart);
                        assoPt = hAsso->Pt();
                        assoPhi = hAsso->Phi();
                        assoEta = hAsso->Eta();
                        assoPdg = hAsso->GetPdgCode();
                        //select just phi mesons in the eta range: |eta| < 0.9 and pt > 150 MeV
                        if((TMath::Abs(assoPdg)==333) && (TMath::Abs(assoEta)< 0.9) && (assoPt > 0.150)){
                            //printf("Got Here! \n");
                            //check that hadron isn't daughter particle of any phi meson
                            Int_t numDaughters = hAsso->GetNDaughters();
                            bool isHadronDaughter = false;
                            /*
                            for(Int_t i=0; i<numDaughters; i++){
                                Int_t daughter = hAsso->GetFirstDaughter()+i; //first we check to make sure the hadron isn't a direct daughter of the meson in hAsso
                                if(daughter == part){
                                    isHadronDaughter = true;
                                }
                            }
                            */
                            Int_t firstMotherIndex = particle->GetMother(0); //we check to make sure the hadron isn't the daughter of any phi meson
                            Int_t secondMotherIndex = particle->GetMother(1);
                            if(firstMotherIndex > 0){
                                TParticle *firstmother = stack->Particle(firstMotherIndex);
                                if(TMath::Abs(firstmother->GetPdgCode()) == 333) isHadronDaughter = true;
                            }
                            if(secondMotherIndex > 0){
                                TParticle *lastmother = stack->Particle(secondMotherIndex);
                                if(TMath::Abs(lastmother->GetPdgCode()) == 333) isHadronDaughter = true;
                            }
                            if(isHadronDaughter) continue;
                            Dphi = phi - assoPhi;
                            if(Dphi >3*TMath::Pi()/2){
                                Dphi = Dphi - 2*TMath::Pi();
                            }
                            if(Dphi < -TMath::Pi()/2){
                                Dphi = Dphi + 2*TMath::Pi();
                            }
                            DphiHPhi->Fill(pt, assoPt, Dphi);
                        }
                    }
                }
            }
        }
    }
//    printf("Total hadrons counted in eta range: %d\n", count);
    rl->~AliRunLoader(); 


    //open new file to contain histograms
    string outputFullPath = outputDir+"/"+outputFile;
    TFile *histoOutput = new TFile(outputFullPath.c_str(), "RECREATE");
    histoOutput->cd();

    //write histograms to file
    hadronPt->Write();
    phiPt->Write();
    DphiHPhi->Write();

    histoOutput->Close();

}
