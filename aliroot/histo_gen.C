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
    TH2F *DphiPhiH = new TH2F("DphiPhiH", "#Delta#phi correlation for Phi-Hadron",200,0,20,64,-1.57,4.71);
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

    //Loop over all events
    for(Int_t nev=0; nev<numEvents; nev++){
        //if(nev%100 == 0)printf("Event: %d\n", nev);
        rl->GetEvent(nev);
        AliStack* stack = rl->Stack();
        Int_t npart = stack->GetNprimary();
        Double_t pt=-1, E=-1, theta=-1, phi=-1, eta=-99, y=-99;
        Double_t assoPt=-1, assoPhi=-1, Dphi=-99;
        Int_t parPdg = 0, assoPdg=0;
        //loop over all particles in stack
        for(Int_t part=0; part<npart; part++){
            TParticle *particle = stack->Particle(part);
            pt = particle->Pt();
            E = particle->Energy();
            theta = particle->Theta();
            phi = particle->Phi();
            eta = particle->Eta();
            parPdg = particle->GetPdgCode();
            //select only hadrons within the eta range -0.9 < eta <0.9
            if((TMath::Abs(parPdg)==2212 || TMath::Abs(parPdg)==211 || TMath::Abs(parPdg)==321 || TMath::Abs(parPdg)==11 || TMath::Abs(parPdg)==333)&&(TMath::Abs(eta)<0.9)){
                hadronPt->Fill(pt);
                if(TMath::Abs(parPdg)==333){                  
                    phiPt->Fill(pt);
                    //loop over all particles that aren't this phi-meson to compute delta-phi
                    for(Int_t apart = 0; apart<npart; apart++){
                        if(apart == part) continue;

                        TParticle *asso = stack->Particle(apart);
                        assoPt = asso->Pt();
                        assoPhi = asso->Phi();
                        assoPdg = asso->GetPdgCode();
                        //select just hadrons:
                        if(TMath::Abs(assoPdg)==2212 ||TMath::Abs(assoPdg)==11 ||TMath::Abs(assoPdg)==333 ||TMath::Abs(assoPdg)==321 || TMath::Abs(assoPdg)==211){
                            Dphi= phi - assoPhi;
                            if(Dphi >3*TMath::Pi()/2){
                                Dphi = Dphi - 2*TMath::Pi();
                            }
                            if(Dphi < -TMath::Pi()/2){
                                Dphi = Dphi + 2*TMath::Pi();
                            }
                            DphiPhiH->Fill(pt, Dphi);
                        }
                    }
                }else{
                    //loop over all particles that aren't this hadron to compute delta-phi
                    for(Int_t apart = 0; apart<npart; apart++){
                        if(apart == part) continue;
                        TParticle *hAsso = stack->Particle(apart);
                        assoPt = hAsso->Pt();
                        assoPhi = hAsso->Phi();
                        assoPdg = hAsso->GetPdgCode();
                        //select just phi mesons in the eta range: |eta| < 0.9
                        if((TMath::Abs(assoPdg)==333) && (TMath::Abs(eta)< 0.9)){
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
                            if(firstMotherIndex != 0){
                                TParticle *firstmother = stack->Particle(firstMotherIndex);
                                if(TMath::Abs(firstmother->GetPdgCode()) == 333) isHadronDaughter = true;
                            }
                            if(lastMotherIndex != 0){
                                TParticle *lastmother = stack->Particle(lastMotherIndex);
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
    rl->~AliRunLoader(); 


    //open new file to contain histograms
    string outputFullPath = outputDir+"/"+outputFile;
    TFile *histoOutput = new TFile(outputFullPath.c_str(), "RECREATE");
    histoOutput->cd();

    //write histograms to file
    hadronPt->Write();
    phiPt->Write();
    DphiPhiH->Write();
    DphiHPhi->Write();

    histoOutput->Close();

}
