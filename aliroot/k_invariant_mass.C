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

    TH1F *invMassHisto = new TH1F("InvMassHisto", "Invariant Mass distribution for K+/K- pairs", 1000, 500, 1500);
    TH2F *invMassDPhi = new TH2F("InvMassDPhi", "Invariant Mass Distribution for K+/K- pairs in delta-phi bins", 1000, 500, 1500, 40, -1.6, 4.7);
    TH1F *invMassPhiOnly = new TH1F("InvMassPhiOnly", "Invariant Mass distribution for K+/K- pairs known to come from phi mesons", 1000, 500, 1500);
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
        Double_t pt_1=-1, E_1=0, m_1=0, px_1=0, py_1=0, pz_1=0 phi_1=-99, eta_1=-99;
        Double_t pt_2=-1, E_2=0, m_2=0, px_2=0, py_2=0, pz_2=0 phi_2=-99, eta_2=-99;
        Double_t invMass = 0;
        Int_t parPdg_1 = 0, parPdg_2 = 0;
        //loop over all particles in stack
        for(Int_t part=0; part<npart; part++){
            TParticle *particle = stack->Particle(part);
            pt_1 = particle->Pt();
            E_1 = particle->Energy();
            m_1 = particle->GetCalcMass();
            px_1 = particle->Px();
            py_1 = particle->Py();
            pz_1 = particle->Pz();
            phi_1 = particle->Phi();
            eta_1 = particle->Eta();
            parPdg_1 = particle->GetPdgCode();
            
            //select only K+ within the eta range -0.9 < eta <0.9
            if(parPdg_1==321 && TMath::Abs(eta_1)<0.9){
                //loop over all particles that aren't this K+ to look for K-
                for(Int_t apart = 0; apart<npart; apart++){
                    if(apart == part) continue;
                    TParticle *secondPart = stack->Particle(apart);
                    pt_2 = secondPart->Pt();
                    px_2 = secondPart->Px();
                    py_2 = secondPart->Py();
                    pz_2 = secondPart->Pz();
                    E_2 = secondPart->Energy;
                    m_2 = secondPart->GetCalcMass();
                    phi_2 = secondPart->Phi();
                    parPdg_2 = secondPart->GetPdgCode();
                    //select just K- in the eta range: |eta| < 0.9
                    if(parPdg_2 == -321 && TMath::Abs(eta_2)< 0.9){
                        //calculate invariant mass
                        invMass = TMath::Sqrt(m_1*m_1 + m_2*m_2 + 2*E_1*E_2 - 2*(px_1*px_2 + py_1*py_2 + pz_1*pz_2));
                        invMassHisto->Fill(invMass);

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
    DphiHPhi->Write();

    histoOutput->Close();

}
