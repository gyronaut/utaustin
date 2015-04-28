using namespace std;

#include <string>
#include <TROOT.h>
#include "TH1F.h"
#include "TTree.h"
#include "TSystem.h"
#include "TClassTable.h"
#include "TObjArray.h"
#include "TStopwatch.h"

void histo_gen(string input_dir, string input_file, string output_dir, string output_file){
	if (gClassTable->GetID("AliRun") < 0) {
		gROOT->LoadMacro("loadlibs.C");
		loadlibs();
    }
    TStopwatch timer_io;

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

    TH2F *k0PhiDist = new TH2F("k0Dist", "Phi distribution for k0", 100, -0.1, 6.29);
    TH3F *DphiHK0 = new TH3F("DphiHK0", "#Delta#phi correlation for Hadron-K^{0}", 1000, 0, 50, 1000, 0, 50, 64, -1.57, 4.71);

    string input_full_path = input_dir+"/"+input_file;
    AliRunLoader* rl = AliRunLoader::Open(input_full_path.c_str());

    if(!rl){
        fprintf(stderr, "Run Loader not initialized, stopping program.");
        return;
    }

    Int_t num_events = rl->GetNumberOfEvents();
    
    TDatabasePDG* DataBase = new TDatabasePDG();

    printf("Number of Events: %d", num_events);

    rl->LoadKinematics();
    rl->LoadHeader();

    Int_t count = 0;
    Int_t npart = 0;
    Double_t pt=-1, E=-1, theta=-1, phi=-1, eta=-99, y=-99, eta_calc=-99, p_tot = 0, p_z = 0;
    Double_t asso_pt=-1, asso_phi=-1, asso_eta=-9, delta_phi=-99;
    Int_t par_pdg = 0, asso_pdg=0;
    TParticle *particle, *hAsso;
    AliStack* stack;

    timer_io.Stop();
    fprintf(stdout, "IO ");
    timer_io.Print("u");

    //Loop over all events
    for(Int_t nev=0; nev<num_events; nev++){
//        if(nev%100 == 0)printf("Event: %d\n", nev);
        TStopwatch timer_getevent;

        rl->GetEvent(nev);
        stack = rl->Stack();
        
        npart = stack->GetNprimary();
        
        //reset all variables before new event
        pt=-1, E=-1, theta=-1, phi=-1, eta=-99, y=-99, eta_calc=-99, p_tot = 0, p_z = 0;
        asso_pt= -1, asso_phi=-1, asso_eta=-9, delta_phi=-99;
        par_pdg = 0, asso_pdg = 0;
        
        timer_getevent.Stop();
        fprintf(stdout, "GETEVENT ");
        timer_getevent.Print("u");
        
        TStopwatch timer_analysis;
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
            par_pdg = particle->GetPdgCode();
            eta_calc = 0.5*TMath::Log((p_tot + p_z)/(p_tot - p_z));
            //select only hadrons within the eta range -0.9 < eta <0.9 and pt > 150 MeV
            //if(nev%1000==0 && TMath::Abs(eta_calc)<1) printf("  Particle eta: %f, Calculated eta: %f, Pt: %f, PDG: %d\n", eta, eta_calc, pt, par_pdg);
            if((TMath::Abs(par_pdg)==2212 || TMath::Abs(par_pdg)==311 || TMath::Abs(par_pdg)==211 || TMath::Abs(par_pdg)==321 || TMath::Abs(par_pdg)==11 || TMath::Abs(par_pdg)==333)&&(TMath::Abs(eta_calc)<0.9)&&(pt > 0.150)){
                hadronPt->Fill(pt);
                count++;
                if(TMath::Abs(par_pdg)==333){                  
                    phiPt->Fill(pt);
                    phiPhiDist->Fill(phi);
                }else if(TMath::Abs(par_pdg)==311){
                    k0PhiDist->Fill(phi);
                }else{
                    hPhiDist->Fill(phi);
                    //loop over all particles that aren't this hadron to compute delta-phi
                    for(Int_t apart = 0; apart<npart; apart++){
                        if(apart == part) continue;
                        hAsso = stack->Particle(apart);
                        asso_pt = hAsso->Pt();
                        asso_phi = hAsso->Phi();
                        asso_eta = hAsso->Eta();
                        asso_pdg = hAsso->GetPdgCode();
                        //select just phi mesons in the eta range: |eta| < 0.9 and pt > 150 MeV
                        if((TMath::Abs(asso_pdg)==333 || TMath::Abs(asso_pdg)==311) && (TMath::Abs(asso_eta)< 0.9) && (asso_pt > 0.150)){
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
                            Int_t firstMotherIndex = particle->GetMother(0); //we check to make sure the hadron isn't the daughter of any phi meson (k0)
                            Int_t secondMotherIndex = particle->GetMother(1);
                            if(TMath::Abs(asso_pdg)==333){
                                if(firstMotherIndex > 0){
                                    TParticle *firstmother = stack->Particle(firstMotherIndex);
                                    if(TMath::Abs(firstmother->GetPdgCode()) == 333) isHadronDaughter = true;
                                }
                                if(secondMotherIndex > 0){
                                    TParticle *lastmother = stack->Particle(secondMotherIndex);
                                    if(TMath::Abs(lastmother->GetPdgCode()) == 333) isHadronDaughter = true;
                                }
                            }else if(TMath::Abs(asso_pdg)==311){
                                if(firstMotherIndex > 0){
                                    TParticle *firstmother = stack->Particle(firstMotherIndex);
                                    if(TMath::Abs(firstmother->GetPdgCode()) == 311) isHadronDaughter = true;
                                }
                                if(secondMotherIndex > 0){
                                    TParticle *lastmother = stack->Particle(secondMotherIndex);
                                    if(TMath::Abs(lastmother->GetPdgCode()) == 311) isHadronDaughter = true;
                                }
                            }
                            if(isHadronDaughter) continue;
                            delta_phi = phi - asso_phi;
                            if(delta_phi >3*TMath::Pi()/2){
                                delta_phi = delta_phi - 2*TMath::Pi();
                            }
                            if(delta_phi < -TMath::Pi()/2){
                                delta_phi = delta_phi + 2*TMath::Pi();
                            }
                            if(TMath::Abs(asso_pdg)==333){
                                DphiHPhi->Fill(pt, asso_pt, delta_phi);
                            }else{
                                DphiHK0->Fill(pt, asso_pt, delta_phi);
                            }
                        }
                    }
                }
            }
        }
        timer_analysis.Stop();
        fprintf(stdout, "ANALYSIS ");
        timer_analysis.Print("u");
    }

    TStopwatch timer_write;

//    printf("Total hadrons counted in eta range: %d\n", count);
    rl->~AliRunLoader(); 


    //open new file to contain histograms
    string outputFullPath = output_dir+"/"+output_file;
    TFile *histoOutput = new TFile(outputFullPath.c_str(), "RECREATE");
    histoOutput->cd();

    //write histograms to file
    hadronPt->Write();
    phiPt->Write();
    DphiHPhi->Write();

    histoOutput->Close();

    timer_write.Stop();
    fprintf(stdout, "WRITE ");
    timer_write.Print("u");
}
