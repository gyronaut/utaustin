using namespace std;

#include <string>
#include <TROOT.h>
#include "TH1F.h"
#include "TTree.h"
#include "TSystem.h"
#include "TClassTable.h"
#include "TObjArray.h"
#include "TStopwatch.h"

void K0_histo_gen(string input_dir, string input_file, string output_dir, string output_file){
	if (gClassTable->GetID("AliRun") < 0) {
		gROOT->LoadMacro("loadlibs.C");
		loadlibs();
    }
    TStopwatch timer_total;
    //    TStopwatch timer_io;

    gSystem->SetIncludePath("-I$ROOTSYS/include -I$ALICE_ROOT/include -I$ALICE_ROOT");
    gSystem->Load("libhijing.so");
    gSystem->Load("libTHijing.so");     // AliGenHijing is defined here
    gSystem->Load("liblhapdf.so");      // Parton density functions
    gSystem->Load("libEGPythia6.so");   // TGenerator interface
    gSystem->Load("libpythia6.so");     // Pythia
    gSystem->Load("libAliPythia6.so");  // ALICE specific implementations

    /* Bin limits for pT, phi, and eta THnSparse histos */
    Int_t dist_bins[3] = {100, 100, 100};
    Double_t dist_min[3] = {0.0, -0.1, -3};
    Double_t dist_max[3] = {20, 6.29, 3};

    /* PT, PHI, AND ETA DISTRIUBTIONS (pT, phi, eta)*/
    THnSparseF *hPtPhiEtaDist = new THnSparseF("hPtPhiEtaDist", "hadron angular phi distribution", 3, dist_bins, dist_min, dist_max);
 
    TH2F *K0daughterPT = new TH2F("K0daughterPT", "PT distribution for #pi^{+/-} daughters of K^{0}_{S}; p_{T}^{#pi^+}; p_{T}^{#pi^-}", 50, 0, 10, 50, 0, 10);

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
    Int_t ntrigger=0;
    Double_t pt=-1, E=-1, theta=-1, phi=-1, eta=-99, y=-99, eta_calc=-99, p_tot = 0, p_z = 0;
    Double_t point[3] = {0,0,0};
    Int_t par_pdg = 0, asso_pdg=0;
    TParticle *particle, *hAsso;
    AliStack* stack;

//    timer_io.Stop();
//    fprintf(stdout, "IO ");
//    timer_io.Print("u");

    /* Loop over all events */
    for(Int_t nev=0; nev<num_events; nev++){
//        if(nev%100 == 0)printf("Event: %d\n", nev);
//        TStopwatch timer_getevent;
        ntrigger=0;
        rl->GetEvent(nev);
        stack = rl->Stack();
        
        npart = stack->GetNprimary();
        
        /* Reset all variables before each event is analyzed */
        pt=-1, E=-1, theta=-1, phi=-1, eta=-99, y=-99, eta_calc=-99, p_tot = 0, p_z = 0;
        par_pdg = 0;
        
//        timer_getevent.Stop();
//        fprintf(stdout, "GETEVENT ");
//        timer_getevent.Print("u");
        
//        TStopwatch timer_analysis;
        /* Loop over all particles in stack */
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
            /* Select only the species of hadrons we're interested in (phi, pi, k0, k+-, p) */
            if(TMath::Abs(par_pdg)==310 || TMath::Abs(par_pdg)==2212 || TMath::Abs(par_pdg)==311 || TMath::Abs(par_pdg)==211 || TMath::Abs(par_pdg)==321 || TMath::Abs(par_pdg)==333){         
                point[0] = pt;
                point[1] = phi;
                point[2] = eta_calc;
                hPtPhiEtaDist->Fill(point);
            
                if(TMath::Abs(particle->GetPdgCode()) == 310 && particle->GetNDaughters()==2){
                            TParticle *daughter1 = stack->Particle(particle->GetFirstDaughter());
                            TParticle *daughter2 = stack->Particle(particle->GetFirstDaughter()+1);
                            if(daughter1->GetPdgCode()==211 && daughter2->GetPdgCode()==-211){
                                K0daughterPT->Fill(daughter1->Pt(), daughter2->Pt());
                            }else if(daughter1->GetPdgCode()==-211 && daughter2->GetPdgCode()==211){
                                K0daughterPT->Fill(daughter2->Pt(), daughter1->Pt());
                            }
                        }
                    default:
                        fprintf(stdout, "How did this particle get through? \n");
                        break;
                }
            }
        }

//    TStopwatch timer_write;
//    printf("Total hadrons counted in eta range: %d\n", count);

    rl->~AliRunLoader(); 

    /* Open new file to contain histograms */
    string outputFullPath = output_dir+"/"+output_file;
    TFile *histoOutput = new TFile(outputFullPath.c_str(), "RECREATE");
    histoOutput->cd();

    /* Write histograms to file */
    hPtPhiEtaDist->Write();
    K0daughterPT->Write();
    histoOutput->Close();

//    timer_write.Stop();
//    fprintf(stdout, "WRITE ");
//    timer_write.Print("u");
    timer_total.Stop();
    fprintf(stdout, "TOTAL: ");
    timer_total.Print("u");
}
