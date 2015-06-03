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
    THnSparseF *phiPtPhiEtaDist = new THnSparseF("phiPtPhiEtaDist", "phi meson angular phi & eta distribution", 3, dist_bins, dist_min, dist_max);
    THnSparseF *hPtPhiEtaDist = new THnSparseF("hPtPhiEtaDist", "hadron angular phi distribution", 3, dist_bins, dist_min, dist_max);
    THnSparseF *k0PtPhiEtaDist = new THnSparseF("k0PtPhiEtaDist", "angular phi distribution for k0", 3, dist_bins, dist_min, dist_max);
    THnSparseF *kPtPhiEtaDist = new THnSparseF("kPtPhiEtaDist", "angular phi distribution for charged kaons", 3, dist_bins, dist_min, dist_max);
    THnSparseF *piPtPhiEtaDist = new THnSparseF("piPtPhiEtaDist", "angular phi distribution for charged pions", 3, dist_bins, dist_min, dist_max);
    THnSparseF *pPtPhiEtaDist = new THnSparseF("pPtPhiEtaDist", "angular phi distribution for protons", 3, dist_bins, dist_min, dist_max);
 
    /* TRIGGER PARTICLES PER EVENT */
    TH1I *triggerHist = new TH1I("triggerHist", "Number of Trigger particles per event", 20, 0.5, 20.5);

    /* Bin limits for Delta-phi THnSparse histos */
    Int_t d_phi_bins[3] = {100,100,256};
    Double_t d_phi_min[3] = {0.0, 0.0, -1.57};
    Double_t d_phi_max[3] = {20.0, 20.0, 4.71};


    /* DELTA-PHI HISTOGRAMS (Trigger pT, correlation pT, delta phi)
     * These are the histograms for delta phi correlation between various particle
     * species.  These are what the main analysis is running on.
     */
    THnSparseF *DphiHPhi = new THnSparseF("DphiHPhi", "#Delta#phi correlation for Hadron-Phi", 3, d_phi_bins, d_phi_min, d_phi_max);
    THnSparseF *DphiPiPhi = new THnSparseF("DphiPiPhi", "#Delta#phi correlation for Pion-Phi", 3, d_phi_bins, d_phi_min, d_phi_max);
    THnSparseF *DphiKPhi = new THnSparseF("DphiKPhi", "#Delta#phi correlation for Kaon-Phi", 3, d_phi_bins, d_phi_min, d_phi_max);

    THnSparseF *DphiHK0 = new THnSparseF("DphiHK0", "#Delta#phi correlation for Hadron-K^{0}", 3, d_phi_bins, d_phi_min, d_phi_max);
    THnSparseF *DphiPiK0 = new THnSparseF("DphiPiK0", "#Delta#phi correlation for Pion-K^{0}", 3, d_phi_bins, d_phi_min, d_phi_max);

    THnSparseF *DphiHp = new THnSparseF("DphiHp", "#Delta#phi correlation for Hadron-proton", 3, d_phi_bins, d_phi_min, d_phi_max);
    
    THnSparseF *DphiHK = new THnSparseF("DphiHK", "#Delta#phi correlation for Hadron-kaon", 3, d_phi_bins, d_phi_min, d_phi_max);

    THnSparseF *DphiHPi = new THnSparseF("DphiHpi", "#Delta#phi correlation for Hadron-pion", 3, d_phi_bins, d_phi_min, d_phi_max);
    
    THnSparseF *DphiHH = new THnSparseF("DphiHH", "#Delta#phi correlation for Hadron-Hadron", 3, d_phi_bins, d_phi_min, d_phi_max);

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
    Double_t asso_pt=-1, asso_phi=-1, asso_eta=-9, delta_phi=-99;
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
        asso_pt= -1, asso_phi=-1, asso_eta=-9, delta_phi=-99;
        par_pdg = 0, asso_pdg = 0;
        
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
            if(TMath::Abs(par_pdg)==2212 || TMath::Abs(par_pdg)==311 || TMath::Abs(par_pdg)==211 || TMath::Abs(par_pdg)==321 || TMath::Abs(par_pdg)==333){         
                point[0] = pt;
                point[1] = phi;
                point[2] = eta_calc;
                hPtPhiEtaDist->Fill(point);
                switch(TMath::Abs(par_pdg)){
                    case 2212:
                        pPtPhiEtaDist->Fill(point);
                        break;
                    case 311:
                        k0PtPhiEtaDist->Fill(point);
                        break;
                    case 321:
                        kPtPhiEtaDist->Fill(point);
                        break;
                    case 333:
                        phiPtPhiEtaDist->Fill(point);
                        break;
                    case 211:
                        piPtPhiEtaDist->Fill(point);
                        break;
                    default:
                        fprintf(stdout, "How did this particle get through? \n");
                        break;
                }
                if((TMath::Abs(eta_calc)>0.9)||(pt < 0.150))continue;
                count++;
                /* Only charged hadrons (pions, kaons, protons) should be considered for the trigger particle, everything else discarded */
                if(!((TMath::Abs(par_pdg))==2212||(TMath::Abs(par_pdg))==211||(TMath::Abs(par_pdg))==321)) continue;
                if(pt > 4.0) ntrigger++;
                /* Loop over all particles that aren't this hadron to compute delta-phi */
                for(Int_t apart = 0; apart<npart; apart++){
                    if(apart == part) continue;
                    hAsso = stack->Particle(apart);
                    asso_pt = hAsso->Pt();
                    asso_phi = hAsso->Phi();
                    asso_eta = hAsso->Eta();
                    asso_pdg = hAsso->GetPdgCode();
                    
                    /* Exclude non hadrons (species we're not interested in) */
                    if(!(TMath::Abs(asso_pdg)==2212 || TMath::Abs(asso_pdg)==311 || TMath::Abs(asso_pdg)==211 || TMath::Abs(asso_pdg)==321 || TMath::Abs(asso_pdg)==333)) continue;
                    /* Exclude any hadron that's not within the acceptance */
                    if(asso_pt < 0.15 || TMath::Abs(asso_eta) > 0.9) continue;

                    /* Check that trigger hadron isn't daughter particle of the associated particle*/
                    Int_t numDaughters = hAsso->GetNDaughters();
                    bool isHadronDaughter = false;
                    bool daughterOutOfBounds = false;
                    
                    for(Int_t i=0; i<numDaughters; i++){
                        /* First we check to make sure the trigger hadron isn't a direct daughter of the associated hadron */
                        Int_t daughter = hAsso->GetFirstDaughter()+i;
                        if(daughter == part){
                            isHadronDaughter = true;
                        }
                        /* Next we check to make sure all decay particles of associated hadron are within the acceptance bounds in pT and eta */
                        if(daughter > 0 && daughter < npart){
                            if(stack->Particle(daughter)->Pt() < 0.15 || TMath::Abs(stack->Particle(daughter)->Eta()) > 0.9){
                                daughterOutOfBounds = true;
                            }
                        }
                    }
                    
                    /*  
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
                    */
                    if(daughterOutOfBounds) continue;
                    if(isHadronDaughter) continue;

                    delta_phi = phi - asso_phi;
                    if(delta_phi >3.0*TMath::Pi()/2.0){
                        delta_phi = delta_phi - 2.0*TMath::Pi();
                    }
                    if(delta_phi < -TMath::Pi()/2.0){
                        delta_phi = delta_phi + 2.0*TMath::Pi();
                    }

                    point[0] = pt;
                    point[1] = asso_pt;
                    point[2] = delta_phi;
                    DphiHH->Fill(point);

                    switch(TMath::Abs(asso_pdg)){
                        case 333:
                            DphiHPhi->Fill(point);
                            if(TMath::Abs(par_pdg)==211){
                                DphiPiPhi->Fill(point);
                            }else if(TMath::Abs(par_pdg)==321){
                                DphiKPhi->Fill(point);
                            }
                            break;
                        case 311:
                            DphiHK0->Fill(point);
                            if(TMath::Abs(par_pdg)==211){
                                DphiPiK0->Fill(point);
                            }
                            break;
                        case 211:
                            DphiHPi->Fill(point);
                            break;
                        case 2212:
                            DphiHp->Fill(point);
                            break;
                        case 321:
                            DphiHK->Fill(point);
                            break;
                        default:
                            fprintf(stdout, "Associated particle that shouldn't have gotten through...\n PDG: %d \n", asso_pdg);
                            break;
                    }
                }
            }
        }
    }
//        timer_analysis.Stop();
//        fprintf(stdout, "ANALYSIS ");
//        timer_analysis.Print("u");
    triggerHist->Fill(ntrigger);

//    TStopwatch timer_write;
//    printf("Total hadrons counted in eta range: %d\n", count);

    rl->~AliRunLoader(); 

    /* Open new file to contain histograms */
    string outputFullPath = output_dir+"/"+output_file;
    TFile *histoOutput = new TFile(outputFullPath.c_str(), "RECREATE");
    histoOutput->cd();

    /* Write histograms to file */
    phiPtPhiEtaDist->Write();
    k0PtPhiEtaDist->Write();
    kPtPhiEtaDist->Write();
    piPtPhiEtaDist->Write();
    hPtPhiEtaDist->Write();
    pPtPhiEtaDist->Write();

    DphiHPhi->Write();
    DphiPiPhi->Write();
    DphiKPhi->Write();
    DphiHK0->Write();
    DphiPiK0->Write();
    DphiHK->Write();
    DphiHPi->Write();
    DphiHp->Write();
    DphiHH->Write();

    triggerHist->Write();

    histoOutput->Close();

//    timer_write.Stop();
//    fprintf(stdout, "WRITE ");
//    timer_write.Print("u");
    timer_total.Stop();
    fprintf(stdout, "TOTAL: ");
    timer_total.Print("u");
}
