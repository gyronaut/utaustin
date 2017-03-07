/**************************************************************************
 * Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 *                                                                        *
 * Author: The ALICE Off-line Project.                                    *
 * Contributors are mentioned in the code where appropriate.              *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          *
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/

#include "TChain.h"
#include "TTree.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TCanvas.h"
#include "THnSparse.h"
#include "TParticle.h"
#include "TStopwatch.h"

#include "AliAnalysisTask.h"
#include "AliAnalysisManager.h"

#include "AliESDEvent.h"
#include "AliESDInputHandler.h"
#include "AliESDtrackCuts.h"
#include "AliAODEvent.h"
#include "AliAODHandler.h"
#include "AliAODMCParticle.h"

#include "AliPID.h"
#include "AliESDpid.h"
#include "AliAODPid.h"
#include "AliPIDResponse.h"

#include "AliMCEvent.h"
#include "AliStack.h"
#include "AliMCEventHandler.h"

#include "AliAnalysisTaskQA.h"

using std::cout;
using std::endl;

ClassImp(AliAnalysisTaskQA)
//________________________________________________________________________
AliAnalysisTaskQA::AliAnalysisTaskQA(const char *name)
: AliAnalysisTaskSE(name),
fVevent(0),
fESD(0),
fAOD(0),
fpidResponse(0),
fOutputList(0),
fNevents(0),
fVtxZ(0),
fVtxX(0),
fVtxY(0),
fTrigMulti(0),
fTrkPt(0),
fTrketa(0),
fTrkphi(0),
fHybridTrkPt(0),
fHybridTrketa(0),
fHybridTrkphi(0),
fHybridGlobalTrkPt(0),
fHybridGlobalTrketa(0),
fHybridGlobalTrkphi(0),
fdEdx(0),
fTPCNpts(0),
fKKUSDist(0),
fKKLSDist(0)
{
    // Constructor
    // Define input and output slots here
    // Input slot #0 works with a TChain
    DefineInput(0, TChain::Class());
    // Output slot #0 id reserved by the base class for AOD
    // Output slot #1 writes into a TH1 container
    DefineOutput(1, TList::Class());
    //printf("\n\n!!!!!!!!!!]\n done with the constructor! \n");
    //fflush(stdout);

}
//________________________________________________________________________
AliAnalysisTaskQA::AliAnalysisTaskQA()
: AliAnalysisTaskSE("DefaultTask_HfeEMCQA"),
fVevent(0),
fESD(0),
fAOD(0),
fpidResponse(0),
fOutputList(0),
fNevents(0),
fVtxZ(0),
fVtxX(0),
fVtxY(0),
fTrigMulti(0),
fTrkPt(0),
fTrketa(0),
fTrkphi(0),
fHybridTrkPt(0),
fHybridTrketa(0),
fHybridTrkphi(0),
fHybridGlobalTrkPt(0),
fHybridGlobalTrketa(0),
fHybridGlobalTrkphi(0),
fdEdx(0),
fTPCNpts(0),
fKKUSDist(0),
fKKLSDist(0)
{
    //Default constructor
    // Define input and output slots here
    // Input slot #0 works with a TChain
    DefineInput(0, TChain::Class());
    // Output slot #0 id reserved by the base class for AOD
    // Output slot #1 writes into a TH1 container
    // DefineOutput(1, TH1I::Class());
    DefineOutput(1, TList::Class());
    //DefineOutput(3, TTree::Class());
}
//________________________________________________________________________
AliAnalysisTaskQA::~AliAnalysisTaskQA()
{
    //Destructor
    delete fOutputList;
}
//________________________________________________________________________
void AliAnalysisTaskQA::UserCreateOutputObjects()
{
    //printf("\n!!!!!\n Starting UserCreateOutputObjects \n\n");
    //fflush(stdout);
    // Create histograms
    // Called once
    AliDebug(3, "Creating Output Objects");
    
    /////////////////////////////////////////////////
    //Automatic determination of the analysis mode//
    ////////////////////////////////////////////////
    AliVEventHandler *inputHandler = dynamic_cast<AliVEventHandler *>(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler());
    if(!TString(inputHandler->IsA()->GetName()).CompareTo("AliAODInputHandler")){
        SetAODAnalysis();
    } else {
        SetESDAnalysis();
    }
    printf("Analysis Mode: %s Analysis\n", IsAODanalysis() ? "AOD" : "ESD");
    
    ////////////////
    //Output list//
    ///////////////
    fOutputList = new TList();
    fOutputList->SetOwner();
    
    fNevents = new TH1F("fNevents","No of events",3,-0.5,2.5);
    fOutputList->Add(fNevents);
    fNevents->GetYaxis()->SetTitle("counts");
    fNevents->GetXaxis()->SetBinLabel(1,"All");
    fNevents->GetXaxis()->SetBinLabel(2,"With >2 Trks");
    fNevents->GetXaxis()->SetBinLabel(3,"Vtx_{z}<10cm");
    
    fVtxZ = new TH1F("fVtxZ","Z vertex position;Vtx_{z};counts",1000,-50,50);
    fOutputList->Add(fVtxZ);
    
    fVtxY = new TH1F("fVtxY","Y vertex position;Vtx_{y};counts",1000,-50,50);
    fOutputList->Add(fVtxY);
    
    fVtxX = new TH1F("fVtxX","X vertex position;Vtx_{x};counts",1000,-50,50);
    fOutputList->Add(fVtxX);
    
    fTrigMulti = new TH2F("fTrigMulti","Multiplicity distribution for different triggers; Trigger type; multiplicity",11,-1,10,2000,0,2000);
    fOutputList->Add(fTrigMulti);
    
    //Global track histos
    fTrkPt = new TH1F("fTrkPt","p_{T} distribution of all tracks;p_{T} (GeV/c);counts",1000,0,100);
    fOutputList->Add(fTrkPt);
    
    fTrketa = new TH1F("fTrketa","All Track #eta distribution;#eta;counts",100,-1.5,1.5);
    fOutputList->Add(fTrketa);
    
    fTrkphi = new TH1F("fTrkphi","All Track #phi distribution;#phi;counts",100,0,6.3);
    fOutputList->Add(fTrkphi);

    //HybridTPC track histos
    fHybridTrkPt = new TH1F("fHybridTrkPt","p_{T} distribution of all hybrid tracks;p_{T} (GeV/c);counts",1000,0,100);
    fOutputList->Add(fHybridTrkPt);
    
    fHybridTrketa = new TH1F("fHybridTrketa","All Hybrid Track #eta distribution;#eta;counts",100,-1.5,1.5);
    fOutputList->Add(fHybridTrketa);
    
    fHybridTrkphi = new TH1F("fHybridTrkphi","All Hybrid Track #phi distribution;#phi;counts",100,0,6.3);
    fOutputList->Add(fHybridTrkphi);
    
    //HybridGlobal track histos
    fHybridGlobalTrkPt = new TH1F("fHybridGlobalTrkPt","p_{T} distribution of all hybrid tracks;p_{T} (GeV/c);counts",1000,0,100);
    fOutputList->Add(fHybridGlobalTrkPt);
    
    fHybridGlobalTrketa = new TH1F("fHybridGlobalTrketa","All HybridGlobal Track #eta distribution;#eta;counts",100,-1.5,1.5);
    fOutputList->Add(fHybridGlobalTrketa);
    
    fHybridGlobalTrkphi = new TH1F("fHybridGlobalTrkphi","All HybridGlobal Track #phi distribution;#phi;counts",100,0,6.3);
    fOutputList->Add(fHybridGlobalTrkphi);
 
    fdEdx = new TH2F("fdEdx","All Track dE/dx distribution;p (GeV/c);dE/dx",200,0,20,500,0,160);
    fOutputList->Add(fdEdx);
    
    fTPCNpts = new TH2F("fTPCNpts","All track TPC Npoints used for dE/dx calculation;p (GeV/c);N points",200,0,20,200,0.,200.);
    fOutputList->Add(fTPCNpts);

    //Histograms for US and LS Kaon pairs for PID studies:
    Int_t bins[6] = {50, 50, 50, 50, 50, 100}; //K1 nsigTPC, K2 nsigTPC, K1 nsigTOF, K2 nsigTOF, KK pT, KK invmass
    Double_t min[6] = {-5, -5, -5, -5, 0.1, 0.98};
    Double_t max[6] = {5, 5, 5, 5, 10.1, 1.1};
 
    fKKUSDist = new THnSparseF("fkkUSDist", "Distribution for all US Kaon pairs", 6, bins, min, max);
    fOutputList->Add(fKKUSDist);

    fKKLSDist = new THnSparseF("fkkLSDist", "Distribution for all LS Kaon pairs", 6, bins, min, max);
    fOutputList->Add(fKKLSDist);
    PostData(1,fOutputList);
}

//________________________________________________________________________
void AliAnalysisTaskQA::UserExec(Option_t *)
{
    //printf("\n!!!!!!!!!!!!!!\nMade it to UserExec! \n");
    //fflush(stdout);
    // Main loop
    // Called for each event
    // Post output data.
    
//    TStopwatch *initTimer = new TStopwatch();
//    initTimer->Start();

    UInt_t evSelMask=((AliInputEventHandler*)(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler()))->IsEventSelected();
    
    fVevent = dynamic_cast<AliVEvent*>(InputEvent());
    if (!fVevent) {
        printf("ERROR: fVEvent not available\n");
        return;
    }
     
    fESD = dynamic_cast<AliESDEvent*>(InputEvent());

    ////////////////////
    //cuts initialised//
    ///////////////////
    AliESDtrackCuts* esdTrackCutsH = AliESDtrackCuts::GetStandardITSTPCTrackCuts2011(kFALSE);
    esdTrackCutsH->SetMaxDCAToVertexXY(2.4);
    esdTrackCutsH->SetMaxDCAToVertexZ(3.2);
    esdTrackCutsH->SetDCAToVertex2D(kTRUE);

    fAOD = dynamic_cast<AliAODEvent*>(InputEvent());
    if (fAOD) {
        // printf("fAOD available\n");
        //return;
    }

    ///////////////////
    //PID initialised//
    //////////////////
    fpidResponse = fInputHandler->GetPIDResponse();

    ////////////////
    //Event vertex//
    ///////////////
    Int_t ntracks = -999;
    ntracks = fVevent->GetNumberOfTracks();
    if(ntracks < 1) printf("There are %d tracks in this event\n",ntracks);

    fNevents->Fill(0); //all events
    Double_t Zvertex = -100, Xvertex = -100, Yvertex = -100;
    const AliVVertex *pVtx = fVevent->GetPrimaryVertex();
    Double_t NcontV = pVtx->GetNContributors();
    if(NcontV<2)return;
    fNevents->Fill(1); //events with 2 tracks

    Zvertex = pVtx->GetZ();
    Yvertex = pVtx->GetY();
    Xvertex = pVtx->GetX();
    fVtxZ->Fill(Zvertex);
    fVtxX->Fill(Xvertex);
    fVtxY->Fill(Yvertex);

    /////////////////
    //trigger check//
    /////////////////
    fVevent->GetFiredTriggerClasses();

    Int_t trigger = -1;
    if (fAOD){
        //Double_t multiplicity=fAOD->GetHeader()->GetRefMultiplicity();
        AliAODHeader *header = dynamic_cast<AliAODHeader*>(fAOD->GetHeader());
        if(!header) AliFatal("Not a standard AOD");
        Double_t multiplicity = header->GetRefMultiplicity();
        
        fTrigMulti->Fill(-0.5, multiplicity);
        if(evSelMask & AliVEvent::kAny) fTrigMulti->Fill(0.5, multiplicity);
        if(evSelMask & AliVEvent::kMB) fTrigMulti->Fill(1.5, multiplicity);
        if(evSelMask & AliVEvent::kINT7) fTrigMulti->Fill(2.5, multiplicity);
        if(evSelMask & AliVEvent::kINT8) fTrigMulti->Fill(3.5, multiplicity);
        if(evSelMask & AliVEvent::kEMC1) fTrigMulti->Fill(4.5, multiplicity);
        if(evSelMask & AliVEvent::kEMC7) fTrigMulti->Fill(5.5, multiplicity);
        if(evSelMask & AliVEvent::kEMC8) fTrigMulti->Fill(6.5, multiplicity);
        if(evSelMask & AliVEvent::kEMCEJE) fTrigMulti->Fill(7.5, multiplicity);
        if(evSelMask & AliVEvent::kEMCEGA) fTrigMulti->Fill(8.5, multiplicity);
    }
    
    ////////////////////
    //event selection//
    ///////////////////
    if(fabs(Zvertex>10.0))return;
    fNevents->Fill(2); //events after z vtx cut

    //Initialize the vectors/points that will be used to fill the histograms
    std::vector<TLorentzVector> phiCandidates;
    std::vector<Int_t> phiDaughterTrackNum;
//    std::vector<TLorentzVector> phiReals;
    std::vector<TLorentzVector> phiLikesignCandidates;
    std::vector<Int_t> likesignDaughterTrackNum;
    std::vector<TLorentzVector> KCandidates;
    TLorentzVector phi;
    TLorentzVector K;

    Double_t distPoint[6] = {0, 0, 0, 0, 0, 0}; //K+ TPCnsigma, K- TPCnsigma, K+ TOFnsigma, K- TOFnsigma, KK pT, KK invmass

//    TParticle *MCFirstDecay = 0x0;
//    AliAODMCParticle* MCFirstDecayTrack = 0x0;
    AliVTrack *firstDecayTrack = 0x0;
    AliESDtrack *eFirstDecayTrack = 0x0;
    AliAODTrack *aFirstDecayTrack = 0x0;
    AliVParticle *vFirstDecayTrack = 0x0;

    /* First Loop - Filling a vector with Lorentz vectors for all
     * Phi Candidates/real Phi's
     */
    for(Int_t i_track = 0; i_track < ntracks; i_track++){
        vFirstDecayTrack = 0x0;
        vFirstDecayTrack = fVevent->GetTrack(i_track);

        if(!vFirstDecayTrack){
            printf("Error: Could not receive track %d\n", i_track);
            continue;
        }
        firstDecayTrack = dynamic_cast<AliVTrack*>(vFirstDecayTrack);
        eFirstDecayTrack = dynamic_cast<AliESDtrack*>(vFirstDecayTrack);
        aFirstDecayTrack = dynamic_cast<AliAODTrack*>(vFirstDecayTrack);

        if(fAOD)
            if(!aFirstDecayTrack->TestFilterMask(AliAODTrack::kTrkGlobalNoDCA)) continue; //mimimum cuts

        if(fESD)
            if(!esdTrackCutsH->AcceptTrack(eFirstDecayTrack))continue;

        if(firstDecayTrack->Pt() > 0.15 && TMath::Abs(firstDecayTrack->Eta()) < 0.8){
            Double_t fposTPCnSigma = -999, fnegTPCnSigma = -999;
            Double_t fposTOFnSigma = -999, fnegTOFnSigma = -999;
            //check for labels
            Int_t label = 0;
            label = firstDecayTrack->GetLabel();

            fposTPCnSigma = fpidResponse->NumberOfSigmasTPC(firstDecayTrack, AliPID::kKaon);
            fposTOFnSigma = fpidResponse->GetNumberOfSigmasTOF(firstDecayTrack, AliPID::kKaon);
            AliVTrack *secondDecayTrack = 0x0;
            AliESDtrack *eSecondDecayTrack = 0x0;
            AliAODTrack *aSecondDecayTrack = 0x0;
            AliVParticle *vSecondDecayTrack = 0x0;

            for(Int_t j_track = 0; j_track < ntracks; j_track++){
                if(i_track == j_track) continue;
                vSecondDecayTrack = 0x0;
                vSecondDecayTrack = fVevent->GetTrack(j_track);

                if(!vSecondDecayTrack){
                    printf("Error: Could not receive track %d\n", j_track);
                    continue;
                }
                secondDecayTrack = dynamic_cast<AliVTrack*>(vSecondDecayTrack);
                eSecondDecayTrack = dynamic_cast<AliESDtrack*>(vSecondDecayTrack);
                aSecondDecayTrack = dynamic_cast<AliAODTrack*>(vSecondDecayTrack);

                if(fAOD)
                    if(!aSecondDecayTrack->TestFilterMask(AliAODTrack::kTrkGlobalNoDCA)) continue; //mimimum cuts

                if(fESD)
                    if(!esdTrackCutsH->AcceptTrack(eSecondDecayTrack))continue;

                if(secondDecayTrack->Pt() > 0.15 && TMath::Abs(secondDecayTrack->Eta()) < 0.8){
                    fnegTPCnSigma = -999;
                    fnegTOFnSigma = -999;
                    fnegTPCnSigma = fpidResponse->NumberOfSigmasTPC(secondDecayTrack, AliPID::kKaon);
                    fnegTOFnSigma = fpidResponse->GetNumberOfSigmasTOF(secondDecayTrack, AliPID::kKaon);
                    Double_t calcPx = 0.0, calcPy = 0.0, calcPz = 0.0;
                    Double_t calcE = 0.0, calcPt = 0.0, calcInvMass = 0.0, calcE1=0.0, calcE2=0.0;
                    calcPx = firstDecayTrack->Px()+secondDecayTrack->Px();
                    calcPy = firstDecayTrack->Py()+secondDecayTrack->Py();
                    calcPz = firstDecayTrack->Pz()+secondDecayTrack->Pz();
                    calcPt = TMath::Sqrt(calcPx*calcPx + calcPy*calcPy);

                    //calcE = firstDecayTrack->E() + secondDecayTrack->E();
                    //calculating energy of each kaon candidate assuming k mass
                    calcE1 = TMath::Sqrt(0.4937*0.4937 + firstDecayTrack->P()*firstDecayTrack->P());
                    calcE2 = TMath::Sqrt(0.4937*0.4937 + secondDecayTrack->P()*secondDecayTrack->P());
                    calcE = calcE1 + calcE2;
                    calcInvMass = TMath::Sqrt(calcE*calcE - (calcPx*calcPx + calcPy*calcPy + calcPz*calcPz));

                    //Unlike sign pairs - create phi inv-mass distribution (and add to vector if invmass between 0.98 and 1.1)
                    if(firstDecayTrack->Charge() == 1 && secondDecayTrack->Charge() == -1){
                        if(calcInvMass >= 0.98 && calcInvMass <= 1.1){
                            //Set-up TLorenztVector (px, py, pz, E), then push to vector
                            phi.SetPx(firstDecayTrack->Px()+secondDecayTrack->Px());
                            phi.SetPy(firstDecayTrack->Py()+secondDecayTrack->Py());
                            phi.SetPz(firstDecayTrack->Pz()+secondDecayTrack->Pz());
                            phi.SetE(calcE);

                            distPoint[0] = fposTPCnSigma;
                            distPoint[1] = fnegTPCnSigma;
                            distPoint[2] = fposTOFnSigma;
                            distPoint[3] = fnegTOFnSigma;
                            distPoint[4] = phi.Pt();
                            distPoint[5] = calcInvMass;
                            fKKUSDist->Fill(distPoint);
                        }
                    }else if(firstDecayTrack->Charge()*secondDecayTrack->Charge() == 1){
                        if(calcInvMass >= 0.98 && calcInvMass <=1.1){
                            phi.SetPx(firstDecayTrack->Px()+secondDecayTrack->Px());
                            phi.SetPy(firstDecayTrack->Py()+secondDecayTrack->Py());
                            phi.SetPz(firstDecayTrack->Pz()+secondDecayTrack->Pz());
                            phi.SetE(calcE);

                            distPoint[0] = fposTPCnSigma;
                            distPoint[1] = fnegTPCnSigma;
                            distPoint[2] = fposTOFnSigma;
                            distPoint[3] = fnegTOFnSigma;                                  
                            distPoint[4] = phi.Pt();
                            distPoint[5] = calcInvMass;
                            fKKLSDist->Fill(distPoint);
                        }

                    }
                }
            }
        }
    } 

    //////////////////
    // Trigger loop // 
    //////////////////
    /*
    TParticle* MCTriggerParticle = 0x0;
    AliAODMCParticle* MCTriggerTrack = 0x0;
    */
    AliVTrack *triggerTrack = 0x0;
    AliESDtrack *etriggerTrack = 0x0;
    AliAODTrack *atriggerTrack = 0x0;
    AliVParticle* VtriggerTrack = 0x0;
    
    ////////////////////////////////////////////
    // Second loop, building d-phi histograms //
    ////////////////////////////////////////////
    for (Int_t i_track = 0; i_track < ntracks; i_track++) {
        
        VtriggerTrack = 0x0;
        VtriggerTrack  = fVevent->GetTrack(i_track);
        
        if (!VtriggerTrack) {
            printf("ERROR: Could not receive track %d\n", i_track);
            continue;
        }
        triggerTrack = dynamic_cast<AliVTrack*>(VtriggerTrack);
        etriggerTrack = dynamic_cast<AliESDtrack*>(VtriggerTrack);
        atriggerTrack = dynamic_cast<AliAODTrack*>(VtriggerTrack);
        
        //fill hybrid track histos if the track is hybridTPC
        if(triggerTrack->Pt() > 0.15 && TMath::Abs(triggerTrack->Eta()) < 0.8 && atriggerTrack->IsHybridTPCConstrainedGlobal()){
            fHybridTrkPt->Fill(triggerTrack->Pt());
            fHybridTrketa->Fill(triggerTrack->Eta());
            fHybridTrkphi->Fill(triggerTrack->Phi());
        }

        //fill hybrid track histos if the track is hybridGlobal
        if(triggerTrack->Pt() > 0.15 && TMath::Abs(triggerTrack->Eta()) < 0.8 && atriggerTrack->IsHybridGlobalConstrainedGlobal()){
            fHybridGlobalTrkPt->Fill(triggerTrack->Pt());
            fHybridGlobalTrketa->Fill(triggerTrack->Eta());
            fHybridGlobalTrkphi->Fill(triggerTrack->Phi());
        }

        //fill global track histos if the track is global
        if( triggerTrack->Pt() > 0.15 && TMath::Abs(triggerTrack->Eta()) < 0.8 && atriggerTrack->TestFilterMask(AliAODTrack::kTrkGlobalNoDCA)){
            fTrkPt->Fill(triggerTrack->Pt());
            fTrketa->Fill(triggerTrack->Eta());
            fTrkphi->Fill(triggerTrack->Phi());
        }



        ////////////////////
        //Apply track cuts//
        ////////////////////
        if(fAOD)
            if(!atriggerTrack->IsHybridGlobalConstrainedGlobal()) continue; //selecting just hybrid-global tracks for trigger, continue otherwise
        
        if(fESD)
            if(!esdTrackCutsH->AcceptTrack(etriggerTrack))continue;
        
    } //track loop
//    analysisTimer->Stop();
//    printf("ANALYSIS: ");
//    analysisTimer->Print("u");
//    printf("\n");  

//    TStopwatch *writeTimer = new TStopwatch();
//    writeTimer->Start();
    PostData(1, fOutputList);
//    printf("WRITE: ");
//    writeTimer->Print("u");
//    printf("\n");

//    delete initTimer;
//    delete analysisTimer;
//    delete writeTimer;
}      
//________________________________________________________________________
void AliAnalysisTaskQA::Terminate(Option_t *) 
{
    // Draw result to the screen
    // Called once at the end of the query
    printf("terminating task... \n");
    fOutputList = dynamic_cast<TList*> (GetOutputData(1));
    if (!fOutputList) {
        printf("ERROR: Output list not available\n");
        return;
    }
    
}
