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
#include "TCanvas.h"
#include "THnSparse.h"
#include "TParticle.h"

#include "AliAnalysisTask.h"
#include "AliAnalysisManager.h"

#include "AliESDEvent.h"
#include "AliESDInputHandler.h"
#include "AliESDtrackCuts.h"
#include "AliAODEvent.h"
#include "AliAODHandler.h"

#include "AliPID.h"
#include "AliESDpid.h"
#include "AliAODPid.h"
#include "AliPIDResponse.h"

#include "AliMCEvent.h"
#include "AliStack.h"

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
fdEdx(0),
fTPCNpts(0),
fTPCnsig(0),
fPDGCodes(0),
fPhiPt(0),
fPhiInvMass(0),
fTPCnsigK(0),
fTPCnsigp(0),
fTPCnsige(0),
fTPCnsigPi(0)
{
    // Constructor
    // Define input and output slots here
    // Input slot #0 works with a TChain
    DefineInput(0, TChain::Class());
    // Output slot #0 id reserved by the base class for AOD
    // Output slot #1 writes into a TH1 container
    DefineOutput(1, TList::Class());
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
fdEdx(0),
fTPCNpts(0),
fTPCnsig(0),
fPDGCodes(0),
fPhiPt(0),
fPhiInvMass(0),
fTPCnsigp(0),
fTPCnsige(0),
fTPCnsigK(0),
fTPCnsigPi(0)
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
    
    fTrkPt = new TH1F("fTrkPt","p_{T} distribution of all tracks;p_{T} (GeV/c);counts",1000,0,100);
    fOutputList->Add(fTrkPt);
    
    fTrketa = new TH1F("fTrketa","All Track #eta distribution;#eta;counts",100,-1.5,1.5);
    fOutputList->Add(fTrketa);
    
    fTrkphi = new TH1F("fTrkphi","All Track #phi distribution;#phi;counts",100,0,6.3);
    fOutputList->Add(fTrkphi);
    
    fdEdx = new TH2F("fdEdx","All Track dE/dx distribution;p (GeV/c);dE/dx",200,0,20,500,0,160);
    fOutputList->Add(fdEdx);
    
    fTPCNpts = new TH2F("fTPCNpts","All track TPC Npoints used for dE/dx calculation;p (GeV/c);N points",200,0,20,200,0.,200.);
    fOutputList->Add(fTPCNpts);
    
    fTPCnsig = new TH2F("fTPCnsig","All Track TPC Nsigma distribution;p (GeV/c);#sigma_{TPC-dE/dx}",1000,0,50,200,-10,10);
    fOutputList->Add(fTPCnsig);
    
    fPDGCodes = new TH1F("fPDGCodes", "PDG codes of tracks", 2000, -1000, 1000);
    fOutputList->Add(fPDGCodes);

    // Additional Histograms for Reconstructed Phi mesons
    fPhiPt = new TH1F("fPhiPt", "p_{T} distribution of all reconstructed #phi mesons; p_{T} (GeV/c);counts", 1000, 0, 100);
    fOutputList->Add(fPhiPt);

    fPhiInvMass = new TH1F("fPhiInvMass", "Invariant mass distribution for all K0 pairs; m (GeV/c^{2}); counts", 1000, 0.5, 2.0);
    fOutputList->Add(fPhiInvMass);

    // Additional TPC Histograms for different particle species
    fTPCnsigK = new TH2F("fTPCnsigK", "Only Kaon TPC Nsigma distribution; p (GeV/c{;#sigma_{TPC-dE/dx}", 1000, 0, 50, 200, -10, 10);
    fOutputList->Add(fTPCnsigK);

    fTPCnsigPi = new TH2F("fTPCnsigPi", "Only Pion TPC Nsigma distribution; p (GeV/c{;#sigma_{TPC-dE/dx}", 1000, 0, 50, 200, -10, 10);
    fOutputList->Add(fTPCnsigPi);

    fTPCnsige = new TH2F("fTPCnsige", "Only Electron TPC Nsigma distribution; p (GeV/c{;#sigma_{TPC-dE/dx}", 1000, 0, 50, 200, -10, 10);
    fOutputList->Add(fTPCnsige);

    fTPCnsigp = new TH2F("fTPCnsigp", "Only Proton TPC Nsigma distribution; p (GeV/c{;#sigma_{TPC-dE/dx}", 1000, 0, 50, 200, -10, 10);
    fOutputList->Add(fTPCnsigp);

    PostData(1,fOutputList);
}

//________________________________________________________________________
void AliAnalysisTaskQA::UserExec(Option_t *)
{
    printf("Made it to UserExec! \n");
    // Main loop
    // Called for each event
    // Post output data.
    
    UInt_t evSelMask=((AliInputEventHandler*)(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler()))->IsEventSelected();
    
    fVevent = dynamic_cast<AliVEvent*>(InputEvent());
    if (!fVevent) {
        printf("ERROR: fVEvent not available\n");
        return;
    }
     
    fESD = dynamic_cast<AliESDEvent*>(InputEvent());
    if (fESD) {
        //   printf("fESD available\n");
        //return;
    }
    ///////////////////////
    // Setting up fStack //
    ///////////////////////
    AliInputEventHandler *fMcHandler = dynamic_cast<AliInputEventHandler*> (AliAnalysisManager::GetAnalysisManager()->GetMCtruthEventHandler());

    AliMCEvent* fMcEvent = 0x0;

    if(fMcHandler){
        fMcEvent = fMcHandler->MCEvent();
        fprintf(stdout, "handler set up! \n");
    }else{
        if(fDebug > 1) printf("AliAnalysisTaskHFDijetMC::UserExec fMcHandler=NULL\n");
        fprintf(stderr, "no MCHandler!\n");
        return;
    }

    if(!fMcEvent){
        if(fDebug > 1) printf("AliAnalysisTaskHFDijetMC::UserExec fMcEvent=NULL \n");
        fprintf(stderr, "No MCEvent!\n");
        return;
    }else{
        fprintf(stdout, "MCEvent Found!\n");
    }

    AliStack *fStack = ((AliMCEvent*)fMcEvent)->Stack();
    if(!fStack){
        fprintf(stderr, "No Stack :(\n");
        return;
    }else{
        fprintf(stdout, "Stack found!\n");
    }

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
    
    
    ///////////////
    //Track loop///
    ///////////////
    for (Int_t iTracks = 0; iTracks < ntracks; iTracks++) {
        
        AliVParticle* Vtrack = 0x0;
        Vtrack  = fVevent->GetTrack(iTracks);
        
        if (!Vtrack) {
            printf("ERROR: Could not receive track %d\n", iTracks);
            continue;
        }
        AliVTrack *track = dynamic_cast<AliVTrack*>(Vtrack);
        AliESDtrack *etrack = dynamic_cast<AliESDtrack*>(Vtrack);
        AliAODTrack *atrack = dynamic_cast<AliAODTrack*>(Vtrack);
        
        ////////////////////
        //Apply track cuts//
        ////////////////////
        if(fAOD)
            if(!atrack->TestFilterMask(AliAODTrack::kTrkGlobalNoDCA)) continue; //mimimum cuts
        
        if(fESD)
            if(!esdTrackCutsH->AcceptTrack(etrack))continue;
        
        
        ////////////////////
        //Track properties//
        ///////////////////
        Double_t dEdx =-999, fTPCnSigma=-999;
        dEdx = track->GetTPCsignal();
        fTPCnSigma = fpidResponse->NumberOfSigmasTPC(track, AliPID::kKaon);
        
        fTrkPt->Fill(track->Pt());
        fTrketa->Fill(track->Eta());
        fTrkphi->Fill(track->Phi());
        fdEdx->Fill(track->P(),dEdx);
        fTPCNpts->Fill(track->P(),track->GetTPCsignalN());
        fTPCnsig->Fill(track->P(),fTPCnSigma);


        //check for PIDs
        Int_t PID = 0;
        PID = track->GetLabel();
        //Int_t PID = 0;
        fPDGCodes->Fill(PID);
        // Using stack to get actual particle PDG codes
        if(PID > 1){ 
            TParticle *MCPart = fStack->Particle(PID);
            if(MCPart){
                Int_t fPDG = MCPart->GetPdgCode();
                if(TMath::Abs(fPDG)==2212)fTPCnsigp->Fill(track->P(), fTPCnSigma);
                if(TMath::Abs(fPDG)==321)fTPCnsigK->Fill(track->P(), fTPCnSigma);
                if(TMath::Abs(fPDG)==211)fTPCnsigPi->Fill(track->P(), fTPCnSigma);
                if(TMath::Abs(fPDG)==11)fTPCnsige->Fill(track->P(), fTPCnSigma);
            }
        }

        ///////////////////////////////
        // Do TPC Cut to identify K+ //
        ///////////////////////////////
        if(fTPCnSigma < 2.0 && track->Pt() > 0.15 && TMath::Abs(track->Eta()) < 0.8 && track->Charge()==1){

            //Do Second Track Loop to find other K0 and calculate their invariant mass
            for(Int_t jTracks = 0; jTracks < ntracks; jTracks++){
                //Exclude double counted particles
                if(jTracks == iTracks) continue;
                
                Double_t fSecondTPCnsig = -999;

                AliVParticle* vSecondTrack = 0x0;
                vSecondTrack = fVevent->GetTrack(jTracks);
                
                if(!vSecondTrack) continue;
                
                AliVTrack *secondTrack = dynamic_cast<AliVTrack*>(vSecondTrack);
                AliESDtrack *eSecondTrack = dynamic_cast<AliESDtrack*>(vSecondTrack);
                AliAODTrack *aSecondTrack = dynamic_cast<AliAODTrack*>(vSecondTrack);

                if(fESD){
                    if(!esdTrackCutsH->AcceptTrack(eSecondTrack))continue;
                }

                if(fAOD){
                   if(!aSecondTrack->TestFilterMask(AliAODTrack::kTrkGlobalNoDCA)) continue; //mimimum cuts
                }
                
                fSecondTPCnsig = fpidResponse->NumberOfSigmasTPC(secondTrack, AliPID::kKaon);

                /////////////////////////////////////////////
                // Do TPC cut to identify second Kaon (K-) //
                ////////////////////////////////////////////
                if(fSecondTPCnsig < 2.0 && secondTrack->Pt() > 0.15 && TMath::Abs(secondTrack->Eta()) < 0.8 && secondTrack->Charge()==-1){

                    Double_t calcPx = 0.0, calcPy = 0.0, calcPz = 0.0;
                    Double_t calcE = 0.0, calcPt = 0.0, calcInvMass = 0.0;

                    calcPx = track->Px()+secondTrack->Px();
                    calcPy = track->Py()+secondTrack->Py();
                    calcPz = track->Pz()+secondTrack->Pz();
                    calcPt = TMath::Sqrt(calcPx*calcPx + calcPy*calcPy);

                    calcE = track->E() + secondTrack->E();

                    fPhiPt->Fill(calcPt);

                    calcInvMass = TMath::Sqrt(calcE*calcE - (calcPx*calcPx + calcPy*calcPy + calcPz*calcPz));

                    fPhiInvMass->Fill(calcInvMass);
                }
            }//inner track loop
        }

    } //track loop
    
    PostData(1, fOutputList);
}      
//________________________________________________________________________
void AliAnalysisTaskQA::Terminate(Option_t *) 
{
    // Draw result to the screen
    // Called once at the end of the query
    
    fOutputList = dynamic_cast<TList*> (GetOutputData(1));
    if (!fOutputList) {
        printf("ERROR: Output list not available\n");
        return;
    }
    
}
