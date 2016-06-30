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
fdEdx(0),
fTPCNpts(0),
fTPCKaonNSig(0),
fPDGCodes(0),
fPhiInvMass(0),
fTruthPhiInvMass(0),
fTruthTracksPhiInvMass(0),
fPhiLikeSignInvMass(0),
fKInvMass(0),
fKLikeSignInvMass(0),
fTPCKaonNSigp(0),
fTPCKaonNSige(0),
fTPCKaonNSigK(0),
fTPCKaonNSigPi(0),
fPhiDaughterPTKept(0),
fPhiDaughterPTCut(0),
fKstarDaughterPTKept(0),
fKstarDaughterPTCut(0),
fK0DaughterPTKept(0),
fK0DaughterPTCut(0),
fDphiHPhi(0),
fDphiHK(0),
fDphiHKstar(0),
fDphiHK0(0),
fDphiHPi(0),
fDphiHp(0)
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
fdEdx(0),
fTPCNpts(0),
fTPCKaonNSig(0),
fPDGCodes(0),
fPhiInvMass(0),
fTruthPhiInvMass(0),
fTruthTracksPhiInvMass(0),
fPhiLikeSignInvMass(0),
fKInvMass(0),
fKLikeSignInvMass(0),
fTPCKaonNSigp(0),
fTPCKaonNSige(0),
fTPCKaonNSigK(0),
fTPCKaonNSigPi(0),
fPhiDaughterPTKept(0),
fPhiDaughterPTCut(0),
fKstarDaughterPTKept(0),
fKstarDaughterPTCut(0),
fK0DaughterPTKept(0),
fK0DaughterPTCut(0),
fDphiHPhi(0),
fDphiHK(0),
fDphiHKstar(0),
fDphiHK0(0),
fDphiHPi(0),
fDphiHp(0)
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
    
    // TPC NSig histogram 
    fTPCKaonNSig = new TH2F("fTPCKaonNSig","All Track TPC Nsigma distribution;p (GeV/c);#sigma_{TPC-dE/dx}",1000,0,50,200,-10,10);
    fOutputList->Add(fTPCKaonNSig);
    
    fPDGCodes = new TH1F("fPDGCodes", "PDG codes of tracks", 2000, -1000, 1000);
    fOutputList->Add(fPDGCodes);

    // Additional Histograms for Reconstructed Phi mesons
    Int_t bins[2] = {100, 1000};
    Double_t min[2] = {0.0, 0.5};
    Double_t max[2] = {10.0, 2.0};
 
    fPhiInvMass = new THnSparseF("fPhiInvMass", "Invariant mass distribution for all K+- pairs per p_{T}", 2, bins, min, max);
    fOutputList->Add(fPhiInvMass);

    fTruthPhiInvMass = new THnSparseF("fTruthPhiInvMass", "Invariant mass distribution for true K+- that come from #phi per p_{T}", 2, bins, min, max);
    fOutputList->Add(fTruthPhiInvMass);

    fTruthTracksPhiInvMass = new THnSparseF("fTruthTracksPhiInvMass", "Invariant mass distribution of reconstructed tracks known to come from #phi per p_{T}", 2, bins, min, max);
    fOutputList->Add(fTruthTracksPhiInvMass);

    fKInvMass = new THnSparseF("fKInvMass", "Invariant mass distribution for all unlike sign Kaon/Pion pairs per p_{T}", 2, bins, min, max);
    fOutputList->Add(fKInvMass);

    Int_t lsbins[3] = {100, 1000, 3};
    Double_t lsmin[3] = {0.0, 0.5, -1.1};
    Double_t lsmax[3] = {10.0, 2.0, 1.1};

    fPhiLikeSignInvMass = new THnSparseF("fPhiLikeSignInvMass", "Invariant mass distribution of like-sign Kaon pairs per p_{T}", 3, lsbins, lsmin, lsmax);
    fOutputList->Add(fPhiLikeSignInvMass);

    fKLikeSignInvMass = new THnSparseF("fKLikeSignInvMass", "Invariant mass distribution of like-sign Kaon/Pion pairs per p_{T}", 3, lsbins, lsmin, lsmax);
    fOutputList->Add(fKLikeSignInvMass);
   
    // Additional TPC Histograms for different particle species (in relation to Kaon)
    fTPCKaonNSigK = new TH2F("fTPCKaonNSigK", "Only Kaon TPC Nsigma distribution; p (GeV/c); #sigma_{TPC-dE/dx}", 1000, 0, 50, 200, -10, 10);
    fOutputList->Add(fTPCKaonNSigK);

    fTPCKaonNSigPi = new TH2F("fTPCKaonNSigPi", "Only Pion TPC Nsigma distribution; p (GeV/c); #sigma_{TPC-dE/dx}", 1000, 0, 50, 200, -10, 10);
    fOutputList->Add(fTPCKaonNSigPi);

    fTPCKaonNSige = new TH2F("fTPCKaonNSige", "Only Electron TPC Nsigma distribution; p (GeV/c); #sigma_{TPC-dE/dx}", 1000, 0, 50, 200, -10, 10);
    fOutputList->Add(fTPCKaonNSige);

    fTPCKaonNSigp = new TH2F("fTPCKaonNSigp", "Only Proton TPC Nsigma distribution; p (GeV/c); #sigma_{TPC-dE/dx}", 1000, 0, 50, 200, -10, 10);
    fOutputList->Add(fTPCKaonNSigp);

    fPhiDaughterPTKept = new TH2F("fPhiDaughterPTKept", "Pt distribution of Phi Meson Daughter Particles, Kept; p_{T}^{K+} (GeV/c); p_{T}^{K-} (GeV/c)", 16,0,4,16,0,4);
    fOutputList->Add(fPhiDaughterPTKept);

    fPhiDaughterPTCut = new TH2F("fPhiDaughterPTCut", "Pt distribution of Phi Meson Daughter Particles, Cut; p_{T}^{K+} (GeV/c); p_{T}^{K-} (GeV/c)", 40,0,10,40,0,10);
    fOutputList->Add(fPhiDaughterPTCut);
 
    fKstarDaughterPTKept = new TH2F("fKstarDaughterPTKept", "Pt distribution of K*(892) Meson Daughter Particles, Kept; p_{T}^{K} (GeV/c); p_{T}^{#pi} (GeV/c)", 16,0,4,16,0,4);
    fOutputList->Add(fKstarDaughterPTKept);

    fKstarDaughterPTCut = new TH2F("fKstarDaughterPTCut", "Pt distribution of K*(892) Meson Daughter Particles, Cut; p_{T}^{K} (GeV/c); p_{T}^{#pi} (GeV/c)", 40,0,10,40,0,10);
    fOutputList->Add(fKstarDaughterPTCut);
 
    fK0DaughterPTKept = new TH2F("fK0DaughterPTKept", "Pt distribution of K0 Meson Daughter Particles, Kept; p_{T}^{#pi^{+}} (GeV/c); p_{T}^{#pi^{-}} (GeV/c)", 16,0,4,16,0,4);
    fOutputList->Add(fK0DaughterPTKept);

    fK0DaughterPTCut = new TH2F("fK0DaughterPTCut", "Pt distribution of K0 Meson Daughter Particles, Cut; p_{T}^{#pi^{+}} (GeV/c); p_{T}^{#pi^{-}} (GeV/c)", 40,0,10,40,0,10);
    fOutputList->Add(fK0DaughterPTCut);
   
    // Delta-phi histograms for different hadron-particle correlations (trigger pT, correlation pT, delta-phi)
    Int_t dphi_bins[3]= {100, 100, 256};
    Double_t dphi_min[3] = {0.0, 0.0, -1.57};
    Double_t dphi_max[3] = {20.0, 20.0, 4.71};

    fDphiHPhi = new THnSparseF("fDphiHPhi", "Hadron-#Phi #Delta#phi correlations", 3, dphi_bins, dphi_min, dphi_max);
    fOutputList->Add(fDphiHPhi);

    fDphiHKstar = new THnSparseF("fDphiHKstar", "Hadron-K*(892) #Delta#phi correlations", 3, dphi_bins, dphi_min, dphi_max);
    fOutputList->Add(fDphiHKstar);

    fDphiHK = new THnSparseF("fDphiHK", "Hadron-Kaon #Delta#phi correlations", 3, dphi_bins, dphi_min, dphi_max);
    fOutputList->Add(fDphiHK);

    fDphiHK0 = new THnSparseF("fDphiHK0", "Hadron-K0 #Delta#phi correlations", 3, dphi_bins, dphi_min, dphi_max);
    fOutputList->Add(fDphiHK0);

    fDphiHPi = new THnSparseF("fDphiHPi", "Hadron-Pi #Delta#phi correlations", 3, dphi_bins, dphi_min, dphi_max);
    fOutputList->Add(fDphiHPi);

    fDphiHp = new THnSparseF("fDphiHp", "Hadron-proton #Delta#phi correlations", 3, dphi_bins, dphi_min, dphi_max);
    fOutputList->Add(fDphiHp);

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

    /////////////////////////////////
    // Setting up fStack (MC Only) //
    /////////////////////////////////
/*    
    AliMCEventHandler *fMCHandler = 0x0;
    AliMCEvent *fMcEvent = 0x0;
    AliStack *fStack = 0x0;
    TClonesArray *mcArray = 0x0;
    
    if(fESD){
        fMCHandler = dynamic_cast<AliMCEventHandler*>(AliAnalysisManager::GetAnalysisManager()->GetMCtruthEventHandler()); 
   
        
        if(fMCHandler){
            fMcEvent = fMCHandler->MCEvent();
            // printf("handler set up! \n");
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
            //fprintf(stdout, "MCEvent Found!\n");
        }

        fStack = ((AliMCEvent*)fMcEvent)->Stack();
        if(!fStack){
            fprintf(stderr, "No Stack :(\n");
            return;
        }else{
            //printf("Stack found!\n");
        }
    }else if(fAOD){
        mcArray = dynamic_cast<TClonesArray*>(fAOD->FindListObject(AliAODMCParticle::StdBranchName()));
        if(!mcArray){
            AliError("Array of MC particles not found");
            return;
        }
    } 
*/
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
//    std::vector<TLorentzVector> phiReals;
    std::vector<TLorentzVector> KCandidates;
//    std::vector<TLorentzVector> KReals;
    std::vector<TLorentzVector> K0Reals;
    TLorentzVector phi;
    TLorentzVector K;

    Double_t point[2] = {0, 0};
    Double_t lspoint[3] = {0, 0, 0};
    Double_t dphi_point[3] = {0, 0, 0};

//    TParticle *MCFirstDecay = 0x0;
//    AliAODMCParticle* MCFirstDecayTrack = 0x0;
    AliVTrack *firstDecayTrack = 0x0;
    AliESDtrack *eFirstDecayTrack = 0x0;
    AliAODTrack *aFirstDecayTrack = 0x0;
    AliVParticle *vFirstDecayTrack = 0x0;

    /* First Loop - Filling a vector with Lorentz vectors for all
     * Phi Candidates/real Phi's, as well as K*(892) candidates/reals
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
            Double_t fTPCnSigma = -999;
            Double_t fpiTPCnSigma = -999;
            //check for labels
            Int_t label = 0;
            label = firstDecayTrack->GetLabel();

            // Using stack to get actual particle PDG codes (MC Only)
            /* 
            Int_t fPDG = 0;
            Int_t motherPDG = 0;
            Int_t motherIndex = 0;
            if(label > 1){             
                if(fESD){
                    MCFirstDecay = fStack->Particle(label);
                    if(MCFirstDecay){
                        fPDG = MCFirstDecay->GetPdgCode();
                        if(fPDG > 0){
                            motherIndex = MCFirstDecay->GetMother(0);
                            if(motherIndex > 0){
                                motherPDG = fStack->Particle(motherIndex)->GetPdgCode();
                            }
                        }
                    }
                }else if(fAOD){
                    MCFirstDecayTrack = (AliAODMCParticle*)mcArray->At(label);
                    if(MCFirstDecayTrack){
                        fPDG = MCFirstDecayTrack->GetPdgCode();
                    }
                }
            }
            */

            fTPCnSigma = fpidResponse->NumberOfSigmasTPC(firstDecayTrack, AliPID::kKaon);

//            TParticle *MCSecondDecay = 0x0;
//            AliAODMCParticle* MCSecondDecayTrack = 0x0;
            AliVTrack *secondDecayTrack = 0x0;
            AliESDtrack *eSecondDecayTrack = 0x0;
            AliAODTrack *aSecondDecayTrack = 0x0;
            AliVParticle *vSecondDecayTrack = 0x0;

            //Cut on kaon candidates
            if(TMath::Abs(fTPCnSigma) < 2.0){
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
                        fTPCnSigma = -999;
                      
                        fTPCnSigma = fpidResponse->NumberOfSigmasTPC(secondDecayTrack, AliPID::kKaon);
                        fpiTPCnSigma = fpidResponse->NumberOfSigmasTPC(secondDecayTrack, AliPID::kPion);
                        Double_t calcPx = 0.0, calcPy = 0.0, calcPz = 0.0;
                        Double_t calcE = 0.0, calcPt = 0.0, calcInvMass = 0.0;
                        if(TMath::Abs(fTPCnSigma) < 2.0){
                            calcPx = firstDecayTrack->Px()+secondDecayTrack->Px();
                            calcPy = firstDecayTrack->Py()+secondDecayTrack->Py();
                            calcPz = firstDecayTrack->Pz()+secondDecayTrack->Pz();
                            calcPt = TMath::Sqrt(calcPx*calcPx + calcPy*calcPy);

                            calcE = firstDecayTrack->E() + secondDecayTrack->E();
                            calcInvMass = TMath::Sqrt(calcE*calcE - (calcPx*calcPx + calcPy*calcPy + calcPz*calcPz));

                            point[0] = calcPt;
                            lspoint[0] = calcPt;
                            point[1] = calcInvMass;
                            lspoint[1] = calcInvMass;
                            lspoint[2] = firstDecayTrack->Charge();
                            //Unlike sign pairs - create actual phi inv-mass distribution
                            if(firstDecayTrack->Charge() == 1 && secondDecayTrack->Charge() == -1){
                                fPhiInvMass->Fill(point);
                                
                                //Set-up TLorenztVector (px, py, pz, E), then push to vector
                                phi.SetPx(firstDecayTrack->Px()+secondDecayTrack->Px());
                                phi.SetPy(firstDecayTrack->Py()+secondDecayTrack->Py());
                                phi.SetPz(firstDecayTrack->Pz()+secondDecayTrack->Pz());
                                phi.SetE(firstDecayTrack->E()+secondDecayTrack->E());
                                phiCandidates.push_back(phi);
                            }else if(firstDecayTrack->Charge()*secondDecayTrack->Charge() == 1){
                                fPhiLikeSignInvMass->Fill(lspoint);
                            }
                        }
                        if(TMath::Abs(fpiTPCnSigma) < 2.0){
                            calcPx = firstDecayTrack->Px()+secondDecayTrack->Px();
                            calcPy = firstDecayTrack->Py()+secondDecayTrack->Py();
                            calcPz = firstDecayTrack->Pz()+secondDecayTrack->Pz();
                            calcPt = TMath::Sqrt(calcPx*calcPx + calcPy*calcPy);

                            calcE = firstDecayTrack->E() + secondDecayTrack->E();
                            calcInvMass = TMath::Sqrt(calcE*calcE - (calcPx*calcPx + calcPy*calcPy + calcPz*calcPz));

                            point[0] = calcPt;
                            lspoint[0] = calcPt;
                            point[1] = calcInvMass;
                            lspoint[1] = calcInvMass;
                            lspoint[2] = firstDecayTrack->Charge();
                            //Unlike sign pairs - create actual K*(892) inv-mass distribution
                            if(firstDecayTrack->Charge()*secondDecayTrack->Charge() == -1){
                                fKInvMass->Fill(point);
                                
                                //Set-up TLorenztVector (px, py, pz, E), then push to vector
                                K.SetPx(firstDecayTrack->Px()+secondDecayTrack->Px());
                                K.SetPy(firstDecayTrack->Py()+secondDecayTrack->Py());
                                K.SetPz(firstDecayTrack->Pz()+secondDecayTrack->Pz());
                                K.SetE(firstDecayTrack->E()+secondDecayTrack->E());
                                KCandidates.push_back(K);
                            }else if(firstDecayTrack->Charge()*secondDecayTrack->Charge() == 1){
                                fKLikeSignInvMass->Fill(lspoint);
                            }
                            
                        }
                    }
                }
            }

            //Cut on just real kaon+ from real phi, or real kaons from real K*(892), or real pi+ from K0 (MC Only)
            /*
            if((fPDG == 321 && TMath::Abs(motherPDG) == 333) || ((TMath::Abs(fPDG) == 321 && TMath::Abs(motherPDG) == 313)) || (fPDG == 211 && TMath::Abs(motherPDG)==310)){
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
                        Double_t fTPCnSigma = -999;
                        //check for labels
                        Int_t label = 0;
                        label = secondDecayTrack->GetLabel();

                        // Using stack to get actual particle PDG codes
                        Int_t fsecondPDG = 0;
                        Int_t secondMotherIndex = 0;
                        if(label > 1){             
                            if(fESD){
                                MCSecondDecay = fStack->Particle(label);
                                if(MCSecondDecay){
                                    fsecondPDG = MCSecondDecay->GetPdgCode();
                                    if(fsecondPDG != 0){
                                        secondMotherIndex = MCSecondDecay->GetMother(0);
                                    }
                                }
                            }else if(fAOD){
                                MCSecondDecayTrack = (AliAODMCParticle*)mcArray->At(label);
                                if(MCSecondDecayTrack){
                                    fsecondPDG = MCSecondDecayTrack->GetPdgCode();
                                }
                            }
                        }
                        if(fPDG == 321 && fsecondPDG == -321 && secondMotherIndex == motherIndex){
                            //Set-up TLorenztVector (px, py, pz, E), then push to vector
                            phi.SetPx(firstDecayTrack->Px()+secondDecayTrack->Px());
                            phi.SetPy(firstDecayTrack->Py()+secondDecayTrack->Py());
                            phi.SetPz(firstDecayTrack->Pz()+secondDecayTrack->Pz());
                            phi.SetE(firstDecayTrack->E()+secondDecayTrack->E());
                            if(firstDecayTrack->Pt() < 4 && secondDecayTrack->Pt() < 4){
                                fPhiDaughterPTKept->Fill(firstDecayTrack->Pt(), secondDecayTrack->Pt());
                                phiReals.push_back(phi);
                            }else{
                                fPhiDaughterPTCut->Fill(firstDecayTrack->Pt(), secondDecayTrack->Pt());
                            }
                            point[0] = TMath::Sqrt(phi.Px()*phi.Px() + phi.Py()*phi.Py());
                            point[1] = TMath::Sqrt(phi.E()*phi.E() - (phi.Px()*phi.Px() + phi.Py()*phi.Py() + phi.Pz()*phi.Pz()));
                            fTruthTracksPhiInvMass->Fill(point); 
                        }
                        if(TMath::Abs(fPDG) == 321 && TMath::Abs(fsecondPDG) == 211 && secondMotherIndex == motherIndex){
                            K.SetPx(firstDecayTrack->Px() + secondDecayTrack->Px());
                            K.SetPy(firstDecayTrack->Py() + secondDecayTrack->Py());
                            K.SetPz(firstDecayTrack->Pz() + secondDecayTrack->Px());
                            K.SetE(firstDecayTrack->E() + secondDecayTrack->E());
                            if(firstDecayTrack->Pt() < 4 && secondDecayTrack->Pt()<4){
                                fKstarDaughterPTKept->Fill(firstDecayTrack->Pt(), secondDecayTrack->Pt());
                                KReals.push_back(K);
                            }else{
                                fKstarDaughterPTCut->Fill(firstDecayTrack->Pt(), secondDecayTrack->Pt());
                            }
                        }
                        if(fPDG == 211 && fsecondPDG == -211 && secondMotherIndex == motherIndex){
                            K.SetPx(firstDecayTrack->Px() + secondDecayTrack->Px());
                            K.SetPy(firstDecayTrack->Py() + secondDecayTrack->Py());
                            K.SetPz(firstDecayTrack->Pz() + secondDecayTrack->Px());
                            K.SetE(firstDecayTrack->E() + secondDecayTrack->E());
                            if(firstDecayTrack->Pt() < 4 && secondDecayTrack->Pt()<4){
                                fK0DaughterPTKept->Fill(firstDecayTrack->Pt(), secondDecayTrack->Pt());
                                K0Reals.push_back(K);
                            }else{
                                fK0DaughterPTCut->Fill(firstDecayTrack->Pt(), secondDecayTrack->Pt());
                            }
                        }
                    }
                }
            }
            */
        }
    } 

    ////////////////
    // Track loop // (MC only right now...) 
    ////////////////
    /*
    TParticle* MCTriggerParticle = 0x0;
    AliAODMCParticle* MCTriggerTrack = 0x0;
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
        
        ////////////////////
        //Apply track cuts//
        ////////////////////
        if(fAOD)
            if(!atriggerTrack->TestFilterMask(AliAODTrack::kTrkGlobalNoDCA)) continue; //mimimum cuts
        
        if(fESD)
            if(!esdTrackCutsH->AcceptTrack(etriggerTrack))continue;
        
        
        ////////////////////
        //Track properties//
        ////////////////////
        Double_t dEdx =-999, fTPCnSigma=-999;
        dEdx = triggerTrack->GetTPCsignal();
       
        //Cut on p_T and eta
        if(triggerTrack->Pt() > 0.15 && TMath::Abs(triggerTrack->Eta()) < 0.8){
            fTrkPt->Fill(triggerTrack->Pt());
            fTrketa->Fill(triggerTrack->Eta());
            fTrkphi->Fill(triggerTrack->Phi());
            fdEdx->Fill(triggerTrack->P(),dEdx);
            fTPCNpts->Fill(triggerTrack->P(),triggerTrack->GetTPCsignalN());
            
            Double_t trigger_phi = triggerTrack->Phi();
            dphi_point[0] = triggerTrack->Pt();
            //check for labels
            Int_t label = 0;
            label = triggerTrack->GetLabel();

            // Using stack to get actual particle PDG codes (MC Only)
            /*
            Int_t fPDG = 0;
            if(label > 1){             
                if(fESD){
                    MCTriggerParticle = fStack->Particle(label);
                    if(MCTriggerParticle){
                        fPDG = MCTriggerParticle->GetPdgCode();
                        fPDGCodes->Fill(fPDG);
                        if(TMath::Abs(fPDG)==2212)fTPCKaonNSigp->Fill(triggerTrack->P(), fTPCnSigma);
                        if(TMath::Abs(fPDG)==321)fTPCKaonNSigK->Fill(triggerTrack->P(), fTPCnSigma);
                        if(TMath::Abs(fPDG)==211)fTPCKaonNSigPi->Fill(triggerTrack->P(), fTPCnSigma);
                        if(TMath::Abs(fPDG)==11)fTPCKaonNSige->Fill(triggerTrack->P(), fTPCnSigma);
                    }
                }else if(fAOD){
                    MCTriggerTrack = (AliAODMCParticle*)mcArray->At(label);
                    if(MCTriggerTrack){
                        fPDG = MCTriggerTrack->GetPdgCode();
                        if(TMath::Abs(fPDG)==2212)fTPCKaonNSigp->Fill(triggerTrack->P(), fTPCnSigma);
                        if(TMath::Abs(fPDG)==321)fTPCKaonNSigK->Fill(triggerTrack->P(), fTPCnSigma);
                        if(TMath::Abs(fPDG)==211)fTPCKaonNSigPi->Fill(triggerTrack->P(), fTPCnSigma);
                        if(TMath::Abs(fPDG)==11)fTPCKaonNSige->Fill(triggerTrack->P(), fTPCnSigma);
                    }
                }
            }
            */
    /*
            TParticle* MCFirstParticle=0x0;
            AliAODMCParticle* MCFirsttrk = 0x0;
            AliVParticle* vFirstTrack = 0x0;
            AliVTrack* firstTrack = 0x0;
            AliESDtrack* eFirstTrack = 0x0;
            AliAODTrack* aFirstTrack = 0x0;

            //Do Correlation Track Loop, finding correlation particles
            for(Int_t j_track = 0; j_track < ntracks; j_track++){
                //Exclude double counted particles
                if(j_track == i_track) continue;

                vFirstTrack = fVevent->GetTrack(j_track);

                if(!vFirstTrack) continue;

                firstTrack = dynamic_cast<AliVTrack*>(vFirstTrack);
                eFirstTrack = dynamic_cast<AliESDtrack*>(vFirstTrack);
                aFirstTrack = dynamic_cast<AliAODTrack*>(vFirstTrack);

                if(fESD){
                    if(!esdTrackCutsH->AcceptTrack(eFirstTrack))continue;
                }
                if(fAOD){
                    if(!aFirstTrack->TestFilterMask(AliAODTrack::kTrkGlobalNoDCA)) continue; //mimimum cuts
                }

                //cut on p_T and eta 
                if(firstTrack->Pt() > 0.15 && TMath::Abs(firstTrack->Eta()) < 0.8){

                    //Get Truth for first particle
                    Int_t firstLabel = 0, firstPDG = 0;
                    firstLabel = firstTrack->GetLabel();
                    if(fESD && firstLabel > 0){
                        MCFirstParticle = fStack->Particle(firstLabel);
                        if(MCFirstParticle){
                            firstPDG = MCFirstParticle->GetPdgCode();
                            if(TMath::Abs(firstPDG) == 2212){
                                dphi_point[1] = firstTrack->Pt();
                                dphi_point[2] = trigger_phi - firstTrack->Phi();
                                if(dphi_point[2] < -TMath::Pi()/2.0){
                                    dphi_point[2] += 2.0*TMath::Pi();
                                }else if(dphi_point[2] > 3.0*TMath::Pi()/2.0){
                                    dphi_point[2] -= 2.0*TMath::Pi();
                                }
                                fDphiHp->Fill(dphi_point);
                            }else if(TMath::Abs(firstPDG) == 211){
                                dphi_point[1] = firstTrack->Pt();
                                dphi_point[2] = trigger_phi - firstTrack->Phi();
                                if(dphi_point[2] < -TMath::Pi()/2.0){
                                    dphi_point[2] += 2.0*TMath::Pi();
                                }else if(dphi_point[2] > 3.0*TMath::Pi()/2.0){
                                    dphi_point[2] -= 2.0*TMath::Pi();
                                }
                                fDphiHPi->Fill(dphi_point);
                            }else if(TMath::Abs(firstPDG) == 321){
                                dphi_point[1] = firstTrack->Pt();
                                dphi_point[2] = trigger_phi - firstTrack->Phi();
                                if(dphi_point[2] < -TMath::Pi()/2.0){
                                    dphi_point[2] += 2.0*TMath::Pi();
                                }else if(dphi_point[2] > 3.0*TMath::Pi()/2.0){
                                    dphi_point[2] -= 2.0*TMath::Pi();
                                }
                                fDphiHK->Fill(dphi_point);
                            }
                        }
                    }
                    if(fAOD && firstLabel > 0){
                        MCFirsttrk = (AliAODMCParticle*)mcArray->At(firstLabel);
                        if(MCFirsttrk)
                            firstPDG = MCFirsttrk->GetPdgCode();
                    }
                }
            }
            
            //Loop over all phi reals
            Double_t cor_phi = 0, cor_pt = 0;

            for(Int_t i_particle = 0; i_particle < phiReals.size(); i_particle++){
                 cor_pt = TMath::Sqrt(phiReals[i_particle].Px()*phiReals[i_particle].Px() + phiReals[i_particle].Py()*phiReals[i_particle].Py());
                 cor_phi = TMath::Pi() + TMath::ATan2(-1*phiReals[i_particle].Py(), -1*phiReals[i_particle].Px());
                 dphi_point[1] = cor_pt;
                 dphi_point[2] = trigger_phi - cor_phi;
                 if(dphi_point[2] < -TMath::Pi()/2.0){
                     dphi_point[2] += 2.0*TMath::Pi();
                 }else if(dphi_point[2] > 3.0*TMath::Pi()/2.0){
                     dphi_point[2] -= 2.0*TMath::Pi();
                 }
                 fDphiHPhi->Fill(dphi_point);
            }
            //Loop over all K* reals
            for(Int_t j_particle = 0; j_particle < KReals.size(); j_particle++){
                 cor_pt = TMath::Sqrt(KReals[j_particle].Px()*KReals[j_particle].Px() + KReals[j_particle].Py()*KReals[j_particle].Py());
                 cor_phi = TMath::Pi() + TMath::ATan2(-1*KReals[j_particle].Py(), -1*KReals[j_particle].Px());
                 dphi_point[1] = cor_pt;
                 dphi_point[2] = trigger_phi - cor_phi;
                 if(dphi_point[2] < -TMath::Pi()/2.0){
                     dphi_point[2] += 2.0*TMath::Pi();
                 }else if(dphi_point[2] > 3.0*TMath::Pi()/2.0){
                     dphi_point[2] -= 2.0*TMath::Pi();
                 }
                 fDphiHKstar->Fill(dphi_point);
            }
            //Loop over all K0 reals
            for(Int_t k_particle = 0; k_particle < K0Reals.size(); k_particle++){
                 cor_pt = TMath::Sqrt(K0Reals[k_particle].Px()*K0Reals[k_particle].Px() + K0Reals[k_particle].Py()*K0Reals[k_particle].Py());
                 cor_phi = TMath::Pi() + TMath::ATan2(-1*K0Reals[k_particle].Py(), -1*K0Reals[k_particle].Px());
                 dphi_point[1] = cor_pt;
                 dphi_point[2] = trigger_phi - cor_phi;
                 if(dphi_point[2] < -TMath::Pi()/2.0){
                     dphi_point[2] += 2.0*TMath::Pi();
                 }else if(dphi_point[2] > 3.0*TMath::Pi()/2.0){
                     dphi_point[2] -= 2.0*TMath::Pi();
                 }
                 fDphiHK0->Fill(dphi_point);
            }
       }
    } //track loop
    */
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
