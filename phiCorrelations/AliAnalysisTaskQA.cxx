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
fKKLSDist(0),
fTrigDist(0),
fDphiHPhi(0),
fDphiHKK(0)
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
fKKLSDist(0),
fTrigDist(0),
fDphiHPhi(0),
fDphiHKK(0)
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
    
    // TPC NSig histogram 
/*    fTPCKaonNSig = new TH2F("fTPCKaonNSig","All Track TPC Nsigma distribution;p (GeV/c);#sigma_{TPC-dE/dx}",1000,0,50,200,-10,10);
    fOutputList->Add(fTPCKaonNSig);
    
    fPDGCodes = new TH1F("fPDGCodes", "PDG codes of tracks", 2000, -1000, 1000);
    fOutputList->Add(fPDGCodes);
*/
    // Histogram for trigger distribution
    Int_t trigBins[3] = {100,100,50};
    Double_t trigMin[3] = {0.1, 0.0, -2.0};
    Double_t trigMax[3] = {10.1, 6.28, 2.0};

    fTrigDist = new THnSparseF("fTrigDist", "Distribution for trigger particles", 3, trigBins, trigMin, trigMax);
    fOutputList->Add(fTrigDist);

    // Additional Histograms for US and LS Kaon pairs:
    Int_t bins[4] = {100, 200, 100, 50}; //pt, invmass, phi, eta
    Double_t min[4] = {0.1, 0.98, 0.0, -2.0};
    Double_t max[4] = {10.1, 1.1, 6.28, 2.0};
 
    fKKUSDist = new THnSparseF("fkkUSDist", "Distribution for all US Kaon pairs", 4, bins, min, max);
    fOutputList->Add(fKKUSDist);

    fKKLSDist = new THnSparseF("fkkLSDist", "Distribution for all LS Kaon pairs", 4, bins, min, max);
    fOutputList->Add(fKKLSDist);
/*
    Int_t lsbins[3] = {100, 1000, 3};
    Double_t lsmin[3] = {0.0, 0.5, -1.1};
    Double_t lsmax[3] = {10.0, 2.0, 1.1};

    fPhiUSInvMass = new THnSparseF("fPhiUSInvMass", "Invariant mass dist. of unlike-sign Kaon pairs per p_{T}", 2, lsbins, lsmin, lsmax);
    fOutputList->Add(fPhiUSInvMass);

    fPhiLSInvMass = new THnSparseF("fPhiLikeSignInvMass", "Invariant mass distribution of like-sign Kaon pairs per p_{T}", 3, lsbins, lsmin, lsmax);
    fOutputList->Add(fPhiLSInvMass);

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
*/   
    // Delta-phi histograms for different hadron-particle correlations (trigger pT, correlation pT, delta-phi, delta-eta, inv mass)
    Int_t dphi_bins[5]=    {85,   98,    64,   24, 120};
    Double_t dphi_min[5] = {3.0,   0.4, -1.57, -3.0, 0.98};
    Double_t dphi_max[5] = {20.0, 20.0,  4.71,  3.0, 1.1};

    fDphiHPhi = new THnSparseF("fDphiHPhi", "Hadron-#Phi #Delta#phi correlations", 5, dphi_bins, dphi_min, dphi_max);
    fOutputList->Add(fDphiHPhi);

    fDphiHKK = new THnSparseF("fDphiHKK", "Hadron-#KK likesign #Delta#phi correlations", 5, dphi_bins, dphi_min, dphi_max);
    fOutputList->Add(fDphiHKK);
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
    std::vector<Int_t> phiDaughterTrackNum;
//    std::vector<TLorentzVector> phiReals;
    std::vector<TLorentzVector> phiLikesignCandidates;
    std::vector<Int_t> likesignDaughterTrackNum;
    std::vector<TLorentzVector> KCandidates;
//    std::vector<TLorentzVector> KReals;
    std::vector<TLorentzVector> K0Reals;
    TLorentzVector phi;
    TLorentzVector K;

    Double_t distPoint[4] = {0, 0, 0, 0}; //pt, invmass, phi, eta
    Double_t trigPoint[3] = {0, 0, 0}; //pt, phi, eta
    Double_t dphi_point[5] = {0, 0, 0, 0, 0}; //trigger pt, phi pt, delta-phi, delta-eta, phi invmass

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
            Double_t fTOFnSigma = -999;
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
            fTOFnSigma = fpidResponse->GetNumberOfSigmasTOF(firstDecayTrack, AliPID::kKaon);
//            TParticle *MCSecondDecay = 0x0;
//            AliAODMCParticle* MCSecondDecayTrack = 0x0;
            AliVTrack *secondDecayTrack = 0x0;
            AliESDtrack *eSecondDecayTrack = 0x0;
            AliAODTrack *aSecondDecayTrack = 0x0;
            AliVParticle *vSecondDecayTrack = 0x0;

            //Cut on kaon candidates
            if((TMath::Abs(fTPCnSigma) < 2.0) && (TMath::Abs(fTOFnSigma) < 2.0)){
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
                        fTOFnSigma = -999;
                        fTPCnSigma = fpidResponse->NumberOfSigmasTPC(secondDecayTrack, AliPID::kKaon);
                        fTOFnSigma = fpidResponse->GetNumberOfSigmasTOF(secondDecayTrack, AliPID::kKaon);
                        fpiTPCnSigma = fpidResponse->NumberOfSigmasTPC(secondDecayTrack, AliPID::kPion);
                        Double_t calcPx = 0.0, calcPy = 0.0, calcPz = 0.0;
                        Double_t calcE = 0.0, calcPt = 0.0, calcInvMass = 0.0, calcE1=0.0, calcE2=0.0;
                        if((TMath::Abs(fTPCnSigma) < 2.0) && (TMath::Abs(fTOFnSigma) < 2.0)){
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
                                    phiCandidates.push_back(phi);
                                    phiDaughterTrackNum.push_back(i_track);
                                    phiDaughterTrackNum.push_back(j_track);

                                    distPoint[0] = phi.Pt();
                                    distPoint[1] = calcInvMass;
                                    distPoint[2] = phi.Phi() + TMath::Pi(); //adding pi to get number in range (0, 2pi)
                                    distPoint[3] = phi.Eta();
                                    fKKUSDist->Fill(distPoint);
                                }
                            }else if(firstDecayTrack->Charge()*secondDecayTrack->Charge() == 1){
                                if(calcInvMass >= 0.98 && calcInvMass <=1.1){
                                    phi.SetPx(firstDecayTrack->Px()+secondDecayTrack->Px());
                                    phi.SetPy(firstDecayTrack->Py()+secondDecayTrack->Py());
                                    phi.SetPz(firstDecayTrack->Pz()+secondDecayTrack->Pz());
                                    phi.SetE(calcE);
                                    phiLikesignCandidates.push_back(phi);
                                    likesignDaughterTrackNum.push_back(i_track);
                                    likesignDaughterTrackNum.push_back(j_track);
                                    
                                    distPoint[0] = phi.Pt();
                                    distPoint[1] = calcInvMass;
                                    distPoint[2] = phi.Phi() + TMath::Pi(); //adding pi to get number in range (0, 2pi)
                                    distPoint[3] = phi.Eta();
                                    fKKLSDist->Fill(distPoint);
                                }

                            }
                        }
                        /*
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
                            
                        }*/
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
        
        
        ////////////////////
        //Track properties//
        ////////////////////
        Double_t dEdx =-999, fTPCnSigma=-999;
        dEdx = triggerTrack->GetTPCsignal();
       
        //Cut on p_T and eta
        if(triggerTrack->Pt() > 0.15 && TMath::Abs(triggerTrack->Eta()) < 0.8){
            //fTrkPt->Fill(triggerTrack->Pt());
            //fTrketa->Fill(triggerTrack->Eta());
            //fTrkphi->Fill(triggerTrack->Phi());
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
            */
            //Do Correlation Track Loop, finding correlation particles
            for(int i_phi = 0; i_phi < phiCandidates.size(); i_phi++){
                if(i_track == phiDaughterTrackNum[2*i_phi] || i_track == phiDaughterTrackNum[2*i_phi + 1]) continue; //skip if hadron is one of the daughter particles
                trigPoint[0] = triggerTrack->Pt();
                trigPoint[1] = triggerTrack->Phi();
                trigPoint[2] = triggerTrack->Eta();
                fTrigDist->Fill(trigPoint);
                dphi_point[1] = phiCandidates[i_phi].Pt();
                dphi_point[2] = trigger_phi - phiCandidates[i_phi].Phi();
                if(dphi_point[2] < -TMath::Pi()/2.0){
                    dphi_point[2] += 2.0*TMath::Pi();
                }else if(dphi_point[2] > 3.0*TMath::Pi()/2.0){
                    dphi_point[2] -= 2.0*TMath::Pi();
                }
                dphi_point[3] = triggerTrack->Eta() - phiCandidates[i_phi].Eta();
                dphi_point[4] = TMath::Sqrt(phiCandidates[i_phi].E()*phiCandidates[i_phi].E() - (phiCandidates[i_phi].Px()*phiCandidates[i_phi].Px() + phiCandidates[i_phi].Py()*phiCandidates[i_phi].Py() + phiCandidates[i_phi].Pz()*phiCandidates[i_phi].Pz()));
                fDphiHPhi->Fill(dphi_point);
            }

             for(int i_KK = 0; i_KK < phiLikesignCandidates.size(); i_KK++){
                if(i_track == likesignDaughterTrackNum[2*i_KK] || i_track == likesignDaughterTrackNum[2*i_KK + 1]) continue; //skip if hadron is one of the daughter particles
                dphi_point[1] = phiLikesignCandidates[i_KK].Pt();
                dphi_point[2] = trigger_phi - phiLikesignCandidates[i_KK].Phi();
                if(dphi_point[2] < -TMath::Pi()/2.0){
                    dphi_point[2] += 2.0*TMath::Pi();
                }else if(dphi_point[2] > 3.0*TMath::Pi()/2.0){
                    dphi_point[2] -= 2.0*TMath::Pi();
                }
                dphi_point[3] = triggerTrack->Eta() - phiLikesignCandidates[i_KK].Eta();
                dphi_point[4] = TMath::Sqrt(phiLikesignCandidates[i_KK].E()*phiLikesignCandidates[i_KK].E() - (phiLikesignCandidates[i_KK].Px()*phiLikesignCandidates[i_KK].Px() + phiLikesignCandidates[i_KK].Py()*phiLikesignCandidates[i_KK].Py() + phiLikesignCandidates[i_KK].Pz()*phiLikesignCandidates[i_KK].Pz()));
                fDphiHKK->Fill(dphi_point);
            }
       }
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
