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
#include "AliCFParticle.h"

#include "AliEventPoolManager.h"
#include "AliMultSelection.h"

#include "AliPID.h"
#include "AliESDpid.h"
#include "AliAODPid.h"
#include "AliPIDResponse.h"

#include "AliMCEvent.h"
#include "AliStack.h"
#include "AliMCEventHandler.h"

#include "AliAnalysisTaskPhiEff.h"

using std::cout;
using std::endl;

ClassImp(AliAnalysisTaskPhiEff)
//________________________________________________________________________
AliAnalysisTaskPhiEff::AliAnalysisTaskPhiEff(const char *name, Float_t multLow, Float_t multHigh)
: AliAnalysisTaskSE(name),
fVevent(0),
fPoolMgr(0x0),
fLSPoolMgr(0x0),
fHHPoolMgr(0x0),
fESD(0),
fAOD(0),
fpidResponse(0),
fOutputList(0),
fNevents(0),
fNumTracks(0),
fVtxZ(0),
fVtxX(0),
fVtxY(0),
fRealKDist(0),
fRealPiDist(0),
fRealpDist(0),
fRealeDist(0),
fRealMuonDist(0),
fRealPhiDist(0),
fRecoPhiDist(0)
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
    MULT_LOW = multLow;
    MULT_HIGH = multHigh;
    CENT_ESTIMATOR = "V0A";
    KAON_TRK_BIT = 1024;
    ASSOC_TRK_BIT = 1024;
    TRIG_TRK_BIT = 768;
    KAON_ETA_CUT = 0.8;

}
//________________________________________________________________________
AliAnalysisTaskPhiEff::AliAnalysisTaskPhiEff()
: AliAnalysisTaskSE("DefaultTask_HfeEMCQA"),
fVevent(0),
fPoolMgr(0x0),
fLSPoolMgr(0x0),
fHHPoolMgr(0x0),
fESD(0),
fAOD(0),
fpidResponse(0),
fOutputList(0),
fNevents(0),
fNumTracks(0),
fVtxZ(0),
fVtxX(0),
fVtxY(0),
fRealPhiDist(0),
fRecoPhiDist(0)
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
    MULT_LOW = 0.0;
    MULT_HIGH = 100.0;
    CENT_ESTIMATOR = "V0A";
    KAON_TRK_BIT = 1024;
    ASSOC_TRK_BIT = 1024;
    TRIG_TRK_BIT = 768;
    KAON_ETA_CUT = 0.8;
}
//________________________________________________________________________
AliAnalysisTaskPhiEff::~AliAnalysisTaskPhiEff()
{
    //Destructor
    delete fOutputList;
}
//________________________________________________________________________
void AliAnalysisTaskPhiEff::UserCreateOutputObjects()
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
   
    ////////////////////////////
    // Set-up for Mixed Event //
    ////////////////////////////


    Int_t numVtxZBins = 10;
    //Double_t vtxZBins[11] = {-10.0, -6.15, -3.90, -2.13, -0.59, 0.86, 2.29, 3.77, 5.39, 7.30, 10.0};
    Double_t vtxZBins[11] = {-10.0, -8.0, -6.0, -4.0, -2.0, 0.0, 2.0, 4.0, 6.0, 8.0, 10.0};
    Int_t numMultBins = 3;
    Double_t multBins[4] = {0.0, 20.0, 50.0, 90.0};


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

    fNumTracks = new TH1F("fNumTracks", "Number of Tracks/evt", 1000, 0, 1000);
    fOutputList->Add(fNumTracks);
    
    fVtxZ = new TH1F("fVtxZ","Z vertex position;Vtx_{z};counts",1000,-50,50);
    fOutputList->Add(fVtxZ);
    
    fVtxY = new TH1F("fVtxY","Y vertex position;Vtx_{y};counts",1000,-50,50);
    fOutputList->Add(fVtxY);
    
    fVtxX = new TH1F("fVtxX","X vertex position;Vtx_{x};counts",1000,-50,50);
    fOutputList->Add(fVtxX);
    
    // Single Particle and inclusive charged hadron histos
    Int_t numbinsSingle[5] = {95, 64, 64, 10, 10};
    Double_t minvalSingle[5] = {0.5, -3.14159, -2.0, -10.0, 0.0};
    Double_t maxvalSingle[5] = {10.0, 3.14159, 2.0, 10.0, 100.0};

    fRealChargedDist = new THnSparseF("fRealChargedDist", "Real Charged Hadron distribution;p_{T};#varphi;#eta;y;Z_{vtx};Multiplicity Percentile", 5, numbinsSingle, minvalSingle, maxvalSingle);
    fOutputList->Add(fRealChargedDist);

    fRealKDist = new THnSparseF("fRealKDist", "Real Kaon distribution;p_{T};#varphi;#eta;y;Z_{vtx};Multiplicity Percentile", 5, numbinsSingle, minvalSingle, maxvalSingle);
    fOutputList->Add(fRealKDist);

    fRealPiDist = new THnSparseF("fRealPiDist", "Real Pion distribution;p_{T};#varphi;#eta;y;Z_{vtx};Multiplicity Percentile", 5, numbinsSingle, minvalSingle, maxvalSingle);
    fOutputList->Add(fRealPiDist);

    fRealpDist = new THnSparseF("fRealpDist", "Real proton distribution;p_{T};#varphi;#eta;y;Z_{vtx};Multiplicity Percentile", 5, numbinsSingle, minvalSingle, maxvalSingle);
    fOutputList->Add(fRealpDist);

    fRealeDist = new THnSparseF("fRealeDist", "Real electron distribution;p_{T};#varphi;#eta;y;Z_{vtx};Multiplicity Percentile", 5, numbinsSingle, minvalSingle, maxvalSingle);
    fOutputList->Add(fRealeDist);

    fRealMuonDist = new THnSparseF("fRealMuonDist", "Real Muon distribution;p_{T};#varphi;#eta;y;Z_{vtx};Multiplicity Percentile", 5, numbinsSingle, minvalSingle, maxvalSingle);
    fOutputList->Add(fRealMuonDist);

    fRecoChargedDist = new THnSparseF("fRecoChargedDist", "Reco Charged Hadron distribution;p_{T};#varphi;#eta;y;Z_{vtx};Multiplicity Percentile", 5, numbinsSingle, minvalSingle, maxvalSingle);
    fOutputList->Add(fRecoChargedDist);

    fRecoKDist = new THnSparseF("fRecoKDist", "Reco Kaon distribution;p_{T};#varphi;#eta;y;Z_{vtx};Multiplicity Percentile", 5, numbinsSingle, minvalSingle, maxvalSingle);
    fOutputList->Add(fRecoKDist);
    
    fTOFKDist = new THnSparseF("fTOFKDist", "TOF Kaon distribution;p_{T};#varphi;#eta;y;Z_{vtx};Multiplicity Percentile", 5, numbinsSingle, minvalSingle, maxvalSingle);
    fOutputList->Add(fTOFKDist);
    
    fRecoPiDist = new THnSparseF("fRecoPiDist", "Reco Pion distribution;p_{T};#varphi;#eta;y;Z_{vtx};Multiplicity Percentile", 5, numbinsSingle, minvalSingle, maxvalSingle);
    fOutputList->Add(fRecoPiDist);

    fRecopDist = new THnSparseF("fRecopDist", "Reco proton distribution;p_{T};#varphi;#eta;y;Z_{vtx};Multiplicity Percentile", 5, numbinsSingle, minvalSingle, maxvalSingle);
    fOutputList->Add(fRecopDist);

    fRecoeDist = new THnSparseF("fRecoeDist", "Reco electron distribution;p_{T};#varphi;#eta;y;Z_{vtx};Multiplicity Percentile", 5, numbinsSingle, minvalSingle, maxvalSingle);
    fOutputList->Add(fRecoeDist);

    fRecoMuonDist = new THnSparseF("fRecoMuonDist", "Reco Muon distribution;p_{T};#varphi;#eta;y;Z_{vtx};Multiplicity Percentile", 5, numbinsSingle, minvalSingle, maxvalSingle);
    fOutputList->Add(fRecoMuonDist);

    fRecoChargedTriggerDist = new THnSparseF("fRecoChargedTriggerDist", "Reco Charged Hadron Trigger distribution;p_{T};#varphi;#eta;y;Z_{vtx};Multiplicity Percentile", 5, numbinsSingle, minvalSingle, maxvalSingle);
    fOutputList->Add(fRecoChargedTriggerDist);

    fRecoKTriggerDist = new THnSparseF("fRecoKTriggerDist", "Reco Kaon Trigger distribution;p_{T};#varphi;#eta;y;Z_{vtx};Multiplicity Percentile", 5, numbinsSingle, minvalSingle, maxvalSingle);
    fOutputList->Add(fRecoKTriggerDist);
    
    fRecoPiTriggerDist = new THnSparseF("fRecoPiTriggerDist", "Reco Pion Trigger distribution;p_{T};#varphi;#eta;y;Z_{vtx};Multiplicity Percentile", 5, numbinsSingle, minvalSingle, maxvalSingle);
    fOutputList->Add(fRecoPiTriggerDist);

    fRecopTriggerDist = new THnSparseF("fRecopTriggerDist", "Reco proton Trigger distribution;p_{T};#varphi;#eta;y;Z_{vtx};Multiplicity Percentile", 5, numbinsSingle, minvalSingle, maxvalSingle);
    fOutputList->Add(fRecopTriggerDist);

    fRecoeTriggerDist = new THnSparseF("fRecoeTriggerDist", "Reco electron Trigger distribution;p_{T};#varphi;#eta;y;Z_{vtx};Multiplicity Percentile", 5, numbinsSingle, minvalSingle, maxvalSingle);
    fOutputList->Add(fRecoeTriggerDist);

    fRecoMuonTriggerDist = new THnSparseF("fRecoMuonTriggerDist", "Reco Muon Trigger distribution;p_{T};#varphi;#eta;y;Z_{vtx};Multiplicity Percentile", 5, numbinsSingle, minvalSingle, maxvalSingle);
    fOutputList->Add(fRecoMuonTriggerDist);

    //Real & MC Phi track histos
    
    Int_t numbins[6] = {95, 64, 64, 10, 80, 10};
    Double_t minval[6] = {0.5, -3.14159, -2., -10, 0.99, 0.0};
    Double_t maxval[6] = {10.0, 3.14159,  2.,  10, 1.07, 100.0};

    fRealPhiDist = new THnSparseF("fRealPhiDist", "Real #phi distribution;p_{T};#varphi;#eta;y;Z_{vtx};m_{KK};Multiplicity Pctl.", 6, numbins, minval, maxval);
    fOutputList->Add(fRealPhiDist);

    fRealNoDecayCutPhiDist = new THnSparseF("fRealNoDecayCutPhiDist", "Real #phi distribution;p_{T};#varphi;#eta;y;Z_{vtx};m_{KK};Multiplicity Pctl.", 6, numbins, minval, maxval);
    fOutputList->Add(fRealNoDecayCutPhiDist);


    fRecoPhiDist = new THnSparseF("fRecoPhiDist", "Reco #phi distribution;p_{T};#varphi;#eta;y;Z_{vtx};m_{KK};Multiplicity Pctl.", 6, numbins, minval, maxval);
    fOutputList->Add(fRecoPhiDist);

    fTrackRecoPhiDist = new THnSparseF("fTrackRecoPhiDist", "Track Cut Reco #phi distribution;p_{T};#varphi;#eta;y;Z_{vtx};m_{KK};Multiplicity Pctl.",6, numbins, minval, maxval);
    fOutputList->Add(fTrackRecoPhiDist);

    fTOFRecoPhiDist = new THnSparseF("fTOFRecoPhiDist","Track Cut & TOF hit Reco #phi distribution;p_{T};#varphi;#eta;y;Z_{vtx};m_{KK};Multiplicity Pctl.",6, numbins, minval, maxval);
    fOutputList->Add(fTOFRecoPhiDist);
    
    fTPCPIDTrackRecoPhiDist = new THnSparseF("fTPCPIDTrackRecoPhiDist", "Track Cut & TPC PID 3#sigma Reco #phi distribution;p_{T};#varphi;#eta;y;Z_{vtx};m_{KK};Multiplicity Pctl.", 6, numbins, minval, maxval);
    fOutputList->Add(fTPCPIDTrackRecoPhiDist);
    
    fTPCPIDRecoPhiDist = new THnSparseF("fTPCPIDRecoPhiDist", "Track Cut & TOF hit & TPC PID 3#sigma Reco #phi distribution;p_{T};#varphi;#eta;y;Z_{vtx};m_{KK};Multiplicity Pctl.", 6, numbins, minval, maxval);
    fOutputList->Add(fTPCPIDRecoPhiDist);

    fPIDRecoPhiDist = new THnSparseF("fPIDRecoPhiDist", "Track Cut & TOF hit & TOF&TPC 3#sigma Reco #phi distribution;p_{T};#varphi;#eta;y;Z_{vtx};m_{KK};Multiplicity Pctl.", 6, numbins, minval, maxval);
    fOutputList->Add(fPIDRecoPhiDist);

    PostData(1,fOutputList);

    fReactionPlane = new TH1D("fReactionPlane", "Reaction Plane Angle; Angle", 64, -3.14159, 3.14159);
    fOutputList->Add(fReactionPlane);
}

//_______________________________________________________________________
UInt_t AliAnalysisTaskPhiEff::PassKaonCuts(AliAODTrack *track, Double_t TPCnSigma, Double_t TOFnSigma){
    //returns the level of cuts that the track passed
    //cutLevel: 1 = track cuts, 2 = TOF Hit, 4 = TPC PID cut, 8 = TOF PID cut
    UInt_t passLevel = 0;
    Bool_t pass = kTRUE;
    pass = pass && (TMath::Abs(track->Eta()) <= 0.8);
    pass = pass && (track->Pt() >= 0.15);
    pass = pass && (track->TestFilterMask(KAON_TRK_BIT));
    pass = pass && (track->GetTPCCrossedRows() > 80);
    if(pass) passLevel |= TRACK_BIT;    
    if(TOFnSigma != -999){ //check if there is a TOF signal, but don't care what the signal is
        passLevel |= TOF_HIT_BIT;
    }
    if(TMath::Abs(TPCnSigma) <= 3.0){ // check that kaon passed the TPC nsigma cut
        passLevel |= TPC_PID_BIT;
    }
    if(TMath::Abs(TOFnSigma) <= 3.0){ // check that kaon passed TOF nsigma cut
        passLevel |= TOF_PID_BIT;
    }

    return passLevel;
}

//_______________________________________________________________________
Bool_t AliAnalysisTaskPhiEff::PassHadronCuts(AliAODTrack *track, Bool_t isTrigger){
    //returns the level of cuts that the track passed
    //cutLevel: 1 = track cuts, 2 = TOF Hit, 4 = TPC PID cut, 8 = TOF PID cut
    Bool_t pass = kTRUE;
    pass = pass && (TMath::Abs(track->Eta()) <= 0.8);
    pass = pass && (track->Pt() >= 0.15);
    if(isTrigger){
        pass = pass && (track->TestBit(TRIG_TRK_BIT));
    }else{
        pass = pass && (track->TestFilterMask(ASSOC_TRK_BIT));
    }
    return pass;
}

//________________________________________________________________________
void AliAnalysisTaskPhiEff::UserExec(Option_t *){

    //masks for the different cut configurations: (track cuts only), (track + TOF hit), (track + TPC PID), (track + TOF hit + TPC PID), (track + TOF hit + TPC PID + TOF PID)
    UInt_t maskTrackOnly = TRACK_BIT;
    UInt_t maskTrackTOF = TRACK_BIT + TOF_HIT_BIT;
    UInt_t maskTrackTPC = TRACK_BIT + TPC_PID_BIT;
    UInt_t maskTrackTOFTPC = TRACK_BIT + TOF_HIT_BIT + TPC_PID_BIT;
    UInt_t maskTrackPID = TRACK_BIT + TOF_HIT_BIT + TPC_PID_BIT + TOF_PID_BIT;

    UInt_t evSelMask=((AliInputEventHandler*)(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler()))->IsEventSelected();
    
    fVevent = dynamic_cast<AliVEvent*>(InputEvent());
    if (!fVevent) {
        printf("ERROR: fVEvent not available\n");
        return;
    }
     
    

    fAOD = dynamic_cast<AliAODEvent*>(InputEvent());
    if (fAOD) {
        //printf("fAOD available\n");
        //return;
    }else{
        return;
    }

    ///////////////////
    //PID initialised//
    //////////////////
   fpidResponse = fInputHandler->GetPIDResponse();

    if(!fpidResponse){
        AliFatal("No PID response loaded!");
    }
    
    ////////////////
    //Event vertex//
    ///////////////
    Int_t ntracks = -999;
    ntracks = fVevent->GetNumberOfTracks();
    
    /////////////////
    //trigger check//
    /////////////////
    fVevent->GetFiredTriggerClasses();
    
    Int_t trigger = -1;
    //Multiplicity stuff
    Double_t multPercentile = 10.0;
    if (fAOD){
        //Double_t multiplicity=fAOD->GetHeader()->GetRefMultiplicity();
        AliAODHeader *header = dynamic_cast<AliAODHeader*>(fAOD->GetHeader());
        if(!header) AliFatal("Not a standard AOD");
        Double_t multiplicity = header->GetRefMultiplicity();
        /*
        fMultSelection = (AliMultSelection*)fAOD->FindListObject("MultSelection");
        if(fMultSelection){
            multPercentile = fMultSelection->GetMultiplicityPercentile(CENT_ESTIMATOR.Data());
        }else{
            return;
        }
        //if(multPercentile < MULT_LOW || multPercentile > MULT_HIGH) return;
        if(multPercentile < 0.0 || multPercentile > 100.0) return;
        */
    }
    
    fNevents->Fill(0); //all events
    Double_t Zvertex = -100, Xvertex = -100, Yvertex = -100;
    const AliVVertex *pVtx = fVevent->GetPrimaryVertex();
    Double_t NcontV = pVtx->GetNContributors();
    if(NcontV<3)return;
    fNevents->Fill(1); //events with 3 tracks

    Zvertex = pVtx->GetZ();
    Yvertex = pVtx->GetY();
    Xvertex = pVtx->GetX();
    fVtxZ->Fill(Zvertex);
    fVtxX->Fill(Xvertex);
    fVtxY->Fill(Yvertex);

    fNumTracks->Fill(ntracks);

    ////////////////////
    //event selection//
    ///////////////////
    if(fabs(Zvertex)>10.0)return;
    fNevents->Fill(2); //events after z vtx cut

    Double_t distPoint[6] = {0., 0., 0., 0., 0., 0.};
    Double_t singledistPoint[5] = {0., 0., 0., 0., 0.};
    //Loop over all particles in stack to get real phi, looking for phi->KK
    //AliMCEvent *fMCEvent = dynamic_cast<AliMCEvent*>(InputEvent());
    fMCArray = dynamic_cast<TClonesArray*>(fAOD->FindListObject(AliAODMCParticle::StdBranchName()));
    if(!fMCArray){
        AliError("Array of MC particles not found");
        return;
    }

    fMCHeader = dynamic_cast<AliAODMCHeader*>(fAOD->GetList()->FindObject(AliAODMCHeader::StdBranchName()));
    if (!fMCHeader) {
        AliError("Could not find MC Header in AOD");
        return;
    }

    Double_t RA = fMCHeader->GetReactionPlaneAngle();
    if(RA > TMath::Pi()){
        RA -= 2.0*TMath::Pi();
    }else if(RA < -1.0*TMath::Pi()){
        RA += 2.0*TMath::Pi();
    }
    fReactionPlane->Fill(RA);

    //TRandom3* randomGen = new TRandom3();
    //Double_t randangle = randomGen->Uniform(0.0, 2.0*TMath::Pi());

    UInt_t negPassCuts = 0;
    UInt_t posPassCuts = 0;

    Int_t negparPDG = 0;
    Int_t posparPDG = 0;
    Float_t recoPx = 0.0;
    Float_t recoPy = 0.0;
    Float_t recoPz = 0.0;
    Float_t recoP = 0.0;
    Float_t recoE = 0.0;
    Float_t recoM = 0.0;
    Float_t recoEta = 0.0;
    Float_t recoY = 0.0;
    Float_t recoPt = 0.0;
    Float_t recoPhi = 0.0;

    Int_t motherIndex = 0;
    Int_t motherPDG = 0;
    Int_t pdgcode = 0;
    
    //loop over tracks to get kaons, find phi daughters and fill Reco Dist
    for(int itrack = 0; itrack < ntracks; itrack++){
        AliVParticle *vnegpart = dynamic_cast<AliVParticle*>(fVevent->GetTrack(itrack));
        AliVTrack *negtrack = dynamic_cast<AliVTrack*>(vnegpart);
        AliAODTrack *aodnegtrack = dynamic_cast<AliAODTrack*>(vnegpart);

        Int_t tracklabel = aodnegtrack->GetLabel();
        if(tracklabel < 0) continue;
    
        Double_t negTPCnSigma = -999;
        Double_t negTOFnSigma = -999;
        if(negtrack->Pt() > 0.15){
            negTPCnSigma = fpidResponse->NumberOfSigmasTPC(negtrack, AliPID::kKaon);
            negTOFnSigma = fpidResponse->NumberOfSigmasTOF(negtrack, AliPID::kKaon);
          /*negTOFnSigma = negtrack->GetTOFsignal();
            if(negTOFnSigma == 0.0){
                negTOFnSigma = -999;
            }
          */  
        }

        //Get all single particle distributions
        AliAODMCParticle* mcpart = (AliAODMCParticle*)fMCArray->At(tracklabel);
        if(mcpart->IsPhysicalPrimary()){
            pdgcode = mcpart->GetPdgCode();
            Bool_t pass = PassHadronCuts(aodnegtrack, kFALSE);
            //Bool_t trigpass = PassHadronCuts(aodnegtrack, kTRUE);
            Bool_t trigpass = kFALSE;
            singledistPoint[0] = aodnegtrack->Pt();
            singledistPoint[1] = aodnegtrack->Phi();
            if(singledistPoint[1] > TMath::Pi()){
                singledistPoint[1] -= 2.0*TMath::Pi();
            }else if(singledistPoint[1] < -1.0*TMath::Pi()){
                singledistPoint[1] += 2.0*TMath::Pi();
            }
            singledistPoint[2] = aodnegtrack->Eta();
            singledistPoint[3] = Zvertex;
            singledistPoint[4] = multPercentile;
            if(pass){ 
                if(TMath::Abs(pdgcode) == 321){
                    Int_t kaonPass = PassKaonCuts(aodnegtrack, negTPCnSigma, negTOFnSigma);
                    fRecoKDist->Fill(singledistPoint);
                    fRecoChargedDist->Fill(singledistPoint);
                    if((kaonPass & maskTrackTOF) == maskTrackTOF){
                        fTOFKDist->Fill(singledistPoint);
                    }
                }else if(TMath::Abs(pdgcode) == 211){
                    fRecoPiDist->Fill(singledistPoint);
                    fRecoChargedDist->Fill(singledistPoint);
                }else if(TMath::Abs(pdgcode) == 2212){
                    fRecopDist->Fill(singledistPoint);
                    fRecoChargedDist->Fill(singledistPoint);
                }else if(TMath::Abs(pdgcode) == 11){
                    fRecoeDist->Fill(singledistPoint);
                    fRecoChargedDist->Fill(singledistPoint);
                }else if(TMath::Abs(pdgcode) == 13){
                    fRecoMuonDist->Fill(singledistPoint);
                    fRecoChargedDist->Fill(singledistPoint);
                }
            }
            if(trigpass){
                if(TMath::Abs(pdgcode) == 321){
                    Int_t kaonPass = PassKaonCuts(aodnegtrack, negTPCnSigma, negTOFnSigma);
                    fRecoKTriggerDist->Fill(singledistPoint);
                    fRecoChargedTriggerDist->Fill(singledistPoint);
                }else if(TMath::Abs(pdgcode) == 211){
                    fRecoPiTriggerDist->Fill(singledistPoint);
                    fRecoChargedTriggerDist->Fill(singledistPoint);
                }else if(TMath::Abs(pdgcode) == 2212){
                    fRecopTriggerDist->Fill(singledistPoint);
                    fRecoChargedTriggerDist->Fill(singledistPoint);
                }else if(TMath::Abs(pdgcode) == 11){
                    fRecoeTriggerDist->Fill(singledistPoint);
                    fRecoChargedTriggerDist->Fill(singledistPoint);
                }else if(TMath::Abs(pdgcode) == 13){
                    fRecoMuonTriggerDist->Fill(singledistPoint);
                    fRecoChargedTriggerDist->Fill(singledistPoint);
                }

            }

        }

        //GetKaons that came from phi for phi reco
        negPassCuts = PassKaonCuts(aodnegtrack,  negTPCnSigma, negTOFnSigma);
        if(negPassCuts == 0) continue;

        AliAODMCParticle* mcnegpart = (AliAODMCParticle*)fMCArray->At(tracklabel);
        negparPDG = mcnegpart->GetPdgCode();
        if(negparPDG != -321) continue;

        motherIndex = mcnegpart->GetMother();
        if(motherIndex < 0) continue;

        AliAODMCParticle* mcmother = (AliAODMCParticle*)fMCArray->At(motherIndex);
        motherPDG = mcmother->GetPdgCode();
        if(motherPDG != 333) continue;

        for(int jtrack = 0; jtrack < ntracks; jtrack++){
            AliVParticle *vpospart = dynamic_cast<AliVParticle*>(fVevent->GetTrack(jtrack));
            AliVTrack *postrack = dynamic_cast<AliVTrack*>(vpospart);
            AliAODTrack *aodpostrack = dynamic_cast<AliAODTrack*>(vpospart);

            Double_t posTPCnSigma = -999;
            Double_t posTOFnSigma = -999;
            if(postrack->Pt() > 0.15){
                posTPCnSigma = fpidResponse->NumberOfSigmasTPC(postrack, AliPID::kKaon);
                posTOFnSigma = fpidResponse->NumberOfSigmasTOF(postrack, AliPID::kKaon);
              /*  posTOFnSigma = postrack->GetTOFsignal();
                if(posTOFnSigma == 0.0){
                    posTOFnSigma = -999;
                }
              */  
            }

            posPassCuts = PassKaonCuts(aodpostrack, posTPCnSigma, posTOFnSigma);
            if(posPassCuts == 0) continue;

            Int_t postracklabel = aodpostrack->GetLabel();
            if(postracklabel < 0) continue;

            AliAODMCParticle* mcpospart = (AliAODMCParticle*)fMCArray->At(postracklabel);
            posparPDG = mcpospart->GetPdgCode();
            if(posparPDG != 321) continue;

            if(mcpospart->GetMother() == motherIndex){
                recoPx = aodnegtrack->Px() + aodpostrack->Px();
                recoPy = aodnegtrack->Py() + aodpostrack->Py();
                recoPz = aodnegtrack->Pz() + aodpostrack->Pz();
                
                recoP = TMath::Sqrt(recoPx*recoPx + recoPy*recoPy + recoPz*recoPz);
                recoE = TMath::Sqrt(aodnegtrack->Px()*aodnegtrack->Px() + aodnegtrack->Py()*aodnegtrack->Py() + aodnegtrack->Pz()*aodnegtrack->Pz() + 0.4937*0.4937) + TMath::Sqrt(aodpostrack->Px()*aodpostrack->Px() + aodpostrack->Py()*aodpostrack->Py() + aodpostrack->Pz()*aodpostrack->Pz() + 0.4937*0.4937);
                recoM = TMath::Sqrt(recoE*recoE - recoP*recoP);
                recoPt = TMath::Sqrt(recoPx*recoPx + recoPy*recoPy);
                recoEta = 0.5*TMath::Log((recoP + recoPz)/(recoP -  recoPz));
                recoY = 0.5*TMath::Log((recoE + recoPz)/(recoE - recoPz));
                recoPhi = TMath::ATan2(recoPy, recoPx);
                if(recoPhi < -TMath::Pi()){
                    recoPhi += 2.0*TMath::Pi();
                }else if(recoPhi > TMath::Pi()){
                    recoPhi -= 2.0*TMath::Pi();
                }
                distPoint[0] = recoPt;
                distPoint[1] = recoPhi;
                distPoint[2] = recoEta;
                distPoint[3] = Zvertex;
                distPoint[4] = recoM;
                distPoint[5] = multPercentile;
                
                fRecoPhiDist->Fill(distPoint);
                //fill with phi's where daughter kaons pass track cuts
                if(((negPassCuts & maskTrackOnly) == maskTrackOnly) && ((posPassCuts & maskTrackOnly)== maskTrackOnly)){
                    fTrackRecoPhiDist->Fill(distPoint);
                }
                //fill with phi's where daughter kaons pass track cuts and have a TOF hit
                if(((negPassCuts & maskTrackTOF) == maskTrackTOF) && ((posPassCuts & maskTrackTOF)== maskTrackTOF)){
                    fTOFRecoPhiDist->Fill(distPoint);
                }
               
                //fill with phi's where daughter kaons pass track cuts and TPC PID cuts
                if(((negPassCuts & maskTrackTPC) == maskTrackTPC) && ((posPassCuts & maskTrackTPC)== maskTrackTPC)){
                    fTPCPIDTrackRecoPhiDist->Fill(distPoint);
                } 
                
                //fill with phi's where daughter kaons pass track cuts and TOF hit and TPC PID cuts 
                if(((negPassCuts & maskTrackTOFTPC) == maskTrackTOFTPC) && ((posPassCuts & maskTrackTOFTPC)== maskTrackTOFTPC)){
                    fTPCPIDRecoPhiDist->Fill(distPoint);
                }

                //fill with phi daughter kaons that pass track cuts and pass TOF&TPC PID cuts
                if(((negPassCuts & maskTrackPID) == maskTrackPID) && ((posPassCuts & maskTrackPID)== maskTrackPID)){
                    fPIDRecoPhiDist->Fill(distPoint);
                }
            }
        }
    }

    //loop over MC particles to get original phi and fill Real dist (and real single particle dist
    
    for(Int_t imcpart=0; imcpart< fMCArray->GetEntries(); imcpart++){
        AliAODMCParticle *AODMCtrack = (AliAODMCParticle*)fMCArray->At(imcpart);
        pdgcode = TMath::Abs(AODMCtrack->GetPdgCode());
        if(AODMCtrack->IsPhysicalPrimary()){
            singledistPoint[0] = AODMCtrack->Pt();
            singledistPoint[1] = AODMCtrack->Phi();
            if(singledistPoint[1] > TMath::Pi()){
                singledistPoint[1] -= 2.0*TMath::Pi();
            }else if(singledistPoint[1] < -1.0*TMath::Pi()){
                singledistPoint[1] += 2.0*TMath::Pi();
            }
            singledistPoint[2] = AODMCtrack->Eta();
            singledistPoint[3] = Zvertex;
            singledistPoint[4] = multPercentile;
            if(TMath::Abs(pdgcode) == 321){
                fRealKDist->Fill(singledistPoint);
                fRealChargedDist->Fill(singledistPoint);
            }else if(TMath::Abs(pdgcode) == 211){
                fRealPiDist->Fill(singledistPoint);
                fRealChargedDist->Fill(singledistPoint);
            }else if(TMath::Abs(pdgcode) == 2212){
                fRealpDist->Fill(singledistPoint);
                fRealChargedDist->Fill(singledistPoint);
            }else if(TMath::Abs(pdgcode) == 11){
                fRealeDist->Fill(singledistPoint);
                fRealChargedDist->Fill(singledistPoint);
            }else if(TMath::Abs(pdgcode) == 13){
                fRealMuonDist->Fill(singledistPoint);
                fRealChargedDist->Fill(singledistPoint);
            }

        }

        //select phis
        if(pdgcode != 333) continue;
        Int_t indexFirstDaughter = 0, indexSecondDaughter = 0;
        indexFirstDaughter = AODMCtrack->GetDaughterFirst();
        indexSecondDaughter = AODMCtrack->GetDaughterLast();

        //indexFirstDaughter = AODMCtrack->GetDaughter(0);
        //indexSecondDaughter = AODMCtrack->GetDaughter(1);

        if(indexFirstDaughter < 0 || indexSecondDaughter < 0) continue;
        AliAODMCParticle* firstDaughter = (AliAODMCParticle*)fMCArray->At(indexFirstDaughter);
        AliAODMCParticle* secondDaughter = (AliAODMCParticle*)fMCArray->At(indexSecondDaughter);

        //select only phi that decay to two kaons
        if(TMath::Abs(firstDaughter->GetPdgCode()) == 321 && TMath::Abs(secondDaughter->GetPdgCode()) == 321 && (firstDaughter->GetPdgCode())*(secondDaughter->GetPdgCode()) <0 && TMath::Abs(firstDaughter->Eta()) <= 0.8 && TMath::Abs(secondDaughter->Eta()) <= 0.8){
            distPoint[0] = AODMCtrack->Pt();
            distPoint[1] = AODMCtrack->Phi();
            if(distPoint[1] > TMath::Pi()){
                distPoint[1] -= 2.0*TMath::Pi();
            }else if(distPoint[1] < -1.0*TMath::Pi()){
                distPoint[1] += 2.0*TMath::Pi();
            }
            distPoint[2] = AODMCtrack->Eta();
            distPoint[3] = Zvertex;
            distPoint[4] = AODMCtrack->GetCalcMass();
            distPoint[5] = multPercentile;
            fRealPhiDist->Fill(distPoint);
        } 

        if(TMath::Abs(firstDaughter->GetPdgCode()) == 321 && TMath::Abs(secondDaughter->GetPdgCode()) == 321 && (firstDaughter->GetPdgCode())*(secondDaughter->GetPdgCode())){
            distPoint[0] = AODMCtrack->Pt();
            distPoint[1] = AODMCtrack->Phi();
            if(distPoint[1] > TMath::Pi()){
                distPoint[1] -= 2.0*TMath::Pi();
            }else if(distPoint[1] < -1.0*TMath::Pi()){
                distPoint[1] += 2.0*TMath::Pi();
            }
            distPoint[2] = AODMCtrack->Eta();
            distPoint[3] = Zvertex;
            distPoint[4] = AODMCtrack->GetCalcMass();
            distPoint[5] = multPercentile;
            fRealNoDecayCutPhiDist->Fill(distPoint);
        } 
    }

    PostData(1, fOutputList);
}    
//________________________________________________________________________
void AliAnalysisTaskPhiEff::Terminate(Option_t *) 
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
