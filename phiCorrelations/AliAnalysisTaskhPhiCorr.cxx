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

#include "AliPID.h"
#include "AliESDpid.h"
#include "AliAODPid.h"
#include "AliPIDResponse.h"

#include "AliMCEvent.h"
#include "AliStack.h"
#include "AliMCEventHandler.h"

#include "AliAnalysisTaskhPhiCorr.h"

using std::cout;
using std::endl;

ClassImp(AliAnalysisTaskhPhiCorr)
//________________________________________________________________________
AliAnalysisTaskhPhiCorr::AliAnalysisTaskhPhiCorr(const char *name)
: AliAnalysisTaskSE(name),
fVevent(0),
fPoolMgr(0x0),
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
fMixStatZVtx(0),
fNoMixEvents(0),
fDphiHPhi(0),
fDphiHKK(0),
fDphiHPhiMixed(0),
fDphiHKKMixed(0)
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
AliAnalysisTaskhPhiCorr::AliAnalysisTaskhPhiCorr()
: AliAnalysisTaskSE("DefaultTask_HfeEMCQA"),
fVevent(0),
fPoolMgr(0),
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
fMixStatZVtx(0),
fNoMixEvents(0),
fDphiHPhi(0),
fDphiHKK(0),
fDphiHPhiMixed(0),
fDphiHKKMixed(0)
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
AliAnalysisTaskhPhiCorr::~AliAnalysisTaskhPhiCorr()
{
    //Destructor
    delete fOutputList;
}
//________________________________________________________________________
void AliAnalysisTaskhPhiCorr::UserCreateOutputObjects()
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

    Int_t poolSize = 1000;
    Int_t trackDepth = 100000;

    Int_t numVtxZBins = 10;
    Double_t vtxZBins[11] = {-10.0, -6.15, -3.90, -2.13, -0.59, 0.86, 2.29, 3.77, 5.39, 7.30, 10.0};

    Int_t numMultBins = 1;
    Double_t multBins[2] = {0.0, 10000.0};

    fPoolMgr = new AliEventPoolManager(poolSize, trackDepth, numMultBins, multBins, numVtxZBins, vtxZBins);


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
    
    // Histogram for trigger distribution
    Int_t trigBins[3] = {100,100,50};
    Double_t trigMin[3] = {0.1, 0.0, -2.0};
    Double_t trigMax[3] = {10.1, 6.28, 2.0};

    fTrigDist = new THnSparseF("fTrigDist", "Distribution for trigger particles", 3, trigBins, trigMin, trigMax);
    fOutputList->Add(fTrigDist);

    //Histograms for Mixed Event Stats
    fMixStatZVtx = new TH2D("fMixStatZVtx", "Mixed Event Statistics;NEvent in pool;Vtx_z", 200, 0, 5000, numVtxZBins, vtxZBins);
    fOutputList->Add(fMixStatZVtx);

    fNoMixEvents = new TH1D("fNoMixEvents", "Number of Mixed Events", 1, -0.5, 0.5);
    fOutputList->Add(fNoMixEvents);

    // Additional Histograms for US and LS Kaon pairs:
    Int_t bins[4] = {100, 200, 100, 50}; //pt, invmass, phi, eta
    Double_t min[4] = {0.1, 0.98, 0.0, -2.0};
    Double_t max[4] = {10.1, 1.1, 6.28, 2.0};
 
    fKKUSDist = new THnSparseF("fkkUSDist", "Distribution for all US Kaon pairs", 4, bins, min, max);
    fOutputList->Add(fKKUSDist);

    fKKLSDist = new THnSparseF("fkkLSDist", "Distribution for all LS Kaon pairs", 4, bins, min, max);
    fOutputList->Add(fKKLSDist);
  
    // Delta-phi histograms for different hadron-particle correlations (trigger pT, correlation pT, delta-phi, delta-eta, inv mass)
    Int_t dphi_bins[5]=    {85,   98,    128,   128, 120};
    Double_t dphi_min[5] = {3.0,   0.4, -1.57, -2.0, 0.98};
    Double_t dphi_max[5] = {20.0, 20.0,  4.71,  2.0, 1.1};

    fDphiHPhi = new THnSparseF("fDphiHPhi", "Hadron-#Phi #Delta#phi correlations", 5, dphi_bins, dphi_min, dphi_max);
    fOutputList->Add(fDphiHPhi);

    fDphiHKK = new THnSparseF("fDphiHKK", "Hadron-#KK likesign #Delta#phi correlations", 5, dphi_bins, dphi_min, dphi_max);
    fOutputList->Add(fDphiHKK);

    fDphiHPhiMixed = new THnSparseF("fDphiHPhiMixed", "Hadron-#Phi #Delta#phi mixed event Correlations", 5, dphi_bins, dphi_min, dphi_max);
    fOutputList->Add(fDphiHPhiMixed);

    fDphiHKKMixed = new THnSparseF("fDphiHKKMixed", "Hadron-#KK likesign #Delta#phi mixed event Correlations", 5, dphi_bins, dphi_min, dphi_max);
    fOutputList->Add(fDphiHKKMixed);

    PostData(1,fOutputList);
}

//___________________________________________________________________________
TObjArray* AliAnalysisTaskhPhiCorr::AddToTracks(){
    TObjArray* fArrayTracksMix = new TObjArray;
    fArrayTracksMix->SetOwner(kTRUE);

    Int_t ntracks = fVevent->GetNumberOfTracks();
    AliVParticle *particle= 0x0;
    AliAODTrack *aparticle= 0x0;
    for(int itrack=0; itrack<ntracks; itrack++){
        particle = fVevent->GetTrack(itrack);
        aparticle = dynamic_cast<AliAODTrack*>(particle);

        //cuts for the hadron track
        if(!aparticle->TestFilterMask(AliAODTrack::kTrkGlobalNoDCA)) continue;
        if(aparticle->Pt() > 0.15 && TMath::Abs(aparticle->Eta()) < 0.8){
            AliCFParticle *cfPart = new AliCFParticle(aparticle->Pt(), aparticle->Eta(), aparticle->Phi(), aparticle->Charge(), 0);
            fArrayTracksMix->Add(cfPart);
        } 
    }
    return fArrayTracksMix;
}

//___________________________________________________________________________
void AliAnalysisTaskhPhiCorr::MakeCorrelations(Int_t triggerIndex, AliVParticle *trigger, std::vector<AliPhiContainer> phiVec, THnSparse *fDphi){

    Double_t dphi_point[5];

    dphi_point[0] = trigger->Pt();
    for(int iphi = 0; iphi < phiVec.size(); iphi++){
        AliPhiContainer phi = phiVec[iphi];
        if(triggerIndex == phi.daughter1TrackNum || triggerIndex == phi.daughter2TrackNum) continue; //skip if hadron is one of the daughter particles
        dphi_point[1] = phi.particle.Pt();
        dphi_point[2] = trigger->Phi() - phi.particle.Phi();
        if(dphi_point[2] < -TMath::Pi()/2.0){
            dphi_point[2] += 2.0*TMath::Pi();
        }else if(dphi_point[2] > 3.0*TMath::Pi()/2.0){
            dphi_point[2] -= 2.0*TMath::Pi();
        }
        dphi_point[3] = trigger->Eta() - phi.particle.Eta();
        dphi_point[4] = phi.particle.M();
        fDphi->Fill(dphi_point);
    }
}

//___________________________________________________________________________
void AliAnalysisTaskhPhiCorr::MakeMixCorrelations(std::vector<AliPhiContainer> phiVec, THnSparse *fDphiMixed, Int_t mult, Double_t zVtx){

    Double_t dphi_point[5];
    AliEventPool* fPool;
    fPool = fPoolMgr->GetEventPool(mult, zVtx); // Get the buffer associated with the current multiplicity and z-vtx
    if (!fPool)
    {
        AliFatal(Form("No pool found for centrality = %f, zVtx = %f", mult, zVtx));
        return;
    }
    
    Int_t nMix = fPool->GetCurrentNEvents();
    fMixStatZVtx->Fill(nMix, zVtx);
    
    if(nMix >=5){
        fNoMixEvents->Fill(0);
        for(int iMix=0; iMix < nMix; iMix++){            
            TObjArray *tracks = fPool->GetEvent(iMix);
            
            for(int ihadron = 0; ihadron < tracks->GetEntriesFast(); ihadron++){
                AliVParticle *hadron = (AliVParticle*) tracks->At(ihadron);
                if(!hadron){
                    printf("ERROR: Could not receive mix pool track %d\n",ihadron);
                    continue;
                }
                for(int iphi = 0; iphi < phiVec.size(); iphi++){
                    AliPhiContainer phi = phiVec[iphi];
                    dphi_point[0] = hadron->Pt();
                    dphi_point[1] = phi.particle.Pt();
                    dphi_point[2] = hadron->Phi() - phi.particle.Phi();
                    if(dphi_point[2] < -TMath::Pi()/2.0){
                        dphi_point[2] += 2.0*TMath::Pi();
                    }else if(dphi_point[2] > 3.0*TMath::Pi()/2.0){
                        dphi_point[2] -= 2.0*TMath::Pi();
                    }
                    dphi_point[3] = hadron->Eta() - phi.particle.Eta();
                    dphi_point[4] = phi.particle.M();
                    fDphiMixed->Fill(dphi_point);
                }
            }
        }
    }
}


//________________________________________________________________________
void AliAnalysisTaskhPhiCorr::UserExec(Option_t *)
{

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
    if(fabs(Zvertex)>10.0)return;
    fNevents->Fill(2); //events after z vtx cut

    //Initialize the vectors/points that will be used to fill the histograms
    std::vector<AliPhiContainer> phiCandidates;
    std::vector<AliPhiContainer> phiLikeSignCandidates;
    std::vector<AliKaonContainer> kPlusList;
    std::vector<AliKaonContainer> kMinusList;

    Double_t distPoint[4] = {0, 0, 0, 0}; //pt, invmass, phi, eta
    Double_t trigPoint[3] = {0, 0, 0}; //pt, phi, eta
    Double_t dphi_point[5] = {0, 0, 0, 0, 0}; //trigger pt, phi pt, delta-phi, delta-eta, phi invmass

    AliVTrack *kaonTrack = 0x0;
    AliESDtrack *eKaonTrack = 0x0;
    AliAODTrack *aKaonTrack = 0x0;
    AliVParticle *vKaonTrack = 0x0;

    /* First Loop - Filling two vector for all Kaons (plus and minus) */
    for(Int_t itrack = 0; itrack < ntracks; itrack++){
        vKaonTrack = 0x0;
        vKaonTrack = fVevent->GetTrack(itrack);

        if(!vKaonTrack){
            printf("Error: Could not receive track %d\n", itrack);
            continue;
        }
        kaonTrack = dynamic_cast<AliVTrack*>(vKaonTrack);
        eKaonTrack = dynamic_cast<AliESDtrack*>(vKaonTrack);
        aKaonTrack = dynamic_cast<AliAODTrack*>(vKaonTrack);

        if(fAOD)
            if(!aKaonTrack->TestFilterMask(AliAODTrack::kTrkGlobalNoDCA)) continue; //mimimum cuts

        if(fESD)
            if(!esdTrackCutsH->AcceptTrack(eKaonTrack))continue;

        // Cut on pT and eta for possible Kaons
        if(kaonTrack->Pt() > 0.15 && TMath::Abs(kaonTrack->Eta()) < 0.8){
            Double_t fTPCnSigma = -999;
            Double_t fTOFnSigma = -999;
            Double_t fpiTPCnSigma = -999;
            //check for labels
            Int_t label = 0;
            label = kaonTrack->GetLabel();

            fTPCnSigma = fpidResponse->NumberOfSigmasTPC(kaonTrack, AliPID::kKaon);
            fTOFnSigma = fpidResponse->GetNumberOfSigmasTOF(kaonTrack, AliPID::kKaon);
            //Cut on kaon candidates
            if((TMath::Abs(fTPCnSigma) < 2.0) && (TMath::Abs(fTOFnSigma) < 2.0)){
                AliKaonContainer kaon;
                kaon.trackNum = itrack;
                kaon.particle.SetPx(kaonTrack->Px());
                kaon.particle.SetPy(kaonTrack->Py());
                kaon.particle.SetPz(kaonTrack->Pz());
                Double_t calcE = TMath::Sqrt(0.4937*0.4937 + kaonTrack->P()*kaonTrack->P());
                kaon.particle.SetE(calcE);

                if(kaonTrack->Charge() == 1){
                    kPlusList.push_back(kaon);
                }else{
                    kMinusList.push_back(kaon);
                }
            }
        }
    }

    // Go through the Kaon lists and create the phi candidates and like sign pairs
    // Also fill in the US and LS K pair distribution histograms
    AliPhiContainer phi;
    for(Int_t i_kplus = 0; i_kplus < kPlusList.size(); i_kplus++){
        for(Int_t j_kplus = i_kplus+1; j_kplus < kPlusList.size(); j_kplus++){
            phi.particle.SetPx(kPlusList[i_kplus].particle.Px() + kPlusList[j_kplus].particle.Px());
            phi.particle.SetPy(kPlusList[i_kplus].particle.Py() + kPlusList[j_kplus].particle.Py());
            phi.particle.SetPz(kPlusList[i_kplus].particle.Pz() + kPlusList[j_kplus].particle.Pz());
            phi.particle.SetE(kPlusList[i_kplus].particle.E() + kPlusList[j_kplus].particle.E());
            phi.daughter1TrackNum = kPlusList[i_kplus].trackNum;
            phi.daughter2TrackNum = kPlusList[j_kplus].trackNum;

            distPoint[0] = phi.particle.Pt();
            distPoint[1] = phi.particle.M();
            distPoint[2] = phi.particle.Phi() + TMath::Pi(); //adding pi to get number in range (0, 2pi)
            distPoint[3] = phi.particle.Eta();
            fKKLSDist->Fill(distPoint);
            phiLikeSignCandidates.push_back(phi);
       }
        for(Int_t i_kminus =0; i_kminus < kMinusList.size(); i_kminus++){
            phi.particle.SetPx(kPlusList[i_kplus].particle.Px() + kMinusList[i_kminus].particle.Px());
            phi.particle.SetPy(kPlusList[i_kplus].particle.Py() + kMinusList[i_kminus].particle.Py());
            phi.particle.SetPz(kPlusList[i_kplus].particle.Pz() + kMinusList[i_kminus].particle.Pz());
            phi.particle.SetE(kPlusList[i_kplus].particle.E() + kMinusList[i_kminus].particle.E());
            phi.daughter1TrackNum = kPlusList[i_kplus].trackNum;
            phi.daughter2TrackNum = kMinusList[i_kminus].trackNum;

            distPoint[0] = phi.particle.Pt();
            distPoint[1] = phi.particle.M();
            distPoint[2] = phi.particle.Phi() + TMath::Pi(); //adding pi to get number in range (0, 2pi)
            distPoint[3] = phi.particle.Eta();
            fKKUSDist->Fill(distPoint);
            phiCandidates.push_back(phi);
       }
    }
    for(Int_t i_kminus =0; i_kminus < kMinusList.size(); i_kminus++){
        for(Int_t j_kminus = i_kminus+1; j_kminus < kMinusList.size(); j_kminus++){
            phi.particle.SetPx(kMinusList[i_kminus].particle.Px() + kMinusList[j_kminus].particle.Px());
            phi.particle.SetPy(kMinusList[i_kminus].particle.Py() + kMinusList[j_kminus].particle.Py());
            phi.particle.SetPz(kMinusList[i_kminus].particle.Pz() + kMinusList[j_kminus].particle.Pz());
            phi.particle.SetE(kMinusList[i_kminus].particle.E() + kMinusList[j_kminus].particle.E());
            phi.daughter1TrackNum = kMinusList[i_kminus].trackNum;
            phi.daughter2TrackNum = kMinusList[j_kminus].trackNum;
 
            distPoint[0] = phi.particle.Pt();
            distPoint[1] = phi.particle.M();
            distPoint[2] = phi.particle.Phi() + TMath::Pi(); //adding pi to get number in range (0, 2pi)
            distPoint[3] = phi.particle.Eta();
            fKKLSDist->Fill(distPoint);
            phiLikeSignCandidates.push_back(phi); 
      }        
    }        

    ///////////////////////////////
    // Building d-phi histograms //
    ///////////////////////////////
    AliVTrack *triggerTrack = 0x0;
    AliESDtrack *etriggerTrack = 0x0;
    AliAODTrack *atriggerTrack = 0x0;
    AliVParticle* VtriggerTrack = 0x0;   
    for (Int_t itrack = 0; itrack < ntracks; itrack++) {

        VtriggerTrack = 0x0;
        VtriggerTrack  = fVevent->GetTrack(itrack);
        
        if (!VtriggerTrack) {
            printf("ERROR: Could not receive track %d\n", itrack);
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

            trigPoint[0] = triggerTrack->Pt();
            trigPoint[1] = triggerTrack->Phi();
            trigPoint[2] = triggerTrack->Eta();
            fTrigDist->Fill(trigPoint);
            
            MakeCorrelations(itrack, VtriggerTrack, phiCandidates, fDphiHPhi);
            MakeCorrelations(itrack, VtriggerTrack, phiLikeSignCandidates, fDphiHKK);
        }
    } //track loop


    ntracks = fVevent->GetNumberOfTracks();

    if(ntracks > 0 && TMath::Abs(Zvertex) < 10.0){
        MakeMixCorrelations(phiCandidates, fDphiHPhiMixed, ntracks, Zvertex);
        MakeMixCorrelations(phiLikeSignCandidates, fDphiHKKMixed, ntracks, Zvertex);
    

        AliEventPool *fPool;
        fPool = fPoolMgr->GetEventPool(ntracks, Zvertex); // Get the buffer associated with the current centrality and z-vtx
        if(!fPool){
            AliFatal(Form("No pool found for multiplicity = %f, zVtx = %f", ntracks, Zvertex));
            return;
        }
        TObjArray * fArrayTracksMix = AddToTracks();
        fArrayTracksMix->SetOwner(kTRUE);
        fPool->UpdatePool(fArrayTracksMix);
    }

    PostData(1, fOutputList);
}      
//________________________________________________________________________
void AliAnalysisTaskhPhiCorr::Terminate(Option_t *) 
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