#ifndef AliAnalysisTaskhPhiCorr_cxx
#define AliAnalysisTaskhPhiCorr_cxx

//QA task for EMCAL electron analysis

class TH1F;
class AliEventPoolManager;
class THnSparse;
class AliESDEvent;
class AliAODEvent;

#include "AliAnalysisTaskSE.h"

class AliAnalysisTaskhPhiCorr : public AliAnalysisTaskSE {
public:
    AliAnalysisTaskhPhiCorr();
    AliAnalysisTaskhPhiCorr(const char *name);
    virtual ~AliAnalysisTaskhPhiCorr();
    
    virtual void   UserCreateOutputObjects();
    virtual void   UserExec(Option_t *option);
    virtual void   Terminate(Option_t *);

    void SetAODAnalysis() { SetBit(kAODanalysis, kTRUE); };
    void SetESDAnalysis() { SetBit(kAODanalysis, kFALSE); };
    Bool_t IsAODanalysis() const { return TestBit(kAODanalysis); };
    
private:
    enum{
        kAODanalysis = BIT(20),
    };
 
    struct AliKaonContainer{
        Int_t trackNum;
        TLorentzVector particle;
    };

    struct AliPhiContainer{
        Int_t daughter1TrackNum;
        Int_t daughter2TrackNum;
        TLorentzVector particle;
    };
   
    TObjArray* AddToTracks();
    void MakeCorrelations(Int_t itrack, AliVParticle *trigger, std::vector<AliPhiContainer> phiVec, THnSparse *fDphi, Float_t mult, Double_t zVtx);
    void MakeMixCorrelations(std::vector<AliPhiContainer> phiVec, THnSparse *fDphiMixed, Float_t mult, Double_t zVtx);
    void MakeHHMixCorrelations(AliCFParticle *cfPart, THnSparse *fDphiMixed, Float_t mult, Double_t zVtx);
  
    AliVEvent   *fVevent;  //!event object
    AliEventPoolManager *fPoolMgr; //! Event pool manager for mixed event
    AliESDEvent *fESD;    //!ESD object
    AliAODEvent *fAOD;    //!AOD object
    AliPIDResponse *fpidResponse; //!pid response
    AliMultSelection *fMultSelection; //!mult selection
    
    TList       *fOutputList; //!Output list
    TH1F        *fNevents;//! no of events
    TH1F        *fNumTracks;//! number of Tracks/evt
    TH1F        *fVtxZ;//!Vertex z
    TH1F        *fVtxX;//!Vertex x
    TH1F        *fVtxY;//!Vertex y
    TH2F        *fTrigMulti;//!trigger multiplicity
    TH1F        *fTrkPt;//!track pt
    TH1F        *fTrketa;//!track eta
    TH1F        *fTrkphi;//!track phi
    TH1F        *fHybridTrkPt;//!hybridTPC track pt
    TH1F        *fHybridTrketa;//!hybridTPC track eta
    TH1F        *fHybridTrkphi;//!hybridTPC track phi
    TH1F        *fHybridGlobalTrkPt;//!hybridGlobal track pt
    TH1F        *fHybridGlobalTrketa;//!hybridGlobal track eta
    TH1F        *fHybridGlobalTrkphi;//!hybridGlobal track phi
    TH2F        *fdEdx;//!dedx vs pt
    TH2F        *fTPCNpts;//!TPC Npoints used for dedx
    TH2F        *fTPCKaonNSig;//!TPC Nsigma

    THnSparseF  *fTrigDist;//! trigger distribution
    TH2D        *fMixStatZVtx;//! stats for mixed events
    TH1D        *fNoMixEvents;//! number of mixed events
    THnSparseF  *fKKUSDist;//! unlike sign kaon distribution
    THnSparseF  *fKKLSDist;//! like sign kaon distribution
    
    THnSparseF  *fDphiHPhi;//! delta-phi distribution with unlike sign kaon pairs
    THnSparseF  *fDphiHKK;//! delta-phi distribution with like sign kaon pairs
    THnSparseF  *fDphiHPhiMixed;//! hadron-US mixed correlation
    THnSparseF  *fDphiHKKMixed;//! hadron-LS mixed correlation
    THnSparseF  *fDphiHH;//! hadron-hadron correlation
    THnSparseF  *fDphiHHMixed;//! hadron-hadron mixed correlation

    AliAnalysisTaskhPhiCorr(const AliAnalysisTaskhPhiCorr&); // not implemented
    AliAnalysisTaskhPhiCorr& operator=(const AliAnalysisTaskhPhiCorr&); // not implemented
   
    ClassDef(AliAnalysisTaskhPhiCorr, 1); // example of analysis
};

#endif

