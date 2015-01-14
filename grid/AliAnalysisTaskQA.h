#ifndef AliAnalysisTaskQA_cxx
#define AliAnalysisTaskQA_cxx

//QA task for EMCAL electron analysis

class TH1F;
class THnSparse;
class AliESDEvent;
class AliAODEvent;

#include "AliAnalysisTaskSE.h"

class AliAnalysisTaskQA : public AliAnalysisTaskSE {
public:
    AliAnalysisTaskQA();
    AliAnalysisTaskQA(const char *name);
    virtual ~AliAnalysisTaskQA();
    
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
    
    AliVEvent   *fVevent;  //!event object
    AliESDEvent *fESD;    //!ESD object
    AliAODEvent *fAOD;    //!AOD object
    AliPIDResponse *fpidResponse; //!pid response
    
    TList       *fOutputList; //!Output list
    TH1F        *fNevents;//! no of events
    TH1F        *fVtxZ;//!Vertex z
    TH1F        *fVtxX;//!Vertex x
    TH1F        *fVtxY;//!Vertex y
    TH2F        *fTrigMulti;//!trigger multiplicity
    TH1F        *fTrkPt;//!track pt
    TH1F        *fTrketa;//!track eta
    TH1F        *fTrkphi;//!track phi
    TH2F        *fdEdx;//!dedx vs pt
    TH2F        *fTPCNpts;//!TPC Npoints used for dedx
    TH2F        *fTPCnsig;//!TPC Nsigma
    
    AliAnalysisTaskQA(const AliAnalysisTaskQA&); // not implemented
    AliAnalysisTaskQA& operator=(const AliAnalysisTaskQA&); // not implemented
    
    ClassDef(AliAnalysisTaskQA, 1); // example of analysis
};

#endif

