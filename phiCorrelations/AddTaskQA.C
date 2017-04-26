AliAnalysisTask *AddTaskQA(){
    //get the current analysis manager
    AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
    if (!mgr) {
        Error("AddTaskHFE", "No analysis manager found.");
        return NULL;
    }
    if (!mgr->GetInputEventHandler()) {
        ::Error("AddTaskHFE", "This task requires an input event handler");
        return NULL;
    }
    TString type = mgr->GetInputEventHandler()->GetDataType(); // can be "ESD" or "AOD"
    
/*    Bool_t MCthere=kFALSE;
    AliMCEventHandler *mcH = dynamic_cast<AliMCEventHandler*>(mgr->GetMCtruthEventHandler());
    if(!mcH){
        MCthere=kFALSE;
    }else{
        MCthere=kTRUE;
    }
  */  
    //char calib[100];
    //    sprintf(calib,"QA");

    printf("\n!!!!!!!!!!!!\nSetting up AliAnalysisTaskQA\n");
    fflush(stdout); 
    AliAnalysisTaskhPhiCorr *hPhiCorr = new AliAnalysisTaskhPhiCorr("hPhiCOrr"); 
    hPhiCorr->SelectCollisionCandidates(AliVEvent::kINT7);
    
    TString containerName7 = mgr->GetCommonFileName();
    containerName7 += ":PhiReconstruction";
    printf("\n!!!!!!!!!!!!!\n Setting up input/output containers\n");
    fflush(stdout);
    AliAnalysisDataContainer *cinput  = mgr->GetCommonInputContainer();
    AliAnalysisDataContainer *coutput1 = mgr->CreateContainer("InvMass", TList::Class(),AliAnalysisManager::kOutputContainer, containerName7.Data());
    printf("\n!!!!!!!!!!!!!\n Connecting input and output containers\n");
    fflush(stdout);
    mgr->ConnectInput(hPhiCorr, 0, cinput);
    mgr->ConnectOutput(hPhiCorr, 1, coutput1);
    mgr->AddTask(hPhiCorr);
    printf("\n!!!!!!!!!!!!\n Done with addtask macro \n");
    fflush(stdout);
    return hPhiCorr;
}
