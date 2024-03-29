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
    AliAnalysisTaskQA *hfecalqa7 = new AliAnalysisTaskQA("emcqa");
    mgr->AddTask(hfecalqa7);
    //hfecalqa7->SelectCollisionCandidates(AliVEvent::kINT7);
    
    TString containerName7 = mgr->GetCommonFileName();
    containerName7 += ":PhiPIDCuts";
    printf("\n!!!!!!!!!!!!!\n Setting up input/output containers\n");
    fflush(stdout);
    AliAnalysisDataContainer *cinput  = mgr->GetCommonInputContainer();
    AliAnalysisDataContainer *coutput1 = mgr->CreateContainer("InvMass", TList::Class(),AliAnalysisManager::kOutputContainer, containerName7.Data());
    printf("\n!!!!!!!!!!!!!\n Connecting input and output containers\n");
    fflush(stdout);
    mgr->ConnectInput(hfecalqa7, 0, cinput);
    mgr->ConnectOutput(hfecalqa7, 1, coutput1);
    printf("\n!!!!!!!!!!!!\n Done with addtask macro \n");
    fflush(stdout);
    return NULL;
}
