AliAnalysisTask *AddTaskQA(Float_t multLow = 0.0, Float_t multHigh = 100.0){
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
    
    stringstream taskStream;
    taskStream << "hPhiCorr_mult_" << multLow << "_" << multHigh;

    TString taskName = taskStream.str();

    AliAnalysisTaskhPhiCorr *hPhiCorr = new AliAnalysisTaskhPhiCorr(taskName.Data(), multLow, multHigh); 
    hPhiCorr->SelectCollisionCandidates(AliVEvent::kINT7);
    
    stringstream fileStream;
    fileStream << "phiCorrelations_mult_" << multLow << "_" << multHigh << ".root";
    TString filename = fileStream.str();

    stringstream containerStream;
    containerStream << "phiCorr_mult_" << multLow << "_" << multHigh;
    TString containerName = containerStream.str();

    printf("\n!!!!!!!!!!!!!\n Setting up input/output containers\n");
    fflush(stdout);
    AliAnalysisDataContainer *cinput  = mgr->GetCommonInputContainer();
    AliAnalysisDataContainer *coutput1 = mgr->CreateContainer(containerName, TList::Class(),AliAnalysisManager::kOutputContainer, filename.Data());
    printf("\n!!!!!!!!!!!!!\n Connecting input and output containers\n");
    fflush(stdout);
    mgr->ConnectInput(hPhiCorr, 0, cinput);
    mgr->ConnectOutput(hPhiCorr, 1, coutput1);
    mgr->AddTask(hPhiCorr);
    printf("\n!!!!!!!!!!!!\n Done with addtask macro \n");
    fflush(stdout);
    return hPhiCorr;
}
