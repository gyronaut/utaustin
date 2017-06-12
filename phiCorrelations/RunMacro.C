/*
 *  
 *
 *
 * */
#include "TRoot.h"
#include "TRint.h"
#include "TSystem.h"

void RunMacro()
{
//     check the function for asymmetric TPC cut in ConfigHFEemcalMod....the rest is still necessary????

   // Firstly, set some variables
   const char* launch = "grid"; // grid, local (if your data is on your local machine, doesn't connect at all)
   const char*  mode = "test"; //test, full, terminate  (test= connect to grid but run locally, full= run on grid, terminate= merge output on grid)
   Bool_t pre_final_stage = kTRUE; //TRUE = merging done on grid, FALSE = merge happens locally   
   Int_t cyclenumber = 1;
   Bool_t debug = kTRUE;
   char* work_dir = "PhiCorrelations_LHC16q";
   char* output_dir = "output_2017_06_12";
   Int_t ttl = 50000;
   Int_t noffiles = 20;
//   Int_t runcycle[]={0,32};
   Int_t runcycle[]={0,1,5,11,18,24,32};
   Bool_t UseParfiles = kFALSE;

// create and customize the alien handler
  AliAnalysisAlien *alienHandler = new AliAnalysisAlien();
      
    alienHandler->AddIncludePath("-I. -I$ROOTSYS/include -I$ALICE_ROOT -I$ALICE_ROOT/EMCAL -I$ALICE_ROOT/PYTHIA6 -I$ALICE_ROOT/ANALYSIS -I$ALICE_PHYSICS/PWGGA -I$ALICE_PHYSICS/PWGHF -I$ALICE_PHYSICS/PWGHF/hfe -I$ALICE_ROOT/CONTAINERS -I$ALICE_ROOT/STEER/STEER -I$ALICE_ROOT/STEER/STEERBase -I$ALICE_ROOT/STEER/ESD -I$ALICE_ROOT/STEER/AOD -I$ALICE_PHYSICS/OADB -I$ALICE_PHYSICS/PWGHF/base  -I$ALICE_ROOT/include -I$ALICE_ROOT/ITS -I$ALICE_ROOT/TPC -I$ALICE_ROOT/CONTAINERS -I$ALICE_ROOT/STEER -I$ALICE_ROOT/TRD -I$ALICE_ROOT/macros -I$ALICE_ROOT/ANALYSIS  -I$ALICE_PHYSICS/OADB/macros -I$ALICE_PHYSICS/PWGCF/Correlations -I$ALICE_PHYSICS/PWGCF -I$ALICE_PHYSICS/PWGCF/Correlations/Base -I$ALICE_PHYSICS/include -g");
    
    alienHandler->SetAdditionalLibs("AliAnalysisTaskhPhiCorr.cxx AliAnalysisTaskhPhiCorr.h AddTaskQA.C libpythia6.so libEGPythia6.so libAliPythia6.so libPWGHFhfe.so libCDB.so libSTEER.so libCORRFW.so libPWGflowBase.so libPWGflowTasks.so libGui.so libProof.so libMinuit.so libXMLParser.so libRAWDatabase.so libRAWDatarec.so libCDB.so libSTEERBase.so libSTEER.so libTPCbase.so libTOFbase.so libTOFrec.so libTRDbase.so libVZERObase.so libVZEROrec.so libT0base.so libT0rec.so libPWGTools.so libPWGCFCorrelationsBase.so");

  // load libraries
   LoadLibraries();
    
  alienHandler->SetAnalysisSource("AliAnalysisTaskhPhiCorr.cxx");
  //alienHandler->SetOverwriteMode();
  alienHandler->SetRunMode(mode);
  alienHandler->SetNtestFiles(2);
  alienHandler->SetAPIVersion("V1.1x");
  alienHandler->SetAliPhysicsVersion("vAN-20170425-1");
  //alienHandler->SetFileForTestMode("File_LHC12dPass1.txt");  //txt file that tells where to look for local files if launch=local
  //alienHandler->SetGridDataDir("/alice/sim/LHC10d4/");
  //alienHandler->SetDataPattern("*ESDs.root");
  alienHandler->SetGridDataDir("//alice/data/2016/LHC16q/");
  alienHandler->SetDataPattern("*/pass1_FAST/AOD/*/*AOD.root");
  //alienHandler->SetDataPattern("*/pass4/AOD/*AOD.root");
  alienHandler->SetRunPrefix("000"); // IMPORTANT! Only need for real data, comment this line out for MC data

   
//LHC12d   
    //Int_t runArray[] = {186320, 186319, 186318, 186229, 186208, 186205, 186200, 186167, 186165, 186164, 186163, 185912, 185909, 185784, 185778, 185776, 185775, 185768, 185765, 185764, 185757, 185756, 185738, 185735, 185734, 185701, 185699, 185698, 185697, 185695, 185687, 185680, 185589, 185588, 185583, 185582, 185581, 185580, 185578};
  //Int_t runArray[] = {186320, 186319, 186318, 186229, 186208, 186205, 186200, 186167, 186165, 186164, 186163, 185912, 185909, 185784, 185778, 185776, 185775, 185768, 185765, 185764, 185757, 185756, 185738, 185735, 185734, 185701, 185699, 185698, 185697, 185695, 185687, 185680, 185589, 185588, 185583, 185582, 185581, 185580, 185578, 185575, 185574, 185565, 185563, 185474, 185465, 185461, 185457, 185375, 185371, 185363, 185362, 185361, 185360, 185359, 185356, 185351, 185350, 185349, 185303, 185302, 185300, 185299, 185296, 185293, 185292, 185291, 185289, 185288, 185284, 185282, 185221, 185217, 185208, 185206, 185203, 185198, 185196, 185189};

 //LHC13b
    //Int_t runArray[] = {195483, 195482, 195481, 195480, 195479, 195478, 195391, 195389, 195351, 195346, 195344}; 

 //LHC13c
    //Int_t runArray[] = {195529, 195531, 195566, 195567, 195568, 195592, 195593, 195596, 195633, 195635, 195644, 195673, 195675, 195677};

//LHC10d4 - MC Data
    //Int_t runArray[] = {119159, 119161, 119163, 119841, 119842, 119844, 119845, 119846, 119849, 119853, 119856, 119859, 119862, 120067, 120069, 120072, 120073, 120076, 120079, 120244, 120503, 120504, 120505, 120616, 120617, 120671, 120741, 120750, 120758, 120820, 120821, 120822, 120823, 120824, 120825, 120829};
   // Int_t runArray[] = {120073}; //for testing why files were being opened but not closed

//LHC16r - 8 TeV pPb data
    //Int_t runArray[] = {266318, 266317, 266316, 266305, 266304, 266300, 266299, 266296, 266208, 266197, 266196, 266193, 266190, 266189, 266187, 266117, 266086, 266085, 266084, 266083, 266081, 266076, 266074, 266034, 265797, 265795, 265789, 265788, 265756, 265754, 265746, 265744, 265742, 265741, 265714, 265713,265709, 265705, 265701, 265700, 265698, 265697, 265696, 265607, 265596, 265594};
   //Int_t runArray[] = {266318, 266317, 266316, 266208, 266197, 266196, 266187, 265754, 265744, 265607, 265596, 265594};

//LHC16q - 5 TeV pPb data
   Int_t runArray[] = {265525, 265521, 265501, 265500, 265499, 265435, 265427, 265426, 265425, 265424, 265422, 265421, 265420, 265419, 265388, 265387, 265385, 265384, 265383, 265381, 265378, 265377, 265344, 265343, 265342, 265339, 265338, 265336, 265335, 265334, 265332, 265309};  
   for (Int_t i =  runcycle[cyclenumber - 1]; i < runcycle[cyclenumber] ; i++)
   {
    if (i == sizeof(runArray) / sizeof(runArray[1])) break;
    alienHandler->AddRunNumber(runArray[i]);
   }

   printf("\n\nSetting Up alienHandler.\n\n");
   alienHandler->SetGridWorkingDir(work_dir);
   alienHandler->SetGridOutputDir(output_dir);
   alienHandler->SetDefaultOutputs();
   alienHandler->SetAnalysisMacro("PhiInvMass.C");
   alienHandler->SetSplitMaxInputFileNumber(noffiles);
   alienHandler->SetExecutable("PhiInvMass.sh");
   alienHandler->SetExecutableCommand("aliroot -b -q");
   alienHandler->SetTTL(ttl); //10000
   alienHandler->SetInputFormat("xml-single");
   alienHandler->SetJDLName("PhiInvMass.jdl");
   alienHandler->SetPrice(1);
   alienHandler->SetSplitMode("se");
   alienHandler->SetMasterResubmitThreshold(10);
   alienHandler->SetMergeExcludes("EventStat_temp.root");
   alienHandler->SetOutputToRunNo(kTRUE);
   alienHandler->SetKeepLogs(kTRUE);
   alienHandler->SetMaxMergeFiles(15);
   alienHandler->SetMaxMergeStages(7);
   alienHandler->SetMergeViaJDL(pre_final_stage);
//    alienHandler->SetOneStageMerging(kFALSE);   //???????????????????????????????-------------------
    if (!alienHandler) return;

    
// load the necessary macros
   //gROOT->LoadMacro("$ALICE_ROOT/ANALYSIS/macros/AddTaskPhysicsSelection.C");
   //gROOT->LoadMacro("$ALICE_ROOT/ANALYSIS/macros/AddTaskPIDResponse.C");
// Use AliRoot includes to compile our task
   gROOT->ProcessLine(".include $ALICE_ROOT/include");
   gROOT->ProcessLine(".include $ALICE_ROOT/EMCAL");
   gROOT->ProcessLine(".include $ALICE_PHYSICS/include");
   gROOT->ProcessLine(".include $ALICE_ROOT/ANALYSIS/");
   gROOT->ProcessLine(".include $PWD/.");

   gSystem->AddIncludePath("-I. -I$ROOTSYS/include -I$ALICE_ROOT -I$ALICE_ROOT/EMCAL -I$ALICE_ROOT/PYTHIA6 -I$ALICE_ROOT/ANALYSIS -I$ALICE_ROOT/CONTAINERS -I$ALICE_ROOT/STEER/STEER -I$ALICE_ROOT/STEER/STEERBase -I$ALICE_ROOT/STEER/ESD -I$ALICE_ROOT/STEER/AOD -I$ALICE_PHYSICS/OADB -I$ALICE_PHYSICS/PWGHF/base  -I$ALICE_ROOT/include -I$ALICE_ROOT/ITS -I$ALICE_ROOT/TPC -I$ALICE_ROOT/CONTAINERS -I$ALICE_ROOT/STEER -I$ALICE_ROOT/TRD -I$ALICE_ROOT/macros -I$ALICE_ROOT/ANALYSIS  -I$ALICE_PHYSICS/OADB -I$ALICE_ROOT/PWG/FLOW/Base -I$ALICE_PHYSICS/PWGCF/Correlations/Base -I$ALICE_PHYSICS/PWGCF/Correlations -g ");
   // gSystem->AddIncludePath("-I. -I$ROOTSYS/include -I$ALICE_ROOT -I$ALICE_ROOT/include -I$ALICE_ROOT/macros -I$ALICE_PHYSICS/include -g");

   //printf("\n!!!!!!!!!!!!!!!!!!!!!!\n AliAnalysis Manager \n\n");
   AliAnalysisManager *mgr = new AliAnalysisManager("PhiAnalysis");
   mgr->SetGridHandler(alienHandler);

   AliAODInputHandler* aodH = new AliAODInputHandler();
   mgr->SetInputEventHandler(aodH);
//   AliESDInputHandler* esdH = new AliESDInputHandler();
//   mgr->SetInputEventHandler(esdH);

//    AliMCEventHandler* mcH = new AliMCEventHandler();
//    mgr->SetMCtruthEventHandler(mcH);   
//    mcH->SetReadTR(kFALSE);

   gROOT->LoadMacro("AddTaskQA.C");
   gROOT->LoadMacro("AliAnalysisTaskhPhiCorr.cxx++g");
   gROOT->LoadMacro("$ALICE_PHYSICS/OADB/macros/AddTaskPhysicsSelection.C");
   gROOT->LoadMacro("$ALICE_PHYSICS/OADB/COMMON/MULTIPLICITY/macros/AddTaskMultSelection.C");
   gROOT->LoadMacro("$ALICE_ROOT/ANALYSIS/macros/AddTaskPIDResponse.C");


    //switch on aliphysicsselection
    AddTaskMultSelection();
    AliPhysicsSelectionTask* physSelTask = AddTaskPhysicsSelection(kFALSE, kTRUE); 
    //Only set true for MC
    Bool_t isMC = kFALSE;
    AddTaskPIDResponse(isMC);

    //create a task
    AliAnalysisTaskhPhiCorr *task1 = AddTaskQA(0.0, 20.0);
    AliAnalysisTaskhPhiCorr *task2 = AddTaskQA(20.0, 50.0);
    AliAnalysisTaskhPhiCorr *task3 = AddTaskQA(50.0, 100.0);

   if (!mgr->InitAnalysis())
     return;

   mgr->PrintStatus();
   //fprintf(stdout, "\n!!!!!!!!!!!!!\nAbout to launch analysis... \n");
   // Start analysis in grid.
   mgr->StartAnalysis(launch);
   //printf("\n!!!!!!!!!!!!!\nDone with StartAnalysis(launch)\n");
   fflush(stdout);
}
//---------------------------------------
void LoadLibraries()
{
  gSystem->Load("libCore");
  gSystem->Load("libTree");
  gSystem->Load("libGeom");
  gSystem->Load("libVMC");
  gSystem->Load("libXMLIO");
  gSystem->Load("libPhysics");
  gSystem->Load("libXMLParser");
  gSystem->Load("libProof");
  gSystem->Load("libMinuit");
  gSystem->Load("libSTEERBase");
  gSystem->Load("libCDB");
  gSystem->Load("libRAWDatabase");
  gSystem->Load("libRAWDatarec");
  gSystem->Load("libESD");
  gSystem->Load("libAOD");
  gSystem->Load("libOADB");
  gSystem->Load("libANALYSIS");
  gSystem->Load("libANALYSISalice");
  gSystem->Load("libOADB");
    gSystem->Load("libCDB");
    gSystem->Load("libRAWDatabase");
    gSystem->Load("libSTEER");
  gSystem->Load("libTPCbase");
  gSystem->Load("libTOFbase");
  gSystem->Load("libTOFrec");
  gSystem->Load("libTRDbase");
  gSystem->Load("libVZERObase");
  gSystem->Load("libVZEROrec");
  gSystem->Load("libT0base");
  gSystem->Load("libT0rec");
  gSystem->Load("libCORRFW");

  gSystem->Load("libPWGTools");
    gSystem->Load("libPWGHFhfe");
  gSystem->Load("libPWGflowBase");
  gSystem->Load("libPWGflowTasks");
    gSystem->Load("libPWGCFCorrelationsBase.so");

  gSystem->Load("libEMCALbase.so");
  gSystem->Load("libEMCALUtils.so");
  gSystem->Load("libEMCALrec.so");
  //  gSystem->Load("libPWG4CaloCalib.so");
  gSystem->Load("libPWGCaloTrackCorrBase.so");

  gSystem->Load("libpythia6.so");

  printf("!!!!!!!!!!!!!! loaded all libraries\n\n");
  //    if(use_parFiles)
  //    {
  // //     AliAnalysisAlien::SetupPar("PWGflowBase");
  //      AliAnalysisAlien::SetupPar("PWGflowTasks");
  //    }
}

