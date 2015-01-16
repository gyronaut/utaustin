/*
 *  
 *
 *
 * */


void RunMacro()
{
//     check the function for asymmetric TPC cut in ConfigHFEemcalMod....the rest is still necessary????

   // Firstly, set some variables
   const char* launch = "grid"; // grid, local (if your data is on your local machine, doesn't connect at all)
   const char*  mode = "terminate"; //test, full, terminate  (test= connect to grid but run locally, full= run on grid, terminate= merge output on grid)
   Bool_t pre_final_stage = kTRUE; //true = merging done on grid, false = merge happens locally
   //Int_t cyclenumber = 10;    
   Int_t cyclenumber = 1;    
   Bool_t debug = kTRUE;
   char* work_dir = "Test-12-12";
   char* output_dir = "output";
   Int_t ttl = 50000;
   Int_t noffiles = 40;
   //Int_t runcycle[] = {0,36,78};
   Int_t runcycle[] = {0,1};
// load libraries
   LoadLibraries();

   Bool_t UseParfiles = kFALSE;

// create and customize the alien handler
  AliAnalysisAlien *alienHandler = new AliAnalysisAlien();
    
    
    alienHandler->AddIncludePath("-I. -I$ROOTSYS/include -I$ALICE_ROOT -I$ALICE_ROOT/EMCAL -I$ALICE_ROOT/ANALYSIS -I$ALICE_ROOT/PWGGA -I$ALICE_ROOT/PWGHF -I$ALICE_ROOT/PWGHF/hfe -I$ALICE_ROOT/CONTAINERS -I$ALICE_ROOT/STEER/STEER -I$ALICE_ROOT/STEER/STEERBase -I$ALICE_ROOT/STEER/ESD -I$ALICE_ROOT/STEER/AOD -I$ALICE_ROOT/OADB -I$ALICE_ROOT/PWGHF/base  -I$ALICE_ROOT/include -I$ALICE_ROOT/ITS -I$ALICE_ROOT/TPC -I$ALICE_ROOT/CONTAINERS -I$ALICE_ROOT/STEER -I$ALICE_ROOT/TRD -I$ALICE_ROOT/macros -I$ALICE_ROOT/ANALYSIS  -I$ALICE_ROOT/OADB  -I$ALICE_ROOT/PWGHF/hfe -I$ALICE_ROOT/PWG/EMCAL -g");
    
    alienHandler->SetAdditionalLibs("AliAnalysisTaskQA.cxx AliAnalysisTaskQA.h AddTaskQA.C libPWGHFhfe.so libCDB.so libSTEER.so libCORRFW.so libPWGflowBase.so libPWGflowTasks.so libGui.so libProof.so libMinuit.so libXMLParser.so libRAWDatabase.so libRAWDatarec.so libCDB.so libSTEERBase.so libSTEER.so libTPCbase.so libTOFbase.so libTOFrec.so libTRDbase.so libVZERObase.so libVZEROrec.so libT0base.so libT0rec.so libTENDER.so libTENDERSupplies.so libPWGTools.so libPWGEMCAL.so");
    
  if(UseParfiles){
    alienHandler->SetupPar("PWGHFhfe");
    alienHandler->EnablePackage("PWGHFhfe.par");  
  }
  alienHandler->SetAnalysisSource("AliAnalysisTaskQA.cxx");
//  alienHandler->SetOverwriteMode();
  alienHandler->SetRunMode(mode);
  alienHandler->SetNtestFiles(5);
  alienHandler->SetAPIVersion("V1.1x");
  alienHandler->SetROOTVersion("v5-34-08-6");
  alienHandler->SetAliROOTVersion("vAN-20141202");
  alienHandler->SetFileForTestMode("File_LHC12dPass1.txt");  //txt file that tells where to look for local files if launch=local
  alienHandler->SetGridDataDir("/alice/data/2012/LHC12d/");
     //alienHandler->SetDataPattern("pass1/*ESDs.root");
  alienHandler->SetDataPattern("*/pass1/*/*AOD.root");
   // alienHandler->SetGridDataDir("//alice/data/2013/LHC13g/");
   // alienHandler->SetDataPattern("*/pass1/*/*AOD.root");
  alienHandler->SetRunPrefix("000"); // IMPORTANT!

   
//LHC12d   
    //Int_t runArray[] = {186320, 186319, 186318, 186229, 186208, 186205, 186200, 186167, 186165, 186164, 186163, 185912, 185909, 185784, 185778, 185776, 185775, 185768, 185765, 185764, 185757, 185756, 185738, 185735, 185734, 185701, 185699, 185698, 185697, 185695, 185687, 185680, 185589, 185588, 185583, 185582, 185581, 185580, 185578};
  Int_t runArray[] = {186320, 186319, 186318, 186229, 186208, 186205, 186200, 186167, 186165, 186164, 186163, 185912, 185909, 185784, 185778, 185776, 185775, 185768, 185765, 185764, 185757, 185756, 185738, 185735, 185734, 185701, 185699, 185698, 185697, 185695, 185687, 185680, 185589, 185588, 185583, 185582, 185581, 185580, 185578, 185575, 185574, 185565, 185563, 185474, 185465, 185461, 185457, 185375, 185371, 185363, 185362, 185361, 185360, 185359, 185356, 185351, 185350, 185349, 185303, 185302, 185300, 185299, 185296, 185293, 185292, 185291, 185289, 185288, 185284, 185282, 185221, 185217, 185208, 185206, 185203, 185198, 185196, 185189};

 //LHC13g
    //Int_t runArray[] = {197606};
    
   for (Int_t i =  runcycle[cyclenumber - 1]; i < runcycle[cyclenumber] ; i++)
   {
    if (i == sizeof(runArray) / sizeof(runArray[1])) break;
    alienHandler->AddRunNumber(runArray[i]);
   }

   alienHandler->SetGridWorkingDir(work_dir);
   alienHandler->SetGridOutputDir(output_dir);
   alienHandler->SetDefaultOutputs();
   alienHandler->SetAnalysisMacro("ElectronV2.C");
   alienHandler->SetSplitMaxInputFileNumber(noffiles);
   alienHandler->SetExecutable("Electron.sh");
   alienHandler->SetExecutableCommand("aliroot -b -q");
   alienHandler->SetTTL(ttl); //10000
   alienHandler->SetInputFormat("xml-single");
   alienHandler->SetJDLName("Electron.jdl");
   alienHandler->SetPrice(1);
   alienHandler->SetSplitMode("se");
   alienHandler->SetMasterResubmitThreshold(10);
   alienHandler->SetMergeExcludes("EventStat_temp.root");
   alienHandler->SetOutputToRunNo(kTRUE);
   alienHandler->SetKeepLogs(kTRUE);
   alienHandler->SetMaxMergeStages(4);
   alienHandler->SetMergeViaJDL(pre_final_stage);
//    alienHandler->SetOneStageMerging(kFALSE);   //???????????????????????????????-------------------
    if (!alienHandler) return;

    
// load the necessary macros
   //gROOT->LoadMacro("$ALICE_ROOT/ANALYSIS/macros/AddTaskPhysicsSelection.C");
   //gROOT->LoadMacro("$ALICE_ROOT/ANALYSIS/macros/AddTaskPIDResponse.C");
// Use AliRoot includes to compile our task
   gROOT->ProcessLine(".include $ALICE_ROOT/include");
   gROOT->ProcessLine(".include $ALICE_ROOT/EMCAL");
   gROOT->ProcessLine(".include $ALICE_ROOT/PWGGA/");
   gROOT->ProcessLine(".include $ALICE_ROOT/ANALYSIS/");
   gROOT->ProcessLine(".include $PWD/.");
   gROOT->ProcessLine(".include $ALICE_ROOT/PWGHF");
   gROOT->ProcessLine(".include $ALICE_ROOT/PWGHF/hfe");

gSystem->AddIncludePath("-I. -I$ROOTSYS/include -I$ALICE_ROOT -I$ALICE_ROOT/EMCAL -I$ALICE_ROOT/ANALYSIS -I$ALICE_ROOT/PWGGA -I$ALICE_ROOT/PWGHF -I$ALICE_ROOT/PWGHF/hfe -I$ALICE_ROOT/CONTAINERS -I$ALICE_ROOT/STEER/STEER -I$ALICE_ROOT/STEER/STEERBase -I$ALICE_ROOT/STEER/ESD -I$ALICE_ROOT/STEER/AOD -I$ALICE_ROOT/OADB -I$ALICE_ROOT/PWGHF/base  -I$ALICE_ROOT/include -I$ALICE_ROOT/ITS -I$ALICE_ROOT/TPC -I$ALICE_ROOT/CONTAINERS -I$ALICE_ROOT/STEER -I$ALICE_ROOT/TRD -I$ALICE_ROOT/macros -I$ALICE_ROOT/ANALYSIS  -I$ALICE_ROOT/OADB -I$ALICE_ROOT/PWG/FLOW/Base -g ");
    
   AliAnalysisManager *mgr = new AliAnalysisManager("ElectronAnalysis");
   mgr->SetGridHandler(alienHandler);

   AliAODInputHandler* aodH = new AliAODInputHandler();
   mgr->SetInputEventHandler(aodH);
   // AliESDInputHandler* esdH = new AliESDInputHandler();
   // mgr->SetInputEventHandler(esdH);
   

   gROOT->LoadMacro("AddTaskQA.C");
   gROOT->LoadMacro("AliAnalysisTaskQA.cxx++g");
   gROOT->LoadMacro("$ALICE_ROOT/ANALYSIS/macros/AddTaskPhysicsSelection.C");
   gROOT->LoadMacro("$ALICE_ROOT/ANALYSIS/macros/AddTaskPIDResponse.C");
   AddTaskPIDResponse(kFALSE);

   //create a task
   AliAnalysisTaskQA *taskQA = AddTaskQA();

   if (!mgr->InitAnalysis())
     return;

   mgr->PrintStatus();
   // Start analysis in grid.
   mgr->StartAnalysis(launch);
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
  gSystem->Load("libTENDER");
  gSystem->Load("libTENDERSupplies");

    gSystem->Load("libPWGTools");
    gSystem->Load("libPWGEMCAL");
    gSystem->Load("libPWGHFhfe");
  gSystem->Load("libPWGflowBase");
  gSystem->Load("libPWGflowTasks");

  gSystem->Load("libEMCALbase.so");
  gSystem->Load("libEMCALUtils.so");
  gSystem->Load("libEMCALrec.so");
  //  gSystem->Load("libPWG4CaloCalib.so");
  gSystem->Load("libPWGCaloTrackCorrBase.so");



  //    if(use_parFiles)
  //    {
  // //     AliAnalysisAlien::SetupPar("PWGflowBase");
  //      AliAnalysisAlien::SetupPar("PWGflowTasks");
  //    }
}

