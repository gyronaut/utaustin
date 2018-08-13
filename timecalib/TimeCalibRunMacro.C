/*
 *  
 *
 *
 * */
#include "TRoot.h"
#include "TRint.h"
#include "TSystem.h"
void LoadLibraries();

void TimeCalibRunMacro()
{
//     check the function for asymmetric TPC cut in ConfigHFEemcalMod....the rest is still necessary????

   // Firstly, set some variables
   const char* launch = "grid"; // grid, local (if your data is on your local machine, doesn't connect at all)
   const char*  mode = "full"; //test, full, terminate  (test= connect to grid but run locally, full= run on grid, terminate= merge output on grid)
   Bool_t pre_final_stage = kTRUE; //TRUE = merging done on grid, FALSE = merge happens locally   
   Int_t cyclenumber = 2;
   Bool_t debug = kTRUE;
   char* work_dir = "TimeCalibWork";
   char* output_dir = "LHC18f_PARtest2";
   Int_t ttl = 50000;
   Int_t noffiles = 40;
   //Int_t runcycle[]={0,20,35,50,65,80,95,110,125,140,158};
   Int_t runcycle[] ={0,5,13};
   Bool_t UseParfiles = kFALSE;

// create and customize the alien handler
  AliAnalysisAlien *alienHandler = new AliAnalysisAlien();
    
 // load libraries
   LoadLibraries();
   
    alienHandler->AddIncludePath("-I. -I$ROOTSYS/include -I$ALICE_ROOT -I$ALICE_ROOT/EMCAL -I$ALICE_ROOT/PYTHIA6 -I$ALICE_ROOT/ANALYSIS -I$ALICE_PHYSICS/PWGGA -I$ALICE_PHYSICS/PWGHF -I$ALICE_PHYSICS/PWGHF/hfe -I$ALICE_ROOT/CONTAINERS -I$ALICE_ROOT/STEER/STEER -I$ALICE_ROOT/STEER/STEERBase -I$ALICE_ROOT/STEER/ESD -I$ALICE_ROOT/STEER/AOD -I$ALICE_PHYSICS/OADB -I$ALICE_PHYSICS/PWGHF/base  -I$ALICE_ROOT/include -I$ALICE_ROOT/ITS -I$ALICE_ROOT/TPC -I$ALICE_ROOT/CONTAINERS -I$ALICE_ROOT/STEER -I$ALICE_ROOT/TRD -I$ALICE_ROOT/macros -I$ALICE_ROOT/ANALYSIS  -I$ALICE_PHYSICS/OADB/macros  -I$ALICE_PHYSICS/PWGHF/hfe -I$ALICE_PHYSICS/PWG/EMCAL -g");
    
    alienHandler->SetAdditionalLibs("AliAnalysisTaskEMCALTimeCalibPAR.cxx AliAnalysisTaskEMCALTimeCalibPAR.h AddTaskEMCALTimeCalibrationPAR.C libpythia6.so libEGPythia6.so libAliPythia6.so libPWGHFhfe.so libCDB.so libSTEER.so libCORRFW.so libPWGflowBase.so libPWGflowTasks.so libGui.so libProof.so libMinuit.so libXMLParser.so libRAWDatabase.so libRAWDatarec.so libCDB.so libSTEERBase.so libSTEER.so libTPCbase.so libTOFbase.so libTOFrec.so libTRDbase.so libVZERObase.so libVZEROrec.so libT0base.so libT0rec.so libTENDER.so libTENDERSupplies.so libPWGTools.so libPWGEMCAL.so");
    
  if(UseParfiles){
    alienHandler->SetupPar("PWGHFhfe");
    alienHandler->EnablePackage("PWGHFhfe.par");  
  }

// Trying to add new PHYSICS package
//  alienHandler->AddExternalPackage("AliPhysics::vAN-20160606-1");

  alienHandler->SetAnalysisSource("AliAnalysisTaskEMCALTimeCalibPAR.cxx");
  //alienHandler->SetOverwriteMode();
  alienHandler->SetRunMode(mode);
  alienHandler->SetNtestFiles(5);
  alienHandler->SetAPIVersion("V1.1x");
  //alienHandler->SetROOTVersion("v5-34-30-alice5-2");
  //alienHandler->SetAliROOTVersion("v5-08-16-1");
  alienHandler->SetAliPhysicsVersion("vAN-20180809-1");
  //alienHandler->SetFileForTestMode("File_LHC12dPass1.txt");  //txt file that tells where to look for local files if launch=local
  //alienHandler->SetGridDataDir("/alice/sim/LHC10d4/");
  //alienHandler->SetDataPattern("*ESDs.root");
  //alienHandler->SetDataPattern("*/pass1/*/*AOD.root");
  alienHandler->SetGridDataDir("//alice/data/2018/LHC18f/");
  alienHandler->SetDataPattern("*/muon_calo_pass1/*/*ESDs.root");
  alienHandler->SetRunPrefix("000"); // IMPORTANT! Only need for real data, comment this line out for MC data

  //LHC18c PAR test
  //Int_t runArray[] = {285892, 285756};
  
  //LHC18d PAR test
  //Int_t runArray[] = {286284, 286313, 286341, 286350};

  //LHC18e PAR test
  //Int_t runArray[] = {286427, 286594, 286877};

  //LHC18f PAR test
  Int_t runArray[] = {287072, 287201, 287203, 287208, 287325, 287349, 287518, 287524, 287578, 287658, 287613, 287783, 287784};

   for (Int_t i =  runcycle[cyclenumber - 1]; i < runcycle[cyclenumber] ; i++)
   {
    if (i == sizeof(runArray) / sizeof(runArray[1])) break;
    alienHandler->AddRunNumber(runArray[i]);
   }

   printf("\n\nSetting Up alienHandler.\n\n");
   alienHandler->SetGridWorkingDir(work_dir);
   alienHandler->SetGridOutputDir(output_dir);
   alienHandler->SetDefaultOutputs(kTRUE);
   alienHandler->SetAnalysisMacro("MacroTimeCalib.C");
   alienHandler->SetSplitMaxInputFileNumber(noffiles);
   alienHandler->SetExecutable("ScriptTimeCalib.sh");
   alienHandler->SetExecutableCommand("aliroot -b -q");
   alienHandler->SetTTL(ttl); //10000
   alienHandler->SetInputFormat("xml-single");
   alienHandler->SetJDLName("TimeCalib.jdl");
   alienHandler->SetPrice(1);
   alienHandler->SetSplitMode("se");
   alienHandler->SetMasterResubmitThreshold(10);
   alienHandler->SetMergeExcludes("EventStat_temp.root");
   alienHandler->SetOutputToRunNo(kTRUE);
   alienHandler->SetKeepLogs(kTRUE);
   alienHandler->SetMaxMergeFiles(20);
//   alienHandler->SetMaxMergeStages(5);
   alienHandler->SetMergeViaJDL(pre_final_stage);
//    alienHandler->SetOneStageMerging(kFALSE);   //???????????????????????????????-------------------
    if (!alienHandler) return;

    
// load the necessary macros
   //gROOT->LoadMacro("$ALICE_ROOT/ANALYSIS/macros/AddTaskPhysicsSelection.C");
   //gROOT->LoadMacro("$ALICE_ROOT/ANALYSIS/macros/AddTaskPIDResponse.C");
// Use AliRoot includes to compile our task
   gROOT->ProcessLine(".include $ALICE_ROOT/include");
   gROOT->ProcessLine(".include $ALICE_ROOT/EMCAL");
   gROOT->ProcessLine(".include $ALICE_PHYSICS/PWGGA/");
   gROOT->ProcessLine(".include $ALICE_ROOT/ANALYSIS/");
   gROOT->ProcessLine(".include $PWD/.");
   gROOT->ProcessLine(".include $ALICE_PHYSICS/PWGHF");
   gROOT->ProcessLine(".include $ALICE_PHYSICS/PWGHF/hfe");

gSystem->AddIncludePath("-I. -I$ROOTSYS/include -I$ALICE_ROOT -I$ALICE_ROOT/EMCAL -I$ALICE_ROOT/PYTHIA6 -I$ALICE_ROOT/ANALYSIS -I$ALICE_PHYSICS/PWGGA -I$ALICE_PHYSICS/PWGHF -I$ALICE_PHYSICS/PWGHF/hfe -I$ALICE_ROOT/CONTAINERS -I$ALICE_ROOT/STEER/STEER -I$ALICE_ROOT/STEER/STEERBase -I$ALICE_ROOT/STEER/ESD -I$ALICE_ROOT/STEER/AOD -I$ALICE_PHYSICS/OADB -I$ALICE_PHYSICS/PWGHF/base  -I$ALICE_ROOT/include -I$ALICE_ROOT/ITS -I$ALICE_ROOT/TPC -I$ALICE_ROOT/CONTAINERS -I$ALICE_ROOT/STEER -I$ALICE_ROOT/TRD -I$ALICE_ROOT/macros -I$ALICE_ROOT/ANALYSIS  -I$ALICE_PHYSICS/OADB -I$ALICE_ROOT/PWG/FLOW/Base -g ");

   //printf("\n!!!!!!!!!!!!!!!!!!!!!!\n AliAnalysis Manager \n\n");
   AliAnalysisManager *mgr = new AliAnalysisManager("TimeCalib");
   mgr->SetGridHandler(alienHandler);

//   AliAODInputHandler* aodH = new AliAODInputHandler();
//   mgr->SetInputEventHandler(aodH);
   AliESDInputHandler* esdH = new AliESDInputHandler();
   mgr->SetInputEventHandler(esdH);

//    AliMCEventHandler* mcH = new AliMCEventHandler();
//    mgr->SetMCtruthEventHandler(mcH);   
//    mcH->SetReadTR(kFALSE);

   //gROOT->LoadMacro("AliAnalysisTaskEMCALTimeCalibJB.cxx++g");
   //gROOT->LoadMacro("AddTaskEMCALTimeCalibrationJB.C");
   //gROOT->LoadMacro("$ALICE_PHYSICS/OADB/macros/AddTaskPhysicsSelection.C");
//   gROOT->LoadMacro("$ALICE_ROOT/ANALYSIS/macros/AddTaskPIDResponse.C");


    //switch on aliphysicsselection
    AliPhysicsSelectionTask* physSelTask = reinterpret_cast<AliPhysicsSelectionTask*>(gInterpreter->ProcessLine(Form(".x %s(kFALSE, kTRUE)", gSystem->ExpandPathName("$ALICE_PHYSICS/OADB/macros/AddTaskPhysicsSelection.C")))); 
 
    //Only set true for MC
    Bool_t isMC = kFALSE;
//   AddTaskPIDResponse(isMC);
    //create a task
    if(!TGrid::Connect("alien://")) return;
    //AliAnalysisTaskEMCALTimeCalibPAR *task = reinterpret_cast<AliAnalysisTaskEMCALTimeCalibPAR*>(gInterpreter->ProcessLine(Form(".x %s(\"AnalysisResults.root\",\"\",0.01,500,2,200,0.01,100.,0.01,100.,0.025,0.01,-20., 20., kFALSE, \"\", \"\",kFALSE, kFALSE, 1, \"\", \"alien::///alice/cern.ch/user/j/jblair/TimeCalibRef/LHC18c_PARs.txt\")", gSystem->ExpandPathName("AddTaskEMCALTimeCalibrationPAR.C"))));
AliAnalysisTaskEMCALTimeCalibPAR *task = reinterpret_cast<AliAnalysisTaskEMCALTimeCalibPAR*>(gInterpreter->ProcessLine(Form(".x %s(\"AnalysisResults.root\",\"\",0.01,500,2,200,0.01,100.,0.01,100.,0.025,0.01,-20., 20., kFALSE, \"\", \"\",kFALSE, kTRUE, 1, \"\", \"LHC18f_PARs.txt\")", gSystem->ExpandPathName("AddTaskEMCALTimeCalibrationPAR.C"))));

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

  gSystem->Load("libpythia6.so");

  //    if(use_parFiles)
  //    {
  // //     AliAnalysisAlien::SetupPar("PWGflowBase");
  //      AliAnalysisAlien::SetupPar("PWGflowTasks");
  //    }
}

