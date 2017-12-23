/*
 *  
 *
 *
 * */
#include "TRoot.h"
#include "TRint.h"
#include "TSystem.h"

void TimeCalibRunMacro()
{
//     check the function for asymmetric TPC cut in ConfigHFEemcalMod....the rest is still necessary????

   // Firstly, set some variables
   const char* launch = "grid"; // grid, local (if your data is on your local machine, doesn't connect at all)
   const char*  mode = "terminate"; //test, full, terminate  (test= connect to grid but run locally, full= run on grid, terminate= merge output on grid)
   Bool_t pre_final_stage = kTRUE; //TRUE = merging done on grid, FALSE = merge happens locally   
   Int_t cyclenumber = 1;
   Bool_t debug = kTRUE;
   char* work_dir = "TimeCalibWork";
   char* output_dir = "LHC16k_check";
   Int_t ttl = 50000;
   Int_t noffiles = 40;
   //Int_t runcycle[]={0,10,40,60,80,95,110,125,140,155,170};
   Int_t runcycle[] ={0,70,140};
   Bool_t UseParfiles = kFALSE;

// create and customize the alien handler
  AliAnalysisAlien *alienHandler = new AliAnalysisAlien();
    
 // load libraries
   LoadLibraries();
   
    alienHandler->AddIncludePath("-I. -I$ROOTSYS/include -I$ALICE_ROOT -I$ALICE_ROOT/EMCAL -I$ALICE_ROOT/PYTHIA6 -I$ALICE_ROOT/ANALYSIS -I$ALICE_PHYSICS/PWGGA -I$ALICE_PHYSICS/PWGHF -I$ALICE_PHYSICS/PWGHF/hfe -I$ALICE_ROOT/CONTAINERS -I$ALICE_ROOT/STEER/STEER -I$ALICE_ROOT/STEER/STEERBase -I$ALICE_ROOT/STEER/ESD -I$ALICE_ROOT/STEER/AOD -I$ALICE_PHYSICS/OADB -I$ALICE_PHYSICS/PWGHF/base  -I$ALICE_ROOT/include -I$ALICE_ROOT/ITS -I$ALICE_ROOT/TPC -I$ALICE_ROOT/CONTAINERS -I$ALICE_ROOT/STEER -I$ALICE_ROOT/TRD -I$ALICE_ROOT/macros -I$ALICE_ROOT/ANALYSIS  -I$ALICE_PHYSICS/OADB/macros  -I$ALICE_PHYSICS/PWGHF/hfe -I$ALICE_PHYSICS/PWG/EMCAL -g");
    
    alienHandler->SetAdditionalLibs("AliAnalysisTaskEMCALTimeCalibJB.cxx AliAnalysisTaskEMCALTimeCalibJB.h AddTaskEMCALTimeCalibrationJB.C libpythia6.so libEGPythia6.so libAliPythia6.so libPWGHFhfe.so libCDB.so libSTEER.so libCORRFW.so libPWGflowBase.so libPWGflowTasks.so libGui.so libProof.so libMinuit.so libXMLParser.so libRAWDatabase.so libRAWDatarec.so libCDB.so libSTEERBase.so libSTEER.so libTPCbase.so libTOFbase.so libTOFrec.so libTRDbase.so libVZERObase.so libVZEROrec.so libT0base.so libT0rec.so libTENDER.so libTENDERSupplies.so libPWGTools.so libPWGEMCAL.so");
    
  if(UseParfiles){
    alienHandler->SetupPar("PWGHFhfe");
    alienHandler->EnablePackage("PWGHFhfe.par");  
  }

// Trying to add new PHYSICS package
//  alienHandler->AddExternalPackage("AliPhysics::vAN-20160606-1");

  alienHandler->SetAnalysisSource("AliAnalysisTaskEMCALTimeCalibJB.cxx");
  //alienHandler->SetOverwriteMode();
  alienHandler->SetRunMode(mode);
  alienHandler->SetNtestFiles(5);
  alienHandler->SetAPIVersion("V1.1x");
  //alienHandler->SetROOTVersion("v5-34-30-alice5-2");
  //alienHandler->SetAliROOTVersion("v5-08-16-1");
  alienHandler->SetAliPhysicsVersion("vAN-20171125-1");
  //alienHandler->SetFileForTestMode("File_LHC12dPass1.txt");  //txt file that tells where to look for local files if launch=local
  //alienHandler->SetGridDataDir("/alice/sim/LHC10d4/");
  //alienHandler->SetDataPattern("*ESDs.root");
  //alienHandler->SetDataPattern("*/pass1/*/*AOD.root");
  alienHandler->SetGridDataDir("//alice/data/2016/LHC16k/");
  alienHandler->SetDataPattern("*/muon_calo_pass1/*/*ESDs.root");
  alienHandler->SetRunPrefix("000"); // IMPORTANT! Only need for real data, comment this line out for MC data

   
//LHC12d   
    //Int_t runArray[] = {186320, 186319, 186318, 186229, 186208, 186205, 186200, 186167, 186165, 186164, 186163, 185912, 185909, 185784, 185778, 185776, 185775, 185768, 185765, 185764, 185757, 185756, 185738, 185735, 185734, 185701, 185699, 185698, 185697, 185695, 185687, 185680, 185589, 185588, 185583, 185582, 185581, 185580, 185578};
  //Int_t runArray[] = {186320, 186319, 186318, 186229, 186208, 186205, 186200, 186167, 186165, 186164, 186163, 185912, 185909, 185784, 185778, 185776, 185775, 185768, 185765, 185764, 185757, 185756, 185738, 185735, 185734, 185701, 185699, 185698, 185697, 185695, 185687, 185680, 185589, 185588, 185583, 185582, 185581, 185580, 185578, 185575, 185574, 185565, 185563, 185474, 185465, 185461, 185457, 185375, 185371, 185363, 185362, 185361, 185360, 185359, 185356, 185351, 185350, 185349, 185303, 185302, 185300, 185299, 185296, 185293, 185292, 185291, 185289, 185288, 185284, 185282, 185221, 185217, 185208, 185206, 185203, 185198, 185196, 185189};

 //LHC15n_mcp1
 //   Int_t runArray[] = {244628, 244627, 244626, 244619, 244618, 244617, 244542, 244540, 244531, 244484, 244483, 244482, 244481, 244480, 244456, 244453, 244421, 244418, 244416, 244411, 244377, 244364, 244359, 244355, 244351, 244340};

//LHC15j_mcp2
    //Int_t runArray[] = {237050};

//LHC10d4 - MC Data
    //Int_t runArray[] = {119159, 119161, 119163, 119841, 119842, 119844, 119845, 119846, 119849, 119853, 119856, 119859, 119862, 120067, 120069, 120072, 120073, 120076, 120079, 120244, 120503, 120504, 120505, 120616, 120617, 120671, 120741, 120750, 120758, 120820, 120821, 120822, 120823, 120824, 120825, 120829};
   // Int_t runArray[] = {120073}; //for testing why files were being opened but not closed

//LHC17n Xe-Xe run
   // Int_t runArray[] = {280234, 280235};

  //LHC16k
  Int_t runArray[] = {258048, 258049, 257026, 257539, 257540, 257541, 257537, 258059, 256514, 258062, 258063, 257560, 257561, 257562, 257563, 257564, 257566, 256506, 256552, 256554, 256556, 256560, 256561, 256562, 256564, 257077, 257590, 256567, 257080, 256692, 257594, 258107, 258108, 258109, 258113, 258114, 257092, 257605, 257606, 257100, 256589, 256591, 256592, 256697, 257635, 257642, 256619, 256620, 257136, 257137, 257138, 257139, 257140, 257141, 257142, 257144, 257145, 258178, 257474, 257682, 258197, 258198, 257687, 257689, 258202, 258203, 258204, 257694, 257697, 256676, 256677, 256681, 256684, 256694, 256691, 257204, 257206, 256695, 257209, 257724, 257733, 257734, 257735, 257224, 257737, 258256, 258257, 258258, 257754, 258270, 258271, 258273, 258274, 257765, 258278, 258280, 257260, 257773, 258299, 258301, 258302, 258303, 258306, 258307, 257797, 257798, 257799, 257800, 257803, 257318, 257320, 257322, 257587, 258359, 257850, 257855, 256941, 258387, 258388, 256512, 258393, 257082, 257083, 257892, 257893, 257084, 256658, 257912, 258426, 256565, 257936, 257937, 257939, 258454, 258456, 257433, 258117, 257691, 257958, 257960, 257692, 257963, 258477, 256942, 256944, 257457, 258498, 258499, 257487, 257490, 257491, 257492, 258012, 257530, 258014, 258017, 258019, 258537, 257021, 257011, 257012, 256504, 257364, 258042, 257531, 258045, 256510};

   for (Int_t i =  runcycle[cyclenumber - 1]; i < runcycle[cyclenumber] ; i++)
   {
    if (i == sizeof(runArray) / sizeof(runArray[1])) break;
    alienHandler->AddRunNumber(runArray[i]);
   }

   printf("\n\nSetting Up alienHandler.\n\n");
   alienHandler->SetGridWorkingDir(work_dir);
   alienHandler->SetGridOutputDir(output_dir);
   alienHandler->SetDefaultOutputs();
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

   gROOT->LoadMacro("AliAnalysisTaskEMCALTimeCalibJB.cxx++g");
   gROOT->LoadMacro("AddTaskEMCALTimeCalibrationJB.C");
   gROOT->LoadMacro("$ALICE_PHYSICS/OADB/macros/AddTaskPhysicsSelection.C");
//   gROOT->LoadMacro("$ALICE_ROOT/ANALYSIS/macros/AddTaskPIDResponse.C");


    //switch on aliphysicsselection
//    AliPhysicsSelectionTask* physSelTask = AddTaskPhysicsSelection(kTRUE); 
 
    //Only set true for MC
    Bool_t isMC = kFALSE;
//   AddTaskPIDResponse(isMC);
    //create a task
    if(!TGrid::Connect("alien://")) return;
   AliAnalysisTaskEMCALTimeCalibJB *task = AddTaskEMCALTimeCalibrationJB("AnalysisResults.root","",0.9,500,2,200,0.1,4.,0.1,4.,0.025,0.4,-20., 20., kFALSE, "alien:///alice/cern.ch/user/a/amatyja/TimeCalibRef/Reference_LHC16k.root", "alien:///alice/cern.ch/user/a/amatyja/TimeCalibRef/ReferenceSM_LHC16k_pass1_v2.root",kFALSE, kFALSE, 1, "");

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

