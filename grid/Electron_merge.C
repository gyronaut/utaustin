void Electron_merge(const char *dir, Int_t stage=0)
{
// Automatically generated merging macro executed in grid subjobs

   TStopwatch timer;
   timer.Start();

// Reset existing include path and add current directory first in the search
   gSystem->SetIncludePath("-I.");
// Load analysis framework libraries
   gSystem->Load("libANALYSIS");
   gSystem->Load("libOADB");
   gSystem->Load("libANALYSISalice");
   gSystem->Load("libCORRFW");

// include path
   TString intPath = gInterpreter->GetIncludePath();
   TObjArray *listpaths = intPath.Tokenize(" ");
   TIter nextpath(listpaths);
   TObjString *pname;
   while ((pname=(TObjString*)nextpath())) {
      TString current = pname->GetName();
      if (current.Contains("AliRoot") || current.Contains("ALICE_ROOT")) continue;
      gSystem->AddIncludePath(current);
   }
   if (listpaths) delete listpaths;
   gSystem->AddIncludePath("-I. -I$ROOTSYS/include -I$ALICE_ROOT -I$ALICE_ROOT/EMCAL -I$ALICE_ROOT/ANALYSIS -I$ALICE_ROOT/PWGGA -I$ALICE_ROOT/PWGHF -I$ALICE_ROOT/PWGHF/hfe -I$ALICE_ROOT/CONTAINERS -I$ALICE_ROOT/STEER/STEER -I$ALICE_ROOT/STEER/STEERBase -I$ALICE_ROOT/STEER/ESD -I$ALICE_ROOT/STEER/AOD -I$ALICE_ROOT/OADB -I$ALICE_ROOT/PWGHF/base  -I$ALICE_ROOT/include -I$ALICE_ROOT/ITS -I$ALICE_ROOT/TPC -I$ALICE_ROOT/CONTAINERS -I$ALICE_ROOT/STEER -I$ALICE_ROOT/TRD -I$ALICE_ROOT/macros -I$ALICE_ROOT/ANALYSIS  -I$ALICE_ROOT/OADB  -I$ALICE_ROOT/PWGHF/hfe -I$ALICE_ROOT/PWG/EMCAL -g ");
   gROOT->ProcessLine(".include $ALICE_ROOT/include");
   printf("Include path: %s\n", gSystem->GetIncludePath());

// Add aditional AliRoot libraries
   gSystem->Load("libPWGHFhfe");
   gSystem->Load("libCDB");
   gSystem->Load("libSTEER");
   gSystem->Load("libCORRFW");
   gSystem->Load("libPWGflowBase");
   gSystem->Load("libPWGflowTasks");
   gSystem->Load("libGui");
   gSystem->Load("libProof");
   gSystem->Load("libMinuit");
   gSystem->Load("libXMLParser");
   gSystem->Load("libRAWDatabase");
   gSystem->Load("libRAWDatarec");
   gSystem->Load("libCDB");
   gSystem->Load("libSTEERBase");
   gSystem->Load("libSTEER");
   gSystem->Load("libTPCbase");
   gSystem->Load("libTOFbase");
   gSystem->Load("libTOFrec");
   gSystem->Load("libTRDbase");
   gSystem->Load("libVZERObase");
   gSystem->Load("libVZEROrec");
   gSystem->Load("libT0base");
   gSystem->Load("libT0rec");
   gSystem->Load("libTENDER");
   gSystem->Load("libTENDERSupplies");
   gSystem->Load("libPWGTools");
   gSystem->Load("libPWGEMCAL");

// Analysis source to be compiled at runtime (if any)
   gROOT->ProcessLine(".L AliAnalysisTaskQA.cxx+g");

// Connect to AliEn
   if (!TGrid::Connect("alien://")) return;
// Set temporary merging directory to current one
   gSystem->Setenv("TMPDIR", gSystem->pwd());

// Set temporary compilation directory to current one
   gSystem->SetBuildDir(gSystem->pwd(), kTRUE);

   TString outputDir = dir;
   TString outputFiles = "AnalysisResults.root";
   TString mergeExcludes = "EventStat_temp.root ";
   TObjArray *list = outputFiles.Tokenize(",");
   TIter *iter = new TIter(list);
   TObjString *str;
   TString outputFile;
   Bool_t merged = kTRUE;
   AliAnalysisManager *mgr = AliAnalysisAlien::LoadAnalysisManager("Electron.root");
   if (!mgr) {
      printf("ERROR: Analysis manager could not be extracted from file ");
      return;
   }
   while((str=(TObjString*)iter->Next())) {
      outputFile = str->GetString();
      if (outputFile.Contains("*")) continue;
      Int_t index = outputFile.Index("@");
      if (index > 0) outputFile.Remove(index);
      // Skip already merged outputs
      if (!gSystem->AccessPathName(outputFile)) {
         printf("Output file <%s> found. Not merging again.\n",outputFile.Data());
         continue;
      }
      if (mergeExcludes.Contains(outputFile.Data())) continue;
      merged = AliAnalysisAlien::MergeOutput(outputFile, outputDir, 100, stage);
      if (!merged) {
         printf("ERROR: Cannot merge %s\n", outputFile.Data());
         return;
      }
   }
   // all outputs merged, validate
   ofstream out;
   out.open("outputs_valid", ios::out);
   out.close();
   // read the analysis manager from file
   if (!outputDir.Contains("Stage")) return;
   mgr->SetRunFromPath(mgr->GetRunFromAlienPath(dir));
   mgr->SetSkipTerminate(kFALSE);
   mgr->PrintStatus();
   AliLog::SetGlobalLogLevel(AliLog::kError);
   TTree *tree = NULL;
   mgr->StartAnalysis("gridterminate", tree);
}

