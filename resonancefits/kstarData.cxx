void kstarData(){
    TFile* input = new TFile("~/Downloads/AnalysisResults_KstarPbPb.root");
    TList* list = (TList*)input->Get("RsnOut_KstarTPC2TOFVeto3SigTrCut2011DefaultPidFix2017");
    
    TH3D* kstar = (TH3D*)list->FindObject("KStarPbPbData_CustomId1_UnlikePM");
    TH3D* kstarmixed = (TH3D*)list->FindObject("KStarPbPbData_CustomId1_MixingPM");
    TH3D* kstarbar = (TH3D*)list->FindObject("KStarPbPbData_CustomId1_UnlikeMP");
    TH3D* kstarbarmixed = (TH3D*)list->FindObject("KStarPbPbData_CustomId1_MixingPM");
}
