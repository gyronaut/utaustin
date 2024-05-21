#include <iostream>
#include <fstream>


TString varheader(TString name){
    TString res(
Form("%s%s%s", R"X(    - header:
        name: )X", name.Data(), "\n"));
   return res; 
}

TString indentstr(int n, string s){
    TString res;
    TString str(s);
    for(int i =0; i < n; i++){
        res += "    ";
    }
    return res+str;
}

void hepdataparser(){

    TString depvar("dependent_variables:\n");
    TString indepvar("independent_variables:\n");

    TString var = varheader("yields");

    TString qualifiers( 
Form(R"X(      qualifiers:
        - name: REACTION
          value: PPB
        - name: SQRT(S)
          units: GEV
          value: 5020.0
        - name: ETARAP
          value: -0.8 to 0.8
        - name: TRIGPT
          value: 4.0 to 8.0 GeV/c
        - name: ASSOCPT
          value: 1.5 to 2.5 GeV/c
)X"));

//    printf("%s%s%s", depvar.Data(), var.Data(), qualifiers.Data());

    TFile* lowptfile = new TFile("./plotcode/ratiocode/fitsyst_cr2lowtest.root");
//  jet-yield graphs   
    TGraphErrors* low_near_phi = (TGraphErrors*)lowptfile->Get("yieldsNear");
    TGraphErrors* low_near_phi_v2 = (TGraphErrors*)lowptfile->Get("yieldsNearv2");
    TGraphErrors* low_away_phi = (TGraphErrors*)lowptfile->Get("yieldsAway");
    TGraphErrors* low_away_phi_v2 = (TGraphErrors*)lowptfile->Get("yieldsAwayv2");
    TGraphErrors* low_near_h = (TGraphErrors*)lowptfile->Get("yieldshhNear");
    TGraphErrors* low_near_h_v2 = (TGraphErrors*)lowptfile->Get("yieldshhNearv2");
    TGraphErrors* low_away_h = (TGraphErrors*)lowptfile->Get("yieldshhAway");
    TGraphErrors* low_away_h_v2 = (TGraphErrors*)lowptfile->Get("yieldshhAwayv2");
// ratio graphs
    TGraphErrors* low_near_ratio = (TGraphErrors*)lowptfile->Get("nchratiosNear");
    TGraphErrors* low_near_ratio_syst = (TGraphErrors*)lowptfile->Get("nchratiosNearSyst");
    TGraphErrors* low_near_ratio_v2 = (TGraphErrors*)lowptfile->Get("nearv2fly");
    
    TGraphErrors* low_away_ratio = (TGraphErrors*)lowptfile->Get("nchratiosAway");
    TGraphErrors* low_away_ratio_syst = (TGraphErrors*)lowptfile->Get("nchratiosAwaySyst");
    TGraphErrors* low_away_ratio_v2 = (TGraphErrors*)lowptfile->Get("awayv2fly");
    
    TGraphErrors* low_ue_ratio = (TGraphErrors*)lowptfile->Get("nchratiosBulk");
    TGraphErrors* low_ue_ratio_syst = (TGraphErrors*)lowptfile->Get("nchratiosBulkSyst");

    TGraphErrors* low_tot_ratio = (TGraphErrors*)lowptfile->Get("nchratiosTot");
    TGraphErrors* low_tot_ratio_syst = (TGraphErrors*)lowptfile->Get("nchratiosTotSyst");

    
    TString low_yields_indep, low_yields_dep;
    low_yields_indep += indepvar + varheader(R"($\langle N_{\mathrm{ch}} \rangle_{|\eta|<0.5}$)");
    
// h-phi jet yields
    
    low_yields_dep += depvar + varheader(R"(Per-trigger h--$\phi$ yields in near-side jet)")+qualifiers;
    TString low_yields_yval = indentstr(1, "  values:\n");
    TString low_yields_xval = indentstr(1, "  values:\n");

    for(int i =  0; i<3; i++){
        TString yval(Form("%s%e\n%s\n%s\n%s%e\n%s\n%s%e}\n",indentstr(2, "- value: ").Data(), low_near_phi->GetPointY(i)/250.0, indentstr(2, "  errors:").Data(), indentstr(3, "- label: stat").Data(), indentstr(3, "  symerror: ").Data(), low_near_phi->GetErrorY(i)/250.0, indentstr(3, "- label: sys,v2").Data(), indentstr(3, "  asymerror: {plus: 0.0, minus: ").Data(), -(low_near_phi->GetPointY(i)-low_near_phi_v2->GetPointY(i))/250.0));
        low_yields_yval += yval;
        TString xval(Form("%s%e\n",indentstr(2, "- value: ").Data(), low_near_phi->GetPointX(i)));
        low_yields_xval += xval;
    }
    low_yields_dep += low_yields_yval;
    low_yields_indep += low_yields_xval;
    
    low_yields_dep += varheader(R"(Per-trigger h--$\phi$ yields in away-side jet)")+qualifiers;
    low_yields_yval = indentstr(1, "  values:\n");
    //low_yields_xval = indentstr(1, "  values:\n");

    for(int i =  0; i<3; i++){
        TString yval(Form("%s%e\n%s\n%s\n%s%e\n%s\n%s%e}\n",indentstr(2, "- value: ").Data(), low_away_phi->GetPointY(i)/250.0, indentstr(2, "  errors:").Data(), indentstr(3, "- label: stat").Data(), indentstr(3, "  symerror: ").Data(), low_away_phi->GetErrorY(i)/250.0, indentstr(3, "- label: sys,v2").Data(), indentstr(3, "  asymerror: {plus: 0.0, minus: ").Data(), -(low_away_phi->GetPointY(i)-low_away_phi_v2->GetPointY(i))/250.0));
        low_yields_yval += yval;
        //TString xval(Form("%s%e\n",indentstr(2, "- value: ").Data(), low_away_phi->GetPointX(i)));
        //low_yields_xval += xval;
    }
    low_yields_dep += low_yields_yval;
    //low_yields_indep += low_yields_xval;

//  h-h jet yields
    low_yields_dep += varheader(R"(Per-trigger h--h yields in near-side jet)")+qualifiers;
    low_yields_yval = indentstr(1, "  values:\n");
    //low_yields_xval = indentstr(1, "  values:\n");

    for(int i =  0; i<3; i++){
        TString yval(Form("%s%e\n%s\n%s\n%s%e\n%s\n%s%e}\n",indentstr(2, "- value: ").Data(), low_near_h->GetPointY(i), indentstr(2, "  errors:").Data(), indentstr(3, "- label: stat").Data(), indentstr(3, "  symerror: ").Data(), low_near_h->GetErrorY(i), indentstr(3, "- label: sys,v2").Data(), indentstr(3, "  asymerror: {plus: 0.0, minus: ").Data(), -(low_near_h->GetPointY(i)-low_near_h_v2->GetPointY(i))));
        low_yields_yval += yval;
        //TString xval(Form("%s%e\n",indentstr(2, "- value: ").Data(), low_near_h->GetPointX(i)));
        //low_yields_xval += xval;
    }
    low_yields_dep += low_yields_yval;
    //low_yields_indep += low_yields_xval;
    
    low_yields_dep += varheader(R"(Per-trigger h--h yields in away-side jet)")+qualifiers;
    low_yields_yval = indentstr(1, "  values:\n");
    //low_yields_xval = indentstr(1, "  values:\n");

    for(int i =  0; i<3; i++){
        TString yval(Form("%s%e\n%s\n%s\n%s%e\n%s\n%s%e}\n",indentstr(2, "- value: ").Data(), low_away_h->GetPointY(i), indentstr(2, "  errors:").Data(), indentstr(3, "- label: stat").Data(), indentstr(3, "  symerror: ").Data(), low_away_h->GetErrorY(i), indentstr(3, "- label: sys,v2").Data(), indentstr(3, "  asymerror: {plus: 0.0, minus: ").Data(), -(low_away_h->GetPointY(i)-low_away_h_v2->GetPointY(i))));
        low_yields_yval += yval;
        //TString xval(Form("%s%e\n",indentstr(2, "- value: ").Data(), low_away_h->GetPointX(i)));
        //low_yields_xval += xval;
    }
    low_yields_dep += low_yields_yval;
    //low_yields_indep += low_yields_xval;

   
    TString low_jet_output = low_yields_indep + low_yields_dep;
 
    // write low pt jet-yield table
    FILE *f = fopen("lowpt_jet_data.yaml", "w+");    
    fprintf(f, "%s", low_jet_output.Data());
    fclose(f);

    
    // low pt ratio data
    TString low_ratio_indep, low_ratio_dep, low_ratio_yval, low_ratio_xval;
    low_ratio_indep += indepvar + varheader(R"($\langle N_{\mathrm{ch}} \rangle_{|\eta|<0.5}$)");
    // near side
    low_ratio_dep += depvar + varheader(R"($1/{N}_{\mathrm{trig}} \mathrm{d}N_{\mathrm{assoc}}/\mathrm{d}{\Delta\varphi} near-side ratio \biggl(\frac{\mathrm{h--}\phi}{\mathrm{h--h}\biggr)$)")+qualifiers;

    low_ratio_yval = indentstr(1, "  values:\n");
    low_ratio_xval = indentstr(1, "  values:\n");


    for(int i =  0; i<3; i++){
        TString yval(Form("%s%e\n%s\n%s\n%s%e\n%s\n%s%e}\n%s\n%s%e\n",indentstr(2, "- value: ").Data(), low_near_ratio->GetPointY(i), indentstr(2, "  errors:").Data(), indentstr(3, "- label: stat").Data(), indentstr(3, "  symerror: ").Data(), low_near_ratio->GetErrorY(i), indentstr(3, "- label: sys,v2").Data(), indentstr(3, "  asymerror: {plus: 0.0, minus: ").Data(), -2.0*(abs(low_near_ratio_v2->GetErrorY(i))), indentstr(3, "- label: sys").Data(), indentstr(3, "  symerror: ").Data(), low_near_ratio_syst->GetErrorY(i)));
        low_ratio_yval += yval;
        TString xval(Form("%s%e\n",indentstr(2, "- value: ").Data(), low_near_ratio->GetPointX(i)));
        low_ratio_xval += xval;
    }
    low_ratio_dep += low_ratio_yval;
    low_ratio_indep += low_ratio_xval;

    // away side
    low_ratio_dep += varheader(R"($1/{N}_{\mathrm{trig}} \mathrm{d}N_{\mathrm{assoc}}/\mathrm{d}{\Delta\varphi} away-side ratio \biggl(\frac{\mathrm{h--}\phi}{\mathrm{h--h}\biggr)$)")+qualifiers;

    low_ratio_yval = indentstr(1, "  values:\n");


    for(int i =  0; i<3; i++){
        TString yval(Form("%s%e\n%s\n%s\n%s%e\n%s\n%s%e}\n%s\n%s%e\n",indentstr(2, "- value: ").Data(), low_away_ratio->GetPointY(i), indentstr(2, "  errors:").Data(), indentstr(3, "- label: stat").Data(), indentstr(3, "  symerror: ").Data(), low_away_ratio->GetErrorY(i), indentstr(3, "- label: sys,v2").Data(), indentstr(3, "  asymerror: {plus: 0.0, minus: ").Data(), -2.0*(abs(low_away_ratio_v2->GetErrorY(i))), indentstr(3, "- label: sys").Data(), indentstr(3, "  symerror: ").Data(), low_away_ratio_syst->GetErrorY(i)));
        low_ratio_yval += yval;
    }
    low_ratio_dep += low_ratio_yval;

    // UE
    low_ratio_dep += varheader(R"($1/{N}_{\mathrm{trig}} \mathrm{d}N_{\mathrm{assoc}}/\mathrm{d}{\Delta\varphi} U.E. ratio \biggl(\frac{\mathrm{h--}\phi}{\mathrm{h--h}\biggr)$)")+qualifiers;

    low_ratio_yval = indentstr(1, "  values:\n");

    for(int i =  0; i<3; i++){
        TString yval(Form("%s%e\n%s\n%s\n%s%e\n%s\n%s%e\n",indentstr(2, "- value: ").Data(), low_ue_ratio->GetPointY(i), indentstr(2, "  errors:").Data(), indentstr(3, "- label: stat").Data(), indentstr(3, "  symerror: ").Data(), low_ue_ratio->GetErrorY(i), indentstr(3, "- label: sys").Data(), indentstr(3, "  symerror: ").Data(), low_ue_ratio_syst->GetErrorY(i)));
        low_ratio_yval += yval;
    }
    low_ratio_dep += low_ratio_yval;

    // total
    low_ratio_dep += varheader(R"($1/{N}_{\mathrm{trig}} \mathrm{d}N_{\mathrm{assoc}}/\mathrm{d}{\Delta\varphi} total ratio \biggl(\frac{\mathrm{h--}\phi}{\mathrm{h--h}\biggr)$)")+qualifiers;

    low_ratio_yval = indentstr(1, "  values:\n");

    for(int i =  0; i<3; i++){
        TString yval(Form("%s%e\n%s\n%s\n%s%e\n%s\n%s%e\n",indentstr(2, "- value: ").Data(), low_tot_ratio->GetPointY(i), indentstr(2, "  errors:").Data(), indentstr(3, "- label: stat").Data(), indentstr(3, "  symerror: ").Data(), low_tot_ratio->GetErrorY(i), indentstr(3, "- label: sys").Data(), indentstr(3, "  symerror: ").Data(), low_tot_ratio_syst->GetErrorY(i)));
        low_ratio_yval += yval;
    }
    low_ratio_dep += low_ratio_yval;

    TString low_ratio_output = low_ratio_indep + low_ratio_dep;
 
    // write low pt jet-yield table
    f = fopen("lowpt_ratio_data.yaml", "w+");    
    fprintf(f, "%s", low_ratio_output.Data());
    fclose(f);


}
