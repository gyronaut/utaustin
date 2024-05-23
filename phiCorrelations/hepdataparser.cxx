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



TString jet_data_to_yaml(TGraphErrors* jetgraph, TGraphErrors* jetv2graph, bool isXaxis, bool isHPhi){
    TString values_string = indentstr(1, "  values:\n");
    int n = jetgraph->GetN();
    float scale = 1.0;
    if(isHPhi) scale = 250.0;
    if(isXaxis){
        for(int i =  0; i<n; i++){
            TString xval(Form("%s%e\n",indentstr(2, "- value: ").Data(), jetgraph->GetPointX(i)));
            values_string += xval;
        }

    }else{
        for(int i =  0; i<n; i++){
            TString yval(Form("%s%e\n%s\n%s\n%s%e\n%s\n%s%e}\n",indentstr(2, "- value: ").Data(), jetgraph->GetPointY(i)/scale, indentstr(2, "  errors:").Data(), indentstr(3, "- label: stat").Data(), indentstr(3, "  symerror: ").Data(), jetgraph->GetErrorY(i)/scale, indentstr(3, "- label: sys,v2").Data(), indentstr(3, "  asymerror: {plus: 0.0, minus: ").Data(), -(jetgraph->GetPointY(i)-jetv2graph->GetPointY(i))/scale));
            values_string += yval;
        }
    }
    return values_string;
}

TString ratio_data_to_yaml(TGraphErrors* ratiograph, TGraphErrors* ratiov2graph, TGraphErrors* ratiosystgraph, bool isXaxis, bool isJet = true){

    TString values_string = indentstr(1, "  values:\n");
    int n = ratiograph->GetN();
    if(isXaxis){
        for(int i =  0; i<n; i++){
            TString xval(Form("%s%e\n",indentstr(2, "- value: ").Data(), ratiograph->GetPointX(i)));
            values_string += xval;
        }
    }else{
        if(isJet){
            for(int i =  0; i<n; i++){
                TString yval(Form("%s%e\n%s\n%s\n%s%e\n%s\n%s%e}\n%s\n%s%e\n",indentstr(2, "- value: ").Data(), ratiograph->GetPointY(i), indentstr(2, "  errors:").Data(), indentstr(3, "- label: stat").Data(), indentstr(3, "  symerror: ").Data(), ratiograph->GetErrorY(i), indentstr(3, "- label: sys,v2").Data(), indentstr(3, "  asymerror: {plus: 0.0, minus: ").Data(), -2.0*(abs(ratiov2graph->GetErrorY(i))), indentstr(3, "- label: sys").Data(), indentstr(3, "  symerror: ").Data(), ratiosystgraph->GetErrorY(i)));
                values_string += yval;
            }
        }else{
            for(int i =  0; i<n; i++){
                TString yval(Form("%s%e\n%s\n%s\n%s%e\n%s\n%s%e\n",indentstr(2, "- value: ").Data(), ratiograph->GetPointY(i), indentstr(2, "  errors:").Data(), indentstr(3, "- label: stat").Data(), indentstr(3, "  symerror: ").Data(), ratiograph->GetErrorY(i), indentstr(3, "- label: sys").Data(), indentstr(3, "  symerror: ").Data(), ratiosystgraph->GetErrorY(i)));
                values_string += yval;
            }
        }

    }
    return values_string;
}


void hepdataparser(bool isLowPT){

    TString qualifiers;
    TFile* datafile;
    TString file_prefix;
    if(isLowPT){
        datafile = new TFile("./plotcode/ratiocode/fitsyst_cr2lowtest.root");
        qualifiers = Form(R"X(      qualifiers:
        - name: RE
          value: P PB --> H PHI X
        - name: SQRT(S)
          units: GEV
          value: 5020.0
        - name: ETARAP
          value: -0.8 to 0.8
        - name: PT(trig)
          units: GEV/C
          value: 4.0 to 8.0
        - name: PT(assoc)
          units: GEV/C
          value: 1.5 to 2.5
)X");
        file_prefix = "lowpt";
    }else{
        datafile = new TFile("./plotcode/ratiocode/fitsyst_cr2hightest.root");
        qualifiers = Form(R"X(      qualifiers:
        - name: RE
          value: P PB --> H PHI X
        - name: SQRT(S)
          units: GEV
          value: 5020.0
        - name: ETARAP
          value: -0.8 to 0.8
        - name: PT(trig)
          units: GEV/C
          value: 4.0 to 8.0
        - name: PT(assoc)
          units: GEV/C
          value: 2.5 to 4.0
)X");
        file_prefix = "highpt";
    }
      

    TString depvar("dependent_variables:\n");
    TString indepvar("independent_variables:\n");

//  get jet-yield graphs   
    TGraphErrors* near_phi = (TGraphErrors*)datafile->Get("yieldsNear");
    TGraphErrors* near_phi_v2 = (TGraphErrors*)datafile->Get("yieldsNearv2");
    TGraphErrors* away_phi = (TGraphErrors*)datafile->Get("yieldsAway");
    TGraphErrors* away_phi_v2 = (TGraphErrors*)datafile->Get("yieldsAwayv2");
    TGraphErrors* near_h = (TGraphErrors*)datafile->Get("yieldshhNear");
    TGraphErrors* near_h_v2 = (TGraphErrors*)datafile->Get("yieldshhNearv2");
    TGraphErrors* away_h = (TGraphErrors*)datafile->Get("yieldshhAway");
    TGraphErrors* away_h_v2 = (TGraphErrors*)datafile->Get("yieldshhAwayv2");
//  get ratio graphs
    TGraphErrors* near_ratio = (TGraphErrors*)datafile->Get("nchratiosNear");
    TGraphErrors* near_ratio_syst = (TGraphErrors*)datafile->Get("nchratiosNearSyst");
    TGraphErrors* near_ratio_v2 = (TGraphErrors*)datafile->Get("nearv2fly");
    
    TGraphErrors* away_ratio = (TGraphErrors*)datafile->Get("nchratiosAway");
    TGraphErrors* away_ratio_syst = (TGraphErrors*)datafile->Get("nchratiosAwaySyst");
    TGraphErrors* away_ratio_v2 = (TGraphErrors*)datafile->Get("awayv2fly");
    
    TGraphErrors* ue_ratio = (TGraphErrors*)datafile->Get("nchratiosBulk");
    TGraphErrors* ue_ratio_syst = (TGraphErrors*)datafile->Get("nchratiosBulkSyst");

    TGraphErrors* tot_ratio = (TGraphErrors*)datafile->Get("nchratiosTot");
    TGraphErrors* tot_ratio_syst = (TGraphErrors*)datafile->Get("nchratiosTotSyst");

    
    TString yields_indep, yields_dep, yields_xval, yields_yval;
    yields_indep += indepvar + varheader(R"($\langle N_{\mathrm{ch}} \rangle_{|\eta|<0.5}$)");
    yields_indep += jet_data_to_yaml(near_phi, near_phi_v2, true, true);

// h-phi jet yields
    yields_dep += depvar + varheader(R"(Per-trigger h--$\phi$ yields in near-side jet)")+qualifiers;
    yields_dep += jet_data_to_yaml(near_phi, near_phi_v2, false, true);

    yields_dep += varheader(R"(Per-trigger h--$\phi$ yields in away-side jet)")+qualifiers;
    yields_dep += jet_data_to_yaml(away_phi, away_phi_v2, false, true);

//  h-h jet yields
    yields_dep += varheader(R"(Per-trigger h--h yields in near-side jet)")+qualifiers;
    yields_dep += jet_data_to_yaml(near_h, near_h_v2, false, false);

    yields_dep += varheader(R"(Per-trigger h--h yields in away-side jet)")+qualifiers;
    yields_dep += jet_data_to_yaml(away_h, away_h_v2, false, false);

    TString jet_output = yields_indep + yields_dep;
 
    // write jet-yield data
    FILE *f = fopen(Form("%s_jet_data.yaml", file_prefix.Data()), "w+");    
    fprintf(f, "%s", jet_output.Data());
    fclose(f);

    
// ratio data
    TString ratio_indep, ratio_dep, ratio_yval, ratio_xval;
    ratio_indep += indepvar + varheader(R"($\langle N_{\mathrm{ch}} \rangle_{|\eta|<0.5}$)");
    ratio_indep += ratio_data_to_yaml(near_ratio, near_ratio_v2, near_ratio_syst, true);
    // near side
    
    ratio_dep += depvar + varheader(R"($1/{N}_{\mathrm{trig}} \mathrm{d}N_{\mathrm{assoc}}/\mathrm{d}{\Delta\varphi} near-side ratio \biggl(\frac{\mathrm{h--}\phi}{\mathrm{h--h}\biggr)$)")+qualifiers;
    ratio_dep += ratio_data_to_yaml(near_ratio, near_ratio_v2, near_ratio_syst, false);
  
    // away side
    ratio_dep += varheader(R"($1/{N}_{\mathrm{trig}} \mathrm{d}N_{\mathrm{assoc}}/\mathrm{d}{\Delta\varphi} away-side ratio \biggl(\frac{\mathrm{h--}\phi}{\mathrm{h--h}\biggr)$)")+qualifiers;
    ratio_dep += ratio_data_to_yaml(away_ratio, away_ratio_v2, away_ratio_syst, false);
    
    // UE
    ratio_dep += varheader(R"($1/{N}_{\mathrm{trig}} \mathrm{d}N_{\mathrm{assoc}}/\mathrm{d}{\Delta\varphi} U.E. ratio \biggl(\frac{\mathrm{h--}\phi}{\mathrm{h--h}\biggr)$)")+qualifiers;
    ratio_dep += ratio_data_to_yaml(ue_ratio, NULL, ue_ratio_syst, false, false);
    
    // total
    ratio_dep += varheader(R"($1/{N}_{\mathrm{trig}} \mathrm{d}N_{\mathrm{assoc}}/\mathrm{d}{\Delta\varphi} total ratio \biggl(\frac{\mathrm{h--}\phi}{\mathrm{h--h}\biggr)$)")+qualifiers;
    ratio_dep += ratio_data_to_yaml(tot_ratio, NULL, tot_ratio_syst, false, false);
    
    TString ratio_output = ratio_indep + ratio_dep;
 
    // write ratio data
    f = fopen(Form("%s_ratio_data.yaml", file_prefix.Data()), "w+");    
    fprintf(f, "%s", ratio_output.Data());
    fclose(f);


}
