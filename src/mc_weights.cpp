#include "../include/mc_weights.h"

using namespace std;

std::string MCWeights::get_mc_selection_str(DataChain* bg_chain, Variable* variable, 
                                            std::vector<Variable*>* variables, std::string mva_cut, bool if_parked)
{
  std::string selection_str;
  if(!if_parked){
    selection_str = variable->build_multicut_selection(false, variables);
  }
  else if(if_parked){
    selection_str = variable->build_parked_selection(variables);
  }
  if (bg_chain->lep_sel != "")
  {
  	 selection_str.insert(selection_str.find("(") + 1, bg_chain->lep_sel);
  }

  if (selection_str.find("(&&", 0) == 0)
  {
  		selection_str.erase(1, 2);
  }
  // This function below checks to see if mva_cut != ""
  selection_str = HistoPlot::add_mva_cut_to_selection(selection_str, mva_cut);
  //std::cout << "Selection str in get_mc_select_str: " << selection_str << std::endl;
  return selection_str;
}

double MCWeights::get_nevents(DataChain* data_chain, Variable* var, bool with_cut, std::vector<Variable*>* variables, 
                              std::string selection, const char* trained_bg_label, bool double_test_bg)
{

TH1F* h = HistoPlot::build_1d_histo(data_chain, var, true, false, "goff", variables, selection);
  double integral = h->Integral();//HistoPlot::get_histo_integral(, with_cut, var);
  //cout<<"integral: "<<integral;

  return integral;
}

double MCWeights::get_all_bg_in_ctrl(std::vector<DataChain*> bg_chains, Variable* var, bool with_cut,
     std::vector<Variable*>* variables, std::string selection, int trained_bg, bool double_test_bg)
{
  double total_integral = 0.0;

  for (int i = 0; i < bg_chains.size(); i++)
  {
    double integral = get_nevents(bg_chains[i], var, with_cut, variables, selection, bg_chains[trained_bg]->label);
      if(double_test_bg){if(i==trained_bg){integral=integral*2;}}
    total_integral += integral;
  }

  return total_integral;
}

double MCWeights::calc_mc_weight(DataChain* data, std::vector<DataChain*> bg_chains, DataChain* bg_chain,Variable* var, 
bool with_cut, std::vector<Variable*>* variables, std::string mva_cut,int trained_bg, bool double_test_bg, bool if_parked)
{
//cout<<"bg to be calculated: !!!!!!!!!!!"<<bg_chain->label<<"!!!!!!!!!!!\n";
  std::string selection;
  if(variables != NULL){selection = get_mc_selection_str(bg_chain, var, variables, mva_cut, if_parked);}
  else{ 
    if(!strcmp(bg_chain->label, "bg_wjets_tauv"))
    {
	selection = "((alljetsmetnomu_mindphi>2.0)&&(classID==0)&&" + bg_chain->lep_sel + ")" + "*total_weight_lepveto";
    }
    	else{selection = "((alljetsmetnomu_mindphi>2.0)&&(classID==0)&&"+bg_chain->lep_sel+")*total_weight_lepveto";//}
    	selection = HistoPlot::add_mva_cut_to_selection(selection, mva_cut);}
   }
  double weight;
  //cout<<"mc weight selection: "<<selection<<"\n";
  double data_in_ctrl     = get_nevents(data, var, with_cut, variables, selection, bg_chains[trained_bg]->label, double_test_bg);
//cout<<"data in control region: "<<data_in_ctrl<<"\n";
  double ctrl_mc_in_ctrl  = get_nevents(bg_chain, var, with_cut, variables, selection, bg_chains[trained_bg]->label, double_test_bg);
//cout<<"desired background in control region: "<<ctrl_mc_in_ctrl<<"\n";
if(ctrl_mc_in_ctrl!=0){ 
  double other_bg_in_ctrl = get_all_bg_in_ctrl(bg_chains, var, with_cut, variables, selection, trained_bg, double_test_bg) - ctrl_mc_in_ctrl;
//cout<<"other backgrounds in control region: "<<other_bg_in_ctrl<<"\n";

 weight = (data_in_ctrl - other_bg_in_ctrl) / ctrl_mc_in_ctrl;
 }
  else{weight=1;}
  return weight;
}

double MCWeights::calc_weight_error(DataChain* data, std::vector<DataChain*> bg_chains, DataChain* bg_chain,
Variable* var, bool with_cut, std::vector<Variable*>* variables, int trained_bg,bool double_test_bg, std::string mva_cut, bool if_parked)
{
  std::string selection;
  if(variables != NULL){selection = get_mc_selection_str(bg_chain, var, variables, mva_cut);}
  else
  { 
    if(!strcmp(bg_chain->label, "bg_wjets_tauv"))
    {
	selection = "((alljetsmetnomu_mindphi>2.0)&&(classID==0)&&" + bg_chain->lep_sel + ")" + "*total_weight_lepveto";
    }
    else
    {selection = "((alljetsmetnomu_mindphi>2.0)&&(classID==0)&&" + bg_chain->lep_sel + ")" + "*total_weight_lepveto";//}
    selection = HistoPlot::add_mva_cut_to_selection(selection, mva_cut);}
  }
  double weight_error;
  double data_N_C       = get_nevents(data, var, with_cut, variables, selection, bg_chains[trained_bg]->label, double_test_bg);
  double MC_N_C         = get_nevents(bg_chain, var, with_cut, variables, selection, bg_chains[trained_bg]->label, double_test_bg);
  double bg_N_C         = get_all_bg_in_ctrl(bg_chains, var, with_cut, variables, selection, trained_bg, double_test_bg) - MC_N_C;
  double sigma_data_N_C = std::pow(data_N_C, 0.5);
  double sigma_MC_N_C   = std::pow(MC_N_C, 0.5);
if(MC_N_C!=0){
  double sigma_bg_N_C   = std::pow(bg_N_C, 0.5);
  double err1           = sigma_data_N_C/MC_N_C;
  double err2           = sigma_bg_N_C/MC_N_C;
  double err3           = (data_N_C- bg_N_C)/(pow(MC_N_C,1.5));
  double error_sq       = std::pow(err1,2) + std::pow(err2,2) + std::pow(err3,2);
  weight_error   = std::pow(error_sq, 0.5);
}
else if(MC_N_C==0){cout<<"Warning!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n"<<bg_chains[trained_bg]->label<<"target background"<<" = "<<MC_N_C;
weight_error = 1000;}
//cout<<bg_chain->label<<"=============================weight error: "<<weight_error<<"\n";
  return weight_error;

}
