#include "../include/mc_weights.h"

using namespace std;

std::string MCWeights::get_mc_selection_str(DataChain* bg_chain, Variable* variable, 
                                            std::vector<Variable*>* variables, std::string mva_cut)
{
  std::string selection_str = variable->build_multicut_selection(false, variables);

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
                              std::string selection)
{

  return HistoPlot::get_histo_integral(HistoPlot::build_1d_histo(data_chain, var, with_cut, false, "goff", variables, 
                                                                 selection), with_cut, var);
}

double MCWeights::get_all_bg_in_ctrl(std::vector<DataChain*> bg_chains, Variable* var, bool with_cut,
                                     std::vector<Variable*>* variables, std::string selection)
{
  double total_integral = 0.0;

  for (int i = 0; i < bg_chains.size(); i++)
  {
    total_integral += get_nevents(bg_chains[i], var, with_cut, variables, selection);
  }

  return total_integral;
}

double MCWeights::calc_mc_weight(DataChain* data, std::vector<DataChain*> bg_chains, DataChain* bg_chain,
                                 Variable* var, bool with_cut, std::vector<Variable*>* variables, std::string mva_cut)
{
//cout<<"bg to be calculated: "<<bg_chain->label<<"\n";
  std::string selection;
  if(variables != NULL){selection = get_mc_selection_str(bg_chain, var, variables, mva_cut);}
  else{ selection = "((classID==0)&&" + bg_chain->lep_sel + ")" + "*total_weight_lepveto";
  selection = HistoPlot::add_mva_cut_to_selection(selection, mva_cut);
   }
  //cout<<"mc weight selection: "<<selection<<"\n";
  double data_in_ctrl     = get_nevents(data, var, with_cut, variables, selection);
//cout<<"data in control region: "<<data_in_ctrl<<"\n";
  double ctrl_mc_in_ctrl  = get_nevents(bg_chain, var, with_cut, variables, selection);
//cout<<"desired background in control region: "<<ctrl_mc_in_ctrl<<"\n";

  double other_bg_in_ctrl = get_all_bg_in_ctrl(bg_chains, var, with_cut, variables, selection) - ctrl_mc_in_ctrl;
//cout<<"other backgrounds in control region: "<<other_bg_in_ctrl<<"\n";


  return (data_in_ctrl - other_bg_in_ctrl) / ctrl_mc_in_ctrl;
}

double MCWeights::calc_weight_error(DataChain* data, std::vector<DataChain*> bg_chains, DataChain* bg_chain,
                                 Variable* var, bool with_cut, std::vector<Variable*>* variables, std::string mva_cut)
{std::string selection;
  if(variables != NULL){selection = get_mc_selection_str(bg_chain, var, variables, mva_cut);}
  else{ selection = "((classID==0)&&" + bg_chain->lep_sel + ")" + "*total_weight_lepveto";}
  double data_N_C       = get_nevents(data, var, with_cut, variables, selection);
  double MC_N_C         = get_nevents(bg_chain, var, with_cut, variables, selection);
  double bg_N_C         = get_all_bg_in_ctrl(bg_chains, var, with_cut, variables, selection) - MC_N_C;
  double sigma_data_N_C = std::pow(data_N_C, 0.5);
  double sigma_MC_N_C   = std::pow(MC_N_C, 0.5);
  double sigma_bg_N_C   = std::pow(bg_N_C, 0.5);
  double err1           = sigma_data_N_C/MC_N_C;
  double err2           = sigma_bg_N_C/MC_N_C;
  double err3           = (data_N_C- bg_N_C)/(pow(MC_N_C,1.5));
  double error_sq       = std::pow(err1,2) + std::pow(err2,2) + std::pow(err3,2);
  double weight_error   = std::pow(error_sq, 0.5);

  return weight_error;

}

