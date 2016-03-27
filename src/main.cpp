#include <initializer_list>
#include <cmath>
#include "../include/mva_analysis.h"

//#include "../include/mlp_analysis.h"

void produce_graphs(bool with_cut, const char* job_ptr) {
  SuperVars* super_vars             = new SuperVars();
  std::vector<Variable*> vars       = super_vars->get_discriminating_vars();
  std::vector<Variable*> cut_vars   = super_vars->get_signal_cut_vars();
  SuperChains* super_chains         = new SuperChains();
  std::vector<DataChain*> bg_chains = super_chains->get_bg_chains();
  DataChain* signal_chain           = super_chains->signal_chain;
  DataChain* data_chain             = super_chains->data_chain;
  std::vector<DataChain*> all_bg_chains = super_chains->get_all_bg_chains();
 
   const char* mva_type = "BDT";
  // the topmost folder for all root files so gitignore ignores properly
  std::string top_folder_name = "analysis";
  //const char* varying_params[] = {"NTrees", "AdaBoostBeta", "nCuts", "SeparationType"};
  const char* varying_params[] = {"HiddenLayers", "NCycles", "LearningRate","NeuronType"};
  // boolean is for whether or not to create separate output app files
  bool unique_output_files = false;
  // boolean is for whether or not to create datacards
  bool create_cards = true;
  std::string job_name = job_ptr;
  std::string mva_cut = "";
  std::string method_name = "BDT";
  int rel_bgs[] = {1, 2, 3, 6};
  std::string sign = ">"; // direction of cut
  int min = -90; // the minimum value you want cuts to be from
  int max = 90; // max value you want cuts to be to
  double digits = 100; // number of digits + 1 of your cuts, e.g. if you put ur min as 40 then put 100 as digits to make it 0.4

  const char* best_bdt_configs[21][5] = {
    // bg, ntrees, adaboostbeta, separation type, ncuts
    {"bg_wjets_ev", "300", "0.2", "MisClassificationError", "15"},
    {"bg_wjets_ev", "200", "0.2", "MisClassificationError", "15"},
    {"bg_wjets_ev", "200", "0.2", "MisClassificationError", "10"},
    {"bg_wjets_ev", "300", "0.2", "GiniIndex", "15"},
    {"bg_wjets_ev", "500", "0.2", "GiniIndex", "10"},
    /////////////////////////////////////////////////////////////
    {"bg_wjets_muv", "500", "0.2", "CrossEntropy", "15"},
    {"bg_wjets_muv", "500", "0.2", "GiniIndex", "15"},
    {"bg_wjets_muv", "300", "0.2", "GiniIndex", "5"},
    {"bg_wjets_muv", "200", "0.2", "CrossEntropy", "15"},
    {"bg_wjets_muv", "100", "0.2", "CrossEntropy", "5"},
    /////////////////////////////////////////////////////////////
    {"bg_wjets_tauv", "200", "0.2", "GiniIndex", "20"},
    {"bg_wjets_tauv", "300", "0.2", "GiniIndex", "15"},
    {"bg_wjets_tauv", "100", "0.5", "MisClassificationError", "5"},
    {"bg_wjets_tauv", "100", "0.5", "GiniIndex", "5"},
    {"bg_wjets_tauv", "200", "0.2", "GiniIndex", "15"},
    /////////////////////////////////////////////////////////////
    {"bg_zjets_vv", "300", "0.2", "CrossEntropy", "10"},
    {"bg_zjets_vv", "300", "0.2", "GiniIndex", "10"},
    {"bg_zjets_vv", "200", "0.2", "CrossEntropy", "10"},
    {"bg_zjets_vv", "80", "0.5", "GiniIndex", "5"},
    {"bg_zjets_vv", "200", "0.2", "GiniIndex", "5"},
    {"bg_zjets_vv", "50", "0.5", "GiniIndex", "10"},
  };

  for (int i = 18; i < 21; i++)
  {  
    int bg_to_train;
    const char* bg_label = best_bdt_configs[i][0];
    if (!strcmp(bg_label, "bg_wjets_ev"))
    {
      bg_to_train = 1;
    }
    else if (!strcmp(bg_label, "bg_wjets_muv"))
    {
      bg_to_train = 2;
    }
    else if (!strcmp(bg_label, "bg_wjets_tauv"))
    {
      bg_to_train = 3;
    }
    else
    {
      bg_to_train = 6;
    }

    std::string bg_label_str = bg_chains[bg_to_train]->label;
    std::string folder_name = "BDT/" + bg_label_str;
       
    MVAAnalysis::get_mva_results(bg_chains, bg_to_train, signal_chain, data_chain, super_vars, folder_name, method_name,
    				  best_bdt_configs[i][1],BoostType[0], best_bdt_configs[i][2], best_bdt_configs[i][3], best_bdt_configs[i][4],
              NeuronType[0], NCycles[0], HiddenLayers[0], LearningRate[0],unique_output_files, create_cards, job_name, 
              mva_cut, sign, min, max, digits);

  }


}

int main(int argc, char** argv) {
  TApplication theApp("tapp", &argc, argv);
  produce_graphs(true, argv[1]);
  theApp.Run();
  return 0;
}
