#include <initializer_list>
#include <cmath>
#include "../include/mva_analysis.h"

//#include "../include/mlp_analysis.h"

void produce_graphs(bool with_cut) {
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
  bool create_cards = false;
  std::string job_name = "1";
  std::string mva_cut = "";
  std::string method_name = "MLP";
  
  std::string sign = ">"; // direction of cut
  int min = 0; // the minimum value you want cuts to be from
  int max = 70; // max value you want cuts to be to
  double digits = 100; // number of digits + 1 of your cuts, e.g. if you put ur min as 40 then put 100 as digits to make it 0.4


  for (int i = 1; i < 2/*bg_chains.size()*/; i++)
  {
  		for (int j = 0; j < 1/*3*/; j++)
  		{
  			 /*MVAAnalysis::get_plots_varying_params(bg_chains, 6, signal_chain, data_chain, super_vars, "MLP", varying_params[0],
																																												NTrees, BoostType, AdaBoostBeta, SeparationType, nCuts, NeuronType, NCycles,
																																												HiddenLayers, LearningRate, unique_output_files, create_cards, "1");*/
		}
  }


MVAAnalysis::get_mva_results(bg_chains, 1, signal_chain, data_chain, super_vars, "test", method_name,
  NTrees[0],BoostType[0], AdaBoostBeta[0], SeparationType[0], nCuts[0], NeuronType[0], 
  NCycles[0], HiddenLayers[0], LearningRate[0],unique_output_files, create_cards, job_name, mva_cut, sign, min, max, digits);

}

int main(int argc, char** argv) {
  TApplication theApp("tapp", &argc, argv);
  produce_graphs(true);
  theApp.Run();
  return 0;
}
