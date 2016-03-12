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
  const char* mva_type = "BDT";  //go into train_and_run_BDT function to change input parameters

  /*MVAAnalysis::get_plots_varying_params(bg_chains, 0, signal_chain, data_chain, super_vars, "BDT", "NTrees", NTrees, BoostType,
  																			                   AdaBoostBeta, SeparationType, nCuts, NeuronType, NCycles, HiddenLayers);*/
  bool unique_output_files = false;
  bool create_cards = false;
  // boolean is for whether or not to create separate output app files
  MVAAnalysis::get_mva_results(bg_chains, 0, signal_chain, data_chain, super_vars, "test", "BDT", NTrees[0],
  																													BoostType[0], AdaBoostBeta[0], SeparationType[0], nCuts[0],
  																													NeuronType[0], NCycles[0], HiddenLayers[0], unique_output_files, create_cards, "1", "output>-0.25");

}

int main(int argc, char** argv) {
  TApplication theApp("tapp", &argc, argv);
  produce_graphs(true);
  theApp.Run();
  return 0;
}
