#include <initializer_list>
#include <cmath>
#include "../include/mva_analysis.h"

//#include "../include/mlp_analysis.h"
using namespace std;
void produce_graphs(bool with_cut, const char* job_ptr, const char* bg_to_train, const char* vary_param) {
  SuperVars* super_vars             = new SuperVars();
  std::vector<Variable*> vars       = super_vars->get_discriminating_vars();
  std::vector<Variable*> cut_vars   = super_vars->get_signal_cut_vars();
  SuperChains* super_chains         = new SuperChains();
  std::vector<DataChain*> bg_chains = super_chains->get_bg_chains();
  DataChain* signal_chain           = super_chains->signal_chain;
  DataChain* data_chain             = super_chains->data_chain;
  const char* mva_type = "BDT";
  // the topmost folder for all root files so gitignore ignores properly
  std::string top_folder_name = "analysis";
  const char* varying_params[] = {"NTrees", "AdaBoostBeta", "nCuts", "SeparationType"};
  // boolean is for whether or not to create separate output app files
  bool unique_output_files = true;
  // boolean is for whether or not to create datacards
  bool create_cards = false;
  std::string job_number = job_ptr;
  int bg_int = atoi(bg_to_train);
  int param_int = atoi(vary_param);
//             0         1             2                3           4          5           6            7
//bg[8] = {"bg_zll","bg_wjets_ev","bg_wjets_muv","bg_wjets_tauv", "bg_top", "bg_VV", "bg_zjets_vv", "bg_qcd"};

      

  			 MVAAnalysis::get_plots_varying_params(bg_chains, bg_int, signal_chain, data_chain, super_vars, "BDT", varying_params[param_int],
																																												NTrees, BoostType, AdaBoostBeta, SeparationType, nCuts, NeuronType, NCycles,
																																												HiddenLayers, unique_output_files, create_cards, job_number);

  /*MVAAnalysis::get_mva_results(bg_chains, 3, signal_chain, data_chain, super_vars, "test", "BDT", NTrees[0],
  																													BoostType[0], AdaBoostBeta[1], SeparationType[2], nCuts[3],
  																													NeuronType[0], NCycles[0], HiddenLayers[0], unique_output_files, create_cards, "1", "");*/

}

int main(int argc, char** argv) {
  TApplication theApp("tapp", &argc, argv);
  produce_graphs(true, argv[1], argv[2], argv[3]);
  theApp.Run();
  return 0;
}

/*
      std::string bg_label = bg_chains[i]->label;
cout<<"bg label :"<<bg_label<<"\n";
						std::string param = varying_params[j];
cout<<"param: "<<param<<"\n";
      std::string param_vary_path = "analysis/BDT_varying_";
      param_vary_path.append(param);
cout<<"remove folder "<<param_vary_path<<"\n";

      param_vary_path.append("/");
cout<<"remove folder "<<param_vary_path<<"\n";

      param_vary_path.append(bg_label);
cout<<"remove folder "<<param_vary_path<<"\n";

      

      param_vary_path.append("/BDT-job_name=2-bg_top-NTrees=50-BoostType=AdaBoost-AdaBoostBeta=0.1-SeparationType=GiniIndex-nCuts=5.root");
for(int trees=0; tress< NTrees.size();trees++)
{
string file = BDTAnalysis::BDT_output_name_str( NTrees[trees], BoostType[0], AdaBoostBeta[0],
																																													 SeparationType[0], nCuts, bg_chains[i]->label,
																																													std::string job_name);
}
cout<<"remove folder "<<param_vary_path<<"\n";
*/

