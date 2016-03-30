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
  bool create_cards = false;
  std::string job_name = job_ptr;
  std::string mva_cut = "";
  std::string method_name = "MLP";
  int rel_bgs[] = {1, 2, 3, 6};
  std::string sign = ">"; // direction of cut
  int min = 20; // the minimum value you want cuts to be from
  int max = 45; // max value you want cuts to be to
  double digits = 100; // number of digits + 1 of your cuts, e.g. if you put ur min as 40 then put 100 as digits to make it 0.4
/*
  for (int i = 0; i < 1i < LearningRate.size(); i++)
  {  
    
    for (int k = 0; k < 1; k++)
    {
      for (int j = 0; j < HiddenLayers.size(); j++)
      {
          int bg_to_train = rel_bgs[k];
          std::string bg_label = bg_chains[bg_to_train]->label;
          std::string folder_name = "MLP/" + bg_label;
*/        

    const char* file_arr[] ={

      "test/MLP-bg_zjets_vv-NeuronType=radial-NCycles=800-HiddenLayers=10,10-LearningRate=0.02-EstimatorType=CE.root",
      "MLP/bg_zjets_vv/MLP-bg_zjets_vv-NeuronType=radial-NCycles=800-HiddenLayers=20,40-LearningRate=0.002-EstimatorType=CE.root"
      //"test/MLP-bg_wjets_muv-NeuronType=radial-NCycles=800-HiddenLayers=5,10-LearningRate=0.002-EstimatorType=CE.root"
    };
    std::vector<const char*> file_paths (file_arr, file_arr +sizeof(file_arr)/sizeof(const char*));
    std::vector<TFile*> files = MVAAnalysis::get_files_from_paths(file_paths);
    RocCurves::get_rocs(files, signal_chain, bg_chains[6], super_vars, "MLP", "test");
    ClassifierOutputs::plot_classifiers_for_all_files(files, "MLP", "test", bg_chains[6]->label);
          /*
    				MVAAnalysis::get_mva_results(bg_chains, 2, signal_chain, data_chain, super_vars, "test", method_name,
    				  NTrees[0],BoostType[0], AdaBoostBeta[0], SeparationType[0], nCuts[0], NeuronType[0],
    				  NCycles[0], HiddenLayers[0], LearningRate[0],unique_output_files, create_cards, job_name, 
              mva_cut, sign, min, max, digits); 
*/
          
/*
    		}
    }
    
  }
*/

}

int main(int argc, char** argv) {
  TApplication theApp("tapp", &argc, argv);
  produce_graphs(true, argv[1]);
  theApp.Run();
  return 0;
}
