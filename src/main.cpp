#include <initializer_list>
#include <cmath>
#include "../include/mva_analysis.h"
#include <vector>


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
 
   //const char* mva_type = "BDT";
  // the topmost folder for all root files so gitignore ignores properly
  std::string top_folder_name = "analysis";
  //const char* varying_params[] = {"NTrees", "AdaBoostBeta", "nCuts", "SeparationType"};
  const char* varying_params[] = {"HiddenLayers", "NCycles", "LearningRate","NeuronType"};
  // boolean is for whether or not to create separate output app files
  bool unique_output_files = false;
  // boolean is for whether or not to create datacards
  bool create_cards = true;
  std::string job_name = "0";
  const char* jn = job_name.c_str();
  int counter = std::atoi(jn);
  std::string mva_cut = "";
  std::string method_name = "MLP";
  
  std::string sign = ">"; // direction of cut
  int min = 0; // the minimum value you want cuts to be from
  int max = 1; // max value you want cuts to be to
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
//TFile* trained_output;

//for(int i =0; i<7;i++){
  MVAAnalysis::get_mva_results(bg_chains, 6, signal_chain, data_chain, super_vars, "test", method_name,
  NTrees[0],BoostType[0], AdaBoostBeta[0], SeparationType[0], nCuts[0], NeuronType[0], 
  NCycles[1], HiddenLayers[counter], LearningRate[counter],unique_output_files, create_cards, job_name, mva_cut, sign, min, max, digits);

   /*const char* train_file_arr[1] = {trained_output->GetName()};
   std::vector<const char*> single_file_vector (train_file_arr,train_file_arr  + sizeof(train_file_arr)/sizeof(const char*));
   MVAAnalysis::get_estimators(single_file_vector);*/
//}


}

int main(int argc, char** argv) {
  TApplication theApp("tapp", &argc, argv);
  produce_graphs(true);
/*const char * files[] = {
"test/MLP-bg_zjets_vv-NeuronType=radial-NCycles=2000-HiddenLayers=2,4,8,16,32,64-LearningRate=0.000001-EstimatorType=CE-PCA.root",
"test/MLP-bg_zjets_vv-NeuronType=radial-NCycles=2000-HiddenLayers=2,4,8,16,32,64-LearningRate=0.00001-EstimatorType=CE-PCA.root",
"test/MLP-bg_zjets_vv-NeuronType=radial-NCycles=2000-HiddenLayers=2,4,8,16,32,64-LearningRate=0.0001-EstimatorType=CE-PCA.root",
"test/MLP-bg_zjets_vv-NeuronType=radial-NCycles=2000-HiddenLayers=2,4,8,16,32,64-LearningRate=0.001-EstimatorType=CE-PCA.root",
"test/MLP-bg_zjets_vv-NeuronType=radial-NCycles=2000-HiddenLayers=2,4,8,16,32,64-LearningRate=0.01-EstimatorType=CE-PCA.root",
"test/MLP-bg_zjets_vv-NeuronType=radial-NCycles=2000-HiddenLayers=2,4,8,16,32,64-LearningRate=0.1-EstimatorType=CE-PCA.root",
"test/MLP-bg_zjets_vv-NeuronType=radial-NCycles=2000-HiddenLayers=2,4,8,16,32,64-LearningRate=0.5-EstimatorType=CE-PCA.root",
"test/MLP-bg_zjets_vv-NeuronType=radial-NCycles=2000-HiddenLayers=2,4,8,16,32-LearningRate=0.000001-EstimatorType=CE-PCA.root",
"test/MLP-bg_zjets_vv-NeuronType=radial-NCycles=2000-HiddenLayers=2,4,8,16,32-LearningRate=0.00001-EstimatorType=CE-PCA.root",
"test/MLP-bg_zjets_vv-NeuronType=radial-NCycles=2000-HiddenLayers=2,4,8,16,32-LearningRate=0.0001-EstimatorType=CE-PCA.root",
"test/MLP-bg_zjets_vv-NeuronType=radial-NCycles=2000-HiddenLayers=2,4,8,16,32-LearningRate=0.001-EstimatorType=CE-PCA.root",
"test/MLP-bg_zjets_vv-NeuronType=radial-NCycles=2000-HiddenLayers=2,4,8,16,32-LearningRate=0.005-EstimatorType=CE-PCA.root",
"test/MLP-bg_zjets_vv-NeuronType=radial-NCycles=2000-HiddenLayers=2,4,8,16,32-LearningRate=0.01-EstimatorType=CE-PCA.root",
"test/MLP-bg_zjets_vv-NeuronType=radial-NCycles=2000-HiddenLayers=2,4,8,16,32-LearningRate=0.1-EstimatorType=CE-PCA.root",
"test/MLP-bg_zjets_vv-NeuronType=radial-NCycles=2000-HiddenLayers=2,4,8,16,32-LearningRate=0.5-EstimatorType=CE-PCA.root",
"test/MLP-bg_zjets_vv-NeuronType=radial-NCycles=2000-HiddenLayers=2,4,8,16-LearningRate=0.000001-EstimatorType=CE-PCA.root",
"test/MLP-bg_zjets_vv-NeuronType=radial-NCycles=2000-HiddenLayers=2,4,8,16-LearningRate=0.00001-EstimatorType=CE-PCA.root",
"test/MLP-bg_zjets_vv-NeuronType=radial-NCycles=2000-HiddenLayers=2,4,8,16-LearningRate=0.0001-EstimatorType=CE-PCA.root",
"test/MLP-bg_zjets_vv-NeuronType=radial-NCycles=2000-HiddenLayers=2,4,8,16-LearningRate=0.001-EstimatorType=CE-PCA.root",
"test/MLP-bg_zjets_vv-NeuronType=radial-NCycles=2000-HiddenLayers=2,4,8,16-LearningRate=0.01-EstimatorType=CE-PCA.root",
"test/MLP-bg_zjets_vv-NeuronType=radial-NCycles=2000-HiddenLayers=2,4,8,16-LearningRate=0.1-EstimatorType=CE-PCA.root",
"test/MLP-bg_zjets_vv-NeuronType=radial-NCycles=2000-HiddenLayers=2,4,8-LearningRate=0.000001-EstimatorType=CE-PCA.root",
"test/MLP-bg_zjets_vv-NeuronType=radial-NCycles=2000-HiddenLayers=2,4,8-LearningRate=0.00001-EstimatorType=CE-PCA.root",
"test/MLP-bg_zjets_vv-NeuronType=radial-NCycles=2000-HiddenLayers=2,4,8-LearningRate=0.0001-EstimatorType=CE-PCA.root",
"test/MLP-bg_zjets_vv-NeuronType=radial-NCycles=2000-HiddenLayers=2,4,8-LearningRate=0.001-EstimatorType=CE-PCA.root",
"test/MLP-bg_zjets_vv-NeuronType=radial-NCycles=2000-HiddenLayers=2,4,8-LearningRate=0.01-EstimatorType=CE-PCA.root",
"test/MLP-bg_zjets_vv-NeuronType=radial-NCycles=2000-HiddenLayers=2,4,8-LearningRate=0.1-EstimatorType=CE-PCA.root",
"test/MLP-bg_zjets_vv-NeuronType=radial-NCycles=2000-HiddenLayers=2,4,8-LearningRate=0.2-EstimatorType=CE-PCA.root",
"test/MLP-bg_zjets_vv-NeuronType=radial-NCycles=2000-HiddenLayers=2,4,8-LearningRate=0.5-EstimatorType=CE-PCA.root",
"test/MLP-bg_zjets_vv-NeuronType=radial-NCycles=2000-HiddenLayers=2,4-LearningRate=0.000001-EstimatorType=CE-PCA.root",
"test/MLP-bg_zjets_vv-NeuronType=radial-NCycles=2000-HiddenLayers=2,4-LearningRate=0.00001-EstimatorType=CE-PCA.root",
"test/MLP-bg_zjets_vv-NeuronType=radial-NCycles=2000-HiddenLayers=2,4-LearningRate=0.0001-EstimatorType=CE-PCA.root",
"test/MLP-bg_zjets_vv-NeuronType=radial-NCycles=2000-HiddenLayers=2,4-LearningRate=0.001-EstimatorType=CE-PCA.root",
"test/MLP-bg_zjets_vv-NeuronType=radial-NCycles=2000-HiddenLayers=2,4-LearningRate=0.002-EstimatorType=CE-PCA.root",
"test/MLP-bg_zjets_vv-NeuronType=radial-NCycles=2000-HiddenLayers=2,4-LearningRate=0.01-EstimatorType=CE-PCA.root",
"test/MLP-bg_zjets_vv-NeuronType=radial-NCycles=2000-HiddenLayers=2,4-LearningRate=0.1-EstimatorType=CE-PCA.root",
"test/MLP-bg_zjets_vv-NeuronType=radial-NCycles=2000-HiddenLayers=2,4-LearningRate=0.5-EstimatorType=CE-PCA.root",
"test/MLP-bg_zjets_vv-NeuronType=radial-NCycles=2000-HiddenLayers=2-LearningRate=0.000001-EstimatorType=CE-PCA.root",
"test/MLP-bg_zjets_vv-NeuronType=radial-NCycles=2000-HiddenLayers=2-LearningRate=0.00001-EstimatorType=CE-PCA.root",
"test/MLP-bg_zjets_vv-NeuronType=radial-NCycles=2000-HiddenLayers=2-LearningRate=0.0001-EstimatorType=CE-PCA.root",
"test/MLP-bg_zjets_vv-NeuronType=radial-NCycles=2000-HiddenLayers=2-LearningRate=0.001-EstimatorType=CE-PCA.root",
"test/MLP-bg_zjets_vv-NeuronType=radial-NCycles=2000-HiddenLayers=2-LearningRate=0.01-EstimatorType=CE-PCA.root",
"test/MLP-bg_zjets_vv-NeuronType=radial-NCycles=2000-HiddenLayers=2-LearningRate=0.1-EstimatorType=CE-PCA.root",
"test/MLP-bg_zjets_vv-NeuronType=radial-NCycles=2000-HiddenLayers=2-LearningRate=0.5-EstimatorType=CE-PCA.root",
};
std::vector<const char*> files_v (files, files + sizeof(files) / sizeof(const char*));
  MVAAnalysis::get_estimators(files_v);*/
  theApp.Run();
  return 0;
}

