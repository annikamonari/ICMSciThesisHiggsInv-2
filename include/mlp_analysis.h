#ifndef Mlp_Analysis_h
#define Mlp_Analysis_h

#include "../include/histo_plot.h"
#include "../include/mc_weights.h"

#include <cstdlib>
#include "TObjString.h"
#include "TSystem.h"
#include "TROOT.h"
#include "TStopwatch.h"
#include "TMVA/Factory.h"
#include "TMVA/Tools.h"
#include "TMVA/Reader.h"
#include "TMVA/MethodCuts.h"

class MLPAnalysis
{
 public:
	 

	 static TFile* create_MLP(DataChain* bg_chain, DataChain* signal_chain, std::vector<Variable*>* variables,
																																									std::string folder_name, const char* NeuronType, const char* NCycles, const char* HiddenLayers,
																																									const char* LearningRate, std::string job_name);

	 static TTree* evaluate_MLP(DataChain* bg_chain, std::vector<Variable*>* variables, std::string output_name,
																																										std::string job_name, bool unique_output_files);

	 static DataChain* get_MLP_results(DataChain* bg_chain, std::vector<Variable*>* variables, std::string output_name,
																																																	std::string job_name, bool unique_output_files);

  static std::string MLP_options_str(const char* NeuronType, const char* NCycles, const char* HiddenLayers, const char* LearningRate);
	
	 static std::string MLP_output_name_str(const char* NeuronType, const char* NCycles,
																																									const char* HiddenLayers, const char* LearningRate,
																																									const char* bg_chain_label, std::string job_name);

	 static std::string MLP_output_file_path(std::string folder_name, std::string job_name, bool is_train_file,
																																										const char* NeuronType, const char* NCycles, const char* HiddenLayers,
																																										const char* LearningRate, const char* bg_chain_label);

};








#endif
