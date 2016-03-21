#ifndef Mva_Analysis_h
#define Mva_Analysis_h

#include "roc_curves.h"
#include "classifier_outputs.h"

class MVAAnalysis
{
 public:
	 static TH1F* plot_output(DataChain* combined_data);

  static std::vector<double> get_categories(TH1F* output_histo);

  static std::vector<std::string> get_category_strs(std::vector<double> categories);

  static TH1F* build_histo(DataChain* combined_output, std::string selection_str, Variable* variable, std::string histo_label);

  static std::string build_output_sel_str(std::string category, std::string final_cuts);

  static TH1F* draw_signal(DataChain* combined_output, std::string category, std::string final_cuts, Variable* variable);

  static TH1F* draw_background(DataChain* combined_output, std::string category, std::string final_cuts, Variable* variable);

  static void style_histo(TH1F* histo);

  static void draw_histo(DataChain* combined_output, std::string final_cuts, Variable* variable);

  static std::vector<const char*> vary_parameters(std::vector<DataChain*> bg_chains, int bg_to_train, DataChain* signal_chain, DataChain* data_chain, SuperVars* super_vars,
																														std::string method_name, std::string dir_name, std::vector<const char*> NTrees, std::vector<const char*> BoostType,
																														std::vector<const char*> AdaBoostBeta, std::vector<const char*> SeparationType, std::vector<const char*> nCuts,
																														std::vector<const char*> NeuronType, std::vector<const char*> NCycles, std::vector<const char*> HiddenLayers,
																														std::vector<const char*> LearningRate, bool unique_output_files = false,
																														bool create_cards = false, std::string job_name = "1");

  static void get_plots_varying_params(std::vector<DataChain*> bg_chains, int bg_to_train, DataChain* signal_chain, DataChain* data_chain, SuperVars* super_vars,
																																							std::string method_name, std::string dir_name, std::vector<const char*> NTrees, std::vector<const char*> BoostType,
																																							std::vector<const char*> AdaBoostBeta, std::vector<const char*> SeparationType, std::vector<const char*> nCuts,
																																							std::vector<const char*> NeuronType, std::vector<const char*> NCycles, std::vector<const char*> HiddenLayers,
																																							std::vector<const char*> LearningRate,
																																							bool unique_output_files = false, bool create_cards = false, std::string job_name = "1");

  static std::vector<TFile*> get_files_from_paths(std::vector<const char*> file_paths);

  static std::string get_app_filename_for_chain(std::string app_output_name, const char* trained_bg_label, const char* app_label);

  static TFile* get_mva_results(std::vector<DataChain*> bg_chains, int bg_to_train, DataChain* signal_chain, DataChain* data_chain,
																														SuperVars* super_vars, std::string folder_name, std::string method_name, const char* NTrees = "800",
																														const char* BoostType = "AdaBoost", const char* AdaBoostBeta = "0.2", const char* SeparationType = "GiniIndex",
																														const char* nCuts = "30", const char* NeuronType = "sigmoid", const char* NCycles = "500",
																														const char* HiddenLayers = "5,5,5,5", const char* LearningRate = "0.02",
																														bool unique_output_files = true, bool create_cards = true,
																														std::string job_name = "1", std::string mva_cut = "", std::string sign = ">", int min = 10, int max = 20);

  static std::string create_auc_line_MLP(const char* bg_label, const char* NeuronType,
																																																						const char* NCycles, const char* HiddenLayers,
																																																						const char* LearningRate, double auc);

  static std::vector<std::string> get_mva_cut_range(std::string sign, int min, int max);

  static void create_datacards(DataChain* output_data_chain, DataChain* output_signal_chain, std::vector<DataChain*> output_bg_chains,
																																												Variable* mva_output, bool with_cut, std::vector<Variable*>* variables, TFile* trained_output,
																																												std::string method_name, std::string sign, int min, int max);

  static std::vector<DataChain*> get_output_bg_chains(std::vector<DataChain*> bg_chains, std::vector<Variable*> vars, std::string method_name,
																																																						std::string app_output_name, std::string job_name, const char* trained_bg_label,
																																																						bool unique_output_files);

  static DataChain* get_output_signal_chain(DataChain* signal_chain, std::vector<Variable*> vars, std::string method_name,
																																												std::string app_output_name, std::string job_name, const char* trained_bg_label,
																																												bool unique_output_files);

  static std::string build_output_graph_name(TFile* trained_output, std::string mva_cut = "");

};








#endif
