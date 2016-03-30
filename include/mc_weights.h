#ifndef Mc_Weights_h
#define Mc_Weights_h

#include "histo_plot.h"

class MCWeights
{
 public:
  static std::string get_mc_selection_str(DataChain* bg_chain, Variable* variable, 
  	                                      std::vector<Variable*>* variables, std::string mva_cut = "",bool if_parked =false);

  static double get_nevents(DataChain* chain_of_data, Variable* var, bool with_cut,
                            std::vector<Variable*>* variables, std::string selection,
																												const char* trained_bg_label = "bg_zjets_vv", bool double_test_bg = false);

  static double get_all_bg_in_ctrl(std::vector<DataChain*> bg_chains, Variable* var, bool with_cut,
                                   std::vector<Variable*>* variables, std::string selection, int trained_bg,
																																			bool double_test_bg = false);

  static double calc_mc_weight(DataChain* data, std::vector<DataChain*> bg_chains, DataChain* bg_chain,
                               Variable* var, bool with_cut, std::vector<Variable*>* variables, std::string mva_cut = "",
int trained_bg = 6, bool double_test_bg = false, bool if_parked =false);

  static double calc_weight_error(DataChain* data, std::vector<DataChain*> bg_chains, DataChain* bg_chain,
                                 Variable* var, bool with_cut, std::vector<Variable*>* variables, int trained_bg,
				bool double_test_bg, std::string mva_cut = "", bool if_parked =false);
};



#endif
