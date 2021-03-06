#include "../include/super_vars.h"

SuperVars::SuperVars()
{
  forward_tag_eta        = new Variable("forward_tag_eta","Forward Tag #eta","-5.0","5.0","1.8",
                                        "5.0","60","5", "", true);
  dijet_deta             = new Variable("dijet_deta","Dijet #Delta#eta","3.5","8.0","4.2","8.0","25","10", "");
  metnomu_significance   = new Variable("metnomu_significance","MET Significance (No Muons)",
                                        "3.0","12.0","3.5","","50","10", "");
  sqrt_ht                = new Variable("sqrt_ht","Square Root HCAL Scalar Sum of Energy","0.0","35.0","9.0",
                                        "18.0","75","10", "GeV^{0.5}");
  alljetsmetnomu_mindphi = new Variable("alljetsmetnomu_mindphi","All Jets - MET Min. #Delta#phi (No Muons)",
                                        "0.0","3.15","2.0","","40","5", "");
  dijet_M                = new Variable("dijet_M","Dijet Mass","0.0","2000.0","800.0","","50","10", "GeV");
  metnomuons             = new Variable("metnomuons","MET (No Muons)","0.0","400.0","120.0","",
                                        "50","5", "GeV");
  jet1_pt                 = new Variable("jet1_pt","Jet1pt","0.0","","50.0","","30","", "GeV");
  jet2_pt                 = new Variable("jet2_pt","Jet2pt","0.0","","45.0","","30","", "GeV");
  jet1_eta = new Variable("jet1_eta","jet1_eta", "-5.0", "5.0", "-5.0", "5.0","50","1","");
  jet2_eta = new Variable("jet2_eta","jet2_eta", "-5.0", "5.0", "-5.0", "5.0","50","1","");
  jet1_phi = new Variable("jet1_phi","jet1_phi", "-4.0", "4.0", "-4.0", "4.0","50","1","");
  jet2_phi = new Variable("jet2_phi","jet2_phi", "-4.0", "4.0", "-4.0", "4.0","50","1","");
  jet_csv1 = new Variable("jet_csv1","jet_csv1", "-1.5", "1.5", "-1.5", "1.5","50","1","");
  jet_csv2 = new Variable("jet_csv2","jet_csv2", "-1.5", "1.5", "-1.5", "1.5","50","1","");
  dijet_dphi = new Variable("dijet_dphi","dijet_dphi", "0.0", "4.0", "0.0", "4.0","50","1","");
  metnomu_x = new Variable("metnomu_x","metnomu_x", "-400.0", "400.0", "-400.0", "400.0","50","1","");
  metnomu_y = new Variable("metnomu_y","metnomu_y", "-400.0", "400.0", "-400.0", "400.0","50","1","");
  sumet = new Variable("sumet","sumet", "0.0", "2400.0", "0.0", "2400.0","50","1","");
  mht = new Variable("mht","mht", "0.0", "3000.0", "0.0", "3000.0","50","1","");
  unclustered_et = new Variable("unclustered_et","unclustered_et", "0.0", "2000.0", "0.0", "2000.0","50","1","");
  jetmetnomu_mindphi = new Variable("jetmetnomu_mindphi","jetmetnomu_mindphi", "0.0", "3.5", "0.0", "3.5","50","1","");
  jetunclet_mindphi = new Variable("jetunclet_mindphi","jetunclet_mindphi", "0.0", "3.5", "0.0", "3.5","50","1","");
  metnomuunclet_dphi = new Variable("metnomuunclet_dphi","metnomuunclet_dphi", "0.0", "3.5", "0.0", "3.5","50","1","");
  dijetmetnomu_vectorialSum_pt = new Variable("dijetmetnomu_vectorialSum_pt","dijetmetnomu_vectorialSum_pt", "0.0", "400.0", "0.0", "400.0","50","1","");
  dijetmetnomu_ptfraction = new Variable("dijetmetnomu_ptfraction","dijetmetnomu_ptfraction", "0.0", "1.0", "0.0", "1.0","50","1","");
  jet1metnomu_scalarprod = new Variable("jet1metnomu_scalarprod","jet1metnomu_scalarprod", "-2000.0", "2000.0", "-2000.0", "2000.0","50","1","");
  jet2metnomu_scalarprod = new Variable("jet2metnomu_scalarprod","jet2metnomu_scalarprod", "-2000.0", "2000.0", "-2000.0", "2000.0","50","1","");
}

std::vector<Variable*> SuperVars::get_discriminating_vars()
{
  Variable* var_arr[] = {
                          alljetsmetnomu_mindphi, forward_tag_eta, dijet_deta, metnomu_significance,
  		  	                   sqrt_ht, dijet_M, metnomuons ,jet1_pt,jet2_pt, jet1_eta,jet2_eta, jet1_phi,jet2_phi,
																										jet_csv1,jet_csv2,dijet_dphi,metnomu_x,metnomu_y,sumet,mht,unclustered_et,
  		                      jetmetnomu_mindphi,jetunclet_mindphi,metnomuunclet_dphi,dijetmetnomu_vectorialSum_pt,
																										dijetmetnomu_ptfraction, jet1metnomu_scalarprod,jet2metnomu_scalarprod
  	                     };

  std::vector<Variable*> vars (var_arr, var_arr + sizeof(var_arr) / sizeof(var_arr[0]));

  return vars;
}

std::vector<Variable*> SuperVars::get_signal_cut_vars()
{
  Variable* var_arr[] =  {alljetsmetnomu_mindphi, metnomu_significance, dijet_deta, jet1_pt, jet2_pt};

  std::vector<Variable*> vars (var_arr, var_arr + sizeof(var_arr) / sizeof(var_arr[0]));

  return vars;
}


std::string SuperVars::get_final_cuts_str()
{
	 std::vector<Variable*> cut_vars = get_signal_cut_vars();
	 std::string cuts_str            = cut_vars[0]->build_multicut_selection(false, &cut_vars);

	 cuts_str.insert(cuts_str.find("(") + 1, HistoPlot::lep_sel_default());

	 return cuts_str;
}

std::string SuperVars::get_cuts_str_for_tmva()
{
  std::string cut_str = HistoPlot::replace_all(
  		                      HistoPlot::replace_all(
  		                     		 HistoPlot::replace_all(get_final_cuts_str(),
																																																			"(", " "),
																																																 "*total_weight_lepveto", ""),
																																														")", " ");

  return cut_str;
}
