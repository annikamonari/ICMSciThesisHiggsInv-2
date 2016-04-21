#include "../include/super_vars.h"

SuperVars::SuperVars()
{
  forward_tag_eta        = new Variable("forward_tag_eta","Forward Tag #eta","1.8","4.5","1.8",
                                        "","20","1", "", true);//why is there upper limit
  dijet_deta             = new Variable("dijet_deta","Dijet #Delta#eta","4.2","","4.2","","20","1", "");//why is there upper limit
  metnomu_significance   = new Variable("metnomu_significance","MET Significance (No Muons)",
                                        "3.0","","3.5","","20","1", "");
  sqrt_ht                = new Variable("sqrt_ht","Square Root HCAL Scalar Sum of Energy","0.0","35.0","9.0",
                                        "18.0","20","1", "GeV^{0.5}");
  alljetsmetnomu_mindphi = new Variable("alljetsmetnomu_mindphi","All Jets - MET Min. #Delta#phi (No Muons)","2.0","","2.0","","20","1", "");
  dijet_M                = new Variable("dijet_M","Dijet Mass","0.0","2000","","","10","1", "GeV");// whyis there upper limit
  metnomuons             = new Variable("metnomuons","MET (No Muons)","0.0","2000","120.0","",//why the upper limit?
                                        "50","50", "GeV");
  jet1_pt                 = new Variable("jet1_pt","Jet1pt","50.0","","50.0","","20","", "GeV");
  jet2_pt                 = new Variable("jet2_pt","Jet2pt","45.0","","45.0","","20","", "GeV");
  jet1_eta = new Variable("jet1_eta","jet1_eta", "-5.0", "5.0", "-5.0", "5.0","50","1","");
  jet2_eta = new Variable("jet2_eta","jet2_eta", "-5.0", "5.0", "-5.0", "5.0","50","1","");
  jet1_phi = new Variable("jet1_phi","jet1_phi", "-4.0", "4.0", "-4.0", "4.0","50","1","");
  jet2_phi = new Variable("jet2_phi","jet2_phi", "-4.0", "4.0", "-4.0", "4.0","50","1","");
  jet_csv1 = new Variable("jet_csv1","jet_csv1", "-1.5", "1.5", "-1.5", "1.5","50","1","");
  jet_csv2 = new Variable("jet_csv2","jet_csv2", "-1.5", "1.5", "-1.5", "1.5","50","1","");
  dijet_dphi = new Variable("dijet_dphi","Dijet dPhi", "0.0", "4.0", "0.0", "4.0","50","1","");
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
  jet1metnomu_scalarprod = new Variable("jet1metnomu_scalarprod","jet1metnomu_scalarprod", "-500", "100", "-500", "100","20","1","");
  jet2metnomu_scalarprod = new Variable("jet2metnomu_scalarprod","jet2metnomu_scalarprod", "-300", "100", "-300", "100","20","1","");
  classID = new Variable("classID","signal_binary", "-1", "2", "-1", "2","2","1","");//0 for background 1 for signal
//CMS parked data variables
  parked_jet1_eta = new Variable("jet1_eta","jet1_eta", "-20", "5", "-4.7", "4.7","50","1","");
  parked_jet2_eta = new Variable("jet2_eta","jet2_eta", "-20", "5", "-4.7", "4.7","50","1","");

  parked_jet1_pt                 = new Variable("jet1_pt","Jet1pt","0.0","","50.0","","50","", "GeV");
  parked_jet2_pt                 = new Variable("jet2_pt","Jet2pt","0.0","","45.0","","50","", "GeV");
  parked_dijet_M                = new Variable("dijet_M","Dijet Mass","0.0","","1200.0","","50","1", "GeV");

  parked_dijet_deta             = new Variable("dijet_deta","Dijet #Delta#eta","3.6","","3.6","","50","1", "");
  parked_metnomuons             = new Variable("metnomuons","MET (No Muons)","0.0","","90","","50","1", "GeV");


  parked_alljetsmetnomu_mindphi = new Variable("alljetsmetnomu_mindphi","All Jets - MET Min. #Delta#phi (No Muons)",
                                        "2.3","3.15","2.3","","50","1", "");
  parked_metnomu_significance   = new Variable("metnomu_significance","MET Significance (No Muons)","0","","4","","50","1", "");


}


std::vector<Variable*> SuperVars::get_discriminating_vars()
{
  Variable* var_arr[] = {
                         dijetmetnomu_ptfraction, 
			dijetmetnomu_vectorialSum_pt, jet_csv2, 
                          dijet_dphi, dijet_M,
                          alljetsmetnomu_mindphi,metnomu_significance/*, forward_tag_eta, dijet_deta, 
  		  	                   sqrt_ht, metnomuons ,jet1_pt,jet2_pt, jet1_eta,jet2_eta, jet1_phi,jet2_phi,
																										jet_csv1, metnomu_x,metnomu_y,sumet,mht,unclustered_et,
  		                      jetmetnomu_mindphi,jetunclet_mindphi,metnomuunclet_dphi, jet1metnomu_scalarprod,jet2metnomu_scalarprod*/
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
std::vector<Variable*> SuperVars::get_parked_vars()
{
  Variable* var_arr[] =  {parked_jet1_pt, parked_jet2_pt, parked_dijet_deta, parked_metnomuons, parked_jet1_eta, parked_jet2_eta,
      parked_dijet_M, parked_alljetsmetnomu_mindphi, parked_metnomu_significance};

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
  HistoPlot::replace_all( HistoPlot::replace_all(get_final_cuts_str(),
       "(", " "), "*total_weight_lepveto", ""),")", " ");

  return cut_str;
}
