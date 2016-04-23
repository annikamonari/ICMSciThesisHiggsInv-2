#include "../include/mc_weights.h"

using namespace std;

std::string MCWeights::get_mc_selection_str(DataChain* bg_chain, Variable* variable, 
                                            std::vector<Variable*>* variables, std::string mva_cut, bool if_parked)
{
  std::string selection_str;
  if(!if_parked){
    selection_str = variable->build_multicut_selection(false, variables);
  }
  else if(if_parked){
	if(!strcmp(bg_chain->label, "bg_wjets_tauv"))
	{
		selection_str = "(&&((jet1_eta<0)&&(jet2_eta>0)||(jet1_eta>0)&&(jet2_eta<0))&&((metnomuons/metnomu_significance)>4)&&(metnomu_significance>4)&&(jetmetnomu_mindphi>1.0)&&(mht>20)&&(dijet_M>1200.0)&&(jet2_eta>-4.7)&&(jet2_eta<4.7)&&(jet1_eta>-4.7)&&(jet1_eta<4.7)&&(metnomuons>90)&&(dijet_deta>3.6)&&(jet2_pt>45.0)&&(jet1_pt>50.0))*total_weight_lepveto";
	}
	else
	{
    		selection_str = "(&&((jet1_eta<0)&&(jet2_eta>0)||(jet1_eta>0)&&(jet2_eta<0))&&((metnomuons/metnomu_significance)>4)&&(metnomu_significance>4)&&(alljetsmetnomu_mindphi>2.3)&&(dijet_M>1200.0)&&(jet2_eta>-4.7)&&(jet2_eta<4.7)&&(jet1_eta>-4.7)&&(jet1_eta<4.7)&&(metnomuons>90)&&(dijet_deta>3.6)&&(jet2_pt>45.0)&&(jet1_pt>50.0))*total_weight_lepveto";//variable->build_parked_selection(variables);
	}
  }
  if (bg_chain->lep_sel != "")
  {
  	 selection_str.insert(selection_str.find("(") + 1, bg_chain->lep_sel);
  }

  if (selection_str.find("(&&", 0) == 0)
  {
  		selection_str.erase(1, 2);
  }

  //cout<<"mc selection string="<<selection_str<<"\n";
  // This function below checks to see if mva_cut != ""
  selection_str = HistoPlot::add_mva_cut_to_selection(selection_str, mva_cut);
  //std::cout << "Selection str in get_mc_select_str: " << selection_str << std::endl;
  return selection_str;
}

double MCWeights::get_nevents(DataChain* data_chain, Variable* var, bool with_cut, std::vector<Variable*>* variables, 
                              std::string selection, const char* trained_bg_label, bool double_test_bg)
{

  TH1F* h = HistoPlot::build_1d_histo(data_chain, var, true, false, "goff", variables, selection);
  double integral = h->Integral();
  //cout<<"integral: "<<integral;

  return integral;
}

double MCWeights::get_other_bg_in_ctrl(int desired_bg_index,std::vector<double> mc_weights,std::vector<DataChain*> bg_chains, Variable* var, bool with_cut,
     std::vector<Variable*>* variables, std::string selection, int trained_bg, bool double_test_bg)
{
  double total_integral = 0.0;
  for (int i = 0; i < bg_chains.size(); i++)
  {
	if(i!=desired_bg_index)
	{
	  	std::string bg_selection = HistoPlot::add_mc_to_selection(bg_chains[i], var, selection, mc_weights[i]);
    		TH1F* bg_h = HistoPlot::build_1d_histo(bg_chains[i], var, true, false, "goff", variables, bg_selection,mc_weights[i]);
    		double integral = bg_h->Integral();     
		if(double_test_bg){if(i==trained_bg){integral=integral*2;}}
    		total_integral += integral;
		//if(i==7){cout<<"QCD bg in control ="<<integral<<"\n";}
	}
  }

  return total_integral;
}
double MCWeights::calc_nunu_mc_weight(DataChain* data, std::vector<DataChain*> bg_chains, DataChain* nunu_chain,Variable* var, 
bool with_cut, std::vector<Variable*>* variables, std::string mva_cut,int trained_bg, bool double_test_bg,bool if_parked)
{
  double mc_arr[8] = {1,1,1,1,1,1,1,1};
  std::vector<double> mc_weights (mc_arr, mc_arr + sizeof(mc_arr)/sizeof(mc_arr[0]));

  string selection;
  if(variables != NULL){selection = get_mc_selection_str(bg_chains[0], var, variables, mva_cut,if_parked);}
  else{  
    selection = "((alljetsmetnomu_mindphi>2.0)&&(classID==0)&&"+bg_chains[0]->lep_sel+")*total_weight_lepveto";
    	selection = HistoPlot::add_mva_cut_to_selection(selection, mva_cut);
   }
  double weight;
  double data_in_ctrl     = get_nevents(data, var, with_cut, variables, selection, bg_chains[trained_bg]->label, double_test_bg);
  double other_bg_in_ctrl = get_other_bg_in_ctrl(0,mc_weights, bg_chains, var, with_cut, variables, selection, trained_bg, double_test_bg);
//cout<<"in histoplot call bg in control, data"<< other_bg_in_ctrl<<","<<data_in_ctrl<<"\n";
  string nunu_selection;

  if(variables != NULL){nunu_selection = "((alljetsmetnomu_mindphi>2.0)&&(nvetomuons==0)&&(nvetoelectrons==0))*total_weight_lepveto";}
  else{ nunu_selection = "((alljetsmetnomu_mindphi>2.0)&&(classID==0)&&(nvetomuons==0)&&(nvetoelectrons==0))*total_weight_lepveto";}
  nunu_selection = HistoPlot::add_mva_cut_to_selection(nunu_selection, mva_cut);
   
  TH1F* nunu_h = HistoPlot::build_1d_histo(nunu_chain, var, true, false, "goff", variables, nunu_selection,1,mva_cut);
  double nunu_unweighted = nunu_h->Integral();
//cout<<"Zvv weighted total: "<<(data_in_ctrl-other_bg_in_ctrl)*5.651*1.513<<"\n";
//cout<<"Zvv unweighted total:"<<nunu_unweighted<<"\n";
  weight = (data_in_ctrl-other_bg_in_ctrl)*5.651/nunu_unweighted;
  return weight;
}
double MCWeights::calc_mc_weight(std::vector<double> mc_weights, int desired_bg_index,DataChain* data, std::vector<DataChain*> bg_chains, DataChain* bg_chain,Variable* var, 
bool with_cut, std::vector<Variable*>* variables, std::string mva_cut,int trained_bg, bool double_test_bg, bool if_parked)
{
//cout<<"bg to be calculated: !!!!!!!!!!!"<<bg_chain->label<<"!!!!!!!!!!!\n";
  std::string selection;
  if(variables != NULL){selection = get_mc_selection_str(bg_chain, var, variables, mva_cut, if_parked);}
  else{  
    if(!strcmp(bg_chain->label, "bg_wjets_tauv"))
    {
	selection = "((jetmetnomu_mindphi>1.0)&&(mht>20)&&(jet_csv2<0.2)&&(classID==0)&&" + bg_chain->lep_sel + ")" + "*total_weight_lepveto";
    }
    	else{selection = "((alljetsmetnomu_mindphi>2.0)&&(classID==0)&&"+bg_chain->lep_sel+")*total_weight_lepveto";}
    	selection = HistoPlot::add_mva_cut_to_selection(selection, mva_cut);//}
   }
  double weight;
  //cout<<"mc weight selection: "<<selection<<"\n";

  double data_in_ctrl     = get_nevents(data, var, with_cut, variables, selection, bg_chains[trained_bg]->label, double_test_bg);
//cout<<"data in control region: "<<data_in_ctrl<<"\n";

  TH1F* desired_h = HistoPlot::build_1d_histo(bg_chains[desired_bg_index], var, true, false, "goff", variables, selection,1);

  double ctrl_mc_in_ctrl  = desired_h->Integral();
//cout<<"desired background in control region: "<<ctrl_mc_in_ctrl<<"\n";

  double other_bg_in_ctrl = get_other_bg_in_ctrl(desired_bg_index,mc_weights, bg_chains, var, with_cut, variables, selection, trained_bg, double_test_bg);
//cout<<"other backgrounds in control region: "<<other_bg_in_ctrl<<"\n";
 if(ctrl_mc_in_ctrl!=0){
 weight = (data_in_ctrl - other_bg_in_ctrl) / ctrl_mc_in_ctrl;
  }
  else{weight=0;}
  return weight;
}

double MCWeights::calc_weight_error(std::vector<double> mc_weights,int desired_bg_index,DataChain* data, std::vector<DataChain*> bg_chains, DataChain* bg_chain,
Variable* var, bool with_cut, std::vector<Variable*>* variables, std::string mva_cut, bool if_parked)
{
  double sigma_MC_N_C;
  double sigma_data_N_C;
  double sigma_bg_N_C;
  std::string selection;
  if(desired_bg_index==0){mc_weights[0]=1;mc_weights[1]=1;mc_weights[2]=1;mc_weights[3]=1;mc_weights[6]=1;}
  else if(desired_bg_index==1){mc_weights[1]=1;mc_weights[2]=1;mc_weights[3]=1;}
  else if(desired_bg_index==2){mc_weights[2]=1;mc_weights[3]=1;}
  else if(desired_bg_index==3){mc_weights[3]=1;}


  if(variables != NULL){selection = get_mc_selection_str(bg_chain, var, variables, mva_cut,if_parked);}
  else
  { 
    if(!strcmp(bg_chain->label, "bg_wjets_tauv"))
    {
	selection = "((jetmetnomu_mindphi>1.0)&&(mht>20)&&(jet_csv2<0.2)&&(classID==0)&&" + bg_chain->lep_sel + ")" + "*total_weight_lepveto";
    }
    else
    {
	selection = "((alljetsmetnomu_mindphi>2.0)&&(classID==0)&&" + bg_chain->lep_sel + ")" + "*total_weight_lepveto";
    }
    selection = HistoPlot::add_mva_cut_to_selection(selection, mva_cut);
  }
  double weight_error;

  TH1F* data_h = HistoPlot::build_1d_histo(data, var, true, false, "goff", variables, selection);
  int nbins = data_h->GetSize();
  //cout<<"n bins of data: "<<nbins<<"\n";
  double data_N_C       = data_h->IntegralAndError(1,nbins,sigma_data_N_C,"");


  TH1F* MC_N_C_h = HistoPlot::build_1d_histo(bg_chain, var, true, false, "goff", variables, selection);
  //nbins = MC_N_C_h->GetSize();//bin no dependant on variable so this should be the same anfd they are
  //cout<<"n bins of bg: "<<nbins<<"\n";
  double MC_N_C = MC_N_C_h->IntegralAndError(1,nbins,sigma_MC_N_C,"");

  double sigma_bg_NC_sq;
  double sigma_bg;
  double bg_N_C;
  for(int i = 0; i<8;i++)
  {
	if(strcmp(bg_chain->label, bg_chains[i]->label))
	{
	  	std:: string bg_selection = HistoPlot::add_mc_to_selection(bg_chains[i], var, selection, mc_weights[i]);
		TH1F* bg_h = HistoPlot::build_1d_histo(bg_chains[i], var, true, false, "goff", variables, bg_selection);
		bg_N_C += bg_h->IntegralAndError(1,nbins,sigma_bg,"");
		sigma_bg_NC_sq += std::pow(sigma_bg,2);
	}
  }
  sigma_bg_N_C = std::pow(sigma_bg_NC_sq,0.5);
if(MC_N_C!=0){
  
  double err1           = sigma_data_N_C/MC_N_C;
  double err2           = sigma_bg_N_C/MC_N_C;
  double err3           = (data_N_C- bg_N_C)/(pow(MC_N_C,1.5));
  double error_sq       = std::pow(err1,2) + std::pow(err2,2) + std::pow(err3,2);
  weight_error   = std::pow(error_sq, 0.5);
}
else if(MC_N_C==0){//cout<<"Warning!!!\n"<<bg_chain->label<<"target background"<<" = "<<MC_N_C;
weight_error = 1000;}
//cout<<bg_chain->label<<"=============================weight error: "<<weight_error<<"\n";
  return weight_error;

}
double MCWeights::calc_nunu_weight_error(DataChain* data, std::vector<DataChain*> bg_chains, DataChain* nunu_chain,
Variable* var, bool with_cut, std::vector<Variable*>* variables, std::string mva_cut, bool if_parked)
{
  string selection;
  double mc_arr[8] = {1,1,1,1,1,1,1,1};
  std::vector<double> mc_weights (mc_arr, mc_arr + sizeof(mc_arr)/sizeof(mc_arr[0]));
  //TH1F* nunu_h = HistoPlot::build_1d_histo(nunu_chain, var, true, false, "goff", variables, "",1,mva_cut);
  //double nunu_unweighted = nunu_h->Integral();

  if(variables != NULL){selection = get_mc_selection_str(bg_chains[0], var, variables, mva_cut, if_parked);}
  else{  
    selection = "((alljetsmetnomu_mindphi>2.0)&&(classID==0)&&"+bg_chains[0]->lep_sel+")*total_weight_lepveto";
    	selection = HistoPlot::add_mva_cut_to_selection(selection, mva_cut);
   }
  double weight;
  double data_in_ctrl     = get_nevents(data, var, with_cut, variables, selection);
  double other_bg_in_ctrl = get_other_bg_in_ctrl(0,mc_weights, bg_chains, var, with_cut, variables, selection,0,false);
  double sigma_data_N_C =pow(data_in_ctrl,0.5);
  double sigma_bg_N_C = pow(other_bg_in_ctrl,0.5);

 /* double err1           = sigma_data_N_C/nunu_unweighted;
  double err2           = sigma_bg_N_C/nunu_unweighted;
  double err3           = (data_in_ctrl- other_bg_in_ctrl)/(pow(nunu_unweighted,1.5));
  double error_sq       = pow(err1,2)+ pow(err2,2) + pow(err3,2);*/
  double error_sq       = data_in_ctrl+other_bg_in_ctrl;

  double weight_error   = (std::pow(error_sq, 0.5)*5.651);//

//cout<<"nnunu mc weight error: "<<weight_error<<"\n";
  return weight_error;
}
