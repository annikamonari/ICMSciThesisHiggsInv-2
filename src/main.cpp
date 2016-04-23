#include <initializer_list>
#include <cmath>
#include "../include/mva_analysis.h"
#include <vector>


//#include "../include/mlp_analysis.h"
using namespace std;
void produce_graphs(bool with_cut, const char* job_ptr) {
  SuperVars* super_vars             = new SuperVars();
  std::vector<Variable*> vars       = super_vars->get_discriminating_vars();
  std::vector<Variable*> cut_vars   = super_vars->get_signal_cut_vars();
  SuperChains* super_chains         = new SuperChains();
  std::vector<DataChain*> bg_chains = super_chains->get_bg_chains();
  DataChain* signal_chain           = super_chains->signal_chain;
  DataChain* data_chain             = super_chains->data_chain;
  DataChain* ewk_chain             = super_chains->bg_zll_ewk;
  DataChain* qcd_chain             = super_chains->bg_zll_qcd;

  std::vector<DataChain*> all_bg_chains = super_chains->get_all_bg_chains();
  std::vector<Variable*> parked_vars = super_vars->get_parked_vars();
   //const char* mva_type = "BDT";
  // the topmost folder for all root files so gitignore ignores properly
  std::string top_folder_name = "analysis";
  //const char* varying_params[] = {"NTrees", "AdaBoostBeta", "nCuts", "SeparationType"};
  const char* varying_params[] = {"HiddenLayers", "NCycles", "LearningRate","NeuronType"};
  // boolean is for whether or not to create separate output app files
  bool unique_output_files = false;
  // boolean is for whether or not to create datacards
  bool create_cards = true;
  std::string job_name = job_ptr;
  const char* jn = job_name.c_str();
  int counter = std::atoi(jn);
  std::string mva_cut = "";
  std::string method_name = "BDT";
  //int rel_bgs[] = {1, 2, 3, 6};
  std::string sign = ">"; // direction of cut
  int min = -80;//1769; // the minimum value you want cuts to be from
  int max =25;//1770; // max value you want cuts to be to
  double digits = 100; // number of digits + 1 of your cuts, e.g. if you put ur min as 40 then put 100 as digits to make it 0.4
/*
for(int i=0; i<vars.size();i++)
{
	string file_name;
	if(with_cut){file_name = "cuts/";}
	else{file_name = "nocuts/";}
	file_name.append(vars[i]->name);
	file_name.append(".png");
	HistoPlot::draw_plot(vars[i], bg_chains,signal_chain, data_chain, with_cut,&cut_vars,  true, file_name,mva_cut);
}*/


//double e_f = HistoPlot::get_efficiency_factor(ewk_chain,qcd_chain,cut_vars[0], &cut_vars, mva_cut);
//double error_on_e_f = HistoPlot::get_error_on_efficiency_factor(ewk_chain,qcd_chain,cut_vars[0], &cut_vars, mva_cut);
//cout<<e_f<<"+/-"<<error_on_e_f<<"\n";

/*for(int i=0;i<61;i++)
{
	DataCard::create_card_from_MC_weights_file("MLP","MLP_extrapolated_weights.csv",i,true);
}*/


/*TFile* trained_file;
trained_file = MVAAnalysis::get_mva_results(bg_chains, 6, signal_chain, data_chain,ewk_chain,qcd_chain, super_vars, "test", method_name,
  "300", BoostType[0], "0.2", "GiniIndex", "-1", NeuronType[0], 
  "800", HiddenLayers[1], LearningRate[1],unique_output_files, create_cards, job_name, mva_cut, sign, min, max, digits);
*/
//  double mc_weights_arr[] = {};
/*std::vector<double> mc_weights_vector = HistoPlot::mc_weights(data_chain, bg_chains, cut_vars[0], true, &cut_vars);
  string selection = "((alljetsmetnomu_mindphi>2.0) &&(nvetomuons==0)&&(nvetoelectrons==0)&&(jet1_pt>50.0)&&(jet2_pt>45.0)&&(metnomu_significance>3.5)&&(dijet_deta>4.2))";

    TH1F* histo = HistoPlot::build_1d_histo(data_chain, cut_vars[0], with_cut, false, "goff", &cut_vars,selection, 1, mva_cut);
    double integral = histo->Integral();
    cout<<data_chain->label<<": "<<integral<<endl;

for(int i=0; i<8;i++){
    TH1F* histo_bg = HistoPlot::build_1d_histo(bg_chains[i], cut_vars[0], with_cut, false, "goff", &cut_vars,"", mc_weights_vector[i], mva_cut);
    integral = histo_bg->Integral();
    cout<<bg_chains[i]->label<<": "<<integral<<endl;
}*/


bool not_tau=true;
int var_index[] ={7,11,12,5,6,8};// forward_tag_eta, jet1_pt, jet2_pt,alljetsmetnomu_mindphi, metnomu_significance, dijet_deta
for(int i=3;i<4;i++)
{
	if(i==3){not_tau=false;}
	string file_name;
	for(int j=0;j<6;j++)
	{
		file_name = vars[var_index[j]]->name;
		file_name.append("_");
		file_name.append(bg_chains[i]->label);
		file_name.append("_control.png");
		HistoPlot::plot_control(not_tau, vars[var_index[j]], data_chain,  bg_chains, &cut_vars, file_name, bg_chains[i]->lep_sel , "");
	}
}
/*for( int i = 6;i < 7;i++)
{
  cout<<vars[i]->name<<endl;
        string plot_name = "cuts/";
	plot_name.append(vars[i]->name);
	plot_name.append(".png");
	HistoPlot::draw_plot(vars[i], bg_chains,signal_chain, data_chain, true, &cut_vars , true, plot_name);
}
cout<<"plotted graph\n";
*/
/*std::vector<double> mc_weights_vector = HistoPlot::mc_weights(data_chain, bg_chains, cut_vars[0], true, &cut_vars);
cout<<"got regular mc weights\n";

std::vector<double> parked_mc_weights_vector = HistoPlot::mc_weights(data_chain, bg_chains, parked_vars[0], true, &parked_vars,"",6,false,true);
cout<<"got park mc weights\n";

double integral;double total=0;
cout<<"===================Gordon and Monari's preselection==============\n";
for(int i=0; i<8;i++){
    TH1F* histo = HistoPlot::build_1d_histo(bg_chains[i], cut_vars[0], with_cut, false, "goff", &cut_vars,"", mc_weights_vector[i], mva_cut);
    integral = HistoPlot::get_histo_integral(histo, with_cut, vars[0]);
    cout<<bg_chains[i]->label<<": "<<integral<<endl;
    total+=integral;
}
 cout << "total background = "<< total <<endl;

  TH1F* signal_histo = HistoPlot::build_1d_histo(signal_chain, cut_vars[0], with_cut, false, "goff", &cut_vars, "");
 cout << "total signal = "<< HistoPlot::get_histo_integral(signal_histo, with_cut, vars[0]) << endl;// taking into account test/train data split

cout<<"===================CMS parked preselection====================\n";
double total,integral;
for(int i=0; i<8;i++){
    TH1F* parked_histo = HistoPlot::build_parked_histo(bg_chains[i], parked_vars[7],&parked_vars, parked_mc_weights_vector[i]);
    integral = HistoPlot::get_histo_integral(parked_histo, with_cut, parked_vars[7]);
    cout<<bg_chains[i]->label<<": "<<integral<<endl;
    total+=integral;
}
 cout << "total background = "<< total <<endl;


// cout << "total signal = "<< HistoPlot::get_histo_integral(parked_signal_histo, with_cut, parked_vars[7]) << endl;// taking into account test/train data s 
*/



  /*TH1F* parked_signal_histo = HistoPlot::build_parked_histo(signal_chain, parked_vars[7], &parked_vars,1);

DataCard::create_datacard(5, data_chain, signal_chain, bg_chains, parked_vars[0], with_cut, &parked_vars, "parked_preselection.png");*/
/*double z = MCWeights::calc_mc_weight(data_chain, bg_chains, bg_chains[0], cut_vars[0], with_cut, &cut_vars, mva_cut, 1, false,false);
double z_error = MCWeights::calc_weight_error( data_chain, bg_chains, bg_chains[0], cut_vars[0], with_cut, &cut_vars, 1, false,  mva_cut);
cout<<"zll: "<<z<<endl;
cout<<"zll error ="<<z_error<<endl;


double enu = MCWeights::calc_mc_weight(data_chain, bg_chains, bg_chains[1], cut_vars[0], with_cut, &cut_vars, mva_cut, 1, false,false);
double enu_error = MCWeights::calc_weight_error( data_chain, bg_chains, bg_chains[1], cut_vars[0], with_cut, &cut_vars, 1, false,  mva_cut);
cout<<"enu: "<<enu<<endl;
cout<<"e error ="<<enu_error<<endl;

double munu = MCWeights::calc_mc_weight(data_chain, bg_chains, bg_chains[2], cut_vars[0], with_cut, &cut_vars, mva_cut, 1, false,false);
double munu_error = MCWeights::calc_weight_error( data_chain, bg_chains, bg_chains[2], cut_vars[0], with_cut, &cut_vars, 1, false,  mva_cut);
cout<<"munu: "<<enu<<endl;
cout<<"mu error ="<<enu_error<<endl;

double taunu = MCWeights::calc_mc_weight(data_chain, bg_chains, bg_chains[3], cut_vars[0], with_cut, &cut_vars, mva_cut, 1, false,false);
double taunu_error = MCWeights::calc_weight_error( data_chain, bg_chains, bg_chains[3], cut_vars[0], with_cut, &cut_vars, 1, false,  mva_cut);
cout<<"taunu: "<<taunu<<endl;
cout<<"tau error ="<<enu_error<<endl;
*/
}

int main(int argc, char** argv) {
TApplication theApp("tapp", &argc, argv);
/*string fline = DataCard::get_line_from_file("weights.txt", 0);
std::vector<string> line_vector = DataCard::get_vector_from_line(fline);
double total = DataCard::get_total_events_from_line(line_vector);
std::vector<double> rates = DataCard::get_rates_from_line(line_vector);
std::vector<double> errors = DataCard::get_errors_from_line(line_vector);
for(int i=0;i<9;i++)
{
	cout<<rates[i]<<","<<errors[i]<<"\n";
}

cout<<"total event: "<<total<<"\n";
cout<<line_vector<<"\n";*/
 produce_graphs(true, argv[1]);
  theApp.Run();
  return 0;
}
/*const char * files[] = {
"test/MLP-bg_zjets_vv-NeuronType=radial-NCycles=500-HiddenLayers=2,4,8,16,32,64-LearningRate=0.0001-EstimatorType=CE-50bins.root",
"test/MLP-bg_zjets_vv-NeuronType=radial-NCycles=500-HiddenLayers=2,4,8,16,32,64-LearningRate=0.001-EstimatorType=CE-50bins.root",
"test/MLP-bg_zjets_vv-NeuronType=radial-NCycles=500-HiddenLayers=2,4,8,16,32,64-LearningRate=0.001-EstimatorType=CE-GDN.root",
"test/MLP-bg_zjets_vv-NeuronType=radial-NCycles=500-HiddenLayers=2,4,8,16,32,64-LearningRate=0.001-EstimatorType=CE-GN.root",
"test/MLP-bg_zjets_vv-NeuronType=radial-NCycles=500-HiddenLayers=2,4,8,16,32,64-LearningRate=0.001-EstimatorType=CE-UDN.root",
"test/MLP-bg_zjets_vv-NeuronType=radial-NCycles=500-HiddenLayers=2,4,8,16,32,64-LearningRate=0.001-EstimatorType=CE-UGDN.root",
"test/MLP-bg_zjets_vv-NeuronType=radial-NCycles=500-HiddenLayers=2,4,8,16,32,64-LearningRate=0.001-EstimatorType=CE-UN.root",
"test/MLP-bg_zjets_vv-NeuronType=radial-NCycles=500-HiddenLayers=2,4,8,16,32,64-LearningRate=0.01-EstimatorType=CE-50bins.root",
"test/MLP-bg_zjets_vv-NeuronType=radial-NCycles=500-HiddenLayers=2,4,8,16,32,64-LearningRate=0.01-EstimatorType=CE-GDN.root",
"test/MLP-bg_zjets_vv-NeuronType=radial-NCycles=500-HiddenLayers=2,4,8,16,32,64-LearningRate=0.01-EstimatorType=CE-GN.root",
"test/MLP-bg_zjets_vv-NeuronType=radial-NCycles=500-HiddenLayers=2,4,8,16,32,64-LearningRate=0.01-EstimatorType=CE-GPN.root",
"test/MLP-bg_zjets_vv-NeuronType=radial-NCycles=500-HiddenLayers=2,4,8,16,32,64-LearningRate=0.01-EstimatorType=CE-UDN.root",
"test/MLP-bg_zjets_vv-NeuronType=radial-NCycles=500-HiddenLayers=2,4,8,16,32,64-LearningRate=0.01-EstimatorType=CE-UGDN.root",
"test/MLP-bg_zjets_vv-NeuronType=radial-NCycles=500-HiddenLayers=2,4,8,16,32,64-LearningRate=0.01-EstimatorType=CE-UGPN.root",
"test/MLP-bg_zjets_vv-NeuronType=radial-NCycles=500-HiddenLayers=2,4,8,16,32,64-LearningRate=0.01-EstimatorType=CE-UN.root",
"test/MLP-bg_zjets_vv-NeuronType=radial-NCycles=500-HiddenLayers=2,4,8,16,32,64-LearningRate=0.01-EstimatorType=CE-UPN.root",
"test/MLP-bg_zjets_vv-NeuronType=radial-NCycles=500-HiddenLayers=2,4,8,16,32-LearningRate=0.0001-EstimatorType=CE-50bins.root",
"test/MLP-bg_zjets_vv-NeuronType=radial-NCycles=500-HiddenLayers=2,4,8,16,32-LearningRate=0.0001-EstimatorType=CE-GDN.root",
"test/MLP-bg_zjets_vv-NeuronType=radial-NCycles=500-HiddenLayers=2,4,8,16,32-LearningRate=0.0001-EstimatorType=CE-GN.root",
"test/MLP-bg_zjets_vv-NeuronType=radial-NCycles=500-HiddenLayers=2,4,8,16,32-LearningRate=0.0001-EstimatorType=CE-UDN.root",
"test/MLP-bg_zjets_vv-NeuronType=radial-NCycles=500-HiddenLayers=2,4,8,16,32-LearningRate=0.001-EstimatorType=CE-50bins.root",
"test/MLP-bg_zjets_vv-NeuronType=radial-NCycles=500-HiddenLayers=2,4,8,16,32-LearningRate=0.001-EstimatorType=CE-GDN.root",
"test/MLP-bg_zjets_vv-NeuronType=radial-NCycles=500-HiddenLayers=2,4,8,16,32-LearningRate=0.001-EstimatorType=CE-GN.root",
"test/MLP-bg_zjets_vv-NeuronType=radial-NCycles=500-HiddenLayers=2,4,8,16,32-LearningRate=0.001-EstimatorType=CE-GPN.root",
"test/MLP-bg_zjets_vv-NeuronType=radial-NCycles=500-HiddenLayers=2,4,8,16,32-LearningRate=0.001-EstimatorType=CE-UDN.root",
"test/MLP-bg_zjets_vv-NeuronType=radial-NCycles=500-HiddenLayers=2,4,8,16,32-LearningRate=0.001-EstimatorType=CE-UGDN.root",
"test/MLP-bg_zjets_vv-NeuronType=radial-NCycles=500-HiddenLayers=2,4,8,16,32-LearningRate=0.001-EstimatorType=CE-UGPN.root",
"test/MLP-bg_zjets_vv-NeuronType=radial-NCycles=500-HiddenLayers=2,4,8,16,32-LearningRate=0.001-EstimatorType=CE-UN.root",
"test/MLP-bg_zjets_vv-NeuronType=radial-NCycles=500-HiddenLayers=2,4,8,16,32-LearningRate=0.001-EstimatorType=CE-UPN.root",
"test/MLP-bg_zjets_vv-NeuronType=radial-NCycles=500-HiddenLayers=2,4,8,16,32-LearningRate=0.01-EstimatorType=CE-50bins.root",
"test/MLP-bg_zjets_vv-NeuronType=radial-NCycles=500-HiddenLayers=2,4,8,16,32-LearningRate=0.01-EstimatorType=CE-GDN.root",
"test/MLP-bg_zjets_vv-NeuronType=radial-NCycles=500-HiddenLayers=2,4,8,16,32-LearningRate=0.01-EstimatorType=CE-GN.root",
"test/MLP-bg_zjets_vv-NeuronType=radial-NCycles=500-HiddenLayers=2,4,8,16,32-LearningRate=0.01-EstimatorType=CE-GPN.root",
"test/MLP-bg_zjets_vv-NeuronType=radial-NCycles=500-HiddenLayers=2,4,8,16,32-LearningRate=0.01-EstimatorType=CE-UDN.root",
"test/MLP-bg_zjets_vv-NeuronType=radial-NCycles=500-HiddenLayers=2,4,8,16,32-LearningRate=0.01-EstimatorType=CE-UGDN.root",
"test/MLP-bg_zjets_vv-NeuronType=radial-NCycles=500-HiddenLayers=2,4,8,16,32-LearningRate=0.01-EstimatorType=CE-UGPN.root",
"test/MLP-bg_zjets_vv-NeuronType=radial-NCycles=500-HiddenLayers=2,4,8,16,32-LearningRate=0.01-EstimatorType=CE-UN.root",
"test/MLP-bg_zjets_vv-NeuronType=radial-NCycles=500-HiddenLayers=2,4,8,16,32-LearningRate=0.01-EstimatorType=CE-UPN.root",
"test/MLP-bg_zjets_vv-NeuronType=radial-NCycles=500-HiddenLayers=2,4,8,16-LearningRate=0.0001-EstimatorType=CE-50bins.root",
"test/MLP-bg_zjets_vv-NeuronType=radial-NCycles=500-HiddenLayers=2,4,8,16-LearningRate=0.0001-EstimatorType=CE-GDN.root",
"test/MLP-bg_zjets_vv-NeuronType=radial-NCycles=500-HiddenLayers=2,4,8,16-LearningRate=0.0001-EstimatorType=CE-GN.root",
"test/MLP-bg_zjets_vv-NeuronType=radial-NCycles=500-HiddenLayers=2,4,8,16-LearningRate=0.0001-EstimatorType=CE-UDN.root",
"test/MLP-bg_zjets_vv-NeuronType=radial-NCycles=500-HiddenLayers=2,4,8,16-LearningRate=0.0001-EstimatorType=CE-UN.root",
"test/MLP-bg_zjets_vv-NeuronType=radial-NCycles=500-HiddenLayers=2,4,8,16-LearningRate=0.001-EstimatorType=CE-50bins.root",
"test/MLP-bg_zjets_vv-NeuronType=radial-NCycles=500-HiddenLayers=2,4,8,16-LearningRate=0.001-EstimatorType=CE-GDN.root",
"test/MLP-bg_zjets_vv-NeuronType=radial-NCycles=500-HiddenLayers=2,4,8,16-LearningRate=0.001-EstimatorType=CE-GN.root",
"test/MLP-bg_zjets_vv-NeuronType=radial-NCycles=500-HiddenLayers=2,4,8,16-LearningRate=0.001-EstimatorType=CE-GPN.root",
"test/MLP-bg_zjets_vv-NeuronType=radial-NCycles=500-HiddenLayers=2,4,8,16-LearningRate=0.001-EstimatorType=CE-UDN.root",
"test/MLP-bg_zjets_vv-NeuronType=radial-NCycles=500-HiddenLayers=2,4,8,16-LearningRate=0.001-EstimatorType=CE-UGDN.root",
"test/MLP-bg_zjets_vv-NeuronType=radial-NCycles=500-HiddenLayers=2,4,8,16-LearningRate=0.001-EstimatorType=CE-UGPN.root",
"test/MLP-bg_zjets_vv-NeuronType=radial-NCycles=500-HiddenLayers=2,4,8,16-LearningRate=0.001-EstimatorType=CE-UN.root",
"test/MLP-bg_zjets_vv-NeuronType=radial-NCycles=500-HiddenLayers=2,4,8,16-LearningRate=0.001-EstimatorType=CE-UPN.root",
"test/MLP-bg_zjets_vv-NeuronType=radial-NCycles=500-HiddenLayers=2,4,8,16-LearningRate=0.01-EstimatorType=CE-50bins.root",
"test/MLP-bg_zjets_vv-NeuronType=radial-NCycles=500-HiddenLayers=2,4,8,16-LearningRate=0.01-EstimatorType=CE-GDN.root",
"test/MLP-bg_zjets_vv-NeuronType=radial-NCycles=500-HiddenLayers=2,4,8,16-LearningRate=0.01-EstimatorType=CE-GN.root",
"test/MLP-bg_zjets_vv-NeuronType=radial-NCycles=500-HiddenLayers=2,4,8,16-LearningRate=0.01-EstimatorType=CE-GPN.root",
"test/MLP-bg_zjets_vv-NeuronType=radial-NCycles=500-HiddenLayers=2,4,8,16-LearningRate=0.01-EstimatorType=CE-UDN.root",
"test/MLP-bg_zjets_vv-NeuronType=radial-NCycles=500-HiddenLayers=2,4,8,16-LearningRate=0.01-EstimatorType=CE-UGDN.root",
"test/MLP-bg_zjets_vv-NeuronType=radial-NCycles=500-HiddenLayers=2,4,8,16-LearningRate=0.01-EstimatorType=CE-UGPN.root",
"test/MLP-bg_zjets_vv-NeuronType=radial-NCycles=500-HiddenLayers=2,4,8,16-LearningRate=0.01-EstimatorType=CE-UN.root",
"test/MLP-bg_zjets_vv-NeuronType=radial-NCycles=500-HiddenLayers=2,4,8,16-LearningRate=0.01-EstimatorType=CE-UPN.root",
"test/MLP-bg_zjets_vv-NeuronType=radial-NCycles=500-HiddenLayers=2,4,8-LearningRate=0.0001-EstimatorType=CE-50bins.root",
"test/MLP-bg_zjets_vv-NeuronType=radial-NCycles=500-HiddenLayers=2,4,8-LearningRate=0.0001-EstimatorType=CE-GDN.root",
"test/MLP-bg_zjets_vv-NeuronType=radial-NCycles=500-HiddenLayers=2,4,8-LearningRate=0.0001-EstimatorType=CE-GN.root",
"test/MLP-bg_zjets_vv-NeuronType=radial-NCycles=500-HiddenLayers=2,4,8-LearningRate=0.0001-EstimatorType=CE-UDN.root",
"test/MLP-bg_zjets_vv-NeuronType=radial-NCycles=500-HiddenLayers=2,4,8-LearningRate=0.0001-EstimatorType=CE-UN.root",
"test/MLP-bg_zjets_vv-NeuronType=radial-NCycles=500-HiddenLayers=2,4,8-LearningRate=0.001-EstimatorType=CE-50bins.root",
"test/MLP-bg_zjets_vv-NeuronType=radial-NCycles=500-HiddenLayers=2,4,8-LearningRate=0.001-EstimatorType=CE-GDN.root",
"test/MLP-bg_zjets_vv-NeuronType=radial-NCycles=500-HiddenLayers=2,4,8-LearningRate=0.001-EstimatorType=CE-GN.root",
"test/MLP-bg_zjets_vv-NeuronType=radial-NCycles=500-HiddenLayers=2,4,8-LearningRate=0.001-EstimatorType=CE-GPN.root",
"test/MLP-bg_zjets_vv-NeuronType=radial-NCycles=500-HiddenLayers=2,4,8-LearningRate=0.001-EstimatorType=CE-UDN.root",
"test/MLP-bg_zjets_vv-NeuronType=radial-NCycles=500-HiddenLayers=2,4,8-LearningRate=0.001-EstimatorType=CE-UGDN.root",
"test/MLP-bg_zjets_vv-NeuronType=radial-NCycles=500-HiddenLayers=2,4,8-LearningRate=0.001-EstimatorType=CE-UGPN.root",
"test/MLP-bg_zjets_vv-NeuronType=radial-NCycles=500-HiddenLayers=2,4,8-LearningRate=0.001-EstimatorType=CE-UN.root",
"test/MLP-bg_zjets_vv-NeuronType=radial-NCycles=500-HiddenLayers=2,4,8-LearningRate=0.001-EstimatorType=CE-UPN.root",
"test/MLP-bg_zjets_vv-NeuronType=radial-NCycles=500-HiddenLayers=2,4,8-LearningRate=0.01-EstimatorType=CE-50bins.root",
"test/MLP-bg_zjets_vv-NeuronType=radial-NCycles=500-HiddenLayers=2,4,8-LearningRate=0.01-EstimatorType=CE-GDN.root",
"test/MLP-bg_zjets_vv-NeuronType=radial-NCycles=500-HiddenLayers=2,4,8-LearningRate=0.01-EstimatorType=CE-GN.root",
"test/MLP-bg_zjets_vv-NeuronType=radial-NCycles=500-HiddenLayers=2,4,8-LearningRate=0.01-EstimatorType=CE-GPN.root",
"test/MLP-bg_zjets_vv-NeuronType=radial-NCycles=500-HiddenLayers=2,4,8-LearningRate=0.01-EstimatorType=CE-UDN.root",
"test/MLP-bg_zjets_vv-NeuronType=radial-NCycles=500-HiddenLayers=2,4,8-LearningRate=0.01-EstimatorType=CE-UGDN.root",
"test/MLP-bg_zjets_vv-NeuronType=radial-NCycles=500-HiddenLayers=2,4,8-LearningRate=0.01-EstimatorType=CE-UGPN.root",
"test/MLP-bg_zjets_vv-NeuronType=radial-NCycles=500-HiddenLayers=2,4,8-LearningRate=0.01-EstimatorType=CE-UN.root",
"test/MLP-bg_zjets_vv-NeuronType=radial-NCycles=500-HiddenLayers=2,4,8-LearningRate=0.01-EstimatorType=CE-UPN.root",
"test/MLP-bg_zjets_vv-NeuronType=radial-NCycles=500-HiddenLayers=2,4-LearningRate=0.0001-EstimatorType=CE-50bins.root",
"test/MLP-bg_zjets_vv-NeuronType=radial-NCycles=500-HiddenLayers=2,4-LearningRate=0.0001-EstimatorType=CE-GDN.root",
"test/MLP-bg_zjets_vv-NeuronType=radial-NCycles=500-HiddenLayers=2,4-LearningRate=0.0001-EstimatorType=CE-GN.root",
"test/MLP-bg_zjets_vv-NeuronType=radial-NCycles=500-HiddenLayers=2,4-LearningRate=0.0001-EstimatorType=CE-UDN.root",
"test/MLP-bg_zjets_vv-NeuronType=radial-NCycles=500-HiddenLayers=2,4-LearningRate=0.0001-EstimatorType=CE-UN.root",
"test/MLP-bg_zjets_vv-NeuronType=radial-NCycles=500-HiddenLayers=2,4-LearningRate=0.001-EstimatorType=CE-50bins.root",
"test/MLP-bg_zjets_vv-NeuronType=radial-NCycles=500-HiddenLayers=2,4-LearningRate=0.001-EstimatorType=CE-GDN.root",
"test/MLP-bg_zjets_vv-NeuronType=radial-NCycles=500-HiddenLayers=2,4-LearningRate=0.001-EstimatorType=CE-GN.root",
"test/MLP-bg_zjets_vv-NeuronType=radial-NCycles=500-HiddenLayers=2,4-LearningRate=0.001-EstimatorType=CE-GPN.root",
"test/MLP-bg_zjets_vv-NeuronType=radial-NCycles=500-HiddenLayers=2,4-LearningRate=0.001-EstimatorType=CE-UDN.root",
"test/MLP-bg_zjets_vv-NeuronType=radial-NCycles=500-HiddenLayers=2,4-LearningRate=0.001-EstimatorType=CE-UGDN.root",
"test/MLP-bg_zjets_vv-NeuronType=radial-NCycles=500-HiddenLayers=2,4-LearningRate=0.001-EstimatorType=CE-UGPN.root",
"test/MLP-bg_zjets_vv-NeuronType=radial-NCycles=500-HiddenLayers=2,4-LearningRate=0.001-EstimatorType=CE-UN.root",
"test/MLP-bg_zjets_vv-NeuronType=radial-NCycles=500-HiddenLayers=2,4-LearningRate=0.001-EstimatorType=CE-UPN.root",
"test/MLP-bg_zjets_vv-NeuronType=radial-NCycles=500-HiddenLayers=2,4-LearningRate=0.01-EstimatorType=CE-50bins.root",
"test/MLP-bg_zjets_vv-NeuronType=radial-NCycles=500-HiddenLayers=2,4-LearningRate=0.01-EstimatorType=CE-GDN.root",
"test/MLP-bg_zjets_vv-NeuronType=radial-NCycles=500-HiddenLayers=2,4-LearningRate=0.01-EstimatorType=CE-GN.root",
"test/MLP-bg_zjets_vv-NeuronType=radial-NCycles=500-HiddenLayers=2,4-LearningRate=0.01-EstimatorType=CE-GPN.root",
"test/MLP-bg_zjets_vv-NeuronType=radial-NCycles=500-HiddenLayers=2,4-LearningRate=0.01-EstimatorType=CE-UDN.root",
"test/MLP-bg_zjets_vv-NeuronType=radial-NCycles=500-HiddenLayers=2,4-LearningRate=0.01-EstimatorType=CE-UGDN.root",
"test/MLP-bg_zjets_vv-NeuronType=radial-NCycles=500-HiddenLayers=2,4-LearningRate=0.01-EstimatorType=CE-UGPN.root",
"test/MLP-bg_zjets_vv-NeuronType=radial-NCycles=500-HiddenLayers=2,4-LearningRate=0.01-EstimatorType=CE-UN.root",
"test/MLP-bg_zjets_vv-NeuronType=radial-NCycles=500-HiddenLayers=2,4-LearningRate=0.01-EstimatorType=CE-UPN.root",
"test/MLP-bg_zjets_vv-NeuronType=radial-NCycles=500-HiddenLayers=2-LearningRate=0.0001-EstimatorType=CE-50bins.root",
"test/MLP-bg_zjets_vv-NeuronType=radial-NCycles=500-HiddenLayers=2-LearningRate=0.0001-EstimatorType=CE-GDN.root",
"test/MLP-bg_zjets_vv-NeuronType=radial-NCycles=500-HiddenLayers=2-LearningRate=0.0001-EstimatorType=CE-GN.root",
"test/MLP-bg_zjets_vv-NeuronType=radial-NCycles=500-HiddenLayers=2-LearningRate=0.0001-EstimatorType=CE-UDN.root",
"test/MLP-bg_zjets_vv-NeuronType=radial-NCycles=500-HiddenLayers=2-LearningRate=0.0001-EstimatorType=CE-UN.root",
"test/MLP-bg_zjets_vv-NeuronType=radial-NCycles=500-HiddenLayers=2-LearningRate=0.001-EstimatorType=CE-50bins.root",
"test/MLP-bg_zjets_vv-NeuronType=radial-NCycles=500-HiddenLayers=2-LearningRate=0.001-EstimatorType=CE-GDN.root",
"test/MLP-bg_zjets_vv-NeuronType=radial-NCycles=500-HiddenLayers=2-LearningRate=0.001-EstimatorType=CE-GN.root",
"test/MLP-bg_zjets_vv-NeuronType=radial-NCycles=500-HiddenLayers=2-LearningRate=0.001-EstimatorType=CE-GPN.root",
"test/MLP-bg_zjets_vv-NeuronType=radial-NCycles=500-HiddenLayers=2-LearningRate=0.001-EstimatorType=CE-UDN.root",
"test/MLP-bg_zjets_vv-NeuronType=radial-NCycles=500-HiddenLayers=2-LearningRate=0.001-EstimatorType=CE-UGDN.root",
"test/MLP-bg_zjets_vv-NeuronType=radial-NCycles=500-HiddenLayers=2-LearningRate=0.001-EstimatorType=CE-UGPN.root",
"test/MLP-bg_zjets_vv-NeuronType=radial-NCycles=500-HiddenLayers=2-LearningRate=0.001-EstimatorType=CE-UN.root",
"test/MLP-bg_zjets_vv-NeuronType=radial-NCycles=500-HiddenLayers=2-LearningRate=0.001-EstimatorType=CE-UPN.root",
"test/MLP-bg_zjets_vv-NeuronType=radial-NCycles=500-HiddenLayers=2-LearningRate=0.01-EstimatorType=CE-50bins.root",
"test/MLP-bg_zjets_vv-NeuronType=radial-NCycles=500-HiddenLayers=2-LearningRate=0.01-EstimatorType=CE-GDN.root",
"test/MLP-bg_zjets_vv-NeuronType=radial-NCycles=500-HiddenLayers=2-LearningRate=0.1-EstimatorType=CE-50bins.root",
"test/MLP-bg_zjets_vv-NeuronType=radial-NCycles=500-HiddenLayers=2-LearningRate=0.1-EstimatorType=CE-GDN.root",
"test/MLP-bg_zjets_vv-NeuronType=radial-NCycles=500-HiddenLayers=2-LearningRate=0.1-EstimatorType=CE-GN.root",
"test/MLP-bg_zjets_vv-NeuronType=radial-NCycles=500-HiddenLayers=2-LearningRate=0.1-EstimatorType=CE-GPN.root",
"test/MLP-bg_zjets_vv-NeuronType=radial-NCycles=500-HiddenLayers=2-LearningRate=0.1-EstimatorType=CE-UDN.root",
"test/MLP-bg_zjets_vv-NeuronType=radial-NCycles=500-HiddenLayers=2-LearningRate=0.1-EstimatorType=CE-UGDN.root",
"test/MLP-bg_zjets_vv-NeuronType=radial-NCycles=500-HiddenLayers=2-LearningRate=0.1-EstimatorType=CE-UGPN.root",
"test/MLP-bg_zjets_vv-NeuronType=radial-NCycles=500-HiddenLayers=2-LearningRate=0.1-EstimatorType=CE-UN.root",
"test/MLP-bg_zjets_vv-NeuronType=radial-NCycles=500-HiddenLayers=2-LearningRate=0.1-EstimatorType=CE-UPN.root",
"test/MLP-bg_zjets_vv-NeuronType=radial-NCycles=500-HiddenLayers=2,4-LearningRate=0.1-EstimatorType=CE-50bins.root",
"test/MLP-bg_zjets_vv-NeuronType=radial-NCycles=500-HiddenLayers=2,4-LearningRate=0.1-EstimatorType=CE-GDN.root",
"test/MLP-bg_zjets_vv-NeuronType=radial-NCycles=500-HiddenLayers=2,4-LearningRate=0.1-EstimatorType=CE-GN.root",
"test/MLP-bg_zjets_vv-NeuronType=radial-NCycles=500-HiddenLayers=2,4-LearningRate=0.1-EstimatorType=CE-GPN.root",
"test/MLP-bg_zjets_vv-NeuronType=radial-NCycles=500-HiddenLayers=2,4-LearningRate=0.1-EstimatorType=CE-UDN.root",
"test/MLP-bg_zjets_vv-NeuronType=radial-NCycles=500-HiddenLayers=2,4-LearningRate=0.1-EstimatorType=CE-UGDN.root",
"test/MLP-bg_zjets_vv-NeuronType=radial-NCycles=500-HiddenLayers=2,4-LearningRate=0.1-EstimatorType=CE-UGPN.root",
"test/MLP-bg_zjets_vv-NeuronType=radial-NCycles=500-HiddenLayers=2,4-LearningRate=0.1-EstimatorType=CE-UN.root",
"test/MLP-bg_zjets_vv-NeuronType=radial-NCycles=500-HiddenLayers=2,4-LearningRate=0.1-EstimatorType=CE-UPN.root",
"test/MLP-bg_zjets_vv-NeuronType=radial-NCycles=500-HiddenLayers=2,4,8-LearningRate=0.1-EstimatorType=CE-50bins.root",
"test/MLP-bg_zjets_vv-NeuronType=radial-NCycles=500-HiddenLayers=2,4,8-LearningRate=0.1-EstimatorType=CE-GDN.root",
"test/MLP-bg_zjets_vv-NeuronType=radial-NCycles=500-HiddenLayers=2,4,8-LearningRate=0.1-EstimatorType=CE-GN.root",
"test/MLP-bg_zjets_vv-NeuronType=radial-NCycles=500-HiddenLayers=2,4,8-LearningRate=0.1-EstimatorType=CE-GPN.root",
"test/MLP-bg_zjets_vv-NeuronType=radial-NCycles=500-HiddenLayers=2,4,8-LearningRate=0.1-EstimatorType=CE-UDN.root",
"test/MLP-bg_zjets_vv-NeuronType=radial-NCycles=500-HiddenLayers=2,4,8-LearningRate=0.1-EstimatorType=CE-UGDN.root",
"test/MLP-bg_zjets_vv-NeuronType=radial-NCycles=500-HiddenLayers=2,4,8-LearningRate=0.1-EstimatorType=CE-UGPN.root",
"test/MLP-bg_zjets_vv-NeuronType=radial-NCycles=500-HiddenLayers=2,4,8-LearningRate=0.1-EstimatorType=CE-UN.root",
"test/MLP-bg_zjets_vv-NeuronType=radial-NCycles=500-HiddenLayers=2,4,8-LearningRate=0.1-EstimatorType=CE-UPN.root",
"test/MLP-bg_zjets_vv-NeuronType=radial-NCycles=500-HiddenLayers=2,4,8,16-LearningRate=0.1-EstimatorType=CE-50bins.root",
"test/MLP-bg_zjets_vv-NeuronType=radial-NCycles=500-HiddenLayers=2,4,8,16-LearningRate=0.1-EstimatorType=CE-GDN.root",
"test/MLP-bg_zjets_vv-NeuronType=radial-NCycles=500-HiddenLayers=2,4,8,16-LearningRate=0.1-EstimatorType=CE-GN.root",
"test/MLP-bg_zjets_vv-NeuronType=radial-NCycles=500-HiddenLayers=2,4,8,16-LearningRate=0.1-EstimatorType=CE-GPN.root",
"test/MLP-bg_zjets_vv-NeuronType=radial-NCycles=500-HiddenLayers=2,4,8,16-LearningRate=0.1-EstimatorType=CE-UDN.root",
"test/MLP-bg_zjets_vv-NeuronType=radial-NCycles=500-HiddenLayers=2,4,8,16-LearningRate=0.1-EstimatorType=CE-UGDN.root",
"test/MLP-bg_zjets_vv-NeuronType=radial-NCycles=500-HiddenLayers=2,4,8,16-LearningRate=0.1-EstimatorType=CE-UGPN.root",
"test/MLP-bg_zjets_vv-NeuronType=radial-NCycles=500-HiddenLayers=2,4,8,16-LearningRate=0.1-EstimatorType=CE-UN.root",
"teMLP-bg_zjets_vv-NeuronType=sigmoid-NCycles=800-HiddenLayers=2,4-LearningRate=0.01-EstimatorType=CE-GN.root
MLP-bg_zjets_vv-NeuronType=tanh-NCycles=800-HiddenLayers=2,4-LearningRate=0.01-EstimatorType=CE-GN.rootst/MLP-bg_zjets_vv-NeuronType=radial-NCycles=500-HiddenLayers=2,4,8,16-LearningRate=0.1-EstimatorType=CE-UPN.root",
"test/MLP-bg_zjets_vv-NeuronType=radial-NCycles=500-HiddenLayers=2,4,8,16,32-LearningRate=0.1-EstimatorType=CE-50bins.root",
"test/MLP-bg_zjets_vv-NeuronType=radial-NCycles=500-HiddenLayers=2,4,8,16,32-LearningRate=0.1-EstimatorType=CE-GDN.root",
"test/MLP-bg_zjets_vv-NeuronType=radial-NCycles=500-HiddenLayers=2,4,8,16,32-LearningRate=0.1-EstimatorType=CE-GN.root",
"test/MLP-bg_zjets_vv-NeuronType=radial-NCycles=500-HiddenLayers=2,4,8,16,32-LearningRate=0.1-EstimatorType=CE-GPN.root",
"test/MLP-bg_zjets_vv-NeuronType=radial-NCycles=500-HiddenLayers=2,4,8,16,32-LearningRate=0.1-EstimatorType=CE-UDN.root",
"test/MLP-bg_zjets_vv-NeuronType=radial-NCycles=500-HiddenLayers=2,4,8,16,32-LearningRate=0.1-EstimatorType=CE-UGDN.root",
"test/MLP-bg_zjets_vv-NeuronType=radial-NCycles=500-HiddenLayers=2,4,8,16,32-LearningRate=0.1-EstimatorType=CE-UGPN.root",
"test/MLP-bg_zjets_vv-NeuronType=radial-NCycles=500-HiddenLayers=2,4,8,16,32-LearningRate=0.1-EstimatorType=CE-UN.root",
"test/MLP-bg_zjets_vv-NeuronType=radial-NCycles=500-HiddenLayers=2,4,8,16,32-LearningRate=0.1-EstimatorType=CE-UPN.root",
"test/MLP-bg_zjets_vv-NeuronType=radial-NCycles=500-HiddenLayers=2,4,8,16,32,64-LearningRate=0.1-EstimatorType=CE-50bins.root",
"test/MLP-bg_zjets_vv-NeuronType=radial-NCycles=800-HiddenLayers=2,4-LearningRate=0.01-EstimatorType=CE-GN.root",
"test/MLP-bg_zjets_vv-NeuronType=sigmoid-NCycles=800-HiddenLayers=2,4-LearningRate=0.01-EstimatorType=CE-GN.root",
"test/MLP-bg_zjets_vv-NeuronType=tanh-NCycles=800-HiddenLayers=2,4-LearningRate=0.01-EstimatorType=CE-GN.root"
};
std::vector<const char*> files_v (files, files + sizeof(files) / sizeof(const char*));
  MVAAnalysis::get_estimators(files_v);
*/
