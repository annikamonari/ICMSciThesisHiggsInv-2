#include "../include/data_card.h"
#include "../include/mlp_analysis.h"
#include "../include/bdt_analysis.h"
#include <sstream>
#include <string>
using namespace std;
void DataCard::create_datacard(int bg_to_train,DataChain* ewk_chain, DataChain* qcd_chain, DataChain* data_chain, DataChain* signal_chain, std::vector<DataChain*> bg_chains,
 Variable* var, bool with_cut, std::vector<Variable*>* variables, std::string output_graph_name,
 std::string mva_cut)
{
  std::vector<DataChain*> bg_chs = bg_chains;
  bg_chs[bg_to_train] = signal_chain;//if parked remove this line

  std::vector<double> mc_weights = HistoPlot::mc_weights(ewk_chain, qcd_chain,data_chain, bg_chains, var, with_cut,variables , mva_cut,1,false, false);//last var is if parked
  //std::fstream fs;
  //std::string data_card_name = get_data_card_name(output_graph_name, mva_cut);
  //std::cout << data_card_name << std::endl;
  //fs.open (data_card_name.c_str(), std::fstream::out | std::fstream::trunc);
  //int size = 1 + bg_chains.size();
  cout<<mva_cut<<"\n";
  //cout<<"opened text file\n";
  //fs << imax_string();
  //fs << jmax_string(size - 1);
  //fs << kmax_string(size);
  //fs << no_shape_line();
  //fs << dashed_line();
  //fs << bin_header_string();
  //cout<<"about to get observation string\n";
  double data = get_total_data_events(data_chain, var,with_cut, variables,mva_cut);
  //fs << bin_observation_string(data);//get_total_nevents(data_chain, bg_to_train, bg_chains, var, with_cut, variables, mc_weights, mva_cut));
//cout<<"data at opt val is: "<<get_total_data_events(data_chain, var,with_cut, variables,mva_cut);
  //cout<<"got observation string\n";  
  //fs << dashed_line();
  //fs << bin_grid_line(size);
  //fs << process_labels(bg_chains, signal_chain);
  //fs << process_2_string(process_line_2(size));
  //std::vector<double> rates_w = get_rates( ewk_chain,  qcd_chain,bg_to_train, data_chain, bg_chains, signal_chain, var, with_cut, variables, mc_weights, mva_cut);
  //string rate_str = rate_string(rates_w);
  //cout<<rate_str<<"\n";
  //fs <<  rate_str;
  //cout<<"got rates string\n";
  //fs << dashed_line();
  //cout<<"got dashed line\n";
  //fs << get_systematic_string( ewk_chain,  qcd_chain,bg_to_train, data_chain, bg_chains,bg_chs, signal_chain, var, with_cut, variables, mc_weights, mva_cut);
  //std::cout << "Data card created" << std::endl;
  //fs.close();  

  std::vector<double> rates_d = get_raw_rates( ewk_chain,  qcd_chain,bg_to_train, data_chain, bg_chains, signal_chain, var, with_cut, variables, mva_cut);
  create_weights_series( ewk_chain,  qcd_chain,data_chain,signal_chain,bg_chains, var, with_cut, variables, output_graph_name,mva_cut,rates_d, data);
  std::cout << "weights and errors added to series\n"; 
}
void DataCard::create_card_from_MC_weights_file(const char* weights_file_name,int cut_number, bool use_data)
{
  string fline = get_line_from_file(weights_file_name, cut_number);
  std::vector<string> line_vector = DataCard::get_vector_from_line(fline);
  double total = DataCard::get_total_events_from_line(line_vector,use_data);
  std::vector<double> rates = DataCard::get_rates_from_line(line_vector);
  std::vector<double> errors = DataCard::get_errors_from_line(line_vector);
  
  string mva_cut_formatted = double_to_str(cut_number);
  

  std::fstream fs;
  string intro ="";
  string cut_num_str;

  if(cut_number<100)
  {
  	intro = "0.";
	if(cut_number<10){ intro += "0";}
	cut_num_str = intro + mva_cut_formatted;
  }
  else
  {
	if(cut_number %100>9)
	{
		cut_num_str = "1."+double_to_str(cut_number % 100);
	}
	else
	{	
		cut_num_str = "1.0"+double_to_str(cut_number % 100);
	}
  }
  std::string data_card_name = "Opt-BDT-Data-output>" + cut_num_str + ".txt";
  std::cout << data_card_name << std::endl;
  fs.open (data_card_name.c_str(), std::fstream::out | std::fstream::trunc);
  int size = 9;
  //cout<<"opened text file\n";
  fs << imax_string();
  fs << jmax_string(size - 1);
  fs << kmax_string(size);
  fs << no_shape_line();
  fs << dashed_line();
  fs << bin_header_string();
  //cout<<"about to get observation string\n";
  fs << bin_observation_string(total);
  //cout<<"got observation string\n";  
  fs << dashed_line();
  fs << bin_grid_line(size);
  fs << "process   Signal   Z_ll   W_enu   W_munu   W_taunu   top   VV   Z_nunu   QCD\n";
  fs << process_2_string(process_line_2(size));
  string rate_str = rate_string(rates);
  //cout<<rate_str<<"\n";
  fs <<  rate_str;
  //cout<<"got rates string\n";
  fs << dashed_line();
  //cout<<"got dashed line\n";
  fs << get_uncertainties_string(get_uncertainty_vectors(1.0,errors));
  std::cout << "Data card created" << std::endl;
  fs.close();  
}

void DataCard::create_weights_series(DataChain* ewk_chain, DataChain* qcd_chain,DataChain* data_chain, DataChain* signal_chain, std::vector<DataChain*> bg_chains,
 Variable* var, bool with_cut, std::vector<Variable*>* variables, std::string output_graph_name,
 std::string mva_cut, std::vector<double> rates_d, double data)
{
	std::vector<double> mc_weights = HistoPlot::mc_weights(ewk_chain, qcd_chain,data_chain, bg_chains, var, with_cut, variables, mva_cut,1,false,false);
	bool if_parked=false;
        bool double_test_bg=true;
  	std::vector<DataChain*> bg_chs = bg_chains;

	double mc_weight_errors[5];
	mc_weight_errors[0]=MCWeights::calc_weight_error(mc_weights, 0, data_chain, bg_chains, bg_chains[0], // zll
			var, with_cut,  variables, mva_cut, if_parked)*1.513;// zll
	mc_weight_errors[1]=MCWeights::calc_nunu_weight_error(ewk_chain, qcd_chain,data_chain, bg_chains, bg_chains[6],
        var, with_cut, variables, mva_cut, if_parked);// znunu

	mc_weight_errors[2]=MCWeights::calc_weight_error(mc_weights, 1,data_chain, bg_chains, bg_chains[1], // W enu
			var, with_cut,  variables, mva_cut, if_parked);// W enu
	mc_weight_errors[3]=MCWeights::calc_weight_error(mc_weights, 2,data_chain, bg_chains, bg_chains[2], // W munu
			var, with_cut,  variables, mva_cut, if_parked);// W munu
	mc_weight_errors[4]=MCWeights::calc_weight_error(mc_weights, 3,data_chain, bg_chains, bg_chains[3], // W taunu
			var, with_cut,  variables, mva_cut, if_parked);// W taunu

	std::fstream fsw;
	fsw.open ("BDTweights.txt", std::fstream::in | std::fstream::out | std::fstream::app);
	fsw<<mva_cut;
	fsw<<","<<rates_d[0]<<",0";
	fsw<<","<<rates_d[1]<<",,"<<mc_weights[0]<<","<<mc_weight_errors[0];
	fsw<<","<<rates_d[2]<<",,"<<mc_weights[1]<<","<<mc_weight_errors[2];
	fsw<<","<<rates_d[3]<<",,"<<mc_weights[2]<<","<<mc_weight_errors[3];
	fsw<<","<<rates_d[4]<<",,"<<mc_weights[3]<<","<<mc_weight_errors[4];
        fsw<<","<<rates_d[5]<<",,1,0";
        fsw<<","<<rates_d[6]<<",,1,0";
        fsw<<","<<rates_d[7]<<",,"<<mc_weights[6]<<","<<mc_weight_errors[1];
        fsw<<","<<rates_d[8]<<",,1,0,"<<data;

	fsw<<"\n";
	fsw.close();  
}

string DataCard::get_line_from_file(const char* weights_file_name,int cut_number)
{
	int line_number = cut_number+ 2; //cut number of 1 corresponds to output>0.01 
        double total;
	std::ifstream ifs;
	ifs.open (weights_file_name, std::fstream::in);
	string fline;
	for(int i=0;i<line_number;i++)
	{
		getline(ifs,fline);
	}
	ifs.close();
	cout<<fline<<"\n";
	return fline;
}

std::vector<string> DataCard::get_vector_from_line(string line)
{
	std::vector<string> line_vector;
	string string_before_comma;
        int comma_position;

	for(int i=0;i<57;i++)
	{
		comma_position =line.find(",");
	//	cout<<"commas position"<<comma_position<<"\n";
		string_before_comma = line.substr(0,comma_position);
	//	cout<<string_before_comma<<"\n";
		line_vector.push_back(string_before_comma);
		line  = line.substr(line.find(",")+1);
	}
	return line_vector;
}
double DataCard::get_total_events_from_line(std::vector<string> line_vector, bool use_data)
{
	double total_events;
  if(use_data)
  {
  	total_events = atof(line_vector[56].c_str());
  }
  else{
	for(int i=1;i<6;i++)
	{
		total_events += atof(line_vector[7*i-4].c_str());
	}
//cout<<"this wimpy vector is has a length of : "<<line_vector.size()<<"\n";
	total_events += atof(line_vector[37].c_str());//VV
	total_events += atof(line_vector[43].c_str());// Z nunu
	total_events += atof(line_vector[50].c_str());// QCD
  }
	return total_events;
}
std::vector<double> DataCard::get_rates_from_line(std::vector<string> line_vector)
{
	double rates_arr[9];
	rates_arr[0] = atof(line_vector[1].c_str());
	for(int i=1;i<6;i++)
	{
		rates_arr[i] = atof(line_vector[7*i-4].c_str());
	}
	rates_arr[6] = atof(line_vector[37].c_str());//VV
	rates_arr[7] = atof(line_vector[43].c_str());// Z nunu
	rates_arr[8] = atof(line_vector[50].c_str());// QCD

	std::vector<double> rates (rates_arr, rates_arr + sizeof(rates_arr)/sizeof(double)); 
	return rates;
}

std::vector<double> DataCard::get_errors_from_line(std::vector<string> line_vector)
{
	double error_arr[8];

		for(int i=0;i<4;i++)
	{
		error_arr[i] = atof(line_vector[7*i+9].c_str());
	}
	error_arr[4] = atof(line_vector[36].c_str());//top
	error_arr[5] = atof(line_vector[42].c_str());//VV
	error_arr[6] = atof(line_vector[49].c_str());// Z nunu
	error_arr[7] = atof(line_vector[55].c_str());// QCD

	std::vector<double> errors (error_arr, error_arr + sizeof(error_arr)/sizeof(double));
	return errors;
}

double DataCard::get_signal_error(DataChain* signal_chain, Variable* var, bool with_cut, std::vector<Variable*>* variables,
  std::string mva_cut)
{
  std::string selection = "((alljetsmetnomu_mindphi>2.0)&&(nvetomuons==0)&&(nvetoelectrons==0))*total_weight_lepveto";
  selection = HistoPlot::add_classID_to_selection(selection, true);

  selection = HistoPlot::add_mva_cut_to_selection(selection, mva_cut);

  TH1F* signalh = HistoPlot::build_1d_histo(signal_chain, var, with_cut, false, "goff", variables, selection, 1, mva_cut);//add selection back in here
  //TH1F* signalh = HistoPlot::build_parked_histo(signal_chain, var, variables,1);
  double total_signal = 2*HistoPlot::get_histo_integral(signalh, with_cut, var);//
  //double total_signal = HistoPlot::get_histo_integral(signalh, with_cut, var);// parked

  double sig_sqrt = std::pow(total_signal, 0.5);
  double sig_err_formatted;
  if(total_signal!=0){sig_err_formatted = 1 + (sig_sqrt/total_signal);}
  else if (total_signal==0){sig_err_formatted=2;}
  return sig_err_formatted;
}

std::vector<double> DataCard::get_bg_errors(DataChain* ewk_chain, DataChain* qcd_chain,int bg_to_train, DataChain* data, std::vector<DataChain*> bg_chains,std::vector<DataChain*> bg_chs_tree, DataChain* signal_chain,
  Variable* var, bool with_cut, std::vector<Variable*>* variables,
  std::vector<double> bg_mc_weights, std::string mva_cut)
{
  double bg_errors_parsed[bg_chains.size()];

  std::vector<double> bg_errors = HistoPlot::get_bg_errors(bg_to_train, ewk_chain,  qcd_chain, data, bg_chs_tree, var, with_cut, variables, bg_mc_weights, mva_cut, 
            false);// last var is if parked
  std::vector<double> rates = get_rates(ewk_chain, qcd_chain,bg_to_train, data, bg_chs_tree, signal_chain, var,with_cut, variables, bg_mc_weights, mva_cut);
  std::cout << "got mc weight errors" << std::endl;
  for(int i = 0; i < bg_chains.size(); i++)
  {
    if(rates[i+1]!=0)
    {
      bg_errors_parsed[i] = 1 + (bg_errors[i] / rates[i+1]);
    }
    else if(rates[i+1]==0)
    {
      bg_errors_parsed[i] = 1;     
    }
  }
  std::vector<double> bg_error_vector (bg_errors_parsed, bg_errors_parsed + sizeof(bg_errors_parsed) / sizeof(bg_errors_parsed[0]));

  return bg_error_vector;
}
//_____________________________________________________________________________________________________________________________
std::vector<double> DataCard::get_raw_rates(DataChain* ewk_chain, DataChain* qcd_chain,int bg_to_train, DataChain* data, std::vector<DataChain*> bg_chains, DataChain* signal_chain,
 Variable* var, bool with_cut, std::vector<Variable*>* variables, std::string mva_cut)
{
  double rates[bg_chains.size() + 1];
  std::string selection1 = "((alljetsmetnomu_mindphi>2.0)&&(nvetomuons==0)&&(nvetoelectrons==0))*total_weight_lepveto";
  selection1 = HistoPlot::add_classID_to_selection(selection1, true);

  selection1 = HistoPlot::add_mva_cut_to_selection(selection1, mva_cut);
 // std::cout << "===signal" << std::endl;
  //std::cout << selection1 << std::endl;
  TH1F* signal_histo = HistoPlot::build_1d_histo(signal_chain, var, with_cut, false, "goff", variables, selection1, 1, mva_cut);
  //TH1F* signal_histo = HistoPlot::build_parked_histo(signal_chain, var, variables,1);
// multiply rates by 2
  rates[0] =2* signal_histo->Integral();// taking into account test/train data split
  //  rates[0] =   HistoPlot::get_histo_integral(signal_histo, with_cut, var);// parked taking into account test/train data split
//std::cout<<rates[0]<<"\n";
  for(int i = 0; i < bg_chains.size();i++)
  {
  if(i!=6){
    //std::cout << "=====" << bg_chains[i]->label << std::endl;
    std::string selection = "((alljetsmetnomu_mindphi>2.0)&&(nvetomuons==0)&&(nvetoelectrons==0))*total_weight_lepveto";
    selection = HistoPlot::add_classID_to_selection(selection, false);
    selection = HistoPlot::add_mva_cut_to_selection(selection, mva_cut);
    //std::cout <<"==========="<< bg_chains[i]->label<<"-mc weight: " << bg_mc_weights[i] << std::endl;
    //std::cout << selection << std::endl;
    TH1F* histo = HistoPlot::build_1d_histo(bg_chains[i], var, with_cut, false, "goff", variables, selection,1 , mva_cut);
    //TH1F* histo = HistoPlot::build_parked_histo(bg_chains[i], var, variables,bg_mc_weights[i]);
    double N = histo->Integral();
    if(!strcmp(bg_chains[i]->label, bg_chains[bg_to_train]->label)){N = 2*N;}// comment outif parked , *2 for taking into account test/train data split
    rates[i + 1]= N;
    //std::cout << bg_chains[i]->label << " -rate: " << N << std::endl;
   // cout<<"rates: "<<i+1<<", "<<N<<"\n";
  }
  else
	{
		double mc_arr[8] = {1,1,1,1,1,1,1,1};
  		std::vector<double> mc_weights (mc_arr, mc_arr + sizeof(mc_arr)/sizeof(mc_arr[0]));
  		string selection;
  		if(variables != NULL){selection = MCWeights::get_mc_selection_str(bg_chains[0], var, variables, mva_cut,false);}
 		else{  
    			selection = "((alljetsmetnomu_mindphi>2.0)&&(classID==0)&&"+bg_chains[0]->lep_sel+")*total_weight_lepveto";
    			selection = HistoPlot::add_mva_cut_to_selection(selection, mva_cut);
   		}
  		double weight;
  		double data_in_ctrl     = MCWeights::get_nevents(data, var, with_cut, variables, selection, bg_chains[bg_to_train]->label, false);
  		double other_bg_in_ctrl = MCWeights::get_other_bg_in_ctrl(0,mc_weights, bg_chains, var, 
		with_cut, variables, selection,bg_to_train , true);
		double e_f = HistoPlot::get_efficiency_factor(ewk_chain, qcd_chain,var, variables, mva_cut);
  		double N = (data_in_ctrl-other_bg_in_ctrl)*5.651*e_f;
		rates[i + 1]= N;
	//	cout<<"rates: "<<i+1<<", "<<N<<"\n";
	}
  }
  std::vector<double> rates_vector (rates, rates + sizeof(rates) / sizeof(rates[0]));
  return rates_vector;
}

//_____________________________________________________________________________________________________________________________


std::vector<double> DataCard::get_rates(DataChain* ewk_chain, DataChain* qcd_chain,int bg_to_train, DataChain* data, std::vector<DataChain*> bg_chains, DataChain* signal_chain,
 Variable* var, bool with_cut, std::vector<Variable*>* variables, std::vector<double> bg_mc_weights,
 std::string mva_cut)
{
  double rates[bg_chains.size() + 1];
  std::string selection1 = "((alljetsmetnomu_mindphi>2.0)&&(nvetomuons==0)&&(nvetoelectrons==0))*total_weight_lepveto";
  selection1 = HistoPlot::add_classID_to_selection(selection1, true);

  selection1 = HistoPlot::add_mva_cut_to_selection(selection1, mva_cut);
 // std::cout << "===signal" << std::endl;
  //std::cout << selection1 << std::endl;
  TH1F* signal_histo = HistoPlot::build_1d_histo(signal_chain, var, with_cut, false, "goff", variables, selection1, 1, mva_cut);
  //TH1F* signal_histo = HistoPlot::build_parked_histo(signal_chain, var, variables,1);
// multiply rates by 2
  rates[0] =2* HistoPlot::get_histo_integral(signal_histo, with_cut, var);// taking into account test/train data split
  //  rates[0] =   HistoPlot::get_histo_integral(signal_histo, with_cut, var);// parked taking into account test/train data split
//std::cout<<rates[0]<<"\n";
  for(int i = 0; i < bg_chains.size();i++)
  {
  if(i!=6){
    //std::cout << "=====" << bg_chains[i]->label << std::endl;
    std::string selection = "((alljetsmetnomu_mindphi>2.0)&&(nvetomuons==0)&&(nvetoelectrons==0))*total_weight_lepveto";
    selection = HistoPlot::add_classID_to_selection(selection, false);
    selection = HistoPlot::add_mc_to_selection(bg_chains[i],var , selection, bg_mc_weights[i]);
    selection = HistoPlot::add_mva_cut_to_selection(selection, mva_cut);
    //std::cout <<"==========="<< bg_chains[i]->label<<"-mc weight: " << bg_mc_weights[i] << std::endl;
    //std::cout << selection << std::endl;
    TH1F* histo = HistoPlot::build_1d_histo(bg_chains[i], var, with_cut, false, "goff", variables, selection, bg_mc_weights[i], mva_cut);
    //TH1F* histo = HistoPlot::build_parked_histo(bg_chains[i], var, variables,bg_mc_weights[i]);
    double N = HistoPlot::get_histo_integral(histo, with_cut, var);
    if(!strcmp(bg_chains[i]->label, bg_chains[bg_to_train]->label)){N = 2*N;}// comment outif parked , *2 for taking into account test/train data split
    rates[i + 1]= N;
    //std::cout << bg_chains[i]->label << " -rate: " << N << std::endl;
   // cout<<"rates: "<<i+1<<", "<<N<<"\n";
  }
  else
	{
		double mc_arr[8] = {1,1,1,1,1,1,1,1};
  		std::vector<double> mc_weights (mc_arr, mc_arr + sizeof(mc_arr)/sizeof(mc_arr[0]));
  		string selection;
  		if(variables != NULL){selection = MCWeights::get_mc_selection_str(bg_chains[0], var, variables, mva_cut,false);}
 		else{  
    			selection = "((alljetsmetnomu_mindphi>2.0)&&(classID==0)&&"+bg_chains[0]->lep_sel+")*total_weight_lepveto";
    			selection = HistoPlot::add_mva_cut_to_selection(selection, mva_cut);
   		}
  		double weight;
  		double data_in_ctrl     = MCWeights::get_nevents(data, var, with_cut, variables, selection, bg_chains[bg_to_train]->label, false);
  		double other_bg_in_ctrl = MCWeights::get_other_bg_in_ctrl(0,mc_weights, bg_chains, var, 
		with_cut, variables, selection,bg_to_train , true);
		double e_f = HistoPlot::get_efficiency_factor(ewk_chain, qcd_chain,var, variables, mva_cut);
  		double N = (data_in_ctrl-other_bg_in_ctrl)*5.651*e_f;
		rates[i + 1]= N;
	//	cout<<"rates: "<<i+1<<", "<<N<<"\n";
	}
  }
  std::vector<double> rates_vector (rates, rates + sizeof(rates) / sizeof(rates[0]));
//std::cout<<"bg zll rate is: "<<rates[1]<<"\n";
  return rates_vector;
}

std::vector<int> DataCard::process_line_2(int size)
{
  int p_line[size];
  for(int i=0;i<size;i++)
  {
   p_line[i] = i; 
 }
 std::vector<int> pl_vector (p_line, p_line + sizeof(p_line) / sizeof(p_line[0]));
 return pl_vector;

}

std::string DataCard::double_to_str(double sint)
{
  std::ostringstream ss;
  ss << sint;
  std::string str(ss.str());
  return str;
}
std::string DataCard::jmax_string(int jmax)
{
  std::string jmax_str = "jmax ";
  jmax_str += double_to_str(jmax);

  jmax_str += "  number of backgrounds \n";

  return jmax_str;
}

std::string DataCard::imax_string()
{
  return "imax 1  number of channels \n";
}

std::string DataCard::kmax_string(int kmax)
{
  std::string kmax_str = "kmax ";
  kmax_str += double_to_str(kmax);
  kmax_str += "  number of nuisance parameters (sources of systematical uncertainties) \n";
  return kmax_str;
}

std::string DataCard::bin_header_string()
{
  return "bin c1 \n";
}

std::string DataCard::bin_observation_string(int nbins)
{
  std::string bin_obs_str;
  bin_obs_str += "observation ";
  bin_obs_str += double_to_str(nbins);
  bin_obs_str += " \n";


  return bin_obs_str;
}

std::string DataCard::bin_grid_line(int cols)
{
  std::string bin_grid_str = "bin";
  for (int i = 0; i < cols; i++)
  {
    bin_grid_str += "   c1";
  }
  bin_grid_str += " \n";

  return bin_grid_str;
}

std::string DataCard::process_labels(std::vector<DataChain*> bg_chains, DataChain* signal_chain)
{
  std::string process_labels_str = "process   ";
  process_labels_str.append("signal");

  for (int i = 0; i < bg_chains.size(); i++)
  {
    process_labels_str += "   ";
    process_labels_str.append(bg_chains[i]->label);
  }
  process_labels_str += " \n";

  return process_labels_str;
}

std::string DataCard::dashed_line()
{
  return "------------ \n";
}

//signal is 0, backgrounds are all 1 onwards
std::string DataCard::process_2_string(std::vector<int> line_2_vals)
{
  std::string line_2 = "process";
  for (int i = 0; i < line_2_vals.size(); i++)
  {
    line_2 += "   ";
    line_2 += double_to_str(line_2_vals[i]);
  }
  line_2 += "\n";

  return line_2;
}

std::string DataCard::rate_string(std::vector<double> rates)
{
  std::string rate_str = "rate";
  for (int i = 0; i < rates.size(); i++)
  {
    rate_str += "   ";
    rate_str += double_to_str(rates[i]);
  }
  rate_str += " \n";
//std::cout<<"bg zll rate in rate_string is: "<<rates[1]<<"\n";
//std::cout<<"rate str:" <<rate_str<<"\n";
  return rate_str;
}

std::vector<std::vector<double> > DataCard::get_uncertainty_vectors(double signal_error, std::vector<double> bg_errors)
{
  int size = bg_errors.size() + 1;
  std::vector<std::vector<double> > error_lines(size, std::vector<double>(size));
  error_lines[0] = get_zeros(size);
  error_lines[0][0] = signal_error;

  for (int i = 0; i < bg_errors.size() ; i++)
  {
    error_lines[i+1] = get_zeros(size);
    error_lines[i+1][i+1] = bg_errors[i];
    //std::cout << error_lines[i+1][i+1] << std::endl;
  }

  return error_lines;
}

std::string DataCard::get_uncertainties_string(std::vector<std::vector<double> > uncertainty_vectors)
{
  std::string uncertainties = "";

  for (int i = 0; i < uncertainty_vectors.size(); i++)
  {
    uncertainties += ("uncertainty" + double_to_str(i) + " lnN  ");
    uncertainties += get_single_uncertainty_str(uncertainty_vectors[i]);
    uncertainties += "\n";
  }

  return uncertainties;
}

std::string DataCard::get_systematic_string(DataChain* ewk_chain, DataChain* qcd_chain,int bg_to_train, DataChain* data, std::vector<DataChain*> bg_chains,std::vector<DataChain*> bg_chs_tree,
  DataChain* signal_chain, Variable* var, bool with_cut, std::vector<Variable*>* variables,
  std::vector<double> bg_mc_weights, std::string mva_cut)
{
  double signal_error = 1.0;//get_signal_error(signal_chain, var, with_cut, variables, mva_cut);
  //std::cout << "got signal error" << std::endl;
  std::vector<double> bg_errors = get_bg_errors(ewk_chain, qcd_chain,bg_to_train, data, bg_chains,bg_chs_tree, signal_chain, var, with_cut, variables, bg_mc_weights, mva_cut);//if parked is last variable
  std::vector<std::vector<double> > uncertainty_vectors = DataCard::get_uncertainty_vectors(signal_error, bg_errors);

  return get_uncertainties_string(uncertainty_vectors);
}

std::string DataCard::get_single_uncertainty_str(std::vector<double> single_uncertainty_vector)
{
  std::string vector_str = "";
  for (int i = 0; i < single_uncertainty_vector.size(); i++)
  {
    vector_str += "   ";
    if (single_uncertainty_vector[i] == 0)
    {
      vector_str += "-";
    }
    else
    {
      vector_str += double_to_str(single_uncertainty_vector[i]);
    }
  }

  return vector_str;
}

std::vector<double> DataCard::get_zeros(int size)
{
  double zero_arr[size];
  for (int i = 0; i < size; i++)
  {
    zero_arr[i] = 0;
  }
  std::vector<double> zeros (zero_arr, zero_arr + sizeof(zero_arr) / sizeof(zero_arr[0]));

  return zeros;
}

std::string DataCard::no_shape_line()
{
  return "shapes *    c1  FAKE \n";
}

double DataCard::get_total_data_events(DataChain* data, Variable* var, bool with_cut, std::vector<Variable*>* variables,std::string mva_cut)
{
  std::string selection = "((alljetsmetnomu_mindphi>2.0)&&(nvetomuons==0)&&(nvetoelectrons==0))*total_weight_lepveto";
    selection = HistoPlot::add_classID_to_selection(selection, false);
//cout<<selection<<"\n";
  selection = HistoPlot::add_mva_cut_to_selection(selection, mva_cut);

    TH1F* histo = HistoPlot::build_1d_histo(data, var, with_cut, false, "goff", variables, selection, 1, mva_cut);

    double integral = histo->Integral();

  return integral;
}

double DataCard::get_total_nevents(DataChain* data, int bg_to_train, std::vector<DataChain*> bg_chains, Variable* var, bool with_cut, std::vector<Variable*>* variables,
 std::vector<double> bg_mc_weights, std::string mva_cut)
{
  std::string selection; 

  double total = 0;
  double integral;

  for (int i = 0; i < bg_chains.size(); i++)
  {
    if(i!=6){
    selection = "((alljetsmetnomu_mindphi>2.0)&&(nvetomuons==0)&&(nvetoelectrons==0))*total_weight_lepveto";
    selection = HistoPlot::add_classID_to_selection(selection, false);
    selection = HistoPlot::add_mc_to_selection(bg_chains[i],var , selection, bg_mc_weights[i]);
    selection = HistoPlot::add_mva_cut_to_selection(selection, mva_cut);

    TH1F* histo = HistoPlot::build_1d_histo(bg_chains[i], var, with_cut, false, "goff", variables, selection, bg_mc_weights[i], mva_cut);
     // TH1F* histo = HistoPlot::build_parked_histo(bg_chains[i], var, variables,bg_mc_weights[i]);

    integral =histo->Integral();

    if(i==bg_to_train) {integral = 2*integral;}// multiply this by 2 //if  parked comment out 
//cout<<integral<<"\n";
	}
	else
	{
		double mc_arr[8] = {1,1,1,1,1,1,1,1};
  		std::vector<double> mc_weights (mc_arr, mc_arr + sizeof(mc_arr)/sizeof(mc_arr[0]));
  		string selection;
  		if(variables != NULL){selection = MCWeights::get_mc_selection_str(bg_chains[0], var, variables, mva_cut,false);}
 		else{  
    			selection = "((alljetsmetnomu_mindphi>2.0)&&(classID==0)&&"+bg_chains[0]->lep_sel+")*total_weight_lepveto";
    			selection = HistoPlot::add_mva_cut_to_selection(selection, mva_cut);
   		}
  		double weight;
  		double data_in_ctrl     = MCWeights::get_nevents(data, var, with_cut, variables, selection, bg_chains[bg_to_train]->label, false);
  		double other_bg_in_ctrl = MCWeights::get_other_bg_in_ctrl(0,mc_weights, bg_chains, var, 
		with_cut, variables, selection,bg_to_train , true);

  		double N = (data_in_ctrl-other_bg_in_ctrl)*5.651*1.513;
		integral= N;
		//cout<<"rates: "<<i+1<<", "<<N<<"\n";
	}
    total += integral;
  }

  return total;
}

std::string DataCard::get_data_card_name(std::string output_graph_name, std::string mva_cut)
{

  std::string card_file_name;

  if (output_graph_name == "")
  {
    card_file_name = "preselection_only.txt";
  }
  else
  {
    card_file_name = HistoPlot::replace_all(output_graph_name, ".png", ".txt");
  }

  return card_file_name;
}

/*
int imax = 1;//bin number
int jmax = bg_chains.size(); //number of backgrounds
int kmax;//number of indepenmdant sources of uncertaintied
bool with_cut = true;*/

