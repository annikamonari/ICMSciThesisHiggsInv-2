#include "../include/data_card.h"
#include "../include/mlp_analysis.h"
#include "../include/bdt_analysis.h"
#include <sstream>
#include <string>
using namespace std;
void DataCard::create_datacard(int bg_to_train, DataChain* data_chain, DataChain* signal_chain, std::vector<DataChain*> bg_chains,
                               Variable* var, bool with_cut, std::vector<Variable*>* variables, std::string output_graph_name,
																															std::string mva_cut)
{
	 std::vector<double> mc_weights = HistoPlot::mc_weights(data_chain, bg_chains, var, with_cut, variables, mva_cut);
	 std::fstream fs;
	 std::string data_card_name = get_data_card_name(output_graph_name, mva_cut);
	 //std::cout << data_card_name << std::endl;
	 fs.open (data_card_name.c_str(), std::fstream::out | std::fstream::trunc);
  int size = 1 + bg_chains.size();
//cout<<"opened text file\n";
  fs << imax_string();
  fs << jmax_string(size - 1);
  fs << kmax_string(size);
  fs << no_shape_line();
  fs << dashed_line();
  fs << bin_header_string();
//cout<<"about to get observation string\n";
  fs << bin_observation_string(get_total_nevents(bg_to_train, bg_chains, var, with_cut, NULL, mc_weights, mva_cut));
//cout<<"got observation string\n";  
  fs << dashed_line();
  fs << bin_grid_line(size);
  fs << process_labels(bg_chains, signal_chain);
  fs << process_2_string(process_line_2(size));
  fs << rate_string(get_rates(bg_to_train, data_chain, bg_chains, signal_chain, var, with_cut, NULL, mc_weights, mva_cut));
//cout<<"got rates string\n";
  fs << dashed_line();
//cout<<"got dashed line\n";
  fs << get_systematic_string(bg_to_train, data_chain, bg_chains, signal_chain, var, with_cut, NULL, mc_weights, mva_cut);
  //std::cout << "Data card created" << std::endl;
	 fs.close();
}

double DataCard::get_signal_error(DataChain* signal_chain, Variable* var, bool with_cut, std::vector<Variable*>* variables,
																																	 std::string mva_cut)
{
  std::string selection = "((alljetsmetnomu_mindphi>2.0)&&(nvetomuons==0)&&(nvetoelectrons==0))*total_weight_lepveto";
  selection = HistoPlot::add_classID_to_selection(selection, true);

  selection = HistoPlot::add_mva_cut_to_selection(selection, mva_cut);

  TH1F* signalh = HistoPlot::build_1d_histo(signal_chain, var, with_cut, false, "goff", variables, "", 1, mva_cut);
  double total_signal = 2 * HistoPlot::get_histo_integral(signalh, with_cut, var);
  double sig_sqrt = std::pow(total_signal, 0.5);

  return 1 + (sig_sqrt/total_signal);
}
//______________________________________________________________________________________________________________

std::vector<double> DataCard::get_bg_errors(int bg_to_train, DataChain* data, std::vector<DataChain*> bg_chains, DataChain* signal_chain,Variable* var, bool with_cut, std::vector<Variable*>* variables,
	std::vector<double> bg_mc_weights, std::string mva_cut)
{
	 double bg_errors_parsed[bg_chains.size()];

	 std::vector<double> bg_errors = HistoPlot::get_mc_weight_errors(bg_to_train, data, bg_chains, var, with_cut, variables, bg_mc_weights, mva_cut);

  std::vector<double> rates = get_rates(bg_to_train, data, bg_chains, signal_chain, var,with_cut, variables, bg_mc_weights, mva_cut);
  //std::cout << "got mc weight errors" << std::endl;
  for(int i = 0; i < bg_chains.size(); i++)
  {
    bg_errors_parsed[i] = 1 + (bg_errors[i] / rates[i+1]);
cout<<"parsed bg error in data for bg: "<<bg_chains[i]->label<<"= "<<bg_errors_parsed[i]<<endl;
  }
  std::vector<double> bg_error_vector (bg_errors_parsed, bg_errors_parsed + sizeof(bg_errors_parsed) / sizeof(bg_errors_parsed[0]));

  return bg_error_vector;
}
//______________________________________________________________________________________________________________
std::vector<double> DataCard::get_rates(int bg_to_train, DataChain* data, std::vector<DataChain*> bg_chains, DataChain* signal_chain,Variable* var, bool with_cut, std::vector<Variable*>* variables, std::vector<double> bg_mc_weights, std::string mva_cut)
{
  double rates[bg_chains.size() + 1];
  std::string selection1 = "((alljetsmetnomu_mindphi>2.0)&&(nvetomuons==0)&&(nvetoelectrons==0))*total_weight_lepveto";
  selection1 = HistoPlot::add_classID_to_selection(selection1, true);

  selection1 = HistoPlot::add_mva_cut_to_selection(selection1, mva_cut);
  //std::cout << "===signal" << std::endl;
  //std::cout << selection1 << std::endl;
  TH1F* signal_histo = HistoPlot::build_1d_histo(signal_chain, var, with_cut, false, "goff", variables, selection1, 1, mva_cut);
  rates[0] = 2*HistoPlot::get_histo_integral(signal_histo, with_cut, var);// taking into account test/train data split

    for(int i = 0; i < bg_chains.size();i++)
  {
    //std::cout << "=====" << bg_chains[i]->label << std::endl;
    std::string selection = "((alljetsmetnomu_mindphi>2.0)&&(nvetomuons==0)&&(nvetoelectrons==0))*total_weight_lepveto";
    selection = HistoPlot::add_classID_to_selection(selection, false);
    selection = HistoPlot::add_mc_to_selection(bg_chains[i],var , selection, bg_mc_weights[i]);
    selection = HistoPlot::add_mva_cut_to_selection(selection, mva_cut);
    //std::cout << bg_mc_weights[i] << std::endl;
    //std::cout << selection << std::endl;
		  TH1F* histo = HistoPlot::build_1d_histo(bg_chains[i], var, with_cut, false, "goff", variables, selection, bg_mc_weights[i], mva_cut);

		  double N = HistoPlot::get_histo_integral(histo, with_cut, var);
    if(!strcmp(bg_chains[i]->label, bg_chains[bg_to_train]->label)){N = 2*N;}// taking into account test/train data split
    rates[i + 1]= N;
    //std::cout << bg_chains[i]->label << " - " << N << std::endl;
  }
  std::vector<double> rates_vector (rates, rates + sizeof(rates) / sizeof(rates[0]));

  return rates_vector;
}
//______________________________________________________________________________________________________________

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
//______________________________________________________________________________________________________________

std::string DataCard::double_to_str(double sint)
{
  std::ostringstream ss;
  ss << sint;
  std::string str(ss.str());
  return str;
}
//______________________________________________________________________________________________________________

std::string DataCard::jmax_string(int jmax)
{
  std::string jmax_str = "jmax ";
  jmax_str += double_to_str(jmax);

  jmax_str += "  number of backgrounds \n";

  return jmax_str;
}
//______________________________________________________________________________________________________________

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

  return rate_str;
}
//______________________________________________________________________________________________________________

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

std::string DataCard::get_systematic_string(int bg_to_train, DataChain* data, std::vector<DataChain*> bg_chains,
																																												DataChain* signal_chain, Variable* var, bool with_cut, std::vector<Variable*>* variables,
																																												std::vector<double> bg_mc_weights, std::string mva_cut)
{
  double signal_error = get_signal_error(signal_chain, var, with_cut, variables, mva_cut);
//std::cout << "got signal erro" << std::endl;
  std::vector<double> bg_errors = get_bg_errors(bg_to_train, data, bg_chains, signal_chain, var, with_cut, variables, bg_mc_weights, mva_cut);
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

double DataCard::get_total_nevents(int bg_to_train, std::vector<DataChain*> bg_chains, Variable* var, bool with_cut, std::vector<Variable*>* variables,
																																			std::vector<double> bg_mc_weights, std::string mva_cut)
{
    std::string selection; 

	 double total = 0;
	 for (int i = 0; i < bg_chains.size(); i++)
	 	{
    selection = "((alljetsmetnomu_mindphi>2.0)&&(nvetomuons==0)&&(nvetoelectrons==0))*total_weight_lepveto";
    selection = HistoPlot::add_classID_to_selection(selection, false);
    selection = HistoPlot::add_mc_to_selection(bg_chains[i],var , selection, bg_mc_weights[i]);
    selection = HistoPlot::add_mva_cut_to_selection(selection, mva_cut);

	 		TH1F* histo = HistoPlot::build_1d_histo(bg_chains[i], var, with_cut, false, "goff", variables, selection, bg_mc_weights[i], mva_cut);
	 		double integral;
			integral = HistoPlot::get_histo_integral(histo, with_cut, var);

			if(!strcmp(bg_chains[i]->label, bg_chains[bg_to_train]->label)) {integral = integral*2;}

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

