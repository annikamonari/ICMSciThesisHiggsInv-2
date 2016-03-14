#include "../include/mva_analysis.h"
#include "../include/bdt_analysis.h"
#include "../include/mlp_analysis.h"

void MVAAnalysis::get_plots_varying_params(std::vector<DataChain*> bg_chains, int bg_to_train, DataChain* signal_chain, DataChain* data_chain,
																																											SuperVars* super_vars, std::string method_name, std::string dir_name,
																																											std::vector<const char*> NTrees, std::vector<const char*> BoostType,
																																											std::vector<const char*> AdaBoostBeta, std::vector<const char*> SeparationType,
																																											std::vector<const char*> nCuts, std::vector<const char*> NeuronType, std::vector<const char*> NCycles,
																																											std::vector<const char*> HiddenLayers, bool unique_output_files, bool create_cards, std::string job_name)
{
  std::vector<const char*> file_paths = vary_parameters(bg_chains, bg_to_train, signal_chain, data_chain, super_vars, method_name, dir_name,
																																															   NTrees, BoostType, AdaBoostBeta, SeparationType, nCuts, NeuronType, NCycles, HiddenLayers,
																																																		unique_output_files, create_cards, job_name);

  std::vector<TFile*> files =  get_files_from_paths(file_paths);

  std::string bg_label = bg_chains[bg_to_train]->label;
  std::string folder_name = "analysis/" + method_name + "_varying_" + dir_name + "/" + bg_label + "/";
  std::cout << "=> Set Folder Name: " << folder_name << std::endl;

  std::vector<Variable*> variables = super_vars->get_signal_cut_vars();
  ClassifierOutputs::plot_classifiers_for_all_files(files, method_name, folder_name, bg_chains[bg_to_train]->label);
  RocCurves::get_rocs(files, signal_chain, bg_chains[bg_to_train], super_vars, method_name, folder_name);

for (int i = 0; i < file_paths.size(); i++)
  	{
    
    std::remove(file_paths[i]);
    //files_arr[i] = TFile::Open(file_paths[i]);
  	}

}
//__________________________________________________________________________________________________________________________________________

std::vector<TFile*> MVAAnalysis::get_files_from_paths(std::vector<const char*> file_paths)
{
  TFile* files_arr[file_paths.size()];
	 for (int i = 0; i < file_paths.size(); i++)
  	{
    
    //std::remove(file_paths[i]);
    files_arr[i] = TFile::Open(file_paths[i]);
  	}

	 std::vector<TFile*> files (files_arr, files_arr + sizeof(files_arr) / sizeof(files_arr[0]));
  std::cout << "File paths size: " << file_paths.size() << std::endl;
  std::cout << "files size: " << files.size() << std::endl;
	 return files;
}
//__________________________________________________________________________________________________________________________________________
TFile* MVAAnalysis::get_mva_results(std::vector<DataChain*> bg_chains, int bg_to_train, DataChain* signal_chain, DataChain* data_chain,
																																			SuperVars* super_vars, std::string folder_name, std::string method_name, const char* NTrees,
																																			const char* BoostType, const char* AdaBoostBeta,const char* SeparationType,const char* nCuts,
																																			const char* NeuronType, const char* NCycles, const char* HiddenLayers, bool unique_output_files,
																																			bool create_cards, std::string job_name, std::string mva_cut)
{
	 std::vector<Variable*> vars = super_vars->get_signal_cut_vars();
	 std::vector<Variable*> vars2 = super_vars->get_discriminating_vars();
	 std::string selection_str = super_vars->get_final_cuts_str();
	 TFile* trained_output;
  const char* trained_bg_label = bg_chains[bg_to_train]->label;

	 std::string app_output_name = BDTAnalysis::BDT_output_file_path(folder_name, job_name, false,
																																																													 NTrees, BoostType, AdaBoostBeta, SeparationType, nCuts,
																																																														trained_bg_label);


  if (method_name == "BDT")
	 {
  		trained_output = BDTAnalysis::create_BDT(bg_chains[bg_to_train], signal_chain, &vars2, folder_name,
																																													NTrees,BoostType,AdaBoostBeta, SeparationType, nCuts, job_name);
	 }
  else if (method_name == "MLP")
  {
  		trained_output = MLPAnalysis::create_MLP(bg_chains[bg_to_train], signal_chain, &vars2, folder_name,
																																													NeuronType, NCycles, HiddenLayers);
  }
  std::cout << "=> Trained method " << method_name << ", output file: " << trained_output->GetName() << std::endl;

  std::vector<DataChain*> output_bg_chains = get_output_bg_chains(bg_chains, vars, method_name, app_output_name, job_name,
																																																																		trained_bg_label, unique_output_files);

  std::cout << "=> All background put through BDT" << std::endl;

  DataChain* output_signal_chain = get_output_signal_chain(signal_chain, vars, method_name, app_output_name, job_name,
																																																											trained_bg_label, unique_output_files);
	 std::cout << "=> Signal put through BDT" << std::endl;

	 DataChain* output_data_chain = get_output_signal_chain(data_chain, vars, method_name, app_output_name, job_name,
																																																									trained_bg_label, unique_output_files);

   std::string a = get_app_filename_for_chain( app_output_name,  bg_chains[bg_to_train]->label, data_chain->label);

	 std::cout << "=> Data put through BDT" << std::endl;

	 Variable* mva_output = new Variable("output","MVA Output","-1.0","1.0","-0.8","0.8","125","1", "", false);

	 std::cout << "=> Declared MVA_Output Variable" << std::endl;

	 std::string output_graph_name = build_output_graph_name(trained_output, mva_cut);

	 std::cout << "mva output graph name: " << output_graph_name << std::endl;

	 HistoPlot::draw_plot(mva_output, output_bg_chains, output_signal_chain, output_data_chain, true, &vars, false,
																							output_graph_name, mva_cut);


  if (create_cards) {create_datacards(output_data_chain, output_signal_chain, output_bg_chains,
																																										mva_output, true, &vars, trained_output, method_name);}

  std::cout << "=> Drew MVA Output plot for all backgrounds and signal" << std::endl;
  std::cout << "Trained output name: "<< trained_output->GetName() << " " << trained_output << std::endl;

  return trained_output;
}
//__________________________________________________________________________________________________________________________________________

// creates datacards for a variety of output values
void MVAAnalysis::create_datacards(DataChain* output_data_chain, DataChain* output_signal_chain, std::vector<DataChain*> output_bg_chains,
																																			Variable* mva_output, bool with_cut, std::vector<Variable*>* variables, TFile* trained_output,
																																			std::string method_name)
{
	if (method_name == "BDT")
	{
			std::string cut_arr[] = {"output>-0.8", "output>-0.7", "output>-0.6", "output>-0.5", "output>-0.4", "output>-0.3", "output>-0.2",
																												"output>-0.1", "output>0.0", "output>0.1", "output>0.2", "output>0.3", "output>0.4", "output>0.5",
																												"output>0.6","output>0.7", "output>0.8", "output>0.9"};

			for (int i = 0; i < sizeof(cut_arr)/sizeof(cut_arr[0]); i++)
			{
					std::string output_graph_name = build_output_graph_name(trained_output, cut_arr[i]);
					DataCard::create_datacard(output_data_chain, output_signal_chain, output_bg_chains,
																															mva_output, true, variables, output_graph_name, cut_arr[i]);
			}
	}
	else
	{
			std::string cut_arr[] = {"output>-1.5", "output>-1.4", "output>-1.3", "output>-1.2", "output>-1.1", "output>-1.0", "output>-0.8",
																												"output>-0.7", "output>-0.6", "output>-0.5", "output>-0.4", "output>-0.3", "output>-0.2",
																												"output>-0.1", "output>0.0", "output>0.1", "output>0.2", "output>0.3", "output>0.4", "output>0.5",
																												"output>0.6","output>0.7", "output>0.8", "output>0.9", "output>1.0", "output>1.1"};

			for (int i = 0; i < sizeof(cut_arr)/sizeof(cut_arr[0]); i++)
			{
					std::string output_graph_name = build_output_graph_name(trained_output, cut_arr[i]);
			  DataCard::create_datacard(output_data_chain, output_signal_chain, output_bg_chains,
																															mva_output, true, variables, output_graph_name, cut_arr[i]);
			}
	}
}
//__________________________________________________________________________________________________________________________________________

// gets the name for the reader output Tfile
std::string MVAAnalysis::get_app_filename_for_chain(std::string app_output_name, const char* trained_bg_label, const char* app_label)
{
	 std::string app_label_str = app_label;
//std::cout<<"original file: "<<app_output_name<<"\n";
	 std::string trained_bg_label_str = trained_bg_label;
  std::string analysis =  app_output_name.substr(0,app_output_name.find("/")+1); //  analysis
  std::string m1 =  app_output_name.substr(app_output_name.find("/")+1,-1); //gives varying/background/file name
//std::cout<<"analysis: "<<analysis<<"\n";
//std::cout<<"m1: "<<m1<<"\n";
  std::string varying_param =  m1.substr(0,m1.find("/")+1); // gives varying param/
//std::cout<<"varying param: "<<varying_param<<"\n";

  std::string m3 = m1.substr(m1.find("/")+1,-1); //gives background/file name
//std::cout<<"m3: "<<m3<<"\n";
  std::string background = m3.substr(0,m3.find("/")+1); // give background/
//std::cout<<"background: "<<background<<"\n";
  std::string file = m3.substr(m3.find("/")+1,-1); //file name
//std::cout<<"file: "<<file<<"\n";
  file = HistoPlot::replace_all(file, trained_bg_label_str, app_label_str);
  std::string file_path = analysis + varying_param + background + file;
	 return file_path;
}

std::vector<DataChain*> MVAAnalysis::get_output_bg_chains(std::vector<DataChain*> bg_chains, std::vector<Variable*> vars,
																																																										std::string method_name, std::string app_output_name, std::string job_name,
																																																										const char* trained_bg_label, bool unique_output_files)
{
  std::vector<DataChain*> output_bg_chains;

	 for (int i = 0; i < bg_chains.size(); i++)
  {
	 		DataChain* combined_output;
	 		std::string real_app_output_name = get_app_filename_for_chain(app_output_name, trained_bg_label, bg_chains[i]->label);

	 		if (method_name == "BDT")
	 		{
	 		  combined_output = BDTAnalysis::get_BDT_results(bg_chains[i], &vars, real_app_output_name, job_name,
																																																					unique_output_files);
	 		}
	 		else if (method_name == "MLP")
	 		{
	 				//needs fixing combined_output = MLPAnalysis::get_MLP_results(bg_chains[i], &vars, app_output_name, job_name);
	 		}

	 		output_bg_chains.push_back(combined_output);
    if(std::ifstream(real_app_output_name.c_str())){remove(real_app_output_name.c_str());} 

  }

	 return output_bg_chains;
}
//__________________________________________________________________________________________________________________________________________

DataChain* MVAAnalysis::get_output_signal_chain(DataChain* signal_chain, std::vector<Variable*> vars, std::string method_name,
																																																std::string app_output_name, std::string job_name, const char* trained_bg_label,
																																																bool unique_output_files)
{
	 std::string real_app_output_name = get_app_filename_for_chain(app_output_name, trained_bg_label, signal_chain->label);
  DataChain* combined_output;
	 if (method_name == "BDT")
		{
		  combined_output =  BDTAnalysis::get_BDT_results(signal_chain, &vars, real_app_output_name, job_name, unique_output_files);
    if(std::ifstream(real_app_output_name.c_str())){remove(real_app_output_name.c_str());} 
    return combined_output;
		}
		else
		{
		 	// needs fixing: return MLPAnalysis::get_MLP_results(signal_chain, &vars, app_output_name, job_name);
		}
}

std::string MVAAnalysis::build_output_graph_name(TFile* trained_output, std::string mva_cut)
{
  std::string file_name = trained_output->GetName();

  if (mva_cut != "")
  {
  		return HistoPlot::replace_all(file_name, ".root", mva_cut + ".png");
  }
  else
  {
  		return HistoPlot::replace_all(file_name, ".root", "output_nocuts.png");
  }
}

//__________________________________________________________________________________________________________________________________________

std::vector<const char*> MVAAnalysis::vary_parameters(std::vector<DataChain*> bg_chains, int bg_to_train, DataChain* signal_chain,
																																																						DataChain* data_chain, SuperVars* super_vars, std::string method_name,
																																																						std::string dir_name, std::vector<const char*> NTrees, std::vector<const char*> BoostType,
																				                                  std::vector<const char*> AdaBoostBeta, std::vector<const char*> SeparationType,
																																																						std::vector<const char*> nCuts, std::vector<const char*> NeuronType,
																																																						std::vector<const char*> NCycles, std::vector<const char*> HiddenLayers,
																																																						bool unique_output_files, bool create_cards, std::string job_name)
{
	 std::string bg_label = bg_chains[bg_to_train]->label;
	 std::string folder_name = "analysis/" + method_name + "_varying_" + dir_name + "/" + bg_label;

	 if (method_name == "BDT")
	 {
    if (dir_name == "NTrees")
    	{
    		const char* files_arr[NTrees.size()];
    		for (int i = 0; i < NTrees.size(); i++)
      	{
    				TFile* file = get_mva_results(bg_chains, bg_to_train, signal_chain, data_chain, super_vars, folder_name,
																																						method_name, NTrees[i], BoostType[0], AdaBoostBeta[0], SeparationType[0], nCuts[0], NeuronType[0],
																																						NCycles[0], HiddenLayers[0], unique_output_files, create_cards, job_name, "");
    				const char* file_path = file->GetName();
    				files_arr[i] = file_path;

      	}
    		 std::vector<const char*> files (files_arr, files_arr + sizeof(files_arr) / sizeof(files_arr[0]));

    		 return files;
    	}
    else if (dir_name == "BoostType")
    	{
    		const char* files_arr[BoostType.size()];
    		for (int i = 0; i < BoostType.size(); i++)
      	{
    				TFile* file = get_mva_results(bg_chains, bg_to_train, signal_chain, data_chain, super_vars,
																																						folder_name, method_name, NTrees[0],
																																						BoostType[i], AdaBoostBeta[0], SeparationType[0], nCuts[0], NeuronType[0],
																																						NCycles[0], HiddenLayers[0], unique_output_files, create_cards, job_name, "");
    				const char* file_path = file->GetName();
    				files_arr[i] = file_path;
      	}
    		std::vector<const char*> files (files_arr, files_arr + sizeof(files_arr) / sizeof(files_arr[0]));

    		return files;
    	}
    else if (dir_name == "AdaBoostBeta")
    	{
    		const char* files_arr[AdaBoostBeta.size()];
    		for (int i = 0; i < AdaBoostBeta.size(); i++)
      	{
    				TFile* file = get_mva_results(bg_chains, bg_to_train, signal_chain, data_chain, super_vars,
																																						folder_name, method_name, NTrees[0],
																																						BoostType[0], AdaBoostBeta[i], SeparationType[0], nCuts[0], NeuronType[0],
																																						NCycles[0], HiddenLayers[0], unique_output_files, create_cards, job_name, "");
    				const char* file_path = file->GetName();
    				files_arr[i] = file_path;
      	}
    		std::vector<const char*> files (files_arr, files_arr + sizeof(files_arr) / sizeof(files_arr[0]));

    		return files;
    	}
    else if (dir_name == "SeparationType")
    	{
    		const char* files_arr[SeparationType.size()];
    		for (int i = 0; i < SeparationType.size(); i++)
      	{
    				TFile* file = get_mva_results(bg_chains, bg_to_train, signal_chain, data_chain, super_vars,
																																						folder_name, method_name, NTrees[0],
																																						BoostType[0], AdaBoostBeta[0], SeparationType[i], nCuts[0], NeuronType[0],
																																						NCycles[0], HiddenLayers[0], unique_output_files, create_cards, job_name, "");
    				const char* file_path = file->GetName();
    				files_arr[i] = file_path;
      	}
    		std::vector<const char*> files (files_arr, files_arr + sizeof(files_arr) / sizeof(files_arr[0]));

    		return files;
    	}
    else if (dir_name == "nCuts")
    	{
    		const char* files_arr[nCuts.size()];
    		for (int i = 0; i < nCuts.size(); i++)
      	{
    				TFile* file = get_mva_results(bg_chains, bg_to_train, signal_chain, data_chain, super_vars,
																																						folder_name, method_name, NTrees[0],
																																						BoostType[0], AdaBoostBeta[0], SeparationType[0], nCuts[i],
																																						NeuronType[0], NCycles[0], HiddenLayers[0], unique_output_files,
																																						create_cards, job_name, "");
    				const char* file_path = file->GetName();
    				files_arr[i] = file_path;
      	}
    		std::vector<const char*> files (files_arr, files_arr + sizeof(files_arr) / sizeof(files_arr[0]));

    		return files;
    	}
	 }
	 else if (method_name == "MLP")
	 {
	 		if (dir_name == "NeuronType")
	 			{
	 				const char* files_arr[NeuronType.size()];

    		for (int i = 0; i < NeuronType.size(); i++)
      	{
    				TFile* file = get_mva_results(bg_chains, bg_to_train, signal_chain, data_chain, super_vars,
																																						folder_name, method_name, NTrees[0], BoostType[0],
																																						AdaBoostBeta[0], SeparationType[0], nCuts[0], NeuronType[i], NCycles[0],
																																						HiddenLayers[0], unique_output_files, create_cards, job_name, "");
    				const char* file_path = file->GetName();
    				std::cout << file_path << std::endl;
    				files_arr[i] = file_path;
    				std::cout << files_arr[i] << std::endl;
      	}

    		std::vector<const char*> files (files_arr, files_arr + sizeof(files_arr) / sizeof(files_arr[0]));

    		return files;
	 			}
	 		else if (dir_name == "NCycles")
	 			{
	 				const char* files_arr[NCycles.size()];
    		for (int i = 0; i < NCycles.size(); i++)
      	{
    				TFile* file = get_mva_results(bg_chains, bg_to_train, signal_chain, data_chain, super_vars,
																																						folder_name, method_name, NTrees[0], BoostType[0],
																																						AdaBoostBeta[0], SeparationType[0], nCuts[0], NeuronType[0],
																																						NCycles[i], HiddenLayers[0], unique_output_files, create_cards, job_name, "");
    				const char* file_path = file->GetName();
    				files_arr[i] = file_path;
      	}
    		std::vector<const char*> files (files_arr, files_arr + sizeof(files_arr) / sizeof(files_arr[0]));

    		return files;
	 			}
	 		else if (dir_name == "HiddenLayers")
	 			{
	 				const char* files_arr[HiddenLayers.size()];
    		for (int i = 0; i < HiddenLayers.size(); i++)
      	{
    				TFile* file = get_mva_results(bg_chains, bg_to_train, signal_chain, data_chain, super_vars,
																																						folder_name, method_name, NTrees[0],
																																						BoostType[0], AdaBoostBeta[0], SeparationType[0], nCuts[0],
																																						NeuronType[0], NCycles[0], HiddenLayers[i], unique_output_files,
																																						create_cards, job_name, "");
    				const char* file_path = file->GetName();
    				files_arr[i] = file_path;
      	}
    		std::vector<const char*> files (files_arr, files_arr + sizeof(files_arr) / sizeof(files_arr[0]));

    		return files;
	 			}
	 }
}

//_________________________________________________________________________________________________________
//__________________ Unused methods from old TMVA classification attempt with categories __________________

TH1F* MVAAnalysis::plot_output(DataChain* combined_data)
{
  combined_data->chain->Draw("output>>output_hist(100, -0.8, 0.8)");
  TH1F* histo = (TH1F*)gDirectory->Get("output_hist");

  return histo;
}

std::vector<double> MVAAnalysis::get_categories(TH1F* output_histo)
{
  int first_bin = output_histo->FindFirstBinAbove(0.0, 1);
  int zero_bin = output_histo->FindBin(0.0);
  int last_bin = output_histo->FindLastBinAbove(0.0, 1);
  double total_integral = output_histo->Integral(0, output_histo->GetNbinsX() + 1);

  double tmp_integral_bg;
  double tmp_integral_sig;
  int bg_bin;
  int sig_bin;
  //bg
  for (int i = zero_bin; i > first_bin; i--)
  	{
    tmp_integral_bg = output_histo->Integral(first_bin, i);
    bg_bin = i;
    if ((tmp_integral_bg / total_integral) < 0.1)
    	{
    		break;
    	}
  	}
  std::cout << "bg" << bg_bin << std::endl;
  //sig
  for (int i = zero_bin; i < last_bin; i++)
  {
    tmp_integral_sig = output_histo->Integral(i, last_bin);
    sig_bin = i;
    if ((tmp_integral_bg / total_integral) < 0.1)
    {
      break;
    }
  }
  std::cout << "sig" << sig_bin << std::endl;
	 double bg_up_lim = output_histo->GetBinCenter(bg_bin);
  double sig_low_lim = output_histo->GetBinCenter(sig_bin);
  double x_arr[] = {-0.8, bg_up_lim, 0, sig_low_lim, 0.8};
  std::cout << "bg_xmax" << bg_up_lim << std::endl;
  std::cout << "sig_xmin" << sig_low_lim << std::endl;
  std::vector<double> categories (x_arr, x_arr + sizeof(x_arr) / sizeof(x_arr[0]));

  return categories;
}

std::vector<std::string> MVAAnalysis::get_category_strs(std::vector<double> categories)
{
  std::string cat_strs[4];

	 for (int i = 0; i < 4; i++)
  {
    cat_strs[i] = "(output > " + HistoPlot::get_string_from_double(categories[i]) + ")&&(output <= " + HistoPlot::get_string_from_double(categories[i+1]) + ")";
  }
	 std::vector<std::string> cat_strs_vector (cat_strs, cat_strs + sizeof(cat_strs) / sizeof(cat_strs[0]));
	 std::cout << "cat_strs_vector done" << cat_strs_vector[0] << std::endl;
	 return cat_strs_vector;
}

TH1F* MVAAnalysis::build_histo(DataChain* combined_output, std::string selection_str, Variable* variable, std::string histo_label)
{
  std::string hist_label = combined_output->extra_label + histo_label;

  std::string var_arg = variable->build_var_string(hist_label.c_str(), true);

  combined_output->chain->Draw(var_arg.c_str(), selection_str.c_str(), "goff");

  TH1F* histo = (TH1F*)gDirectory->Get(hist_label.c_str());

  return histo;
}

std::string MVAAnalysis::build_output_sel_str(std::string category, std::string final_cuts)
{
	 std::string selection_str = final_cuts;
	 selection_str.insert(selection_str.find("(") + 1, category + "&&");

	 return selection_str;
}

TH1F* MVAAnalysis::draw_signal(DataChain* combined_output, std::string category, std::string final_cuts, Variable* variable)
{
  combined_output->chain->SetLineColor(2);
  combined_output->chain->SetLineWidth(3);
  combined_output->chain->SetFillColor(0);

  TH1F* histo = build_histo(combined_output, build_output_sel_str(category, final_cuts), variable, "signal");

  return histo;
}

TH1F* MVAAnalysis::draw_background(DataChain* combined_output, std::string category, std::string final_cuts, Variable* variable)
{
  combined_output->chain->SetLineColor(1);
  combined_output->chain->SetFillColor(40);

  TH1F* hist = build_histo(combined_output, build_output_sel_str(category, final_cuts), variable, "bg");

  return hist;
}

void MVAAnalysis::style_histo(TH1F* histo)
{
	 histo->GetYaxis()->SetTitle("Events");
	 histo->GetYaxis()->SetLabelSize(0.035);
	 histo->GetYaxis()->SetTitleOffset(1.55);
	 histo->GetXaxis()->SetLabelSize(0);
	 histo->SetStats(kFALSE);
	 histo->SetTitle("");
}

void MVAAnalysis::draw_histo(DataChain* combined_output, std::string final_cuts, Variable* variable)
{
	 std::vector<std::string> category_strs = get_category_strs(get_categories(plot_output(combined_output)));

  TLegend* legend = new TLegend(0.0, 0.5, 0.0, 0.88);
  TCanvas* c1     = new TCanvas("c1", variable->name_styled, 800, 800);
  TPad* p1        = new TPad("p1", "p1", 0.0, 0.95, 1.0, 1.0);
  TPad* p2        = new TPad("p2", "p2", 0.0, 0.2, 1.0, 0.95);
  TPad* p3        = new TPad("p3", "p3", 0.0, 0.0, 1.0, 0.2);

  p1->SetLeftMargin(0.102);
  p2->SetBottomMargin(0.012);
  p2->SetLeftMargin(0.105);
  p3->SetBottomMargin(0.3);
  p3->SetLeftMargin(0.102);
  p1->Draw();
  p2->Draw();
  p3->Draw();
  p2->cd();

  TH1F* very_sig = draw_signal(combined_output, category_strs[3], final_cuts, variable);
  TH1F* very_bg  = draw_background(combined_output, category_strs[0], final_cuts, variable);

  legend->AddEntry(very_sig, "Signal (Top 10% of output)", "f");
  legend->AddEntry(very_bg, "Background (Bottom 10% of output)", "f");
  HistoPlot::style_legend(legend);

  very_bg->Draw();
  very_sig->Draw("SAME");
  legend->Draw("SAME");
  style_histo(very_bg);
  style_histo(very_sig);
  TH1F* plot_histos[2] = {very_sig, very_bg};
  std::vector<TH1F*> plot_histos_vector (plot_histos, plot_histos + sizeof(plot_histos) / sizeof(plot_histos[0]));
  TH1F* max_histo = HistoPlot::get_max_histo(plot_histos_vector);

  very_bg->SetMaximum(HistoPlot::get_histo_y_max(max_histo)*1.1);

  HistoPlot::draw_subtitle(variable, NULL, true, combined_output, final_cuts);

  p3->cd();
  TH1F* data_bg_ratio_histo = HistoPlot::data_to_bg_ratio_histo(plot_histos[0], plot_histos[1]);
  data_bg_ratio_histo->Draw("e1");
  HistoPlot::style_ratio_histo(data_bg_ratio_histo, variable->name_styled);
  HistoPlot::draw_yline_on_plot(variable, true, 1.0);
  p1->cd();
  HistoPlot::draw_title((combined_output->extra_label).c_str());
  c1->SaveAs("output_test.png");
  c1->Close();
}
