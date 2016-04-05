#include "../include/mva_analysis.h"
#include "../include/bdt_analysis.h"
#include "../include/mlp_analysis.h"
#include <sstream>

using namespace std;
void MVAAnalysis::get_plots_varying_params(std::vector<DataChain*> bg_chains, int bg_to_train, DataChain* signal_chain, DataChain* data_chain,
	SuperVars* super_vars, std::string method_name, std::string dir_name,
	std::vector<const char*> NTrees, std::vector<const char*> BoostType,
	std::vector<const char*> AdaBoostBeta, std::vector<const char*> SeparationType,
	std::vector<const char*> nCuts, std::vector<const char*> NeuronType, std::vector<const char*> NCycles,
	std::vector<const char*> HiddenLayers, std::vector<const char*> LearningRate,
	bool unique_output_files, bool create_cards, std::string job_name)
{
	/*std::vector<const char*> file_paths = vary_parameters(bg_chains, bg_to_train, signal_chain, data_chain, super_vars, method_name, dir_name,
		NTrees, BoostType, AdaBoostBeta, SeparationType, nCuts, NeuronType, NCycles, HiddenLayers,
		LearningRate, unique_output_files, create_cards, job_name);*/

		const char* file_arr[] ={

			"test/MLP-bg_zjets_vv-NeuronType=radial-NCycles=2000-HiddenLayers=2,4,8,16,32-LearningRate=0.0000001-EstimatorType=CE.root",
			"test/MLP-bg_zjets_vv-NeuronType=radial-NCycles=2000-HiddenLayers=2,4,8,16,32-LearningRate=0.0000002-EstimatorType=CE.root",
			"test/MLP-bg_zjets_vv-NeuronType=radial-NCycles=2000-HiddenLayers=2,4,8,16,32-LearningRate=0.00000007-EstimatorType=CE.root",
			"test/MLP-bg_zjets_vv-NeuronType=radial-NCycles=2000-HiddenLayers=2,4,8,16,32-LearningRate=0.00000008-EstimatorType=CE.root"
		};
		std::vector<const char*> file_paths (file_arr, file_arr +sizeof(file_arr)/sizeof(const char*));
		std::vector<TFile*> files = get_files_from_paths(file_paths);
		std::string bg_label = bg_chains[bg_to_train]->label;
		std::string folder_name = "analysis/" + method_name + "_varying_" + dir_name + "/" + bg_label + "/";
		std::cout << "=> Set Folder Name: " << folder_name << std::endl;

		std::vector<Variable*> variables = super_vars->get_signal_cut_vars();
  //ClassifierOutputs::plot_classifiers_for_all_files(files, method_name, folder_name, bg_chains[bg_to_train]->label);
		RocCurves::get_rocs(files, signal_chain, bg_chains[bg_to_train], super_vars, method_name, folder_name);
	}
//_________________________________________________________________________________________________________________________________________________

	std::vector<TFile*> MVAAnalysis::get_files_from_paths(std::vector<const char*> file_paths)
	{
		TFile* files_arr[file_paths.size()];
		for (int i = 0; i < file_paths.size(); i++)
		{
			files_arr[i] = TFile::Open(file_paths[i]);
		}

		std::vector<TFile*> files (files_arr, files_arr + sizeof(files_arr) / sizeof(files_arr[0]));
		std::cout << "File paths size: " << file_paths.size() << std::endl;
		std::cout << "files size: " << files.size() << std::endl;
		return files;
	}
//_________________________________________________________________________________________________________________________________________________
	TFile* MVAAnalysis::get_mva_results(std::vector<DataChain*> bg_chains, int bg_to_train, DataChain* signal_chain, DataChain* data_chain,
		SuperVars* super_vars, std::string folder_name, std::string method_name, const char* NTrees,
		const char* BoostType, const char* AdaBoostBeta,const char* SeparationType,const char* nCuts,
		const char* NeuronType, const char* NCycles, const char* HiddenLayers, const char* LearningRate,
		bool unique_output_files,
		bool create_cards, std::string job_name, std::string mva_cut, std::string sign, int min, int max, double digits)
	{
		std::vector<Variable*> vars = super_vars->get_signal_cut_vars();
		std::vector<Variable*> vars2 = super_vars->get_discriminating_vars();
		std::string selection_str = super_vars->get_final_cuts_str();
		TFile* trained_output;
//for(int i=0; i<8;i++){
		const char* trained_bg_label = bg_chains[bg_to_train]->label;
//step 1 get output name and train MVa
//________________________________________________________________________________________________________________________________________________
		std::string app_output_name; 
		if (method_name == "BDT")
		{
			app_output_name = BDTAnalysis::BDT_output_file_path(folder_name, job_name, false,
				NTrees, BoostType, AdaBoostBeta, SeparationType, nCuts,
				trained_bg_label);
			trained_output = BDTAnalysis::create_BDT(bg_chains[bg_to_train], signal_chain, &vars2, folder_name,
				NTrees,BoostType,AdaBoostBeta, SeparationType, nCuts, job_name);
		}
		else if (method_name == "MLP")
		{
			app_output_name = MLPAnalysis::MLP_output_file_path(folder_name, job_name, false,
				NeuronType, NCycles, HiddenLayers, LearningRate, trained_bg_label);

			trained_output = /*TFile::Open("test/MLP-all_bg-NeuronType=radial-NCycles=500-HiddenLayers=2-LearningRate=0.01-EstimatorType=CE-50bins.root");*/ MLPAnalysis::create_MLP(bg_chains[bg_to_train], signal_chain, &vars2, folder_name,
				NeuronType, NCycles, HiddenLayers, LearningRate, job_name);
		}
//MLPAnalysis::create_MLP(data_chain, signal_chain, &vars2, folder_name,
//				NeuronType, NCycles, HiddenLayers, LearningRate, job_name);
		std::cout << "=> Trained method " << method_name << ", output file: " << trained_output->GetName() << std::endl;
//}
 if (create_cards) {		


	//step 2 evaluate MVA's
//________________________________________________________________________________________________________________________________________________	

//step 2.1 get test tree 
        //step 2.1.1 get trained output file name
		const char* t_arr[] = {trained_output->GetName()};
	//step 2.1.2 intialise DataChain paramters for test set
		const char* t_label = bg_chains[bg_to_train]->label;
		const char* t_legend ="test_set";
		std::string lep_sel = bg_chains[bg_to_train]->lep_sel;
		std::vector<const char*> t_vector (t_arr, t_arr + 
			sizeof(t_arr)/sizeof(const char*));
        //step 2.1.3 open trained_output TFile
		TFile* f = new TFile(trained_output->GetName(),"update");
        //step 2.1.4 copy test tree(called lighttree) from trained_output tFile
		TTree* tree = (TTree*)f->Get("TestTree");
	tree->SetName("LightTree"); //renames tree for reader input
	tree->Write();

        //step 2.1.5 create DataChain of the test set tree
	DataChain* test_chain     = new DataChain(t_vector,t_label,t_legend, lep_sel, "");
	f->Close();

//step 2.2 get data tree
	const char* d_arr[] = {"dataTrees/data_chain.root"};
	//step 2.2.2 intialise DataChain paramters for test set
	std::vector<const char*> data_vector (d_arr, d_arr + sizeof(d_arr)/sizeof(const char*));

        //step 2.2.5 create DataChain of the test set tree
	DataChain* data_ch = new DataChain(data_vector,data_label,data_legend,"");

//step 2.2 get output chains for test,data and other BGs
//////////qcd
	const char* qcd_arr[] = {"dataTrees/bg_qcd.root"};

	std::vector<const char*> qcd_vector (qcd_arr, qcd_arr + sizeof(qcd_arr)/sizeof(const char*));

	DataChain* qcd_ch = new DataChain(qcd_vector,qcd_label,qcd_legend,"");
//////////VV
	const char* VV_arr[] = {"dataTrees/bg_vv.root"};

	std::vector<const char*> VV_vector (VV_arr, VV_arr + sizeof(VV_arr)/sizeof(const char*));

	DataChain* VV_ch = new DataChain(VV_vector,vv_label,vv_legend,"");

//////////wjets->enu
	const char* ev_arr[] = {"dataTrees/bg_wjets_ev.root"};

	std::vector<const char*> ev_vector (ev_arr, ev_arr + sizeof(ev_arr)/sizeof(const char*));

	DataChain* wjets_ev_ch = new DataChain(ev_vector,wjets_ev_label,wjets_ev_legend,"(nselelectrons == 1)");

//////////wjets->munu
	const char* muv_arr[] = {"dataTrees/bg_wjets_muv.root"};

	std::vector<const char*> muv_vector (muv_arr, muv_arr + sizeof(muv_arr)/sizeof(const char*));

	DataChain* wjets_muv_ch = new DataChain(muv_vector,wjets_muv_label,wjets_muv_legend,"(nselmuons == 1)");

//////////wjets->taunu
	const char* tauv_arr[] = {"dataTrees/bg_wjets_tauv.root"};

	std::vector<const char*> tauv_vector (tauv_arr, tauv_arr + sizeof(tauv_arr)/sizeof(const char*));

	DataChain* wjets_tauv_ch = new DataChain(tauv_vector,wjets_tauv_label,wjets_tauv_legend,"(nvetomuons==0)&&(nvetoelectrons==0)&&(ntaus == 1)");

//////////zjets->nunu
	const char* zjets_vv_arr[] = {"dataTrees/bg_zjets_vv.root"};

	std::vector<const char*> zjets_vv_vector (zjets_vv_arr, zjets_vv_arr + sizeof(zjets_vv_arr)/sizeof(const char*));

	DataChain* zjets_vv_ch = new DataChain(zjets_vv_vector,zjets_vv_label,zjets_vv_legend,"");

//////////top
	const char* top_arr[] = {"dataTrees/bg_top.root"};

	std::vector<const char*> top_vector (top_arr, top_arr + sizeof(top_arr)/sizeof(const char*));

	DataChain* top_ch = new DataChain(top_vector,top_label,top_legend,"");

//////////Z->ll
	const char* zll_arr[] = {"dataTrees/bg_zll.root"};

	std::vector<const char*> zll_vector (zll_arr, zll_arr + sizeof(zll_arr)/sizeof(const char*));

	DataChain* zll_ch = new DataChain(zll_vector,z_ll_label,z_ll_legend,"(nselmuons == 2)&&(m_mumu>60)&&(m_mumu<120)");

	DataChain* bg_ch_arr[] = {zll_ch, wjets_ev_ch, wjets_muv_ch, wjets_tauv_ch, top_ch,VV_ch, zjets_vv_ch, qcd_ch };

	std::vector<DataChain*> bg_chs (bg_ch_arr, bg_ch_arr + sizeof(bg_ch_arr) / sizeof(bg_ch_arr[0]));


//_________________________

	DataChain* mva_output_test_chain = evaluate_test_data(test_chain,  vars,  method_name,
		app_output_name, job_name, trained_bg_label, unique_output_files);

	std::cout << "=> test_tree put through MVA" << std::endl;

///other backgrounds and data passed through tree for MC weight calculation

	std::vector<DataChain*> output_bg_chains = get_output_bg_chains(bg_chs, vars, method_name, app_output_name, job_name, trained_bg_label,
																																																																	unique_output_files);

	std::cout << "=> All background put through MVA" << std::endl;



	DataChain* output_data_chain = get_output_signal_chain(data_ch, vars, method_name, app_output_name, job_name,                                 		trained_bg_label, unique_output_files);

	std::cout << "=> Data put through MVA" << std::endl;

//step 3 initalise output variable
//________________________________________________________________________________________________________________

	Variable* mva_output;

	if (method_name == "BDT")
	{
		mva_output = new Variable("output","MVA Output","-1.0","1.0","-0.8","0.8","125","1", "", false);
	}
	else if (method_name == "MLP")
	{
		mva_output = new Variable("output","MVA Output","0.0","1.0","","","50","1", "", false);
	}
	std::cout << "=> Declared MVA_Output Variable" << std::endl;

	std::string output_graph_name = build_output_graph_name(trained_output, mva_cut);

	std::cout << "mva output graph name: " << output_graph_name << std::endl;
//step 4 draw plot
//_________________________
cout<<"step 4 in mva analysis"<<endl;
HistoPlot::plot_evaluated_zjets_vv_testTree(bg_to_train, mva_output, mva_output_test_chain,output_data_chain, output_bg_chains,&vars, output_graph_name, mva_cut);


//output_bg_chains[1]->chain->Draw("output>>test(100,-1.25,1.5)", "((output>0.1)&&(classID==0)&&(nselelectrons == 1))*total_weight_lepveto");
//step 5 create datacards
//create array of test file bg and all other bgs remebering to halve the other mc weights later...



	if (create_cards)
		{
			create_datacards(bg_to_train, output_data_chain,mva_output_test_chain , output_bg_chains,
							mva_output, true, NULL, trained_output, method_name, sign, min, max, digits);
			
		}

//roc curve plots plot train set so useless
//std::cout<<"Area unde ROC: "<<RocCurves::get_auc( method_name, trained_output->GetName());



//step 6 get estimators
 //if (method_name == "MLP"){get_estimators(trained_output->GetName());}
  std::cout << "=> Drew MVA Output plot for all backgrounds and signal" << std::endl;
  std::cout << "Trained output name: "<< trained_output->GetName() << " " << trained_output << std::endl;
  //get_estimators(trained_output->GetName());

}
		return trained_output;
	}

//________________________________________________________________________________________
void MVAAnalysis::get_estimators(std::vector<const char*> training_file_paths)
{
  TFile *_file0;
//get get params from file name
//"test/MLP-bg_wjets_ev-NeuronType=radial-NCycles=2000-HiddenLayers=2,4,8,16-LearningRate=0.00001-EstimatorType=CE.root";
  string file_path;
  string background;
  string NeuronType;
  string NCycles;
  string HiddenLayers;
  string LearningRate;
  int len = training_file_paths.size();
  Double_t test_estimator[len];
  Double_t train_estimator[len];
 
//get arr of estimators;
  for(int i=0; i< len;i++)
  {
    _file0 = TFile::Open(training_file_paths[i]);
    //go into estimator hiustogram directory
    _file0->cd("Method_MLP;1/MLP;1");
    //declare histogram pointers
    TH1F* test;
    TH1F * train;
    //get histograms from subdirectory of file
    gDirectory->GetObject("estimatorHistTest;1",test);
    gDirectory->GetObject("estimatorHistTrain;1",train);
    test_estimator[i] = test->GetBinContent(100);
    train_estimator[i] = train->GetBinContent(100);
    _file0->Close();
  cout<<"test_estimator= "<<test_estimator[i]<<" train_estimator= "
  <<train_estimator[i]<<" train/test= "<<train_estimator[i]/test_estimator[i]<<"\n";

  }
  int bglen;
  int bgpos;
  int NTpos;
  int NTlen;
  int NCpos;
  int NClen;
  int HLpos;
  int HLlen;
  int LRpos;
  int LRlen;
  bool human_readable=false;
  std::fstream fs;
  fs.open ("Estimator_statistics.txt", std::fstream::out | std::fstream::trunc);
  fs << "Estimator statistics for last epoch of MLP training with G,P,N transform\n" ;
  fs << "Background NeuronType  NCycles  HiddenLayers   LearningRate  test_estimator train_estimator train/test \n";
  //done title lines 
  for(int i=0; i< len;i++)
  {
    file_path= training_file_paths[i];

    bgpos = file_path.find("bg");
    bglen = file_path.find("-NeuronType") - (file_path.find("bg"));
    background = file_path.substr(bgpos, bglen);

    NTpos = file_path.find("NeuronType=")+11;
    NTlen = file_path.find("-NCycles") - (11+file_path.find("NeuronType="));
    NeuronType = file_path.substr(NTpos, NTlen);

    NCpos = file_path.find("NCycles=")+8;
    NClen = file_path.find("-HiddenLayers") - (8+file_path.find("NCycles="));
    NCycles = file_path.substr(NCpos, NClen);

    HLpos = file_path.find("HiddenLayers=")+13;
    HLlen = file_path.find("-LearningRate") - (13+file_path.find("HiddenLayers="));
    HiddenLayers = file_path.substr(HLpos, HLlen);

    LRpos = file_path.find("LearningRate=")+13;
    LRlen = file_path.find("-EstimatorType") - (13+file_path.find("LearningRate="));
    LearningRate= file_path.substr(LRpos, LRlen);
    if(human_readable){
      fs<<background;
      for(int i =15-bglen; i>0;i--){fs<<" ";}
      fs << NeuronType <<"      "<<NCycles<<"     "<<HiddenLayers;
      for(int i =15-HLlen; i>0;i--){fs<<" ";}
      fs<<LearningRate;
      for(int i =15-LRlen; i>0;i--){fs<<" ";}
      fs<<test_estimator[i]<<"        "<<
         train_estimator[i]<<"        "<<
        train_estimator[i]/test_estimator[i]<<"\n";
    }
    else if (!human_readable){fs<<background<<" "<<NeuronType <<" "<<NCycles<<" "<<HiddenLayers<<" "<<LearningRate<<" ";
    fs<<test_estimator[i]<<" "<<train_estimator[i]<<" "<<train_estimator[i]/test_estimator[i]<<"\n";}
  }
  fs.close();

}

//________________________________________________________________________________________
std::string MVAAnalysis::create_auc_line_MLP(const char* bg_label, const char* NeuronType,
		const char* NCycles, const char* HiddenLayers,
		const char* LearningRate, double auc)
	{
		std::string line = bg_label;
		line += ",";
		line.append(NeuronType);
		line += ",";
		line.append(NCycles);
		line += ",";
		std::string hidden_layers = HiddenLayers;
		std::string hidden_layers_str = HistoPlot::replace_all(hidden_layers, ",", " ");
		line.append(hidden_layers_str);
		line += ",";
		line.append(LearningRate);
		line += ",";
		line.append(DataCard::double_to_str(auc));

		return line;
	}

//_________________________________________________________________________________________________________________
	std::vector<std::string> MVAAnalysis::get_mva_cut_range(std::string sign, int min, int max, double digits)
	{
		int arr_length = max - min + 1;
		const char* cut_arr[2] = {"2"};
		std::vector<std::string> cuts(cut_arr, cut_arr + 1);
		int counter = 0;
		cuts.erase(cuts.begin());

		for (int i = min; i < (max + 1); i++)
		{
			std::string cut = "output" + sign;
			double cut_val = double(i) / digits;
			std::string cut_val_str = DataCard::double_to_str(cut_val);
			cut += cut_val_str;
			cuts.push_back(cut);
			counter += 1;
		}

		return cuts;
	}
//________________________________________________________________________________________________________________________________________________
// creates datacards for a variety of output values
	void MVAAnalysis::create_datacards(int bg_to_train, DataChain* output_data_chain, DataChain* output_signal_chain, std::vector<DataChain*> output_bg_chains,
		Variable* mva_output, bool with_cut, std::vector<Variable*>* variables, TFile* trained_output,
		std::string method_name, std::string sign, int min, int max, double digits)
	{
		std::string trained_output_str = trained_output->GetName();
		std::string folder_name = trained_output_str;

		if (sign == ">") {folder_name += "greater_than";}
		else {folder_name += "less_than";}

		folder_name = HistoPlot::replace_all(folder_name, ".root", "");

		if (!opendir(folder_name.c_str()))
			  {
			    mkdir(folder_name.c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
			  }
	    int output_name_idx = trained_output_str.rfind("/");
		folder_name += "/" + trained_output_str.substr(output_name_idx + 1, -1);
		
		std::vector<std::string> cut_arr = get_mva_cut_range(sign, min, max, digits);

		for (int i = 0; i < cut_arr.size(); i++)
		{
			std::string output_graph_name = HistoPlot::replace_all(folder_name, ".root", cut_arr[i] + ".png");
			//build_output_graph_name(trained_output, cut_arr[i]);

			DataCard::create_datacard(bg_to_train, output_data_chain, output_signal_chain, output_bg_chains,
				mva_output, true, variables, output_graph_name, cut_arr[i]);
		}

	}

//________________________________________________________________________________________________________________________________________________
// gets the name for the reader output Tfile
	std::string MVAAnalysis::get_app_filename_for_chain(std::string app_output_name, const char* trained_bg_label, const char* app_label)
	{
		std::string app_label_str = app_label;
		std::string trained_bg_label_str = trained_bg_label;
		int last_trained_idx = app_output_name.rfind(trained_bg_label_str);

		return app_output_name.replace(last_trained_idx, trained_bg_label_str.length(), app_label_str);
	}
//________________________________________________________________________________________________________________________________________________

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
				combined_output = MLPAnalysis::get_MLP_results(bg_chains[i], &vars, real_app_output_name, job_name,
					unique_output_files);
			}

			output_bg_chains.push_back(combined_output);
		}

		return output_bg_chains;
	}
//________________________________________________________________________________________________________________________________________________

	DataChain* MVAAnalysis::evaluate_test_data(DataChain* test_chain, std::vector<Variable*> vars, std::string method_name,
		std::string app_output_name, std::string job_name, const char* trained_bg_label,
		bool unique_output_files)
	{
		std::string real_app_output_name = get_app_filename_for_chain(app_output_name, trained_bg_label, test_chain->label);

		if (method_name == "BDT")
		{
			return BDTAnalysis::get_BDT_results(test_chain, &vars, real_app_output_name, job_name, unique_output_files);
		}
		else
		{
			return MLPAnalysis::get_MLP_results(test_chain, &vars, real_app_output_name, job_name, unique_output_files);
		}
	}


//________________________________________________________________________________________________________________________________________________
	DataChain* MVAAnalysis::get_output_signal_chain(DataChain* signal_chain, std::vector<Variable*> vars, std::string method_name,
		std::string app_output_name, std::string job_name, const char* trained_bg_label,
		bool unique_output_files)
	{
		std::string real_app_output_name = get_app_filename_for_chain(app_output_name, trained_bg_label, signal_chain->label);

		if (method_name == "BDT")
		{
			return BDTAnalysis::get_BDT_results(signal_chain, &vars, real_app_output_name, job_name, unique_output_files);
		}
		else
		{
			return MLPAnalysis::get_MLP_results(signal_chain, &vars, real_app_output_name, job_name, unique_output_files);
		}
	}
//________________________________________________________________________________________________________________________________________________

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

//________________________________________________________________________________________________________________________________________________
	std::vector<const char*> MVAAnalysis::vary_parameters(std::vector<DataChain*> bg_chains, int bg_to_train, DataChain* signal_chain,
		DataChain* data_chain, SuperVars* super_vars, std::string method_name,
		std::string dir_name, std::vector<const char*> NTrees, std::vector<const char*> BoostType,
		std::vector<const char*> AdaBoostBeta, std::vector<const char*> SeparationType,
		std::vector<const char*> nCuts, std::vector<const char*> NeuronType,
		std::vector<const char*> NCycles, std::vector<const char*> HiddenLayers,
		std::vector<const char*> LearningRate,
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
						NCycles[0], HiddenLayers[0], LearningRate[0],
						unique_output_files, create_cards, job_name, "");
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
						NCycles[0], HiddenLayers[0], LearningRate[0],
						unique_output_files, create_cards, job_name, "");
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
						NCycles[0], HiddenLayers[0], LearningRate[0],
						unique_output_files, create_cards, job_name, "");
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
						NCycles[0], HiddenLayers[0], LearningRate[0], unique_output_files,
						create_cards, job_name, "");
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
						NeuronType[0], NCycles[0], HiddenLayers[0], LearningRate[0], unique_output_files,
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
						HiddenLayers[0], LearningRate[0], unique_output_files, create_cards, job_name, "");
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
						NCycles[i], HiddenLayers[0],LearningRate[0], unique_output_files, create_cards, job_name, "");
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
						NeuronType[0], NCycles[0], HiddenLayers[i], LearningRate[0], unique_output_files,
						create_cards, job_name, "");
					const char* file_path = file->GetName();
					files_arr[i] = file_path;
				}
				std::vector<const char*> files (files_arr, files_arr + sizeof(files_arr) / sizeof(files_arr[0]));

				return files;
			}
			else if (dir_name == "LearningRate")
			{
				const char* files_arr[LearningRate.size()];
				for (int i = 0; i < LearningRate.size(); i++)
				{
					TFile* file = get_mva_results(bg_chains, bg_to_train, signal_chain, data_chain, super_vars,
						folder_name, method_name, NTrees[0],
						BoostType[0], AdaBoostBeta[0], SeparationType[0], nCuts[0],
						NeuronType[0], NCycles[0], HiddenLayers[0], LearningRate[i], unique_output_files,
						create_cards, job_name, "");
					const char* file_path = file->GetName();
					files_arr[i] = file_path;
				}
				std::vector<const char*> files (files_arr, files_arr + sizeof(files_arr) / sizeof(files_arr[0]));

				return files;
			}
		}
	}
//--------------------------------------------------------------------------------------------------------------------------------
//--------------------------------------------------------------------------------------------------------------------------------



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
