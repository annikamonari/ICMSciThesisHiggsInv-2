#include "../include/mlp_analysis.h"

TFile* MLPAnalysis::create_MLP(DataChain* bg_chain, DataChain* signal_chain, std::vector<Variable*>* variables,
																															std::string folder_name, const char* NeuronType, const char* NCycles, const char* HiddenLayers,
																															const char* LearningRate, std::string job_name)
{
	  if (!opendir(folder_name.c_str()))
	  {
	    mkdir(folder_name.c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
	  }

	  std::string output_path = MLP_output_file_path(folder_name, job_name, true, NeuronType, NCycles,
																																																		HiddenLayers, LearningRate, bg_chain->label);

	  TFile* output_tmva = TFile::Open(output_path.c_str(),"RECREATE");

   std::string weight_file_name = "TMVAClassification" + job_name;

	  TMVA::Factory* factory = new TMVA::Factory(weight_file_name.c_str(), output_tmva,
	                                             "!V:!Silent:Color:DrawProgressBar:Transformations=I;D;P;G,D:AnalysisType=Classification");

	  for (int i = 0; i < variables->size(); i++)
	  {
	    factory->AddVariable((*variables)[i]->name, (*variables)[i]->name_styled, (*variables)[i]->units, 'F');
	  }

	  // Background
	    double background_weight = 1.0;
	    factory->AddBackgroundTree(bg_chain->chain,background_weight);
	    factory->SetBackgroundWeightExpression("total_weight_lepveto");

	    // Signal
	    double signal_weight = 1.0;
	    factory->AddSignalTree(signal_chain->chain, signal_weight);
	    factory->SetSignalWeightExpression("total_weight_lepveto");

	    // Apply additional cuts on the signal and background samples (can be different)
	    TCut signal_cuts = "alljetsmetnomu_mindphi>2.0 && jet1_pt>50.0 && jet2_pt>45.0 && metnomu_significance>3.5 && dijet_deta>4.2 && dijet_deta<8.0 && nvetomuons==0 && nvetoelectrons==0"; // for example: TCut mycuts = "abs(var1)<0.5 && abs(var2-0.5)<1";
	    TCut bg_cuts = signal_cuts; // for example: TCut mycutb = "abs(var1)<0.5";

	    factory->PrepareTrainingAndTestTree(signal_cuts, bg_cuts,
	    				       "SplitMode=Random:NormMode=EqualNumEvents:!V" );

	    factory->BookMethod(TMVA::Types::kMLP, "MLP", MLP_options_str(NeuronType, NCycles, HiddenLayers, LearningRate) );

  // Train MVAs using the set of training events
  factory->TrainAllMethods();

  // ---- Evaluate all MVAs using the set of test events
  factory->TestAllMethods();

  // ----- Evaluate and compare performance of all configured MVAs
  factory->EvaluateAllMethods();

  // --------------------------------------------------------------

  // Save the output
  output_tmva->Close();

  std::cout << "==> Wrote root file: " << output_tmva->GetName() << std::endl;
  std::cout << "==> TMVAClassification is done!" << std::endl;
  std::cout << std::endl;
  std::cout << "==> To view the results, launch the GUI: \"root -l ./TMVAGui.C\"" << std::endl;
  std::cout << std::endl;

  delete factory;

  return output_tmva;
}

TTree* MLPAnalysis::evaluate_MLP(DataChain* bg_chain, std::vector<Variable*>* variables, std::string output_name,
																																	std::string job_name, bool unique_output_files)
{
	TMVA::Reader* reader = new TMVA::Reader( "!Color:!Silent" );

		   Float_t dijet_deta;
		   Float_t forward_tag_eta;
		   Float_t metnomu_significance;
		   Float_t sqrt_ht;
		   Float_t alljetsmetnomu_mindphi;
					Float_t dijet_M;
		   Float_t metnomuons;

		   Float_t jet1_pt;
		   Float_t jet2_pt;

		   Float_t jet1_eta;
		   Float_t jet2_eta;
		   Float_t jet1_phi;
		   Float_t jet2_phi;
		   Float_t jet_csv1;
		   Float_t jet_csv2;
		   Float_t dijet_dphi;
		   Float_t metnomu_x;
		   Float_t metnomu_y;
		   Float_t sumet;
		   Float_t mht;
		   Float_t unclustered_et;
		   Float_t jetmetnomu_mindphi;
		   Float_t jetunclet_mindphi;
		   Float_t metnomuunclet_dphi;
		   Float_t dijetmetnomu_vectorialSum_pt;
		   Float_t dijetmetnomu_ptfraction;
		   Float_t jet1metnomu_scalarprod;
		   Float_t jet2metnomu_scalarprod;

     reader->AddVariable("dijetmetnomu_ptfraction", &dijetmetnomu_ptfraction);
     reader->AddVariable("dijetmetnomu_vectorialSum_pt", &dijetmetnomu_vectorialSum_pt);
     reader->AddVariable("jet_csv2", &jet_csv2);
     reader->AddVariable("dijet_dphi", &dijet_dphi);
     reader->AddVariable("dijet_M", &dijet_M);

/*

		   reader->AddVariable("alljetsmetnomu_mindphi", &alljetsmetnomu_mindphi);
		   reader->AddVariable("forward_tag_eta", &forward_tag_eta);
		   reader->AddVariable("dijet_deta", &dijet_deta);
		   reader->AddVariable("metnomu_significance", &metnomu_significance);
		   //reader->AddVariable("sqrt_ht", &sqrt_ht);
		   //reader->AddVariable("dijet_M", &dijet_M);
		   reader->AddVariable("metnomuons", &metnomuons);

		   reader->AddVariable("jet1_pt", &jet1_pt);
		   reader->AddVariable("jet2_pt", &jet2_pt);

		   reader->AddVariable("jet1_eta", &jet1_eta);
		   reader->AddVariable("jet2_eta", &jet2_eta);
		   reader->AddVariable("jet1_phi", &jet1_phi);
		   reader->AddVariable("jet2_phi", &jet2_phi);
		   reader->AddVariable("jet_csv1", &jet_csv1);
		   reader->AddVariable("jet_csv2", &jet_csv2);
		   reader->AddVariable("dijet_dphi", &dijet_dphi);
		   reader->AddVariable("metnomu_x", &metnomu_x);
		   reader->AddVariable("metnomu_y", &metnomu_y);
		   reader->AddVariable("sumet", &sumet);
		   reader->AddVariable("mht", &mht);
		   reader->AddVariable("unclustered_et", &unclustered_et);
		   reader->AddVariable("jetmetnomu_mindphi", &jetmetnomu_mindphi);
		   reader->AddVariable("jetunclet_mindphi", &jetunclet_mindphi);
		   reader->AddVariable("metnomuunclet_dphi", &metnomuunclet_dphi);
		   //reader->AddVariable("dijetmetnomu_vectorialSum_pt", &dijetmetnomu_vectorialSum_pt);
		   reader->AddVariable("dijetmetnomu_ptfraction", &dijetmetnomu_ptfraction);
		   //reader->AddVariable("jet1metnomu_scalarprod", &jet1metnomu_scalarprod);
		   //reader->AddVariable("jet2metnomu_scalarprod", &jet2metnomu_scalarprod);
	*/

		   std::string weight_file_path = "weights/TMVAClassification" + job_name + "_MLP.weights.xml";

	   // Book method(s)
 	   reader->BookMVA( "MLP method", weight_file_path.c_str() );
	   // Book output histograms
	   TH1F* histNn     = new TH1F( "MVA_MLP", "MVA_MLP", 100, -1.25, 1.5 );
	   // --- Event loop

	   	   TChain* data = (TChain*) bg_chain->chain->Clone();

	   	   data->SetBranchAddress("dijet_deta", &dijet_deta);
	   	   data->SetBranchAddress("forward_tag_eta", &forward_tag_eta);
	   	   data->SetBranchAddress("metnomu_significance", &metnomu_significance);
	   	   data->SetBranchAddress("sqrt_ht", &sqrt_ht);
	   	   data->SetBranchAddress("alljetsmetnomu_mindphi", &alljetsmetnomu_mindphi);
	   	   data->SetBranchAddress("dijet_M", &dijet_M);
	   	   data->SetBranchAddress("metnomuons", &metnomuons);

	   	   data->SetBranchAddress("jet1_pt", &jet1_pt);
	   	   data->SetBranchAddress("jet2_pt", &jet2_pt);
	   	   data->SetBranchAddress("jet1_eta", &jet1_eta);
	   	   data->SetBranchAddress("jet2_eta", &jet2_eta);
	   	   data->SetBranchAddress("jet1_phi", &jet1_phi);
	   	   data->SetBranchAddress("jet2_phi", &jet2_phi);
	   	   data->SetBranchAddress("jet_csv1", &jet_csv1);
	   	   data->SetBranchAddress("jet_csv2", &jet_csv2);
	   	   data->SetBranchAddress("dijet_dphi", &dijet_dphi);
	   	   data->SetBranchAddress("metnomu_x", &metnomu_x);
	   	   data->SetBranchAddress("metnomu_y", &metnomu_y);
	   	   data->SetBranchAddress("sumet", &sumet);
	   	   data->SetBranchAddress("mht", &mht);
	   	   data->SetBranchAddress("unclustered_et", &unclustered_et);
	   	   data->SetBranchAddress("jetmetnomu_mindphi", &jetmetnomu_mindphi);
	   	   data->SetBranchAddress("jetunclet_mindphi", &jetunclet_mindphi);
	   	   data->SetBranchAddress("metnomuunclet_dphi", &metnomuunclet_dphi);
	   	   data->SetBranchAddress("dijetmetnomu_vectorialSum_pt", &dijetmetnomu_vectorialSum_pt);
	   	   data->SetBranchAddress("dijetmetnomu_ptfraction", &dijetmetnomu_ptfraction);
	   	   data->SetBranchAddress("jet1metnomu_scalarprod", &jet1metnomu_scalarprod);
	   	   data->SetBranchAddress("jet2metnomu_scalarprod", &jet2metnomu_scalarprod);

	   	   // Efficiency calculator for cut method
	   	   Int_t    nSelCutsGA = 0;
	   	   Double_t effS       = 0.7;
	   	   std::vector<Float_t> vecVar(9); // vector for EvaluateMVA tests

	   	   Float_t output;
	   	   TTree* output_tree = new TTree("MVAtree","Tree with classifier outputs");
	   	   output_tree -> Branch("output", &output, "output");

	   	   std::cout << "--- Processing: " << data->GetEntries() << " events" << std::endl;
	   	   TStopwatch sw;
	   	   sw.Start();
	   	   for (Long64_t ievt=0; ievt<data->GetEntries(); ievt++) {

	   	      if (ievt%10000 == 0) std::cout << "--- ... Processing event: " << ievt << std::endl;

	   	      data->GetEntry(ievt);

	   	      output = reader->EvaluateMVA( "MLP method" );
	   	      output_tree->Fill();
	   	      histNn->Fill(output);
	   	   }

	   	   // Get elapsed time
	   	   sw.Stop();
	   	   std::cout << "--- End of event loop: "; sw.Print();

	   	   // --- Write histograms

	       TTree* output_tree_clone = output_tree->CloneTree();

	       std::string target_name;

	       if (unique_output_files)
	       {
	       		target_name = output_name;
	       }
	       else
	       {
	       		target_name = job_name + "TMVApp.root";
	       }

	   	   TFile* target  = new TFile(target_name.c_str(),"RECREATE" );
	   	   target->cd();
	   	   output_tree_clone->Write();
	   	   histNn->Write();

	   	   target->Close();

	   	   std::cout << "--- Created root file: " << output_name << "containing the MVA output histograms" << std::endl;

	   	   delete reader;

	   	   std::cout << "==> TMVAClassificationApplication is done!" << std::endl;

	   	   return output_tree_clone;
}


//note before calling this method you must call create_bdt to update the xml weight file:
DataChain* MLPAnalysis::get_MLP_results(DataChain* bg_chain, std::vector<Variable*>* variables, std::string output_name,
																																								std::string job_name, bool unique_output_files)
{
	 TTree* output_weight = MLPAnalysis::evaluate_MLP(bg_chain,variables, output_name, job_name, unique_output_files);
	 TTree* output_weight_clone = (TTree*) output_weight->Clone();
	 output_weight_clone->SetBranchStatus("*",1);
	 TChain* bg_clone     = (TChain*) bg_chain->chain->Clone();

	 bg_clone->AddFriend(output_weight_clone);

	 std::string label(bg_chain->label);
	 label += "_w_mva_output";

  DataChain* output_data = new DataChain(z_ll, bg_chain->label, bg_chain->legend, bg_chain->lep_sel, label, bg_clone);

	 return output_data;
}

std::string MLPAnalysis::MLP_output_file_path(std::string folder_name, std::string job_name, bool is_train_file,
																																														const char* NeuronType, const char* NCycles, const char* HiddenLayers,
																																														const char* LearningRate, const char* bg_chain_label)
{
  std::string file_name = MLP_output_name_str(NeuronType, NCycles, HiddenLayers, LearningRate, bg_chain_label, job_name);

  if (!is_train_file)
  {
  		file_name.insert(file_name.find("-") + 1, "App-");
  }

  if (!opendir(folder_name.c_str()))
  {
  		mkdir(folder_name.c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
  }

  std::string output_path = folder_name + "/" + file_name;

  return output_path;
}


std::string MLPAnalysis::MLP_options_str(const char* NeuronType, const char* NCycles, const char* HiddenLayers, const char* LearningRate)
{
	std::string MLP_options = "H:!V:NeuronType=";
	std::string nt = NeuronType;
	std::string nc = NCycles;
	std::string hl = HiddenLayers;
	std::string lr = LearningRate;

	MLP_options.append(nt);
	MLP_options += ":VarTransform=N:NCycles=";
	MLP_options.append(nc);
	MLP_options += ":HiddenLayers=";
	MLP_options.append(hl);
	MLP_options += ":LearningRate=";
	MLP_options.append(lr);
	MLP_options += ":TestRate=5:UseRegulator:EstimatorType=CE";//:ConvergenceTests=10"; ConvergenceTests=100:

 return MLP_options;
}

std::string MLPAnalysis::MLP_output_name_str(const char* NeuronType, const char* NCycles, const char* HiddenLayers,
																																													const char* LearningRate, const char* bg_chain_label, std::string job_name)
{
	std::string nt = NeuronType;
	std::string nc = NCycles;
	std::string hl = HiddenLayers;
	std::string lr = LearningRate;
 std::string bg = bg_chain_label;

	std::string out_nam = "MLP-" + bg + "-NeuronType=";//
	out_nam.append(nt);
	out_nam += "-NCycles=";
	out_nam.append(nc);
	out_nam += "-HiddenLayers=";
	out_nam.append(hl);
	out_nam += "-LearningRate=";
	out_nam.append(lr);
        out_nam += "-EstimatorType=CE";
 out_nam += ".root";

	return out_nam;
}

