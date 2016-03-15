#include "../include/bdt_analysis.h"

void BDTAnalysis::add_trees_to_factory(DataChain* any_chain, bool is_bg_chain, TMVA::Factory* factory)
{
  std::vector<const char*> file_paths = any_chain->all_file_paths;
  std::string selection1 = "(alljetsmetnomu_mindphi>2.0)&&(jet1_pt>50.0)&&(jet2_pt>45.0)&&(metnomu_significance>3.5)&&";
  std::string selection2 = "(dijet_deta>4.2)&&(dijet_deta<8.0)&&(nvetomuons==0)&&(nvetoelectrons==0)";
  std::string selection = selection1 + selection2;

  for (int i = 0; i < file_paths.size(); i++)
  {
  		TFile* file = TFile::Open(file_paths[i]);
  		TTree* tree = (TTree*) file->Get("LightTree");

  		TFile* out_file = new TFile("output.root");

  		TTree* tree_copy = tree->CopyTree(selection.c_str(), "Copy");

				if (is_bg_chain == true)
				{
						factory->AddBackgroundTree(tree_copy, 1.0);
				}
				else
				{
						factory->AddSignalTree(tree_copy, 1.0);
				}

				out_file->Write();
				out_file->Close();
				std::remove(out_file->GetName());
  }
}

const char* BDTAnalysis::create_BDT(DataChain* bg_chain, DataChain* signal_chain, std::vector<Variable*>* variables,
																																				std::string folder_name, const char* NTrees, const char* BoostType,
																																				const char* AdaBoostBeta, const char* SeparationType, const char* nCuts,
																																				std::string job_name)
{
	  if (!opendir(folder_name.c_str()))
	  {
	    mkdir(folder_name.c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
	  }

	  std::string output_path = BDT_output_file_path(folder_name, job_name, true, NTrees, BoostType,
																																																		AdaBoostBeta, SeparationType, nCuts, bg_chain->label);

	  TFile* output_tmva = TFile::Open(output_path.c_str(),"RECREATE");

   std::string weight_file_name = "TMVAClassification" + job_name;

	  TMVA::Factory* factory = new TMVA::Factory(weight_file_name.c_str(), output_tmva,
	                                             "!V:!Silent:Color:DrawProgressBar:Transformations=I;D;P;G,D:AnalysisType=Classification");

	  for (int i = 0; i < variables->size(); i++)
	  {
	    factory->AddVariable((*variables)[i]->name, (*variables)[i]->name_styled, (*variables)[i]->units, 'F');
	  }

	  // Background
	  add_trees_to_factory(bg_chain, true, factory);
	  factory->SetBackgroundWeightExpression("total_weight_lepveto");

	  // Signal
	  add_trees_to_factory(signal_chain, false, factory);
	  factory->SetSignalWeightExpression("total_weight_lepveto");

	  // Apply additional cuts on the signal and background samples (can be different)
	  TCut signal_cuts = "alljetsmetnomu_mindphi>2.0 && jet1_pt>50.0 && jet2_pt>45.0 && metnomu_significance>3.5 && dijet_deta>4.2 && dijet_deta<8.0 && nvetomuons==0 && nvetoelectrons==0"; // for example: TCut mycuts = "abs(var1)<0.5 && abs(var2-0.5)<1";
	  TCut bg_cuts = signal_cuts; // for example: TCut mycutb = "abs(var1)<0.5";

	  factory->PrepareTrainingAndTestTree("", "", "SplitMode=Random:NormMode=NumEvents:!V" );
 
   factory->BookMethod(TMVA::Types::kBDT, "BDT", BDT_options_str(NTrees, BoostType, AdaBoostBeta, SeparationType, nCuts));

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

			return output_tmva->GetName();
}

const char* BDTAnalysis::evaluate_BDT(DataChain* bg_chain, std::vector<Variable*>* variables, std::string output_name,
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

	   /*Float_t jet1_eta;
	   Float_t jet2_eta;
	   Float_t jet1_phi;
	   Float_t jet2_phi;*/
	   Float_t jet_csv1;
	   Float_t jet_csv2;
	   Float_t dijet_dphi;
	   //Float_t metnomu_x;
	   //Float_t metnomu_y;
	   //Float_t sumet;
	   Float_t mht;
	   Float_t unclustered_et;
	   Float_t jetmetnomu_mindphi;
	   //Float_t jetunclet_mindphi;
	   Float_t metnomuunclet_dphi;
	   //Float_t dijetmetnomu_vectorialSum_pt;
	   Float_t dijetmetnomu_ptfraction;
	   Float_t jet1metnomu_scalarprod;
	   Float_t jet2metnomu_scalarprod;

	   reader->AddVariable("alljetsmetnomu_mindphi", &alljetsmetnomu_mindphi);
	   reader->AddVariable("forward_tag_eta", &forward_tag_eta);
	   reader->AddVariable("dijet_deta", &dijet_deta);
	   reader->AddVariable("metnomu_significance", &metnomu_significance);
	   reader->AddVariable("sqrt_ht", &sqrt_ht);
	   reader->AddVariable("dijet_M", &dijet_M);
	   reader->AddVariable("metnomuons", &metnomuons);

	   reader->AddVariable("jet1_pt", &jet1_pt);
	   reader->AddVariable("jet2_pt", &jet2_pt);

	   //reader->AddVariable("jet1_eta", &jet1_eta);
	   //reader->AddVariable("jet2_eta", &jet2_eta);
	   //reader->AddVariable("jet1_phi", &jet1_phi);
	   //reader->AddVariable("jet2_phi", &jet2_phi);
	   reader->AddVariable("jet_csv1", &jet_csv1);
	   reader->AddVariable("jet_csv2", &jet_csv2);
	   reader->AddVariable("dijet_dphi", &dijet_dphi);
	   //reader->AddVariable("metnomu_x", &metnomu_x);
	   //reader->AddVariable("metnomu_y", &metnomu_y);
	   //reader->AddVariable("sumet", &sumet);
	   reader->AddVariable("mht", &mht);
	   reader->AddVariable("unclustered_et", &unclustered_et);
	   reader->AddVariable("jetmetnomu_mindphi", &jetmetnomu_mindphi);
	   //reader->AddVariable("jetunclet_mindphi", &jetunclet_mindphi);
	   reader->AddVariable("metnomuunclet_dphi", &metnomuunclet_dphi);
	   //reader->AddVariable("dijetmetnomu_vectorialSum_pt", &dijetmetnomu_vectorialSum_pt);
	   reader->AddVariable("dijetmetnomu_ptfraction", &dijetmetnomu_ptfraction);
	   reader->AddVariable("jet1metnomu_scalarprod", &jet1metnomu_scalarprod);
	   reader->AddVariable("jet2metnomu_scalarprod", &jet2metnomu_scalarprod);

	   std::string weight_file_path = "weights/TMVAClassification" + job_name + "_BDT.weights.xml";

	   // Book method(s)
	   reader->BookMVA( "BDT method", weight_file_path.c_str() );

	   // Book output histograms
	   TH1F* histBdt     = new TH1F( "MVA_BDT", "MVA_BDT", 100, -0.8, 0.8 );

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
	   /*data->SetBranchAddress("jet1_eta", &jet1_eta);
	   data->SetBranchAddress("jet2_eta", &jet2_eta);
	   data->SetBranchAddress("jet1_phi", &jet1_phi);
	   data->SetBranchAddress("jet2_phi", &jet2_phi);*/
	   data->SetBranchAddress("jet_csv1", &jet_csv1);
	   data->SetBranchAddress("jet_csv2", &jet_csv2);
	   data->SetBranchAddress("dijet_dphi", &dijet_dphi);
	   //data->SetBranchAddress("metnomu_x", &metnomu_x);
	   //data->SetBranchAddress("metnomu_y", &metnomu_y);
	   //data->SetBranchAddress("sumet", &sumet);
	   data->SetBranchAddress("mht", &mht);
	   data->SetBranchAddress("unclustered_et", &unclustered_et);
	   data->SetBranchAddress("jetmetnomu_mindphi", &jetmetnomu_mindphi);
	   //data->SetBranchAddress("jetunclet_mindphi", &jetunclet_mindphi);
	   data->SetBranchAddress("metnomuunclet_dphi", &metnomuunclet_dphi);
	   //data->SetBranchAddress("dijetmetnomu_vectorialSum_pt", &dijetmetnomu_vectorialSum_pt);
	   data->SetBranchAddress("dijetmetnomu_ptfraction", &dijetmetnomu_ptfraction);
	   data->SetBranchAddress("jet1metnomu_scalarprod", &jet1metnomu_scalarprod);
	   data->SetBranchAddress("jet2metnomu_scalarprod", &jet2metnomu_scalarprod);

	   std::string target_name;

	   if (unique_output_files)
	   {
	     target_name = output_name;
	   }
	   else
	   {
	     target_name = "TMVApp.root";
	   }

	   TFile* target = TFile::Open(target_name.c_str(),"RECREATE");

	   Float_t output;
	   TTree* output_tree = new TTree("MVAtree","Tree with classifier outputs");
	   output_tree->Branch("output", &output, "output");

	   std::cout << "--- Processing: " << data->GetEntries() << " events" << std::endl;
	   TStopwatch sw;
	   sw.Start();
	   for (Long64_t ievt=0; ievt<data->GetEntries(); ievt++)
	   {

	      if (ievt%1000 == 0) std::cout << "--- ... Processing event: " << ievt << std::endl;

	      data->GetEntry(ievt);

	      output = reader->EvaluateMVA( "BDT method" );
	      output_tree->Fill();
	      histBdt->Fill(output);
	   }

	   // Get elapsed time
	   sw.Stop();
	   std::cout << "--- End of event loop: "; sw.Print();

	   // --- Write histograms
	   target->cd();
	   output_tree->Write();
	   histBdt->Write();

	   target->Close();

	   std::cout << "--- Created root file: " << output_name << "containing the MVA output histograms" << std::endl;

	   delete reader;

	   std::cout << "==> TMVAClassificationApplication is done!" << std::endl;

	   return target->GetName();
}

//note before calling this method you must call create_bdt to update the xml weight file:
DataChain* BDTAnalysis::get_BDT_results(DataChain* bg_chain, std::vector<Variable*>* variables, std::string output_name,
																																								std::string job_name, bool unique_output_files)
{
	 const char* output_path = BDTAnalysis::evaluate_BDT(bg_chain, variables, output_name, job_name, unique_output_files);
	 TFile* output_file = new TFile(output_path);
	 TTree* output_weight = (TTree*) output_file->Get("MVAtree");
	 TChain* bg_clone     = (TChain*) bg_chain->chain->Clone();

	 bg_clone->AddFriend(output_weight);

	 std::string label(bg_chain->label);
	 label += "_w_mva_output";

  DataChain* output_data = new DataChain(z_ll, bg_chain->label, bg_chain->legend, bg_chain->lep_sel, label, bg_clone);

  output_file->Close();
  //std::remove(output_path);

	 return output_data;
}

std::string BDTAnalysis::BDT_options_str(const char* NTrees, const char* BoostType,
																																									const char* AdaBoostBeta, const char* SeparationType, const char* nCuts)
{
	std::string BDT_options = "!H:!V:NTrees=";
	BDT_options.append(NTrees);
	BDT_options += ":MinNodeSize=2.5%:MaxDepth=3:BoostType=";
	BDT_options.append(BoostType);
	BDT_options += ":AdaBoostBeta=";
	BDT_options.append(AdaBoostBeta);
	BDT_options += ":UseBaggedBoost:BaggedSampleFraction=";
	BDT_options.append(AdaBoostBeta);
	BDT_options += ":SeparationType=";
	BDT_options.append(SeparationType);
	BDT_options += ":nCuts=";
	BDT_options.append(nCuts);

 return BDT_options;
}

std::string BDTAnalysis::BDT_output_file_path(std::string folder_name, std::string job_name, bool is_train_file,
																																													 const char* NTrees, const char* BoostType,
																																														const char* AdaBoostBeta, const char* SeparationType, const char* nCuts,
																																														const char* bg_chain_label)
{
  std::string file_name = BDT_output_name_str(NTrees, BoostType, AdaBoostBeta, SeparationType, nCuts, bg_chain_label, job_name);

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

std::string BDTAnalysis::BDT_output_name_str(const char* NTrees, const char* BoostType,const char* AdaBoostBeta,
																																													const char* SeparationType, const char* nCuts, const char* bg_chain_label,
																																													std::string job_name)
{
 std::string bg = bg_chain_label;
 std::string out_nam = "BDT-job_name=" + job_name + "-" + bg + "-NTrees=";
	out_nam.append(NTrees);
	out_nam += "-BoostType=";
	out_nam.append(BoostType);
	out_nam += "-AdaBoostBeta=";
	out_nam.append(AdaBoostBeta);
	out_nam += "-SeparationType=";
	out_nam.append(SeparationType);
	out_nam += "-nCuts=";
	out_nam.append(nCuts);
	out_nam += ".root";

	return out_nam;
}
