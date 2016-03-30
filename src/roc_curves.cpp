#include "../include/roc_curves.h"

double RocCurves::get_auc(std::string method_name, const char* training_output_path)
{
	  TFile* training_output = TFile::Open(training_output_path);
	  TH1D* roc = plot_rejBvsS(training_output, method_name);
   double auc = 0;

   if (roc != 0x0)
   {
     auc = roc->Integral();
   }

   return auc;
}

void RocCurves::get_rocs(std::vector<TFile*> training_outputs, DataChain* signal,
	DataChain* bg, SuperVars* super_vars, std::string method_name, std::string dir)
{
	std::vector<Variable*> vars      = super_vars->get_discriminating_vars();
	std::string preselection         = super_vars->get_final_cuts_str();
	std::string var_cut_str_tmva     = super_vars->get_cuts_str_for_tmva();

	plot_all_rejBvsS(training_outputs, signal, bg, preselection, vars[0], &vars, method_name, dir);
}

void RocCurves::plot_all_rejBvsS(std::vector<TFile*> training_outputs, DataChain* signal,
	DataChain* bg, std::string preselection, Variable* var, std::vector<Variable*>* variables,
	std::string method_name, std::string dir)
{
	std::string plot_name = "ROC Curves";
	TCanvas* c1     = new TCanvas("c1", plot_name.c_str());
	TLegend* legend = new TLegend(0.15, 0.6, 0.65, 0.15);
	TGraph* point   = parked_data_point(signal, bg, preselection, var, variables, legend);
	std::cout << "=> Parked Data Point done" << std::endl;
	for (int i = 0; i < training_outputs.size(); i++)
	{
		std::string training_name(training_outputs[i]->GetName());
		std::string file_name = training_name.substr(training_name.find("/") + 1 , -1);
		std::string legend_str = HistoPlot::replace_all(file_name, ".root", "");
		std::string legend_text = "#splitline{" + legend_str.substr(0, 45) + "}{" + legend_str.substr(46, -1) + "}";

		TH1D* rejBvsS_histo    = plot_rejBvsS(training_outputs[i], method_name);
		rejBvsS_histo->SetLineColor(HistoPlot::colours()[i+2]);
		rejBvsS_histo->SetStats(false);
		rejBvsS_histo->SetLineWidth(3);
		rejBvsS_histo->SetTitle("ROC Curves");

		if (i == 0)
		{
			rejBvsS_histo->Draw();
		}
		else
		{
			rejBvsS_histo->Draw("SAME");
		}

		legend->AddEntry(rejBvsS_histo, legend_text.c_str(), "f");
	}
	style_legend(legend);
	point->Draw("pSAME");
	legend->Draw();
	std::string bg_name = bg->label;
	std::string img_name = dir + "/" + bg_name + "_roc_curves.png";
	c1->SaveAs(img_name.c_str());
	c1->Close();
}

TH1D* RocCurves::plot_rejBvsS(TFile* training_output, std::string method_name)
{
	std::string file_path;
	if (method_name == "MLP")
	{
		file_path = "Method_MLP/MLP/MVA_MLP_rejBvsS";
	}
	else if (method_name == "BDT")
	{
		file_path = "Method_BDT/BDT/MVA_BDT_rejBvsS";
	}
	TH1D* rejBvsS_histo = (TH1D*) training_output->Get(file_path.c_str());
	std::cout << "Training output name " << training_output->GetName() << std::endl;
	std::cout << "test " << ", " << (TH2F*) training_output->Get("CorrelationMatrixS;1") << std::endl;
	return rejBvsS_histo;
}

double RocCurves::get_presel_effy(DataChain* data_chain, std::string preselection, Variable* var, std::vector<Variable*>* variables)
{
	std::string parked_data = "((nvetomuons==0)&&(nvetoelectrons==0)&&(jet1_eta<4.7)&&(jet2_eta<4.7)&&(jet1_pt > 50)&&(jet2_pt > 45)&&(dijet_M > 1200)&&(metnomuons > 90)&&(metnomu_significance > 4.0)&&(alljetsmetnomu_mindphi > 2.3)&&(dijet_deta>4.2))*total_weight_lepveto";
	TH1F* remaining_histo = histo(data_chain, var, parked_data);
	double remaining = remaining_histo->Integral();
	std::cout << "remaining integral = " << remaining << std::endl;

	TH1F* total_histo = histo(data_chain, var, preselection);
	double total = total_histo->Integral();
	std::cout << "total integral" <<  total << std::endl;

	double effy = remaining / total;
	std::cout << effy << std::endl;
	return effy;
}

TGraph* RocCurves::parked_data_point(DataChain* signal, DataChain* bg, std::string preselection, Variable* var,
	std::vector<Variable*>* variables, TLegend* legend)
{
	double rejB = 1 - get_presel_effy(bg, preselection, var, variables);
	double sigeff = get_presel_effy(signal, preselection, var, variables);
	std::cout << "here\n";
	TGraph* point = new TGraph(1);
	point->SetPoint(1, sigeff, rejB);
	point->SetMarkerStyle(20);
	point->SetMarkerColor(kRed);
	legend->AddEntry(point, "Original Preselection", "p");

	return point;
}

TH1F* RocCurves::histo(DataChain* data, Variable* variable, std::string preselection)
{
	data->chain->Draw((variable->build_var_string("test", true)).c_str(), preselection.c_str());
	TH1F* histo = (TH1F*) gDirectory->Get("test");
	return histo;
}

void RocCurves::style_legend(TLegend* legend)
{
	legend->SetTextSize(0.03);
	legend->SetBorderSize(0);
	legend->SetFillStyle(0);
}
/*
THStack HistoPlot::draw_parked_stack(Variable* var, std::vector<DataChain*> bg_chains,
									std::vector<Variable*>* parked_variables, DataChain* data, std::string parked_selection)
{
  THStack stack("parked", "");
  std::vector<double> mc_weights_vector = mc_weights(data, bg_chains, var, true, parkedvariables, mva_cut);

  for(int i = 0; i < bg_chains.size(); i++) {
  	std::string selection = parked_selection + "*" + HistoPlot::get_string_from_double(mc_weights_vector[i]);
	bg_chains[i]->chain->Draw((variable->build_var_string(bg_chains[i]->label, true)).c_str(), selection.c_str());
	TH1F* histo = (TH1F*) gDirectory->Get(bg_chains[i]->label);																																						
    stack.Add(histo);
  }
  return stack;
}
*/
