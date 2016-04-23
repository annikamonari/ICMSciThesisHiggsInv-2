#include "../include/histo_plot.h"
#include <algorithm>
using namespace std;
void HistoPlot::draw_plot(Variable* var, std::vector<DataChain*> bg_chains,
                          DataChain* signal_chain, DataChain* data, bool with_cut,
                          std::vector<Variable*>* variables, bool plot_data, std::string file_name,																					std::string mva_cut)
{
  TCanvas* c1     = new TCanvas("c1", var->name_styled, 800, 800);
  TPad* p1        = new TPad("p1", "p1", 0.0, 0.95, 1.0, 1.0);
  TPad* p2        = new TPad("p2", "p2", 0.0, 0.2, 1.0, 0.95);
  TPad* p3        = new TPad("p3", "p3", 0.0, 0.0, 1.0, 0.2);
  TLegend* legend = new TLegend(0.0, 0.5, 0.0, 0.88);

  p1->SetLeftMargin(0.102);
  p2->SetBottomMargin(0.012);
  p2->SetLeftMargin(0.105);
  p3->SetBottomMargin(0.3);
  p3->SetLeftMargin(0.102);

  p1->Draw();
  p2->Draw();
  p3->Draw();
  p2->cd();
  //std::cout << "mva cut: " << mva_cut << std::endl;
  THStack stack      = draw_stacked_histo(legend, var, bg_chains, with_cut, variables, data, mva_cut);
  std::cout << "drew stack" << std::endl;
  TH1F* signal_histo = draw_signal(signal_chain, var, with_cut, legend, variables, mva_cut);
  std::cout << "drew signal" << std::endl;
  //std::cout << signal_histo << std::endl;
  TH1F* data_histo;
  if (plot_data) {data_histo = draw_data(data, var, with_cut, legend, variables, mva_cut);}
  //stack.Add(signal_histo);
  stack.Draw();
  signal_histo->Draw("SAME");

  cout<<"total data ="<<data_histo->Integral()<<"\n";
  if (plot_data) {data_histo->Draw("SAME");}

  std::cout << "drew all" << std::endl;
  style_stacked_histo(&stack, var->name_styled);
//cout<<"histo styled"<<endl;
  TH1F* plot_histos[3] = {(TH1F*)(stack.GetStack()->Last()), data_histo, signal_histo};
  std::vector<TH1F*> plot_histos_vector (plot_histos, plot_histos + sizeof(plot_histos) / sizeof(plot_histos[0]));
  TH1F* max_histo      = get_max_histo(plot_histos_vector);
  std::cout << "got max" << max_histo << std::endl;
  stack.SetMaximum(get_histo_y_max(max_histo)*1.5);

  build_legend(legend, signal_histo, var, with_cut);
cout<<"legend built"<<endl;
//  draw_subtitle(var, variables, with_cut, data, "", mva_cut);
  std::cout << "drew subtitle" << std::endl;
  p3->cd();
  TH1F* data_bg_ratio_histo;

  /*if (plot_data)
  {
  		std::cout << "data histo isnt null" << std::endl;
  		data_bg_ratio_histo = data_to_bg_ratio_histo(data_histo, (TH1F*)(stack.GetStack()->Last()));
  }
  else
  {*/
  //		std::cout << "data_histo appaz is null" << std::endl;
  		data_bg_ratio_histo = data_to_bg_ratio_histo(signal_histo, (TH1F*)(stack.GetStack()->Last()));
 // }

  data_bg_ratio_histo->Draw("e1");
  style_ratio_histo(data_bg_ratio_histo, var->name_styled);
  draw_yline_on_plot(var, with_cut, 1.0);
  std::cout << "ratio histo done" << std::endl;
  std::string img_name;

  if (file_name == "")
  {
  		img_name = build_file_name(var, with_cut);
  }
  else
  {
  		img_name = file_name;
  }
  //std::cout << "file name" << std::endl;
  p1->cd();
  draw_title(var->name_styled);
  c1->SaveAs(img_name.c_str());
  c1->Close();
}
//////////////////////////////////////////////////////////////////////
////////////Plots a histogram of the evaluated test data for zjets vv with signal 
/////////// in red overlaid on top of grey background, als includes additional signal/nbackground plot below
///////////
void HistoPlot::plot_evaluated_zjets_vv_testTree(int bg_trained, Variable* mva_output, DataChain* testTree_chain,
DataChain* data, std::vector<DataChain*> bg_chains,std::vector<Variable*>* variables, std::string file_name, std::string mva_cut)
{
  //////////////////////////////////////////////////////////////////
  /*cout<<"---------------------Started HistoPlot::plot_evaluated_zjets_vv_testTree --------------------------------"<<endl;
  //step 1: initialise TCanvas
  TCanvas* c1     = new TCanvas("c1", mva_output->name_styled, 800, 800);
  TPad* p1        = new TPad("p1", "p1", 0.0, 0.95, 1.0, 1.0);
  TPad* p2        = new TPad("p2", "p2", 0.0, 0.2, 1.0, 0.95);
  TPad* p3        = new TPad("p3", "p3", 0.0, 0.0, 1.0, 0.2);
  TLegend* legend = new TLegend(0.0, 0.5, 0.0, 0.88);

  p1->SetLeftMargin(0.102);
  p2->SetBottomMargin(0.012);
  p2->SetLeftMargin(0.105);
  p3->SetBottomMargin(0.3);
  p3->SetLeftMargin(0.102);

  p1->Draw();
  p2->Draw();
  p3->Draw();
  p2->cd();
 // std::cout << "mva cut: " << mva_cut << std::endl;
   //step 1.2 clone data chain
	 TChain* test_clone     = (TChain*) testTree_chain->chain->Clone();
 DataChain* clone_chain  = new DataChain(mc_signal_data, mc_signal_label, mc_signal_legend, "","",test_clone);
   std::cout << "step 1 done" << std::endl;


  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //step 2: create background histo
  //step 2.1: get zjets mc weight
   // don't need to double the trained bg in mc weights becauase this takes in test + train trees
  std::vector<double> mc_weights_vector = mc_weights(data, bg_chains, mva_output, true, NULL, mva_cut,bg_trained, true);

  double trained_mc_weight = mc_weights_vector[bg_trained];

  //cout<<"mc_weight of trained bg= "<< trained_mc_weight <<endl;
//mc weights have been tested and work, now need to edit the draw stack function to not stack the zjets in the datachain vector and instead add the test_set th1
  //step 2.2 create output selection string
  std::string selection = "((alljetsmetnomu_mindphi>2.0)&&(nvetomuons==0)&&(nvetoelectrons==0))*2*total_weight_lepveto";
  selection = add_classID_to_selection(selection, false);
  //selection  += "*total_weight_lepveto";
  //step 2.3 add mva cut to selection string
  //step 2.3.1 reformat mva_cut string for MLPs

  selection = HistoPlot::add_mva_cut_to_selection(selection, mva_cut);
  selection = add_mc_to_selection(testTree_chain,mva_output, selection, trained_mc_weight);

  //step 2.4 add mc weight to selection  */
  /*std::string mc_weight_str = get_string_from_double(mc_weight);
  selection = selection + "*" + mc_weight_str; */
  //std::cout << "selection for test tree background =============" << std::endl;
  //std::cout<<"final selection: "<< selection <<"\n";
/*
  //step 2.5: get background histo  problem area caused by signal using the same chain, try cloning the tree and running agsain
  testTree_chain->chain->SetLineColor(1);
  testTree_chain->chain->SetFillColor(colours()[bg_trained]);
  TH1F* trained_histo = build_1d_histo(testTree_chain, mva_output, true, true, "goff", NULL, selection, trained_mc_weight, mva_cut);
  TH1F bg_histo_inThe_memory_Stack = *trained_histo;// save to stack emory fromn the heap
  //std::cout << "events in zjets vv: ==========" << trained_histo->Integral() << std::endl;
  //step 2.6 add legend entry
  std::string legend_str(testTree_chain->label);
  legend_str += (" #font[12]{(MC weight: " + get_string_from_double(trained_mc_weight) + ")}");// get string from double fails so mc weight is added manually
  legend->AddEntry(trained_histo, legend_str.c_str(), "f");

   //cout<<"legend str: "<<legend_str<<"\n";
   //step 2.7 get stack for all non zjets_vv backgrounds here

   THStack stack = draw_stacked_histo_no_zjets(legend, mva_output, bg_chains, true, variables, data, mc_weights_vector, testTree_chain, mva_cut);

   //step 2.8 add zjets_vv to stack
   stack.Add(trained_histo);
   std::cout << "step 2 done" << std::endl;


  ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //step 3: create signal histo 
  //step 3.1 get regular selection string
  string sig_selection = "((alljetsmetnomu_mindphi>2.0)&&(nvetomuons==0)&&(nvetoelectrons==0))*2*total_weight_lepveto";
  sig_selection = add_classID_to_selection(sig_selection, true);
  //step 3.2: add class ID to selection
  sig_selection = HistoPlot::add_mva_cut_to_selection(sig_selection, mva_cut);

  //std::cout<<"final sig_selection: "<< sig_selection<<"\n";

  //step 3.3 get signal his0-HiddenLayers=2,2,2,2-LearningRate=0.02-EstimatorType=CEoutput_nocuts.pngto 
  test_clone->SetLineColor(2);
  test_clone->SetLineWidth(3);
  test_clone->SetFillColor(0);
  TH1F* signal_histo = build_1d_histo(clone_chain, mva_output, true, true, "goff", NULL , sig_selection, 1, mva_cut);
  TH1F signal_histo_inTheStack = *signal_histo;
  //step 3.4 add legend entry for signal
  std::string signal_leg_str = "Signal";
  signal_leg_str += " (x";
  signal_leg_str.append(mva_output->signal_multiplier);
  signal_leg_str += ")";

  legend->AddEntry(signal_histo, (signal_leg_str).c_str(), "l");
  std::cout << "step 3 done" << std::endl;
  ////////////////////////////////////////////////////////
  //step 3.5 add data histo
  TH1F* data_histo;
  string data_selection = "((alljetsmetnomu_mindphi>2.0)&&(nvetomuons==0)&&(nvetoelectrons==0))";
  data_selection = HistoPlot::add_mva_cut_to_selection(data_selection, mva_cut);
  data_histo = draw_data(data, mva_output, true, legend, NULL, mva_cut,data_selection);
  /////////////////////////////////////////////////////////
  //step 4: draw and style primary histogram
  //step 4.1 draw background histogram
  stack.Draw();
//trained_histo->Draw();
  //std::cout << "step 4.1 done" << std::endl;

  //step 4.2 draw signal histogram
  signal_histo->Draw("SAME");
 // data_histo->Draw("SAME");

  //std::cout << "step 4.2 done" << std::endl;

  //step 4.3 style histogram
  style_stacked_histo(&stack, mva_output->name_styled);
  std::cout << "step 4 done" << std::endl;

  ////////////////////////////////////////////////////////////
  //step 5: build legend and draw substitle
  //step 5.1 build legend
  build_legend(legend, signal_histo, mva_output, true);

  //step 5.2 draw subtitles
  draw_subtitle(mva_output, variables, true, data, "", mva_cut);
  std::cout << "step 5 done" << std::endl;


  /////////////////////////////////////////////////////////////////////////////
  //step 6: draw signal/background secondary histogram
  //step 6.1 change pad then create signal/background histogram
  p3->cd();
  TH1F* signal_bg_ratio_histo = data_to_bg_ratio_histo(signal_histo,*//*trained_histo *//*(TH1F*)(stack.GetStack()->Last()));

  //step 6.2 draw signal/background histo
  signal_bg_ratio_histo->Draw("e1");

  //step 6.3 style signal/background histo
  style_ratio_histo(signal_bg_ratio_histo, mva_output->name_styled);
  //std::cout << "ratio histo done" << std::endl;
  //step 6.4 draw line, y=1 on signal/background histogram
  draw_yline_on_plot(mva_output, true, 1.0);
  draw_yline_on_plot(mva_output, true, 1.0);
  std::cout << "step 6 done" << std::endl;

  ////////////////////////////////////////////////
  //step 7: add title and save histogram to file
  //step 7.1 create file name
  std::string img_name = file_name;
 // std::cout << "file name" << std::endl;
  
  //step 7.2 change pad and add title
  p1->cd();
  draw_title(mva_output->name_styled);
  //
  //step 7.3 close and save the file
  c1->SaveAs(img_name.c_str());
//cout<<"file saved\n";
  c1->Close();
  std::cout << "step 7 done" << std::endl;
cout<<"---------------------Finished HistoPlot::plot_evaluated_zjets_vv_testTree --------------------------------"<<endl;*/
}


//_______________________________________________________________________________________________________________________

std::string HistoPlot::add_classID_to_selection(std::string selection, bool is_signal)
{
  std::string classID;
  if (is_signal)
  {
    classID = "(classID==1)&&";
  }
  else
  {
    classID = "(classID==0)&&";
  }
  selection.insert(selection.find("(") + 1, classID);
  
  return selection;
}


//_______________________________________________________________________________________________________________________

THStack HistoPlot::draw_stacked_histo_no_zjets(TLegend* legend, Variable* var, std::vector<DataChain*> bg_chains,
                                      bool with_cut, std::vector<Variable*>* variables, DataChain* data, std::vector<double> mc_weights_vector,
																																						DataChain* testTree_chain, std::string mva_cut)
{
  THStack stack(var->name_styled, "");
  std::string selection = "((alljetsmetnomu_mindphi>2.0)&&(nvetomuons==0)&&(nvetoelectrons==0))*total_weight_lepveto";
  selection = add_classID_to_selection(selection, false);
  selection = HistoPlot::add_mva_cut_to_selection(selection, mva_cut);
 // std::cout << "selection in draw stacked histo " << selection << std::endl;
  for(int i = 0; i < bg_chains.size(); i++) {
   if(strcmp(bg_chains[i]->label, testTree_chain->label)) {
       //cout<<bg_chains[i]->label<<" === background, mc weight: "<<mc_weights_vector[i]<<endl;
      TH1F* single_bg_histo = draw_background_from_trees(bg_chains[i], var, colours()[i], selection, mc_weights_vector[i], mva_cut);
      stack.Add(single_bg_histo);
      std::string legend_str(bg_chains[i]->legend);
      legend_str += (" #font[12]{(MC weight: " + get_string_from_double(mc_weights_vector[i]) + ")}");
      legend->AddEntry(single_bg_histo, legend_str.c_str(), "f");
      //cout<<"drew background: "<<bg_chains[i]->label<<endl;
    }
  }
  return stack;
}
//_______________________________________________________________________________________________________________________


void HistoPlot::draw_yline_on_plot(Variable* var, bool with_cut, double y)
{
  double x_min = 0.0;
  double x_max = 1.0;

  if (with_cut)
  {
    x_min = atof(var->x_min_cut);
    x_max = atof(var->x_max_cut);
  }
  else
  {
    x_min = atof(var->x_min_nocut);
    x_max = atof(var->x_max_nocut);
  }

  TLine *line = new TLine(x_min, y, x_max, y);
  line->SetLineColor(13);
  line->SetLineStyle(2);
  line->Draw("SAME");
}

void HistoPlot::draw_title(const char* title)
{
  TPaveText* pt = new TPaveText(0.1, 0.1, 0.9, 1.0, "blNDC");
  pt->SetBorderSize(0);
  pt->SetFillColor(0);
  pt->AddText(title);
  pt->SetAllWith(title, "size", 0.8);
  pt->Draw();
}
//_______________________________________________________________________________________________________________________

std::string HistoPlot::get_selection(Variable* variable, std::vector<Variable*>* variables,
                                     bool with_cut, bool is_signal, DataChain* bg_chain, double mc_weight,
																																					std::string mva_cut)
{
  std::string selection;

  if ((variables != NULL) && (with_cut))
  {
    selection = variable->build_multicut_selection(is_signal, variables);
  }
  else
  {
    selection = variable->build_selection_string(with_cut, is_signal);
  }
//std::cout << selection << std::endl;

  selection.insert(selection.find("(") + 1, lep_sel_default());
  selection = HistoPlot::add_mva_cut_to_selection(selection, mva_cut);
  
  return add_mc_to_selection(bg_chain, variable, selection, mc_weight);
}
//_______________________________________________________________________________________________________________________

std::string HistoPlot::get_parked_selection(Variable* variable, std::vector<Variable*>* variables,  DataChain* bg_chain, double mc_weight)
{
  std::string selection;

  selection = variable->build_parked_selection( variables);

  selection.insert(selection.find("(") + 1, lep_sel_default());
  //std::cout << selection << std::endl;

  return add_mc_to_selection(bg_chain, variable, selection, mc_weight);

}
//_______________________________________________________________________________________________________________________

std::string HistoPlot::add_mc_to_selection(DataChain* bg_chain, Variable* variable, std::string selection, double mc_weight)
{
  std::string mc_weight_str = get_string_from_double(mc_weight);
  std::string sel_new = selection += "*" + mc_weight_str; //selection.insert(selection.find("*") + 1, mc_weight_str + "*");

  return sel_new;
}

std::string HistoPlot::add_mva_cut_to_selection(std::string selection, std::string mva_cut_str)
{
  if (mva_cut_str != "")
  {
    std::string mva_cut = "(";
    mva_cut.append(mva_cut_str);
    mva_cut.append(")&&");
    selection.insert(selection.find("(") + 1, mva_cut);
  }

  return selection;
}
//_______________________________________________________________________________________________________________________

std::vector<double> HistoPlot::mc_weights(DataChain* ewk_chain, DataChain* qcd_chain,DataChain* data, std::vector<DataChain*> bg_chains, Variable* var, bool with_cut, 
std::vector<Variable*>* variables,std::string mva_cut, int trained_bg, bool double_test_bg, bool if_parked)
{//cout<<"in HistoPlot::mc_weights\n";
  double mc_weight[8];
  double e_f = HistoPlot::get_efficiency_factor( ewk_chain,  qcd_chain,var, variables, mva_cut);
  std::vector<double> mc_weights_vector (mc_weight, mc_weight + sizeof(mc_weight) / sizeof(mc_weight[0]));

  for(int i =0; i< 8;i++){mc_weights_vector[i]=1;}
  mc_weights_vector[0] = MCWeights::calc_mc_weight(mc_weights_vector,0, data, bg_chains, bg_chains[0], var, with_cut,
						variables, mva_cut, trained_bg, double_test_bg,if_parked)*e_f;
//cout<<"mc weights:"<<mc_weights_vector[0]<<"\n";

  mc_weights_vector[6] = MCWeights::calc_nunu_mc_weight(data, bg_chains, bg_chains[6], var, 
 with_cut, variables, mva_cut, trained_bg, double_test_bg,if_parked);

  for(int i =1; i < 4;i++)
  {
  	mc_weights_vector[i] = MCWeights::calc_mc_weight(mc_weights_vector,i, data, bg_chains, bg_chains[i], var, with_cut,
						variables, mva_cut, trained_bg, double_test_bg,if_parked);
//cout<<","<<mc_weights_vector[1];
//cout<<","<<mc_weights_vector[2];
//cout<<","<<mc_weights_vector[3];
//cout<<","<<mc_weights_vector[6]<<"\n";
  }
 /* for(int i =0; i < 8;i++)
  {
    std::cout << bg_chains[i]->label <<": "<<mc_weights_vector[i]<<"\n";
  }*/

  return mc_weights_vector;
}
//_______________________________________________________________________________________________________________________

// new function written just like the one above: HistoPlot::mc_weights, which calculates the right error for the bgs without a control
// region (its just sqrt(unweighted mc events in signal) / unweighted mc events in signal)
// note:: changed so that if mc weight is 1 then dont calculate the mc weight error
// TODO make error use sumw2 and integralanderror
std::vector<double> HistoPlot::get_bg_errors(int bg_to_train,DataChain* ewk_chain, DataChain* qcd_chain, DataChain* data, std::vector<DataChain*> bg_chains,
Variable* var, bool with_cut, std::vector<Variable*>* variables,std::vector<double> bg_mc_weights, std::string mva_cut, bool if_parked)
{
string selection="";
	 double bg_errors[bg_chains.size()];
	double integral[8];
	 for(int i = 0; i < bg_chains.size();i++)
	 {
		if(variables==NULL)
		{
    			selection = "((alljetsmetnomu_mindphi>2.0)&&(nvetomuons==0)&&(nvetoelectrons==0))*total_weight_lepveto";
    			selection = HistoPlot::add_classID_to_selection(selection, false);
    			selection = HistoPlot::add_mva_cut_to_selection(selection, mva_cut);
//std::cout << "in histoplot mc weight errors" << selection << std::endl;
		}
//cout<<"selection :"<<selection<<"\n";
	   TH1F* histo = build_1d_histo(bg_chains[i], var,with_cut, false, "goff", variables, selection, 1, mva_cut);
//cout<<"got histo\n";
	   integral[i] = get_histo_integral(histo, with_cut, var);
	 		bg_errors[i] = std::pow(integral[i], 0.5);
//cout<<"got integral: "<<integral<<"\n";

	   if (bg_mc_weights[i] != 1) 
	   {
					if(strcmp(bg_chains[i]->label, "bg_zjets_vv"))
					{
			//cout<<"mc wiehgt in if loop: "<<bg_mc_weights[i]<<endl;
							bg_errors[i] = single_bg_error(ewk_chain, qcd_chain,bg_mc_weights,i, data, bg_chains, bg_chains[i], var,
							with_cut, variables, bg_mc_weights[i], mva_cut,selection,if_parked);
	//cout<<"got bg error\n";

					}
           }

	   if (!strcmp(bg_chains[i]->label, "bg_zjets_vv"))
	   {
             double nunu_weight_error =  MCWeights::calc_nunu_weight_error(ewk_chain, qcd_chain,data, bg_chains, bg_chains[i],
 var,  with_cut,  variables,  mva_cut, if_parked);

	     bg_errors[i] = nunu_weight_error;
// cout << bg_chains[i]->label <<" ========ERROR VALUE========== "<<bg_errors[i]<<endl;

	   }
// cout << bg_chains[i]->label <<" ========ERROR VALUE========== "<<bg_errors[i]<<endl;
	 }
	 std::vector<double> bg_error_vector (bg_errors, bg_errors + sizeof(bg_errors) / sizeof(bg_errors[0]));

	 return bg_error_vector;
}
//_______________________________________________________________________________________________________________________

// problem: TH1F* bg doesn't plot with the mc weight? in the function above this, whenever we call this we pass through
// the mc weight (see last arg: double weight), so if you realise we need it then just put it onto the end of the build_1d_histo call
double HistoPlot::single_bg_error(DataChain* ewk_chain, DataChain* qcd_chain, std::vector<double> bg_mc_weights, int desired_bg_index, DataChain* data, std::vector<DataChain*> bg_chains, DataChain* bg_chain,Variable* var, bool with_cut, std::vector<Variable*>* variables, double weight,
std::string mva_cut, std::string selection, bool if_parked)
{//cout<<"selection: "<<selection<<"\n";
  TH1F* bg = build_1d_histo(bg_chain, var, with_cut, false, "goff", variables, selection, 1, mva_cut);
  //cout << "histo" << bg << endl;
  double MC_N_S = get_histo_integral(bg, with_cut, var);
//cout <<"unweighted background: "<< MC_N_S << endl;
  double sigma_N = std::pow(MC_N_S, 0.5);
//cout << "got sigma_N" << endl;
  double sigma_w = MCWeights::calc_weight_error(bg_mc_weights,desired_bg_index,data, bg_chains, bg_chain, var, with_cut, variables, mva_cut,if_parked);
//cout << "got weight error:"<< sigma_w << endl;
  double sigma_total_sq = std::pow(sigma_w*MC_N_S,2)+std::pow(sigma_N*weight,2);
  double sigma_total = std::pow(sigma_total_sq,0.5);
  //std::cout << bg_chain->label << " - single bg error: " << sigma_total << std::endl;
  if(desired_bg_index==0)
  {
	double e_f = HistoPlot::get_efficiency_factor(ewk_chain,  qcd_chain,var, variables, mva_cut);
	double error_on_e_f = HistoPlot::get_error_on_efficiency_factor( ewk_chain,  qcd_chain,var, variables, mva_cut);
	sigma_total =pow((pow((sigma_total*e_f),2)+pow((MC_N_S*error_on_e_f),2)),0.5);
  }
  return sigma_total;
}
//_______________________________________________________________________________________________________________________

std::string HistoPlot::get_string_from_double(double num)
{
  std::ostringstream num_ss;
  num_ss << num;
  std::string num_str(num_ss.str());

  return num_str;
}
//_______________________________________________________________________________________________________________________

double HistoPlot::sig_to_bg_ratio(Variable* var, TH1F* bg,
                                  TH1F* signal_histo, bool with_cut)
{
  double bg_integral  = get_histo_integral(bg, with_cut, var);
  double sig_integral = get_histo_integral(signal_histo, with_cut, var);
  float signal_mult   = atof(var->signal_multiplier);
  float sig_to_bg     = sig_integral / bg_integral / signal_mult;

  return sig_to_bg;
}
//_______________________________________________________________________________________________________________________

double HistoPlot::get_histo_integral(TH1F* histo, bool with_cut, Variable* var)
{
  int nbins;
  if (with_cut)
  {
    nbins = (int) (atof(var->bins_cut.c_str()) + 0.5);
  }
  else
  {
    nbins = (int) (atof(var->bins_nocut) + 0.5);
  }
  double integral = histo->Integral(0, nbins + 1);

  return integral;
}
//_______________________________________________________________________________________________________________________

std::string HistoPlot::replace_all(std::string str, const std::string& from, const std::string& to)
{

  size_t start_pos = 0;

  while((start_pos = str.find(from, start_pos)) != std::string::npos) {
    str.replace(start_pos, from.length(), to);
    start_pos += to.length(); // Handles case where 'to' is a substring of 'from'
  }
  return str;
}

std::string HistoPlot::style_selection(std::string selection)
{

  std::string sele = replace_all(replace_all(replace_all(replace_all(selection, ")", ""), ">", " > "), "==", " = "), "&&", ", ");

  return replace_all(replace_all(replace_all(replace_all(replace_all(sele, "_", " "), "))", ""), "(", ""), "((", ""), "<", " < ");
}
//_______________________________________________________________________________________________________________________

void HistoPlot::draw_subtitle(Variable* variable, std::vector<Variable*>* variables,
                              bool with_cut, DataChain* data, std::string supervar_selection,
				std::string mva_cut)
{
  std::string sel;
	 if (variables == NULL)
  {
    sel = supervar_selection;
  }
	 else
	 {

			 sel = style_selection(get_selection(variable, variables, with_cut, false, data, 1.0, mva_cut));
	 }
 // std::cout << sel << std::endl;
	 std::string selection = "Selection: " + sel;
TPaveText* pts;
  if(with_cut){
  std::string l1        = "#font[12]{" + selection.substr(0, 90) + "-}";//selection.substr(0, 90)
  std::string l2        = "#font[12]{" + selection.substr(88, 90) + "-}";
  std::string l3        = "#font[12]{" + selection.substr(178, 88) + "}";

  pts        = new TPaveText(0.1, 1.0, 0.9, 0.9, "blNDC");

  pts->SetBorderSize(0);
  pts->SetFillColor(0);
  pts->AddText(l1.c_str());
  pts->AddText(l2.c_str());
  pts->AddText(l3.c_str());
  pts->SetAllWith(l1.c_str(), "size", 0.03);
  pts->SetAllWith(l2.c_str(), "size", 0.03);
  pts->SetAllWith(l3.c_str(), "size", 0.03);
}
else{
  std::string l1        = "#font[12]{" + selection + "-}";//;

  pts        = new TPaveText(0.1, 1.0, 0.9, 0.9, "blNDC");

  pts->SetBorderSize(0);
  pts->SetFillColor(0);
  pts->AddText(l1.c_str());
  pts->SetAllWith(l1.c_str(), "size", 0.03);
}
  pts->Draw();
}

//_______________________________________________________________________________________________________________________

THStack HistoPlot::draw_stacked_histo(TLegend* legend, Variable* var, std::vector<DataChain*> bg_chains,
                                      bool with_cut, std::vector<Variable*>* variables, DataChain* data, std::string mva_cut)
{
  THStack stack(var->name_styled, "");
  double mc_weights_arr[] = { 0.91, 0.73, 0.32, 1.24, 1, 1,1.70, 1};
  if(!with_cut)
  {
  	for(int i=0;i<8;i++){mc_weights_arr[i] = 1;}
  }
  for(int i = 0; i < bg_chains.size(); i++) {
    TH1F* single_bg_histo = draw_background(bg_chains[i], var, colours()[i], with_cut, variables, mc_weights_arr[i],mva_cut);
    stack.Add(single_bg_histo);
    std::string legend_str(bg_chains[i]->legend);
    legend_str += (" #font[12]{(MC weight: " + get_string_from_double(mc_weights_arr[i]) + ")}");
    legend->AddEntry(single_bg_histo, legend_str.c_str(), "f");
  }
  return stack;
}
//_______________________________________________________________________________________________________________________

TH1F* HistoPlot::get_max_histo(std::vector<TH1F*> plot_histos)
{
  double plot_max = 0.0;
  TH1F* histo_max = NULL;
  for (int i = 0; i < plot_histos.size(); i++)
  {
    if (plot_histos[i] != NULL)
    {
					double y_max = get_histo_y_max(plot_histos[i]);

					if (y_max > plot_max)
					{
							plot_max = y_max;
							histo_max = plot_histos[i];
					}
    }
  }

  return histo_max;
}
//_______________________________________________________________________________________________________________________

double HistoPlot::get_histo_y_max(TH1F* histo)
{
  return histo->GetBinContent(histo->GetMaximumBin());
}
//_______________________________________________________________________________________________________________________

void HistoPlot::build_legend(TLegend* legend, TH1F* max_histo, Variable* var, bool with_cut)
{
  double x1 = position_legend_x1(max_histo, var, with_cut);
  double x2 = x1 + 0.26;

  style_legend(legend);

  legend->SetX1(x1);
  legend->SetX2(x2);
  legend->Draw();

}
//_______________________________________________________________________________________________________________________

double HistoPlot::position_legend_x1(TH1F* max_histo, Variable* var, bool with_cut)
{
  int max_bin         = max_histo->GetMaximumBin();
  double nbins        = var->get_bins(with_cut);
  double max_bin_x1   = get_x1_from_bin(max_bin, nbins);

  if (max_bin_x1 > 0.5)
  {
    return 0.12;
  }
  else
  {
    return 0.56;
  }
}
//_______________________________________________________________________________________________________________________

double HistoPlot::get_x1_from_bin(double max_bin, double nbins)
{
  return max_bin * 0.8 / nbins + 0.1;
}
//_______________________________________________________________________________________________________________________

void HistoPlot::style_stacked_histo(THStack* hs, const char* x_label)
{
  hs->GetYaxis()->SetTitle("Events");
  hs->GetYaxis()->SetLabelSize(0.035);
  hs->GetYaxis()->SetTitleOffset(1.55);
  hs->GetXaxis()->SetLabelSize(0);
}
//_______________________________________________________________________________________________________________________

void HistoPlot::style_ratio_histo(TH1F* single_histo, const char* x_label)
{
  single_histo->GetYaxis()->SetTitle("Data/Background"); //when not tmva was data/MC
  single_histo->GetYaxis()->SetLabelSize(0.12);
  single_histo->GetYaxis()->SetTitleOffset(0.45);
  single_histo->GetYaxis()->SetTitleSize(0.12);
  single_histo->GetXaxis()->SetLabelSize(0.12);
  single_histo->GetXaxis()->SetTitle(x_label);
  single_histo->GetXaxis()->SetTitleSize(0.12);
  single_histo->GetXaxis()->SetTitleOffset(1.1);

  single_histo->SetTitle("");
  single_histo->SetStats(false);
  single_histo->GetYaxis()->SetNdivisions(5, 5, 0);
}
//_______________________________________________________________________________________________________________________

void HistoPlot::style_legend(TLegend* legend)
{
  legend->SetTextSize(0.025);
  legend->SetBorderSize(0);
  legend->SetFillStyle(0);
}
//_______________________________________________________________________________________________________________________

TH1F* HistoPlot::build_1d_histo(DataChain* data_chain, Variable* variable, bool with_cut, bool is_signal,
                                const char* option, std::vector<Variable*>* variables, std::string selection,
				double mc_weight, std::string mva_cut)
{
  std::string var_arg = variable->build_var_string(data_chain->label, with_cut);

  std::string selection_str;

  if (selection == "")
  {
    selection_str = get_selection(variable, variables, with_cut, is_signal, data_chain, mc_weight, mva_cut);
  }
  else
  {
    selection_str = selection;
  }
//cout<<"var string in build 1d histo: "<<var_arg<<endl;
//cout<<"selection string in build 1d histo: "<<selection_str<<"\n";
  data_chain->chain->Draw(var_arg.c_str(), selection_str.c_str(), option);

  TH1F* histo = (TH1F*)gDirectory->Get(data_chain->label);

  return histo;
}
//_______________________________________________________________________________________________________________________

TH1F* HistoPlot::build_parked_histo(DataChain* data_chain, Variable* variable,std::vector<Variable*>* variables,double mc_weight)
{
  std::string var_arg = variable->build_var_string(data_chain->label, true);

  std::string selection_str;

  selection_str = get_parked_selection(variable, variables, data_chain, mc_weight);
  
//cout<<"selection string in build parked histo: "<<selection_str<<"\n";
  data_chain->chain->Draw(var_arg.c_str(), selection_str.c_str());

  TH1F* histo = (TH1F*)gDirectory->Get(data_chain->label);

  return histo;
}
//_______________________________________________________________________________________________________________________


TH1F* HistoPlot::draw_data(DataChain* data_chain, Variable* variable, bool with_cut, TLegend* legend,
                           std::vector<Variable*>* variables, std::string mva_cut, std::string selection)
{
  data_chain->chain->SetMarkerStyle(7);
  data_chain->chain->SetMarkerColor(1);
  data_chain->chain->SetLineColor(1);
  TH1F* data_histo = set_error_bars(
  build_1d_histo(data_chain, variable, with_cut, false, "E1", variables, selection, 1,mva_cut)
  //build_parked_histo(data_chain, variable, variables, 1)
);
  legend->AddEntry(data_histo, data_chain->legend, "lep");

  return data_histo;
}
//_______________________________________________________________________________________________________________________

TH1F* HistoPlot::draw_signal(DataChain* data_chain, Variable* variable, bool with_cut, TLegend* legend,
                             std::vector<Variable*>* variables, std::string mva_cut)
{
  data_chain->chain->SetLineColor(2);
  data_chain->chain->SetLineWidth(3);
  data_chain->chain->SetFillColor(0);
  TH1F* signal_histo = build_1d_histo(data_chain, variable, with_cut, true, "goff", variables, "", 1, mva_cut);
  //TH1F* signal_histo = build_parked_histo(data_chain, variable, variables,1);
  legend->AddEntry(signal_histo, (build_signal_leg_entry(variable, data_chain)).c_str(), "l");

  return signal_histo;
}
//_______________________________________________________________________________________________________________________

TH1F* HistoPlot::draw_background(DataChain* data_chain, Variable* variable, int fill_colour, bool with_cut,
																																	std::vector<Variable*>* variables, double mc_weight, std::string mva_cut)
{
  data_chain->chain->SetLineColor(1);
  data_chain->chain->SetFillColor(fill_colour);
  //std::cout << "in draw background: " << mva_cut << std::endl;
  return build_1d_histo(data_chain, variable, with_cut, false, "goff", variables, "", mc_weight, mva_cut);
  //return build_parked_histo(data_chain, variable, variables, mc_weight);
}
//_______________________________________________________________________________________________________________________

TH1F* HistoPlot::draw_background_from_trees(DataChain* data_chain, Variable* variable, int fill_colour, string selection,
 double mc_weight, std::string mva_cut)
{
  data_chain->chain->SetLineColor(1);
  data_chain->chain->SetFillColor(fill_colour);
  //std::cout << "in draw background: " << mva_cut << std::endl;
  selection = add_mc_to_selection(data_chain, variable, selection, mc_weight);
  //std::cout << "final selection in draw background from trees " << selection << std::endl;

  return build_1d_histo(data_chain, variable, true, false, "goff", NULL,selection, mc_weight, mva_cut);
}
//_______________________________________________________________________________________________________________________


TH1F* HistoPlot::data_to_bg_ratio_histo(TH1F* data_histo, TH1F* bg_histo)
{
  TH1F* ratio_histo = (TH1F*) data_histo->Clone();
  ratio_histo->Divide(bg_histo);
  ratio_histo->SetMarkerColor(1);

  return set_ratio_error_bars(ratio_histo, data_histo, bg_histo);
}
//_______________________________________________________________________________________________________________________

TH1F* HistoPlot::set_ratio_error_bars(TH1F* ratio_histo, TH1F* data_histo, TH1F* bg_histo)
{
  int nbins = ratio_histo->GetNbinsX();

  for(int i = 0; i < nbins; i++) 
  {
    double error_val = std::pow(get_data_error(data_histo, i), 0.5) / get_data_error(bg_histo, i);
    ratio_histo->SetBinError(i, error_val);
  }

  return ratio_histo;
}
//_______________________________________________________________________________________________________________________

TH1F* HistoPlot::set_error_bars(TH1F* histo) 
{
  int nbins = histo->GetNbinsX();
  
  for(int i = 0; i < nbins; i++) {
    double error_val = std::pow(get_data_error(histo, i), 0.5);
    histo->SetBinError(i, error_val);
  }

  return histo;
}
//_______________________________________________________________________________________________________________________

void HistoPlot::set_th1d_error_bars(TH1D* histo)
{
  int nbins = histo->GetNbinsX();

  for(int i = 0; i < nbins; i++) {
    double error_val = std::pow(get_th1d_data_error(histo, i), 0.5);
    histo->SetBinError(i, error_val);
  }
}
//_______________________________________________________________________________________________________________________

float HistoPlot::get_data_error(TH1F* histo, int bin) 
{
  return histo->Integral(bin, bin + 1);
}

double HistoPlot::get_th1d_data_error(TH1D* histo, int bin)
{
  return histo->Integral(bin, bin + 1);
}

std::string HistoPlot::build_file_name(Variable* variable, bool with_cut) 
{
  std::string file_name(variable->name);

  if (with_cut)
  {
    file_name += "_cut_";
    file_name.append((variable->bins_cut).c_str());
    file_name += "_";
    file_name.append(variable->x_min_cut);
    file_name += "_";
    file_name.append(variable->x_max_cut);
  }
  else
  {
    file_name += "_nocut_";
    file_name.append(variable->bins_nocut);
    file_name += "_";
    file_name.append(variable->x_min_nocut);
    file_name += "_";
    file_name.append(variable->x_max_nocut);
  }
  file_name += ".png";

  return file_name;
}
//_______________________________________________________________________________________________________________________

std::string HistoPlot::build_signal_leg_entry(Variable* var, DataChain* signal_chain)
{
  std::string signal_leg_str(signal_chain->legend);
  signal_leg_str += " (x";
  signal_leg_str.append(var->signal_multiplier);
  signal_leg_str += ")";

  return signal_leg_str;
}
//_______________________________________________________________________________________________________________________
void HistoPlot::plot_control(bool not_tau, Variable* mva_output, DataChain* data, std::vector<DataChain*> bg_chains,
                                      std::vector<Variable*>* variables, std::string file_name, 
                                      std::string control, std::string mva_cut)
{
  //////////////////////////////////////////////////////////////////
  cout<<"---------------------Started HistoPlot::plot_control --------------------------------"<<endl;
  //step 1: initialise TCanvas
  TCanvas* c1     = new TCanvas("c1", mva_output->name_styled, 800, 800);
  TPad* p1        = new TPad("p1", "p1", 0.0, 0.95, 1.0, 1.0);
  TPad* p2        = new TPad("p2", "p2", 0.0, 0.2, 1.0, 0.95);
  TPad* p3        = new TPad("p3", "p3", 0.0, 0.0, 1.0, 0.2);
  TLegend* legend = new TLegend(0.0, 0.5, 0.0, 0.88);

  p1->SetLeftMargin(0.102);
  p2->SetBottomMargin(0.012);
  p2->SetLeftMargin(0.105);
  p3->SetBottomMargin(0.3);
  p3->SetLeftMargin(0.102);

  p1->Draw();
  p2->Draw();
  p3->Draw();
  p2->cd();
  std::cout << "mva cut: " << mva_cut << std::endl;
   //step 1.2 clone data chain

  string selection;
  if(not_tau)
  { 
  	selection = "((jet1_pt>50.0)&&(jet2_pt>45.0)&&(metnomu_significance>3.5)&&(dijet_deta>4.2)&&(alljetsmetnomu_mindphi>2.0)&&" + control + ")*total_weight_lepveto";
  }
  else
  {
	selection = "((jet1_pt>50.0)&&(jet2_pt>45.0)&&(metnomu_significance>3.5)&&(dijet_deta>4.2)&&(jetmetnomu_mindphi>1.0)&&(mht>20)&&(jet_csv2<0.2)&&" + control + ")*total_weight_lepveto";
  }

  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //step 2: create background histo
  //step 2.1: get zjets mc weight
  double mc_weights_arr[] = { 0.22, 0.73, 0.32, 1.24, 1, 1,1.70, 1};
  std::vector<double> mc_weights_vector (mc_weights_arr, mc_weights_arr + sizeof(mc_weights_arr) / sizeof(mc_weights_arr[0]));  


   THStack stack = draw_stacked_control(legend, mva_output, bg_chains, true, variables, data, mc_weights_vector, mva_cut,
			selection);

    std::cout << "step 2 done" << std::endl;


  ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //step 3: create data histo 
  
  TH1F* data_histo = draw_data(data, mva_output, false, legend, variables, "", selection);
  std::cout<<"final sig_selection: "<< selection<<"\n";

  std::cout << "step 3 done" << std::endl;

  /////////////////////////////////////////////////////////
  //step 4: draw and style primary histogram
  //step 4.1 draw background histogram
  stack.Draw();
  data_histo->Draw("same");

  //step 4.3 style histogram
  style_stacked_histo(&stack, mva_output->name_styled);
  std::cout << "step 4 done" << std::endl;

  TH1F* plot_histos[3] = {(TH1F*)(stack.GetStack()->Last()), data_histo};

  std::vector<TH1F*> plot_histos_vector (plot_histos, plot_histos + sizeof(plot_histos) / sizeof(plot_histos[0]));
  TH1F* max_histo      = get_max_histo(plot_histos_vector);
  stack.SetMaximum(get_histo_y_max(max_histo)*1.15);

  ////////////////////////////////////////////////////////////
  //step 5: build legend and draw substitle
  //step 5.1 build legend
  build_legend(legend, data_histo, mva_output, true);

  //step 5.2 draw subtitles
  //draw_subtitle(mva_output, variables, true, data, "", mva_cut);
  std::cout << "step 5 done" << std::endl;


  /////////////////////////////////////////////////////////////////////////////
  //step 6: draw signal/background secondary histogram
  //step 6.1 change pad then create signal/background histogram
  p3->cd();
  TH1F* signal_bg_ratio_histo = data_to_bg_ratio_histo(data_histo, (TH1F*)(stack.GetStack()->Last()));

  //step 6.2 draw signal/background histo
  signal_bg_ratio_histo->Draw("e1");

  //step 6.3 style signal/background histo
  style_ratio_histo(signal_bg_ratio_histo, mva_output->name_styled);
  std::cout << "ratio histo done" << std::endl;
  //step 6.4 draw line, y=1 on signal/background histogram
  draw_yline_on_plot(mva_output, true, 1.0);
  draw_yline_on_plot(mva_output, true, 1.0);

  std::cout << "step 6 done" << std::endl;

  ////////////////////////////////////////////////
  //step 7: add title and save histogram to file
  //step 7.1 create file name
  std::string img_name = file_name;
  std::cout << "file name" << std::endl;
  
  //step 7.2 change pad and add title
  p1->cd();
  draw_title(mva_output->name_styled);
  //
  //step 7.3 close and save the file
  c1->SaveAs(img_name.c_str());
cout<<"file saved\n";
  c1->Close();
  std::cout << "step 7 done" << std::endl;
cout<<"---------------------Finished HistoPlot::plot_control --------------------------------"<<endl;
}
//_______________________________________________________________________________________________________________________
THStack HistoPlot::draw_stacked_control(TLegend* legend, Variable* var, std::vector<DataChain*> bg_chains,
                                      bool with_cut, std::vector<Variable*>* variables, DataChain* data, std::vector<double> mc_weights_vector,
				 std::string mva_cut, string selection)
{
  THStack stack(var->name_styled, "");
  
  //std::cout << "selection in draw stacked control histo " << selection << std::endl;
  for(int i = 0; i < 8; i++) {
      //cout<<bg_chains[i]->label<<" === background, mc weight: "<<mc_weights_vector[i]<<endl;
      TH1F* single_bg_histo = draw_background_from_trees(bg_chains[i], var, colours()[i], selection, mc_weights_vector[i], mva_cut);
      stack.Add(single_bg_histo);
      std::string legend_str(bg_chains[i]->legend);
      legend_str += (" #font[12]{(MC weight: " + get_string_from_double(mc_weights_vector[i]) + ")}");
      legend->AddEntry(single_bg_histo, legend_str.c_str(), "f");
      //cout<<"drew background: "<<bg_chains[i]->label<<endl;
  }
  return stack;
}
//_______________________________________________________________________________________________________________________

double HistoPlot::get_efficiency_factor(DataChain* ewk_chain,DataChain* qcd_chain, Variable* var, std::vector<Variable*>* variables,std::string mva_cut)
{
	double e_S,e_C;

	e_S = get_e_S(ewk_chain,qcd_chain, var, variables, mva_cut);
	e_C = get_e_C(ewk_chain,qcd_chain, var, variables, mva_cut);

	//cout<<e_S/e_C<<"\n";
	return e_S/e_C;
}
//_______________________________________________________________________________________________________________________

double HistoPlot::get_e_S(DataChain* ewk_chain, DataChain* qcd_chain,Variable* var, std::vector<Variable*>* variables,std::string mva_cut)
{
	double X_Zvv_ewk = 1380; //cross section in pb
	double X_Zvv_qcd = 6600;

	double R_S_ewk = get_R_S_ewk(ewk_chain, var, variables, mva_cut);
	double R_S_qcd = get_R_S_qcd(qcd_chain, var, variables, mva_cut);

	//cout<<"e_S:"<<(R_S_ewk*X_Zvv_ewk + R_S_qcd*X_Zvv_qcd)/(X_Zvv_ewk+X_Zvv_qcd)<<"\n";
	return (R_S_ewk*X_Zvv_ewk + R_S_qcd*X_Zvv_qcd)/(X_Zvv_ewk+X_Zvv_qcd);

}
//_______________________________________________________________________________________________________________________

double HistoPlot::get_e_C(DataChain* ewk_chain, DataChain* qcd_chain,Variable* var, std::vector<Variable*>* variables,std::string mva_cut)
{
	double X_Zmm_ewk = 0.303; //cross section in pb there ius msitake in the CMS note
	double X_Zmm_qcd = 1168;

	double R_C_ewk = get_R_C_ewk(ewk_chain, var, variables, mva_cut);
	double R_C_qcd = get_R_C_qcd(qcd_chain, var, variables, mva_cut);

	//cout<<"e_R: "<<(R_C_ewk*X_Zmm_ewk + R_C_qcd*X_Zmm_qcd)/(X_Zmm_ewk+X_Zmm_qcd)<<"\n";
	return (R_C_ewk*X_Zmm_ewk + R_C_qcd*X_Zmm_qcd)/(X_Zmm_ewk+X_Zmm_qcd);
}
//_______________________________________________________________________________________________________________________

double HistoPlot::get_R_S_ewk(DataChain* ewk_chain, Variable* var, std::vector<Variable*>* variables, std::string mva_cut)
{
	double N_S_ewk; // the number of Z ll EWK events in signal region
	double N_g_ewk_Zm = 4226.5; //number of generated EWK z to ll events with a z mass cut applied
	bool with_cut = true;

    	std::string selection ="";
 	if(variables ==NULL)
	{
		selection = "((alljetsmetnomu_mindphi>2.0)&&(nvetomuons==0)&&(nvetoelectrons==0))*total_weight_lepveto";
    		selection = HistoPlot::add_classID_to_selection(selection, false);
	}
	else
	{
		selection = "((nvetomuons==0)&&(nvetoelectrons==0)&&(jet2_pt>45.0)&&(jet1_pt>50.0)&&(dijet_deta>4.2)&&(metnomu_significance>3.5)&&(alljetsmetnomu_mindphi>2.0))*total_weight_lepveto*1";
	}
    	selection = HistoPlot::add_mva_cut_to_selection(selection, mva_cut);

	//cout<<"selection in get _r_S_ewk: "<<selection<<"\n";
    	TH1F* histo = HistoPlot::build_1d_histo(ewk_chain, var, with_cut, false, "goff", variables, selection,1 , mva_cut);
    	N_S_ewk = histo->Integral();

	//cout<<"N_S_ewk: "<<N_S_ewk<<"\n"; 
	//cout<<"R_S_ewk: "<<N_S_ewk/N_g_ewk_Zm<<"\n";
	return N_S_ewk/N_g_ewk_Zm;
}
//_______________________________________________________________________________________________________________________

double HistoPlot::get_R_S_qcd(DataChain* qcd_chain, Variable* var, std::vector<Variable*>* variables, std::string mva_cut)
{
	double N_g_qcd_Zm = 20334900; //number of generated QCD z to ll events with a z mass cut applied
	double N_S_qcd;  // the number of Z ll QCD events in signal region
	bool with_cut = true;

    	std::string selection ="";
 	if(variables ==NULL)
	{
		selection = "((alljetsmetnomu_mindphi>2.0)&&(nvetomuons==0)&&(nvetoelectrons==0))*total_weight_lepveto";
    		selection = HistoPlot::add_classID_to_selection(selection, false);
	}
	else
	{
		selection = "((nvetomuons==0)&&(nvetoelectrons==0)&&(jet2_pt>45.0)&&(jet1_pt>50.0)&&(dijet_deta>4.2)&&(metnomu_significance>3.5)&&(alljetsmetnomu_mindphi>2.0))*total_weight_lepveto";
	}
    	selection = HistoPlot::add_mva_cut_to_selection(selection, mva_cut);

	//cout<<"selection in get _r_S_qcd: "<<selection<<"\n";
    	TH1F* histo = HistoPlot::build_1d_histo(qcd_chain, var, with_cut, false, "goff", variables, selection,1 , mva_cut);
    	N_S_qcd = histo->Integral();

	//cout<<"N_S_qcd: "<<N_S_qcd<<"\n";
	//cout<<"R_S_qcd: "<<N_S_qcd/N_g_qcd_Zm<<"\n";
	return N_S_qcd/N_g_qcd_Zm;
}
//_______________________________________________________________________________________________________________________

double HistoPlot::get_R_C_ewk(DataChain* ewk_chain, Variable* var, std::vector<Variable*>* variables, std::string mva_cut)
{
	double N_g_ewk = 5781.9; //number of generated EWK z to ll events
	double N_C_ewk; //the number of Z ll EWK events in control region
	bool with_cut = true;

    	std::string selection ="";
 	if(variables ==NULL)
	{
		selection = "((alljetsmetnomu_mindphi>2.0)&&(nselmuons == 2)&&(m_mumu>60)&&(m_mumu<120))*total_weight_lepveto";
    		selection = HistoPlot::add_classID_to_selection(selection, false);
	}
	else
	{
		selection ="((nvetoelectrons==0)&&(ntaus == 0)&&(nselmuons == 2)&&(m_mumu>60)&&(m_mumu<120)&&(jet2_pt>45.0)&&(jet1_pt>50.0)&&(dijet_deta>4.2)&&(metnomu_significance>3.5)&&(alljetsmetnomu_mindphi>2.0))*total_weight_lepveto*1";//(&&((jet1_eta<0)&&(jet2_eta>0)||(jet1_eta>0)&&(jet2_eta<0))&&((metnomuons/metnomu_significance)>4)&&(metnomu_significance>4)&&(jetmetnomu_mindphi>1.0)&&(mht>20)&&(dijet_M>1200.0)&&(jet2_eta>-4.7)&&(jet2_eta<4.7)&&(jet1_eta>-4.7)&&(jet1_eta<4.7)&&(metnomuons>90)&&(dijet_deta>3.6)&&(jet2_pt>45.0)&&(jet1_pt>50.0))*total_weight_lepveto  
	}
	selection = HistoPlot::add_mva_cut_to_selection(selection, mva_cut);
    	TH1F* histo = HistoPlot::build_1d_histo(ewk_chain, var, with_cut, false, "goff", variables, selection,1 , mva_cut);

	N_C_ewk = histo->Integral();
//cout<<"N_C_ewk:"<<N_C_ewk<<"\n";
//	cout<<"R_C_ewk: "<<N_C_ewk/N_g_ewk<<"\n";
	return N_C_ewk/N_g_ewk;
}
//_______________________________________________________________________________________________________________________

double HistoPlot::get_R_C_qcd(DataChain* qcd_chain, Variable* var, std::vector<Variable*>* variables, std::string mva_cut)
{
	double N_g_qcd = 22789300; //number of generated QCD z to ll events
	double N_C_qcd;// the number of Z ll QCD events in control region
	bool with_cut = true;

    	std::string selection ="";
 	if(variables ==NULL)
	{
		selection = "((alljetsmetnomu_mindphi>2.0)&&(nselmuons == 2)&&(m_mumu>60)&&(m_mumu<120))*total_weight_lepveto";
    		selection = HistoPlot::add_classID_to_selection(selection, false);
	}
	else
	{
		selection ="((nvetoelectrons==0)&&(ntaus == 0)&&(nselmuons == 2)&&(m_mumu>60)&&(m_mumu<120)&&(jet2_pt>45.0)&&(jet1_pt>50.0)&&(dijet_deta>4.2)&&(metnomu_significance>3.5)&&(alljetsmetnomu_mindphi>2.0))*total_weight_lepveto*1";
	}
	selection = HistoPlot::add_mva_cut_to_selection(selection, mva_cut);
    	TH1F* histo = HistoPlot::build_1d_histo(qcd_chain, var, with_cut, false, "goff", variables, selection,1 , mva_cut);

	N_C_qcd = histo->Integral();

	//cout<<"N_C_qcd: "<<N_C_qcd<<"\n";
	//cout<<"R_C_qcd: "<<N_C_qcd/N_g_qcd<<"\n";
	return N_C_qcd/N_g_qcd;
}
//_______________________________________________________________________________________________________________________

double HistoPlot::get_error_on_efficiency_factor(DataChain* ewk_chain,DataChain* qcd_chain, Variable* var, std::vector<Variable*>* variables,std::string mva_cut)
{
	double e_S,e_C,error_e_S,error_e_C,erp1,erp2,error_on_efficiency_factor;

	e_S = get_e_S(ewk_chain,qcd_chain, var, variables, mva_cut);
	e_C = get_e_C(ewk_chain,qcd_chain, var, variables, mva_cut);
	error_e_S = get_error_on_e_S(ewk_chain,qcd_chain, var, variables, mva_cut);
	error_e_C = get_error_on_e_C(ewk_chain,qcd_chain, var, variables, mva_cut);
	erp1 = pow((error_e_S/e_C),2);
	erp2 = pow((error_e_C*e_S/pow(e_C,2)),2);

	error_on_efficiency_factor = pow((erp1+erp2),0.5);
	
	return error_on_efficiency_factor;
}
//_______________________________________________________________________________________________________________________

double HistoPlot::get_error_on_e_S(DataChain* ewk_chain, DataChain* qcd_chain,Variable* var, std::vector<Variable*>* variables,std::string mva_cut)
{
	double X_Zvv_ewk = 1380; //cross section in pb
	double X_Zvv_qcd = 6600;
	double N_g_ewk_Zm = 4226.5;
	double N_g_qcd_Zm = 20334900;
	
	double R_S_ewk = get_R_S_ewk(ewk_chain, var, variables, mva_cut);
	double R_S_qcd = get_R_S_qcd(qcd_chain, var, variables, mva_cut);
	double errorp1 = R_S_ewk*pow((X_Zvv_ewk/N_g_ewk_Zm),2);
	double errorp2 = R_S_qcd*pow((X_Zvv_qcd/N_g_qcd_Zm),2);
	
	double error_on_e_S = pow((errorp1 + errorp2),0.5)/(X_Zvv_ewk+X_Zvv_qcd);
	
	return error_on_e_S;

}
//_______________________________________________________________________________________________________________________

double HistoPlot::get_error_on_e_C(DataChain* ewk_chain, DataChain* qcd_chain,Variable* var, std::vector<Variable*>* variables,std::string mva_cut)
{
	double X_Zmm_ewk = 0.303; //cross section in pb there ius msitake in the CMS note
	double X_Zmm_qcd = 1168;
	double N_g_ewk = 5781.9;
	double N_g_qcd = 22789300; 

	double R_C_ewk = get_R_C_ewk(ewk_chain, var, variables, mva_cut);
	double R_C_qcd = get_R_C_qcd(qcd_chain, var, variables, mva_cut);
	double errorp1 = R_C_ewk*pow((X_Zmm_ewk/N_g_ewk),2);
	double errorp2 = R_C_qcd*pow((X_Zmm_qcd/N_g_qcd),2);

	double error_on_e_C = pow((errorp1 + errorp2),0.5)/(X_Zmm_ewk+X_Zmm_qcd);

	return  error_on_e_C;
}
//_______________________________________________________________________________________________________________________


