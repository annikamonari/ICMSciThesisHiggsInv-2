#include "../include/histo_plot.h"
#include <algorithm>
using namespace std;
void HistoPlot::draw_plot(Variable* var, std::vector<DataChain*> bg_chains,
                          DataChain* signal_chain, DataChain* data, bool with_cut,
                          std::vector<Variable*>* variables, bool plot_data, std::string file_name,
																										std::string mva_cut)
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
  std::cout << "mva cut: " << mva_cut << std::endl;
  THStack stack      = draw_stacked_histo(legend, var, bg_chains, with_cut, variables, data, mva_cut);
  std::cout << "drew stack" << std::endl;
  TH1F* signal_histo = draw_signal(signal_chain, var, with_cut, legend, variables, mva_cut);
  std::cout << "drew signal" << std::endl;
  std::cout << signal_histo << std::endl;
  TH1F* data_histo;
  if (plot_data) {data_histo = draw_data(data, var, with_cut, legend, variables, mva_cut);}

  stack.Draw();
  signal_histo->Draw("SAME");

  if (plot_data) {data_histo->Draw("SAME");}

  std::cout << "drew all" << std::endl;
  style_stacked_histo(&stack, var->name_styled);

  /*TH1F* plot_histos[3] = {(TH1F*)(stack.GetStack()->Last()), data_histo, signal_histo};
  if (!plot_data) {plot_histos[1] = NULL;}
  std::vector<TH1F*> plot_histos_vector (plot_histos, plot_histos + sizeof(plot_histos) / sizeof(plot_histos[0]));
  TH1F* max_histo      = get_max_histo(plot_histos_vector);
  std::cout << "got max" << max_histo << std::endl;
  stack.SetMaximum(get_histo_y_max(max_histo)*1.15);*/

  build_legend(legend, signal_histo, var, with_cut);

  draw_subtitle(var, variables, with_cut, data, "", mva_cut);
  std::cout << "drew subtitle" << std::endl;
  p3->cd();
  TH1F* data_bg_ratio_histo;

  if (plot_data)
  {
  		std::cout << "data histo isnt null" << std::endl;
  		data_bg_ratio_histo = data_to_bg_ratio_histo(data_histo, (TH1F*)(stack.GetStack()->Last()));
  }
  else
  {
  		std::cout << "data_histo appaz is null" << std::endl;
  		data_bg_ratio_histo = data_to_bg_ratio_histo(signal_histo, (TH1F*)(stack.GetStack()->Last()));
  }

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
  std::cout << "file name" << std::endl;
  p1->cd();
  draw_title(var->name_styled);
  c1->SaveAs(img_name.c_str());
  c1->Close();
}
//////////////////////////////////////////////////////////////////////
////////////Plots a histogram of the evaluated test data for zjets vv with signal 
/////////// in red overlaid on top of grey background, als includes additional signal/nbackground plot below
///////////
void HistoPlot::plot_evaluated_zjets_vv_testTree(Variable* mva_output, DataChain* testTree_chain, DataChain* data, std::vector<DataChain*> bg_chains, std::vector<Variable*>* variables, std::string file_name, std::string mva_cut)
{
  //////////////////////////////////////////////////////////////////
  cout<<"---------------------Started HistoPlot::plot_evaluated_zjets_vv_testTree --------------------------------"<<endl;
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
	 TChain* test_clone     = (TChain*) testTree_chain->chain->Clone();
 DataChain* clone_chain  = new DataChain(mc_signal_data, mc_signal_label, mc_signal_legend, "","",test_clone);
   std::cout << "step 1 done" << std::endl;


  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //step 2: create background histo
  //step 2.1: get zjets mc weight
  std::vector<double> mc_weights_vector = mc_weights(data, bg_chains, mva_output, true, NULL, mva_cut);

  double zjets_weight = mc_weights_vector[6];

  cout<<"mc_weight = "<<zjets_weight<<endl;
//mc weights have been tested and work, now need to edit the draw stack function to not stack the zjets in the datachain vector and instead add the test_set th1
  //step 2.2 create output selection string
  std::string selection = "((nvetomuons==0)&&(nvetoelectrons==0))*total_weight_lepveto";
  selection = add_classID_to_selection(selection, false);
  //selection  += "*total_weight_lepveto";
  //step 2.3 add mva cut to selection string
  //step 2.3.1 reformat mva_cut string for MLPs

  selection = HistoPlot::add_mva_cut_to_selection(selection, mva_cut);
  selection = add_mc_to_selection(testTree_chain,mva_output , selection, zjets_weight);

  //step 2.4 add mc weight to selection  
  /*std::string mc_weight_str = get_string_from_double(mc_weight);
  selection = selection + "*" + mc_weight_str; */
  std::cout<<"final selection: "<<selection<<"\n";

  //step 2.5: get background histo  problem area caused by signal using the same chain, try cloning the tree and running agsain
  testTree_chain->chain->SetLineColor(1);
  testTree_chain->chain->SetFillColor(15);
  TH1F* zjets_histo = build_1d_histo(testTree_chain, mva_output, true, true, "goff", NULL, selection, zjets_weight, mva_cut);  
  TH1F bg_histo_inThe_memory_Stack = *zjets_histo;// save to stack emory fromn the heap

  //step 2.6 add legend entry
  std::string legend_str("Zjets->vv");
  legend_str += (" #font[12]{(MC weight: " + get_string_from_double(zjets_weight) + ")}");// get string from double fails so mc weight is added manually
  legend->AddEntry(zjets_histo, legend_str.c_str(), "f");

   //cout<<"legend str: "<<legend_str<<"\n";
   //step 2.7 get stack for all non zjets_vv backgrounds here
   THStack stack = draw_stacked_histo_no_zjets(legend, mva_output, bg_chains, true, variables, data, mc_weights_vector, mva_cut);

   //step 2.8 add zjets_vv to stack
   stack.Add(zjets_histo);
   std::cout << "step 2 done" << std::endl;


  ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //step 3: create signal histo 
  //step 3.1 get regular selection string
  string sig_selection = "((nvetomuons==0)&&(nvetoelectrons==0))*total_weight_lepveto";
  sig_selection = add_classID_to_selection(sig_selection, true);
  //step 3.2: add class ID to selection
  sig_selection = HistoPlot::add_mva_cut_to_selection(sig_selection, mva_cut);

std::cout<<"final sig_selection: "<<sig_selection<<"\n";

  //step 3.3 get signal histo 
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

  /////////////////////////////////////////////////////////
  //step 4: draw and style primary histogram
  //step 4.1 draw background histogram
  stack.Draw();
  std::cout << "step 4.1 done" << std::endl;

  //step 4.2 draw signal histogram
  signal_histo->Draw("SAME");

  std::cout << "step 4.2 done" << std::endl;

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
  TH1F* signal_bg_ratio_histo = data_to_bg_ratio_histo(signal_histo, (TH1F*)(stack.GetStack()->Last()));

  //step 6.2 draw signal/background histo
  signal_bg_ratio_histo->Draw("e1");

  //step 6.3 style signal/background histo
  style_ratio_histo(signal_bg_ratio_histo, mva_output->name_styled);
  std::cout << "ratio histo done" << std::endl;
  //step 6.4 draw line, y=1 on signal/background histogram
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
cout<<"---------------------Finished HistoPlot::plot_evaluated_zjets_vv_testTree --------------------------------"<<endl;
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
                                      bool with_cut, std::vector<Variable*>* variables, DataChain* data, std::vector<double> mc_weights_vector, std::string mva_cut)
{
  THStack stack(var->name_styled, "");
  std::string selection = "((nvetomuons==0)&&(nvetoelectrons==0))*total_weight_lepveto";
  selection = add_classID_to_selection(selection, false);
  selection = HistoPlot::add_mva_cut_to_selection(selection, mva_cut);

  for(int i = 0; i < bg_chains.size(); i++) {
   if(strcmp(bg_chains[i]->label, "bg_zjets_vv")) {
       cout<<bg_chains[i]->label<<" background, mc weight: "<<mc_weights_vector[i]<<endl;
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

  selection.insert(selection.find("(") + 1, lep_sel_default());
  selection = HistoPlot::add_mva_cut_to_selection(selection, mva_cut);
  //std::cout << selection << std::endl;

  return add_mc_to_selection(bg_chain, variable, selection, mc_weight);
}

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

std::vector<double> HistoPlot::mc_weights(DataChain* data, std::vector<DataChain*> bg_chains,
                                          Variable* var, bool with_cut, std::vector<Variable*>* variables,
					std::string mva_cut, bool double_zjets)
{//cout<<"in HistoPlot::mc_weights\n";
  double mc_weight[bg_chains.size()];
  double zll_weight;

  for(int i=0; i<bg_chains.size();i++)
  {
    mc_weight[i] = 1;

    if (bg_chains[i]->lep_sel != "")
    {
    		double mc_weight_val = MCWeights::calc_mc_weight(data, bg_chains, bg_chains[i], var, with_cut,
																																																			    variables, mva_cut, double_zjets);
      if (mc_weight_val > 0)
    		{
      		mc_weight[i] = mc_weight_val;
    		}

    		if(!strcmp(bg_chains[i]->label, "bg_zll"))
	     {
	     		zll_weight = mc_weight[i];
	     }
    }

    if (!strcmp(bg_chains[i]->label, "bg_zjets_vv"))
    {
    	 if (zll_weight != 1)
    	 {
    	 		mc_weight[i] = zll_weight*5.651*1.513;
    	 }
    }
    //std::cout<<i<<": "<<mc_weight[i]<<"\n";
  }
  std::vector<double> mc_weights_vector (mc_weight, mc_weight + sizeof(mc_weight) / sizeof(mc_weight[0]));

  return mc_weights_vector;
}

// new function written just like the one above: HistoPlot::mc_weights, which calculates the right error for the bgs without a control
// region (its just sqrt(unweighted mc events in signal) / unweighted mc events in signal)
// note:: changed so that if mc weight is 1 then dont calculate the mc weight error
// TODO make error use sumw2 and integralanderror
std::vector<double> HistoPlot::get_mc_weight_errors(DataChain* data, std::vector<DataChain*> bg_chains, Variable* var, bool with_cut, std::vector<Variable*>* variables, std::vector<double> bg_mc_weights, std::string mva_cut)
{
string selection;
	 double mc_weight_errors[bg_chains.size()];
	 double zll_weight_error;

	 for(int i = 0; i < bg_chains.size();i++)
	 {
		if(variables==NULL)
		{
    			selection = "((nvetomuons==0)&&(nvetoelectrons==0))*total_weight_lepveto";
    			selection = HistoPlot::add_classID_to_selection(selection, false);
    			selection = HistoPlot::add_mc_to_selection(bg_chains[i],var , selection, bg_mc_weights[i]);
    			selection = HistoPlot::add_mva_cut_to_selection(selection, mva_cut);
//std::cout << "in histoplot mc weight errors" << selection << std::endl;
		}
//cout<<"selection :"<<selection<<"\n";
	   TH1F* histo = build_1d_histo(bg_chains[i], var,with_cut, false, "goff", variables, selection, 1, mva_cut);
//cout<<"got histo\n";
	   double integral = get_histo_integral(histo, with_cut, var);
	 		mc_weight_errors[i] = std::pow(integral, 0.5);
//cout<<"got integral: "<<integral<<"\n";

	   if (bg_mc_weights[i] != 1) 
	   {
		if(strcmp(bg_chains[i]->label, "bg_zjets_vv")){
		//cout<<"mc wiehgt in if loop: "<<bg_mc_weights[i]<<endl;
	     mc_weight_errors[i] = single_bg_error(data, bg_chains, bg_chains[i], var, with_cut, variables, bg_mc_weights[i], mva_cut,selection);
//cout<<"got bg error\n";
		    if(!strcmp(bg_chains[i]->label, "bg_zll"))
		    {
		     	zll_weight_error = mc_weight_errors[i];
		    }
	         }
           }

	   if (!strcmp(bg_chains[i]->label, "bg_zjets_vv"))
	   {
	     mc_weight_errors[i] = zll_weight_error * 5.651 * 1.513*2;
             if(variables==NULL){mc_weight_errors[i] = mc_weight_errors[i]*2;}//*2 for tes/training split

	   }
 cout << bg_chains[i]->label << endl;
	 }
	 std::vector<double> mc_weights_vector (mc_weight_errors, mc_weight_errors + sizeof(mc_weight_errors) / sizeof(mc_weight_errors[0]));

	 return mc_weights_vector;
}

// problem: TH1F* bg doesn't plot with the mc weight? in the function above this, whenever we call this we pass through
// the mc weight (see last arg: double weight), so if you realise we need it then just put it onto the end of the build_1d_histo call
double HistoPlot::single_bg_error(DataChain* data, std::vector<DataChain*> bg_chains, DataChain* bg_chain,
                                 Variable* var, bool with_cut, std::vector<Variable*>* variables, double weight,
std::string mva_cut, std::string selection)
{//cout<<"selection: "<<selection<<"\n";
  TH1F* bg = build_1d_histo(bg_chain, var, with_cut, false, "goff", variables, selection, 1, mva_cut);
  //cout << "histo" << bg << endl;
  double MC_N_S = get_histo_integral(bg, with_cut, var);
//cout << MC_N_S << endl;
  double sigma_N = std::pow(MC_N_S, 0.5);
//cout << "got sigma_N" << endl;
  double sigma_w = MCWeights::calc_weight_error(data, bg_chains, bg_chain, var, with_cut, variables, mva_cut);
//cout << "got weight error" << endl;
  double sigma_total_sq = std::pow(sigma_w*MC_N_S,2)+std::pow(sigma_N*weight,2);
  double sigma_total = std::pow(sigma_total_sq,0.5);
  //std::cout << bg_chain->label << " - single bg error: " << sigma_total << std::endl;
  return sigma_total;
}

std::string HistoPlot::get_string_from_double(double num)
{
  std::ostringstream num_ss;
  num_ss << num;
  std::string num_str(num_ss.str());

  return num_str;
}

double HistoPlot::sig_to_bg_ratio(Variable* var, TH1F* bg,
                                  TH1F* signal_histo, bool with_cut)
{
  double bg_integral  = get_histo_integral(bg, with_cut, var);
  double sig_integral = get_histo_integral(signal_histo, with_cut, var);
  float signal_mult   = atof(var->signal_multiplier);
  float sig_to_bg     = sig_integral / bg_integral / signal_mult;

  return sig_to_bg;
}

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
  std::cout << sel << std::endl;
	 std::string selection = "Selection: " + sel;
  std::string l1        = "#font[12]{" + selection.substr(0, 90) + "-}";
  std::string l2        = "#font[12]{" + selection.substr(88, 90) + "-}";
  std::string l3        = "#font[12]{" + selection.substr(178, 88) + "}";

  TPaveText* pts        = new TPaveText(0.1, 1.0, 0.9, 0.9, "blNDC");

  pts->SetBorderSize(0);
  pts->SetFillColor(0);
  pts->AddText(l1.c_str());
  pts->AddText(l2.c_str());
  pts->AddText(l3.c_str());
  pts->SetAllWith(l1.c_str(), "size", 0.03);
  pts->SetAllWith(l2.c_str(), "size", 0.03);
  pts->SetAllWith(l3.c_str(), "size", 0.03);
  pts->Draw();
}


THStack HistoPlot::draw_stacked_histo(TLegend* legend, Variable* var, std::vector<DataChain*> bg_chains,
                                      bool with_cut, std::vector<Variable*>* variables, DataChain* data, std::string mva_cut)
{
  THStack stack(var->name_styled, "");
  std::vector<double> mc_weights_vector = mc_weights(data, bg_chains, var, with_cut, variables, mva_cut);

  for(int i = 0; i < bg_chains.size(); i++) {
    TH1F* single_bg_histo = draw_background(bg_chains[i], var, colours()[i], with_cut, variables, mc_weights_vector[i],
																																												mva_cut);
    stack.Add(single_bg_histo);
    std::string legend_str(bg_chains[i]->legend);
    legend_str += (" #font[12]{(MC weight: " + get_string_from_double(mc_weights_vector[i]) + ")}");
    legend->AddEntry(single_bg_histo, legend_str.c_str(), "f");
  }
  return stack;
}

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

double HistoPlot::get_histo_y_max(TH1F* histo)
{
  return histo->GetBinContent(histo->GetMaximumBin());
}

void HistoPlot::build_legend(TLegend* legend, TH1F* max_histo, Variable* var, bool with_cut)
{
  double x1 = position_legend_x1(max_histo, var, with_cut);
  double x2 = x1 + 0.26;

  style_legend(legend);

  legend->SetX1(x1);
  legend->SetX2(x2);
  legend->Draw();

}

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

double HistoPlot::get_x1_from_bin(double max_bin, double nbins)
{
  return max_bin * 0.8 / nbins + 0.1;
}

void HistoPlot::style_stacked_histo(THStack* hs, const char* x_label)
{
  hs->GetYaxis()->SetTitle("Events");
  hs->GetYaxis()->SetLabelSize(0.035);
  hs->GetYaxis()->SetTitleOffset(1.55);
  hs->GetXaxis()->SetLabelSize(0);
}

void HistoPlot::style_ratio_histo(TH1F* single_histo, const char* x_label)
{
  single_histo->GetYaxis()->SetTitle("Signal/Background"); //when not tmva was data/MC
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

void HistoPlot::style_legend(TLegend* legend)
{
  legend->SetTextSize(0.025);
  legend->SetBorderSize(0);
  legend->SetFillStyle(0);
}

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
//cout<<"selection string in build 1d histo: "<<selection_str<<"\n";
  data_chain->chain->Draw(var_arg.c_str(), selection_str.c_str(), option);

  TH1F* histo = (TH1F*)gDirectory->Get(data_chain->label);

  return histo;
}

TH1F* HistoPlot::draw_data(DataChain* data_chain, Variable* variable, bool with_cut, TLegend* legend,
                           std::vector<Variable*>* variables, std::string mva_cut)
{
  data_chain->chain->SetMarkerStyle(7);
  data_chain->chain->SetMarkerColor(1);
  data_chain->chain->SetLineColor(1);
  TH1F* data_histo = set_error_bars(build_1d_histo(data_chain, variable, with_cut, false, "E1", variables, "", 1,
																																																			mva_cut));
  legend->AddEntry(data_histo, data_chain->legend, "lep");

  return data_histo;
}

TH1F* HistoPlot::draw_signal(DataChain* data_chain, Variable* variable, bool with_cut, TLegend* legend,
                             std::vector<Variable*>* variables, std::string mva_cut)
{
  data_chain->chain->SetLineColor(2);
  data_chain->chain->SetLineWidth(3);
  data_chain->chain->SetFillColor(0);
  TH1F* signal_histo = build_1d_histo(data_chain, variable, with_cut, true, "goff", variables, "", 1, mva_cut);
  legend->AddEntry(signal_histo, (build_signal_leg_entry(variable, data_chain)).c_str(), "l");

  return signal_histo;
}

TH1F* HistoPlot::draw_background(DataChain* data_chain, Variable* variable, int fill_colour, bool with_cut,
																																	std::vector<Variable*>* variables, double mc_weight, std::string mva_cut)
{
  data_chain->chain->SetLineColor(1);
  data_chain->chain->SetFillColor(fill_colour);
  //std::cout << "in draw background: " << mva_cut << std::endl;
  return build_1d_histo(data_chain, variable, with_cut, false, "goff", variables, "", mc_weight, mva_cut);
}

TH1F* HistoPlot::draw_background_from_trees(DataChain* data_chain, Variable* variable, int fill_colour, string selection,
 double mc_weight, std::string mva_cut)
{
  data_chain->chain->SetLineColor(1);
  data_chain->chain->SetFillColor(fill_colour);
  //std::cout << "in draw background: " << mva_cut << std::endl;
  selection = add_mc_to_selection(data_chain, variable, selection, mc_weight);
  return build_1d_histo(data_chain, variable, true, false, "goff", NULL,selection, mc_weight, mva_cut);
}


TH1F* HistoPlot::data_to_bg_ratio_histo(TH1F* data_histo, TH1F* bg_histo)
{
  TH1F* ratio_histo = (TH1F*) data_histo->Clone();
  ratio_histo->Divide(bg_histo);
  ratio_histo->SetMarkerColor(1);

  return set_ratio_error_bars(ratio_histo, data_histo, bg_histo);
}

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

TH1F* HistoPlot::set_error_bars(TH1F* histo) 
{
  int nbins = histo->GetNbinsX();
  
  for(int i = 0; i < nbins; i++) {
    double error_val = std::pow(get_data_error(histo, i), 0.5);
    histo->SetBinError(i, error_val);
  }

  return histo;
}

void HistoPlot::set_th1d_error_bars(TH1D* histo)
{
  int nbins = histo->GetNbinsX();

  for(int i = 0; i < nbins; i++) {
    double error_val = std::pow(get_th1d_data_error(histo, i), 0.5);
    histo->SetBinError(i, error_val);
  }
}

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

std::string HistoPlot::build_signal_leg_entry(Variable* var, DataChain* signal_chain)
{
  std::string signal_leg_str(signal_chain->legend);
  signal_leg_str += " (x";
  signal_leg_str.append(var->signal_multiplier);
  signal_leg_str += ")";

  return signal_leg_str;
}


