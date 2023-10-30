//////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                 	//
//    This code creates workspaces needed to be used as input for parametric statistical analysis  	//
//                                                                                                 	//    
//    For any questions contact:																	//
//			    Sandra Consuegra Rodríguez (sandra.consuegra.rodriguez@desy.de)						//
//     			Daina Leyva Pernia (daina.leyva.pernia@desy.de)										//
//				Chayanit Asawatangtrakuldee (chayanit.asawatangtrakuldee@cern.ch)					//
//    for full Run2 Legacy Analysis see https://github.com/consuegs/MSSM_Htobb-combine				//
//																									//
//////////////////////////////////////////////////////////////////////////////////////////////////////

#include <iostream>
#include <cstring>
#include <string>
#include <cmath>
#include "TH1.h"
#include "TStyle.h"
#include "TMath.h"
#include <iostream>
#include "RooWorkspace.h"
using namespace std;
using namespace RooFit;

int AnalysisWorkspaceSR2()
{

	std::ofstream textout("figs/AnalysisWorkspaceSR2.txt");
	TString dir("/afs/desy.de/user/l/leyvaped/public/rootfiles_july2023/");
	TString dir_param("/nfs/dust/cms/user/leyvaped/Analyses/MSSM/FullRun2/Combine/April_2023/CMSSW_11_3_4/src/Analysis/Combine/Run2018/input_doubleCB/");

	int rebin = 1;

	// As usual, load the combine library to get access to the RooParametricHist
	gSystem->Load("libHiggsAnalysisCombinedLimit.so");
	
	vector<double> lumiscalefactors = { 28.75, 23.26, 29.35, 36.57 };	//SR2
	vector<double> prefiringweights = { 0.99939, 0.99932, 0.99924, 0.99914};	//from prefiring weigth studies, saved in csv format in github, e.g.: https://github.com/desy-cms/analysis-calibrations/blob/master/2018/extras/2018_Full_Hadronic_L1Prefiring.csv
	vector<string> srmasses = { "400", "450", "500", "600" };	//SR2
	TString Tsrmasses[4] 	= { "400", "450", "500", "600" };	//SR2
	vector<float> JER_sigma_shift = { 0.00144, 0.00173, 0.00135, 0.00398 }; //JER nuisance shift (sigma)
	vector<float> JER_norm_shift  = { 0.00770, 0.00571, 0.00426, 0.00250 }; //JES nuisance shift (norm)
	vector<float> JES_mean_shift  = { 0.00891, 0.00933, 0.00946, 0.00893 }; //JES nuisance shift (mean)
	vector<float> JES_norm_shift  = { 0.03589, 0.02645, 0.01981, 0.01299 }; //JES nuisance shift (norm)
	
	//some warnings
	if (!(lumiscalefactors.size() == srmasses.size())) 	{   cout << "WARNING: Number of mass points and lumi scale factors does not agree. Please check what you provided." << endl; 		return -1; 	}
	if (!(prefiringweights.size() == srmasses.size())) 	{   cout << "WARNING: Number of mass points and prefiringweights does not agree. Please check what you provided." << endl; 		return -1; 	}
	if (!(JER_sigma_shift.size() == srmasses.size())) 	{   cout << "WARNING: Number of mass points and JER_sigma_shift does not agree. Please check what you provided." << endl; 		return -1; 	}
	if (!(JES_mean_shift.size() == srmasses.size())) 	{   cout << "WARNING: Number of mass points and JES_mean_shift does not agree. Please check what you provided." << endl; 		return -1; 	}
	if (!(JES_norm_shift.size() == srmasses.size())) 	{   cout << "WARNING: Number of mass points and JES_norm_shift does not agree. Please check what you provided." << endl; 		return -1; 	}
	if (!(JER_norm_shift.size() == srmasses.size())) 	{   cout << "WARNING: Number of mass points and JER_norm_shift does not agree. Please check what you provided." << endl; 		return -1; 	}

	//maps relevant for signal scaling (lumi and prefiring)
	map<string, double> assignedlumisf;
	map<string, double> assignedprefiringweights;
	
	for (unsigned int massvalue = 0; massvalue < srmasses.size(); massvalue++)
	{
		assignedlumisf[srmasses[massvalue]] = 1. / lumiscalefactors[massvalue];
		assignedprefiringweights[srmasses[massvalue]] = prefiringweights[massvalue];
	}

	// A search in a mbb tail, define mbb as our variable
	int m12_min = 320, m12_max = 800;
	RooRealVar mbb("mbb", "m_{12}", m12_min, m12_max);	//SR 2: 400/450/500/600
	RooArgList vars(mbb);

	//loop over mass-points in the mass range
	for (unsigned int mass = 0; mass < srmasses.size(); mass++)
	{
		cout << endl;
		cout << endl;
		cout << "mass " << srmasses[mass];

		///
		/// GET SIG NORMALIZATION 
		/// 

		TFile *f_signal_in = new TFile(dir + "/mssmHbb_2018_FH_" + Tsrmasses[mass] + "_sr.root", "READ");		
		TH1F *h_signal_in = (TH1F*) f_signal_in->Get("mbb");
		double lumisf = assignedlumisf[srmasses[mass]];
		double prefiringweight = assignedprefiringweights[srmasses[mass]];
		cout << "  lumi sf = " << lumisf;
		cout << "  prefiring weight = " << prefiringweight;
		double normSignal = h_signal_in->GetSum() *lumisf *prefiringweight;
		cout << "  norm signal = " << normSignal << std::endl;
		h_signal_in->Scale(lumisf*prefiringweight);
		RooDataHist sigHist("sigHist", "sigHist", mbb, h_signal_in);

		///
		/// GET DATA_OBS HISTS FOR CR/SR 
		///
		
		TFile *f_cr_in = new TFile(dir + "/mssmHbb_2018_FH_Run2018ABCD_cr.root", "READ");
		TH1F *h_cr_in = (TH1F*) f_cr_in->Get("mbb");
		h_cr_in->SetName("h_cr_in");
		//int normCR = h_cr_in->GetEntries();
		int normCR = h_cr_in->Integral(h_cr_in->FindBin(m12_min),h_cr_in->FindBin(m12_max));
		cout << "normCR: " << normCR << endl;
		RooDataHist RDHCR("RDHCR", "CR", vars, h_cr_in);


		TFile *f_sr_in = new TFile(dir + "/mssmHbb_2018_FH_Run2018ABCD_sr.root", "READ");
		TH1F *SRHist_norm = (TH1F*) f_sr_in->Get("mbb"); //blinded -> used to get normalization in SR
		TH1F *SRHist = (TH1F*) f_cr_in->Get("mbb");	//data_obs CR -> now using the data in CR with normalization from SR
		SRHist->SetName("SRHist");
		double scaleSRtoys = SRHist_norm->GetEntries()/h_cr_in->GetEntries();
		SRHist->Scale(scaleSRtoys);
		//int normSR = SRHist->Integral(SRHist->FindBin(m12_min),SRHist->FindBin(m12_max));
		int normSR = 289419;
		cout << "normSR: " << normSR << endl;
		cout << "estimation from CR and TF: " <<  SRHist->Integral(SRHist->FindBin(m12_min),SRHist->FindBin(m12_max)) << endl;
		cout << "percent. difference: " <<  (normSR - SRHist->Integral(SRHist->FindBin(m12_min),SRHist->FindBin(m12_max)))/normSR*100 << endl;
		RooDataHist RDHSR("RDHSR", "SR", vars, SRHist);


		h_cr_in->Rebin(rebin);
		SRHist->Rebin(rebin);

		///
		/// GET BG PARAMETRIZATION FROM ROOFIT
		///

		TFile *f_bgfit = new TFile(dir + "/workspaces_bkg_CR/4FRs/FR2/320to800/extnovosibirsk/workspace/FitContainer_workspace.root", "READ");
		RooWorkspace *w_bgfit = (RooWorkspace*) f_bgfit->Get("workspace");
		RooAbsPdf *background = w_bgfit->pdf("background");
		RooRealVar background_norm("background_norm", "Number of background events", normCR, 0.9 *normCR, 1.1 *normCR);
		
		// Change the names of the parameters
		w_bgfit->var("peak")->SetName("peak_2018");
		w_bgfit->var("tail")->SetName("tail_2018");
		w_bgfit->var("par4")->SetName("par4_2018");
		w_bgfit->var("width")->SetName("width_2018");

		///
		/// GET SIG PARAMETRIZATION FROM ROOFIT
		///

		TFile *f_signal_in_unbinned = new TFile(dir_param + "signal_m" + Tsrmasses[mass] + "_SR2.root", "READ");
		RooWorkspace *w_signalfit = (RooWorkspace*) f_signal_in_unbinned->Get("w");
		RooAbsPdf *signalx = w_signalfit->pdf("signal_dcb");
		signalx->SetName("signal");

		RooRealVar *mean_ws = (RooRealVar*) w_signalfit->var("mean");
		RooRealVar *sigma_ws = (RooRealVar*) w_signalfit->var("sigma");
		RooRealVar *alpha1_ws = (RooRealVar*) w_signalfit->var("alpha1");
		RooRealVar *alpha2_ws = (RooRealVar*) w_signalfit->var("alpha2");
		RooRealVar *n1_ws = (RooRealVar*) w_signalfit->var("n1");
		RooRealVar *n2_ws = (RooRealVar*) w_signalfit->var("n2");
		mean_ws->setConstant(true);
		sigma_ws->setConstant(true);
		alpha1_ws->setConstant(true);
		alpha2_ws->setConstant(true);
		n1_ws->setConstant(true);
		n2_ws->setConstant(true);
		cout << "mean       = " << mean_ws->getVal() << endl;
		cout << "sigma     = " << sigma_ws->getVal() << endl;
		cout << "alpha1     = " << alpha1_ws->getVal() << endl;
		cout << "alpha2 = " << alpha2_ws->getVal() << endl;
		cout << "n1 = " << n1_ws->getVal() << endl;
		cout << "n2 = " << n2_ws->getVal() << endl;
		
		RooPlot *xframe = mbb.frame();
		sigHist.plotOn(xframe, LineColor(1), MarkerColor(1));
		signalx->plotOn(xframe, LineColor(2));
		TCanvas *c1 = new TCanvas("c1", "c1", 600, 600);
		xframe->Draw();
		c1->Update();
		c1->Print("figs/sig_SR2_" + Tsrmasses[mass] + ".png");
		c1->Print("figs/sig_SR2_" + Tsrmasses[mass] + ".pdf");
		delete c1;

	    ///
		/// SYSTEMATIC VARIATIONS OF THE SIGNAL SHAPE PARAMETERS
		///

		RooRealVar theta_JES("CMS_JES_2018", "CMS_JES_2018", 0., -5., 5.);
		RooRealVar theta_JER("CMS_JER_2018", "CMS_JER_2018", 0., -5., 5.);

		RooRealVar shift_JES_mean("shift_JES_mean", "shift_JES_mean", JES_mean_shift[mass]);
		RooRealVar shift_JES_norm("shift_JES_norm", "shift_JES_norm", JES_norm_shift[mass]);
		RooRealVar shift_JER_Sigma("shift_JER_Sigma", "shift_JER_Sigma", JER_sigma_shift[mass]);
		RooRealVar shift_JER_norm("shift_JER_norm", "shift_JER_norm", JER_norm_shift[mass]);

		RooRealVar mean("mean", "mean", 300, 280, 2000);
		RooRealVar sigma("sigma", "sigma", 30, 20, 90);
		RooRealVar alpha1("alpha1", "alpha1", 0.8, 0.1, 20);
		RooRealVar alpha2("alpha2", "alpha2", 20, 0.5, 35);
		RooRealVar n1("n1", "n1", 130, 1, 200);
		RooRealVar n2("n2", "n2", 15, 1, 25);

		//mean.setVal(Mean);
		//sigma.setVal(Sigma);
		mean.setVal(mean_ws->getVal());
		sigma.setVal(sigma_ws->getVal());
		alpha1.setVal(alpha1_ws->getVal());
		alpha2.setVal(alpha2_ws->getVal());
		n1.setVal(n1_ws->getVal());
		n2.setVal(n2_ws->getVal());

		mean.setConstant(true);
		sigma.setConstant(true);
		alpha1.setConstant(true);
		alpha2.setConstant(true);
		n1.setConstant(true);
		n2.setConstant(true);

        RooFormulaVar mean_shifted("mean_shifted",
			"@0*(1+@1*@2)",
			RooArgList(mean, shift_JES_mean, theta_JES));
		RooFormulaVar sigma_shifted("sigma_shifted",
			"@0*(1+@1*@2)",
			RooArgList(sigma, shift_JER_Sigma, theta_JER));

		double normSR_JES_Up, normSR_JES_Down, normSR_JER_Up, normSR_JER_Down;

		TFile fileJesUp(dir + "/mssmHbb_2018_FH_" + Tsrmasses[mass] + "_JES_1sup.root");
		TH1F *SRHist_JES_Up = (TH1F*) fileJesUp.Get("mbb");
		normSR_JES_Up = SRHist_JES_Up->GetSumOfWeights() *lumisf;

		TFile fileJesDown(dir + "/mssmHbb_2018_FH_" + Tsrmasses[mass] + "_JES_1sdown.root");
		TH1F *SRHist_JES_Down = (TH1F*) fileJesDown.Get("mbb");
		normSR_JES_Down = SRHist_JES_Down->GetSumOfWeights() *lumisf;

		TFile fileJerUp(dir + "/mssmHbb_2018_FH_" + Tsrmasses[mass] + "_JER_1sup.root");
		TH1F *SRHist_JER_Up = (TH1F*) fileJerUp.Get("mbb");
		normSR_JER_Up = SRHist_JER_Up->GetSumOfWeights() *lumisf;

		TFile fileJerDown(dir + "/mssmHbb_2018_FH_" + Tsrmasses[mass] + "_JER_1sdown.root");
		TH1F *SRHist_JER_Down = (TH1F*) fileJerDown.Get("mbb");
		normSR_JER_Down = SRHist_JER_Down->GetSumOfWeights() *lumisf;

		RooDoubleCB signal("signal", "DoubleCB", mbb, mean_shifted, sigma_shifted, alpha1, n1, alpha2, n2);

		double d_norm_JES_Up = normSR_JES_Up - normSignal;
		double d_norm_JES_Down = normSignal - normSR_JES_Down;
		double d_norm_JER_Up = normSR_JER_Up - normSignal;
		double d_norm_JER_Down = normSignal - normSR_JER_Down;

		RooRealVar norm("norm", "norm", normSignal, 0., 2. *normSignal);
		RooRealVar norm_JES_Up("norm_JES_Up", "norm_JES_Up", d_norm_JES_Up, -100., 100.);
		RooRealVar norm_JES_Down("norm_JES_Down", "norm_JES_Down", d_norm_JES_Down, -100., 100.);
		RooRealVar norm_JER_Up("norm_JER_Up", "norm_JER_Up", d_norm_JER_Up, -100., 100.);
		RooRealVar norm_JER_Down("norm_JER_Down", "norm_JER_Down", d_norm_JER_Down, -100., 100.);

		RooFormulaVar signal_norm("signal_norm",
			"@0 + @1*@2 + @3*@4",
			RooArgList(norm, shift_JES_norm, theta_JES, shift_JER_norm, theta_JER));

		norm.setVal(normSignal);
		norm.setConstant(true);

		textout << "mass = " << Tsrmasses[mass] << endl;
		textout << "signal normalization -> " << endl;
		textout << endl;

		///
		/// DEFINE TRANSFER FACTOR PDF
		///		
		
		///
		/// DEFINE TRANSFER FACTOR PDF
		///		

		// pdf index 0 -> power-based linear TF
		double TF_pol1_linear_centralValue =  5.8445e-03; //qcd mc linear function
		RooRealVar TF_pol1_linear_2018("TF_pol1_linear_2018", "TF_pol1_linear_2018", TF_pol1_linear_centralValue, 0, 0.1);
		RooArgList varsTF_pol1(mbb, TF_pol1_linear_2018);
		RooPolynomial TF_pol1("TF_pol1", "TF_pol1", mbb, RooArgList(TF_pol1_linear_2018), 1);
		cout << "TF_pol1_linear_2018     = " << TF_pol1_linear_2018.getVal() << endl;

		//pdf index 1 -> power-based quadratic TF
		double TF_pol2_linear_centralValue = 3.3172e-03 ; //qcd mc quadratic function
		double TF_pol2_quad_centralValue   = 9.5524e-07 ; //qcd mc quadratic function
		RooRealVar TF_pol2_quad_2018("TF_pol2_quad_2018", "TF_pol2_quad_2018", TF_pol2_quad_centralValue, -0.001, 0.001);
		RooRealVar TF_pol2_linear_2018("TF_pol2_linear_2018", "TF_pol2_linear_2018", TF_pol2_linear_centralValue, -0.1, 0.1);
		RooArgList varsTF(mbb, TF_pol2_quad_2018, TF_pol2_linear_2018);
		RooPolynomial TF_pol2("TF_pol2", "TF_pol2", mbb, RooArgList(TF_pol2_linear_2018, TF_pol2_quad_2018), 1);

		cout << "TF_pol2_quad_2018       = " << TF_pol2_quad_2018.getVal() << endl;
		cout << "TF_pol2_linear_2018     = " << TF_pol2_linear_2018.getVal() << endl;
		

		
		//RooRealVar signalregion_norm("signalregion_norm", "Signal normalization", normSR, 0.9 *normSR, 1.1 *normSR);
		RooRealVar alternative0_norm("alternative0_norm", "Signal alternative0_norm", normSR, 0.1 *normSR, 1.9 *normSR);
		RooRealVar alternative1_norm("alternative1_norm", "Signal alternative1_norm", normSR, 0.1 *normSR, 1.9 *normSR);
		RooRealVar multipdf_norm("multipdf_norm", "Signal multipdf_norm", normSR, 0.1 *normSR, 1.9 *normSR);

		//Output file
		TFile *fOut = new TFile("input_2018_FH/signal_workspace_" + Tsrmasses[mass] + "_SR2.root", "RECREATE");
		RooWorkspace wspace("wspace", "wspace");

		wspace.import(RDHCR);
		wspace.import(RDHSR);
		wspace.import(signal);
		wspace.import(*background);
		wspace.import(background_norm);
		wspace.import(signal_norm);
		//wspace.import(TF);

		wspace.import(TF_pol1);
		wspace.import(TF_pol2);

		wspace.factory("PROD::alternative0(background,TF_pol1)"); //index 0
		wspace.import(alternative0_norm);

		wspace.factory("PROD::alternative1(background,TF_pol2)"); //index 1
		wspace.import(alternative1_norm);

		RooAbsPdf *alternative0 = wspace.pdf("alternative0");
		RooAbsPdf *alternative1 = wspace.pdf("alternative1");

		RooCategory cat("pdf_index", "Index of Pdf which is active");
		RooArgList mypdfs;
		mypdfs.add(*alternative0);
		mypdfs.add(*alternative1);

		RooMultiPdf multipdf("multipdf", "all pdfs", cat, mypdfs);

		wspace.import(cat);
		wspace.import(multipdf);
		wspace.import(multipdf_norm);

		wspace.Write();
		cout << "File created: signal_workspace_" + Tsrmasses[mass] + "_SR2.root" << endl;
		fOut->Close();
	}
	return 0;
}
