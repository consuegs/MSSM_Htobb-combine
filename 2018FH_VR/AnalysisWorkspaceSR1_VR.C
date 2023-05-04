/////////////////////////////////////////////////////////////////////////////////////////
//
//    Original code from Paul Asmuss (https://github.com/pasmuss/analysis-combine)
//    Edited by Sandra Consuegra Rodríguez (sandra.consuegra.rodriguez@desy.de)
//    Edited by Daina Leyva Pernia (daina.leyva.pernia@desy.de)
//
//    This code creates workspaces needed to to be used as input for statistical analysis 
//
/////////////////////////////////////////////////////////////////////////////////////////
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

int AnalysisWorkspaceSR1_VR()
{

	std::ofstream textout("figs/AnalysisWorkspaceSR1.txt");
	TString dir("/afs/desy.de/user/l/leyvaped/public/for_sandra/rootfiles_2018FH_Feb2023/");

	int rebin = 1;

	// As usual, load the combine library to get access to the RooParametricHist
	gSystem->Load("libHiggsAnalysisCombinedLimit.so");

	vector<double> lumiscalefactors = { 29.41, 29.82, 28.75 };	//SR1
	vector<string> srmasses = { "300", "350", "400" };	//SR1

	TString Tsrmasses[3] = { "300", "350", "400" };	//SR1

	if (!(lumiscalefactors.size() == srmasses.size()))
	{
		cout << "Number of mass points and lumi scale factors does not agree. Please check what you provided." << endl;
		return -1;
	}
	map<string, double> assignedlumisf;
	for (unsigned int massvalue = 0; massvalue < srmasses.size(); massvalue++)
	{
		assignedlumisf[srmasses[massvalue]] = 1. / lumiscalefactors[massvalue];
	}

	// A search in a mbb tail, define mbb as our variable
	int m12_min = 270, m12_max = 560;
	RooRealVar mbb("mbb", "m_{12}", m12_min, m12_max);	//SR 1: 300/350/400
	RooArgList vars(mbb);

	for (unsigned int mass = 0; mass < srmasses.size(); mass++)
	{
		cout << endl;
		cout << endl;
		cout << "mass " << srmasses[mass];

		///
		/// GET SIG NORMALIZATION 
		///

		TFile *f_signal_in = new TFile(dir + "/signal_2018FH_nominal/mssmHbb_2018_FH_" + Tsrmasses[mass] + "_sr.root", "READ");		
		TH1F *h_signal_in = (TH1F*) f_signal_in->Get("mbb");
		double lumisf = assignedlumisf[srmasses[mass]];
		cout << "  lumi sf = " << lumisf;
		double normSignal = h_signal_in->GetSum() *lumisf;
		cout << "  norm signal = " << normSignal << std::endl;
		h_signal_in->Scale(lumisf);
		RooDataHist sigHist("sigHist", "sigHist", mbb, h_signal_in);

		///
		/// GET DATA_OBS HISTS FOR CR/SR 
		///
		
		TFile *f_cr_in = new TFile(dir + "/mssmHbb_2018_FH_Run2018ABCD_cr.root", "READ");
		TH1F *h_cr_in = (TH1F*) f_cr_in->Get("mbb");
		h_cr_in->SetName("h_cr_in");
		h_cr_in->Rebin(rebin);
		//int normCR = h_cr_in->GetEntries();
		int normCR = h_cr_in->Integral(h_cr_in->FindBin(m12_min),h_cr_in->FindBin(m12_max));
		cout << "normCR: " << normCR << endl;
		RooDataHist RDHCR("RDHCR", "CR", vars, h_cr_in);

		TFile *f_vr_in = new TFile(dir + "/mssmHbb_2018_FH_Run2018ABCD_vr.root", "READ");
		TH1F *SRHist = (TH1F*) f_vr_in->Get("mbb");	//data_obs VR -> now using the data in VR with normalization from VR
		SRHist->SetName("SRHist");
		SRHist->Rebin(rebin);
		//int normSR = SRHist->GetEntries();
		int normSR = SRHist->Integral(SRHist->FindBin(m12_min),SRHist->FindBin(m12_max));
		cout << "normSR: " << normSR << endl;
		RooDataHist RDHSR("RDHSR", "SR", vars, SRHist);

		///
		/// GET BG PARAMETRIZATION FROM ROOFIT
		///

		TFile *f_bgfit = new TFile(dir + "/workspaces_bkg_CR/4FRs/FR1/270to560/extnovosibirsk/workspace/FitContainer_workspace.root", "READ");
		RooWorkspace *w_bgfit = (RooWorkspace*) f_bgfit->Get("workspace");
		RooAbsPdf *background = w_bgfit->pdf("background");
		RooRealVar background_norm("background_norm", "Number of background events", normCR, 0.1 *normCR, 1.9 *normCR);

		///
		/// GET SIG PARAMETRIZATION FROM ROOFIT
		///

		TFile *f_signal_in_unbinned = new TFile("/nfs/dust/cms/user/leyvaped/Analyses/MSSM/FullRun2/Combine/April_2023/CMSSW_11_3_4/src/Analysis/Combine/Run2018/input_doubleCB_Feb23/signal_m" + Tsrmasses[mass] + "_SR1.root", "READ");
		RooWorkspace *w_signalfit = (RooWorkspace*) f_signal_in_unbinned->Get("w");
		RooAbsPdf *signal = w_signalfit->pdf("signal_dcb");
		signal->SetName("signal");
		RooRealVar signal_norm("signal_norm", "signal_norm", normSignal);

		RooRealVar * mean = (RooRealVar*)w_signalfit->var("mean");
		RooRealVar * sigma = (RooRealVar*)w_signalfit->var("sigma");
		RooRealVar * alpha1 = (RooRealVar*)w_signalfit->var("alpha1");
		RooRealVar * alpha2 = (RooRealVar*)w_signalfit->var("alpha2");
		RooRealVar * n1 = (RooRealVar*)w_signalfit->var("n1");
		RooRealVar * n2 = (RooRealVar*)w_signalfit->var("n2");
		mean->setConstant(true);
		sigma->setConstant(true);
		alpha1->setConstant(true);
		alpha2->setConstant(true);
		n1->setConstant(true);
		n2->setConstant(true);
		cout << "mean       = " << mean->getVal() << endl;
		cout << "sigma     = " << sigma->getVal() << endl;
		cout << "alpha1     = " << alpha1->getVal() << endl;
		cout << "alpha2 = " << alpha2->getVal() << endl;
		cout << "n1 = " << n1->getVal() << endl;
		cout << "n2 = " << n2->getVal() << endl;

		RooPlot *xframe = mbb.frame();
		sigHist.plotOn(xframe, LineColor(1), MarkerColor(1));
		signal->plotOn(xframe, LineColor(2));
		TCanvas *c1 = new TCanvas("c1", "c1", 600, 600);
		xframe->Draw();
		c1->Update();
		c1->Print("figs/sig_SR1_" + Tsrmasses[mass] + ".png");
		c1->Print("figs/sig_SR1_" + Tsrmasses[mass] + ".pdf");
		delete c1;

		///
		/// DEFINE TRANSFER FACTOR PDF
		///
		
		

		double TF_pol2_quad_centralValue = -2.85726e-06; //mc
		double TF_pol2_linear_centralValue = 0.000353419;
		
		RooRealVar TF_pol2_quad("TF_pol2_quad", "TF_pol2_quad", TF_pol2_quad_centralValue, -0.000001, 0.000001);
		RooRealVar TF_pol2_linear("TF_pol2_linear", "TF_pol2_linear", TF_pol2_linear_centralValue, -0.001, 0.001);
		RooArgList varsTF(mbb, TF_pol2_quad, TF_pol2_linear);

		RooGenericPdf transferfactor("transferfactor", "transferfactor", "TF_pol2_quad*(mbb-500)*(mbb-500)+TF_pol2_linear*(mbb-500)+1", varsTF);	// 

		RooRealVar signalregion_norm("signalregion_norm", "Signal normalization", normSR, 0.9 *normSR, 1.1 *normSR);

		cout << "TF_pol2_quad       = " << TF_pol2_quad.getVal() << endl;
		cout << "TF_pol2_linear     = " << TF_pol2_linear.getVal() << endl;

		//Output file
		TFile *fOut = new TFile("input_2018_FH_VR/signal_workspace_" + Tsrmasses[mass] + "_SR2.root", "RECREATE");
		RooWorkspace wspace("wspace", "wspace");

		wspace.import(RDHCR);
		wspace.import(RDHSR);
		wspace.import(*signal);
		wspace.import(*background);
		wspace.import(background_norm);
		wspace.import(signal_norm);
		wspace.import(transferfactor);
		wspace.factory("PROD::signalregion(background,transferfactor)");
		wspace.import(signalregion_norm);
		wspace.Write();
		cout << "File created: signal_workspace_" + Tsrmasses[mass] + "_SR1.root" << endl;
		fOut->Close();
	}
	return 0;
}
