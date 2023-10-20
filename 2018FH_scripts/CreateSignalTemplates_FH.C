#include "CMS_lumi.C"
using namespace RooFit;

// SR1 : 270-560
// SR2 : 320-800
// SR3 : 370-1270
// SR4 : 500-2000

map<int, double> lumi_sf = {
    {300,29.41},
    {350,29.82},
    {400,28.75},
    {450,23.26},
    {500,29.35},
    {600,36.57},
    {700,36.43},
    {800,35.68},
    {900,35.98},
    {1000,35.77},
    {1200,34.61},
    {1400,34.05},
    {1600,33.24},
    {1800,31.45},
};

map<int, double> l1prefiring_sf = {
    {300, 0.99951},
    {350, 0.99944},
    {400, 0.99939},
    {450, 0.99932},
    {500, 0.99924},
    {600, 0.99914},
    {700, 0.99904},
    {800, 0.99894},
    {900, 0.99888},
    {1000, 0.99883},
    {1200, 0.99871},
    {1400, 0.99869},
    {1600, 0.99857},
    {1800, 0.99841},
};

map<int, int> mass_region = {
    {300,1}, 
    {350,1},  
    {400,1},  
//    {400,2},
    {450,2},
    {500,2},
//    {600,2},
    {600,3},
    {700,3},
    {800,3},
    {900,3},
    {1000,3},
//    {1000,4},
    {1200,4},
    {1400,4},
    {1600,4},
    {1800,4},
};

map<int, int> mass_binning = {
    {300,1},
    {350,1},
    {400,1},
    {450,1},
    {500,1},
    {600,1},
    {700,5},
    {800,5},
    {900,5},
    {1000,5},
    {1200,10},
    {1400,10},
    {1600,10},
    {1800,10},
};

map<int, double> mbb_low = {
    {1,270.},
    {2,320.},
//    {3,390.},
    {3,370.},
    {4,500.},
};

map<int, double> mbb_high = {
    {1,560.},
    {2,800.},
    {3,1270.},
    {4,2000.},
};

TString dir("/afs/desy.de/user/l/leyvaped/public/rootfiles_july2023/");

map<TString, TString> histName_suffix = {
		{
		"nominal", "" },
};

TString IntToStr(int i)
{

	char c[10];
	if (i < 10)
		sprintf(c, "%1i", i);
	else if (i < 100)
		sprintf(c, "%2i", i);
	else if (i < 1000)
		sprintf(c, "%3i", i);
	else
		sprintf(c, "%4i", i);

	TString ST(c);

	return ST;
}
// initial fit parameters for double-sided crystall ball
// mass, mean, sigma, alpha1, alpha2, n1, n2
map<int, vector < double>> mass_fit = {
    { 300,{ 290.,30.,1.8,1.3, 12.,1.1}},
    { 350,{ 340.,35.,1.4,1.7,125.,1.0}},
    { 400,{ 380.,37.,0.9,1.6,132.,1.4}},
    { 450,{ 425.,41.,0.8,1.7,136.,1.4}},
    { 500,{ 470.,44.,0.7,1.8,133.,1.4}},
    { 600,{ 570.,48.,0.5,2.0,116.,0.8}},
    { 700,{ 670.,49.,0.5,1.6,118.,2.0}},
    { 800,{ 770.,55.,0.5,1.8,124.,1.9}},
    { 900,{ 860.,60.,0.5,1.8,5.10,1.9}},
    {1000,{ 960.,61.,0.5,1.9,4.20,1.3}},
    {1200,{1165.,66.,0.4,1.9,4.90,2.0}},
    {1400,{1360.,77.,0.4,1.9,5.00,2.2}},
    {1600,{1550.,83.,0.4,1.7,3.00,4.0}},
    {1800,{1750.,83.,0.4,1.7,3.00,4.0}},
};

vector<int> masses = { 300, 350, 400, 450, 500, 600, 700, 800, 900, 1000, 1200, 1400, 1600, 1800 };

std::ofstream outtext("signal_out_FH.txt");

void CreateSignalPDF(int mass,
	int region,
	vector<TString> histNames,
	map<TString, TH1D*> Hists,
	RooWorkspace *w,
	double mbb_min,
	double mbb_max)
{

	TH1::SetDefaultSumw2(true);

	std::map<TString, double> mapMean;
	std::map<TString, double> mapSigma;
	std::map<TString, double> mapAlpha1;
	std::map<TString, double> mapAlpha2;
	std::map<TString, double> mapN1;
	std::map<TString, double> mapN2;
	std::map<TString, double> mapNorm;

	double Mean = mass_fit[mass][0];
	double Sigma = mass_fit[mass][1];
	double Alpha1 = mass_fit[mass][2];
	double Alpha2 = mass_fit[mass][3];
	double N1 = mass_fit[mass][4];
	double N2 = mass_fit[mass][5];

	outtext << "++++++++++++++++++++++++++++++++++" << std::endl;
	outtext << std::endl;
	outtext << "Mass : " << mass << "  SR" << region << std::endl;

	for (auto histName: histNames)
	{

		outtext << std::endl;
		outtext << std::endl;
		outtext << "Histogram name : " << histName << std::endl;
		outtext << std::endl;

		TH1D *hist = Hists[histName];
		if (mass < 700)
		{
			hist->Rebin(5);
		}
		else
		{
			hist->Rebin(10);
		}

		// Bin width
		double nBinsFullRange = hist->GetNbinsX();
		double nBins = 0;
		for (int i = 1; i <= hist->GetNbinsX(); ++i)
		{
			double center = hist->GetBinCenter(i);
			if (center > mbb_min && center < mbb_max)
				++nBins;
		}

		double binWidth = (mbb_max - mbb_min) / nBins;

		for (int iB = 1; iB <= nBinsFullRange; ++iB)
		{
			if (hist->GetBinContent(iB) < 0)
			{
				hist->SetBinContent(iB, 0.);
			}
		}
		double norm = hist->GetSumOfWeights() / lumi_sf[mass] * l1prefiring_sf[mass];

		outtext << "signal normalization = " << norm << endl;
		outtext << endl;

		RooRealVar mbb("mbb", "M_{12}", mbb_min, mbb_max);

		RooRealVar meanx("mean", "Mean", Mean, 0.5 *Mean, 1.5 *Mean);
		RooRealVar sigmax("sigma", "Width", Sigma, 0., 2.0 *Sigma);

		RooRealVar alpha1x("alpha1", "Alpha1", Alpha1, 0, 50 *Alpha1);
		RooRealVar n1x("n1", "N1", N1, 0, 50 *N1);

		RooRealVar alpha2x("alpha2", "Alpha2", Alpha2, 0, 50 *Alpha2);
		RooRealVar n2x("n2", "N2", N2, 0, 50 *N2);

		RooDoubleCB cbx("cb", "CBshape", mbb, meanx, sigmax, alpha1x, n1x, alpha2x, n2x);

		RooDataHist sig("sig", "sig", mbb, hist);   
		w->import(sig);

		RooFitResult *res = cbx.fitTo(sig, Save(), SumW2Error(kTRUE));
		res->Print();
		outtext << "Fit result with DoubleCB[" << mbb_min << "," << mbb_max << "] -> " << std::endl;
		outtext << "Mean   = " << meanx.getVal() << " +/- " << meanx.getError() << std::endl;
		outtext << "Sigma  = " << sigmax.getVal() << " +/- " << sigmax.getError() << std::endl;
		outtext << "Alpha1 = " << alpha1x.getVal() << " +/- " << alpha1x.getError() << std::endl;
		outtext << "N1     = " << n1x.getVal() << " +/- " << n1x.getError() << std::endl;
		outtext << "Alpha2 = " << alpha2x.getVal() << " +/- " << alpha2x.getError() << std::endl;
		outtext << "N2     = " << n2x.getVal() << " +/- " << n2x.getError() << std::endl;

		TString Mass = IntToStr(mass);
		TString Region = IntToStr(region);

		TCanvas *c1 = new TCanvas("c1", "", 700, 700);
		RooPlot *frame1 = mbb.frame(Name("MMSM H->bb analysis"), Title("   "));
		frame1->GetXaxis()->SetTitle("m_{12}[GeV]");
		frame1->GetYaxis()->SetTitle(Form("Events / %.0f GeV", binWidth));
		sig.plotOn(frame1, MarkerSize(1.1));
		cbx.plotOn(frame1, LineColor(kRed), LineWidth(3));
		c1->SetLeftMargin(0.1);
		c1->SetBottomMargin(0.1);
		gPad->SetLeftMargin(0.13);
		frame1->Draw();

		//// Selecting legend names for 	////	
		for (int i = 0; i < frame1->numItems(); i++)
		{
			TString obj_name = frame1->nameOf(i);
			if (obj_name == "") continue;
		}
		TString names[] = { "h_sig",
			"PDF_Norm[x]",
			"" };
		TString signs[] = { "Simulated Signal",
			//  "Fit Error",
			"Global Fit",
		};
		//// Legend	////
		float xpad1Leg;
		if (/*mass == 600 || */ mass == 1600 || mass == 1800)
		{
			xpad1Leg = 0.42;
		}
		else
		{
			xpad1Leg = 0.98;
		}

		int fitRangeMin, fitRangeMax;
		if (mass < 450)
		{
			fitRangeMin = 270;
			fitRangeMax = 560;
		}
		else if (mass >= 450 & mass < 800)
		{
			fitRangeMin = 320;
			fitRangeMax = 800;
		}
		else if (mass >= 800 & mass < 1200)
		{
			fitRangeMin = 390;
			fitRangeMax = 1270;
		}
		else
		{
			fitRangeMin = 500;
			fitRangeMax = 2000;
		}
     
		int nParameters = res->floatParsFinal().getSize();
		double chi = frame1->chiSquare(nParameters);
		double ndof = nBins - nParameters;
		double p_value = TMath::Prob(chi *ndof, ndof);
		cout << "M_{12} = " + Mass + " GeV: " << endl;
		cout << "Bin width: " << binWidth << endl;
		cout << "Number of parameters of double sided crystal ball = " << nParameters << endl;
		cout << "chi^2 / ndof = " << Form("%.2f", chi *ndof) << " / " << nBins - nParameters << " = " << chi << endl;
		cout << "p-value = " << Form("%.2f", p_value) << endl;
		cout << "//////////////////////////////////////////" << endl;

		TPad * pad1;
		pad1 = new TPad("pad1", "", 0, 0.1, 1, 1);
		pad1->SetBottomMargin(0.2);
		pad1->SetRightMargin(0.05);
		pad1->SetLeftMargin(0.16);
		pad1->Draw();
		pad1->cd();
		frame1->GetXaxis()->SetTitleOffset(999);	//Effectively turn off x axis title on main plot
		frame1->GetXaxis()->SetLabelOffset(999);	//Effectively turn off x axis label on main plot
		frame1->GetYaxis()->SetTitleSize(0.038);
		frame1->GetYaxis()->SetTitleOffset(1.6);
		frame1->GetYaxis()->SetLabelSize(0.033);
		frame1->Draw();

		float lumi_2018_FH = 54.5;

		std::string lumistr(Form("%.1f", lumi_2018_FH));
		std::string chi2str(Form("#chi^{2}/ndof = %.2f ", chi));
		std::string probstr(Form("p-value= %.2f", p_value));

		TLatex latex;
		latex.SetTextFont(61);
		latex.SetTextSize(0.06);
		latex.SetTextAlign(11);
		latex.DrawLatexNDC(pad1->GetLeftMargin(), 1.02 - c1->GetTopMargin(),
			std::string("CMS").c_str());
		latex.SetTextFont(52);
		latex.SetTextSize(0.04);
		latex.DrawLatexNDC(pad1->GetLeftMargin() + 0.12, 1.02 - c1->GetTopMargin(),
			std::string("Work in progress").c_str());
		latex.SetTextFont(43);
		latex.SetTextSize(22);
		latex.DrawLatexNDC(0.71, 1.02 - c1->GetTopMargin(),
			(lumistr + std::string(" fb^{-1} (13 TeV)")).c_str());
		latex.SetTextAlign(33);
		latex.SetTextSize(18);
		latex.SetTextColor(kBlack);
		latex.DrawLatexNDC(xpad1Leg - pad1->GetRightMargin(), 0.98 - pad1->GetTopMargin(),
			(std::string("M_{12} = ").c_str() + Mass + std::string(" GeV").c_str()));
		latex.SetTextColor(kBlue + 2);
		latex.DrawLatexNDC(xpad1Leg - pad1->GetRightMargin() + 0.01, 0.93 - pad1->GetTopMargin(),
			chi2str.c_str());
		latex.SetTextColor(kGreen + 2);
		latex.DrawLatexNDC(xpad1Leg - pad1->GetRightMargin(), 0.88 - pad1->GetTopMargin(),
			probstr.c_str());
		latex.SetTextColor(kOrange + 2);
		std::string minstr(Form("%.0d", fitRangeMin));
		std::string maxstr(Form("%.0d", fitRangeMax));
		latex.DrawLatexNDC(xpad1Leg - pad1->GetRightMargin(), 0.83 - pad1->GetTopMargin(),
			(minstr + std::string("<M_{12}<") + maxstr).c_str());

		// Construct a histogram with the pulls of the data w.r.t the curve (('curve' - histogram) / err_histogram)
		RooHist *hpull = frame1->pullHist();
		hpull->SetMarkerSize(0.8);	
		std::unique_ptr<RooPlot> frame2(mbb.frame());
		frame2->addPlotable(hpull, "P");

		c1->cd();
		TPad *pad2 = new TPad("pad2", "", 0, 0.0, 1, 0.265);
		pad2->SetTopMargin(1);
		pad2->SetBottomMargin(0.33);
		pad2->SetLeftMargin(pad1->GetLeftMargin());
		pad2->SetRightMargin(pad1->GetRightMargin());
		pad2->SetGridy();
		pad2->Draw();
		pad2->cd();
		frame2->SetTitle("");
		frame2->GetXaxis()->SetTitleSize(0.15);
		frame2->GetXaxis()->SetTitleOffset(0.9);
		frame2->GetXaxis()->SetLabelSize(0.115);
		frame2->GetXaxis()->SetLabelOffset(0.010);
		frame2->SetYTitle("Pulls");
		frame2->GetYaxis()->CenterTitle(kTRUE);
		frame2->GetYaxis()->SetTitleSize(0.14);
		frame2->GetYaxis()->SetTitleOffset(0.4);
		frame2->GetYaxis()->SetNdivisions(3, 5, 0);
		frame2->GetYaxis()->SetLabelSize(0.115);
		frame2->GetYaxis()->SetLabelOffset(0.011);
		frame2->SetMinimum(-5.);
		frame2->SetMaximum(+5.);
		frame2->Draw();

		c1->Print("Figs/Mass" + Mass + "_SR" + Region + "_" + histName + "_doubleCB.png");
		c1->Print("Figs/Mass" + Mass + "_SR" + Region + "_" + histName + "_doubleCB.pdf");
		delete c1;

		mapMean[histName] = meanx.getVal();
		mapSigma[histName] = sigmax.getVal();
		mapAlpha1[histName] = alpha1x.getVal();
		mapAlpha2[histName] = alpha2x.getVal();
		mapN1[histName] = n1x.getVal();
		mapN2[histName] = n2x.getVal();
		mapNorm[histName] = norm;
	}

	double norm_Val = mapNorm["nominal"];
	double mean_Val = mapMean["nominal"];
	double sigma_Val = mapSigma["nominal"];
	double alpha1_Val = mapAlpha1["nominal"];
	double alpha2_Val = mapAlpha2["nominal"];
	double n1_Val = mapN1["nominal"];
	double n2_Val = mapN2["nominal"];

	RooRealVar mbb("mbb", "M_{12}", mbb_min, mbb_max);
	RooRealVar mean("mean", "mean", 0.5 *mean_Val, 2 *mean_Val);
	RooRealVar sigma("sigma", "sigma", 0.5 *sigma_Val, 2 *sigma_Val);
	RooRealVar alpha1("alpha1", "alpha1", 0.5 *alpha1_Val, 2 *alpha1_Val);
	RooRealVar alpha2("alpha2", "alpha2", 0.5 *alpha2_Val, 2 *alpha2_Val);
	RooRealVar n1("n1", "n1", 0.5 *n1_Val, 2 *n1_Val);
	RooRealVar n2("n2", "n2", 0.5 *n2_Val, 2 *n2_Val);

	RooRealVar norm_dcb("norm_dcb", "norm_dcb", 0.5 *norm_Val, 2 *norm_Val);

	mean.setVal(mean_Val);
	sigma.setVal(sigma_Val);
	alpha1.setVal(alpha1_Val);
	alpha2.setVal(alpha2_Val);
	n1.setVal(n1_Val);
	n2.setVal(n2_Val);

	norm_dcb.setVal(norm_Val);

	mean.setConstant(true);
	sigma.setConstant(true);
	alpha1.setConstant(true);
	alpha2.setConstant(true);
	n1.setConstant(true);
	n2.setConstant(true);

	norm_dcb.setConstant(true);

	RooDoubleCB signal_dcb("signal_dcb", "DoubleCB", mbb, mean, sigma, alpha1, n1, alpha2, n2);

	outtext << endl;
	outtext << endl;
	outtext << "norm   = " << norm_Val << endl;
	outtext << "mean   = " << mean_Val << endl;
	outtext << "sigma  = " << sigma_Val << endl;
	outtext << endl;
	outtext << "Evaluation function within fitted range" << endl;
	int nbins = Hists["nominal"]->GetNbinsX();
	outtext << "number of bins = " << nbins << endl;

	for (int ib = 1; ib <= nbins; ++ib)
	{
		double xbin = Hists["nominal"]->GetXaxis()->GetBinCenter(ib);
		char xbinCh[10];
		sprintf(xbinCh, "%6.1f", xbin);
		TString xbinStr(xbinCh);
		mbb.setVal(xbin);
		double pdf = signal_dcb.getVal();
		//outtext << xbinStr << " : " << pdf << endl;
	}

	outtext << endl;
	outtext << endl;

	w->import(signal_dcb);
	w->import(norm_dcb);

}

void CreateSignalTemplates_FH()
{
	gROOT->SetBatch();
	// Silence INFO messages
	RooMsgService::instance().setGlobalKillBelow(RooFit::FATAL);
	// Silence additional MINUIT output 
	RooMsgService::instance().setSilentMode(true);
	// Silence Info in TCanvas::Print: messages
	gErrorIgnoreLevel = kWarning;

	vector<TString> histNames = { "nominal" };

	for (unsigned int i = 0; i < masses.size(); ++i)
	{

		RooWorkspace *w = new RooWorkspace("wsignal", "signal");

		int mass = masses[i];
		int region = mass_region[mass];
		int bin = mass_binning[mass];
		TString Mass = IntToStr(mass);
		TString Region = IntToStr(region);
		TString Bin = IntToStr(bin);
		map<TString, TH1D*> Hists;
		double mbb_min = mbb_low[region];
		double mbb_max = mbb_high[region];
		for (auto histName: histNames)
		{
			TFile *file = new TFile(dir + "/mssmHbb_2018_FH_" + Mass + "_sr.root", "READ");
			TH1D *hist = (TH1D*) file->Get("mbb");
			Hists[histName] = hist;				
		}
		CreateSignalPDF(mass, region, histNames, Hists, w, mbb_min, mbb_max);
		TFile *fileOutput = new TFile("input_doubleCB/signal_m" + Mass + "_SR" + Region + ".root", "recreate");
		fileOutput->cd("");
		w->Write("w");
		fileOutput->Write();
		fileOutput->Close();
	}
}
