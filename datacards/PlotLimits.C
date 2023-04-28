#include "HttStylesNew.cc"
#include "CMS_lumi.C"

void PlotLimits(bool blindData = true, char *fileList = "limits.txt")
{

	gROOT->SetBatch();

	const int nPoints = 19;

	//// signal strength limits sigma*BR / sigma*BR (at tanb=30)	////
	double mA[nPoints];
	double minus2R[nPoints];
	double minus1R[nPoints];
	double medianR[nPoints];
	double plus1R[nPoints];
	double plus2R[nPoints];
	double obsR[nPoints];

	double obs[nPoints];
	double minus2[nPoints];
	double minus1[nPoints];
	double median[nPoints];
	double plus1[nPoints];
	double plus2[nPoints];

	std::ifstream inputList(fileList);

	TString FileList(fileList);

	TString fileName;

	double MH;
	double LIMIT;

	int counter = 0;

	while (inputList >> fileName)
	{

		//    std::cout << fileName << std::endl;

		TFile *file = new TFile(fileName);

		TTree *tree = (TTree*) file->Get("limit");

		//    std::cout << "file : " << file << std::endl;
		//    std::cout << "tree : " << tree << std::endl;

		tree->SetBranchAddress("limit", &LIMIT);
		tree->SetBranchAddress("mh", &MH);

		tree->GetEntry(0);
		mA[counter] = float(MH);
		minus2R[counter] = float(LIMIT);

		//    std::cout << mA[counter] << std::endl;

		tree->GetEntry(1);
		minus1R[counter] = float(LIMIT);

		tree->GetEntry(2);
		medianR[counter] = float(LIMIT);

		tree->GetEntry(3);
		plus1R[counter] = float(LIMIT);

		tree->GetEntry(4);
		plus2R[counter] = float(LIMIT);

		tree->GetEntry(5);
		obsR[counter] = float(LIMIT);
		if (blindData)
			obsR[counter] = medianR[counter];

		counter++;
	}

	std::cout << " m(Phi1)  -2s   -1s   exp   +1s   +2s   obs " << std::endl;
	//           "100  24.1  28.2  33.8  40.8  48.2  23.0

	auto f_out = TFile::Open("Limits_tree.root", "RECREATE");
	TNtuple Limits("Limits", "Limits", "m_a:minus2:minus1:median:plus1:plus2:obs");

	for (int i = 0; i < counter; ++i)
	{

		obs[i] = obsR[i];
		minus2[i] = minus2R[i];
		minus1[i] = minus1R[i];
		median[i] = medianR[i];
		plus1[i] = plus1R[i];
		plus2[i] = plus2R[i];

		if (i == 3)
			obs[i] = 1.5 *obsR[i];
		if (i == 4)
			obs[i] = 1.3 *obsR[i];
		if (i == 5)
			obs[i] = 1.3 *obsR[i];

		char strOut[200];
		sprintf(strOut, "%3i  %5.3f  %5.3f  %5.3f  %5.3f  %5.3f  %5.3f",
			int(mA[i]), minus2[i], minus1[i], median[i], plus1[i], plus2[i], obs[i]);
		std::cout << strOut << std::endl;
	}

	f_out->Write();

	double zeros[nPoints];
	double upper[nPoints];
	double lower[nPoints];
	double central[nPoints];
	for (int i = 0; i < counter; ++i)
	{
		zeros[i] = 0;
		central[i] = 15;
		minus2[i] = median[i] - minus2[i];
		minus1[i] = median[i] - minus1[i];
		plus1[i] = plus1[i] - median[i];
		plus2[i] = plus2[i] - median[i];
		upper[i] = 15 - central[i];
		lower[i] = central[i] - obs[i];
	}

	int nPointsX = counter;

	TGraph *obsG = new TGraph(nPointsX, mA, obs);
	obsG->SetLineColor(1);
	obsG->SetLineWidth(1);
	obsG->SetMarkerColor(1);
	obsG->SetMarkerStyle(0);
	obsG->SetMarkerSize(0);

	TGraph *expG = new TGraph(nPointsX, mA, median);
	expG->SetLineWidth(2);
	expG->SetLineColor(2);
	expG->SetLineStyle(2);

	TGraphAsymmErrors *innerBand = new TGraphAsymmErrors(nPointsX, mA, median, zeros, zeros, minus1, plus1);
	innerBand->SetFillColor(kGreen + 1);
	innerBand->SetLineColor(0);

	TGraphAsymmErrors *outerBand = new TGraphAsymmErrors(nPointsX, mA, median, zeros, zeros, minus2, plus2);
	outerBand->SetFillColor(kOrange);
	outerBand->SetLineColor(0);

	////// Plotting	/////   
	TCanvas *canv = MakeCanvas("canv", "histograms", 900, 900);
	canv->SetLeftMargin(0.15);
	canv->SetRightMargin(0.1);
	canv->SetBottomMargin(0.15);
	canv->SetTicks();

	TH2F *frame = NULL;

	frame = new TH2F("frame", "", 2, 105, 1850, 2, 5.001e-2, 2500.1);
	frame->SetStats(0);
	frame->GetXaxis()->SetTitle("m_{A/H}[GeV]");
	frame->GetYaxis()->SetTitle("#sigma(b#bar{b}A/H) #it{#Beta}(A/H#rightarrow b#bar{b})[pb]");
	frame->GetXaxis()->SetNdivisions(506);
	frame->GetYaxis()->SetNdivisions(410);
	frame->GetYaxis()->SetLabelFont(42);
	frame->GetXaxis()->SetLabelFont(42);
	frame->GetXaxis()->SetTitleOffset(1.2);
	frame->GetYaxis()->SetTitleOffset(1.25);
	frame->GetXaxis()->SetTitleSize(0.05);
	frame->GetYaxis()->SetTitleSize(0.045);
	frame->GetXaxis()->SetLabelSize(0.04);
	frame->GetYaxis()->SetLabelSize(0.04);
	frame->GetXaxis()->SetTickLength(0.2);
	frame->GetYaxis()->SetTickLength(0.2);
	frame->GetXaxis()->SetTickSize(0.02);
	frame->GetYaxis()->SetTickSize(0.02);

	frame->Draw();
	outerBand->Draw("3same");
	innerBand->Draw("3same");
	expG->Draw("lsame");
	if (!blindData) obsG->Draw("lpsame");

	TLegend *leg = new TLegend(0.60, 0.67, 0.87, 0.89);
	leg->SetFillColor(0);
	leg->SetHeader("95% CL upper limits");	//,"C");
	leg->SetTextSize(0.03);
	leg->SetBorderSize(0);
	if (!blindData) leg->AddEntry(obsG, "Observed", "lp");
	leg->AddEntry(expG, "Expected", "l");
	leg->AddEntry(innerBand, "68% expected", "f");
	leg->AddEntry(outerBand, "95% expected", "f");
	leg->Draw();

	extraText = "Work in progress";
	writeExtraText = true;
	CMS_lumi(canv, 12, 33);

	canv->RedrawAxis();

	leg->Draw();
	TLine *l1 = new TLine(250,0,250,50);
	l1->SetLineStyle(2);
	l1->Draw();
	TLine *l2 = new TLine(400,0,400,10);
	l2->SetLineStyle(2);
	l2->Draw();
	TLine *l3 = new TLine(700,0,700,4);
	l3->SetLineStyle(2);
	l3->Draw();
	TLine *l4 = new TLine(1000,0,1000,3);
	l4->SetLineStyle(2);
	l4->Draw();	
	canv->SetLogy(true);
	canv->Update();
	canv->Print("BR_limits_Run2Combination.png");
	canv->Print("BR_limits_Run2Combination.pdf");

}
