void GetLumiSF()
{

	float lumi_2018 = 54540; // pb^-1

	vector<int> masses = { 300, 350, 400, 450, 500, 600, 700, 800, 900, 1000, 1200, 1400, 1600, 1800 };
	vector<string> srmasses = { "300", "350", "400", "450", "500", "600", "700", "800", "900", "1000", "1200", "1400", "1600", "1800" };

	std::ofstream outtext("LumiSF_2018.txt");
	outtext << "Lumi Scale Factors:" << endl;
	outtext << endl;

	TString dir("/afs/desy.de/user/l/leyvaped/public/for_sandra/rootfiles_2018FH_Feb2023/signal_2018FH_nominal/");

	for (unsigned int i = 0; i < masses.size(); ++i)
	{

		int mass = masses[i];

		TString srmass = srmasses[i];

		TFile *file = new TFile(dir + "/mssmHbb_2018_FH_" + srmass + "_sr.root", "READ");

		TH1D *hist = (TH1D*) file->Get("nentries");

		//float lumiSF = (hist->GetBinContent(2) - hist->GetBinContent(1)) / lumi_2018;
		//cout<<endl<<"Nentries from  (hist->GetBinContent(2) - hist->GetBinContent(1)) "<< int(hist->GetBinContent(2) - hist->GetBinContent(1));

		TH1F* workflow = (TH1F*) file -> Get("workflow");
		double genweightEvents = workflow -> GetBinContent(workflow->GetXaxis()->FindBin("Generated weighted events (sign of weights)"));
		float lumiSF = genweightEvents / lumi_2018;
		//cout<<endl<<"Nentries workflow "<< int(genweightEvents);

		outtext << endl;
		outtext << mass << " GeV: " << setprecision(4) << lumiSF << endl;
		outtext << endl;
	}
}
