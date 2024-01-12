void plotMassVariations() {

  TH1::SetDefaultSumw2();
  // List of mass points
  const int numMassPoints = 14;
  int massPoints[numMassPoints] = {300,350,400,450,500,600,700,800,900,1000,1200,1400,1600,1800};
  std::string syst = "JES";
  //Signal region map
  std::map<int, int> signalregion = {
        {300, 1},
        {350, 1},
        {400, 2},
        {450, 2},
        {500, 2},
        {600, 3},
        {700, 3},
        {800, 3},
        {900, 3},
        {1000, 4},
        {1200, 4},
        {1400, 4},
        {1600, 4},
        {1800, 4}
    };
    
  std::map<int, float> xmin = {
        {1, 270},
        {2, 320},
        {3, 390},
        {4, 500}
    };
    
  std::map<int, float> xmax = {
        {1, 560},
        {2, 800},
        {3, 1270},
        {4, 2000}
    };

  map<int, int> samples_signal_mc;
  samples_signal_mc[300]={1};
  samples_signal_mc[350]={1};
  samples_signal_mc[400]={2};

  // Loop over mass points
  for (int i = 0; i < numMassPoints; ++i) {
    int mass = massPoints[i];

    // Create canvas for plotting
    TCanvas* canvas = new TCanvas("canvas", "Signal Mass Variations", 600, 600);
    canvas->SetLeftMargin(0.15);
    canvas->SetRightMargin(0.05);
    // Open the root files for the three variations
    TString rootFileNominal = TString::Format("signal_m%d_SR%d.root", mass, signalregion[mass]);
    TString rootFileUp = TString::Format("signal_m%d_SR%d_%s_1sup.root", mass, signalregion[mass],  syst.c_str());
    TString rootFileDown = TString::Format("signal_m%d_SR%d_%s_1sdown.root", mass, signalregion[mass],  syst.c_str());

    TFile* fileNominal = TFile::Open(rootFileNominal);
    TFile* fileUp = TFile::Open(rootFileUp);
    TFile* fileDown = TFile::Open(rootFileDown);

    if (!fileNominal || !fileUp || !fileDown) {
      std::cerr << "Error: Could not open root files for mass " << mass << std::endl;
      continue;
    }

    // Get the RooWorkspace objects from the root files
    RooWorkspace* wsNominal = (RooWorkspace*)fileNominal->Get("w");
    RooWorkspace* wsUp = (RooWorkspace*)fileUp->Get("w");
    RooWorkspace* wsDown = (RooWorkspace*)fileDown->Get("w");

    if (!wsNominal || !wsUp || !wsDown) {
      std::cerr << "Error: Could not retrieve RooWorkspace for mass " << mass << std::endl;
      continue;
    }

    // Get the RooDataHist objects
    RooDataHist* dataNominal = (RooDataHist*)wsNominal->data("sig");
    RooDataHist* dataUp = (RooDataHist*)wsUp->data("sig");
    RooDataHist* dataDown = (RooDataHist*)wsDown->data("sig");

    // Get the RooAbsPdf objects
    RooAbsPdf* pdfNominal = wsNominal->pdf("signal_dcb");
    RooAbsPdf* pdfUp = wsUp->pdf("signal_dcb");
    RooAbsPdf* pdfDown = wsDown->pdf("signal_dcb");

    // Plot the nominal data and pdf in black
    RooRealVar mbb("mbb", "mbb", xmin[signalregion[mass]], xmax[signalregion[mass]]); 
    RooPlot * frame1 = mbb.frame(RooFit::Title(" "));
    

    dataNominal->plotOn(frame1, RooFit::Name("dataNominal_plot"), RooFit::MarkerColor(kBlack), RooFit::LineColor(kBlack), RooFit::MarkerStyle(20), RooFit::LineWidth(2));
    pdfNominal->plotOn(frame1, RooFit::LineColor(kBlack), RooFit::LineWidth(2));

    // Plot the up variation data and pdf in blue
    dataUp->plotOn(frame1, RooFit::Name("dataUp_plot"),RooFit::MarkerColor(kBlue), RooFit::LineColor(kBlue), RooFit::MarkerStyle(21), RooFit::LineWidth(2));
    pdfUp->plotOn(frame1, RooFit::LineColor(kBlue), RooFit::LineWidth(2));

    // Plot the down variation data and pdf in red
    dataDown->plotOn(frame1, RooFit::Name("dataDown_plot"), RooFit::MarkerColor(kRed), RooFit::LineColor(kRed), RooFit::MarkerStyle(22), RooFit::LineWidth(2));
    pdfDown->plotOn(frame1, RooFit::LineColor(kRed), RooFit::LineWidth(2));

    frame1->Draw();

    //some cosmetics
    TLatex latex;
    latex.SetTextFont(43);
    latex.SetTextSize(20);
    latex.SetTextAlign(11);
    latex.DrawLatexNDC(gPad->GetLeftMargin(), 1.02 - canvas->GetTopMargin(), "CMS Simulation      Work in progress");
    latex.DrawLatexNDC(gPad->GetLeftMargin() + 0.7, 1.02 - canvas->GetTopMargin(), "(13 TeV)");


    TLegend *  mbb_legend;
    //if (massPoints[i] == 300 || massPoints[i] == 350 ) 
    mbb_legend = new TLegend(0.6,0.7,0.8,.85);
    //else
    //mbb_legend = new TLegend(0.15,0.7,0.35,0.85);
    mbb_legend -> SetBorderSize(0);
    mbb_legend -> AddEntry((TObject*)0, Form("m_{A/H} = %d GeV,", massPoints[i]), "");
    mbb_legend -> AddEntry(frame1->findObject("dataNominal_plot"),"central","lep");
    mbb_legend -> AddEntry(frame1->findObject("dataUp_plot"),Form("%s +1#sigma variation", syst.c_str()),"lep");
    mbb_legend -> AddEntry(frame1->findObject("dataDown_plot"),Form("%s -1#sigma variation", syst.c_str()),"lep");
    mbb_legend -> Draw();

    // Save the plot in pdf and png formats
    TString outputFileName = TString::Format("plots/%s_m%d.pdf", syst.c_str(), mass);
    canvas->SaveAs(outputFileName);
    outputFileName.ReplaceAll(".pdf", ".png");
    canvas->SaveAs(outputFileName);

    // Close the root files
    fileNominal->Close();
    fileUp->Close();
    fileDown->Close();
    delete fileNominal;
    delete fileUp;
    delete fileDown;

    // Delete the canvas
    delete canvas;
  }
}
