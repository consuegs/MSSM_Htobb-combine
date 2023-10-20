
#include <iostream>
#include <map>

void compare_mbb_syst_jkte( int year = 2018, string channel = "FH")
{

   TH1::SetDefaultSumw2();  // proper treatment of errors when scaling histograms
   
   gStyle->SetOptStat(0);
         
   map<int, string> type = 
   {
    {2017,Form("POWHEG NLO signal %d%s", year, channel.c_str())},
    {2018,Form("POWHEG NLO signal %d%s", year, channel.c_str())},
   };
   
   map< int,  map<string, double>> physics_lumi;
   physics_lumi[2017]["SL"] = 36.67;
   physics_lumi[2017]["FH"] = 36.26;
   physics_lumi[2018]["FH"] = 54.54;
   
   
   map<string, vector <int>> samples_signal_mc;
   samples_signal_mc["SL"]={300,350,400,450,500,600,700,800,900,1000,1200,1400,1600,1800};
   samples_signal_mc["FH"]={300,350,400,450,500,600,700,800,900,1000,1200,1400,1600,1800};
   //samples_signal_mc["FH"]={600};
   
   //File and histograms
   TFile * file_nominal;
   TFile * file_up;
   TFile * file_down;
   vector <string> regions  = {"SR", "CR"};
   vector <string> Regions  = {"bbb", "bbnb"};
   double ratiomin = 0.7, ratiomax = 1.3;

   string outdir = Form("plots_mbb_%d%s_syst_jkte", year, channel.c_str());
   gSystem->Exec(Form("mkdir %s", outdir.c_str())); 
   
   
   //Binning section
   int bineo_option = 0; //for regular rebin use any number but 0 or 1. To use variable binning set to -1. To use bin map set to 0
   int bineo_mbb = 10; //for regular rebin. To use variable binning set to -1. To use bin map set to 0
   vector <double> xbins;
   double xbins_array[64];
   xbins.push_back(0);
   double bin;
   for (int i = 1; i <= 12; i++) {bin = 10*i; xbins.push_back(bin);}
   for (int i = 13; i <= 31; i++) {bin = 120+20*(i-12); xbins.push_back(bin);}
   for (int i = 32; i <= 42; i++) {bin = 500+30*(i-31); xbins.push_back(bin);}
   for (int i = 43; i <= 50; i++) {bin = 800+50*(i-42); xbins.push_back(bin);}
   for (int i = 51; i <= 58; i++) {bin = 1200+100*(i-50); xbins.push_back(bin);}
   for (int i = 59; i <= 62; i++) {bin = 2000+150*(i-58); xbins.push_back(bin);}
   for (int i = 63; i <= 63; i++) {bin = 2600+400*(i-62); xbins.push_back(bin);}
   for (int bin = 0; bin<= 63; bin++) {xbins_array[bin] = xbins[bin];}
   //for (int i = 0; i<=63; i++) {cout<<endl<<"bin "<<i<<" xbin "<<xbins_array[i];}
   //cout<<endl;
   
   map <int, int> bineo_map;
   bineo_map[300] = 10;
   bineo_map[350] = 10;
   bineo_map[400] = 10;
   bineo_map[450] = 10;
   bineo_map[500] = 10;
   bineo_map[600] = 10;
   bineo_map[700] = 10;
   bineo_map[800] = 10;
   bineo_map[900] = 20;
   bineo_map[1000] = 20;
   bineo_map[1200] = 20;
   bineo_map[1400] = 20;
   bineo_map[1600] = 20;
   bineo_map[1800] = 20;
   
   //end of Binning section
   
      
   //Histogram ranges section
   map <int, int> xmin_mbb, xmax_mbb;
   
   xmin_mbb[300] = 260;
   xmin_mbb[350] = 260;
   xmin_mbb[400] = 260;
   xmin_mbb[450] = 260;
   xmin_mbb[500] = 320;
   xmin_mbb[600] = 320;
   xmin_mbb[700] = 380;
   xmin_mbb[800] = 400;
   xmin_mbb[900] = 400;
   xmin_mbb[1000] = 450;
   xmin_mbb[1200] = 550;
   xmin_mbb[1400] = 550;
   xmin_mbb[1600] = 550;
   xmin_mbb[1800] = 550;
   
   xmax_mbb[300] = 500;
   xmax_mbb[350] = 500;
   xmax_mbb[400] = 500;
   xmax_mbb[450] = 600;
   xmax_mbb[500] = 650;
   xmax_mbb[600] = 800;
   xmax_mbb[700] = 900;
   xmax_mbb[800] = 1000;
   xmax_mbb[900] = 1200;
   xmax_mbb[1000] = 1300;
   xmax_mbb[1200] = 1500;
   xmax_mbb[1400] = 1800;
   xmax_mbb[1600] = 2000;
   xmax_mbb[1800] = 2000;
   
   int xmin_pT = 0, xmax_pT;
   
     
   // end of Histogram ranges section
   
   //the mbb histograms
   map <string, TH1F *> mbb; //mbb for each strategy
   map <string, TH1F *> mbb_rebinned; //mbb for each strategy
   map <string, TH1F *> mbb_normalized; //mbb for each strategy
   
   map <string, EColor> hist_color  = 
   {
    {"nominal", kBlack},
    {"up", kBlue},
    {"down", kRed}
   };
   
   //map <string, TLegend * > mbb_legend;
   TLegend *  mbb_legend;
   
   auto C = new TCanvas("C","Canvas",800,800);
    
   //CMS Labeling
   TLatex latex;
   latex.SetTextFont(43);
   latex.SetTextSize(30);
   latex.SetTextAlign(11);
   latex.DrawLatexNDC(gPad->GetLeftMargin()+0.3, 1.02-C->GetTopMargin(),
                Form("CMS Work in progress #sqrt{s} = 13 TeV, L = %.2f fb^{-1}",physics_lumi[year][channel]));
                
                
   std::ofstream outFile("2018FH_systematics_jkte.csv");
     if (!outFile.is_open())
     {
        std::cout << "Error opening the output file." << std::endl;
        return 1; // Exit with an error code
     }
     
   outFile << "Higgs_mass_GeV,systematics_down,systematics_up" << std::endl;

   
   //loop over mass-points
   for (int mp = 0; mp < samples_signal_mc["FH"].size() ; mp++)
   {
      TPad *pad1 = new TPad("pad1", "pad1", 0, 0.3, 1, 1.0);
      pad1->SetBottomMargin(0.02); // Upper and lower plot are joined
      pad1->Draw();             // Draw the upper pad: pad1
      pad1->cd();               // pad1 becomes the current pad
      
      //TFiles 
      file_nominal = new TFile(Form("/home/leyvaped/Work2023/20230725_onlinebtagSF/systematics/rootfiles_signal_p1/mssmHbb_%d_%s_%d_sr.root", year,channel.c_str(), samples_signal_mc["FH"][mp] ));
      file_up      = new TFile(Form("rootfiles_signal/mssmHbb_%d_%s_%d_JKTE_1sup.root", year,channel.c_str(), samples_signal_mc["FH"][mp] ));
      file_down    = new TFile(Form("rootfiles_signal/mssmHbb_%d_%s_%d_JKTE_1sdown.root", year,channel.c_str(), samples_signal_mc["FH"][mp] ));
      
      //get mbb histograms
      mbb["nominal"]  = (TH1F*) file_nominal -> Get("mbb"); 
      mbb["up"]       = (TH1F*) file_up -> Get("mbb"); 
      mbb["down"]     = (TH1F*) file_down -> Get("mbb");      
      
      //Cosmetics section
      mbb["nominal"] -> SetMarkerStyle(20);
      mbb["nominal"] -> SetMarkerColor(hist_color["nominal"]);
      mbb["nominal"] -> SetLineColor(hist_color["nominal"]);
      mbb["up"] -> SetMarkerStyle(20);
      mbb["up"] -> SetMarkerColor(hist_color["up"]);
      mbb["up"] -> SetLineColor(hist_color["up"]);
      mbb["down"] -> SetMarkerStyle(20);
      mbb["down"] -> SetMarkerColor(hist_color["down"]);
      mbb["down"] -> SetLineColor(hist_color["down"]);
      
      //mbb[regions[sr]] -> SetTitle(Form("mbb %s",regions[sr].c_str()));
      mbb["nominal"] -> SetTitle("");
      mbb["up"] -> SetTitle("");
      mbb["down"] -> SetTitle("");
      
      if(bineo_option != 0)
      {
         if (fabs(bineo_option) != 1)
         {
            mbb["nominal"] -> GetYaxis() -> SetTitle(Form("Entries/%d [GeV]", bineo_mbb));
            mbb["up"] -> GetYaxis() -> SetTitle(Form("Entries/%d [GeV]", bineo_mbb));
            mbb["down"] -> GetYaxis() -> SetTitle(Form("Entries/%d [GeV]", bineo_mbb));
         }
         else
         {
            mbb["nominal"] -> GetYaxis() -> SetTitle("Entries" );
            mbb["up"] -> GetYaxis() -> SetTitle("Entries" );
            mbb["down"] -> GetYaxis() -> SetTitle("Entries" );
         }
      }
      else
      {
         mbb["nominal"] -> GetYaxis() -> SetTitle(Form("Entries/%d [GeV]", bineo_map[samples_signal_mc["FH"][mp]]));
         mbb["up"] -> GetYaxis() -> SetTitle(Form("Entries/%d [GeV]", bineo_map[samples_signal_mc["FH"][mp]]));
         mbb["down"] -> GetYaxis() -> SetTitle(Form("Entries/%d [GeV]", bineo_map[samples_signal_mc["FH"][mp]]));
      }
      
      
      mbb["nominal"] -> GetXaxis() -> SetTitle("M_{12}");
      mbb["up"]      -> GetXaxis() -> SetTitle("M_{12}");
      mbb["down"]    -> GetXaxis() -> SetTitle("M_{12}");
      
      
      
      //end of Cosmetics section
      
      
      if(bineo_option != 0)
      {
         if(bineo_option == -1)
         {
            mbb_rebinned["nominal"] = (TH1F*)mbb["nominal"]->Rebin(63,"mbb_rebin_nominal");	
            mbb_rebinned["up"] = (TH1F*)mbb["up"]->Rebin(63,"mbb_rebin_up");	
            mbb_rebinned["down"] = (TH1F*)mbb["down"]->Rebin(63,"mbb_rebin_down");	
         }
         else
         {
            mbb_rebinned["nominal"] = (TH1F*)mbb["nominal"]->Clone("mbb_rebinned_nominal");
            mbb_rebinned["up"] = (TH1F*)mbb["up"]->Clone("mbb_rebinned_up");
            mbb_rebinned["down"] = (TH1F*)mbb["down"]->Clone("mbb_rebinned_down");
            mbb_rebinned["nominal"] -> Rebin(bineo_mbb);
            mbb_rebinned["up"] -> Rebin(bineo_mbb);
            mbb_rebinned["down"] -> Rebin(bineo_mbb);
         } 
      }
      else
      {
         mbb_rebinned["nominal"] = (TH1F*)mbb["nominal"]->Clone("mbb_rebinned_nominal");
         mbb_rebinned["up"] = (TH1F*)mbb["up"]->Clone("mbb_rebinned_up");
         mbb_rebinned["down"] = (TH1F*)mbb["down"]->Clone("mbb_rebinned_down");

         mbb_rebinned["nominal"] -> Rebin(bineo_map[samples_signal_mc["FH"][mp]]);
         mbb_rebinned["up"] -> Rebin(bineo_map[samples_signal_mc["FH"][mp]]);
         mbb_rebinned["down"] -> Rebin(bineo_map[samples_signal_mc["FH"][mp]]);
      }
      
      mbb_rebinned["nominal"] -> GetXaxis() -> SetRangeUser(xmin_mbb[samples_signal_mc["FH"][mp]],xmax_mbb[samples_signal_mc["FH"][mp]]);
      mbb_rebinned["up"] -> GetXaxis() -> SetRangeUser(xmin_mbb[samples_signal_mc["FH"][mp]],xmax_mbb[samples_signal_mc["FH"][mp]]);
      mbb_rebinned["down"] -> GetXaxis() -> SetRangeUser(xmin_mbb[samples_signal_mc["FH"][mp]],xmax_mbb[samples_signal_mc["FH"][mp]]);
      //gPad->SetLogy();
          
      mbb_rebinned["nominal"]  -> SetFillColor(kRed-10);
      mbb_rebinned["nominal"]->GetXaxis()->SetLabelOffset(999);
      mbb_rebinned["up"]->GetXaxis()->SetLabelOffset(999);
      mbb_rebinned["down"]->GetXaxis()->SetLabelOffset(999);
      
      mbb_rebinned["nominal"] -> Draw("hist");
      mbb_rebinned["up"] -> Draw("same");
      mbb_rebinned["down"] -> Draw("same");
      
      //mbb legend
      if (samples_signal_mc["FH"][mp] == 300 || samples_signal_mc["FH"][mp] == 350 ) //ok
      mbb_legend = new TLegend(0.65,0.7,0.85,0.85);
      else
      mbb_legend = new TLegend(0.15,0.7,0.35,0.85);
      mbb_legend -> SetBorderSize(0);
      mbb_legend -> AddEntry((TObject*)0, Form("m_{A/H} = %d GeV,", samples_signal_mc["FH"][mp]), "");
      //mbb_legend -> AddEntry((TObject*)0, "HEM correction", "");
      mbb_legend -> AddEntry(mbb_rebinned["nominal"],"central","f");      
      mbb_legend -> AddEntry(mbb_rebinned["up"],"jet kin +1#sigma variation","lep");
      mbb_legend -> AddEntry(mbb_rebinned["down"],"jet kin -1#sigma variation","lep");
      mbb_legend -> Draw();
      
      
      //plot mbb and legend  
      latex.SetTextFont(43);
      latex.SetTextSize(20);
      latex.SetTextAlign(11);
      latex.DrawLatexNDC(gPad->GetLeftMargin(), 1.02-C->GetTopMargin(),
                   (std::string("CMS Simulation      Work in progress")).c_str());
      latex.DrawLatexNDC(gPad->GetLeftMargin()+0.7, 1.02-C->GetTopMargin(),
                   (std::string("(13 TeV)")).c_str());
  
  
      mbb["nominal"] -> GetXaxis() -> SetRangeUser(xmin_mbb[samples_signal_mc["FH"][mp]],xmax_mbb[samples_signal_mc["FH"][mp]]);
      mbb["up"] -> GetXaxis() -> SetRangeUser(xmin_mbb[samples_signal_mc["FH"][mp]],xmax_mbb[samples_signal_mc["FH"][mp]]); 
      mbb["down"] -> GetXaxis() -> SetRangeUser(xmin_mbb[samples_signal_mc["FH"][mp]],xmax_mbb[samples_signal_mc["FH"][mp]]); 
      mbb["nominal"] -> SetMinimum(0);
      mbb["up"] -> SetMinimum(0); 
      mbb["down"] -> SetMinimum(0); 
      
      //plot ratio
      C->cd();          // Go back to the main canvas before defining pad2
      TPad *pad2 = new TPad("pad2", "pad2", 0, 0.05, 1, 0.3);
      pad2->SetTopMargin(0.05);
      pad2->SetBottomMargin(0.3);
      pad2->SetGridy();   // horizontal grid
      pad2->Draw();
      pad2->cd();       // pad2 becomes the current pad
      
      TH1F *ratio_up  = (TH1F*)mbb_rebinned["up"]->Clone("ratio");
      ratio_up->Divide(mbb_rebinned["nominal"]);
      ratio_up -> SetMarkerColor(kBlue);
      ratio_up -> SetLineColor(kBlack);
      ratio_up -> GetYaxis() -> SetRangeUser(ratiomin,ratiomax);
      ratio_up -> GetXaxis() -> SetRangeUser(xmin_mbb[samples_signal_mc["FH"][mp]],xmax_mbb[samples_signal_mc["FH"][mp]]);
   
      ratio_up -> GetYaxis() -> SetTitle("Down(Up)/Central");
      //ratio_up -> GetYaxis()->SetNdivisions(505);
      ratio_up -> SetStats(0);      // No statistics on lower plot
      ratio_up -> SetMarkerStyle(20);
      ratio_up -> SetMarkerSize(0.8);
      
      ratio_up ->GetYaxis()->SetTitleSize(20);
      ratio_up ->GetYaxis()->SetTitleFont(43);
      ratio_up ->GetYaxis()->SetTitleOffset(1.55);
      ratio_up ->GetYaxis()->SetLabelFont(43); // Absolute font size in pixel (precision 3)
      ratio_up ->GetYaxis()->SetLabelSize(15);
 
   // X axis ratio plot settings
      ratio_up ->GetXaxis()->SetLabelOffset(0);
      ratio_up ->GetXaxis()->SetTitleSize(20);
      ratio_up ->GetXaxis()->SetTitleFont(43);
      ratio_up ->GetXaxis()->SetTitleOffset(1);
      ratio_up ->GetXaxis()->SetLabelFont(43); // Absolute font size in pixel (precision 3)
      ratio_up ->GetXaxis()->SetLabelSize(15);
      
      ratio_up -> Draw();
      
      TH1F *ratio_down  = (TH1F*)mbb_rebinned["down"]->Clone("ratio");
      ratio_down->Divide(mbb_rebinned["nominal"]);
      ratio_down -> SetMarkerColor(kRed);
      ratio_down -> SetLineColor(kBlack);
      ratio_down -> GetYaxis() -> SetRangeUser(ratiomin,ratiomax);
      ratio_down -> GetXaxis() -> SetRangeUser(xmin_mbb[samples_signal_mc["FH"][mp]],xmax_mbb[samples_signal_mc["FH"][mp]]);
   
      ratio_down -> GetYaxis() -> SetTitle("Down(Up)/Central");
      //ratio_down -> GetYaxis()->SetNdivisions(505);
      ratio_down -> SetStats(0);      // No statistics on lower plot
      ratio_down -> SetMarkerStyle(20);
      ratio_down -> SetMarkerSize(0.8);
      
      ratio_down ->GetYaxis()->SetTitleSize(20);
      ratio_down ->GetYaxis()->SetTitleFont(43);
      ratio_down ->GetYaxis()->SetTitleOffset(1.55);
      ratio_down ->GetYaxis()->SetLabelFont(43); // Absolute font size in pixel (precision 3)
      ratio_down ->GetYaxis()->SetLabelSize(15);
 
   // X axis ratio plot settings
      ratio_down ->GetXaxis()->SetLabelOffset(0);
      ratio_down ->GetXaxis()->SetTitleSize(20);
      ratio_down ->GetXaxis()->SetTitleFont(43);
      ratio_down ->GetXaxis()->SetTitleOffset(1);
      ratio_down ->GetXaxis()->SetLabelFont(43); // Absolute font size in pixel (precision 3)
      ratio_down ->GetXaxis()->SetLabelSize(15);
      
      ratio_down -> Draw("same");
      
      C -> SaveAs(Form("%s/jkte_%d_GeV.png", outdir.c_str(), samples_signal_mc["FH"][mp])); 
      C -> SaveAs(Form("%s/jkte_%d_GeV.pdf", outdir.c_str(), samples_signal_mc["FH"][mp]));
      C -> SaveAs(Form("%s/jkte_%d_GeV.root", outdir.c_str(), samples_signal_mc["FH"][mp]));  
      
      //Implementing the fit function
      TF1 *f1 ;
      
      if (samples_signal_mc["FH"][mp] == 300 || samples_signal_mc["FH"][mp] == 350 ) 
      f1 = new TF1("constant","[0]",mbb["nominal"] -> GetMean() - 75 ,mbb["nominal"] -> GetMean() + 75);
      
      else if (samples_signal_mc["FH"][mp] >= 400 && samples_signal_mc["FH"][mp] < 600 ) 
      f1 = new TF1("constant","[0]",mbb["nominal"] -> GetMean() - 100 ,mbb["nominal"] -> GetMean() + 100);
      
      else if (samples_signal_mc["FH"][mp] >= 600 && samples_signal_mc["FH"][mp] < 800 )
      f1 = new TF1("constant","[0]",mbb["nominal"] -> GetMean() - 150 ,mbb["nominal"] -> GetMean() + 150);
    
      else if (samples_signal_mc["FH"][mp] >= 800 && samples_signal_mc["FH"][mp] < 1000 ) 
      f1 = new TF1("constant","[0]",mbb["nominal"] -> GetMean() - 200 ,mbb["nominal"] -> GetMean() + 200);
      
      else if (samples_signal_mc["FH"][mp] >= 1000 && samples_signal_mc["FH"][mp] < 1400 ) 
      f1 = new TF1("constant","[0]",mbb["nominal"] -> GetMean() - 300 ,mbb["nominal"] -> GetMean() + 300);
      
      else if (samples_signal_mc["FH"][mp] >= 1400 && samples_signal_mc["FH"][mp] < 1800 ) 
      f1 = new TF1("constant","[0]",mbb["nominal"] -> GetMean() - 350 ,mbb["nominal"] -> GetMean() + 350);
      
      else if (samples_signal_mc["FH"][mp] >= 1800 ) 
      f1 = new TF1("constant","[0]",mbb["nominal"] -> GetMean() - 350 ,mbb["nominal"] -> GetMean() + 500);
      
      f1->SetParName(0, "C");
      f1->SetParameters(0);
      f1->SetLineColor(kBlue);
      
      
      cout << "--------------------------------------------------------------------------------------"<< endl;
      cout << "masspoint: " << samples_signal_mc["FH"][mp] << " GeV"<< endl;
      TFitResultPtr rup = ratio_up->Fit(f1,"RS");
      cout << "Fitting Up/Central" << endl;
      cout << "chi2/ndf = " << rup->Chi2() << "/" << rup->Ndf() << " = " << rup->Chi2()/rup->Ndf() << ", prob = " << rup->Prob() << endl;
      //double c_up = rup->Parameter(0);
      double c_up = mbb_rebinned["up"]->GetSumOfWeights() / mbb_rebinned["nominal"]->GetSumOfWeights();
      cout << "--------------------------------------------------------------------------------------"<< endl;
      
      f1->SetLineColor(kRed);
      cout << "--------------------------------------------------------------------------------------"<< endl;
      cout << "masspoint: " << samples_signal_mc["FH"][mp] << " GeV"<< endl;
      TFitResultPtr rdown = ratio_down->Fit(f1,"RS");
      cout << "Fitting Down/Central" << endl;
      cout << "chi2/ndf = " << rdown->Chi2() << "/" << rdown->Ndf() << " = " << rdown->Chi2()/rdown->Ndf() << ", prob = " << rdown->Prob() << endl;
      //double c_down = rdown->Parameter(0);      
      double c_down = mbb_rebinned["down"]->GetSumOfWeights() / mbb_rebinned["nominal"]->GetSumOfWeights();
      cout << "--------------------------------------------------------------------------------------"<< endl;

      C -> SaveAs(Form("%s/jkte_%d_GeV_param.png", outdir.c_str(), samples_signal_mc["FH"][mp])); 
      C -> SaveAs(Form("%s/jkte_%d_GeV_param.pdf", outdir.c_str(), samples_signal_mc["FH"][mp]));
      C -> SaveAs(Form("%s/jkte_%d_GeV_param.root", outdir.c_str(), samples_signal_mc["FH"][mp]));  
      
      C->Clear();
      outFile << samples_signal_mc["FH"][mp] << "," << c_down << "," << c_up << std::endl;
   }
   outFile.close();
}
