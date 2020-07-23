#include "wrapper.hpp"
#include "MeanRMSCalc.hpp"
#include "ErrorBarAnalysis.hpp"
#include "WaveformFitResult.hpp"
#include "ScanPoint.hpp"
#include "PTFAnalysis.hpp"
#include "PTFQEAnalysis.hpp"
#include "Utilities.hpp"
#include <string>
#include <iostream>
#include <ostream>
#include <fstream>
#include <sstream>
#include <exception>
#include <stdio.h>
#include <stdlib.h>
#include <vector>
#include <algorithm>
#include <cmath>
#include <unordered_set>

#include "TH1D.h"
#include "TH2D.h"
#include "TFile.h"
#include "TF1.h"
#include "TVirtualFitter.h"
#include "TFitter.h"
#include "TMath.h"
#include "TStyle.h"
#include "TDirectory.h"
#include "TFile.h"
#include "TH2D.h"
#include "TStyle.h"
#include "TCanvas.h"
#include "TGraph.h"//scatter plot
#include "TLegend.h"



int main( int argc, char* argv[] ) {

  if ( argc != 3 ){
    std::cerr<<"Usage: ptf_qe_analysis.app ptf_analysis.root run_number\n";
    exit(0);
  }

  // Get utilities
  Utilities utils;

  // Set style
  utils.set_style();

  //TFile * fin = new TFile( argv[1], "read" );
  
  
  //std::vector< ScanPoint > scanpoints = ReadScanPoints( fin );

  //vector< double > xbins = utils.get_bins( scanpoints, 'x' ); // how to use the template correctly
  //vector< double > ybins = utils.get_bins( scanpoints, 'y' );
  //TH2D * pmt0_qe = new TH2D("pmt0_qe", "Detection efficiency", xbins.size()-1, &xbins[0], ybins.size()-1, &ybins[0]);
  
  vector<int> phidgets = {0, 1, 3, 4};
  vector<PTF::PMTChannel> activeChannels = {};
 
  PTF::Wrapper wrapper = PTF::Wrapper(16384, 70, activeChannels, phidgets);  
  
  wrapper.openFile(argv[1], "scan_tree");
  
  vector <double> temperature_int;
  vector <double> temperature_ext1;
  vector <double> temperature_ext2;
  vector <double> Time; 
  
  for (unsigned int i = 0; i < wrapper.getNumEntries(); i++) {
        wrapper.setCurrentEntry(i);
        //auto location = wrapper.getDataForCurrentEntry(PTF::Gantry1);                                                                                                      
        auto reading  = wrapper.getReadingTemperature();
        auto time_value=wrapper.getReadingTime();
		int n1,n2;
		n1=wrapper.getCurrentEntry();
		n2=wrapper.getNumEntries();
		std::cout << "Num entries: " <<  n1 <<"  "<<n2 <<"\t\r" << std::flush;
		
        wrapper.setCurrentEntry(0);
        auto time_s=wrapper.getReadingTime();
	//temperature_int.push_back(reading.int_1);
        //temperature_ext1.push_back(reading.ext_1);
        temperature_ext2.push_back(reading.ext_2);
        Time.push_back(time_value.time_c-time_s.time_c);


   
   }
  //TGraph* temperature_data=new TGraph(Time.size(),&Time[0],&temperature_int[0]);
  //TGraph* temperature_data_2=new TGraph(Time.size(),&Time[0],&temperature_ext1[0]);		
   TGraph* temperature_data_3=new TGraph(Time.size(),&Time[0],&temperature_ext2[0]);
   temperature_data_3->GetYaxis()->SetTitle("Temperature _reading");//Create normal plot for all the corrected efficiency                                                                     
   temperature_data_3->GetXaxis()->SetTitle("Time (sec)");
   //temperature_data->SetTitle("Temperature int (degC)");
   //temperature_data_2->SetTitle("Temperature ext1 (degC)");
   temperature_data_3->SetTitle("Temperature ext2 (degC)");
   //temperature_data->SetMarkerColor(2);
   //temperature_data_2->SetMarkerColor(8);
   temperature_data_3->SetMarkerColor(8);
   temperature_data_3->SetMarkerStyle(1);
     //pmt_correlation->GetXaxis()->SetRangeUser(0.1,0.3);
   temperature_data_3->GetYaxis()->SetRangeUser(18,24);


   TCanvas* c = new TCanvas("canvas");
   string plotname;
   //temperature_data->Draw("ap1");
   //temperature_data_2->Draw("p1 same");
   temperature_data_3->Draw("ap1");
   
//auto leg = new TLegend(0.2,0.5,0.5,0.7);
// leg->SetFillColor(2);
//  leg->AddEntry("temperature_data","Int","p");
//  leg->AddEntry("temperature_data_2","ext_1","p");
//  leg->AddEntry("temperature_data_3","ext_2","p");
//  leg->Draw();
   c->BuildLegend();
   gPad->Modified();
   gPad->Update();
   plotname = string("~/projects/def-pdeperio/vgousy/ptf/ptf-analysis-2/plot_dir/1par/temperature_data")+argv[2]+".pdf";
   c->SaveAs(plotname.c_str(),"pdf");
  } 
