#include "WaveformFitResult.hpp"
#include "ScanPoint.hpp"
#include "Utilities.hpp"
#include "FindCircle.hpp"
#include "TDirectory.h"
#include "TFile.h"
#include "TH2D.h"
#include "TStyle.h"
#include "TCanvas.h"
#include "TGraph.h"//scatter plot
#include "TMinuit.h"//minimisation 
#include "TVirtualFitter.h"
#include <iostream>
#include <vector>
#include <algorithm>
#include <string>
#include <fstream>


int main( int argc, char* argv[] ) {

  if ( argc != 3 ){
    std::cerr<<"Usage: ptf_qe_analysis.app ptf_analysis.root run_number\n";
    exit(0);
  }

  // Get utilities
  Utilities utils;

  // Set style
  utils.set_style();

  TFile * fin = new TFile( argv[1], "read" );
  
  
  std::vector< ScanPoint > scanpoints = ReadScanPoints( fin );

  vector< double > xbins = utils.get_bins( scanpoints, 'x' ); // how to use the template correctly
  vector< double > ybins = utils.get_bins( scanpoints, 'y' );
  TH2D * pmt0_qe = new TH2D("pmt0_qe", "Detection efficiency", xbins.size()-1, &xbins[0], ybins.size()-1, &ybins[0]);
  
  
  
  wrapper.openFile(argv[1], "scanpoints");
  
  vector <double> temperature; 
  vector <double> Time; 
  
  for (unsigned int i = 0; i < wrapper.getNumEntries(); i++) {
        wrapper.setCurrentEntry(i);
        //auto location = wrapper.getDataForCurrentEntry(PTF::Gantry1);                                                                                                      
        auto reading  = wrapper.getReadingTemperature();
        auto time_value=wrapper.getReadingTime();
        wrapper.setCurrentEntry(0);
        auto time_s=wrapper.getReadingTime();
		temperature.push_back(reading.int_1);
		Time.push_back(time_value.time_c-time_s.time_c)
	}
   
   TGraph* temperature_data=new TGraph(temperature.size(),&temperature[0],&Time[0]);
		
   temperature_data->GetXaxis()->SetTitle("Temperature _reading");//Create normal plot for all the corrected efficiency                                                                     
   temperature_data->GetYaxis()->SetTitle("Time");
   temperature_data->SetTitle("Temperature");
   temperature_data->SetMarkerStyle(1);
   temperature_data->SetMarkerColor(2);
	    //pmt_correlation->GetXaxis()->SetRangeUser(0.1,0.3);
	    //pmt_correlation->GetYaxis()->SetRangeUser(0.4,0.6);
  
   temperature_data->Draw("colz0");
   gPad->Modified();
   gPad->Update();
   plotname = string("~/projects/def-pdeperio/vgousy/ptf/ptf-analysis-2/plot_dir/1par/temperature_data")+argv[2]+".pdf";
   c->SaveAs(plotname.c_str(),"pdf");
}  