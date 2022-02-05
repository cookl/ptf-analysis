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


using namespace std;


int main(int argc, char** argv) {
	
	
  if (argc != 3) {
      cerr << "give path to file to read" << endl;
      cerr << "usage: temperature_reading filename.root run_number" << endl;
      return 0;
    }

    // Get utilities
  Utilities utils;

    // Set style
  utils.set_style();
   
  //TFile * fin = new TFile( argv[1], "read" );
  //string outname = string("output_0000") + argv[2] + ".root";
  //TFile * outFile = new TFile(outname.c_str(), "NEW");

  vector<int> phidgets = {0, 1, 3, 4};
  vector<PTF::PMT> activePMTs = {};
  vector<PTF::Gantry> gantries = {PTF::Gantry0, PTF::Gantry1}; 
  Wrapper wrapper = Wrapper(16384, 70, activePMTs, phidgets, gantries, PTF_CAEN_V1730);

  wrapper.openFile(argv[1], "scan_tree");

  cerr << "Num entries: " << wrapper.getNumEntries() << endl << endl;
  
  for (unsigned int i = 0; i < wrapper.getNumEntries(); i++) {
      wrapper.setCurrentEntry(i);
      //auto location = wrapper.getDataForCurrentEntry(PTF::Gantry1);
      auto reading  = wrapper.getReadingTemperature();
      auto time_value=wrapper.getReadingTime();
      wrapper.setCurrentEntry(0);
      auto time_s=wrapper.getReadingTime();
      auto reading_2  = wrapper.getReadingTemperature();
      cout <<wrapper.getCurrentEntry()<<endl;
      //cerr << location.x << " " << location.y << " " << location.z << " "<<endl;
      cerr <<" "<<reading.ext_2<< " "<<time_value.time_c-time_s.time_c <<endl;
      //cout<<setprecision(5)<<time_value.time_c<<endl;
  }
}


