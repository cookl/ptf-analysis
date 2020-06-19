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
   
  TFile * fin = new TFile( argv[1], "read" );
  //string outname = string("output_0000") + argv[2] + ".root";
  //TFile * outFile = new TFile(outname.c_str(), "NEW");

  vector<int> phidgets = {0, 1, 3, 4};
  vector<PTF::PMTChannel> activeChannels = {};
 
  PTF::Wrapper wrapper = PTF::Wrapper(16384, 70, activeChannels, phidgets);

  unordered_set<int> skipLines = {};// {962,1923,2884,5240,6201,9611,10572,11533,12494,13455,15811,16771};

  wrapper.openFile(argv[1], "scan_tree");

  cerr << "Num entries: " << wrapper.getNumEntries() << endl << endl;
  
  for (unsigned int i = 0; i < wrapper.getNumEntries(); i++) {

      auto reading  = wrapper.getReadingTemperature(PTF::T);

      cerr << reading.int_1 << reading.ext_1 << reading.ext_2<<endl;
      
  }
}

