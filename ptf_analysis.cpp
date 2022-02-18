/// Analysis of PTF data
///
/// Definitions of terminology used, and associated structures:
///
/// Each input root file has a number of "ScanPoint" events, where each event has
/// several waveforms at a given location.
///
/// class  MeanRMSCalc           For calculating mean and rms
///
/// class  PTFErrorBarAnalysis   For calculating the error bar size to use on the waveforms
///
/// struct WaveformFitResult     Structure to hold one waveform fit result
///
/// class  ScanPoint             Holds location of scan point, first entry number in TTree of scan point, and number of waveforms
///
/// class  PTFAnalysis           For doing analysis of all of the waveforms, and keep track of scan points, store results in TTree
///
///
/// This program takes a PTF scan root file as input.
/// The analysis proceeds in these steps:
/// 1. Determine the uncertainties on the "collected charge" by studying the pedestal
///    of the waveforms in bins 1-20.  This part of the analysis produces two kinds of
///    pedestal histograms.
/// 2. Fit the waveforms to gaussians.  Each waveform is read in once into a Waveform object,
///    it is then fit to a gaussian using a fitter.  Fit results are stored in a TTee that
///    has one entry per scan-point
/// 3. The now filled TTree of fitted waveform parameters is analyzed to make histograms
///     
/// Author: Blair Jamieson (Sep. 2019) 


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
  if (argc != 4) {
    cerr << "give path to file to read" << endl;
    cerr << "usage: ptf_analysis filename.root run_number config_file" << endl;
    return 0;
  }

  // Get utilities
  Utilities utils;

  // Set style
  utils.set_style();

  // Opening the output root file
  string outname = string("ptf_analysis_run0") + argv[2] + ".root";
  TFile * outFile = new TFile(outname.c_str(), "NEW");
  //TFile * outFile = new TFile("ptf_analysis.root", "NEW");

  // Set up PTF Wrapper
  vector<int> phidgets = {0, 1, 3};
  PTF::PMT PMT0 = {0,0,PTF::Hamamatsu_R3600_PMT}; // only looking at one PMT at a time
  PTF::PMT PMT1 = {1,5,PTF::PTF_Monitor_PMT}; // only looking at one PMT at a time
  PTF::PMT REF = {2,1,PTF::Reference}; // only looking at one PMT at a time
  vector<PTF::PMT> activePMTs = { PMT0, PMT1, REF }; // must be ordered {main,monitor}
  vector<PTF::Gantry> gantries = {PTF::Gantry0, PTF::Gantry1};
  Wrapper wrapper = Wrapper(6000, 70, activePMTs, phidgets, gantries, PTF_CAEN_V1730);
  wrapper.openFile( string(argv[1]), "scan_tree");
  cerr << "Num entries: " << wrapper.getNumEntries() << endl << endl;

  // Determine error bars to use on waveforms
  // Commented to save time
  // Approximate error bars for waveforms are sufficient for fitting
  //ErrorBarAnalysis * errbars0 = new ErrorBarAnalysis( outFile, wrapper, PMT0 );
  //std::cout << "Using PMT0 errorbar size " << errbars0->get_errorbar() << std::endl;
  //ErrorBarAnalysis * errbars1 = new ErrorBarAnalysis( outFile, wrapper, PMT1 );
  //std::cout << "Using PMT1 errorbar size " << errbars1->get_errorbar() << std::endl;
  
  // Do analysis of waveforms for each scanpoint
  PTFAnalysis *analysis0 = new PTFAnalysis( outFile, wrapper, 4.4/*errbars0->get_errorbar()*/, PMT0, string(argv[3]), true );
  analysis0->write_scanpoints();

  // Switch PMT to monitor PMT
  
  // Do analysis of waveforms for each scanpoint
  PTFAnalysis *analysis1 = new PTFAnalysis( outFile, wrapper, 4.4/*errbars1->get_errorbar()*/, PMT1, string(argv[3]), true );

  // Switch to reference waveform
  
  // Do analysis of waveforms for each scanpoint
  PTFAnalysis *analysis2 = new PTFAnalysis( outFile, wrapper, 4.4/*errbars2->get_errorbar()*/, REF, string(argv[3]), true );
  
  // Do quantum efficiency analysis
  // This is now also done in a separate analysis script (including temperature corrections)
  //PTFQEAnalysis *qeanalysis = new PTFQEAnalysis( outFile, analysis0, analysis1 );

  outFile->Write();
  outFile->Close();
    
  cout << "Done" << endl; 

  return 0;
}

