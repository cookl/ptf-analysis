#include "wrapper.hpp"
#include "ScanPoint.hpp"
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

#include "TH1D.h"
// #include "TH2D.h"
#include "TFile.h"
#include "TF1.h"
#include "TMath.h"
#include "Configuration.hpp"

using namespace std;

TH1D* hwaveform{nullptr}; // current waveform

int main(int argc, char** argv) {

  if (argc != 5) {
    if (argc < 2) cerr << "Give path to file to read" << endl;
    cerr << "Give waveform number to display" << endl;
    cerr << "usage: waveform_plotting filename.root run_number config_file event_num" << endl;
    return 0;
  }

  // Get utilities & set syle
  Utilities utils;
  utils.set_style();

  // Opening the output root file
  string outname = string("run0") + argv[2] + "_wfdisplay_" + argv[4] + ".root";
  TFile * outFile = new TFile(outname.c_str(), "NEW");

  // Open config file and retrieve active channels
  Configuration config;
  config.Load(argv[3]);
  vector<int> active_channels;
  if( !config.Get("mpmt_channel_list", active_channels) ){
    cout << "Missing terminal_output parameter from config file." << endl;
    exit( EXIT_FAILURE );
  }

  // Print active channels
  cout << "Number of active channels: " << active_channels.size() << "\nActive channels are: ";
  for(unsigned int i = 0; i < active_channels.size(); i++){ cout << active_channels[i]<< " ";}
  cout << endl;

  // Loop over the active channels to do setup.
  vector<int> phidgets = {0, 1, 3};
  vector<PTF::PMT> activePMTs;
  for(unsigned int i = 0; i < active_channels.size(); i++){ 
    int ch = active_channels[i];
    PTF::PMT PMT = {ch,ch,PTF::mPMT_REV0_PMT};
    activePMTs.push_back(PMT);
  }
  vector<PTF::Gantry> gantries = {PTF::Gantry0, PTF::Gantry1};
  Wrapper wrapper = Wrapper(1, 1024, activePMTs, phidgets, gantries,mPMT_DIGITIZER);
  wrapper.openFile( string(argv[1]), "scan_tree");

  // Get waveform display for specified event in each channel
  unsigned event_num = 2 + atoi(argv[4]);

  for(unsigned int i = 0; i < active_channels.size(); i++){
    PTF::PMT pmt = activePMTs[i];
    // Get digitizer settings
    Digitizer digi = wrapper.getDigitizerSettings();
    double digiCounts = pow(2.0, digi.resolution);

    // Get length of waveforms
    wrapper.setCurrentEntry(0);
    int  numTimeBins= wrapper.getSampleLength();
    
    // Build the waveform histogram
    std::string hname =  "Chan " + std::to_string(pmt.channel) + " - Event " + to_string(event_num-2);
    std::string htitle = hname + "; Time (ns); Voltage (V)";
    outFile->cd();
    hwaveform = new TH1D( hname.c_str(), htitle.c_str(), numTimeBins, 0., float(numTimeBins)*1000/digi.samplingRate );
    outFile->cd();
    wrapper.setCurrentEntry(event_num);
    
    int numWaveforms = wrapper.getNumSamples();
    // for ( int j=0; j<numWaveforms; j++) { 
      double* pmtsample=wrapper.getPmtSample( pmt.pmt, 0 );
      // set the contents of the histogram
      hwaveform->Reset();
      for ( int ibin=1; ibin <= numTimeBins; ++ibin ){
        hwaveform->SetBinContent( ibin, pmtsample[ibin-1] );
        hwaveform->SetBinError( ibin, 2.1e-3 );
      }
      hwaveform->Scale(digi.fullScaleRange/digiCounts);
    // }
  }

  outFile->Write();
  outFile->Close();
  return 0;
}


