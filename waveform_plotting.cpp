#include "wrapper.hpp"
#include "ScanPoint.hpp"
#include "Utilities.hpp"
#include <iostream>

#include "TCanvas.h"
#include "TH1D.h"

using namespace std;

TH1D* hwaveform{nullptr}; // current waveform

int main(int argc, char** argv) {  
  // Last two arguments required:
  //  - zoom_range: use "-1" to view full window or specify a single time value (ns) 
  //  - channel: use "-1" for channels 4-15 or specify a single channel 0-16
  if (argc < 5) {
    if (argc < 2) cerr << "Give path to file to read" << endl;
    cerr << "usage: waveform_plotting filename.root event_num zoom_range channels" << endl;
    return 0;
  }
  
  // Extract run number, time window range, and active channels
  string run_num = string(argv[1]).substr(11,3);

  int zoom = -1;
  if (argc>3) zoom = atoi(argv[3]); 

  int chan[17]= {0};
  if (argc>4 && atoi(argv[4]) != -1) {
    for (int k=4; k<argc; k++) chan[atoi(argv[k])] = 1;
  } 

  // Get utilities & set syle
  Utilities utils;
  utils.set_style();

  // Opening the output root file
  string outname = "run0" + run_num + "_wfdisplay_" + argv[2] + ".root";
  TFile * outFile = new TFile(outname.c_str(), "NEW");

  // Retrieve active channels and set up PMT
  vector<int> active_channels;
  vector<PTF::PMT> activePMTs;

  if (argc>4 && atoi(argv[4])!=-1) {
    for(int i = 0; i < 17; i++){ 
      if (chan[i]==1) {
        PTF::PMT PMT = {i,i,PTF::mPMT_REV0_PMT};
        activePMTs.push_back(PMT);
      }      
    }
  } else {
    for(int i = 4; i <= 15; i++){
      PTF::PMT PMT = {i,i,PTF::mPMT_REV0_PMT};
      activePMTs.push_back(PMT);
    }
  }
  
  vector<int> phidgets = {0, 1, 3};
  vector<PTF::Gantry> gantries = {PTF::Gantry0, PTF::Gantry1};
  Wrapper wrapper = Wrapper(1, 1024, activePMTs, phidgets, gantries,mPMT_DIGITIZER);
  wrapper.openFile( string(argv[1]), "scan_tree");

  // Retrieve waveform display for specified event number in each channel
  unsigned event_num = 2 + atoi(argv[2]);

  for(unsigned int i = 0; i < activePMTs.size(); i++){
    PTF::PMT pmt = activePMTs[i];
    // Get digitizer settings
    Digitizer digi = wrapper.getDigitizerSettings();
    double digiCounts = pow(2.0, digi.resolution);

    // Get length of waveforms
    wrapper.setCurrentEntry(0);
    int  numTimeBins= wrapper.getSampleLength();
    
    // Build the waveform histogram
    string hname = "run0" + run_num + " ch" + to_string(pmt.channel) + " - " + to_string(event_num-2);
    string htitle = hname + "; Time (ns); Voltage (V)";
    outFile->cd();
    hwaveform = new TH1D( hname.c_str(), htitle.c_str(), numTimeBins, 0., float(numTimeBins)*1000/digi.samplingRate );

    wrapper.setCurrentEntry(event_num);
    double* pmtsample=wrapper.getPmtSample( pmt.pmt, 0 );

    // set the contents of the histogram
    hwaveform->Reset();
    for ( int ibin=1; ibin <= numTimeBins; ++ibin ){
      hwaveform->SetBinContent( ibin, pmtsample[ibin-1] );
      hwaveform->SetBinError( ibin, 2.1e-3 );
    }
    hwaveform->Scale(digi.fullScaleRange/digiCounts);

    // Print zoomed in waveform display image
    TCanvas *c1 = new TCanvas("c1", "C1", 800,600);
    hwaveform->Draw();
    if (zoom!=-1) hwaveform->GetXaxis()->SetRangeUser(zoom-300, zoom+300);
    string filename= hname + ".png";
    c1->SaveAs(filename.c_str());
  }

  outFile->Write();
  outFile->Close();
  return 0;
}


