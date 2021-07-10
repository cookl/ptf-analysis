#include "wrapper.hpp"
#include "ScanPoint.hpp"
#include "Utilities.hpp"
#include <iostream>
#include <sstream>

#include "TCanvas.h"
#include <TPaveStats.h>
#include "TStyle.h"
#include "TH1D.h"
using namespace std;


/*
This program takes a mPMT ROOT data file and a specified event number as inputs. 

Other optional inputs include:
  - specified channels (0-20) provided as a comma separated list
  - specified time range in ns (0-8192) in the format (lower_range,higher_range)
Channel and time range options can be specified separately or at the same time 
(order does not matter).

The default settings are:
  - channels 0-20
  - full time range 0-8192 ns

Each waveform for the specified (or otherwise default) channels is outputted 
separately as a png file. It is also collectively stored and outputted in a 
ROOT file under 'runnum_wfdisplay_eventnum.root'.

Additionally, an overlay plot of the waveforms of all the specified (or 
otherwise default) channels are also outputted.
*/

TH1D* hwaveform{nullptr}; // current waveform
int color[20]={8,5,41,28,9,30,46,38,1,920,632,416,600,400,616,432,800,820,880,860};

int main(int argc, char** argv) {  

  // Error and usage message
  if (argc<3 || argc==4 || argc==6 || argc>7) {
    cerr << "Give path to file to read, and event number:\n" << 
            "   usage: './bin/waveform_plotting.app output00000838.root 510940'\n" << 
            "Optionally, provide channels (0-19) in a comma separated list:\n" <<
            "   usage: './bin/waveform_plotting.app output00000838.root 510940 --channels 4,9,12'\n" <<
            "Optionally, provide time range in ns(0-8192) to zoom in on in the format (lower_range,higher_range):\n" <<
            "   usage: './bin/waveform_plotting.app output00000838.root 510940 --time 1290,1320'\n" <<
            "Or combine --channel and --time options." << endl;
    return 0;
  }
 
 // Declare variables
  bool specifyTime = false;   // time range selection option
  bool specifyCh = false;     // channel selection option

  int time[2] = {0,8192};     // default time range
  int chan[20] = {0};         // default channels to read from, given specifyCh==true

  TCanvas *c2 = new TCanvas("overlay", "Overlay",1000,800);

  // Extract run number, time window range, and active channels based on command line inputs
  string run_num = string(argv[1]).substr(11,3);
  unsigned event_num = 2 + atoi(argv[2]);       //offset of 2 included to remain consistent with PTFAnalysis.cpp line 661
  
  for (int arg=3; arg<argc; arg=arg+2){

    // Extract specified channels (if any)
    if (string(argv[arg])=="--channels" && !specifyCh) {
      string str_ch = string(argv[arg+1]);
      stringstream ss_ch(str_ch);
      while (ss_ch.good()) {
        string substr_ch;
        getline(ss_ch, substr_ch, ',');
        chan[stoi(substr_ch)]=1;
      }
      specifyCh = true;
    }

    // Extract specified time range (if any)
    if (string(argv[arg])=="--time" && !specifyTime) {
      string str_time = string(argv[arg+1]);
      stringstream ss_time(str_time);
      int index = 0;
      while (ss_time.good()) {
        string substr_time;
        getline(ss_time, substr_time, ',');
        time[index] = stoi(substr_time);
        index++;
        if (index>1) break;
      }
      if (time[0]>time[1]) cerr << "Please input time in the format (lower_range, higher_range)." << endl;
      specifyTime = true;
    }
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

  if (specifyCh) {
    for(int i = 0; i < 20; i++){ 
      if (chan[i]==1) {
        PTF::PMT PMT = {i,i,PTF::mPMT_REV0_PMT};
        activePMTs.push_back(PMT);
      }      
    }
  } else {
    for(int i = 0; i < 20; i++){
      PTF::PMT PMT = {i,i,PTF::mPMT_REV0_PMT};
      activePMTs.push_back(PMT);
    }
  }
  
  vector<int> phidgets = {0, 1, 3};
  vector<PTF::Gantry> gantries = {PTF::Gantry0, PTF::Gantry1};
  Wrapper wrapper = Wrapper(1, 1024, activePMTs, phidgets, gantries,mPMT_DIGITIZER);
  wrapper.openFile( string(argv[1]), "scan_tree");

  // Retrieve waveform display for specified event number in each channel
  for(unsigned int i = 0; i < activePMTs.size(); i++){
    PTF::PMT pmt = activePMTs[i];
    // Get digitizer settings
    Digitizer digi = wrapper.getDigitizerSettings();
    double digiCounts = pow(2.0, digi.resolution);

    // Get length of waveforms
    wrapper.setCurrentEntry(0);
    int  numTimeBins= wrapper.getSampleLength();
    
    // Build the waveform histogram
    string hname = "run0" + run_num + " - event" + to_string(event_num-2) + " - ch" + to_string(pmt.channel);
    string htitle = hname + "; Time (ns); Voltage (V)";
    outFile->cd();
    hwaveform = new TH1D( hname.c_str(), htitle.c_str(), numTimeBins, 0., float(numTimeBins)*1000/digi.samplingRate );

    wrapper.setCurrentEntry(event_num);
    double* pmtsample=wrapper.getPmtSample( pmt.pmt, 0 );

    // Set the contents of the histogram
    hwaveform->Reset();
    for ( int ibin=1; ibin <= numTimeBins; ++ibin ){
      hwaveform->SetBinContent( ibin, pmtsample[ibin-1] );
      hwaveform->SetBinError( ibin, 2.1e-3 );
    }
    hwaveform->Scale(digi.fullScaleRange/digiCounts);

    // Print waveform separately
    TCanvas *c1 = new TCanvas("c1", "C1",1000,800);
    c1->cd();
    hwaveform->Draw();
    if (specifyTime) hwaveform->GetXaxis()->SetRangeUser(time[0], time[1]);
    string filename= hname + ".png";
    c1->SaveAs(filename.c_str());

    // Overlay current waveform
    c2->cd();
    hwaveform->SetLineColor(color[i]);
    hwaveform->SetMarkerStyle(7);
    hwaveform->SetMarkerColor(color[i]);
    hwaveform->Draw("SAMES");
    c2->Update();
  }

  // Customize and print overlay plot of all specified/default channels
  string title = "run0" + run_num + " - event" + to_string(event_num-2) + " - overlaid";
  string filename2 = "run0" + run_num + " - event" + to_string(event_num-2) + " - overlay.png";
  TPaveText* pt  = (TPaveText*)gPad->FindObject("title");
  TPaveText* ptn = new TPaveText(pt->GetX1(),pt->GetY1(), pt->GetX2(),pt->GetY2(),"");
  ptn->AddText(title.c_str());
  ptn->SetBorderSize(0);
  ptn->SetFillColor(0);
  ptn->SetTextFont(42);
  ptn->Draw();
  gPad->BuildLegend();
  c2->SaveAs(filename2.c_str());

  outFile->Write();
  outFile->Close();
  return 0;
}


