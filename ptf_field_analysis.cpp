#include "wrapper.hpp"
#include "Utilities.hpp"

#include "TFile.h"
#include "TGraph.h"
#include "TMultiGraph.h"
#include "TStyle.h"
#include "TCanvas.h"

#include <string>
#include <iostream>
#include <ostream>
#include <fstream>
#include <unordered_set>
#include <math.h>

using namespace std;

int main(int argc, char** argv) {
  if (argc != 3) {
    cerr << "give path to file to read" << endl;
    cerr << "usage: ptf_field_analysis.app /data/directory run_number" << endl;
    return 0;
  }

  // Get utilities
  Utilities utils;

  // Set style
  utils.set_style();

  // Opening the output root file
  string outname = string("ptf_field_analysis_run0") + argv[2] + ".root";
  TFile * fout = new TFile(outname.c_str(), "NEW");
  //TFile * outFile = new TFile("ptf_analysis.root", "NEW");
  
  // Create the output graphs
  TMultiGraph *mg = new TMultiGraph();
  TGraph *gr_bx = new TGraph();
  TGraph *gr_by = new TGraph();
  TGraph *gr_bz = new TGraph();

  // Set up PTF Wrapper
  vector<int> phidgets = {4};
  vector<PTF::PMT> activePMTs = {}; //Not looking at PMT data
  vector<PTF::Gantry> gantries = {PTF::Gantry0, PTF::Gantry1};
  PTF::Wrapper wrapper = PTF::Wrapper(6000, 70, activePMTs, phidgets, gantries);
  wrapper.openFile( string(argv[1])+"/out_run0"+argv[2]+".root", "scan_tree");
  cerr << "Num entries: " << wrapper.getNumEntries() << endl << endl;

  unordered_set<int> skipLines = {};// {962,1923,2884,5240,6201,9611,10572,11533,12494,13455,15811,16771};

  uint32_t lines = 0;
  const uint32_t freq = 100;
  uint32_t gr_point = 0;
  uint32_t gr_freq = 0;
  gr_freq = ceil( (float)wrapper.getNumEntries() / 1000.0 );

  for (unsigned int i = 0; i < wrapper.getNumEntries(); i++) {
    // cerr << "Entry " << i;
    if (i % freq == 0 || i == wrapper.getNumEntries() - 1) {
      cerr << "Entry " << i << "/" << wrapper.getNumEntries() << "\u001b[34;1m (" << (((double)i)/wrapper.getNumEntries()*100) << "%)\u001b[0m\033[K";
      if (skipLines.find(i) != skipLines.end()) {
        cerr << "\u001b[31;1m Skipping...\u001b[0m\r";
        continue;
      } else {
        cerr << "\r";
      }
    }

    if (skipLines.find(i) != skipLines.end()) continue;
    lines++;

    wrapper.setCurrentEntry(i);

    //auto location = wrapper.getDataForCurrentEntry(PTF::Gantry1);

    if (i % gr_freq == 0){
      for (int phidget : phidgets) {
        auto reading  = wrapper.getReadingForPhidget(phidget);
        gr_bx->SetPoint(gr_point, gr_point, reading.Bx[0]);
        gr_by->SetPoint(gr_point, gr_point, reading.By[0]);
        gr_bz->SetPoint(gr_point, gr_point, reading.Bz[0]);
      }
      gr_point++;
    }

  }

  mg->Add(gr_bx); gr_bx->SetTitle("Bx"); gr_bx->SetLineWidth(3);
  mg->Add(gr_by); gr_by->SetTitle("By"); gr_by->SetLineWidth(3);
  mg->Add(gr_bz); gr_bz->SetTitle("Bz"); gr_bz->SetLineWidth(3);

  mg->SetTitle("Magnetic field strength; Data point; Field strength [G]");

  //Set directories
  //gr_bx->SetDirectory( fout );
  //gr_by->SetDirectory( fout );
  //gr_bz->SetDirectory( fout );
  //mg->SetDirectory( fout );

  //Make plots
  TCanvas* c = new TCanvas("canvas");
  mg->Draw("A pmc plc");
  c->BuildLegend();
  string plotname;
  plotname = string("ptf_field_analysis_run0")+argv[2]+".pdf";
  c->SaveAs(plotname.c_str(),"pdf");

  //Write and close output file
  fout->cd();
  gr_bx->Write();
  gr_by->Write();
  gr_bz->Write();
  mg->Write();
  fout->Write();
  fout->Close();

  cout << "Done. Number of data points: " << lines << endl;

  return 0;
}

