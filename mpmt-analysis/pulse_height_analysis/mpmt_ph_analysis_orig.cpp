#include "WaveformFitResult.hpp"
#include "ScanPoint.hpp"
#include "TFile.h"
#include "TCanvas.h"
#include "TH1D.h"
#include "TH2F.h"
#include "TGraph.h"
#include "TF1.h"
#include "TStyle.h"
#include "TLegend.h"
#include "TProfile.h"

#include <iostream>
#include <vector>



int main( int argc, char* argv[] ) {

  if ( argc != 2 ){
    std::cerr<<"Usage: ptf_ttree_analysis.app ptf_analysis.root\n";
    exit(0);
  }

  TFile * fin = new TFile( argv[1], "read" );


  // get the waveform fit TTree
  std::cout << "TTree 0" << std::endl;
  TTree * tt0 = (TTree*)fin->Get("ptfanalysis0");
  WaveformFitResult * wf0 = new WaveformFitResult;
  if(tt0) wf0->SetBranchAddresses( tt0 );


  
  TH1F *pulse_height     = new TH1F("PH0","Pulse Height ",200,0,0.000488*200*1000);


  // Loop the first scan point and print something 
  for(int i = 0; i < tt0->GetEntries()-1; i++){

    
    tt0->GetEvent(i );

    // Loop over the pulses we found

    
    for(int k = 0; k < wf0->numPulses; k++){
      std::cout << "Pulse found. Time= " << wf0->pulseTimes[k]
                << " Pulse Height = " << wf0->pulseCharges[k] << std::endl;
      pulse_height->Fill(wf0->pulseCharges[k] * 1000.0);  // Convert to mV

      
    }

  } 


      
  TCanvas *c = new TCanvas("C");
  pulse_height->Draw();
    
  c->SaveAs("mpmt_pulse_height.png");



  
  fin->Close();
  return 0;
}

