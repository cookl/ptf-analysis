#include "WaveformFitResult.hpp"
#include "ScanPoint.hpp"
#include "TFile.h"

#include <iostream>
#include <vector>


int main( int argc, char* argv[] ) {

  if ( argc != 2 ){
    std::cerr<<"Usage: ptf_ttree_analysis.app ptf_analysis.root\n";
    exit(0);
  }

  TFile * fin = new TFile( argv[1], "read" );


  // get the waveform fit TTree
  TTree * tt0 = (TTree*)fin->Get("ptfanalysis0");
  WaveformFitResult * wf0 = new WaveformFitResult;
  wf0->SetBranchAddresses( tt0 );
  TTree * tt1 = (TTree*)fin->Get("ptfanalysis1");
  WaveformFitResult * wf1 = new WaveformFitResult;
  wf1->SetBranchAddresses( tt1 );
  
  std::cout << "Looping tree " << tt0->GetEntries() << " " << tt1->GetEntries() << std::endl;
  // Loop the first scan point and print something
  for(int i = 0; i < tt0->GetEntries(); i++){
    tt0->GetEvent(i );
    tt1->GetEvent(i );
    std::cout  <<" number of found pulses="<<wf1->numPulses << std::endl;

  }
      
  return 0;
}
