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
  TTree * tt = (TTree*)fin->Get("ptfanalysis0");
  WaveformFitResult * wf = new WaveformFitResult;
  wf->SetBranchAddresses( tt );
  
  std::cout << "Looping tree " << std::endl;
  // Loop the first scan point and print something
  for(int i = 0; i < tt->GetEntries(); i++){
    tt->GetEvent(i );
    std::cout  <<" number of found pulses="<<wf->numPulses << std::endl;
    if(wf->numPulses > 0){
      for(int i = 0; i < wf->numPulses; i++){
        std::cout <<" pulse " << i << " has time = "<<wf->pulseTimes[i]
                  << " ns " <<  std::endl;
      }

    }

  }
      
  return 0;
}
