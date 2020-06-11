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

  // get the scanpoints information
  std::vector< ScanPoint > scanpoints = ReadScanPoints( fin );
  for ( const ScanPoint& sp : scanpoints ){
    std::cout << sp;
  }

  // get the waveform fit TTree
  TTree * tt = (TTree*)fin->Get("ptfanalysis");
  WaveformFitResult * wf = new WaveformFitResult;
  wf->SetBranchAddresses( tt );

  // loop over the first scan point and print something
  for ( unsigned iev = 0; iev < scanpoints[0].nentries(); ++iev ){
    tt->GetEvent( scanpoints[0].get_entry() + iev );
    std::cout <<"iev = "<<iev
	      <<" scanpt="<<wf->scanpt
	      <<" amp="<<wf->amp
	      <<" iswf="<<wf->haswf<<std::endl;

  }
      
  return 0;
}
