#include "ErrorBarAnalysis.hpp"

using namespace std;

ErrorBarAnalysis::ErrorBarAnalysis( TFile* outfile, PTF::Wrapper & wrapper, PTF::PMT & pmt ){
  static unsigned instance_count = 0;
  ++instance_count;
  
  // Do two error bar analyses:
  // one for absolute pedestal value, and one for difference between first two time bins
  MeanRMSCalc absrmscalc;
  MeanRMSCalc diffrmscalc;

  // First lets determine mean and rms to use for defining histogram ranges
  // Loop over scan points (index i)
  for (unsigned long long i = 0; i < wrapper.getNumEntries(); ++i) {
    if ( i>100 ) continue;
    std::cout<< "ErrorBarAnalysis First Pass scan point " << i << " / " << wrapper.getNumEntries() << std::endl;
    wrapper.setCurrentEntry(i);
    //get number of samples
    int  numSamples= wrapper.getNumSamples();

    // loop over the number of waveforms at this ScanPoint (index j) 
    for ( int j=0; j<numSamples; ++j ) {
      //if ( j>50 ) continue;
      double* pmtsample=wrapper.getPmtSample( pmt.pmt, j );
      diffrmscalc.add( pmtsample[ 1 ] - pmtsample[ 0 ] );
      // only use first 20 time-bins of the waveform (time-bin index k)
      for ( int k = 0; k < std::min( 20, wrapper.getSampleLength() ); ++k ){
        absrmscalc.add( pmtsample[k] );
      }
    }
  }

  cout<<"ABS RMS "<<absrmscalc;
  cout<<"DIF RMS "<<diffrmscalc;
  
  // Build the histograms using known mean and rms
  std::ostringstream hnabs, htabs, hndiff, htdiff;
  hnabs << "habs_" << instance_count;
  htabs << "Absolute Pedestal " << instance_count << " ; Charge collected ; count/bin";
  hndiff << "hdiff_" << instance_count;
  htdiff << "Differential Pedestal " << instance_count << " ; #Delta Charge collected ; count/bin";
  outfile->cd();
  habs = new TH1D( hnabs.str().c_str(), htabs.str().c_str(), 50, absrmscalc.mean()-4*absrmscalc.rms(), absrmscalc.mean()+4*absrmscalc.rms() ); 
  hdiff = new TH1D( hndiff.str().c_str(), htdiff.str().c_str(), 50, diffrmscalc.mean()-4*diffrmscalc.rms(), diffrmscalc.mean()+4*diffrmscalc.rms() ); 

  // Now fill the histograms
  // Loop over scan points (index i)
  for (unsigned long long i = 0; i < wrapper.getNumEntries(); ++i ) {
    //if ( i>100 ) continue;
    std::cout << "ErrorBarAnalysis Second Pass scan point " << i << " / " << wrapper.getNumEntries() << std::endl;
    wrapper.setCurrentEntry(i);
    //get number of samples
    int  numSamples= wrapper.getNumSamples();

    // loop over the number of waveforms at this ScanPoint (index j) 
    for ( int j=0; j<numSamples; ++j) {
      //if ( j>50 ) continue;
      double* pmtsample=wrapper.getPmtSample( pmt.pmt, j );
      hdiff->Fill( pmtsample[1] - pmtsample[ 0 ] );
      // only use first 20 time-bins of the waveform (time-bin index k)
      for ( int k = 0; k < std::min( 20, wrapper.getSampleLength() ); ++k ){
        habs->Fill( pmtsample[k] );
      }
    }
  }
}

std::ostream& operator << ( std::ostream& os, const MeanRMSCalc & calc ){
  os<<"N="<<calc.ndatapts()<<" Mean="<<calc.mean()<<" RMS="<<calc.rms()<<endl;
  return os;
}
