#ifndef __ERRORBARANALYSIS__
#define __ERRORBARANALYSIS__

#include "wrapper.hpp"
#include "MeanRMSCalc.hpp"

#include "TFile.h"
#include "TH1D.h"
#include <ostream>
#include <sstream>

/// Class to calculate the error bars on the PTF waveforms
/// Uses time bins 1-20 of the waveforms
/// Produces two histograms:
///  1) that is the pedestal using all of the data in bins 1-20
///  2) that uses difference between adjacent time-bins 1 and 2 collected charge

//const int ntimebin = 20;
std::ostream& operator << ( std::ostream& os, const MeanRMSCalc & calc );

class ErrorBarAnalysis {

public:
  ErrorBarAnalysis( TFile * outfile, PTF::Wrapper & ptf, PTF::PMTChannel & channel );
  double get_errorbar( ) const { return hdiff->GetRMS(); }

private:
  TH1D * habs; //absolute pedestal
  TH1D * hdiff; //differential pedestals

};

#endif // __ERRORBARANALYSIS__
