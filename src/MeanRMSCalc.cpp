#include "MeanRMSCalc.hpp"

MeanRMSCalc::MeanRMSCalc( const std::vector<double>& invec ){
  for ( double val : invec ){
    add( val );
  }
}

void MeanRMSCalc::add( double val ){
  ++npts;
  sum += val;
  sum2 += val*val;
  fMean = sum / npts;
  if ( npts > 1 ) { 
    fRMS = sum2/(npts-1) - fMean*fMean*npts/(npts-1);
    if ( fRMS < 0 ) fRMS=0.0;
    else fRMS = sqrt( fRMS  );
  }
}
