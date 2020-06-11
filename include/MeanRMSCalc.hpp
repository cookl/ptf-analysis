#ifndef __MEANRMSCALC__
#define __MEANRMSCALC__

#include <vector>
#include <cmath>

/// Class for calculating and keeping track of the mean and RMS of a collection of data

class MeanRMSCalc {

public:
  MeanRMSCalc(){};//default constructor takes no input
  MeanRMSCalc( const std::vector< double > & invec );
  void   add( double val );                  //< add a data point
  double mean() const { return fMean; }  //< get current mean of points added
  double rms()  const { return fRMS; };  //< get current RMS of points added
  unsigned long long ndatapts() const { return npts; }
  
private:
  double sum{0.};
  double sum2{0.};
  double fMean{0.};
  double fRMS{0.};
  unsigned long long npts{0};

};

#endif // __MEANRMSCALC__
