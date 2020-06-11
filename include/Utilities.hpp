#ifndef __UTILITIES__
#define __UTILITIES__

#include "ScanPoint.hpp"
#include "WaveformFitResult.hpp"

#include <vector>
#include <cmath>

/// Class for providing utility functions for analyses

class Utilities {

public:
  //default constructor takes no input
  Utilities(){};
  //get vector of bin edges from scan points
  const std::vector<double> get_bins( std::vector< ScanPoint > scanpoints, char dim);
  //set a style for plots
  void set_style();
  //does the PMT signal contain a waveform?
  bool HasWaveform( WaveformFitResult *wf, int pmt );

private:
  static bool comparison (double i, double j){ return (fabs( i-j ) < 1e-5); }

};

#endif // __UTILITIES__
