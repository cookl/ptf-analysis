#ifndef __PTFANALYSIS__
#define __PTFANALYSIS__

#include "TFile.h"
#include "TF1.h"
#include "TH1D.h"
#include "TTree.h"
#include "TDirectory.h"
#include "TMath.h"
#include <vector>
#include <string>


#include "wrapper.hpp"
#include "ScanPoint.hpp"
#include "WaveformFitResult.hpp"

using namespace std;

/// This class takes a PTFWrapper reference, charge errorbar, and PMT as input
/// It then does main analysis to fill a TTree of WaveformFitResults
/// Has methods to later read back entries of the TTree
/// Keeps track of number of Scan Points, and locations used find entries in TTree
class PTFAnalysis {
public:
  PTFAnalysis( TFile * outfile,Wrapper & ptf, double errorbar, PTF::PMT & pmt, string config_file, bool savewf=false );
  ~PTFAnalysis(){
    if ( fitresult ) delete fitresult;
  }

  // Access fit results
  const std::vector< ScanPoint > & get_scanpoints()             const { return scanpoints; };
  const unsigned                   get_nscanpoints()            const { return scanpoints.size(); }
  const unsigned long long         get_firstentry( int scanpt ) const { return scanpoints[ scanpt ].get_entry(); }
  const unsigned long long         get_nentries( int scanpt )   const { return scanpoints[ scanpt ].nentries(); }
  
  const WaveformFitResult &        get_fitresult( int scanpt, unsigned long long wavenum );

  // write scanpoints information into a separate TTree
  void                             write_scanpoints(){ WriteScanPoints( scanpoints ); }

  // Post-fitresult analysis
  const std::vector< double >      get_bins( char dim );
  
private:
  void ChargeSum( float ped, int bin_low=1, int bin_high=0 ); // Charge sum relative to ped
  bool MonitorCut( float cut ); // Cut if no monitor PMT pulse
  bool FFTCut(); // Do FFT and check if waveform present
  bool PulseLocationCut( int cut ); // Cut on pulse in first or last bins
  void InitializeFitResult( int wavenum, int nwaves  );

  void FitWaveform( int wavenum, int nwaves, PTF::PMT pmt );
  static double pmt0_gaussian(double *x, double *par);
  static double pmt1_gaussian(double *x, double *par);
  static double funcEMG(double* x, double* p);
  static double pmt2_piecewise(double *x, double *par);
  static double bessel(double *x, double *p);
  static bool comparison(double i, double j){ return (fabs( i-j ) < 1e-5); }

  std::vector< ScanPoint > scanpoints;

  //std::vector< ScanPoint > Temperature;
  //TF1* fmygauss{nullptr};  // gaussian function used to fit waveform
  TF1* ffitfunc{nullptr};  // function used to fit waveform

  TH1D* hwaveform{nullptr}; // current waveform
  TH1* hfftm{nullptr}; // fast fourier transform magnitude
  WaveformFitResult * fitresult{nullptr};
  TTree* ptf_tree{nullptr};
  bool  save_waveforms{false};
  TDirectory* wfdir{nullptr};
  TDirectory* nowfdir{nullptr};
  TDirectory* wfdir_fft{nullptr};
  TDirectory* nowfdir_fft{nullptr};

};

#endif // __PTFANALYSIS__
