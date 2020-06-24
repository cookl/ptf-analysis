#ifndef __WAVEFORMFITRESULT__
#define __WAVEFORMFITRESULT__

#include "TTree.h"
#include <string>

#define MAX_PULSES 10

/// Structure to hold information about one waveform fit
/// It is of the form that it can be used to set up a TTree readout
class WaveformFitResult {
public:
  WaveformFitResult(){ Init(); } // default constructor
  ~WaveformFitResult(){}; // default destructor
  void Init();
  const char * GetRootString();
  void MakeTTreeBranches(TTree * t);
  void SetBranchAddresses(TTree * t);

  int scanpt;       //< scan point number
  int wavenum;      //< waveform number in scan
  int nwaves;       //< number of waveforms in this scan point
  float x;          //< x location of this scan point
  float y;          //< y location of this scan point
  float z;          //< z location of this scan point
  float ped;        //< pedestal value of waveform
  float ped_err;    //< fit uncertainty in pedestal value
  float mean;       //< mean time of waveform
  float mean_err;   //< fit uncertainty in mean time of waveform
  float sigma;      //< sigma in time of this waveform
  float sigma_err;  //< fit uncertainty in sigma of time of waveform
  float amp;        //< amplitude of waveform
  float amp_err;    //< fit uncertainty in fitted amplitude of waveform
  float sinamp;     //< ringing amplitude
  float sinamp_err;
  float sinw;       //< ringing frequency (rad/s)
  float sinw_err; 
  float sinphi;     //< ringing phase offset (rad)
  float sinphi_err;
  float chi2;       //< fit chi2
  float ndof;       //< fit ndof
  float prob;       //< fit p-value
  int   fitstat;    //< 0 is good, any other number is an error in fit
  int   fftmaxbin;  //< maximum bin of FFT
  float fftmaxval;  //< value of FFT maximum bin
  int   haswf;      //< 1 if has a pulse, 0 if not
  float qped;       //< fixed pedestal
  float qsum;       //< charge sum

  int numPulses; //< number of pulses found in waveform
  float pulseTimes[MAX_PULSES]; // Pulse times
  float pulseTimeErr[MAX_PULSES]; // Pulse time errors
  float pulseCharges[MAX_PULSES]; // Pulse charges
  float pulseChargeErr[MAX_PULSES]; // Pulse charge errors


private:
  // hold copy of string needed for making TTree
  std::string rootstring;

};

#endif // __WAVEFORMFITRESULT__
