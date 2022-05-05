#ifndef __PULSEFINDING__
#define __PULSEFINDING__

#include "TH1D.h"
#include "WaveformFitResult.hpp"
#include "wrapper.hpp"

// Find pulses in a given waveform
// arguments:
// int algo_type : which pulse finding algorithm to use?
// TH1D *hwaveform : the input waveform
// WaveformFitResult *fitresult : store the list of pulses in WaveformFitResult
void find_pulses(int algo_type, TH1D *hwaveform, WaveformFitResult *fitresult, PTF::PMT pmt );

// Do simple comparison to fixed threshold to find pulses
void simple_threshold_technique(TH1D *hwaveform, WaveformFitResult *fitresult, PTF::PMT pmt );

// Another function to find peaks and compare the results to the first one
void another_pulse_finding_function(TH1D *hwaveform, WaveformFitResult *fitresult, PTF::PMT pmt );

#endif // __PULSEFINDING__

