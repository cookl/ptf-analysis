#ifndef __PULSEFINDING__
#define __PULSEFINDING__

#include "TH1D.h"
#include "WaveformFitResult.hpp"

// Find pulses in a given waveform
// arguments:
// int algo_type : which pulse finding algorithm to use?
// TH1D *hwaveform : the input waveform
// WaveformFitResult *fitresult : store the list of pulses in WaveformFitResult
void find_pulses(int algo_type, TH1D *hwaveform, WaveformFitResult *fitresult );

// Do simple comparison to fixed threshold to find pulses
void simple_threshold_technique(TH1D *hwaveform, WaveformFitResult *fitresult );

#endif // __PULSEFINDING__

