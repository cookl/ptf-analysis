#include "PulseFinding.hpp"
#include <iostream>

#include "BrbSettingsTree.hxx"

#include <vector>



void find_pulses(int algo_type, TH1D *hwaveform, WaveformFitResult *fitresult, PTF::PMT pmt){

  // Reset the number of pulses
  fitresult->numPulses = 0;

  if(algo_type == 0){
    simple_threshold_technique(hwaveform, fitresult, pmt);
  }else{
    std::cerr << "Invalid pulse finding algorithm = " << algo_type
              << ". Exiting"<< std::endl;
    exit(0);
  }
  
  
}



void simple_threshold_technique(TH1D *hwaveform, WaveformFitResult *fitresult, PTF::PMT pmt){
  
  // Loop over waveform; look for every case of waveform going below fixed threshold
  
  double baseline = 1.0;

  // If it is an mPMT channel, the use the baseline from BRB settings tree.
  if(pmt.type == PTF::mPMT_REV0_PMT){
    baseline = BrbSettingsTree::Get()->GetBaseline(pmt.channel);
  }

  double threshold = baseline - 0.004;
    
  //if (pmt.channel == 1) threshold = baseline - 0.2;

  int nsamples = hwaveform->GetNbinsX();

  double min_bin = 9999, min_value = 9999; // Find the minimum bin and value for each pulse.
  bool in_pulse = false; // are we currently in a pulse?
  int start_bin;
  int end_bin;
  for(int ib = 1; ib < nsamples+1; ib++){
    double sample = hwaveform->GetBinContent(ib);

    if(sample < min_value){
      min_value = sample;
      min_bin = ib;
    }
    
    if(sample < threshold && !in_pulse){ // found a pulse
      in_pulse = true;
      start_bin = ib;
    }
    
    if(sample >= threshold && in_pulse){ // finished this pulse
      in_pulse = false;
      end_bin = ib;
      double left_bin;
      double right_bin;
      double left_value;
      double right_value;
      if(fitresult->numPulses < MAX_PULSES){
	double mid_value = baseline - (baseline - min_value) / 2;
	for(int bin = start_bin; bin < end_bin; bin++){
	  double left_check = hwaveform->GetBinContent(bin);
	  double right_check = hwaveform->GetBinContent(bin+1);
	  if((left_check >= mid_value) && (right_check <= mid_value)){
	    left_bin = bin;
	    right_bin = bin+1;
	    left_value = left_check;
	    right_value = right_check;
	    break;
	  }
	}
	double slope = (right_value - left_value) / (right_bin - left_bin);
        double y_int = left_value - slope*left_bin;
        double mid_bin = (mid_value - y_int) / slope;
	fitresult->pulseTimesCFD[fitresult->numPulses] = mid_bin * 8.0;
        fitresult->pulseTimes[fitresult->numPulses] = min_bin * 8.0;
        fitresult->pulseCharges[fitresult->numPulses] = baseline - min_value;

          //if(0)std::cout << "Pulse found : " << fitresult->numPulses << " " << min_bin
          //          << " " << min_value << " " << std::endl;
          fitresult->numPulses++;
      }
      min_bin = 9999, min_value = 9999;
    }
    
  }

  //  if(fitresult->numPulses){ std::cout << "number of pulses found: " << fitresult->numPulses << " " << std::endl;}
  
}

