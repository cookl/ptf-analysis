#include "PulseFinding.hpp"
#include <iostream>
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
  if(pmt.channel == 0){ baseline = 0.9985; }
  if(pmt.channel == 1){ baseline = 1.006; }
  
  double threshold = baseline - 0.004;

  int nsamples = hwaveform->GetNbinsX();

  double min_bin = 9999, min_value = 9999; // Find the minimum bin and value for each pulse.
  bool in_pulse = false; // are we currently in a pulse?
  for(int ib = 1; ib < nsamples+1; ib++){
    double sample = hwaveform->GetBinContent(ib);

    if(sample < min_value){
      min_value = sample;
      min_bin = ib;
    }
    
    if(sample < threshold && !in_pulse){ // found a pulse
      in_pulse = true;      
    }
    
    if(sample >= threshold && in_pulse){ /// finished this pulse
      in_pulse = false;
      if(fitresult->numPulses < MAX_PULSES){
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

void another_pulse_finding_function(TH1D *hwaveform, WaveformFitResult *fitresult, PTF::PMT pmt) {
  
  // initializing the baseline
  double baseline = 1.0;
  if(pmt.channel == 0){ baseline = 0.9985; }
  if(pmt.channel == 1){ baseline = 1.006; }

  // a threshold to determine if peak is above it
  double threshold = baseline - 0.004;

  // number of samples in a waveform
  int nsamples = hwaveform->GetNbinsX();

  // initializing the extremum counter
  int numExtrema = 0;
  int numMinima  = 0;
  int numMaxima  = 0;
  int numDontKnow= 0;

  // extrema contains all minima and maxima
  std::vector<int> extrema;
  std::vector<int> maxima;
  std::vector<int> minima;
  std::vector<int> dontknow;

  for(int ib = 1; ib < nsamples; ib++){

    double sample = hwaveform->GetBinContent(ib);
    double sample_plus = hwaveform->GetBinContent(ib+1);
    double sample_minus = hwaveform->GetBinContent(ib-1);

    if ( abs(sample_minus) <= abs(sample) && abs(sample_plus) <= abs(sample) && abs(sample) >= threshold) {
      // local extremum found
      extrema.push_back(ib);
      numExtrema++;
      if (sample_minus - sample > 0.) {
        // local minimum
        minima.push_back(ib);
        numMaxima++;
      } else if (sample_minus - sample < 0.) {
        // local maximum
        maxima.push_back(ib);
        numMinima++;
      } else {
        // not an extremum. what happened?
        dontknow.push_back(ib);
        numDontKnow++;
      }
    } 
  }
}


