#include "WaveformFitResult.hpp"

void WaveformFitResult::Init(){
  scanpt=-1; wavenum=0; nwaves=0;
  x=0.; y=0.; z=0.;
  ped=0.; ped_err=0.;
  mean=0.; mean_err=0.;
  sigma=0.; sigma_err=0.;
  amp=0.; amp_err=0.;
  sinw=0.; sinw_err=0.;
  sinphi=0.; sinphi_err=0.;
  chi2=0.; ndof=0.; prob=0.;
  fitstat=-1;
  fftmaxbin=-1; fftmaxval=0.;
  haswf=0; qped=0.; qsum=-999.;
  numPulses = 0;
  for(int i = 0; i < MAX_PULSES; i++){
    pulseTimes[i] = 0.0;
    pulseTimeErr[i] = 0.0;
    pulseCharges[i] = 0.0;
    pulseChargeErr[i] = 0.0;
  }
}

const char * WaveformFitResult::GetRootString(){
  if (rootstring != "" ) return rootstring.c_str();
  rootstring ="scanpt/I";
  rootstring+=":wavenum/I";
  rootstring+=":nwaves/I";
  rootstring+=":x/F";
  rootstring+=":y/I";
  rootstring+=":z/I";
  rootstring+=":ped/I";
  rootstring+=":ped_err/I";
  rootstring+=":mean/F";
  rootstring+=":mean_err/F";
  rootstring+=":sigma/F";
  rootstring+=":sigma_err/F";
  rootstring+=":amp/F";
  rootstring+=":amp_err/F";
  rootstring+=":sinamp/F";     //< ringing amplitude
  rootstring+=":sinamp_err/F";
  rootstring+=":sinw/F";       //< ringing frequency (rad/s)
  rootstring+=":sinw_err/F"; 
  rootstring+=":sinphi/F";     //< ringing phase offset (rad)
  rootstring+=":sinphi_err/F";
  rootstring+=":chi2/F";
  rootstring+=":ndof/F";
  rootstring+=":prob/F";
  rootstring+=":fitstat/I";
  rootstring+=":fftmaxbin/I";
  rootstring+=":fftmaxval/F";
  rootstring+=":haswf/I";
  rootstring+=":qped/F";
  rootstring+=":qsum/F";
  rootstring+=":numPulses/I";
  rootstring+=":pulseTimes[numPulses]/F";
  rootstring+=":pulseTimeErr[numPulses]/F";
  rootstring+=":pulseCharge[numPulses]/F";
  rootstring+=":pulseChargeErr[numPulses]/F";
  return rootstring.c_str(); 
}

void WaveformFitResult::MakeTTreeBranches(TTree * t){
  t->Branch( "scanpt",    &scanpt,    "scanpt/I" );
  t->Branch( "wavenum",   &wavenum,   "wavenum/I" );
  t->Branch( "nwaves",    &nwaves,    "nwaves/I" );
  t->Branch( "x",         &x,         "x/F" );
  t->Branch( "y",         &y,         "y/F" );
  t->Branch( "z",         &z,         "z/F" );
  t->Branch( "ped",       &ped,       "ped/F" );
  t->Branch( "ped_err",   &ped_err,   "ped_err/F" );
  t->Branch( "mean",      &mean,      "mean/F" );
  t->Branch( "mean_err",  &mean_err,  "mean_err/F" );
  t->Branch( "sigma",     &sigma,     "sigma/F" );
  t->Branch( "sigma_err", &sigma_err, "sigma_err/F" );
  t->Branch( "amp",       &amp,       "amp/F" );
  t->Branch( "amp_err",   &amp_err,   "amp_err/F" );
  t->Branch( "sinamp",    &sinamp,    "sinamp/F" );     //< ringing amplitude
  t->Branch( "sinamp_err",&sinamp_err,"sinamp_err/F");
  t->Branch( "sinw",      &sinw,      "sinw/F");       //< ringing frequency (rad/s)
  t->Branch( "sinw_err",  &sinw_err,  "siw_err/F"); 
  t->Branch( "sinphi",    &sinphi,    "sinphi/F");     //< ringing phase offset (rad)
  t->Branch( "sinphi_err",&sinphi_err,"sinphi_err/F");
  t->Branch( "chi2",      &chi2,      "chi2/F" );
  t->Branch( "ndof",      &ndof,      "ndof/F" );
  t->Branch( "prob",      &prob,      "prob/F" );
  t->Branch( "fitstat",   &fitstat,   "fitstat/I" );
  t->Branch( "fftmaxbin", &fftmaxbin, "fftmaxbin/I" );
  t->Branch( "fftmaxval", &fftmaxval, "fftmaxval/F" );
  t->Branch( "haswf",     &haswf,     "haswf/I" );
  t->Branch( "qped",      &qped,      "qped/F" );
  t->Branch( "qsum",      &qsum,      "qsum/F" );
  t->Branch( "numPulses", &numPulses, "numPulses/I" );
  t->Branch( "pulseTimes",pulseTimes,"pulseTimes[numPulses]/F" );
  t->Branch( "pulseTimeErr",pulseTimeErr,"pulseTimeErr[numPulses]/F" );
  t->Branch( "pulseCharges",pulseCharges,"pulseCharges[numPulses]/F" );
  t->Branch( "pulseChargeErr",pulseChargeErr,"pulseChargeErr[numPulses]/F" );
  
  return;
}

void WaveformFitResult::SetBranchAddresses(TTree * t){
  t->SetBranchAddress( "scanpt",    &scanpt );
  t->SetBranchAddress( "wavenum",   &wavenum );
  t->SetBranchAddress( "nwaves",    &nwaves );
  t->SetBranchAddress( "x",         &x );
  t->SetBranchAddress( "y",         &y );
  t->SetBranchAddress( "z",         &z );
  t->SetBranchAddress( "ped",       &ped );
  t->SetBranchAddress( "ped_err",   &ped_err );
  t->SetBranchAddress( "mean",      &mean );
  t->SetBranchAddress( "mean_err",  &mean_err );
  t->SetBranchAddress( "sigma",     &sigma );
  t->SetBranchAddress( "sigma_err", &sigma_err );
  t->SetBranchAddress( "amp",       &amp );
  t->SetBranchAddress( "amp_err",   &amp_err );
  t->SetBranchAddress( "sinamp",    &sinamp );     //< ringing amplitude
  t->SetBranchAddress( "sinamp_err",&sinamp_err );
  t->SetBranchAddress( "sinw",      &sinw );       //< ringing frequency (rad/s)
  t->SetBranchAddress( "sinw_err",  &sinw_err ); 
  t->SetBranchAddress( "sinphi",    &sinphi );     //< ringing phase offset (rad)
  t->SetBranchAddress( "sinphi_err",&sinphi_err );
  t->SetBranchAddress( "chi2",      &chi2 );
  t->SetBranchAddress( "ndof",      &ndof );
  t->SetBranchAddress( "prob",      &prob );
  t->SetBranchAddress( "fitstat",   &fitstat );
  t->SetBranchAddress( "fftmaxbin", &fftmaxbin );
  t->SetBranchAddress( "fftmaxval", &fftmaxval );
  t->SetBranchAddress( "haswf",     &haswf );
  t->SetBranchAddress( "qped",      &qped );
  t->SetBranchAddress( "qsum",      &qsum );
  t->SetBranchAddress( "numPulses", &numPulses );
  t->SetBranchAddress( "pulseTimes", pulseTimes );
  t->SetBranchAddress( "pulseTimeErr", pulseTimeErr );
  t->SetBranchAddress( "pulseCharges", pulseCharges );
  t->SetBranchAddress( "pulseChargeErr", pulseChargeErr );
  return;
}
