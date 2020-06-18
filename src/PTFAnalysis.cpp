#include "PTFAnalysis.hpp"
#include "Configuration.hpp"
#include "Utilities.hpp"
#include "TVirtualFFT.h"

#include "TH2D.h"

#include <iostream>
#include <ostream>
#include <fstream>

void PTFAnalysis::ChargeSum( float ped ){
  fitresult->qped = ped;
  float sum = 0.;
  for( int ibin = 1; ibin<=hwaveform->GetNbinsX(); ibin++ ){
    sum += ped - hwaveform->GetBinContent( ibin );
  }
  fitresult->qsum = sum;
}

bool PTFAnalysis::MonitorCut( float cut ){
  float ped = 0.0;
  int nbins = 20;
  for( int ibin = 1; ibin<=nbins; ibin++ ){
      ped += hwaveform->GetBinContent( ibin )/(float)nbins;
  }
  fitresult->ped = ped;
  float amp = 0.0;
  float mean = 0.0;
  for( int ibin = 1; ibin<=hwaveform->GetNbinsX(); ibin++ ){
      if( ped - hwaveform->GetBinContent( ibin ) > amp ){
          amp = ped - hwaveform->GetBinContent( ibin );
          mean = (float)ibin - 1.;
      }
  }
  fitresult->amp = amp;
  fitresult->mean = mean;
  if( amp < cut ){
    return false;
  }
  else{
    return true;
  }
}

bool PTFAnalysis::FFTCut(){
  // Compute the transform
  TVirtualFFT::SetTransform(0);
  //cout << "bin contents: ";
  //for( int i=1; i<=hwaveform->GetNbinsX(); i++ ){
  //  cout << hwaveform->GetBinContent(i) << " ";
  //}
  //cout << endl;
  //hwaveform->Print();
  hfftm = hwaveform->FFT( hfftm, "MAG" ); // Magnitude
  TVirtualFFT *fft = TVirtualFFT::GetCurrentTransform();
  delete fft;
  hfftm->SetBinContent(1, 0.0); // Remove pedestal
  // Cut if max bin not not near 0 Hz or below threshold
  int nbins = hwaveform->GetNbinsX();
  int maxBin = hfftm->GetMaximumBin();
  double maxValue = hfftm->GetBinContent(maxBin);
  fitresult->fftmaxbin = maxBin;
  fitresult->fftmaxval = maxValue;
  if( maxValue > 80.0 && ((maxBin > 1 && maxBin <= 4) || (maxBin >= nbins-3 && maxBin <= nbins)) ){
    return true;
  }
  else{
    return false;
  }
}

bool PTFAnalysis::PulseLocationCut( int cut ){
  // Cut if min bin in first or last bins
  int nbins = hwaveform->GetNbinsX();
  int minBin = hwaveform->GetMinimumBin();
  if( (minBin >= 1 && minBin <= cut) || (minBin >= nbins-cut+1 && minBin <= nbins) ){
    return false;
  }
  else{
    return true;
  }
}

//bool PTFAnalysis::HasWaveform(int pmt){
//  if( pmt == 0 ){
//    // check if fit is valid
//    if ( fitresult->fitstat != 0 ) return false;
//    // check if ringing is bigger than waveform
//    if ( fitresult->amp - fitresult->sinamp < 20.0 ) return false; 
//    // check if pulse width is reasonable
//    if ( fitresult->sigma < 1.0 || fitresult->sigma > 7.0 ) return false;
//    // check if mean pulse time is reasonable
//    if ( fitresult->mean < 28.0 || fitresult->mean > 45.0 ) return false;
//    // cut on chi2
//    //if ( fitresult->chi2 > 200.0 ) return false;
//  }
//  if( pmt == 1 ){
//    // Checks if monitor (PMT1) is being fitted
//    // check if fit is valid
//    //if ( fitresult->fitstat != 0 ) return false;
//    // check if ringing is bigger than waveform
//    //if ( fitresult->amp - fitresult->sinamp < 40.0 ) return false; 
//    // check if pulse width is reasonable
//    //if ( fitresult->sigma < 0.25 || fitresult->sigma > 5.0 ) return false;
//    // check if mean pulse time is reasonable
//    //if ( fitresult->mean < 25.0 || fitresult->mean > 55.0 ) return false;
//    // cut on chi2
//    //if ( fitresult->chi2 > 200.0 ) return false;
//
//    // Checks if using simpler analysis of bin furthest from pedestal
//    if ( fitresult->amp < 25.0 ) return false;
//  }
//  return true;
//}
  
// Returns constant reference to a particular fit result
// it is useable until next entry in TTree is added, or another get_fitresult
// is called again
const WaveformFitResult & PTFAnalysis::get_fitresult( int scanpt, unsigned long long wavenum ){
  ScanPoint scanpoint = scanpoints[ scanpt ];
  long long entry = scanpoint.get_entry();
  entry += wavenum;
  if ( entry > ptf_tree->GetEntries() ) throw std::runtime_error("WaveformFitResult::get_fitresult beyond end of TTree");
  ptf_tree->GetEvent( entry );
  return *fitresult;
}

//Gaussian fitting function
// x[0] is x value
// par[0] = amplitude
// par[1] = mean
// par[2] = sigma
// par[3] = offset
// par[4] = sin amplitude
// par[5] = sin frequency (rad/s)
// par[6] = sin phase (rad)
double PTFAnalysis::pmt0_gaussian(double *x, double *par) {
  double arg=0;
  if(par[2]!=0) arg=(x[0]-par[1])/par[2];
  double gfunc=par[3] - par[0] * TMath::Exp( -0.5*arg*arg ) + par[4]*TMath::Sin( par[5]*x[0] + par[6] );
  return gfunc;
}

//Gaussian fitting function
// x[0] is x value
// par[0] = amplitude
// par[1] = mean
// par[2] = sigma
// par[3] = offset
double PTFAnalysis::pmt1_gaussian(double *x, double *par) {
  double arg=0;
  if(par[2]!=0) arg=(x[0]-par[1])/par[2];
  double gfunc=par[3] - par[0] * TMath::Exp( -0.5*arg*arg );
  return gfunc;
}

void PTFAnalysis::InitializeFitResult( int wavenum, int nwaves  ) {
  fitresult->Init();
  ScanPoint & scanpoint = scanpoints[ scanpoints.size()-1 ];
  fitresult->scanpt    = scanpoints.size()-1;
  fitresult->wavenum   = wavenum;
  fitresult->nwaves    = nwaves;
  fitresult->x         = scanpoint.x();
  fitresult->y         = scanpoint.y();
  fitresult->z         = scanpoint.z();
}

void PTFAnalysis::FitWaveform( int wavenum, int nwaves, int pmt ) {
  // assumes hwaveform already defined and filled
  // assumes fit result structure already setup
  // Fit waveform for main PMT
  if( pmt == 0 ){
    // check if we need to build the function to fit
    if( fmygauss == nullptr ) fmygauss = new TF1("mygauss",pmt0_gaussian,0,70,7);
    fmygauss->SetParameters( 1.0, 35.0, 3.6, 8135.0, 10.0, 0.5, 0.0 );
    fmygauss->SetParNames( "Amplitude", "Mean", "Sigma", "Offset",
      		 "Sine-Amp",  "Sin-Freq", "Sin-Phase" );

    fmygauss->SetParLimits(0, 0.0, 8500.0);
    fmygauss->SetParLimits(1, 1.0, 69.0 );
    fmygauss->SetParLimits(2, 0.25, 10.0 );
    fmygauss->SetParLimits(3, 7500.0, 8500.0 );
    fmygauss->SetParLimits(4, 0.0, 8500.0);
    fmygauss->SetParLimits(5, 0.4, 0.7);
    fmygauss->SetParLimits(6, -TMath::Pi(), TMath::Pi() );
 
    // first fit for sine wave:
    fmygauss->FixParameter(0,1.0);
    fmygauss->FixParameter(1,35.0);
    fmygauss->FixParameter(2,3.6);
    hwaveform->Fit( fmygauss, "Q", "", 0,30.0);

    // then fit gaussian
    fmygauss->ReleaseParameter(0);
    fmygauss->ReleaseParameter(1);
    fmygauss->ReleaseParameter(2);
    fmygauss->SetParLimits(0, 0.0, 8500.0);
    fmygauss->SetParLimits(1, 1.0, 69.0 );
    fmygauss->SetParLimits(2, 0.25, 10.0 );
    fmygauss->FixParameter(3, fmygauss->GetParameter(3) );
    fmygauss->FixParameter(4, fmygauss->GetParameter(4));
    fmygauss->FixParameter(5, fmygauss->GetParameter(5));
    fmygauss->FixParameter(6, fmygauss->GetParameter(6));
    hwaveform->Fit( fmygauss, "Q", "", 20.0, 50.0);

    // then fit sine and gaussian together
    fmygauss->ReleaseParameter(3);
    fmygauss->ReleaseParameter(4);
    fmygauss->ReleaseParameter(5);
    fmygauss->ReleaseParameter(6);
    fmygauss->SetParLimits(0, 0.0, 8500.0);
    fmygauss->SetParLimits(1, 1.0, 69.0 );
    fmygauss->SetParLimits(2, 0.25, 10.0 );
    fmygauss->SetParLimits(3, 7500.0, 8500.0 );
    fmygauss->SetParLimits(4, 0.0, 8500.0);
    fmygauss->SetParLimits(5, 0.4, 0.7);
    fmygauss->SetParLimits(6, -TMath::Pi(), TMath::Pi() );
    int fitstat = hwaveform->Fit( fmygauss, "Q", "", 0, 70);
    // collect fit results
    fitresult->ped       = fmygauss->GetParameter(3);
    fitresult->mean      = fmygauss->GetParameter(1);
    fitresult->sigma     = fmygauss->GetParameter(2);
    fitresult->amp       = fmygauss->GetParameter(0);
    fitresult->sinamp    = fmygauss->GetParameter(4);
    fitresult->sinw      = fmygauss->GetParameter(5);
    fitresult->sinphi    = fmygauss->GetParameter(6);
    fitresult->ped_err   = fmygauss->GetParError(3);
    fitresult->mean_err  = fmygauss->GetParError(1);
    fitresult->sigma_err = fmygauss->GetParError(2);
    fitresult->amp_err   = fmygauss->GetParError(0);
    fitresult->sinamp_err= fmygauss->GetParError(4);
    fitresult->sinw_err  = fmygauss->GetParError(5);
    fitresult->sinphi_err= fmygauss->GetParError(6);
    fitresult->chi2      = fmygauss->GetChisquare();
    fitresult->ndof      = 30-4;
    fitresult->prob      = TMath::Prob( fmygauss->GetChisquare(), 30-4 );
    fitresult->fitstat   = fitstat;
  }
  // Simpler analysis for monitor PMT
  // Fit with simple gaussian
  // OR find bin furthest from pedestal
  //else if( pmt == 1 ){
  //  if( fmygauss == nullptr ) fmygauss = new TF1("mygauss",pmt1_gaussian,0,70,4);
  //  fmygauss->SetParameters( fitresult->amp, fitresult->mean, 1.0, fitresult->ped );
  //  fmygauss->SetParNames( "Amplitude", "Mean", "Sigma", "Offset" );

  //  fmygauss->SetParLimits(0, 0.0, 8500.0);
  //  fmygauss->SetParLimits(1, 0.0, 70.0 );
  //  fmygauss->SetParLimits(2, 0.01, 3.0 );
  //  fmygauss->SetParLimits(3, 7500.0, 9000.0 );

  //  // then fit gaussian
  //  int fitstat = hwaveform->Fit( fmygauss, "Q", "", 30.0, 50.0);

  //  // collect fit results
  //  fitresult->ped       = fmygauss->GetParameter(3);
  //  fitresult->mean      = fmygauss->GetParameter(1);
  //  fitresult->sigma     = fmygauss->GetParameter(2);
  //  fitresult->amp       = fmygauss->GetParameter(0);
  //  fitresult->chi2      = fmygauss->GetChisquare();
  //  fitresult->ndof      = 30-4;
  //  fitresult->prob      = TMath::Prob( fmygauss->GetChisquare(), 30-4 );
  //  fitresult->fitstat   = fitstat;
  //}
  else if( pmt == 1 ){
    float ped = 0.0;
    int nbins = 20;
    for( int ibin = 1; ibin<=nbins; ibin++ ){
      ped += hwaveform->GetBinContent( ibin )/(float)nbins;
    }
    fitresult->ped = ped;
    float amp = 0.0;
    float mean = 0.0;
    for( int ibin = 1; ibin<=hwaveform->GetNbinsX(); ibin++ ){
      if( ped - hwaveform->GetBinContent( ibin ) > amp ){
        amp = ped - hwaveform->GetBinContent( ibin );
        mean = (float)ibin - 1.;
      }
    }
    fitresult->amp = amp;
    fitresult->mean = mean;
  }
  else if( pmt == 2 ){ /// Add a new PMT type for the mPMT analysis.
    if( fmygauss == nullptr ) fmygauss = new TF1("mygauss",pmt1_gaussian,255,290,4);
    fmygauss->SetParameters( fitresult->amp, fitresult->mean, 1.0, fitresult->ped );
    fmygauss->SetParNames( "Amplitude", "Mean", "Sigma", "Offset" );

    fmygauss->SetParameter(0, 100.0);
    fmygauss->SetParameter(1, 270.0 );
    fmygauss->SetParameter(2, 1.4 );
    fmygauss->SetParameter(3, 2040.0 );
    fmygauss->SetParLimits(0, 0.0, 500.0);
    fmygauss->SetParLimits(1, 265.0, 275.0 );
    fmygauss->SetParLimits(2, 0.8, 2.0 );
    fmygauss->SetParLimits(3, 2035.0, 2045.0 );

    // then fit gaussian
    int fitstat = hwaveform->Fit( fmygauss, "Q", "", 255, 290.0);

    // collect fit results
    fitresult->ped       = fmygauss->GetParameter(3);
    fitresult->mean      = fmygauss->GetParameter(1);
    fitresult->sigma     = fmygauss->GetParameter(2);
    fitresult->amp       = fmygauss->GetParameter(0);
    fitresult->chi2      = fmygauss->GetChisquare();
    fitresult->ndof      = 30-4;
    fitresult->prob      = TMath::Prob( fmygauss->GetChisquare(), 30-4 );
    fitresult->fitstat   = fitstat;
  }

  else{
    cout << "PTFAnalysis::FitWaveForm Error: PMT not 0 or 1!" << endl;
    exit( EXIT_FAILURE );
  }
}

PTFAnalysis::PTFAnalysis( TFile* outfile, PTF::Wrapper & wrapper, double errorbar, PTF::PMTChannel & channel, string config_file, bool savewf ){

  // Load config file
  Configuration config;
  bool terminal_output;
  bool pulse_location_cut;
  bool fft_cut;

  config.Load(config_file);
  if( !config.Get("terminal_output", terminal_output) ){
    cout << "Missing terminal_output parameter from config file." << endl;
    exit( EXIT_FAILURE );
  }
  if( !config.Get("pulse_location_cut", pulse_location_cut) ){
    cout << "Missing pulse_location_cut parameter from config file." << endl;
    exit( EXIT_FAILURE );
  }
  if( !config.Get("fft_cut", fft_cut) ){
    cout << "Missing fft_cut parameter from config file." << endl;
    exit( EXIT_FAILURE );
  }

  static int instance_count =0;
  static int savewf_count =0;
  static int savenowf_count = 0;
  ++instance_count;
  save_waveforms = savewf;
  
  // Get utilities
  Utilities utils;

  // get length of waveforms
  wrapper.setCurrentEntry(0);
  int  numTimeBins= wrapper.getSampleLength();
  
  // build the waveform histogram
  std::string hname = "hwaveform" + std::to_string(instance_count);
  std::string hname_fft = "hfftm" + std::to_string(instance_count);
  outfile->cd();
  hwaveform = new TH1D( hname.c_str(), "Pulse waveform; Time bin (tdc clock ticks); Charge (adc units)", numTimeBins, 0., float(numTimeBins) );
  hfftm = new TH1D( hname_fft.c_str(), "Fast Fourier Transform; Frequency; Coefficient", numTimeBins, -5.0e8, 5.0e8 );
  
  // set up the output TTree
  string ptf_tree_name = "ptfanalysis" + std::to_string(channel.pmt);
  ptf_tree = new TTree(ptf_tree_name.c_str(), ptf_tree_name.c_str());
  fitresult = new WaveformFitResult();
  fitresult->MakeTTreeBranches( ptf_tree );

  // Create output directories
  // Directories for waveforms
  string wfdir_name = "PMT" + std::to_string(channel.pmt) + "_Waveforms";
  string nowfdir_name = "PMT" + std::to_string(channel.pmt) + "_NoWaveforms";
  if ( save_waveforms && wfdir==nullptr ) wfdir = outfile->mkdir(wfdir_name.c_str());
  if ( save_waveforms && nowfdir==nullptr ) nowfdir = outfile->mkdir(nowfdir_name.c_str());
  // Directories for FFTs
  string wfdir_fft_name = "FFT" + std::to_string(channel.pmt) + "_Waveforms";
  string nowfdir_fft_name = "FFT" + std::to_string(channel.pmt) + "_NoWaveforms";
  if ( save_waveforms && wfdir_fft==nullptr ) wfdir_fft = outfile->mkdir(wfdir_fft_name.c_str());
  if ( save_waveforms && nowfdir_fft==nullptr ) nowfdir_fft = outfile->mkdir(nowfdir_fft_name.c_str());
  outfile->cd();
  
  // Loop over scan points (index i)
  unsigned long long nfilled = 0;// number of TTree entries so far
  for (unsigned i = 2; i < wrapper.getNumEntries(); i++) {
    //if ( i>2000 ) continue;
    if( terminal_output ){
      cerr << "PTFAnalysis scan point " << i << " / " << wrapper.getNumEntries() << "\u001b[34;1m (" << (((double)i)/wrapper.getNumEntries()*100) << "%)\u001b[0m\033[K";
      cerr << "\r";
    }
    else{
      if ( i % 10 == 0 ){
        std::cout << "PTFAnalysis scan point " << i << " / " << wrapper.getNumEntries() << std::endl;
      }
    }
    wrapper.setCurrentEntry(i);
    auto location = wrapper.getDataForCurrentEntry(PTF::Gantry1);
    
    scanpoints.push_back( ScanPoint( location.x, location.y, location.z, nfilled ) );
    ScanPoint& curscanpoint = scanpoints[ scanpoints.size()-1 ];
    
    // loop over the number of waveforms at this ScanPoint (index j)
    int numWaveforms = wrapper.getNumSamples();
    for ( int j=0; j<numWaveforms; j++) {
      //if( j>20 ) continue;
      double* pmtsample=wrapper.getPmtSample( channel.pmt, j );
      // set the contents of the histogram
      hwaveform->Reset();
      for ( int ibin=1; ibin <= numTimeBins; ++ibin ){
        hwaveform->SetBinContent( ibin, pmtsample[ibin-1] );
	    hwaveform->SetBinError( ibin, errorbar );
      }
      InitializeFitResult( j, numWaveforms );
      // Do simple charge sum calculation
      if( channel.pmt == 0 ) ChargeSum(8135.4);
      // For main PMT do FFT and check if there is a waveform
      // If a waveform present then fit it
      bool dofit = true;
      if( dofit && pulse_location_cut && channel.pmt == 0 ) dofit = PulseLocationCut(10);
      if( dofit && fft_cut && channel.pmt == 0 ) dofit = FFTCut();
      //if( dofit && channel.pmt == 1 ) dofit = MonitorCut( 25. );
      if( dofit ){
        FitWaveform( j, numWaveforms, channel.pmt ); // Fit waveform and copy fit results into TTree
      }
      fitresult->haswf = utils.HasWaveform( fitresult, channel.pmt );
      ptf_tree->Fill();
      
      // check if we should clone waveform histograms
      if ( save_waveforms && savewf_count<500 && savenowf_count<500 ){
	    if  ( fabs( curscanpoint.x() - 0.46 ) < 0.0005 && 
	      fabs( curscanpoint.y() - 0.38 ) < 0.0005 ) {
              
          std::string hwfname = "hwf_" + std::to_string( nfilled );
          std::string hfftmname = "hfftm_" + std::to_string( nfilled );
          if ( fitresult->haswf && savewf_count<500 ) {
            wfdir->cd();
            TH1D* hwf = (TH1D*) hwaveform->Clone( hwfname.c_str() );
            hwf->SetName( hwfname.c_str() );
            hwf->SetTitle("HAS a pulse; Time bin (tdc clock ticks); Charge (adc units)");
            hwf->SetDirectory( wfdir );
            wfdir_fft->cd();
            TH1D* hfftm_tmp = (TH1D*) hfftm->Clone( hfftmname.c_str() );
            hfftm_tmp->SetName( hfftmname.c_str() );
            hfftm_tmp->SetTitle("HAS a pulse; Frequency; Coefficient");
            hfftm_tmp->SetDirectory( wfdir_fft );
            ++savewf_count;	  
          } else if ( !fitresult->haswf && savenowf_count<500 ){
            nowfdir->cd();
            TH1D* hwf = (TH1D*) hwaveform->Clone( hwfname.c_str() );
            hwf->SetName( hwfname.c_str() );
            hwf->SetTitle("Noise pulse; Time bin (tdc clock ticks); Charge (adc units)");
            hwf->SetDirectory( nowfdir );
            nowfdir_fft->cd();
            TH1D* hfftm_tmp = (TH1D*) hfftm->Clone( hfftmname.c_str() );
            hfftm_tmp->SetName( hfftmname.c_str() );
            hfftm_tmp->SetTitle("Noise pulse; Frequency; Coefficient");
            hfftm_tmp->SetDirectory( nowfdir_fft );
            ++savenowf_count;	  
          }

	    outfile->cd();


        }

      }
      ++curscanpoint;  // increment counters
      ++nfilled;
    }
  }
  //cout << endl;
  // Done.

}

const std::vector< double > PTFAnalysis::get_bins( char dim ){

  vector< double > positions;
  for(unsigned int iscan=0; iscan<get_nscanpoints(); iscan++){
    ScanPoint scanpoint = scanpoints[ iscan ];
    if( scanpoint.x() < 1e-5 ) continue; // Ignore position (0,0,0)
    if( dim == 'x' ){
      positions.push_back( scanpoint.x() );
    }
    else if( dim == 'y' ){
      positions.push_back( scanpoint.y() ); 
    }
    else if( dim == 'z' ){
      positions.push_back( scanpoint.z() );
    }
    else{
      cout << "PTFAnalysis::get_bins Error: input must be x, y or z!" << endl;
      exit( EXIT_FAILURE );
    }
  }

  sort( positions.begin(), positions.end() );
  positions.erase( unique( positions.begin(), positions.end(), comparison ), positions.end() ); // comparison function may not be necessary but was concerned about the comparison of doubles

  //std::cout << "positions contains:";
  //vector<double>::iterator it;
  //for (it=positions.begin(); it!=positions.end(); ++it)
  //  std::cout << ' ' << *it;
  //std::cout << '\n';

  if( positions.size() < 3 ){
    cout << "PTFAnalysis::get_bins Error: fewer than 3 positions a problem for automatic binning!" << endl;
    exit( EXIT_FAILURE );
  }

  vector< double > bins;
  for(unsigned int i=0; i<positions.size()-1; i++){
    bins.push_back( (positions[i]+positions[i+1])/2. );
  }
  bins.push_back( 2.*bins[ bins.size()-1 ] - bins[ bins.size()-2 ] );
  bins.insert( bins.begin(), 2.*bins[0] - bins[1] );

  //std::cout << "bins contains:";
  //vector<double>::iterator it;
  //for (it=bins.begin(); it!=bins.end(); ++it)
  //  std::cout << ' ' << *it;
  //std::cout << '\n';

  return bins;
}
