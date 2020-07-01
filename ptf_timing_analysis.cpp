/* ptf_charge_analysis.cpp
   Analysis of fitted waveforms to determine the PMT timing response.
   Produces histogram of the PMT pulse location.
   Currently taken relative to digitised trigger pulse.
    - Could change this to monitor PMT pulse, but narrowness results in larger error.
   Relative transit time (RTT) is the mean minus scan point with earliest time.
   Transit time spread (TTS) is the std dev.
   
   Author: John Walker (Jan 2020)
 */

#include "WaveformFitResult.hpp"
#include "ScanPoint.hpp"
#include "Utilities.hpp"
#include "FindCircle.hpp"
#include "TFile.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TStyle.h"
#include "TCanvas.h"
#include "TMath.h"
#include "TF1.h"
#include <iomanip>

#include <iostream>
#include <vector>
#include <algorithm>
#include <string>
#include <iomanip>

using namespace std;

//Gaussian fitting function
// x[0] is x value
// par[0] = amplitude
// par[1] = mean
// par[2] = sigma
double my_gaussian(double *x, double *par) {
  double arg=0;
  if(par[2]!=0) arg=(x[0]-par[1])/par[2];
  double gfunc=par[0] * TMath::Exp( -0.5*arg*arg );
  return gfunc;
}

int main( int argc, char* argv[] ) {

  if ( argc != 3 ){
    std::cerr<<"Usage: ptf_timing_analysis.app ptf_analysis.root run_number\n";
    exit(0);
  }

  // Get utilities
  Utilities utils;

  // Set style
  utils.set_style();

  std::cout<<"Input file "<<argv[1]<<std::endl;
  TFile * fin = new TFile( argv[1], "read" );

  // Opening the output root file
  string outname = string("ptf_timing_analysis_run0") + argv[2] + ".root";
  TFile * fout = new TFile(outname.c_str(), "NEW");
  std::cout<<"Output file "<<outname<<std::endl;

  // get the scanpoints information
  std::vector< ScanPoint > scanpoints = ReadScanPoints( fin );

  vector< double > xbins = utils.get_bins( scanpoints, 'x' );
  vector< double > ybins = utils.get_bins( scanpoints, 'y' );

  //Create one charge histogram per scan point
  std::vector< TH1D* > h_pmt0_tscanpt;
  std::vector< TH1D* > h_pmt1_tscanpt;

  // Loop through scanpoints
  TDirectory * dirbyscanpt = fout->mkdir( "dirbyscanpt" );
  TDirectory * curdir;
  for(unsigned int iscan=0; iscan<scanpoints.size(); iscan++){
    if ( iscan%1000 == 0 ){
      dirbyscanpt->cd();
      std::ostringstream os_curdir;
      os_curdir << iscan;
      curdir = dirbyscanpt->mkdir( os_curdir.str().c_str() );
    }
    curdir->cd();
    std::ostringstream os_dir;
    os_dir<<"scanpt_"<<iscan<<"_x="<<fixed<<setprecision(3)<<scanpoints[iscan].x()<<",y="<<scanpoints[iscan].y();
    TDirectory* dirscanpt = curdir->mkdir(os_dir.str().c_str());
    dirscanpt->cd();
    std::ostringstream os_name, os_title;
    os_name<<"h_pmt0_tscanpt_"<<iscan;
    os_title<<"PMT0: Time for signal, scan point "<<iscan;
    os_title<<"; T (ns) for X="<<fixed<<setprecision(3)<<scanpoints[iscan].x()<<" Y="<<scanpoints[iscan].y();
    os_title<<"; Counts/bin";
    if (iscan%1000==0) std::cout<<"Scan point "<<iscan<<" creating histogram: "<<os_title.str()<<std::endl;
    TH1D* htmp = new TH1D( os_name.str().c_str(), os_title.str().c_str(), 15, 40., 100. ) ;
    htmp->SetDirectory( dirscanpt );
    h_pmt0_tscanpt.push_back( htmp );

    os_name.str(""); os_name.clear();
    os_title.str(""); os_title.clear();
    os_name<<"h_pmt1_tscanpt_"<<iscan;
    os_title<<"PMT1: Time for signal, scan point "<<iscan;
    os_title<<"; T (ns) for X="<<fixed<<setprecision(3)<<scanpoints[iscan].x()<<" Y="<<scanpoints[iscan].y();
    os_title<<"; Counts/bin";
    htmp = new TH1D( os_name.str().c_str(), os_title.str().c_str(), 15, 40., 100. ) ;
    htmp->SetDirectory( dirscanpt );
    h_pmt1_tscanpt.push_back( htmp );
  }
  
  fout->cd("/");

  // Make 2d histograms
  // Transit time, transit time spread
  TH2D* h_rtt = new TH2D("h_rtt", "Transit time; X position [m]; Y position [m]", xbins.size()-1, &xbins[0], ybins.size()-1,&ybins[0]);
  TH2D* h_tts = new TH2D("h_tts", "Transit time spread; X position [m]; Y position [m]", xbins.size()-1, &xbins[0], ybins.size()-1,&ybins[0]);

  std::cout<<"Finished booking histograms"<<std::endl;

  // get the waveform fit TTree for PMT0 (The signal pmt)
  TTree * tt0 = (TTree*)fin->Get("ptfanalysis0");
  if ( !tt0 ){
    std::cerr<<"Failed to read TTree called ptfanalysis0, exiting"<<std::endl;
    return 0;
  }
  WaveformFitResult * wf0 = new WaveformFitResult;
  wf0->SetBranchAddresses( tt0 );

  // get the waveform fit TTree for PMT1 (The reference pmt)
  TTree * tt1 = (TTree*)fin->Get("ptfanalysis1");
  if ( !tt1 ){
    std::cerr<<"Failed to read TTree called ptfanalysis1, exiting"<<std::endl;
    return 0;
  }
  WaveformFitResult * wf1 = new WaveformFitResult;
  wf1->SetBranchAddresses( tt1 );
  
  //Loop through scanpoints to fill histograms
  for(unsigned int iscan=0; iscan<scanpoints.size(); iscan++){
    if (iscan%1000==0) std::cout<<"Filling histograms for iscan = "<<iscan<<" / "<<scanpoints.size()<<std::endl;
    ScanPoint scanpoint = scanpoints[ iscan ];
    //Loop over scanpoint
    for ( unsigned iev = 0; iev < scanpoint.nentries(); ++iev ){
      
      tt0->GetEvent( scanpoint.get_entry() + iev );
      if ( wf0->haswf ){
        h_pmt0_tscanpt[ iscan ]->Fill( wf0->mean );
      }
      
      tt1->GetEvent( scanpoint.get_entry() + iev );
      if ( wf1->haswf ){
        h_pmt1_tscanpt[ iscan ]->Fill( wf1->mean );
      }
    }
  }

  // Fit for each scanpoint
  std::vector< TF1* > vecpmtresponse;
  //double min_time = 70.;
  for ( unsigned iscan=0; iscan<scanpoints.size(); ++iscan ){
    if ( h_pmt0_tscanpt[iscan]->GetEntries() < 100 ) {
      vecpmtresponse.push_back( nullptr );
      //std::cout<<"Skip "<<iscan
	  //     <<" with "<<h_pmt0_tscanpt[iscan]->GetEntries()
	  //     <<" entries"<<std::endl;
      continue;
    }
    std::ostringstream fname;
    fname << "pmt_response_" << iscan;
    //std::cout<<"Fitting "<<fname.str()
	//     <<" with "<<h_pmt0_tscanpt[iscan]->GetEntries()
	//     <<std::endl;
    TF1* ftmp = new TF1( fname.str().c_str(), my_gaussian, 40., 100., 3 );
    ftmp->SetParNames("Amplitude","Mean","Sigma");
    ftmp->SetParameter(0, 20.);
    ftmp->SetParameter(1, 70.);
    ftmp->SetParameter(2, 10.);
    ftmp->SetParLimits(0, 0.0, 5000.0);
    ftmp->SetParLimits(1, 56.0, 90.0 );
    ftmp->SetParLimits(2, 0.0, 200.0 );
    h_pmt0_tscanpt[iscan]->Fit( ftmp, "Q", "", 40., 100. );
    vecpmtresponse.push_back( ftmp );
    //if( ftmp->GetParameter(1) > 25. &&
    //  ftmp->GetParameter(1) < min_time ) min_time = ftmp->GetParameter(1);
  }

  //Now fill 2d plots
  for ( unsigned iscan=0; iscan<scanpoints.size(); ++iscan){ 
    ScanPoint scanpoint = scanpoints[iscan];
    if ( h_pmt0_tscanpt[iscan]->GetEntries() >= 100 ){
      if ( vecpmtresponse[ iscan ] != nullptr ){
	    TF1* ftmp = vecpmtresponse[ iscan ];
	    h_rtt->Fill( scanpoint.x(), scanpoint.y(), ftmp->GetParameter(1) );
	    h_tts->Fill( scanpoint.x(), scanpoint.y(), ftmp->GetParameter(2) );
      }
    }
  }
  
  //Remove data outside circle
  TH2D* h_rtt_grad;
  Circle_st circ = find_circle_max_grad( h_rtt, h_rtt_grad, 0.5 );
  zero_outside_circle( h_rtt, circ );
  zero_outside_circle( h_tts, circ );

  //Set plot ranges
  h_rtt->SetMinimum(69.0);
  h_rtt->SetMaximum(80.0);
  h_tts->SetMinimum(4.4);
  //h_tts->SetMaximum(8.8);
  h_tts->SetMaximum(6.4);

  //Make plots
  TCanvas* c = new TCanvas("canvas");
  string plotname;
  h_rtt->Draw("colz0");
  plotname = string("ptf_timing_analysis_run0")+argv[2]+"_rtt.pdf";
  c->SaveAs(plotname.c_str(),"pdf");
  //TCanvas* c2 = new TCanvas();
  h_tts->Draw("colz0");
  gPad->Modified();
  gPad->Update();
  plotname = string("ptf_timing_analysis_run0")+argv[2]+"_tts.pdf";
  c->SaveAs(plotname.c_str(),"pdf");

  //Write and close output file
  fout->Write();
  fout->Close();

  cout << "Done" << endl; 

  return 0;
}
