/* ptf_charge_analysis.cpp
   Analysis of fitted waveforms to determine the charge of the events.
   Produces histograms of total charge in the waveform using the gaussian fit results.
   The Charge histogram is made once per scan point, and a sum one for all events.
   The Charge histogram is fit to a model to determine the locations and strengths of the 1pe, 2pe, ... 

   Charge of 1pe is plotted as function of scan point (TH2D with x,y as axes, charge as entry).
   Difference in charge between 1 pe and 2 pe is also plotted as function of scan point (TH2D with x,y as axes, diff as entry).
   Finally, plot the average number of pe.
   
   Author: Blair Jamieson (Jan 2020)
 */

#include "pmt_response_function.hpp"

#include "WaveformFitResult.hpp"
#include "ScanPoint.hpp"
#include "Utilities.hpp"
#include "TFile.h"
#include "TH2D.h"
#include "TStyle.h"
#include "TCanvas.h"
#include "TFitResult.h"
#include "TFitResultPtr.h"
#include "TMatrixDSym.h"
#include "FindCircle.hpp"

#include <iostream>
#include <vector>
#include <algorithm>
#include <string>

#include <sstream>

using namespace std;

const double sqrt2pi = std::sqrt( 2 * std::acos(-1) );

double calculate_charge( const WaveformFitResult& wf ){
  // integral of gaussian a * exp( -(x-b)^2 / 2 c^2 )  = sqrt(2*pi) a * c
  return sqrt2pi * fabs( wf.amp * wf.sigma );
} 

int main( int argc, char* argv[] ) {

  if ( argc != 3 && argc != 4){
    std::cerr<<"Usage: ptf_charge_analysis.app ptf_analysis.root run_number [T/F/I]"<<std::endl;
    std::cerr<<"Where the T/F/I is for True to do/not do circle fit to find PMT, I to cut inside circle (default T)\n"<<std::endl;
    exit(0);
  }

  // Get utilities
  Utilities utils;

  // Set style
  utils.set_style();

  std::cout<<"Input file "<<argv[1]<<std::endl;
  TFile * fin = new TFile( argv[1], "read" );

  // Opening the output root file
  string outname = string("ptf_charge_analysis_run0") + argv[2] + ".root";
  TFile * fout = new TFile(outname.c_str(), "NEW");
  std::cout<<"Output file "<<outname<<std::endl;

  // get the scanpoints information
  std::vector< ScanPoint > scanpoints = ReadScanPoints( fin );

  vector< double > xbins = utils.get_bins( scanpoints, 'x' );
  vector< double > ybins = utils.get_bins( scanpoints, 'y' );

  // Create overall sum charge histogram
  TH1D* hqall = new TH1D("hqall", "Charge deposit for all events; Q (ADC); Counts/bin", 50, 0., 5000. );
  TH1D* hqallfine = new TH1D("hqallfine", "Charge deposit for all events; Q (ADC); Counts/bin", 100, 0., 200. );
  TH1D* hqsum = new TH1D("hqsum", "Charge deposit for good fit events; Q (ADC); Counts/bin", 100, 0., 5000. );
  TH1D* hqped = new TH1D("hqped", "Charge deposit for pedestal events; Q (ADC); Counts/bin", 100, 0., 5000. );

  TH1D* hchi2 = new TH1D("hchi2", "#chi^{2} of scan-point fits; #chi^{2}; scan-points / bin", 100, 0., 50. );

  // Pulse height spectrum
  TH1D* hphall = new TH1D("hphall", "Pulse height for all events; Pulse Height (ADC); Counts/bin", 100, 0., 500. );
  // Create one charge histogram per scan point
  std::vector< TH1D* > hqallscanpt;
  std::vector< TH1D* > hqallfinescanpt;
  std::vector< TH1D* > hqscanpt;
  std::vector< TH1D* > hpedscanpt;
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
    os_dir<<"scanpt_"<<iscan;
    TDirectory* dirscanpt = curdir->mkdir(os_dir.str().c_str());
    dirscanpt->cd();
    std::ostringstream os_name, os_title;
    os_name<<"hqscanpt_"<<iscan;
    os_title<<"Charge for signal, scan point "<<iscan;
    os_title<<"; Q (ADC) for X="<<scanpoints[iscan].x()<<" Y="<<scanpoints[iscan].y();
    os_title<<"; Counts/bin";
    if (iscan%100==0) std::cout<<"Scan point "<<iscan<<" creating histogram: "<<os_title.str()<<std::endl;
    TH1D* htmp = new TH1D( os_name.str().c_str(), os_title.str().c_str(), 50, 0., 5000. ) ;
    htmp->SetDirectory( dirscanpt );
    hqscanpt.push_back( htmp );

    os_name.str(""); os_name.clear();
    os_title.str(""); os_title.clear();
    os_name<<"hpedscanpt_"<<iscan;
    os_title<<"Charge for Pedestal, scan point "<<iscan;
    os_title<<"; Q (ADC) for X="<<scanpoints[iscan].x()<<" Y="<<scanpoints[iscan].y();
    os_title<<"; Counts/bin";
    htmp = new TH1D( os_name.str().c_str(), os_title.str().c_str(), 50, 0., 5000. ) ;
    htmp->SetDirectory( dirscanpt );
    hpedscanpt.push_back( htmp );

    os_name.str(""); os_name.clear();
    os_title.str(""); os_title.clear();
    os_name<<"hqallscanpt_"<<iscan;
    os_title<<"Charge for all, scan point "<<iscan;
    os_title<<"; Q (ADC) for X="<<scanpoints[iscan].x()<<" Y="<<scanpoints[iscan].y();
    os_title<<"; Counts/bin";
    htmp = new TH1D( os_name.str().c_str(), os_title.str().c_str(), 50, 0., 5000. ) ;
    htmp->SetDirectory( dirscanpt );
    hqallscanpt.push_back( htmp );

    os_name.str(""); os_name.clear();
    os_title.str(""); os_title.clear();
    os_name<<"hqallfinescanpt_"<<iscan;
    os_title<<"Charge for all, scan point "<<iscan;
    os_title<<"; Q (ADC) for X="<<scanpoints[iscan].x()<<" Y="<<scanpoints[iscan].y();
    os_title<<"; Counts/bin";
    htmp = new TH1D( os_name.str().c_str(), os_title.str().c_str(), 50, 0., 200. ) ;
    htmp->SetDirectory( dirscanpt );
    hqallfinescanpt.push_back( htmp );
  }
  
  fout->cd("/");

  // Make 2d histograms
  // average charge, 1pe charge, 2pe charge, average number of pe, 2 pe - 1 pe
  TH2D* hscanpt = new TH2D("hscanpt", "Scan point; X position [m]; Y position [m]", xbins.size()-1, &xbins[0], ybins.size()-1,&ybins[0]);
  TH2D* hqavg = new TH2D("hqavg", "Average charge; X position [m]; Y position [m]", xbins.size()-1, &xbins[0], ybins.size()-1,&ybins[0]);
  TH2D* hq1pe = new TH2D("hq1pe", "Q_{1} pe charge; X position [m]; Y position [m]", xbins.size()-1, &xbins[0], ybins.size()-1,&ybins[0]);
  TH2D* hq2pe = new TH2D("hq1sig", "#sigma_{1} charge width; X position [m]; Y position [m]", xbins.size()-1, &xbins[0], ybins.size()-1,&ybins[0]);
  TH2D* hq21pe = new TH2D("hqmupe", "#mu average npe; X position [m]; Y position [m]", xbins.size()-1, &xbins[0], ybins.size()-1,&ybins[0]);

  TH2D* hngoodfits = new TH2D("hngoodfits", "N signals; X position [m]; Y position [m]", xbins.size()-1, &xbins[0], ybins.size()-1,&ybins[0]);

  TH2D* hq1pemu = new TH2D("hq1pemu", "Q_{1} x #mu ; X position [m]; Y position [m]", xbins.size()-1, &xbins[0], ybins.size()-1,&ybins[0]);


  // Histogram of correlation matrix
  TH2D* hcormat = new TH2D("hcormat","Bellamy fit correlations", 6, -0.5, 5.5, 6, -0.5, 5.5 );
  TAxis * xax = hcormat->GetXaxis();
  TAxis * yax = hcormat->GetYaxis();
  //ff->SetParNames("N","Q_{1}","#sigma_{1}", "#mu", "w", "#alpha" );

  xax->SetBinLabel( 1, "N" );
  xax->SetBinLabel( 2, "Q_{1}" );
  xax->SetBinLabel( 3, "#sigma_{1}" );
  xax->SetBinLabel( 4, "#mu" );
  xax->SetBinLabel( 5, "w" );
  xax->SetBinLabel( 6, "#alpha" );

  yax->SetBinLabel( 1, "N" );
  yax->SetBinLabel( 2, "Q_{1}" );
  yax->SetBinLabel( 3, "#sigma_{1}" );
  yax->SetBinLabel( 4, "#mu" );
  yax->SetBinLabel( 5, "w" );
  yax->SetBinLabel( 6, "#alpha" );

  // Histogram of correlation matrix
  TH2D* hcormat2 = new TH2D("hcormat2","Bellamy fit correlations", 6, -0.5, 5.5, 6, -0.5, 5.5 );
  xax = hcormat2->GetXaxis();
  yax = hcormat2->GetYaxis();

  xax->SetBinLabel( 1, "N" );
  xax->SetBinLabel( 2, "Q_{1}" );
  xax->SetBinLabel( 3, "#sigma_{1}" );
  xax->SetBinLabel( 4, "#mu" );
  xax->SetBinLabel( 5, "w" );
  xax->SetBinLabel( 6, "#alpha" );

  yax->SetBinLabel( 1, "N" );
  yax->SetBinLabel( 2, "Q_{1}" );
  yax->SetBinLabel( 3, "#sigma_{1}" );
  yax->SetBinLabel( 4, "#mu" );
  yax->SetBinLabel( 5, "w" );
  yax->SetBinLabel( 6, "#alpha" );

  
  std::cout<<"Finished booking histograms"<<std::endl;

  // get the waveform fit TTree for PMT1 (The signal pmt)
  TTree * tt1 = (TTree*)fin->Get("ptfanalysis0");
  if ( !tt1 ){
    std::cerr<<"Failed to read TTree called ptfanalysis0, exiting"<<std::endl;
    return 0;
  }
  WaveformFitResult * wf = new WaveformFitResult;
  wf->SetBranchAddresses( tt1 );
  
  // First loop through scanpoints to fill a few histograms
  for(unsigned int iscan=0; iscan<scanpoints.size(); iscan++){
    if (iscan%100==0) std::cout<<"pass 1: Filling histograms for iscan = "<<iscan<<" / "<<scanpoints.size()<<std::endl;
    ScanPoint scanpoint = scanpoints[ iscan ];
    //Loop over scanpoint
    for ( unsigned iev = 0; iev < scanpoint.nentries(); ++iev ){
      
      tt1->GetEvent( scanpoint.get_entry() + iev );
      double charge = calculate_charge( *wf );
      hqallscanpt[ iscan ]->Fill( charge );
    }
  }

  // Build 2d average charge plot from hqallscanpt histograms
  for ( unsigned iscan=0; iscan<scanpoints.size(); ++iscan){ 
    ScanPoint scanpoint = scanpoints[iscan];
    hscanpt->Fill( scanpoint.x(), scanpoint.y(), float(iscan) );
    if ( hqallscanpt[iscan]->Integral(2,50) > 200 ){
      hqavg->Fill( scanpoint.x(), scanpoint.y(), hqallscanpt[iscan]->GetMean() );
    }
  }
  
  // Find PMT location from hqavg histogram
  hqavg->Write();
  TH2D* hgrad;
  Circle_st circ;

  bool docirclefit = true;
  if ( argc == 4 && argv[3][0] == 'F' ) docirclefit = false;
  if ( docirclefit ){
    circ = find_circle_max_grad( hqavg, hgrad, 0.25 );
  } else {
    circ.r  = 0.25;
    circ.xc = 0.40;
    circ.yc = 0.36;
  }
  std::cout<<"==================================="<<std::endl;
  std::cout<<"Using circle fit? "<<docirclefit<<std::endl;
  std::cout<<"PMT Located at : "<<std::endl;
  std::cout<<"\t (xc,yc) = ( "<<circ.xc<<", "<<circ.yc<<" )"<<std::endl;
  std::cout<<"\t R       =   "<<circ.r<<std::endl;
  std::cout<<"==================================="<<std::endl;


  //Second pass through scanpoints to fill histograms inside circle of PMT
  for(unsigned int iscan=0; iscan<scanpoints.size(); iscan++){
    if (iscan%100==0) std::cout<<"pass 2 filling histograms for iscan = "<<iscan<<" / "<<scanpoints.size()<<std::endl;

    ScanPoint scanpoint = scanpoints[ iscan ];

    if ( argc == 4 && argv[3][0] == 'I' ) {// cut inside instead of outside
      if ( circ.is_inside( scanpoint.x(), scanpoint.y() ) ) {
	std::cout<<"Skip scan point "<<iscan<<" inside circle " <<std::endl;
	continue;
      }
    } else {
      if ( !circ.is_inside( scanpoint.x(), scanpoint.y() ) ) {
	std::cout<<"Skip scan point "<<iscan<<" outside circle " <<std::endl;
	continue;
      }
    }
    
    //Loop over scanpoint
    for ( unsigned iev = 0; iev < scanpoint.nentries(); ++iev ){
      
      tt1->GetEvent( scanpoint.get_entry() + iev );
      double charge = calculate_charge( *wf );
      if ( wf->haswf ) {
	hqscanpt[ iscan ]->Fill( charge );
	hqsum->Fill( charge );
      } else {
	hqped->Fill( charge );
	hpedscanpt[ iscan ]->Fill( charge );
      }

      hphall->Fill( wf->amp );
      hqall->Fill( charge );
      hqallfine->Fill( charge );
      hqallfinescanpt [ iscan ]->Fill(charge);

    }
  }

  // init binwidth
  PMTResponsePed::set_binwid( hqall->GetBinWidth(1) );
  //PMTResponsePed::set_pedestal( "nofftcut_pedestal_4554.root", "hqall_nofftcut" );

  // Do initial fit to the overall charge histogram
  // With everything freely floating
  TF1* ff = new TF1( "pmt_response", model1, 0., 5000., 6 );
  double Nfix  = hqall->Integral( 1, hqall->GetNbinsX()+1, "width" );
  double N0    = hqall->Integral( 1, 1 );
  double Nrest = hqall->Integral( 2, hqall->GetNbinsX()+1 );
  double mufix = Nrest / N0;
  mufix = log( mufix + 1 ); // corrected to get estimated mu
  ff->SetNpx(1000);
  ff->SetParNames("N","Q_{1}","#sigma_{1}", "#mu", "w", "#alpha" );
  ff->SetParameters( Nfix, 400.0, 152.0, mufix, 0.003, 0.000835  );
  ff->FixParameter( 0, Nfix );
  ff->FixParameter( 1, 400.0 );
  ff->FixParameter( 2, 148. );
  ff->FixParameter( 3, mufix );
  hqall->Fit(ff, "", "", 2000., 5000.0 );

  for ( unsigned ipar=0; ipar < 6; ++ipar ) ff->ReleaseParameter(ipar);
  ff->FixParameter( 0, Nfix );
  ff->FixParameter( 4, ff->GetParameter(4) ) ;
  ff->FixParameter( 5, ff->GetParameter(5) );

  ff->SetLineWidth(3);
  ff->SetLineColor(kRed+2);

  hqall->Fit(ff, "", "", 0., 2000. );

  for ( unsigned ipar=0; ipar < 6; ++ipar ) ff->ReleaseParameter(ipar);
  TCanvas * cpmtres = new TCanvas();
  cpmtres->cd();
  TFitResultPtr  result = hqall->Fit(ff, "S", "", 0, 5000.0 );
  TMatrixDSym cov = result->GetCovarianceMatrix();
  TMatrixDSym cor = result->GetCorrelationMatrix();

  std::cout<<"Fit covariance matrix is"<<std::endl;
  cov.Print();
  std::cout<<"Fit correlation matrix is"<<std::endl;
  cor.Print();

  for (unsigned i=0; i<6 ;++i){
    for (unsigned j=0; j<6 ;++j){
      hcormat->SetBinContent( i+1, j+1, cor[i][j] );
    }
  }


  // Repeat initial fit to the overall charge histogram
  // This time fixing what needs to be fixed

  TF1* ff2 = new TF1( "pmt_response2", model1, 0., 5000., 6 );
  ff2->SetParNames("N","Q_{1}","#sigma_{1}", "#mu", "w", "#alpha" );
  ff2->SetParameters( Nfix, 400.0, 152.0, mufix, 0.003, 0.000835  );
  ff2->SetNpx(1000);
  ff2->SetLineWidth(3);
  ff2->SetLineColor(kRed+2);
  ff2->FixParameter(0, Nfix );
  ff2->FixParameter(1, 400.0 );
  ff2->FixParameter(2, 148.0 );
  ff2->FixParameter(3, mufix );

  hqall->Fit(ff2, "", "", 2000., 5000.0 );

  for ( unsigned ipar=0; ipar < 6; ++ipar ) ff2->ReleaseParameter(ipar);
  ff2->FixParameter( 0, Nfix );
  ff2->FixParameter( 4, ff2->GetParameter(4) ) ;
  ff2->FixParameter( 5, ff2->GetParameter(5) );

  ff2->SetLineWidth(3);
  ff2->SetLineColor(kRed+2);

  hqall->Fit(ff2, "", "", 0., 2000. );

  // Initial fit to pin the parameters we want
  TCanvas * cpmtres2 = new TCanvas();
  cpmtres2->cd();

  // Fix up the fit to the background
  for(unsigned i=0; i<6; ++i ) ff2->ReleaseParameter(i);
  ff2->FixParameter( 0, Nfix  );
  ff2->FixParameter( 4, ff2->GetParameter(4) );
  ff2->FixParameter( 5, ff2->GetParameter(5) );

  TFitResultPtr  result2 = hqall->Fit(ff2, "S", "", 0, 5000.0 );

  TMatrixDSym cov2 = result2->GetCovarianceMatrix();
  TMatrixDSym cor2 = result2->GetCorrelationMatrix();

  std::cout<<"All Bellamy function fit covariance matrix is"<<std::endl;
  cov2.Print();
  std::cout<<"All Bellamy function fit correlation matrix is"<<std::endl;
  cor2.Print();
  
  for (unsigned i=0; i<8 ;++i){
    for (unsigned j=0; j<8 ;++j){
      hcormat2->SetBinContent( i+1, j+1, cor2[i][j] );
    }
  }
  

  std::cout<<"Done global fit"<<std::endl;

  // add fit components to plot
  std::vector< TF1* > fcmp = get_model1_components( ff2->GetParameters() );
  for ( TF1* f : fcmp ){
    hqall->GetListOfFunctions()->Add( f );
    f->Draw("same");
  }
  hqall->Write();


  std::cout<<"Now fit each scanpoint"<<std::endl;
  //PMTResponsePed::set_pedestal( "nofftcut_pedestal_4554.root", "nofftcut" );

  // Fit for each scanpoint
  std::vector< TF1* > vecpmtresponse;
  for ( unsigned iscan=0; iscan<scanpoints.size(); ++iscan ){

    ScanPoint scanpoint = scanpoints[ iscan ];

    if ( argc == 4 && argv[3][0] == 'I' ) {// cut inside instead of outside
      if ( circ.is_inside( scanpoint.x(), scanpoint.y() ) ) {
	vecpmtresponse.push_back( nullptr ); 
	std::cout<<"Skip scan point "<<iscan<<" inside circle " <<std::endl;
	continue;
      }
    } else {
      if ( !circ.is_inside( scanpoint.x(), scanpoint.y() ) ) {
	vecpmtresponse.push_back( nullptr ); 
	std::cout<<"Skip scan point "<<iscan<<" outside circle " <<std::endl;
	continue;
      }
    }

    std::ostringstream fname;
    fname << "pmt_response_" << iscan;
    std::cout<<"Fitting "<<fname.str()
	     <<" with "<<hqallscanpt[iscan]->GetEntries()
	     <<std::endl;

    PMTResponsePed::set_binwid( hqallscanpt[iscan]->GetBinWidth(1) );

    double curNfix  = hqallscanpt[iscan]->Integral(1, hqallscanpt[iscan]->GetNbinsX()+1, "width" );
    double curN0    = hqallscanpt[iscan]->Integral(1,1);
    double curNrest = hqallscanpt[iscan]->Integral(2, hqallscanpt[iscan]->GetNbinsX()+1 );
    double curmufix = curNrest / curN0;
    curmufix = log( curmufix + 1 );

    TF1* ftmp = new TF1( fname.str().c_str(), model1, 0., 5000., 6 );
    ftmp->SetNpx(1000);
    ftmp->SetParNames("N","Q_{1}","#sigma_{1}", "#mu", "w", "#alpha" );
    ftmp->FixParameter(0, curNfix );
    ftmp->SetParameter(1, ff2->GetParameter(1) );
    ftmp->SetParameter(2, ff2->GetParameter(2) );
    ftmp->FixParameter(3, curmufix );
    ftmp->FixParameter(4, ff2->GetParameter(4) );
    ftmp->FixParameter(5, ff2->GetParameter(5) );
    ftmp->SetParLimits(1, 0., 1000. );
    ftmp->SetParLimits(2, 0., 1000. );
    hqallscanpt[iscan]->Fit( ftmp, "", "", 0., 2000. );
    vecpmtresponse.push_back( ftmp );
    hchi2->Fill( ftmp->GetChisquare() );
  }


  //Now fill 2d plots
  for ( unsigned iscan=0; iscan<scanpoints.size(); ++iscan){ 
    ScanPoint scanpoint = scanpoints[iscan];
    hngoodfits->Fill( scanpoint.x(), scanpoint.y(), hqscanpt[iscan]->GetEntries() );
    if ( vecpmtresponse[ iscan ] != nullptr ){
      TF1* ftmp = vecpmtresponse[ iscan ];
      hq1pe->Fill( scanpoint.x(), scanpoint.y(), ftmp->GetParameter(1) );
      hq2pe->Fill( scanpoint.x(), scanpoint.y(), ftmp->GetParameter(2) );
      hq21pe->Fill( scanpoint.x(), scanpoint.y(), ftmp->GetParameter(3) );
      hq1pemu->Fill( scanpoint.x(), scanpoint.y(), ftmp->GetParameter(1)*ftmp->GetParameter(3) );
    }
  }

  //Write and close output file
  //PMTResponsePed::set_pedestal( "nofftcut_pedestal_4554.root", "hqall_nofftcut" );
  fout->Write();
  fout->Close();

  cout << "Done" << endl; 

  return 0;
}
