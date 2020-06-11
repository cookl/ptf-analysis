/// PMT Response Function
/// Based on:
/// E.H. Bellamy, et. al., NIM A 339 (1994) 468.
/// and "Fitting Single Photo-electron peak," Qing He, Princeton, 2010

#include "pmt_response_function.hpp"
#include <cmath>
#include <iostream>
#include <sstream>

#include "TFile.h"
#include "TCanvas.h"
#include "TH1D.h"
#include "TMath.h"

const double pi = acos(-1.0);



/// factorial using TMath::Factorial
double factorial( unsigned n ){
  return TMath::Factorial( n );
}


/// Poisson term
double poisson_term( double mu, unsigned n ){
  return std::pow( mu, n ) * std::exp( -mu ) / factorial(n);
}

/// Gaussian part of signal for n'th photoelectron
/// This gaussian is defined only for x>0 and normalized accordingly
/// such that
///  integra( 0, infty, gaussian ) = 1
double gaussian( double x, double Qn, double sigman ){
  //double norm = 1.0 / ( sigman * std::sqrt( 2*pi ) );
  double norm = sqrt(2/pi) / ( sigman * (1 + std::erf( Qn / ( sqrt(2)*sigman ) ) ) );
  
  return norm * std::exp( -std::pow( (x-Qn)/(sqrt(2)*sigman), 2 ) );
}

double sign( double x ){
  if (x > 0.) return 1.;
  if (x < 0.) return -1.;
  return 0.;
}

/// PMT response Sreal(x) from
/// E.H. Bellamy, et. al., NIM A 339 (1994) 468.
///
/// *** Exclude Pedestal term. ***
/// x[0] is charge
/// p[0] is N, overall normalization
/// p[1] is Q1, charge of single photoelectron
/// p[2] is sigma1, width of single photoelectron gaussian
/// p[3] is mu, the poisson mean number of photoelectrons
/// p[4] is w, the fraction of the background signal that is exponential
/// p[5] is alpha, the exponential constant
double pmtresponse( double * xx, double * p ){
  double x     = xx[0];
  double N     = p[0];
  double Q1    = p[1];
  double s1    = p[2];
  double mu    = p[3];
  double w     = p[4];
  double alpha = p[5];
  
  /// how many terms to include?
  /// go to 3 sigma above mean 
  double approx_sigma = std::sqrt( mu );  
  unsigned nmax = std::max(5, int(approx_sigma*3+mu) );

  double signal = 0;
  for (unsigned n=1; n<=nmax; ++n){
    double Qn = n * Q1;
    double sigman = std::sqrt(n) * s1;
    double poisson = poisson_term( mu, n );
    double exp_term = alpha / 2 * std::exp( -alpha * ( x - Qn - alpha*sigman*sigman/2 ) ); 
    double gaus_term = gaussian( x, Qn, sigman );
    double erf1 = std::erf( fabs(n*Q1 + sigman*sigman*alpha) / ( sigman*sqrt(2.) ) );
    double erfarg = x - Qn - sigman*sigman*alpha;
    double erf2 = sign(erfarg)*std::erf( fabs(erfarg) / ( sigman*sqrt(2.) ) );
    double IGn_term = exp_term*( erf1 + erf2 );
    signal += poisson*( (1-w)*gaus_term + w*IGn_term );
  }
  return N*signal;
}

///  Background response
/// Just the background portionof the function
double pmtbackgroundresponse( double * xx, double *p ){
  double x     = xx[0];
  double N     = p[0];
  double Q1    = p[1];
  double s1    = p[2];
  double mu    = p[3];
  double w     = p[4];
  double alpha = p[5];
  
  /// how many terms to include?
  /// go to 3 sigma above mean 
  double approx_sigma = std::sqrt( mu );  
  unsigned nmax = std::max(5, int(approx_sigma*3+mu) );

  double signal = 0;
  for (unsigned n=1; n<=nmax; ++n){
    double Qn = n * Q1;
    double sigman = std::sqrt(n) * s1;
    double poisson = poisson_term( mu, n );
    double exp_term = alpha / 2 * std::exp( -alpha * ( x - Qn - alpha*sigman*sigman/2 ) ); 
    double gaus_term = gaussian( x, Qn, sigman );
    double erf1 = std::erf( fabs(n*Q1 + sigman*sigman*alpha) / ( sigman*sqrt(2.) ) );
    double erfarg = x - Qn - sigman*sigman*alpha;
    double erf2 = sign(erfarg)*std::erf( fabs(erfarg) / ( sigman*sqrt(2.) ) );
    double IGn_term = exp_term*( erf1 + erf2 );
    signal += poisson* w * IGn_term ;
  }
  return N*signal;
}

/// N'th photoelectron ideal response
/// make seventh parameter N (the photoelectron number)
double pmtG_n( double * xx, double *p ){
  double x     = xx[0];
  double N     = p[0];
  double Q1    = p[1];
  double s1    = p[2];
  double mu    = p[3];
  double w     = p[4];
  double alpha = p[5];
  double n     = p[6];
  
  double signal = 0;

  double Qn = n * Q1;
  double sigman = std::sqrt(n) * s1;
  double poisson = poisson_term( mu, n );
  double exp_term = alpha / 2 * std::exp( -alpha * ( x - Qn - alpha*sigman*sigman/2 ) ); 
  double gaus_term = gaussian( x, Qn, sigman );
  signal += poisson*( (1-w)*gaus_term );
  
  return N*signal;
}

/// Return a vector of functions that sums to be approxametely the Bellamy function
/// First term is the background, then 1 pe, 2 pe, ... etc
std::vector< TF1* > get_pmt_response_components( double * p ){
  double pars[ 7 ];
  for ( unsigned i = 0; i < 6; ++i){
    pars[i] = p[i];
  }
  pars[6]=0.;

  std::vector< TF1* > fcns;
  TF1* bg = new TF1("Background", pmtbackgroundresponse, 0., 5000., 6 );
  bg->SetParameters( pars );
  //bg->SetLineColor( 52 );
  fcns.push_back( bg );

  double mu=p[3];
  double approx_sigma = std::sqrt( mu );  
  unsigned nmax = std::max(5, int(approx_sigma*3+mu) );
  
  for (unsigned n=1; n<=nmax; ++n ){
    std::ostringstream fname;
    fname <<"fnpe_"<<n;
    TF1* ftmp = new TF1( fname.str().c_str(), pmtG_n, 0., 5000., 7 );
    pars[6] = n;
    ftmp->SetNpx(1000);
    ftmp->SetLineWidth(3);
    ftmp->SetParameters( pars );
    ftmp->SetLineColor( 52 + 5*n );
    fcns.push_back( ftmp );
  }
  return fcns;
}



/////****************************************************************
///// Repeat all of above functions but including pedestal
/////****************************************************************  
///
/// PMT response Sreal(x) from
/// E.H. Bellamy, et. al., NIM A 339 (1994) 468.
///
/// *** INclude Pedestal term. ***
/// x[0] is charge
/// p[0] is N, overall normalization
/// p[1] is Q0, charge of pedestal
/// p[2] is sigma0, width of pedestal
/// p[3] is Q1, charge of single photoelectron
/// p[4] is sigma1, width of single photoelectron gaussian
/// p[5] is mu, the poisson mean number of photoelectrons
/// p[6] is w, the fraction of the background signal that is exponential
/// p[7] is alpha, the exponential constant
double pmtresponseped( double * xx, double * p ){
  double x     = xx[0];
  double N     = p[0];
  double Q0    = p[1];
  double s0    = p[2];
  double Q1    = p[3];
  double s1    = p[4];
  double mu    = p[5];
  double w     = p[6];
  double alpha = p[7];
  
  /// how many terms to include?
  /// go to 3 sigma above mean 
  double approx_sigma = std::sqrt( mu );  
  unsigned nmax = std::max(5, int(approx_sigma*3+mu) );

  double signal = 0;

  // Use approx equation from eq. 10 of Bellamy
  double exp0term = w * std::exp( -alpha * ( x - Q0 ) );
  if ( x < Q0 ) exp0term = 0.;
  double gaus0term = (1-w) * gaussian( x, Q0, s0 );
  signal = poisson_term( mu, 0 ) * ( gaus0term + exp0term );


  for (unsigned n=1; n<=nmax; ++n){
    double gaus_term = gaussian( x, Q0 + n*Q1 + w/alpha, std::sqrt(n) * s1 );
    signal += poisson_term(mu,n) * gaus_term ;
  }
  return N*signal;
}

///  Background response *** with pedesal ***
/// Just the background portion of the function, excludes pedestal term
double pmtbackgroundresponseped( double * xx, double *p ){
  double x     = xx[0];
  double N     = p[0];
  double Q0    = p[1];
  double s0    = p[2];
  double Q1    = p[3];
  double s1    = p[4];
  double mu    = p[5];
  double w     = p[6];
  double alpha = p[7];
  
  double signal = 0;

  // Use approx equation from eq. 10 of Bellamy
  double exp0term = w * std::exp( -alpha * ( x - Q0 ) );
  if ( x < Q0 ) exp0term = 0.;
  double gaus0term = (1-w) * gaussian( x, Q0, s0 );
  signal = poisson_term( mu, 0 ) * ( gaus0term + exp0term );

  return N*signal;
}


/// N'th photoelectron ideal response
/// make ninth parameter n (the photoelectron number)
/// 0 photoelectron term is 1pe
double pmtG_n_ped( double * xx, double *p ){
  double x     = xx[0];
  double N     = p[0];
  double Q0    = p[1];
  double s0    = p[2];
  double Q1    = p[3];
  double s1    = p[4];
  double mu    = p[5];
  double w     = p[6];
  double alpha = p[7];
  unsigned n   = unsigned(p[8]);
  
  double signal = 0;

  double gaus_term = gaussian( x, Q0 + n*Q1 + w/alpha, std::sqrt(n) * s1 );
  signal += poisson_term(mu,n) * gaus_term ;

  return N*signal;
}

/// Return a vector of functions that sums to be approxametely the Bellamy function
/// First term is the background, then 0 pe, 1 pe, 2 pe, ... etc
std::vector< TF1* > get_pmt_response_components_ped( double * p ){
  double pars[ 9 ];
  for ( unsigned i = 0; i < 8; ++i){
    pars[i] = p[i];
  }
  pars[8]=0.;

  std::vector< TF1* > fcns;
  TF1* bg = new TF1("Background", pmtbackgroundresponseped, 0., 5000., 8 );
  bg->SetParameters( pars );
  //bg->SetLineColor( 52 );
  fcns.push_back( bg );

  double mu=p[5];
  double approx_sigma = std::sqrt( mu );  
  unsigned nmax = std::max(5, int(approx_sigma*3+mu) );
  
  for (unsigned n=1; n<=nmax; ++n ){
    std::ostringstream fname;
    fname <<"fnpe_"<<n;
    TF1* ftmp = new TF1( fname.str().c_str(), pmtG_n_ped, 0., 5000., 9 );
    pars[8] = n;
    ftmp->SetNpx(1000);
    ftmp->SetLineWidth(3);
    ftmp->SetParameters( pars );
    ftmp->SetLineColor( 52 + 5*n );
    fcns.push_back( ftmp );
  }
  return fcns;
}


PMTResponsePed * PMTResponsePed::fInstance = nullptr;
TH1D* PMTResponsePed::fPDF = nullptr;
double PMTResponsePed::fWid = 0;
TFile* PMTResponsePed::fin = nullptr;

void PMTResponsePed::set_pedestal( const std::string& fname, const std::string& histname ){
  if ( fInstance == nullptr ) fInstance = new PMTResponsePed( 0.0 );
  fWid = 0.0;
  if ( fin == nullptr ) fin = new TFile( fname.c_str(), "read" );
  if ( !fin ) {
    std::cout<<"Could not find background root file: "<<fname<<std::endl;
    exit(0);
  }
  fPDF = (TH1D*) fin->Get( histname.c_str() );
  if ( !fPDF ){
    std::cout<<"Could not find histogram: "<<histname<<" in root file: "<<fname<<std::endl;
    exit(0);
  } 
}


/*
Model 1:
x[0] = charge
p[0] = normalization (count)
p[1] = q1 single pe mean charge
p[2] = s1 single pe charge width
p[3] = mu mean number of pe
p[4] = w weight of exponential bg (0-1)
p[5] = alpha exponential constant bg

Change background to be part of pedestal signal!
 */
double model1( double * x, double *p ){
  double N = p[0];
  double q1 = p[1];
  double s1 = p[2];
  double mu = p[3];
  double w = p[4];
  double alpha = p[5];

  double retval = 0.0;

  retval = w * alpha * exp( -alpha*x[0] );
  
  retval += (1-w)*poisson_term( mu, 0 ) * PMTResponsePed::get_prob_density( x[0] );
  
  for ( unsigned npe =1; npe<6; ++npe ){
    retval += (1-w)*poisson_term( mu, npe ) * gaussian( x[0], npe*q1, sqrt( npe )*s1 )
      ;
  }

  return retval * N;
}


double model1bg( double * x, double *p ){
  double N = p[0];
  double q1 = p[1];
  double s1 = p[2];
  double mu = p[3];
  double w = p[4];
  double alpha = p[5];
  
  double retval = 0.0;

  //double binwid = PMTResponsePed_BinWid::GetBinWid();
  //double expnorm = alpha / (1.0 - exp( -alpha*binwid ) );
  //if ( x[0] < binwid ) retval = (1-w)*poisson_term( mu, 0 )/binwid;


  retval = w * alpha * exp( -alpha*x[0] ); 
  retval += (1-w)* poisson_term( mu, 0 ) * PMTResponsePed::get_prob_density( x[0] );

  
  return retval * N;
}


double model1npe( double * x, double *p ){
  double N = p[0];
  double q1 = p[1];
  double s1 = p[2];
  double mu = p[3];
  double w = p[4];
  double alpha = p[5];
  
  unsigned npe = p[6];
  
  double retval = 0.0;

  retval += (1-w) * poisson_term( mu, npe ) * gaussian( x[0], npe*q1, sqrt( npe ) * s1 );

  return retval * N;
}



/// Return a vector of functions that sums to be approxametely model1
/// First term is 0 pe + bg, 1 pe, 2 pe, ... etc
std::vector< TF1* > get_model1_components( double * p ){
  double pars[ 7 ];
  for ( unsigned i = 0; i < 6; ++i){
    pars[i] = p[i];
  }
  pars[6]=0.;

  std::vector< TF1* > fcns;
  TF1* bg = new TF1("Background", model1bg, 0., 5000., 6 );
  bg->SetParameters( pars );
  //bg->SetLineColor( 52 );
  fcns.push_back( bg );
  
  for (unsigned n=1; n<6; ++n ){
    std::ostringstream fname;
    fname <<"fnpe_"<<n;
    TF1* ftmp = new TF1( fname.str().c_str(), model1npe, 0., 5000., 7 );
    pars[6] = n;
    ftmp->SetLineWidth(3);
    ftmp->SetParameters( pars );
    ftmp->SetLineColor( 52 + 5*n );
    fcns.push_back( ftmp );
  }
  for ( TF1 * fcn : fcns ) fcn->SetNpx(1000);
  return fcns;
}




/// Test case code to test the fitting ***with pedestal***
void pmt_response_function(){
  
  // try it out
  TFile* fin = new TFile( "ptf_charge_analysis_run04525.root", "read" );

  TH1D* h = (TH1D*) fin->Get("hqall");//get_sum_histogram( fin, "hsum_4525" );
  
  //TH1D* h=(TH1D*)fin->Get("hqall");
  h->SetLineWidth(3);
  //h->Rebin(2);

  PMTResponsePed::set_binwid( h->GetBinWidth(1) );
  //PMTResponsePed_BinWid::SetBinWid( h->GetBinWidth(1) );

  TCanvas *c=new TCanvas();
  c->cd();
  h->Draw("e2");
  h->SetMarkerStyle(20);
  h->SetMarkerSize(0.5);

  double Nfix = h->Integral( 1, h->GetNbinsX()+1, "width" );

  TF1* ff = new TF1( "pmt_response", model1, 0., 5000., 6 );

  double N0 = h->Integral(1,1);
  double Nrest = h->Integral(2, 15 ); // beyond bin 15 dominated by exp background
  double mufix = Nrest / N0;
  mufix = log( mufix + 1 ); // correction factor to go from estimated mu, to true

  //cout << "Nfix = "<<Nfix<<" Nrest="<<Nrest<<" / N0 = "<<N0<<" = "<<mufix<<std::endl;
  //cout << "mufix = " << mufix << std::endl;
                           
  ff->SetNpx(1000);
  ff->SetParNames("N","Q_{1}","#sigma_{1}", "#mu", "w", "#alpha" );
  ff->SetParameters( Nfix, 397., 148., mufix, 0.003, 0.000835  );
  ff->FixParameter( 0, Nfix );
  ff->FixParameter( 1, 400.0 );
  ff->FixParameter( 2, 148. );
  ff->FixParameter( 3, mufix );

  h->Fit(ff, "", "", 2000., 5000.0 );

  for ( unsigned ipar=0; ipar < 6; ++ipar ) ff->ReleaseParameter(ipar);
  //ff->FixParameter( 0, Nfix );
  //ff->FixParameter( 4, ff->GetParameter(4) ) ;
  //ff->FixParameter( 5, ff->GetParameter(5) );

  ff->SetLineWidth(3);
  ff->SetLineColor(kRed+2);


  h->Fit(ff, "", "", 0., 5000. );
  //h->Draw();
  //ff->Draw("same");


  //ff->SetParameter(1,400);
  //ff->SetParameter(2,80);
  
  //h->Draw("same");
  std::vector< TF1* > fcmp = get_model1_components( ff->GetParameters() );

  for ( TF1* f : fcmp ){
    f->Draw("same");
  }
  
  return;
}
