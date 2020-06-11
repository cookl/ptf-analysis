#ifndef _PMTResponseFunction_hpp_
#define _PMTResponseFunction_hpp_

#include <vector>
#include "TF1.h"
#include "TH1D.h"
#include "TFile.h"

/// Functions to allow fitting to PMT Response function of
/// E.H. Bellamy, et. al., NIM A 339 (1994) 468.
/// and "Fitting Single Photo-electron peak," Qing He, Princeton, 2010

/// PMT response Sreal(x) from
/// E.H. Bellamy, et. al., NIM A 339 (1994) 468.
///
/// Exclude Pedestal term.
/// x[0] is charge
/// p[0] is N, overall normalization
/// p[1] is Q1, charge of single photoelectron
/// p[2] is sigma1, width of single photoelectron gaussian
/// p[3] is mu, the poisson mean number of photoelectrons
/// p[4] is w, the fraction of the background signal that is exponential
/// p[5] is alpha, the exponential constant
double pmtresponse( double * xx, double * p );

///  Background response
/// Just the background portionof the function
/// Same parameter arguments as pmtresponse
double pmtbackgroundresponse( double * xx, double *p );

/// N'th photoelectron ideal response
/// make seventh parameter N (the photoelectron number)
/// Same parameter arguments as pmtresponse
double pmtG_n( double * xx, double *p );

/// Return a vector of functions that sums to be the Bellamy function
/// First term is the background, then 1 pe, 2 pe, ... etc
std::vector< TF1* > get_pmt_response_components( double * p );


/// Repeat all of above functions but include pedestal
/// E.H. Bellamy, et. al., NIM A 339 (1994) 468.
/// and "Fitting Single Photo-electron peak," Qing He, Princeton, 2010

/// PMT response Sreal(x) from
/// E.H. Bellamy, et. al., NIM A 339 (1994) 468.
///
/// INclude Pedestal term.
/// x[0] is charge
/// p[0] is N, overall normalization
/// p[1] is Q0, charge of pedestal
/// p[2] is sigma0, width of pedestal
/// p[3] is Q1, charge of single photoelectron
/// p[4] is sigma1, width of single photoelectron gaussian
/// p[5] is mu, the poisson mean number of photoelectrons
/// p[6] is w, the fraction of the background signal that is exponential
/// p[7] is alpha, the exponential constant
double pmtresponseped( double * xx, double * p );

///  Background response
/// Just the background portionof the function
/// Same parameter arguments as pmtresponseped
double pmtbackgroundresponseped( double * xx, double *p );

/// N'th photoelectron ideal response
/// make ninth parameter N (the photoelectron number)
/// Same parameter arguments as pmtresponseped
double pmtG_n_ped( double * xx, double *p );

/// Return a vector of functions that sums to be the Bellamy function
/// First term is the background, then 0 pe, 1 pe, 2 pe, ... etc
std::vector< TF1* > get_pmt_response_components_ped( double * p );


/// For pedestal and background from dark run histogram
/// Set nonzero binwidth to use bin width instead
class PMTResponsePed {
public:
  static double get_prob_density( double x ) {
    if ( fInstance == nullptr ) fInstance = new PMTResponsePed( "pedestal_from_run4554.root" );
    if ( fWid > 0.0 ){
      if ( x < fWid ) return 1.0/fWid;
      else return 0.0;
    } else {
      int ibin = fPDF->FindBin( x );
      return fPDF->GetBinContent( ibin );
    }
  }

  static void set_pedestal( const std::string& fname, const std::string& histname );

  static void set_binwid( double wid ){
    if ( fInstance == nullptr ) fInstance = new PMTResponsePed( wid );
    else fWid = wid;
  }
      
protected:
  PMTResponsePed( double wid){
    fWid = wid;
  }
  
  PMTResponsePed( std::string filename ) {
    fWid = 0.0;
    fin = new TFile( filename.c_str(), "read" );
    if ( !fin ) {
      std::cout<<"Could not find background root file: "<<filename<<std::endl;
      exit(0);
    }
    fPDF = (TH1D*) fin->Get("nofftcut");
    if ( !fPDF ){
      std::cout<<"Could not find nofftcut in root file: "<<filename<<std::endl;
      exit(0);
    }
  }
  static double fWid;
  static PMTResponsePed * fInstance;
  static TH1D * fPDF;
  static TFile* fin;
};


///
/// Model 1: arXiv 1710.00219
/// Gaussians for each photoelectron, assumed to dissappear below x=0
/// Poisson distribution decides strength of each gaussian
/// exponential background starts after pedestal
/// pedestal is Poisson 0 term, but is uniform distr. in first bin
///
/// Uses PMTResponsPed_BinWid::GetBinWid() to determine width of pedestal bin 
///
/// x[0] = charge
/// p[0] = normalization (count)
/// p[1] = q1 single pe mean charge
/// p[2] = s1 single pe charge width
/// p[3] = mu mean number of pe
/// p[4] = w weight of exponential
/// p[5] = alpha decay constant of exponential
/// 
double model1( double * x, double *p );


/// Return a vector of functions that sums to be approxametely model1
/// First term is the background + 0 pe, 1 pe, 2 pe, ... etc
std::vector< TF1* > get_model1_components( double * p );



#endif
