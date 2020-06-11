#include "FindCircle.hpp"
#include "Hough.hpp"
#include "HoughDisplay.hpp"

#include <sstream>

#include "TAxis.h"
#include "TCanvas.h"
#include "TFile.h"
#include "TArc.h"

#include <iostream>

using namespace std;

bool Circle_st::is_inside( double x, double y ) const{
  return ( (x-xc)*(x-xc) + (y-yc)*(y-yc) < r*r );
}



void th2d_to_xypointvector( const TH2D* in, std::vector<xypoint>& out ){
  // find binning
  const TAxis* xax = in->GetXaxis();
  unsigned nbinsx = xax->GetNbins();
  const TAxis* yax = in->GetYaxis();
  unsigned nbinsy = yax->GetNbins();

  // loop over bins
  for (unsigned ix=1; ix<nbinsx; ++ix ){
    double x = xax->GetBinCenter( ix );
    for (unsigned iy=1; iy<nbinsy; ++iy ){
      double y = yax->GetBinCenter( iy );
      if ( in->GetBinContent(ix,iy ) > 0.0 ) {
	out.push_back( xypoint( x,y ) );
      }
    }
  }   
}



// Take an input TH2D histogram, and compute the gradient (largest change in value between 8 nearest neighbours)
// and fill that in the gradhist that is built in this function.  Set any gradient values that are below
// cut_below * maximum_gradient to zero, then find circle using hough transform
Circle_st find_circle_max_grad( const TH2D* in, TH2D*& grad, double cut_below ){
  Circle_st circle;
  
  // build output gradient histogram
  ostringstream grad_name;
  grad_name << in->GetName() <<"_grad";
  grad = (TH2D*)in->Clone( grad_name.str().c_str() );
  grad->SetName( grad_name.str().c_str() );
  grad->Reset();
  grad->SetTitle( "Gradient; x (m); y(m) ");
  
  // find binning
  const TAxis* xax = in->GetXaxis();
  unsigned nbinsx = xax->GetNbins();
  const TAxis* yax = in->GetYaxis();
  unsigned nbinsy = yax->GetNbins();
  
  // compute the gradient for each bin
  for (unsigned ix=1; ix<nbinsx; ++ix){
    for (unsigned iy=1; iy<nbinsy; ++iy){
      double curval = in->GetBinContent( ix, iy );
      double maxgrad = 0.0;
      for ( int dx = -1; dx <= 1 ; ++dx ){
	for ( int dy = -1; dy <= 1 ; ++dy ){
	  if (dx ==0 && dy == 0 ) continue;
	  double val = in->GetBinContent( ix+dx, iy+dy );
	  double diff = fabs( val - curval );
	  if ( diff > maxgrad ) {
	    maxgrad = diff;
	  }
	}
      }
      grad->SetBinContent( ix, iy, maxgrad );
    }
  }
  
  grad->Write();
  // build graph, ignoring values below cutval
  std::vector<double> xx;
  std::vector<double> yy;
  double maxval = grad->GetMaximum();
  for (unsigned ix=1; ix<nbinsx; ++ix ){
    double x = xax->GetBinCenter( ix );
    for (unsigned iy=1; iy<nbinsy; ++iy ){
      double y = yax->GetBinCenter( iy );
      double curval = grad->GetBinContent( ix, iy );
      if ( curval > cut_below * maxval ){
	xx.push_back( x );
	yy.push_back( y );
	//grad->SetBinContent( ix, iy, 1.0 );
      } else {
	grad->SetBinContent( ix, iy, 0.0 );
      }
    }
  }
  grad->Write();

  // copy grad into vector of xypoints
  std::vector< xypoint > data;
  th2d_to_xypointvector( grad, data );
  CircleHough hough( 50, 0.1, 0.3, 140, 0.0, 0.7, 140, 0.0, 0.7 );
  hough.set_distance_factor( 2 );
  hough.set_minhits( 25 );
  HoughResults res = hough.find_circles( data );

  std::cout<<res<<std::endl;
  hough_display( hough, res );

  if ( res.size() < 1 ) {
    std::cout<<"Failed to find a circle!"<<std::endl;
    std::cout<<"Set Huge Circle"<<std::endl;
    circle.r=1.0;
    circle.xc=0.25;
    circle.yc=0.25;
  } else {
    // take first circle as answer
    circle.r  = res[0].rc;
    circle.xc = res[0].xyc.x;
    circle.yc = res[0].xyc.y;
  }
  
  return circle;
}

void zero_outside_circle( TH2D* in, const Circle_st& c ){
  // find binning
  const TAxis* xax = in->GetXaxis();
  unsigned nbinsx = xax->GetNbins();
  const TAxis* yax = in->GetYaxis();
  unsigned nbinsy = yax->GetNbins();

  // loop over bins
  for (unsigned ix=1; ix<nbinsx+1; ++ix ){
    double x = xax->GetBinCenter( ix );
    for (unsigned iy=1; iy<nbinsy+1; ++iy ){
      double y = yax->GetBinCenter( iy );

      if ( ! c.is_inside( x, y ) ){
	in->SetBinContent( ix, iy, 0.0 );
      }
    }
  }
}

void zero_inside_circle( TH2D* in, const Circle_st& c ){
  // find binning
  const TAxis* xax = in->GetXaxis();
  unsigned nbinsx = xax->GetNbins();
  const TAxis* yax = in->GetYaxis();
  unsigned nbinsy = yax->GetNbins();

  // loop over bins
  for (unsigned ix=1; ix<nbinsx+1; ++ix ){
    double x = xax->GetBinCenter( ix );
    for (unsigned iy=1; iy<nbinsy+1; ++iy ){
      double y = yax->GetBinCenter( iy );

      if ( c.is_inside( x, y ) ){
	in->SetBinContent( ix, iy, 0.0 );
      }
    }
  }
}


// test case code
void find_circle(){

  TFile* fin = new TFile("ptf_charge_analysis_run04511.root");
  //TFile* fin = new TFile("ptf_charge_analysis_run04562.root");
  //TFile* fin = new TFile("ptf_charge_analysis_run04525.root");
  
  TH2D* hq1_2d = (TH2D*)fin->Get("hq1pe");
  
  TH2D* hq1_2d_grad;
  
  Circle_st circ = find_circle_max_grad( hq1_2d, hq1_2d_grad, 0.5 );


  zero_outside_circle( hq1_2d, circ );
  
  TCanvas* tc1 = new TCanvas("tc1","tc1", 600, 600 );
  hq1_2d->Draw("colz");

  TArc *arc1 = new TArc(circ.xc, circ.yc, circ.r );
  arc1->SetLineColor(kRed);
  arc1->SetLineWidth(3);
  arc1->SetFillStyle(0);
  arc1->Draw();


  
  TCanvas* tc2 = new TCanvas("tc2","tc2", 600, 600 );
  hq1_2d_grad->Draw("colz");
  arc1->Draw();


  
}
