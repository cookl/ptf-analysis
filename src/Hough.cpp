#include "Hough.hpp"

#include <sstream>

#include "TH2D.h"
#include "TDirectory.h"
#include "TEllipse.h"
#include "TMarker.h"

const double pi = std::acos(-1);

unsigned num_circles( const HoughResults& hrs) {
  unsigned nc=0;
  for ( const HoughResult& hr : hrs ){
    if ( hr.type == HoughCircle ) ++nc;
  }
  return nc;
}

std::ostream& operator<<( std::ostream& os, const HoughResults& hrs ){
  for ( const HoughResult& hr : hrs ){
    os << hr;
  }
  return os;
}

std::ostream& operator<<( std::ostream& os, const HoughResult& hr ){
  if ( hr.type == HoughCircle ){
    os<<"Hough Circle: (xc,yc)= ( "<<hr.xyc.x<<", "<<hr.xyc.y
      <<" )  radius="<<hr.rc<<" hough-peak="<<hr.peakval<<std::endl;
  } else {
    os<<"Hough Unused Hits::"<<std::endl;
  }
  os << "  Nhits="<<hr.data.size()<<std::endl;
  for ( const xypoint& xy : hr.data ){
    os << "    (x,y)= ( "<<xy.x<<", "<<xy.y<<" )"<<std::endl;
  }
  return os;
}



CircleHough::CircleHough(  unsigned nbins_radius, double rmin, double  rmax,
			   unsigned nbins_x, double xmin, double xmax,
			   unsigned nbins_y, double ymin, double ymax ){
  static unsigned instance_count=0;
  ++instance_count;
  std::ostringstream os;

  TDirectory * curdir = gDirectory;
  houghdir = gDirectory->mkdir( (std::string("circlehough_")+std::to_string(instance_count)).c_str() );
  houghdir->cd();

  double dr = (rmax-rmin)/nbins_radius;
  double rbmin = rmin;
  double rbmax = rmin+dr;
  for (int i = 0; i< nbins_radius; ++i){
    std::string hname = std::string( "ht_") + std::to_string(instance_count) 
      + "_rbin_"+ std::to_string(i);
    std::string htitle = std::string( "R = ") + std::to_string( (rbmin+rbmax)/2 ) + " ; xc; yc ";
    fRbins.push_back( std::pair<double,double>( rbmin, rbmax ) );
    fTransformed.push_back( new TH2D( hname.c_str(), htitle.c_str(),
				      nbins_x, xmin, xmax,
				      nbins_y, ymin, ymax )  );
    rbmin+=dr;
    rbmax+=dr;
  }
  curdir->cd();
}

CircleHough::~CircleHough(){
   for ( unsigned i=0; i<fTransformed.size(); ++i ){
    fTransformed[i]->SetDirectory(0);
    delete fTransformed[i];
  }
}


const HoughResults& CircleHough::find_circles( const std::vector< xypoint >& data ){
  ++nfind_circles;
  fresults.clear();
  std::vector< xypoint > unused_hits = data;
  bool done = false;
  while ( unused_hits.size() > minhits && !done ){
    hough_transform( unused_hits );
    HoughResult hr = find_maximum( unused_hits );
    if ( hr.peakval > threshold && hr.data.size() > minhits ) {
      hr.type = HoughCircle;
    } else {
      hr.data.insert( hr.data.end(), unused_hits.begin(), unused_hits.end() );
      hr.type = HoughUnusedPoints;
      done = true;
    }

    fresults.push_back( hr );
    save_hough_histo( fresults.size(), fTransformed[ hr.rbin ] );
    //plot_candidate( fresults.size(), hr );
  }

  // cleanup
  for ( unsigned i=0; i<fTransformed.size(); ++i ){
    fTransformed[i]->SetDirectory(0);
    //delete fTransformed[i];
  }

  return fresults;
}

void CircleHough::hough_transform( const std::vector< xypoint >& data ){
  for ( TH2D* h : fTransformed) h->Reset();
  TH2D*h = fTransformed[0];
  double xbwid = h->GetXaxis()->GetBinWidth(1);
  double ybwid = h->GetYaxis()->GetBinWidth(1);
  for ( const xypoint& xy : data ){
    for ( unsigned ibin=0; ibin < fRbins.size(); ++ibin ){ 
      std::pair< double, double > rbin = fRbins[ibin];
      double rc = (rbin.first+rbin.second)/2;
      // pick number of angles based on radius
      unsigned nang = unsigned( 2 * rc / xbwid );
      double   dtheta = 2*pi/nang;
      for ( unsigned itheta = 0; itheta<nang; ++itheta ){
	double a = xy.x - rc*std::cos( itheta * dtheta );
	double b = xy.y - rc*std::sin( itheta * dtheta );
	fTransformed[ ibin ]->Fill( a, b );
	fTransformed[ ibin ]->Fill( a+xbwid, b, 0.5 );
	fTransformed[ ibin ]->Fill( a-xbwid, b, 0.5 );
	fTransformed[ ibin ]->Fill( a, b+ybwid, 0.5 );
	fTransformed[ ibin ]->Fill( a, b-ybwid, 0.5 );
	fTransformed[ ibin ]->Fill( a+xbwid, b+ybwid, 0.25 );
	fTransformed[ ibin ]->Fill( a-xbwid, b-ybwid, 0.25 );
	fTransformed[ ibin ]->Fill( a-xbwid, b+ybwid, 0.25 );
	fTransformed[ ibin ]->Fill( a+xbwid, b-ybwid, 0.25 );
      }
    }
  }
}


HoughResult CircleHough::find_maximum( std::vector< xypoint >& hits ){
  HoughResult curbest( 0., xypoint(0., 0.), 0 );
  // loop over the hough transform histograms to find the peak
  // and store the "best" circle center and radius
  for ( unsigned ibin=0; ibin < fTransformed.size(); ++ibin ){
    std::pair< double, double > rbin = fRbins[ibin];
    double rc = (rbin.first+rbin.second)/2;
    TH2D* h = fTransformed[ibin];
    for ( int ix=1; ix<=h->GetNbinsX(); ++ix ){
      float x = h->GetXaxis()->GetBinCenter(ix);
      for ( int iy=1; iy<=h->GetNbinsY(); ++iy ){
	if ( h->GetBinContent( ix, iy ) > curbest.peakval ) {
	  float y = h->GetYaxis()->GetBinCenter(iy);
	  curbest.rc = rc;
	  curbest.xyc = xypoint( x, y );
	  curbest.peakval = h->GetBinContent( ix, iy );
	  curbest.rbin = ibin;
	}
      }
    }
  }
  // find the hits that are associated with the circle and add them
  // to the result, otherwise add them to list of unused_hits
  // use bin sizes in x,y as threshold distance for hit to be from circle
  TH2D*h = fTransformed[0];
  double xbwid = h->GetXaxis()->GetBinWidth(1);
  double ybwid = h->GetYaxis()->GetBinWidth(1);
  double rthres = drscaling * std::sqrt( xbwid*xbwid + ybwid*ybwid );
  std::vector< xypoint > unused_hits;
  for ( xypoint xy : hits ){
    double dr = std::sqrt( (xy.x-curbest.xyc.x)*(xy.x-curbest.xyc.x) + (xy.y-curbest.xyc.y)*(xy.y-curbest.xyc.y) );
    if ( fabs( dr - curbest.rc ) < rthres ){
      curbest.data.push_back( xy );
    } else {
      unused_hits.push_back( xy );
    }
  }

  hits = unused_hits;
  return curbest;
}

void CircleHough::save_hough_histo( unsigned num, TH2D* histo ){
  if ( nfind_circles > 10 ) return;
  TDirectory* curdir = gDirectory;
  houghdir->cd();

  std::string hname = std::string( "ev_" ) + std::to_string(nfind_circles) + "_cand_" + std::to_string( num ) + histo->GetName() ;
  TH2D* savehist = (TH2D*)histo->Clone( hname.c_str() );
  savehist->SetName( hname.c_str() );
  savehist->SetDirectory( houghdir );
  curdir->cd();
}

void CircleHough::plot_candidate( unsigned num, const HoughResult & hr ){
  TDirectory* curdir = gDirectory;
  houghdir->cd();
  
  std::string hname = std::string("hcircle_")+std::to_string(num);
  TH2D* hplot = (TH2D*)fTransformed[ hr.rbin ]->Clone( hname.c_str() );
  hplot->SetName( hname.c_str() );
  hplot->SetDirectory( houghdir );
  hplot->Reset();

  for ( const xypoint& xy : hr.data ){
    hplot->Fill( xy.x, xy.y );
  }

  TEllipse *el = new TEllipse( hr.xyc.x, hr.xyc.y, hr.rc );
  el->SetFillStyle(0);
  hplot->GetListOfFunctions()->Add( el );
  TMarker *mark = new TMarker( hr.xyc.x, hr.xyc.y, 20 );
  hplot->GetListOfFunctions()->Add( mark );

  curdir->cd();
}
  
