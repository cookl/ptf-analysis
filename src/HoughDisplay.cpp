#include "HoughDisplay.hpp"

#include <TH2D.h>
#include <TDirectory.h>
#include <TGraph.h>
#include <TCanvas.h>
#include <TBox.h>
#include <TEllipse.h>

// new event_display
void hough_display( CircleHough & ch,
		    const HoughResults & hcr ){
		    
  const std::vector<int> colors = { kRed   , kBlue, kGreen,   kOrange,
				    kViolet, kCyan, kMagenta, kPink+6 };   

  // build histogram range from transformed histogram
  std::vector< TH2D* > vcht = ch.get_transform();
  if ( vcht.size() < 1 ) return;
  unsigned nbins_x = vcht[0]->GetNbinsX();
  unsigned nbins_y = vcht[0]->GetNbinsX();
  double   xmin    = vcht[0]->GetXaxis()->GetXmin();
  double   xmax    = vcht[0]->GetXaxis()->GetXmax();
  double   ymin    = vcht[0]->GetYaxis()->GetXmin();
  double   ymax    = vcht[0]->GetYaxis()->GetXmax();
  // build histogram name
  std::string hname = std::string{"hough_display"};
  TH2D* hist = new TH2D( hname.c_str(), " ; X (cm); Y (cm) ",
			 nbins_x, xmin, xmax,
			 nbins_y, ymin, ymax );

  // put the histogram on a canvas
  std::string canname = std::string{"tc_"} + hname;
  TCanvas * tc = new TCanvas( canname.c_str(), canname.c_str() );
  tc->cd();
  hist->Draw();

  // Now we will add each of the datapoints color coded for each ellipse
  float xbwid = (xmax-xmin)/nbins_x;
  float ybwid = (ymax-ymin)/nbins_y;
  for ( unsigned iel = 0; iel < hcr.size(); ++iel ){
    int curcol = 0;
    if ( iel<colors.size() ) curcol = colors[iel];
    if ( hcr[iel].type == HoughUnusedPoints ) curcol=kGray;
    
    // Add datapoints as boxes to the plot
    for ( const xypoint & pt : hcr[iel].data ){
      TBox *b =new TBox( pt.x-xbwid, pt.y-ybwid, pt.x+xbwid, pt.y+ybwid );
      b->SetFillStyle(1001);
      b->SetFillColor( curcol );
      b->Draw();
      b->SetBit( kCanDelete );
      b->SetBit( kMustCleanup );
    }
    if ( hcr[iel].type == HoughUnusedPoints ) continue;
    
    // Draw circle from hough transform
    TEllipse *tel= new TEllipse( hcr[iel].xyc.x, hcr[iel].xyc.y, hcr[iel].rc);
    tel->SetFillStyle(0);
    tel->SetLineColor( curcol );
    tel->Draw();
    hist->GetListOfFunctions()->Add( tel );
    

  }

  //  tc->Write();
  //curdir->cd();
}
