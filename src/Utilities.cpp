#include "Utilities.hpp"
#include "PTFStyle.hpp"

#include "TROOT.h"
#include "TStyle.h"
#include "TColor.h"

using namespace std;

const std::vector< double > Utilities::get_bins( std::vector< ScanPoint > scanpoints, char dim ){

  vector< double > positions;

  for(unsigned int iscan=0; iscan<scanpoints.size(); iscan++){
    ScanPoint scanpoint = scanpoints[ iscan ];
    if( scanpoint.x() < 1e-5 && scanpoint.y() < 1e-5 ) continue; // Ignore position (0,0,0)
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
      cout << "Utilities::get_bins Error: input must be x, y or z!" << endl;
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
    cout << "Utilities::get_bins Error: fewer than 3 positions a problem for automatic binning!" << endl;
    exit( EXIT_FAILURE );
  }

  vector< double > bins;
  for(unsigned int i=0; i<positions.size()-1; i++){
    bins.push_back( (positions[i]+positions[i+1])/2. );
  }
  //bins.push_back( 2.*bins[ bins.size()-1 ] - bins[ bins.size()-2 ] );
  //bins.insert( bins.begin(), 2.*bins[0] - bins[1] );
  bins.push_back( 0.75 );
  bins.insert( bins.begin(), 0.0 );

  //std::cout << "bins contains:";
  //vector<double>::iterator it;
  //for (it=bins.begin(); it!=bins.end(); ++it)
  //  std::cout << ' ' << *it;
  //std::cout << '\n';

  return bins;
}

void Utilities::set_style(){ // PTF style
  TString localStyleName = "PTF";
  // -- WhichStyle --
  // 1 = presentation large fonts
  // 2 = presentation small fonts
  // 3 = publication/paper
  Int_t localWhichStyle = 3;
  TStyle* ptfstyle = SetPTFStyle(localWhichStyle, localStyleName);
  gROOT->SetStyle(ptfstyle->GetName());

  // -- Commands to overide the PTF style

  // -- margin --
  gStyle->SetPadTopMargin(0.15);
  gStyle->SetPadBottomMargin(0.15);
  gStyle->SetPadRightMargin(0.15); //more space if no colz option
  gStyle->SetPadLeftMargin(0.15);
  // -- title/lable offset --
  // gStyle->SetTitleOffset(1.5, "x");
  // gStyle->SetTitleOffset(1.7, "y");
  // gStyle->SetLabelOffset(0.02, "x");
  // gStyle->SetLabelOffset(0.02, "y");
  // -- title/label size --
  // gStyle->SetTitleSize(0.04, "x");
  // gStyle->SetTitleSize(0.04, "y");
  // gStyle->SetTitleSize(0.04, "z");
  // gStyle->SetLabelSize(0.034,"x");
  // gStyle->SetLabelSize(0.034,"y");
  // gStyle->SetLabelSize(0.034,"z");
  // -- statistic and title info --
  // gStyle->SetOptTitle(1);
  // gStyle->SetOptStat(1111);
  // -- lines --
  // gStyle->SetLineWidth(4);
  // -- fills/color --
  // gStyle->SetFrameFillColor(0);
  // white color for backgroud
  // gStyle->SetFillColor(1);
  // -- color scheme --
  // - "warm"/red-ish -
  // const Int_t NRGBs = 3;
  // const Int_t NCont = 500;
  // Double_t red[] = {1.00, 1.00, 0.25 };
  // Double_t green[] = {1.00, 0.00, 0.00 };
  // Double_t blue[] = {0.00, 0.00, 0.00 };
  // Double_t stops[] = { 0.25, 0.75, 1.00 };
  // TColor::CreateGradientColorTable(NRGBs, stops, red, green, blue, NCont);
  // gStyle->SetNumberContours(NCont);
  // - inbuilt color schemes -
  // gStyle->SetPalette(1); // use the rainbow color set
  // gStyle->SetPalette(kViridis); // use the viridis color set
  // TColor::InvertPalette(); // invert color palette
  // -- horizontal error bars back --
  // gStyle->SetErrorX(0.5);
  // -- transparent stuff --
  // gStyle->SetFillStyle(0);
  // -- canvas size --
  gStyle->SetCanvasDefW(500);
  gStyle->SetCanvasDefH(500);
  // -- grid lines --
  gStyle->SetPadGridX(0);
  gStyle->SetPadGridY(0);

}

//void Utilities::set_style(){
//  TStyle* ptf_style = new TStyle("ptf_style", "PTF Style");
//  ptf_style->SetPalette(kInvertedDarkBodyRadiator); //change color scheme
//  ptf_style->SetOptStat(0);
//  //ptf_style->SetOptTitle(0);
//  ptf_style->SetOptDate(0);
//  ptf_style->SetLabelSize(0.03, "xyz"); //size of axis value font
//  ptf_style->SetTitleSize(0.035, "xyz"); //size of axis title font
//  ptf_style->SetTitleFont(22, "xyz"); //font option
//  ptf_style->SetLabelFont(22, "xyz");
//  ptf_style->SetTitleOffset(1.2, "y");
//  ptf_style->SetCanvasDefW(500);
//  ptf_style->SetCanvasDefH(500);
//  ptf_style->SetCanvasColor(0);
//  ptf_style->SetCanvasBorderMode(0);
//  ptf_style->SetCanvasBorderSize(0);
//  ptf_style->SetPadTopMargin(0.1);
//  ptf_style->SetPadBottomMargin(0.1);
//  ptf_style->SetPadLeftMargin(0.1);
//  ptf_style->SetPadRightMargin(0.1);
//  ptf_style->SetPadGridX(0);
//  ptf_style->SetPadGridY(0);
//  ptf_style->SetPadTickX(1);
//  ptf_style->SetPadTickY(1);
//  ptf_style->SetFrameBorderMode(0);
//  ptf_style->SetPaperSize(20,24); //US letter size
//
//  gROOT->SetStyle("ptf_style");
//}

bool Utilities::HasWaveform( WaveformFitResult *wf, int pmt ){
  if( pmt == 0 ){
    // check if fit is valid
    if ( wf->fitstat != 0 ) return false;
    // check if ringing is bigger than waveform
    if ( wf->amp - wf->sinamp < 20.0 ) return false; 
    // check if pulse width is reasonable
    if ( wf->sigma < 1.0 || wf->sigma > 7.0 ) return false;
    // check if mean pulse time is reasonable
    if ( wf->mean < 28.0 || wf->mean > 45.0 ) return false;
    // cut on chi2
    //if ( wf->chi2 > 200.0 ) return false;
  }
  if( pmt == 1 ){
    // Checks if using simpler analysis of bin furthest from pedestal
    if ( wf->amp < 25.0 ) return false; 
  }
  return true;
}
