/// collect_timing_plots
/// John Walker (Apr. 2020)
///
/// Make uniform set of plots from each PTF run
/// No inputs.  All files hard coded.
/// Input files are from output of ptf_timing_analysis.app

/// make_1d_from_2d
/// Take z-axis values of TH2D and use them to fill a newly constructed TH1D
/// Last four arguments of function are used to build the TH1D
/// Ignores values that are less than 0.001
/// Returns a pointer to the newly created TH1D

#include "TFile.h"
#include "TCanvas.h"
#include "TH2D.h"
#include "TH1D.h"
#include "TROOT.h"
#include "TStyle.h"
#include "TColor.h"
#include "TPad.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TMultiGraph.h"
#include "TVector2.h"
#include "TAxis.h"

#include <stdlib.h>
#include <string>
#include <iostream>
#include <sstream>

using namespace std;

// --- PTF style ---
TStyle* SetPTFStyle(Int_t WhichStyle = 1, TString styleName = "PTF") {
  TStyle *ptfStyle= new TStyle(styleName, "PTF approved plots style");
  // -- WhichStyle --
  // 1 = presentation large fonts
  // 2 = presentation small fonts
  // 3 = publication/paper
  Int_t FontStyle = 22;
  Float_t FontSizeLabel = 0.035;
  Float_t FontSizeTitle = 0.05;
  Float_t YOffsetTitle = 1.3;
  switch(WhichStyle) {
    case 1:
      FontStyle = 42;
      FontSizeLabel = 0.05;
      FontSizeTitle = 0.065;
      YOffsetTitle = 1.19;
      break;
    case 2:
      FontStyle = 42;
      FontSizeLabel = 0.035;
      FontSizeTitle = 0.05;
      YOffsetTitle = 1.6;
      break;
    case 3:
      FontStyle = 132;
      FontSizeLabel = 0.03;
      FontSizeTitle = 0.035;
      YOffsetTitle = 1.3;
      break;
  }
  // use plain black on white colors
  ptfStyle->SetFrameBorderMode(0);
  ptfStyle->SetCanvasBorderMode(0);
  ptfStyle->SetCanvasBorderSize(0);
  ptfStyle->SetPadBorderMode(0);
  ptfStyle->SetPadColor(0);
  ptfStyle->SetCanvasColor(0);
  ptfStyle->SetStatColor(0);
  ptfStyle->SetFillColor(0);
  ptfStyle->SetEndErrorSize(4);
  ptfStyle->SetStripDecimals(kFALSE);
  ptfStyle->SetLegendBorderSize(0);
  ptfStyle->SetLegendFont(FontStyle);
  // set the paper & margin sizes
  ptfStyle->SetPaperSize(20, 26);
  ptfStyle->SetPadTopMargin(0.1);
  ptfStyle->SetPadBottomMargin(0.15);
  ptfStyle->SetPadRightMargin(0.13);
  // 0.075 -> 0.13 for colz option
  ptfStyle->SetPadLeftMargin(0.16);
  //to include both large/small font options
  // Fonts, sizes, offsets
  ptfStyle->SetTextFont(FontStyle);
  ptfStyle->SetTextSize(0.08);
  ptfStyle->SetLabelFont(FontStyle, "x");
  ptfStyle->SetLabelFont(FontStyle, "y");
  ptfStyle->SetLabelFont(FontStyle, "z");
  ptfStyle->SetLabelFont(FontStyle, "t");
  ptfStyle->SetLabelSize(FontSizeLabel, "x");
  ptfStyle->SetLabelSize(FontSizeLabel, "y");
  ptfStyle->SetLabelSize(FontSizeLabel, "z");
  ptfStyle->SetLabelOffset(0.015, "x");
  ptfStyle->SetLabelOffset(0.015, "y");
  ptfStyle->SetLabelOffset(0.015, "z");
  ptfStyle->SetTitleFont(FontStyle, "x");
  ptfStyle->SetTitleFont(FontStyle, "y");
  ptfStyle->SetTitleFont(FontStyle, "z");
  ptfStyle->SetTitleFont(FontStyle, "t");
  ptfStyle->SetTitleSize(FontSizeTitle, "y");
  ptfStyle->SetTitleSize(FontSizeTitle, "x");
  ptfStyle->SetTitleSize(FontSizeTitle, "z");
  ptfStyle->SetTitleSize(0.05, "t");
  ptfStyle->SetTitleOffset(1.14, "x");
  ptfStyle->SetTitleOffset(YOffsetTitle, "y");
  ptfStyle->SetTitleOffset(1.2, "z");
  ptfStyle->SetTitleStyle(0);
  //ptfStyle->SetTitleFontSize(0.05);
  //0.08
  ptfStyle->SetTitleFont(FontStyle, "pad");
  ptfStyle->SetTitleBorderSize(0);
  ptfStyle->SetTitleX(0.1f);
  ptfStyle->SetTitleW(0.8f);
  ptfStyle->SetTitleY(0.93f);
  // use bold lines and markers
  ptfStyle->SetMarkerStyle(20);
  ptfStyle->SetHistLineWidth( Width_t(2.5) );
  ptfStyle->SetLineStyleString(2, "[12 12]");
  // postscript dashes
  // get rid of X error bars and y error bar caps
  ptfStyle->SetErrorX(0.001);
  // do not display any of the standard histogram decorations
  //ptfStyle->SetOptTitle(0);
  ptfStyle->SetOptStat(0);
  ptfStyle->SetOptFit(0);
  // put tick marks on top and RHS of plots
  ptfStyle->SetPadTickX(1);
  ptfStyle->SetPadTickY(1);
  // -- color --
  ptfStyle->SetFillColor(1);
  // make color fillings (not white)
  // - color setup for 2D -
  // - "cold"/ blue-ish -
  Double_t red[] = { 0.00, 0.00, 0.00 };
  Double_t green[] = { 1.00, 0.00, 0.00 };
  Double_t blue[] = { 1.00, 1.00, 0.25 };
  // - "warm" red-ish colors -
  // Double_t red[] = {1.00, 1.00, 0.25 };
  // Double_t green[] = {1.00, 0.00, 0.00 };
  // Double_t blue[] = {0.00, 0.00, 0.00 };
  Double_t stops[] = { 0.25, 0.75, 1.00 };
  const Int_t NRGBs = 3;
  const Int_t NCont = 500;
  TColor::CreateGradientColorTable(NRGBs, stops, red, green, blue, NCont);
  ptfStyle->SetNumberContours(NCont);
  // - inbuilt color schemes -
  // ptfStyle->SetPalette(1); // use the rainbow color set
  //ptfStyle->SetPalette(kViridis); // use the viridis color set
  //ptfStyle->SetPalette(kSunset); // use the sunset color set
  //ptfStyle->SetPalette(kColorPrintableOnGrey); // use the colorPrintableOnGrey color set
  //ptfStyle->SetPalette(kCubehelix); // use the cubehelix color set
  ptfStyle->SetPalette(kLightTemperature);
  //TColor::InvertPalette(); // invert color palette

  return ptfStyle;
}

TH1D* make_1d_from_2d( TH2D* hin, std::string title, int nbins, double xmin, double xmax ){
  static int count=0;
  ++count;
  std::ostringstream os;
  os<<"h1d_from_2d_"<<count;
  TH1D* h = new TH1D(os.str().c_str(), title.c_str(), nbins, xmin, xmax );

  for (unsigned ix = 1; ix <= hin->GetNbinsX(); ++ix ){
    for (unsigned iy = 1; iy <= hin->GetNbinsY(); ++iy ){
      if ( fabs( hin->GetBinContent( ix, iy) ) > 0.02 ){
	    h->Fill( hin->GetBinContent( ix, iy ) );
      }
    }
  } 
  return h;
}


void collect_timing_plots(){

  TString localStyleName = "PTF";
  // -- WhichStyle --
  // 1 = presentation large fonts
  // 2 = presentation small fonts
  // 3 = publication/paper
  Int_t localWhichStyle = 3;
  TStyle* ptfstyle = SetPTFStyle(localWhichStyle, localStyleName);
  gROOT->SetStyle(ptfstyle->GetName());

  // -- margin --
  gStyle->SetPadTopMargin(0.15);
  gStyle->SetPadBottomMargin(0.15);
  gStyle->SetPadRightMargin(0.15); //more space if no colz option
  gStyle->SetPadLeftMargin(0.15);
  // -- canvas size --
  gStyle->SetCanvasDefW(500);
  gStyle->SetCanvasDefH(500);
  // -- grid lines --
  gStyle->SetPadGridX(1);
  gStyle->SetPadGridY(1);
  
  std::vector< std::string > bscan2_tags = {
    "Run 4427 | 0mG | Air | Cover off              ",
    "Run 4431 | 0mG | Air | Cover off              ",
    "Run 4446 | X=-100mG | Air | Cover off         ",
    "Run 4447 | X=-50mG | Air | Cover off          ",
    "Run 4434 | X=+50mG | Air | Cover off          ",
    "Run 4433 | X=+100mG | Air | Cover off         ",
    "Run 4458 | Y=-100mG | Air | Cover off         ",
    "Run 4459 | Y=-50mG | Air | Cover off          ",
    "Run 4457 | Y=+50mG | Air | Cover off          ",
    "Run 4448 | Y=+100mG | Air | Cover off         ",
    "Run 4444 | Z=-100mG | Air | Cover off         ",
    "Run 4445 | Z=-50mG | Air | Cover off          ",
    "Run 4443 | Z=+50mG | Air | Cover off          ",
    "Run 4442 | Z=+100mG | Air | Cover off         ",

    "Run 4525 | 0mG | Water | Cover off            ",
    "Run 4526 | 0mG | Water | Cover off | Rotated  ",
    "Run 4533 | X=-100mG | Water | Cover off       ",
    "Run 4534 | X=-50mG | Water | Cover off        ",
    "Run 4531 | X=+50mG | Water | Cover off        ",
    "Run 4530 | X=+100mG | Water | Cover off       ",
    "Run 4538 | Y=-100mG | Water | Cover off       ",
    "Run 4540 | Y=-50mG | Water | Cover off        ",
    "Run 4537 | Y=+50mG | Water | Cover off        ",
    "Run 4536 | Y=+100mG | Water | Cover off       ",

    "Run 4562 | No compensation | Water | Cover off",
    "Run 4551 | No compensation | Water | Cover on ",
    "Run 4584 | No compensation | Water | Cover off",
    "Run 4581 | No compensation | Water | Cover on "
  };
   
  std::vector< TFile* > bscan2_files = {
    new TFile( "ptf_timing_analysis_run04427.root" ),
    new TFile( "ptf_timing_analysis_run04431.root" ),
    new TFile( "ptf_timing_analysis_run04446.root" ),
    new TFile( "ptf_timing_analysis_run04447.root" ),
    new TFile( "ptf_timing_analysis_run04434.root" ),
    new TFile( "ptf_timing_analysis_run04433.root" ),
    new TFile( "ptf_timing_analysis_run04458.root" ),
    new TFile( "ptf_timing_analysis_run04459.root" ),
    new TFile( "ptf_timing_analysis_run04457.root" ),
    new TFile( "ptf_timing_analysis_run04448.root" ),
    new TFile( "ptf_timing_analysis_run04444.root" ),
    new TFile( "ptf_timing_analysis_run04445.root" ),
    new TFile( "ptf_timing_analysis_run04443.root" ),
    new TFile( "ptf_timing_analysis_run04442.root" ),

    new TFile( "ptf_timing_analysis_run04525.root" ),
    new TFile( "ptf_timing_analysis_run04526.root" ),
    new TFile( "ptf_timing_analysis_run04533.root" ),
    new TFile( "ptf_timing_analysis_run04534.root" ),
    new TFile( "ptf_timing_analysis_run04531.root" ), 
    new TFile( "ptf_timing_analysis_run04530.root" ),
    new TFile( "ptf_timing_analysis_run04538.root" ),
    new TFile( "ptf_timing_analysis_run04540.root" ),
    new TFile( "ptf_timing_analysis_run04537.root" ),
    new TFile( "ptf_timing_analysis_run04536.root" ),

    new TFile( "ptf_timing_analysis_run04562.root" ),
    new TFile( "ptf_timing_analysis_run04551.root" ),
    new TFile( "ptf_timing_analysis_run04584.root" ),
    new TFile( "ptf_timing_analysis_run04581.root" )
  };


  //KEY: TH2D h_rtt; Transit time
  //KEY: TH2D  h_tts; Transit time spread
  std::vector< TH1D* > vec_h_rtt;
  std::vector< TH1D* > vec_h_tts;
  for ( unsigned i=0 ; i < bscan2_files.size() ; ++i ){
    TH2D* h_rtt_2d = (TH2D*)bscan2_files[i]->Get("h_rtt");
    TH2D* h_tts_2d = (TH2D*)bscan2_files[i]->Get("h_tts");

    TH1D* h_rtt = make_1d_from_2d( h_rtt_2d, " ; RTT [ns] ; scan points / bin", 100, 28., 40. );
    double scale = h_rtt->GetMaximum();
    h_rtt->Scale( 1.0/scale );
    TH1D* h_tts = make_1d_from_2d( h_tts_2d, " ; TTS [ns] ; scan points / bin", 100, 1.0, 4.5 );
    scale = h_tts->GetMaximum();
    h_tts->Scale( 1.0/scale );

    vec_h_rtt.push_back( h_rtt );
    vec_h_tts.push_back( h_tts );

    ostringstream os;
    os << "tc_" << i;
    TCanvas * tc = new TCanvas( os.str().c_str(), bscan2_tags[i].c_str() );
    h_rtt->Draw("e1");
    string fname(bscan2_files[i]->GetName());
    fname.erase(fname.length()-5);
    ostringstream os2;
    os2 << fname << "_rtt_1d.pdf";
    tc->Print( os2.str().c_str() );

    h_tts->Draw("e1");
    gPad->Modified();
    gPad->Update();
    ostringstream os3;
    os3 << fname << "_tts_1d.pdf";
    tc->Print( os3.str().c_str() );
  }


  // print table of results from the TH1D's
  cout<< "\\begin{center}" << endl;
  cout<< "\\begin{tabular}{ l | c c }" << endl;
  cout << "Run and Description\t& RTT RMS\t& TTS mean \\\\ \\hline" <<endl;
  for ( unsigned i=0 ; i < bscan2_files.size() ; ++i ){

    cout<< bscan2_tags[i] <<"\t& ";
    cout<< fixed << setprecision(4) << setw(12 )<< vec_h_rtt[i]->GetRMS() << " $\\pm$ " << vec_h_rtt[i]->GetRMSError() <<"\t& ";
    cout<< fixed << setprecision(4) << setw(12 )<< vec_h_tts[i]->GetMean() << " $\\pm$ " << vec_h_tts[i]->GetMeanError() ;
    cout<< "\\\\ "<<endl;
  }
  cout<< "\\end{tabular}" <<endl;
  cout<< "\\end{center}" <<endl;


  // Plot RTT spread vs magnetic field
  // x-direction, air, cover off
  TCanvas * tc = new TCanvas( "timing_vs_bfield", "timing_vs_bfield" );
  //tc->SetGrid();
  const int n = 6;
  double x_rxa0[n] = { -100., -50., 0., 0., 50., 100. };
  double y_rxa0[n] = { vec_h_rtt[2]->GetRMS(),
                  vec_h_rtt[3]->GetRMS(),
                  vec_h_rtt[0]->GetRMS(),
                  vec_h_rtt[1]->GetRMS(),
                  vec_h_rtt[4]->GetRMS(),
                  vec_h_rtt[5]->GetRMS() };
  double ex_rxa0[n] = { 10., 10., 10., 10., 10., 10. };
  double ey_rxa0[n] = { vec_h_rtt[2]->GetRMSError(),
                   vec_h_rtt[3]->GetRMSError(),
                   vec_h_rtt[0]->GetRMSError(),
                   vec_h_rtt[4]->GetRMSError(),
                   vec_h_rtt[5]->GetRMSError() };
  TGraphErrors* gr_rxa0 = new TGraphErrors(n,x_rxa0,y_rxa0,ex_rxa0,ey_rxa0);
  gr_rxa0->SetTitle( "Air | Cover off; Bx [mG]; RTT RMS [ns]" );
  gr_rxa0->SetMarkerColor(kMagenta+3);
  gr_rxa0->SetMarkerStyle(21);
  gr_rxa0->Draw("AP");
  string plotname = string("ptf_timing_analysis_rtt_vs_bfield_air_coveroff_x.pdf");
  tc->SaveAs(plotname.c_str(),"pdf");

  // y-direction, air, cover off
  double x_rya0[n] = { -100., -50., 0., 0., 50., 100. };
  double y_rya0[n] = { vec_h_rtt[6]->GetRMS(),
                  vec_h_rtt[7]->GetRMS(),
                  vec_h_rtt[0]->GetRMS(),
                  vec_h_rtt[1]->GetRMS(),
                  vec_h_rtt[8]->GetRMS(),
                  vec_h_rtt[9]->GetRMS() };
  double ex_rya0[n] = { 10., 10., 10., 10., 10., 10. };
  double ey_rya0[n] = { vec_h_rtt[6]->GetRMSError(),
                   vec_h_rtt[7]->GetRMSError(),
                   vec_h_rtt[0]->GetRMSError(),
                   vec_h_rtt[1]->GetRMSError(),
                   vec_h_rtt[8]->GetRMSError(),
                   vec_h_rtt[9]->GetRMSError() };
  TGraphErrors* gr_rya0 = new TGraphErrors(n,x_rya0,y_rya0,ex_rya0,ey_rya0);
  gr_rya0->SetTitle( "Air | Cover off; By [mG]; RTT RMS [ns]" );
  gr_rya0->SetMarkerColor(kMagenta+3);
  gr_rya0->SetMarkerStyle(21);
  gr_rya0->Draw("AP");
  gPad->Modified();
  gPad->Update();
  plotname = string("ptf_timing_analysis_rtt_vs_bfield_air_coveroff_y.pdf");
  tc->SaveAs(plotname.c_str(),"pdf");

  // z-direction, air, cover off
  double x_rza0[n] = { -100., -50., 0., 0., 50., 100. };
  double y_rza0[n] = { vec_h_rtt[10]->GetRMS(),
                  vec_h_rtt[11]->GetRMS(),
                  vec_h_rtt[0]->GetRMS(),
                  vec_h_rtt[1]->GetRMS(),
                  vec_h_rtt[12]->GetRMS(),
                  vec_h_rtt[13]->GetRMS() };
  double ex_rza0[n] = { 40., 40., 40., 40., 40., 40. };
  double ey_rza0[n] = { vec_h_rtt[10]->GetRMSError(),
                   vec_h_rtt[11]->GetRMSError(),
                   vec_h_rtt[0]->GetRMSError(),
                   vec_h_rtt[1]->GetRMSError(),
                   vec_h_rtt[12]->GetRMSError(),
                   vec_h_rtt[13]->GetRMSError() };
  TGraphErrors* gr_rza0 = new TGraphErrors(n,x_rza0,y_rza0,ex_rza0,ey_rza0);
  gr_rza0->SetTitle( "Air | Cover off; Bz [mG]; RTT RMS [ns]" );
  gr_rza0->SetMarkerColor(kMagenta+3);
  gr_rza0->SetMarkerStyle(21);
  gr_rza0->Draw("AP");
  gPad->Modified();
  gPad->Update();
  plotname = string("ptf_timing_analysis_rtt_vs_bfield_air_coveroff_z.pdf");
  tc->SaveAs(plotname.c_str(),"pdf");

  // x-direction, water, cover off
  const int nwx = 5;
  double x_rxw0[nwx] = { -100., -50., 0., 0., 50. };
  double y_rxw0[nwx] = { vec_h_rtt[16]->GetRMS(),
                  vec_h_rtt[17]->GetRMS(),
                  vec_h_rtt[14]->GetRMS(),
                  vec_h_rtt[15]->GetRMS(),
                  vec_h_rtt[18]->GetRMS() };
  double ex_rxw0[nwx] = { 10., 10., 10., 10., 10. };
  double ey_rxw0[nwx] = { vec_h_rtt[16]->GetRMSError(),
                   vec_h_rtt[17]->GetRMSError(),
                   vec_h_rtt[14]->GetRMSError(),
                   vec_h_rtt[15]->GetRMSError(),
                   vec_h_rtt[18]->GetRMSError() };
  TGraphErrors* gr_rxw0 = new TGraphErrors(nwx,x_rxw0,y_rxw0,ex_rxw0,ey_rxw0);
  gr_rxw0->SetTitle( "Water | Cover off; Bx [mG]; RTT RMS [ns]" );
  gr_rxw0->SetMarkerColor(kOrange-3);
  gr_rxw0->SetMarkerStyle(21);
  gr_rxw0->Draw("AP");
  gPad->Modified();
  gPad->Update();
  plotname = string("ptf_timing_analysis_rtt_vs_bfield_water_coveroff_x.pdf");
  tc->SaveAs(plotname.c_str(),"pdf");

  // y-direction, water, cover off
  double x_ryw0[n] = { -100., -50., 0., 0., 50., 100. };
  double y_ryw0[n] = { vec_h_rtt[20]->GetRMS(),
                  vec_h_rtt[21]->GetRMS(),
                  vec_h_rtt[14]->GetRMS(),
                  vec_h_rtt[15]->GetRMS(),
                  vec_h_rtt[22]->GetRMS(),
                  vec_h_rtt[23]->GetRMS() };
  double ex_ryw0[n] = { 10., 10., 10., 10., 10., 10. };
  double ey_ryw0[n] = { vec_h_rtt[20]->GetRMSError(),
                   vec_h_rtt[21]->GetRMSError(),
                   vec_h_rtt[14]->GetRMSError(),
                   vec_h_rtt[15]->GetRMSError(),
                   vec_h_rtt[22]->GetRMSError(),
                   vec_h_rtt[23]->GetRMSError() };
  TGraphErrors* gr_ryw0 = new TGraphErrors(n,x_ryw0,y_ryw0,ex_ryw0,ey_ryw0);
  gr_ryw0->SetTitle( "Water | Cover off; By [mG]; RTT RMS [ns]" );
  gr_ryw0->SetMarkerColor(kOrange-3);
  gr_ryw0->SetMarkerStyle(21);
  gr_ryw0->Draw("AP");
  gPad->Modified();
  gPad->Update();
  plotname = string("ptf_timing_analysis_rtt_vs_bfield_water_coveroff_y.pdf");
  tc->SaveAs(plotname.c_str(),"pdf");


  // Plot average TTS vs magnetic field
  // x-direction, air, cover off
  double x_sxa0[n] = { -100., -50., 0., 0., 50., 100. };
  double y_sxa0[n] = { vec_h_tts[2]->GetMean(),
                  vec_h_tts[3]->GetMean(),
                  vec_h_tts[0]->GetMean(),
                  vec_h_tts[1]->GetMean(),
                  vec_h_tts[4]->GetMean(),
                  vec_h_tts[5]->GetMean() };
  double ex_sxa0[n] = { 10., 10., 10., 10., 10., 10. };
  double ey_sxa0[n] = { vec_h_tts[2]->GetMeanError(),
                   vec_h_tts[3]->GetMeanError(),
                   vec_h_tts[0]->GetMeanError(),
                   vec_h_tts[1]->GetMeanError(),
                   vec_h_tts[4]->GetMeanError(),
                   vec_h_tts[5]->GetMeanError() };
  TGraphErrors* gr_sxa0 = new TGraphErrors(n,x_sxa0,y_sxa0,ex_sxa0,ey_sxa0);
  gr_sxa0->SetTitle( "Air | Cover off; Bx [mG]; Mean TTS [ns]" );
  gr_sxa0->SetMarkerColor(kMagenta+3);
  gr_sxa0->SetMarkerStyle(21);
  gr_sxa0->Draw("AP");
  gPad->Modified();
  gPad->Update();
  plotname = string("ptf_timing_analysis_tts_vs_bfield_air_coveroff_x.pdf");
  tc->SaveAs(plotname.c_str(),"pdf");

  // y-direction, air, cover off
  double x_sya0[n] = { -100., -50., 0., 0., 50., 100. };
  double y_sya0[n] = { vec_h_tts[6]->GetMean(),
                  vec_h_tts[7]->GetMean(),
                  vec_h_tts[0]->GetMean(),
                  vec_h_tts[1]->GetMean(),
                  vec_h_tts[8]->GetMean(),
                  vec_h_tts[9]->GetMean() };
  double ex_sya0[n] = { 10., 10., 10., 10., 10., 10. };
  double ey_sya0[n] = { vec_h_tts[6]->GetMeanError(),
                   vec_h_tts[7]->GetMeanError(),
                   vec_h_tts[0]->GetMeanError(),
                   vec_h_tts[1]->GetMeanError(),
                   vec_h_tts[8]->GetMeanError(),
                   vec_h_tts[9]->GetMeanError() };
  TGraphErrors* gr_sya0 = new TGraphErrors(n,x_sya0,y_sya0,ex_sya0,ey_sya0);
  gr_sya0->SetTitle( "Air | Cover off; By [mG]; Mean TTS [ns]" );
  gr_sya0->SetMarkerColor(kMagenta+3);
  gr_sya0->SetMarkerStyle(21);
  gr_sya0->Draw("AP");
  gPad->Modified();
  gPad->Update();
  plotname = string("ptf_timing_analysis_tts_vs_bfield_air_coveroff_y.pdf");
  tc->SaveAs(plotname.c_str(),"pdf");

  // z-direction, air, cover off
  double x_sza0[n] = { -100., -50., 0., 0., 50., 100. };
  double y_sza0[n] = { vec_h_tts[10]->GetMean(),
                  vec_h_tts[11]->GetMean(),
                  vec_h_tts[0]->GetMean(),
                  vec_h_tts[1]->GetMean(),
                  vec_h_tts[12]->GetMean(),
                  vec_h_tts[13]->GetMean() };
  double ex_sza0[n] = { 40., 40., 40., 40., 40., 40. };
  double ey_sza0[n] = { vec_h_tts[10]->GetMeanError(),
                   vec_h_tts[11]->GetMeanError(),
                   vec_h_tts[0]->GetMeanError(),
                   vec_h_tts[1]->GetMeanError(),
                   vec_h_tts[12]->GetMeanError(),
                   vec_h_tts[13]->GetMeanError() };
  TGraphErrors* gr_sza0 = new TGraphErrors(n,x_sza0,y_sza0,ex_sza0,ey_sza0);
  gr_sza0->SetTitle( "Air | Cover off; Bz [mG]; Mean TTS [ns]" );
  gr_sza0->SetMarkerColor(kMagenta+3);
  gr_sza0->SetMarkerStyle(21);
  gr_sza0->Draw("AP");
  gPad->Modified();
  gPad->Update();
  plotname = string("ptf_timing_analysis_tts_vs_bfield_air_coveroff_z.pdf");
  tc->SaveAs(plotname.c_str(),"pdf");

  // x-direction, water, cover off
  double x_sxw0[nwx] = { -100., -50., 0., 0., 50. };
  double y_sxw0[nwx] = { vec_h_tts[16]->GetMean(),
                  vec_h_tts[17]->GetMean(),
                  vec_h_tts[14]->GetMean(),
                  vec_h_tts[15]->GetMean(),
                  vec_h_tts[18]->GetMean() };
  double ex_sxw0[nwx] = { 10., 10., 10., 10., 10. };
  double ey_sxw0[nwx] = { vec_h_tts[16]->GetMeanError(),
                   vec_h_tts[17]->GetMeanError(),
                   vec_h_tts[14]->GetMeanError(),
                   vec_h_tts[15]->GetMeanError(),
                   vec_h_tts[18]->GetMeanError() };
  TGraphErrors* gr_sxw0 = new TGraphErrors(nwx,x_sxw0,y_sxw0,ex_sxw0,ey_sxw0);
  gr_sxw0->SetTitle( "Water | Cover off; Bx [mG]; Mean TTS [ns]" );
  gr_sxw0->SetMarkerColor(kOrange-3);
  gr_sxw0->SetMarkerStyle(21);
  gr_sxw0->Draw("AP");
  gPad->Modified();
  gPad->Update();
  plotname = string("ptf_timing_analysis_tts_vs_bfield_water_coveroff_x.pdf");
  tc->SaveAs(plotname.c_str(),"pdf");

  // y-direction, water, cover off
  double x_syw0[n] = { -100., -50., 0., 0., 50., 100. };
  double y_syw0[n] = { vec_h_tts[20]->GetMean(),
                  vec_h_tts[21]->GetMean(),
                  vec_h_tts[14]->GetMean(),
                  vec_h_tts[15]->GetMean(),
                  vec_h_tts[22]->GetMean(),
                  vec_h_tts[23]->GetMean() };
  double ex_syw0[n] = { 10., 10., 10., 10., 10., 10. };
  double ey_syw0[n] = { vec_h_tts[20]->GetMeanError(),
                   vec_h_tts[21]->GetMeanError(),
                   vec_h_tts[14]->GetMeanError(),
                   vec_h_tts[15]->GetMeanError(),
                   vec_h_tts[22]->GetMeanError(),
                   vec_h_tts[23]->GetMeanError() };
  TGraphErrors* gr_syw0 = new TGraphErrors(n,x_syw0,y_syw0,ex_syw0,ey_syw0);
  gr_syw0->SetTitle( "Water | Cover off; By [mG]; Mean TTS [ns]" );
  gr_syw0->SetMarkerColor(kOrange-3);
  gr_syw0->SetMarkerStyle(21);
  gr_syw0->Draw("AP");
  gPad->Modified();
  gPad->Update();
  plotname = string("ptf_timing_analysis_tts_vs_bfield_water_coveroff_y.pdf");
  tc->SaveAs(plotname.c_str(),"pdf");

  // TMultiGraph
  // RTT:
  // Air X vs Water X
  // Air Y vs Water Y
  // TTS:
  // Air X vs Water X
  // Air Y vs Water Y

  TMultiGraph *mg_rx = new TMultiGraph();
  gr_rxa0->SetTitle( "Air | Cover off" );
  gr_rxw0->SetTitle( "Water | Cover off" );
  mg_rx->Add(gr_rxa0);
  mg_rx->Add(gr_rxw0);
  mg_rx->SetTitle( "; Bx [mG]; RTT RMS [ns]" );
  mg_rx->Draw("AP");
  tc->BuildLegend(0.2, 0.2, 0.45, 0.35);
  gPad->Modified();
  gPad->Update();
  plotname = string("ptf_timing_analysis_rtt_vs_bfield_coveroff_x.pdf");
  tc->SaveAs(plotname.c_str(),"pdf");

  TMultiGraph *mg_ry = new TMultiGraph();
  gr_rya0->SetTitle( "Air | Cover off" );
  gr_ryw0->SetTitle( "Water | Cover off" );
  mg_ry->Add(gr_rya0);
  mg_ry->Add(gr_ryw0);
  mg_ry->SetTitle( "; By [mG]; RTT RMS [ns]" );
  mg_ry->Draw("AP");
  tc->BuildLegend(0.55, 0.65, 0.8, 0.8);
  gPad->Modified();
  gPad->Update();
  plotname = string("ptf_timing_analysis_rtt_vs_bfield_coveroff_y.pdf");
  tc->SaveAs(plotname.c_str(),"pdf");

  TMultiGraph *mg_sx = new TMultiGraph();
  gr_sxa0->SetTitle( "Air | Cover off" );
  gr_sxw0->SetTitle( "Water | Cover off" );
  mg_sx->Add(gr_sxa0);
  mg_sx->Add(gr_sxw0);
  mg_sx->SetTitle( "; Bx [mG]; Mean TTS [ns]" );
  mg_sx->Draw("AP");
  tc->BuildLegend(0.2, 0.65, 0.45, 0.8);
  gPad->Modified();
  gPad->Update();
  plotname = string("ptf_timing_analysis_tts_vs_bfield_coveroff_x.pdf");
  tc->SaveAs(plotname.c_str(),"pdf");

  TMultiGraph *mg_sy = new TMultiGraph();
  gr_sya0->SetTitle( "Air | Cover off" );
  gr_syw0->SetTitle( "Water | Cover off" );
  mg_sy->Add(gr_sya0);
  mg_sy->Add(gr_syw0);
  mg_sy->SetTitle( "; By [mG]; Mean TTS [ns]" );
  mg_sy->Draw("AP");
  tc->BuildLegend(0.2, 0.65, 0.45, 0.8);
  gPad->Modified();
  gPad->Update();
  plotname = string("ptf_timing_analysis_tts_vs_bfield_coveroff_y.pdf");
  tc->SaveAs(plotname.c_str(),"pdf");
}
