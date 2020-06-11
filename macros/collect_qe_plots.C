/// collect_qe_plots
/// John Walker (Apr. 2020)
///
/// Make uniform set of plots from each PTF run
/// No inputs.  All files hard coded.
/// Input files are from output of ptf_qe_analysis.app

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
  ptfStyle->SetPalette(kSunset); // use the sunset color set
  //ptfStyle->SetPalette(kColorPrintableOnGrey); // use the colorPrintableOnGrey color set
  //ptfStyle->SetPalette(kCubehelix); // use the cubehelix color set
  //ptfStyle->SetPalette(kLightTemperature);
  TColor::InvertPalette(); // invert color palette

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


void collect_qe_plots(){

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
    new TFile( "ptf_qe_analysis_run04427.root" ),
    new TFile( "ptf_qe_analysis_run04431.root" ),
    new TFile( "ptf_qe_analysis_run04446.root" ),
    new TFile( "ptf_qe_analysis_run04447.root" ),
    new TFile( "ptf_qe_analysis_run04434.root" ),
    new TFile( "ptf_qe_analysis_run04433.root" ),
    new TFile( "ptf_qe_analysis_run04458.root" ),
    new TFile( "ptf_qe_analysis_run04459.root" ),
    new TFile( "ptf_qe_analysis_run04457.root" ),
    new TFile( "ptf_qe_analysis_run04448.root" ),
    new TFile( "ptf_qe_analysis_run04444.root" ),
    new TFile( "ptf_qe_analysis_run04445.root" ),
    new TFile( "ptf_qe_analysis_run04443.root" ),
    new TFile( "ptf_qe_analysis_run04442.root" ),

    new TFile( "ptf_qe_analysis_run04525.root" ),
    new TFile( "ptf_qe_analysis_run04526.root" ),
    new TFile( "ptf_qe_analysis_run04533.root" ),
    new TFile( "ptf_qe_analysis_run04534.root" ),
    new TFile( "ptf_qe_analysis_run04531.root" ), 
    new TFile( "ptf_qe_analysis_run04530.root" ),
    new TFile( "ptf_qe_analysis_run04538.root" ),
    new TFile( "ptf_qe_analysis_run04540.root" ),
    new TFile( "ptf_qe_analysis_run04537.root" ),
    new TFile( "ptf_qe_analysis_run04536.root" ),

    new TFile( "ptf_qe_analysis_run04562.root" ),
    new TFile( "ptf_qe_analysis_run04551.root" ),
    new TFile( "ptf_qe_analysis_run04584.root" ),
    new TFile( "ptf_qe_analysis_run04581.root" )
  };


  //KEY: TH2D pmt0_qe_corr;  Detection efficiency (corrected)
  std::vector< TH1D* > vec_pmt0_qe_corr;
  std::vector< TH1D* > vec_pmt1_qe;
  for ( unsigned i=0 ; i < bscan2_files.size() ; ++i ){
    TH2D* pmt0_qe_corr_2d = (TH2D*)bscan2_files[i]->Get("pmt0_qe_corr");
    TH2D* pmt1_qe_2d = (TH2D*)bscan2_files[i]->Get("pmt1_qe");

    TH1D* pmt0_qe_corr = make_1d_from_2d( pmt0_qe_corr_2d, " ; QE ; scan points / bin", 100, 0.00, 0.4 );
    TH1D* pmt1_qe = make_1d_from_2d( pmt1_qe_2d, " ; QE ; scan points / bin", 30, 0.3, 0.7 );
    double scale = pmt0_qe_corr->GetMaximum();
    double scale_pmt1 = pmt1_qe->GetMaximum();
    pmt0_qe_corr->Scale( 1.0/scale );
    pmt1_qe->Scale( 1.0/scale_pmt1 );

    vec_pmt0_qe_corr.push_back( pmt0_qe_corr );
    vec_pmt1_qe.push_back( pmt1_qe );

    ostringstream os;
    os << "tc_" << i;
    TCanvas * tc = new TCanvas( os.str().c_str(), bscan2_tags[i].c_str() );
    pmt0_qe_corr->Draw("e1");

    string fname(bscan2_files[i]->GetName());
    fname.erase(fname.length()-5);
    ostringstream os2;
    os2 << fname << "_pmt0corr_1d.pdf";
    tc->Print( os2.str().c_str() );
  }

  // Plot to compare the 1D plots for 0mG
  TCanvas * tc_comp = new TCanvas( "tc_comp", "tc_comp" );
  vec_pmt0_qe_corr[14]->Draw("e1");
  vec_pmt0_qe_corr[15]->SetMarkerColor(kMagenta+3);
  vec_pmt0_qe_corr[15]->SetLineColor(kMagenta+3);
  vec_pmt0_qe_corr[15]->Draw("e1 same");
  string plotname_comp = string("ptf_qe_analysis_0mg_1d_comparison_pmt0.pdf");
  tc_comp->SaveAs(plotname_comp.c_str(),"pdf");
  TCanvas * tc_comp_pmt1 = new TCanvas( "tc_comp_pmt1", "tc_comp_pmt1" );
  vec_pmt1_qe[14]->Draw("le1");
  vec_pmt1_qe[15]->SetMarkerColor(kMagenta+3);
  vec_pmt1_qe[15]->SetLineColor(kMagenta+3);
  vec_pmt1_qe[15]->Draw("le1 same");
  string plotname_comp_pmt1 = string("ptf_qe_analysis_0mg_1d_comparison_pmt1.pdf");
  tc_comp_pmt1->SaveAs(plotname_comp_pmt1.c_str(),"pdf");

  // print table of results from the TH1D's
  cout<< "\\begin{center}" << endl;
  cout<< "\\begin{tabular}{ l | c }" << endl;
  cout << "Run and Description\t& QE \\\\ \\hline" <<endl;
  for ( unsigned i=0 ; i < bscan2_files.size() ; ++i ){

    cout<< bscan2_tags[i] <<"\t& ";
    cout<< fixed << setprecision(4) << setw(12 )<< vec_pmt0_qe_corr[i]->GetMean() << " $\\pm$ " << vec_pmt0_qe_corr[i]->GetMeanError() ;
    cout<< "\\\\ "<<endl;
  }
  cout<< "\\end{tabular}" <<endl;
  cout<< "\\end{center}" <<endl;


  // Plot average QE vs magnetic field
  // x-direction, air, cover off
  TCanvas * tc = new TCanvas( "av_qe_vs_bfield", "av_qe_vs_bfield" );
  //tc->SetGrid();
  const int n = 6;
  double x_xa0[n] = { -100., -50., 0., 0., 50., 100. };
  double y_xa0[n] = { vec_pmt0_qe_corr[2]->GetMean(),
                  vec_pmt0_qe_corr[3]->GetMean(),
                  vec_pmt0_qe_corr[0]->GetMean(),
                  vec_pmt0_qe_corr[1]->GetMean(),
                  vec_pmt0_qe_corr[4]->GetMean(),
                  vec_pmt0_qe_corr[5]->GetMean() };
  double ex_xa0[n] = { 10., 10., 10., 10., 10., 10. };
  double ey_xa0[n] = { vec_pmt0_qe_corr[2]->GetMeanError(),
                   vec_pmt0_qe_corr[3]->GetMeanError(),
                   vec_pmt0_qe_corr[0]->GetMeanError(),
                   vec_pmt0_qe_corr[1]->GetMeanError(),
                   vec_pmt0_qe_corr[4]->GetMeanError(),
                   vec_pmt0_qe_corr[5]->GetMeanError() };
  TGraphErrors* gr_xa0 = new TGraphErrors(n,x_xa0,y_xa0,ex_xa0,ey_xa0);
  gr_xa0->SetTitle( "Air | Cover off; Bx [mG]; Mean QE" );
  gr_xa0->SetMarkerColor(kMagenta+3);
  gr_xa0->SetMarkerStyle(21);
  gr_xa0->Draw("AP");
  string plotname = string("ptf_qe_analysis_vs_bfield_air_coveroff_x.pdf");
  tc->SaveAs(plotname.c_str(),"pdf");

  // y-direction, air, cover off
  double x_ya0[n] = { -100., -50., 0., 0., 50., 100. };
  double y_ya0[n] = { vec_pmt0_qe_corr[6]->GetMean(),
                  vec_pmt0_qe_corr[7]->GetMean(),
                  vec_pmt0_qe_corr[0]->GetMean(),
                  vec_pmt0_qe_corr[1]->GetMean(),
                  vec_pmt0_qe_corr[8]->GetMean(),
                  vec_pmt0_qe_corr[9]->GetMean() };
  double ex_ya0[n] = { 10., 10., 10., 10., 10., 10. };
  double ey_ya0[n] = { vec_pmt0_qe_corr[6]->GetMeanError(),
                   vec_pmt0_qe_corr[7]->GetMeanError(),
                   vec_pmt0_qe_corr[0]->GetMeanError(),
                   vec_pmt0_qe_corr[1]->GetMeanError(),
                   vec_pmt0_qe_corr[8]->GetMeanError(),
                   vec_pmt0_qe_corr[9]->GetMeanError() };
  TGraphErrors* gr_ya0 = new TGraphErrors(n,x_ya0,y_ya0,ex_ya0,ey_ya0);
  gr_ya0->SetTitle( "Air | Cover off; By [mG]; Mean QE" );
  gr_ya0->SetMarkerColor(kMagenta+3);
  gr_ya0->SetMarkerStyle(21);
  gr_ya0->Draw("AP");
  gPad->Modified();
  gPad->Update();
  plotname = string("ptf_qe_analysis_vs_bfield_air_coveroff_y.pdf");
  tc->SaveAs(plotname.c_str(),"pdf");

  // z-direction, air, cover off
  double x_za0[n] = { -100., -50., 0., 0., 50., 100. };
  double y_za0[n] = { vec_pmt0_qe_corr[10]->GetMean(),
                  vec_pmt0_qe_corr[11]->GetMean(),
                  vec_pmt0_qe_corr[0]->GetMean(),
                  vec_pmt0_qe_corr[1]->GetMean(),
                  vec_pmt0_qe_corr[12]->GetMean(),
                  vec_pmt0_qe_corr[13]->GetMean() };
  double ex_za0[n] = { 40., 40., 40., 40., 40., 40. };
  double ey_za0[n] = { vec_pmt0_qe_corr[10]->GetMeanError(),
                   vec_pmt0_qe_corr[11]->GetMeanError(),
                   vec_pmt0_qe_corr[0]->GetMeanError(),
                   vec_pmt0_qe_corr[1]->GetMeanError(),
                   vec_pmt0_qe_corr[12]->GetMeanError(),
                   vec_pmt0_qe_corr[13]->GetMeanError() };
  TGraphErrors* gr_za0 = new TGraphErrors(n,x_za0,y_za0,ex_za0,ey_za0);
  gr_za0->SetTitle( "Air | Cover off; Bz [mG]; Mean QE" );
  gr_za0->SetMarkerColor(kMagenta+3);
  gr_za0->SetMarkerStyle(21);
  gr_za0->Draw("AP");
  gPad->Modified();
  gPad->Update();
  plotname = string("ptf_qe_analysis_vs_bfield_air_coveroff_z.pdf");
  tc->SaveAs(plotname.c_str(),"pdf");

  // x-direction, water, cover off
  const int nwx = 5;
  double x_xw0[nwx] = { -100., -50., 0., 0., 50. };
  double y_xw0[nwx] = { vec_pmt0_qe_corr[16]->GetMean(),
                  vec_pmt0_qe_corr[17]->GetMean(),
                  vec_pmt0_qe_corr[14]->GetMean(),
                  vec_pmt0_qe_corr[15]->GetMean(),
                  vec_pmt0_qe_corr[18]->GetMean() };
  double ex_xw0[nwx] = { 10., 10., 10., 10., 10. };
  double ey_xw0[nwx] = { vec_pmt0_qe_corr[16]->GetMeanError(),
                   vec_pmt0_qe_corr[17]->GetMeanError(),
                   vec_pmt0_qe_corr[14]->GetMeanError(),
                   vec_pmt0_qe_corr[15]->GetMeanError(),
                   vec_pmt0_qe_corr[18]->GetMeanError() };
  TGraphErrors* gr_xw0 = new TGraphErrors(nwx,x_xw0,y_xw0,ex_xw0,ey_xw0);
  gr_xw0->SetTitle( "Water | Cover off; Bx [mG]; Mean QE" );
  gr_xw0->SetMarkerColor(kOrange-3);
  gr_xw0->SetMarkerStyle(21);
  gr_xw0->Draw("AP");
  gPad->Modified();
  gPad->Update();
  plotname = string("ptf_qe_analysis_vs_bfield_water_coveroff_x.pdf");
  tc->SaveAs(plotname.c_str(),"pdf");

  // y-direction, water, cover off
  double x_yw0[n] = { -100., -50., 0., 0., 50., 100. };
  double y_yw0[n] = { vec_pmt0_qe_corr[20]->GetMean(),
                  vec_pmt0_qe_corr[21]->GetMean(),
                  vec_pmt0_qe_corr[14]->GetMean(),
                  vec_pmt0_qe_corr[15]->GetMean(),
                  vec_pmt0_qe_corr[22]->GetMean(),
                  vec_pmt0_qe_corr[23]->GetMean() };
  double ex_yw0[n] = { 10., 10., 10., 10., 10., 10. };
  double ey_yw0[n] = { vec_pmt0_qe_corr[20]->GetMeanError(),
                   vec_pmt0_qe_corr[21]->GetMeanError(),
                   vec_pmt0_qe_corr[14]->GetMeanError(),
                   vec_pmt0_qe_corr[15]->GetMeanError(),
                   vec_pmt0_qe_corr[22]->GetMeanError(),
                   vec_pmt0_qe_corr[23]->GetMeanError() };
  TGraphErrors* gr_yw0 = new TGraphErrors(n,x_yw0,y_yw0,ex_yw0,ey_yw0);
  gr_yw0->SetTitle( "Water | Cover off; By [mG]; Mean QE" );
  gr_yw0->SetMarkerColor(kOrange-3);
  gr_yw0->SetMarkerStyle(21);
  gr_yw0->Draw("AP");
  gPad->Modified();
  gPad->Update();
  plotname = string("ptf_qe_analysis_vs_bfield_water_coveroff_y.pdf");
  tc->SaveAs(plotname.c_str(),"pdf");

  // TMultiGraph
  // RTT:
  // Air X vs Water X
  // Air Y vs Water Y
  // TTS:
  // Air X vs Water X
  // Air Y vs Water Y

  TMultiGraph *mg_x = new TMultiGraph();
  gr_xa0->SetTitle( "Air | Cover off" );
  gr_xw0->SetTitle( "Water | Cover off" );
  mg_x->Add(gr_xa0);
  mg_x->Add(gr_xw0);
  mg_x->SetTitle( "; Bx [mG]; Mean QE" );
  mg_x->Draw("AP");
  tc->BuildLegend(0.2, 0.5, 0.45, 0.65);
  gPad->Modified();
  gPad->Update();
  plotname = string("ptf_qe_analysis_vs_bfield_coveroff_x.pdf");
  tc->SaveAs(plotname.c_str(),"pdf");

  TMultiGraph *mg_y = new TMultiGraph();
  gr_ya0->SetTitle( "Air | Cover off" );
  gr_yw0->SetTitle( "Water | Cover off" );
  mg_y->Add(gr_ya0);
  mg_y->Add(gr_yw0);
  mg_y->SetTitle( "; By [mG]; Mean QE" );
  mg_y->Draw("AP");
  tc->BuildLegend(0.2, 0.5, 0.45, 0.65);
  gPad->Modified();
  gPad->Update();
  plotname = string("ptf_qe_analysis_vs_bfield_coveroff_y.pdf");
  tc->SaveAs(plotname.c_str(),"pdf");

}
