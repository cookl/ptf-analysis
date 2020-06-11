#include "TFile.h"
#include "TCanvas.h"
#include "TH2D.h"
#include "TROOT.h"
#include "TStyle.h"
#include "TColor.h"
#include "TPad.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TVector2.h"
#include "TAxis.h"

#include <iostream>
#include <stdlib.h>

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

// Calculate mean of vector
double fmean(vector<double> v){
  double sum = 0.;
  for( unsigned int i=0; i<v.size(); i++ ){
    sum += v[i];
  }
  if( v.size() > 0 )
    return sum / (double)v.size();
  else
    return sum;
}

// Calculate variance of vector
double fvar(vector<double> v){
  double var = 0.;
  double mean = fmean(v);
  for( unsigned int i=0; i<v.size(); i++ ){
    var += (v[i]-mean) * (v[i]-mean);
  }
  if( v.size() > 0 )
    return var / (double)v.size();
  else
    return var;
}

void ratio_plot(){

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
  gStyle->SetPadGridX(0);
  gStyle->SetPadGridY(0);

  TFile *f_on = new TFile("ptf_qe_analysis_run04581.root"); //4551 or 4581
  TH2D *h_ratio_plot = (TH2D*)f_on->Get("pmt0_qe_corr");
  
  //Translate plot down 7 bins in y
  //int offset = 4; // for 4562 vs 4551
  int offset = 0;
  for( int i=1; i<=h_ratio_plot->GetNbinsX(); i++ ){
    for( int j=1; j<=h_ratio_plot->GetNbinsY(); j++ ){
      if( j<=h_ratio_plot->GetNbinsY()-offset ){
        h_ratio_plot->SetBinContent(i,j,h_ratio_plot->GetBinContent(i,j+offset));
      }
      else{
        h_ratio_plot->SetBinContent(i,j,0.);
      }
    }
  }

  TFile *f_off = new TFile("ptf_qe_analysis_run04584.root"); //4562 or 4584
  TH2D *h_pmt0_cover_off = (TH2D*)f_off->Get("pmt0_qe_corr");
  //h_ratio_plot->Divide(h_pmt0_cover_off);
  
  //Loop through cover on bins
  //Get x and y min and max
  //Loop through cover off bins
  //Add all bins inside range
  //Divide by number of bins
  //Divide cover on by cover off
  for( int i=1; i<=h_ratio_plot->GetNbinsX(); i++ ){
    for( int j=1; j<=h_ratio_plot->GetNbinsY(); j++ ){
      double xl = ((TAxis*)h_ratio_plot->GetXaxis())->GetBinLowEdge(i);
      double xu = ((TAxis*)h_ratio_plot->GetXaxis())->GetBinUpEdge(i);
      double yl = ((TAxis*)h_ratio_plot->GetYaxis())->GetBinLowEdge(j);
      double yu = ((TAxis*)h_ratio_plot->GetYaxis())->GetBinUpEdge(j);
      double av = 0.;
      int nbins = 0;
      for( int k=1; k<=h_pmt0_cover_off->GetNbinsX(); k++ ){
        for( int l=1; l<=h_pmt0_cover_off->GetNbinsY(); l++ ){
          double bcx = ((TAxis*)h_pmt0_cover_off->GetXaxis())->GetBinCenter(k);
          double bcy = ((TAxis*)h_pmt0_cover_off->GetYaxis())->GetBinCenter(l);
          if( bcx > xl && bcx < xu && bcy > yl && bcy < yu ){
            av += h_pmt0_cover_off->GetBinContent(k,l);
            nbins++;
          }
        }
      }
      if( nbins > 0 ) av = av/(double)nbins;
      double num = h_ratio_plot->GetBinContent(i,j);
      if( av > 0.001 ){
        h_ratio_plot->SetBinContent(i,j,num/av);
      }
      else{
        h_ratio_plot->SetBinContent(i,j,0.);
      }
    }
  }
  
  h_ratio_plot->SetTitle("Cover on:off detection efficiency (corrected)");
  h_ratio_plot->SetMinimum(0.8);
  h_ratio_plot->SetMaximum(1.2);

  // x=0.39, y=0.35
  // length 0.29
  // for 0->0.29
  // loop through bins
  // sum content if bin centre within radius
 
  const Int_t n = 30;
  Double_t x[n], y[n], ex[n], ey[n];
  for( int r=0; r<n; r++ ){
    //if( r > 1 ) exit (EXIT_FAILURE);
    double radius = 0.01 + 0.22*(double)r/(double)(n-1);
    //cout << "radius " << radius << endl;
    //double sum = 0.;
    //int bins = 0;
    vector<double> bins;
    for( int i=1; i<=h_ratio_plot->GetNbinsX(); i++ ){
      for( int j=1; j<=h_ratio_plot->GetNbinsY(); j++ ){
        double bcx = ((TAxis*)h_ratio_plot->GetXaxis())->GetBinCenter(i);
        double bcy = ((TAxis*)h_ratio_plot->GetYaxis())->GetBinCenter(j);
        TVector2 binVec(bcx-0.39,bcy-0.35);
        if( binVec.Mod() < radius ){
          //cout << "bcx: " << bcx << endl;
          //cout << "bcy: " << bcy << endl;
          //cout << "mod: " << binVec.Mod() << endl;
          //cout << "bin content: " << h_ratio_plot->GetBinContent(i,j) << endl;
          //sum += h_ratio_plot->GetBinContent(i,j);
          //bins++;
          bins.push_back( h_ratio_plot->GetBinContent(i,j) );
        }
      }
    }
    x[r] = radius;
    //y[r] = sum/(double)bins;
    y[r] = fmean(bins);
    ex[r] = 0.;
    if( bins.size() > 0 )
      ey[r] = sqrt( fvar(bins) / bins.size() );
    else
      ey[r] = 0.;
    //cout << "sum: " << sum << endl;
    //cout << "bins: " << bins << endl;
  }

  TGraphErrors *g_average_ratio = new TGraphErrors(n,x,y,ex,ey);
  //g_average_ratio->SetLineColor(kMagenta+3);
  //g_average_ratio->SetLineWidth(4);
  //g_average_ratio->SetMarkerColor(kMagenta+3);
  //g_average_ratio->SetMarkerStyle(21);
  g_average_ratio->SetTitle("Average ratio cover on:off");
  g_average_ratio->GetXaxis()->SetTitle("Radius [m]");
  g_average_ratio->GetYaxis()->SetTitle("Ratio");

  //Make plots
  TCanvas* c = new TCanvas("canvas");
  h_ratio_plot->Draw("colz0");
  string plotname = string("ptf_qe_analysis_ratio.pdf");
  c->SaveAs(plotname.c_str(),"pdf");

  //g_average_ratio->Draw("AC");
  g_average_ratio->Draw("AP");
  gPad->Modified();
  gPad->Update();
  plotname = string("ptf_qe_analysis_ratio_average.pdf");
  c->SaveAs(plotname.c_str(),"pdf");

}
