#include "TStyle.h"
#include "TColor.h"
#include "TH1.h"
#include "TPad.h"

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
  TColor::InvertPalette(); // invert color palette

  return ptfStyle;
}

void CenterHistoTitles(TH1 *thisHisto){
  thisHisto->GetXaxis()->CenterTitle();
  thisHisto->GetYaxis()->CenterTitle();
  thisHisto->GetZaxis()->CenterTitle();
}

int AddGridLinesToPad(TPad *thisPad) {
  thisPad->SetGridx();
  thisPad->SetGridy();
  return 0;
}

