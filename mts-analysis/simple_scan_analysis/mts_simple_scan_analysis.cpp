// Ashley/Thomas/Angela 2D histograms
// 2022-05-27
#include "WaveformFitResult.hpp"
#include "ScanPoint.hpp"
#include "TCanvas.h"
#include "TFile.h"
#include "TF1.h"
#include "TF2.h"
#include "TMath.h"
#include "TH1D.h"
#include "TH2F.h"
#include "THStack.h"
#include "TGraph.h"
#include "TGaxis.h"
#include "TStyle.h"
#include "TLatex.h"
#include "TLine.h"
#include "TLegend.h"
#include "TPaveStats.h"
#include "TPad.h"
#include "TProfile.h"
#include "TH2F.h"
#include "TEllipse.h"
#include "TPaletteAxis.h"
#include <iostream>
#include <vector>
#include <string>
#include <cstring>
//using namespace std;

int main( int argc, char* argv[] )
{//START MAIN

  if ( argc != 2 ){
    std::cerr<<"Usage: ptf_ttree_analysis.app ptf_analysis.root\n";
    exit(0); }
  TFile * fin = new TFile( argv[1], "read" );
  TTree * tt0;    
  WaveformFitResult * wf0; 
  // Getting the waveform fit TTree for desired channel
  tt0 = (TTree*)fin->Get("ptfanalysis2");//select channel here
  wf0 = new WaveformFitResult;
  if(tt0) wf0->SetBranchAddresses( tt0 );
  
  //Define bin paramaters for scan. (m)
  /*---------------------------------------------------------------------------------*/
  float step_size = 0.005;
  float xstart = 0.42;
  float ystart = 0.45;
  float x_scan_dist = 0.12;
  float y_scan_dist = 0.105;
  float scnsize = step_size/2.0;//half of the x/y square you want to look over
  
  int num_bins_x = static_cast<int>(x_scan_dist/step_size)+1 ;
  int num_bins_y = static_cast<int>(y_scan_dist/step_size)+1 ;
  float x_low = xstart - step_size/2;
  float x_high = xstart + x_scan_dist + step_size/2;
  float y_low = ystart - step_size/2;;
  float y_high = ystart + y_scan_dist + step_size/2;

  std::cout << "There are approx" << num_bins_x * num_bins_y << "scan points" << std::endl;
  std::cout << "We expect" << 3500 * num_bins_x * num_bins_y << "events in total" << std::endl;
  
  /*----------------------------------------------------------------------------------*/
  //Initialize
  /* --------------------------------------------------------------------------------- */
  //2D Histograms; 
  TH2F *h_mph = new TH2F("mean pulse height distribution","Mean Pulse Height",num_bins_x,x_low,x_high,num_bins_y,y_low,y_high);
  TH2F *h_p = new TH2F("Number of Pulses","Number of Pulses",num_bins_x,x_low,x_high,num_bins_y,y_low,y_high);
  TH2F *h_RMS = new TH2F("RMS", "Standard dev. dist.",num_bins_x,x_low,x_high,num_bins_y,y_low,y_high);
  TH2F *h_scan_pt = new TH2F("scan point","scan point",num_bins_x,x_low,x_high,num_bins_y,y_low,y_high);
  TH2F *h_events = new TH2F("Events","Number of events",num_bins_x,x_low,x_high,num_bins_y,y_low,y_high);
  TH2F *h_eff = new TH2F("Efficiency","pulses/events",num_bins_x,x_low,x_high,num_bins_y,y_low,y_high);
  
  //2D array bins for constructing the 2D histograms;
  TH1F *h[num_bins_x][num_bins_y];
  Double_t events_bin[num_bins_x][num_bins_y];
  for (int x=0;x<num_bins_x;x++) {
    for (int y=0;y<num_bins_y;y++){
      char name1[100];
      sprintf(name1,"h bin %i %i",x,y);
      h[x][y] = new TH1F(name1,"pulse heights",200,0,200*0.48828125);
      events_bin[x][y]=0;
    }
  }
  /* --------------------------------------------------------------------------------- */

  
  std::cout << "Analyzing " << tt0->GetEntries() << " waveforms" << std::endl;
  for(int i = 0; i < tt0->GetEntries()-1 ; i++)
    {//START WAVEFORM LOOP
      
      tt0->GetEvent(i);
      //std::cout << "Number of pulses found: " << wf0->numPulses << std::endl;

      for(int ypoint = 0; ypoint < num_bins_y ; ypoint++)
	{//START YLOOP

	  float ycenter = step_size*ypoint + ystart;
          float y_l = ycenter - scnsize;
          float y_r = ycenter + scnsize;

	  for(int xpoint = 0; xpoint < num_bins_x ; xpoint++)
	    {//START XLOOP
	      float xcenter = step_size*xpoint + xstart;
	      float x_l = xcenter - scnsize;
	      float x_r = xcenter + scnsize;

	      if((wf0->x>x_l && wf0->x < x_r) and (wf0->y >y_l && wf0->y < y_r))
		{//FILTER XY COORDS

		  events_bin[xpoint][ypoint] += 1;

		  for(int k = 0; k < wf0->numPulses; k++ )
		    {//START PULSE LOOP

		      if(wf0->pulseTimes[k] > 2200 and wf0->pulseTimes[k] < 2400 and wf0->pulseCharges[k]*1000.0 > 2.0)
			{//FILTER PULSE TIME
		      
			  h[xpoint][ypoint]->Fill(wf0->pulseCharges[k]*1000.0);
			  h_scan_pt->SetBinContent(xpoint+1,ypoint+1,wf0->scanpt);
			  			 
			}//DONE FILTER PULSE TIME

		    }//DONE PULSE LOOP
		  
	      	}//DONE FILTER XY COORDS
		 
	    }//DONE XLOOP
		    
	}//DONE YLOOP

    }//DONE WAVEFORM LOOP
  
 	  
  //Constructing 2D histograms from h (array of histograms)
   /* --------------------------------------------------------------------------------- */

  for(int x = 0; x < num_bins_x ; x++){
      for(int y = 0; y < num_bins_y ; y++){
       	  h_p->SetBinContent(x+1,y+1,h[x][y]->GetEntries());
	  h_mph->SetBinContent(x+1,y+1,h[x][y]->GetMean());
	  h_RMS->SetBinContent(x+1,y+1,h[x][y]->GetRMS());
	  h_events->SetBinContent(x+1,y+1,events_bin[x][y]);
      }} 
  /* --------------------------------------------------------------------------------- */

  
  //Plotting
  /* --------------------------------------------------------------------------------- */
  //Style;
  gStyle->SetPalette(1);
  gStyle->SetOptTitle(1); gStyle->SetOptStat(0);
  gStyle->SetOptFit(1111); gStyle->SetStatBorderSize(0);
  gStyle->SetStatX(.89); gStyle->SetStatY(.89);

  //Draw circle to outline PMT around maximum h_p value;
  h_eff = (TH2F*)h_p->Clone();
  h_eff->Divide(h_events);
  h_eff->Scale(100.0);
  Int_t MaxBin = h_eff->GetMaximumBin();
  Int_t x,y,z;
  h_eff->GetBinXYZ(MaxBin,x,y,z);
  std::cout << "MAXBIN  = (" << x << "," << y << "," << z << ")" << std::endl;
  double x_max = (static_cast<float>(x)-1.0)*step_size+xstart;
  double y_max = (static_cast<float>(y)-1.0)*step_size+ystart;
  std::cout << "MAXcoords  = (" << x_max << "," << y_max << ")" << std::endl;
  TEllipse *pmt = new TEllipse(0.474,0.493, 0.04,0.04);
  pmt->SetLineColor(1);
  pmt->SetFillStyle(0);
  
  /* --------------------------------------------------------------------------------- */

  //h_scan_pt histogram;
  TCanvas *c0 = new TCanvas("C0","",2000,1000);
  c0->SetRightMargin(0.2);
  c0->SetLeftMargin(0.15);
  c0->SetBottomMargin(0.10);
  h_scan_pt->Draw("TEXT");
  h_scan_pt->GetXaxis()->SetTitle("X (m)");
  h_scan_pt->GetYaxis()->SetTitle("Y (m)");
  h_scan_pt->GetZaxis()->SetTitle("Scan point");
  h_scan_pt->SetTitleSize(50,"t");
  c0->SaveAs("scan_points.png");
  

  //h_p histogram;
  TCanvas *c1 = new TCanvas("C1","",1200,500);
  c1->SetRightMargin(0.2);
  c1->SetLeftMargin(0.15);
  c1->SetBottomMargin(0.10);
  h_p->Draw("COLZ");
  h_p->GetXaxis()->SetTitle("X (m)");
  h_p->GetYaxis()->SetTitle("Y (m)");
  h_p->GetZaxis()->SetTitle("Pulse counts");
  h_p->SetTitleSize(50,"t");   
  pmt->Draw();
  c1->SetRealAspectRatio();
  c1->Modified();
  c1->Update();
  c1->SaveAs("2Dhist_p.png");
  std::cout << "p histogram  mean X  = " << h_p->GetMean(1) << std::endl;
  std::cout << "p histogram  mean Y  = " << h_p->GetMean(2) << std::endl;

  
  //h_mph histogram
  TCanvas *c2 = new TCanvas("C2","");
  c2->SetRightMargin(0.2);
  c2->SetLeftMargin(0.15);
  c2->SetBottomMargin(0.10); 
  h_mph->Draw("COLZ");
  h_mph->GetXaxis()->SetTitle("X (m)");
  h_mph->GetYaxis()->SetTitle("Y (m)");
  h_mph->GetZaxis()->SetTitle("Mean Pulse height (mV)");
  pmt->Draw();
  c2->SetRealAspectRatio();
  c2->Modified();
  c2->Update();
  c2->SaveAs("2Dhist_mph.png"); 
  std::cout << "mph histogram mean X = " << h_mph->GetMean(1) << std::endl;
  std::cout << "mph histogram mean Y = " << h_mph->GetMean(2) << std::endl;

 
  //h_RMS histogram;
  TCanvas *c3 = new TCanvas("C3","",1200,500);
  c3->SetRightMargin(0.2);
  c3->SetLeftMargin(0.15);
  c3->SetBottomMargin(0.10);
  h_RMS->Draw("COLZ");
  h_RMS->GetXaxis()->SetTitle("X (m)");
  h_RMS->GetYaxis()->SetTitle("Y (m)");
  h_RMS->GetZaxis()->SetTitle("RMS (mV)");
  h_RMS->SetTitleSize(50,"t"); 
  pmt->Draw();
  c3->SetRealAspectRatio();
  c3->Modified();
  c3->Update();
  c3->SaveAs("2Dhist_RMS.png");

  //h_eff histogram;
  TCanvas *c4 = new TCanvas("C4","",1200,500);
  c4->SetRightMargin(0.2);
  c4->SetLeftMargin(0.15);
  c4->SetBottomMargin(0.10);
  h_eff->GetXaxis()->SetTitle("X (m)");
  h_eff->GetYaxis()->SetTitle("Y (m)");
  h_eff->GetZaxis()->SetTitle("efficiency (%)");
  h_eff->SetTitle("Efficiency");
  h_eff->Draw("COLZ");
  h_eff->SetTitleSize(50,"t"); 
  pmt->Draw();
  c4->SetRealAspectRatio();
  c4->Modified();
  c4->Update();
  c4->SaveAs("2D_Efficiency.png");

  //h_events histogram;
  TCanvas *c5 = new TCanvas("C5","",1200,500);
  c5->SetRightMargin(0.2);
  c5->SetLeftMargin(0.15);
  c5->SetBottomMargin(0.10);
  h_events->GetXaxis()->SetTitle("X (m)");
  h_events->GetYaxis()->SetTitle("Y (m)");
  h_events->GetZaxis()->SetTitle("number of events");
  h_events->SetTitle("Events");
  h_events->Draw("COLZ");
  h_events->SetTitleSize(50,"t"); 
  pmt->Draw();
  c5->SetRealAspectRatio();
  c5->Modified();
  c5->Update();
  c5->SaveAs("2D_Events.png");
  
  
  return 0;

}//DONE MAIN

