// Ashley/Thomas simple pulse height histogram
// 2020-01-20


#include "WaveformFitResult.hpp"
#include "ScanPoint.hpp"
#include "TCanvas.h"
#include "TFile.h"
#include "TF1.h"
#include "TH1D.h"
#include "TH2F.h"
#include "THStack.h"
#include "TGraph.h"
#include "TGaxis.h"


#include "TF1.h"
#include "TStyle.h"
#include "TLatex.h"
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
   
  // Get the waveform fit TTree for desired channel
  tt0 = (TTree*)fin->Get("ptfanalysis2");
  wf0 = new WaveformFitResult;
  if(tt0) wf0->SetBranchAddresses( tt0 );

  // Define bin paramaters for scan
  float step_size = 0.01;//(m)
  float x_scan_dist = 0.2;//(m)
  float y_scan_dist = 0.24;//(m)
  
  int num_bins_x = static_cast<int>(x_scan_dist/step_size)+1 ;//number of bins for x
  int num_bins_y = static_cast<int>(y_scan_dist/step_size)+1 ;//number of bins for y
  float x_low = 0.375 - step_size/2;
  float x_high = 0.375 + x_scan_dist + step_size/2;
  float y_low = x_low;
  float y_high = 0.375 + y_scan_dist + step_size/2;
  float scnsize = step_size/2.0;//half of the x/y square you want to look over
  
  //Initialize 2d histogram for mean pulse height and pulse counts
  TH2F *h_mean = new TH2F("mean pulse height distribution","Mean Pulse Height",num_bins_x,x_low,x_high,num_bins_y,y_low,y_high);
  TH2F *h_counts = new TH2F("h3","Number of Pulses",num_bins_x,x_low,x_high,num_bins_y,y_low,y_high);

  //Initialize histogram for x/y scan points
  TH1F *xpos = new TH1F("xpos","xpos",5000,0.35,0.6);
  TH1F *ypos = new TH1F("ypos","ypos",5000,0.35,0.65);


  //Initialize bin for waveform pulse heights histogram
  TH1F *h1 = new TH1F("pulse_height","Pulse Height",200,0,200*0.48828125);
  int last_x = -1;
  int last_y = -1;

  std::cout << "Analyzing " << tt0->GetEntries() << " waveforms" << std::endl;
  for(int i = 0; i < tt0->GetEntries()-1 ; i++)
    {//START WAVEFORM LOOP
      
      tt0->GetEvent(i);
      //std::cout << "Number of pulses found: " << wf0->numPulses << std::endl;
	for(int k = 0; k < wf0->numPulses; k++ )
	  {//START PULSE LOOP

	    if(wf0->pulseTimes[k] > 2200 and wf0->pulseTimes[k] < 2400 and wf0->pulseCharges[k]*1000.0 > 4.5)
	      {//FILTER PULSE TIME

		xpos->Fill(wf0->x);
		ypos->Fill(wf0->y);
				
		for(int xpoint = 0; xpoint < num_bins_x ; xpoint++)
		  {//START XLOOP

		    float xcenter = step_size*xpoint + 0.375;
		    float x_bin_edge_l = xcenter - scnsize;
		    float x_bin_edge_r = xcenter + scnsize;
		    
		    for(int ypoint = 0; ypoint < num_bins_y ; ypoint++)
		      {//START YLOOP

			float ycenter = step_size*ypoint + 0.375;
			float y_bin_edge_l = ycenter - scnsize;
			float y_bin_edge_r = ycenter + scnsize;
						
			if((wf0->x>x_bin_edge_l && wf0->x < x_bin_edge_r) && (wf0->y >y_bin_edge_l && wf0->y < y_bin_edge_r))
			  {//FILTER XY COORDS
			    
			    if(xpoint == last_x and ypoint == last_y)
			      {h1->Fill(wf0->pulseCharges[k]*1000.0);}
			    else
			      {
				std::cout << "SCAN POINT " << wf0->scanpt << " : " << std::endl;
				std::cout << "X = " << x_bin_edge_l << " - " << x_bin_edge_r << std::endl;
				std::cout << "Y = " << y_bin_edge_l << " - " << y_bin_edge_r << std::endl;
				int Num_Entries = h1->GetEntries();
				float Mean_Entries = h1->GetMean();
				std::cout << "Num Entries " << Num_Entries << "|  Mean " << Mean_Entries << std::endl;
        			                                                          
				h_mean->SetBinContent(last_x+1,last_y+1,Mean_Entries);
				h_counts->SetBinContent(last_x+1,last_y+1,Num_Entries); 
		      	       	h1->Reset();
			      }
		        
			    last_x = xpoint;
			    last_y = ypoint;
	        

			  }//DONE FILTER XY COORDS

		      }//DONE YLOOP
		    
		  }//DONE XLOOP

	      }//DONE FILTER PULSE TIME

	  }//DONE PULSE LOOP

    }//DONE WAVEFORM LOOP

  
  //Style
  gStyle->SetPalette(1);
  gStyle->SetOptTitle(1); gStyle->SetOptStat(0);
  gStyle->SetOptFit(1111); gStyle->SetStatBorderSize(0);
  gStyle->SetStatX(.89); gStyle->SetStatY(.89);

  //Mean histogram
  TCanvas *c1 = new TCanvas("C1","",1200,500);
  c1->SetRightMargin(0.2);
  c1->SetLeftMargin(0.15);
  c1->SetBottomMargin(0.10);

  h_mean->Draw("COLZ");
  h_mean->GetXaxis()->SetTitle("X (m)");
  h_mean->GetYaxis()->SetTitle("Y (m)");
  h_mean->GetZaxis()->SetTitle("Pulse height (mV)");
  h_mean->SetTitleSize(50,"t");
  
  
  //Initialize circle for outlining pmt
  TEllipse *el_mean = new TEllipse(h_mean->GetMean(1), h_mean->GetMean(2), 0.04,0.04);
  el_mean->SetLineColor(4);
  el_mean->SetFillStyle(0);
  el_mean->Draw();
 
  c1->SetRealAspectRatio();
  c1->Modified();
  c1->Update();
  c1->SaveAs("mpmt_mean_colz_2d.png");
  
  std::cout << "mean counts histogram  mean x  = " << h_mean->GetMean(1) << std::endl;
  std::cout << "mean counts histogram  mean y  = " << h_mean->GetMean(2) << std::endl;

  //Counts histogram
  TCanvas *c2 = new TCanvas("C2","");
  c2->SetRightMargin(0.2);
  c2->SetLeftMargin(0.15);
  c2->SetBottomMargin(0.10);
  
  h_counts->Draw("COLZ");
  h_counts->GetXaxis()->SetTitle("X (m)");
  h_counts->GetYaxis()->SetTitle("Y (m)");
  h_counts->GetZaxis()->SetTitle("counts");
  
  el_mean->Draw();
  c2->SetRealAspectRatio();
  c2->Modified();
  c2->Update();
  c2->SaveAs("mpmt_entries_colz_2d.png");

  
  std::cout << "counts histogram mean x = " << h_counts->GetMean(1) << std::endl;
  std::cout << "counts histogram mean y = " << h_counts->GetMean(2) << std::endl;

  //Look at how we divide the bins

  Double_t x_bin_edge_l[num_bins_x],x_bin_edge_r[num_bins_x],x_border[num_bins_x],x_bins[num_bins_x];
  Double_t y_bin_edge_l[num_bins_y],y_bin_edge_r[num_bins_y],y_border[num_bins_y],y_bins[num_bins_y];

  
  for(Int_t xpoint = 0; xpoint < num_bins_x ; xpoint++)
    {//START XLOOP
      float xcenter = step_size*xpoint + 0.375;
      x_bins[xpoint] = xcenter;
      x_bin_edge_l[xpoint] = xcenter - scnsize;
      x_bin_edge_r[xpoint] = xcenter + scnsize;
      x_border[xpoint]=35000;
    }//END XLOOP
  TGraph *grx_l = new TGraph (num_bins_x, x_bin_edge_l, x_border);
  TGraph *grx_r = new TGraph (num_bins_x, x_bin_edge_r, x_border);
  TGraph *grx = new TGraph (num_bins_x, x_bins, x_border);
  
  for(Int_t ypoint = 0; ypoint < num_bins_y ; ypoint++)
    {//START YLOOP
      float ycenter = step_size*ypoint + 0.375;
      y_bins[ypoint] = ycenter;
      y_bin_edge_l[ypoint] = ycenter - scnsize;
      y_bin_edge_r[ypoint] = ycenter + scnsize;
      y_border[ypoint]= 80000;
    }//END YLOOP
  TGraph *gry_l = new TGraph (num_bins_y, y_bin_edge_l, y_border);
  TGraph *gry_r = new TGraph (num_bins_y, y_bin_edge_r, y_border);
  TGraph *gry = new TGraph (num_bins_y, y_bins, y_border);

  //X position histogram
  TCanvas *c3 = new TCanvas("C3","",1200,500);
  xpos->SetLineColor(2);
  xpos->SetLineWidth(1);
  xpos->Draw();
  gStyle->SetBarWidth(0.02);
  grx_l->SetFillColor(3);
  grx_l->Draw("B");
  grx_r->SetFillColor(3);
  grx_r->Draw("B");
  grx->SetFillColorAlpha(6,0.1);
  grx->Draw("B");
  c3->SaveAs("mts_xpos.png");

  //Y position histogram
  TCanvas *c4 = new TCanvas("C4","",1200,500);
  ypos->SetLineColor(2);
  ypos->SetLineWidth(1);
  ypos->Draw();
  gStyle->SetBarWidth(0.02);
  gry_l->SetFillColor(3);
  gry_l->Draw("B");
  gry_r->SetFillColor(3);
  gry_r->Draw("B");
  gry->SetFillColorAlpha(6,0.1);
  gry->Draw("B");
  c4->SaveAs("mts_ypos.png");
  
  return 0;

}//DONE MAIN

