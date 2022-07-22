// Ashley/Thomas/Angela 2D histograms
//Modified by Mohit Gola
// 2022-07-04
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
#include <math.h> 
//using namespace std;


//Functions for drawing PMTs and plotting;
/*---------------------------------------------------------------------------------*/

void drawPMT(Double_t x, Double_t y){
  TEllipse pmt(x,y,0.03,0.03);
  pmt.SetLineColor(1);
  pmt.SetFillStyle(0);
  pmt.SetLineWidth(2);
  pmt.DrawClone();
  
  }

void plot2D( TCanvas *c,bool plotting_on,TH2F *hist,char *title,const char *zt,char *name){

  if(plotting_on == true){
    c->SetRightMargin(0.2);
    c->SetLeftMargin(0.15);
    c->SetBottomMargin(0.10);
    hist->GetXaxis()->SetTitle("X (m)");
    //    hist->GetXaxis()->SetRangeUser(0.,1.);
    //    hist->GetYaxis()->SetRangeUser(0.,1.);
    hist->GetYaxis()->SetTitle("Y (m)");
    hist->GetZaxis()->SetTitle(zt);
    hist->SetTitle(title);
    hist->Draw("COLZ");
    hist->SetTitleSize(50,"t");
    c->SetRealAspectRatio();
    c->Modified();
    c->Update();
    drawPMT(0.29,0.125); //PMT2 Location
    drawPMT(0.305,0.215); //PMT3 Location
    c->SaveAs(name);}
  
}

void plot1D( TCanvas *c,bool plotting_on,TH1F *hist,char *title,char *name){

  if(plotting_on == true){
    hist->Draw();
    hist->GetXaxis()->SetTitle("Pulse height (mV)");
    hist->GetYaxis()->SetTitle("Number of events");
    hist->SetTitle(title);
    gStyle->SetOptFit(11);
    c->SaveAs(name);}
  
}


int main( int argc, char* argv[] )
{//START MAIN 

   //Specify desired plots
  /*---------------------------------------------------------------------------------*/
  bool plot_mph = true;
  bool plot_p = true;
  bool plot_RMS = false;
  bool plot_scan_pt = false;
  bool plot_events = true;
  bool plot_eff = true;
  bool plot_pP = true;
  bool plot_Pw = true;
  bool plot_h = true;
 
  //Define paramaters for scan. (m)
  /*---------------------------------------------------------------------------------*/
  float step_size = 0.01;
  float xstart = 0.15;
  float ystart = 0.05;
  float x_scan_dist = 0.32;
  float y_scan_dist = 0.3;
  float dwelltime = 3500;
  int num_ch = 2; //number of active channels
  int f_ch = 2; //first channel
  
  int num_bins_x = static_cast<int>(x_scan_dist/step_size)+1 ;
  //  int num_bins_x = 6;  
  //  int num_bins_y = 11;
  int num_bins_y = static_cast<int>(y_scan_dist/step_size)+1 ;

  float x_low = xstart - step_size/2;
  //float x_low = 0.;
  float x_high = xstart + x_scan_dist + step_size/2;
  //float x_high = 1.;
  float y_low = ystart - step_size/2;;
  //float y_low = 0.;
  float y_high = ystart + y_scan_dist + step_size/2;
  //float y_high = 1.;
  
  std::cout << "There are approx" << num_bins_x * num_bins_y << "scan points" << std::endl;
  std::cout << "We expect" << dwelltime * num_bins_x * num_bins_y << "events in total" << std::endl;

  //Specify the coordinate you want to look at for the 1D histogram;
  /*---------------------------------------------------------------------------------*/
  //  float POI_y = 0.1; //scan point of interest y
  //  float POI_x = 0.1;
  //scan point of interest x
  float POI_y = 0.385; //scan point of interest y
  float POI_x = 0.37; //scan point of interest x
  //int BOI_x = 0.5;
  //  int BOI_y = 0.5;
  int BOI_x = static_cast<int>((POI_x-xstart)/step_size);
  int BOI_y = static_cast<int>((POI_y-ystart)/step_size);

  //  TH1::AddDirectory(false);
  //  TH2::AddDirectory(false);
  //Initialize histograms and bins of combined channels
  /*---------------------------------------------------------------------------------*/
  TH2F *h_eff_sum = new TH2F("Total Efficiency","LED pulses/waveform",num_bins_x,x_low,x_high,num_bins_y,y_low,y_high);
  TH2F *h_pP_sum = new TH2F(" LED pulses / Raw Pulse Count","",num_bins_x,x_low,x_high,num_bins_y,y_low,y_high);
  TH2F *h_Pw_sum = new TH2F("Pulse / waveforms","",num_bins_x,x_low,x_high,num_bins_y,y_low,y_high) ;
  TH2F *h_p_sum = new TH2F("pulses","",num_bins_x,x_low,x_high,num_bins_y,y_low,y_high) ;
  TH2F *h_mph_sum = new TH2F("mean pulse height","",num_bins_x,x_low,x_high,num_bins_y,y_low,y_high) ;
  TH2F *h_eff_bin[num_ch];
  TH2F *h_pP_bin[num_ch];
  TH2F *h_Pw_bin[num_ch];
  TH2F *h_p_bin[num_ch];
  TH2F *h_mph_bin[num_ch];

  
  //Read ROOT file;
  if ( argc != 2 ){
    std::cerr<<"Usage: ptf_ttree_analysis.app ptf_analysis.root\n";
    exit(0); }
  TFile * fin = new TFile( argv[1], "read" );
  TTree * tt0;    
  WaveformFitResult * wf0;
 

 //Start the loops
  /*---------------------------------------------------------------------------------*/
  
  for (int ch=0;ch<num_ch;ch++)
    {//START CHANNEL LOOP
      int ch_name = ch + f_ch;
      std::cout << "CHANNEL " << ch_name << std::endl;
      //      if (ch_name == 2) continue;
      //      if (ch_name == 4) continue;
      
      /*      
      if (ch_name == 4) continue;
      if (ch_name == 5) continue;
      if (ch_name == 7) continue;
      if (ch_name == 8) continue;
      if (ch_name == 9) continue;
      if (ch_name == 10) continue;
      if (ch_name == 11) continue;
      if (ch_name == 13) continue;
      if (ch_name == 14) continue;
      */
      
  
      // Getting the waveform fit TTree for desired channel;
      char channel[50];
      sprintf(channel,"ptfanalysis%i",ch_name);
      tt0 = (TTree*)fin->Get(channel);
      wf0 = new WaveformFitResult;
      if(tt0) wf0->SetBranchAddresses( tt0 );
  
      /*----------------------------------------------------------------------------------*/
      //Initialize
      /* --------------------------------------------------------------------------------- */
      //2D Histograms for separate channels; 
      TH2F *h_mph = new TH2F("mean pulse height distribution","Mean Pulse Height",num_bins_x,x_low,x_high,num_bins_y,y_low,y_high);
      TH2F *h_p = new TH2F("Number of Pulses","Number of Pulses",num_bins_x,x_low,x_high,num_bins_y,y_low,y_high);
      //      TH2F *h_p = new TH2F("Number of Pulses","Number of Pulses",11,0.,1.,11,0.,1.);
      TH2F *h_RMS = new TH2F("RMS", "Standard dev. dist.",num_bins_x,x_low,x_high,num_bins_y,y_low,y_high);
      TH2F *h_scan_pt = new TH2F("scan point","scan point",num_bins_x,x_low,x_high,num_bins_y,y_low,y_high);
      TH2F *h_events = new TH2F("Events","Number of events",num_bins_x,x_low,x_high,num_bins_y,y_low,y_high);
      TH2F *h_P = new TH2F("Raw Pulse Count","Number of Pulses",num_bins_x,x_low,x_high,num_bins_y,y_low,y_high);
      TH2F *h_eff = new TH2F("Efficiency","laser pulses/waveform",num_bins_x,x_low,x_high,num_bins_y,y_low,y_high);
      //TH2F *h_eff = new TH2F("Efficiency","laser pulses/waveform",11,0.,1.,11,0.,1.);
      TH2F *h_pP = new TH2F("laser pulse / Raw Pulse Count","",num_bins_x,x_low,x_high,num_bins_y,y_low,y_high);
      TH2F *h_Pw = new TH2F("Pulse / waveforms","",num_bins_x,x_low,x_high,num_bins_y,y_low,y_high) ;
  

  
      //2D array bins for constructing the 2D histograms;
      TH1F *h[num_bins_x][num_bins_y];
      Double_t events_bin[num_bins_x][num_bins_y];
      Double_t pulse_bin[num_bins_x][num_bins_y];
      
      for (int x=0;x<num_bins_x;x++) {
	for (int y=0;y<num_bins_y;y++){
	  char name1[100];
	  sprintf(name1,"h bin %i %i %i",x,y,ch);
	  h[x][y] = new TH1F(name1,"pulse heights",200,0,200*0.48828125);
	  events_bin[x][y]=0;
	  pulse_bin[x][y]=0;
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
	      //	      float ycenter = 0.38;
	      float y_l = ycenter - step_size/2.0;
	      float y_r = ycenter + step_size/2.0;
											  
	      for(int xpoint = 0; xpoint < num_bins_x ; xpoint++)
		{//START XLOOP
		  float xcenter = step_size*xpoint + xstart;
		  //		  float xcenter = 0.46;
		  float x_l = xcenter - step_size/2.0;
		  float x_r = xcenter + step_size/2.0;

		  if((wf0->x>x_l && wf0->x < x_r) and (wf0->y >y_l && wf0->y < y_r))
		    {//FILTER XY COORDS
		      
		      events_bin[xpoint][ypoint] += 1;
		      pulse_bin[xpoint][ypoint] += wf0->numPulses;

		      for(int k = 0; k < wf0->numPulses; k++ )
			{//START PULSE LOOP

			  //if(wf0->pulseTimes[k] > 2300 and wf0->pulseTimes[k] < 2420 and wf0->pulseCharges[k]*1000.0 > 2.0)
			  if(wf0->pulseTimes[k] > 2300 and wf0->pulseTimes[k] < 2420 and wf0->pulseCharges[k]*1000.0 > 2.0 )
			    {//FILTER PULSE TIME

			      //			      std::cout << "Number of pulses found: " << wf0->pulseCharges[k]*1000 << std::endl;
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
	  h_P->SetBinContent(x+1,y+1,pulse_bin[x][y]);
	  //	  std::cout << "Mean Values: " << h[x][y]->GetMean() << std::endl; 
	}}

      h_p_bin[ch] = h_p;
      h_p_sum->Add(h_p_bin[ch]);

      h_mph_bin[ch] = h_mph;
      h_mph_sum->Add(h_mph_bin[ch]);
      
      h_eff = (TH2F*)h_p->Clone();
      h_eff->Divide(h_events);
      h_eff->Scale(100.0);
      h_eff_bin[ch] = h_eff;
      h_eff_sum->Add(h_eff_bin[ch]);

      h_Pw = (TH2F*)h_P->Clone();
      h_Pw->Divide(h_events);
      h_Pw->Scale(100.0);
      h_Pw_bin[ch]=h_Pw;
      h_Pw_sum->Add(h_Pw_bin[ch]);

      h_pP = (TH2F*)h_p->Clone();
      h_pP->Divide(h_P);
      h_pP->Scale(100.0);
      h_pP_bin[ch]=h_pP;
      h_pP_sum->Add(h_pP_bin[ch]);
      /* --------------------------------------------------------------------------------- */

  
      //Plotting
      /* --------------------------------------------------------------------------------- */

      //      TCanvas *c0 = new TCanvas("C0");
      //      char *title0 = Form("pulse height distribution channel %d (%f,%f) bin (%d,%d)",ch_name,POI_x,POI_y,BOI_x,BOI_y);
      //      char *name0 = Form("pulse_height_distribution%d.png",ch_name);
      //      plot1D(c0,plot_h,h[BOI_x][BOI_y],title0, name0);

      //      std::cout<<"check"<<std::endl; 
	    

      
      //Style;
      gStyle->SetPalette(1);
      gStyle->SetOptTitle(1); gStyle->SetOptStat(0);
      gStyle->SetOptFit(1111); gStyle->SetStatBorderSize(0);
      gStyle->SetStatX(.89); gStyle->SetStatY(.89);	      
  
      TCanvas *c1 = new TCanvas("C1","",1200,500);
      char *title1 = Form("Scan points channel %d",ch_name);
      char *name1 = Form("scan_points_%d.png",ch_name);
      plot2D(c1,plot_scan_pt,h_scan_pt,title1,"scan point", name1);
 
      TCanvas *c2 = new TCanvas("C2","",1200,500);
      char *title2 = Form("Pulse counts channel %d",ch_name);
      char *name2 = Form("2D_p_%d.png",ch_name);
      plot2D(c2,plot_p,h_p,title2,"pulses", name2);
 
      TCanvas *c3 = new TCanvas("C3","",1200,500);
      char *title3 = Form("Mean Pulse height (mV) channel %d",ch_name);
      char *name3 = Form("2D_mph_%d.png",ch_name);
      plot2D(c3,plot_mph,h_mph,title3,"mV", name3);

      TCanvas *c4 = new TCanvas("C4","",1200,500);
      char *title4 = Form("RMS (mV) channel %d",ch_name);
      char *name4 = Form("2D_RMS_%d.png",ch_name);
      plot2D(c4,plot_RMS,h_RMS,title4,"mV", name4);

      TCanvas *c5 = new TCanvas("C5","",1200,500);
      char *title5 = Form("Efficiency (laser pulses / waveforms) channel %d",ch_name);
      char *name5 = Form("2D_Efficiency_%d.png",ch_name);
      plot2D(c5,plot_eff,h_eff,title5,"%", name5);

      TCanvas *c6 = new TCanvas("C6","",1200,500);
      char *title6 = Form("Number of waveforms channel %d",ch_name);
      char *name6 = Form("2D_Events_%d.png",ch_name);
      plot2D(c6,plot_events,h_events,title6,"events", name6);

      TCanvas *c7 = new TCanvas("C7","",1200,500);
      char *title7 = Form("laser pulses / Pulses channel %d",ch_name);
      char *name7 = Form("2D_pulse.Pulse_%d.png",ch_name);
      plot2D(c7,plot_pP,h_pP,title7,"%", name7);

      TCanvas *c8 = new TCanvas("C8","",1200,500);
      char *title8 = Form("Pulses / waveforms channel %d",ch_name);
      char *name8 = Form("2D_Pulse.waveform_%d.png",ch_name);
      plot2D(c8,plot_Pw,h_Pw,title8,"%", name8);		

      TCanvas *c14 = new TCanvas("C14","",1200,500);
      char *title14 = Form("Pulses / waveforms channel %d",ch_name);
      char *name14 = Form("2D_Pulse.waveform_%d.png",ch_name);
      plot2D(c14,plot_Pw,h_Pw,title14,"%", name14);

      TCanvas *c15 = new TCanvas("C15","",1200,500);
      char *title15 = Form("Pulses / waveforms channel %d",ch_name);
      char *name15 = Form("2D_Pulse.waveform_%d.png",ch_name);
      plot2D(c15,plot_Pw,h_Pw,title15,"%", name15);

      TCanvas *c16 = new TCanvas("C16","",1200,500);
      char *title16 = Form("Pulses / waveforms channel %d",ch_name);
      char *name16 = Form("2D_Pulse.waveform_%d.png",ch_name);
      plot2D(c16,plot_Pw,h_Pw,title16,"%", name16);

      TCanvas *c17 = new TCanvas("C17","",1200,500);
      char *title17 = Form("Pulses / waveforms channel %d",ch_name);
      char *name17 = Form("2D_Pulse.waveform_%d.png",ch_name);
      plot2D(c17,plot_Pw,h_Pw,title17,"%", name17);

      TCanvas *c18 = new TCanvas("C18","",1200,500);
      char *title18 = Form("Pulses / waveforms channel %d",ch_name);
      char *name18 = Form("2D_Pulse.waveform_%d.png",ch_name);
      plot2D(c18,plot_Pw,h_Pw,title18,"%", name18);

      TCanvas *c19 = new TCanvas("C19","",1200,500);
      char *title19 = Form("Pulses / waveforms channel %d",ch_name);
      char *name19 = Form("2D_Pulse.waveform_%d.png",ch_name);
      plot2D(c19,plot_Pw,h_Pw,title19,"%", name19);

      TCanvas *c20 = new TCanvas("C20","",1200,500);
      char *title20 = Form("Pulses / waveforms channel %d",ch_name);
      char *name20 = Form("2D_Pulse.waveform_%d.png",ch_name);
      plot2D(c20,plot_Pw,h_Pw,title20,"%", name20);

      TCanvas *c21 = new TCanvas("C21","",1200,500);
      char *title21 = Form("Pulses / waveforms channel %d",ch_name);
      char *name21 = Form("2D_Pulse.waveform_%d.png",ch_name);
      plot2D(c21,plot_Pw,h_Pw,title21,"%", name21);

      TCanvas *c22 = new TCanvas("C22","",1200,500);
      char *title22 = Form("Pulses / waveforms channel %d",ch_name);
      char *name22 = Form("2D_Pulse.waveform_%d.png",ch_name);
      plot2D(c22,plot_Pw,h_Pw,title22,"%", name22);

      TCanvas *c23 = new TCanvas("C23","",1200,500);
      char *title23 = Form("Pulses / waveforms channel %d",ch_name);
      char *name23 = Form("2D_Pulse.waveform_%d.png",ch_name);
      plot2D(c23,plot_Pw,h_Pw,title23,"%", name23);

      TCanvas *c24 = new TCanvas("C24","",1200,500);
      char *title24 = Form("Pulses / waveforms channel %d",ch_name);
      char *name24 = Form("2D_Pulse.waveform_%d.png",ch_name);
      plot2D(c24,plot_Pw,h_Pw,title24,"%", name24);
      
    }//END CHANNEL LOOP

  /*Plotting over sum of all channels*/
  TCanvas *c9 = new TCanvas("C9","",1200,500);
  char *title9 = Form("Efficiency sum over %d channels",num_ch);
  char *name9 = Form("2D_SUM_Efficiency_%d.png",num_ch);
  plot2D(c9,true,h_eff_sum,title9,"%",name9);

  TCanvas *c10 = new TCanvas("C10","",1200,500);
  char *title10 = Form("LED pulse / total (sum over %d channels)",num_ch);
  char *name10 = Form("2D_SUM_pP_%d.png",num_ch);
  plot2D(c10,true,h_pP_sum,title10,"%",name10);

  TCanvas *c11 = new TCanvas("C11","",1200,500);
  char *title11 = Form("Pulses / waveforms (sum over %d channels)",num_ch);
  char *name11 = Form("2D_SUM_Pw_%d.png",num_ch);
  plot2D(c11,true,h_Pw_sum,title11,"%",name11);

  TCanvas *c12 = new TCanvas("C12","",1200,500);
  char *title12 = Form("number of pulses (sum over %d channels)",num_ch);
  char *name12 = Form("2D_SUM_p_%d.png",num_ch);
  plot2D(c12,true,h_p_sum,title12,"%",name12);

  TCanvas *c13 = new TCanvas("C13","",1200,500);
  char *title13 = Form("mean pulse height (sum over %d channels)",num_ch);
  char *name13 = Form("2D_SUM_mph_%d.png",num_ch);
  plot2D(c13,true,h_mph_sum,title13,"%",name13);

  return 0;

}//DONE MAIN

