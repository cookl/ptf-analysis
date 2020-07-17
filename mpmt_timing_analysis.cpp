#include "WaveformFitResult.hpp"
#include "ScanPoint.hpp"
#include "TFile.h"
#include "TCanvas.h"
#include "TH1D.h"
#include "TGraph.h"
#include "TF1.h"
#include "TStyle.h"

#include <iostream>
#include <vector>


int main( int argc, char* argv[] ) {

  if ( argc != 2 ){
    std::cerr<<"Usage: ptf_ttree_analysis.app ptf_analysis.root\n";
    exit(0);
  }

  TFile * fin = new TFile( argv[1], "read" );


  // get the waveform fit TTree
  TTree * tt0 = (TTree*)fin->Get("ptfanalysis0");
  WaveformFitResult * wf0 = new WaveformFitResult;
  wf0->SetBranchAddresses( tt0 );
  TTree * tt1 = (TTree*)fin->Get("ptfanalysis1");
  WaveformFitResult * wf1 = new WaveformFitResult;
  wf1->SetBranchAddresses( tt1 );


  TH1F *tdiffs[15];
  for(int i =0; i < 15; i++){
    char name[100];
    char title[100];
    sprintf(name,"timediff_%i_pe",i+1);
    sprintf(title,"Time diff %i PE",i+1);
    tdiffs[i] = new TH1F(name,title,100,-5,5);
  }
  

  TH1F *tdiff = new TH1F("time diff","timediff",100,-5,5);
  TH1F *tdiff2 = new TH1F("time diff2","timediff2",100,-5,5);
  std::cout << "Looping tree " << tt0->GetEntries() << " " << tt1->GetEntries() << std::endl;
  // Loop the first scan point and print something
  for(int i = 0; i < tt0->GetEntries(); i++){
  //for(int i = 0; i < 50000; i++){
    tt0->GetEvent(i );
    tt1->GetEvent(i );

    // Find the pulses in list of pulses
    double pulse_time[2] = {-1,-1};
    double pulse_height[2] {-1,-1};

    double baseline[2] = {0.9985,1.004};
    for(int j =0; j < 2; j++){
      WaveformFitResult *wf;
      if(j==0) wf = wf0;
      if(j==1) wf = wf1;
      for(int k = 0; k < wf->numPulses; k++){
        if(wf->pulseTimes[k] > 2150 && wf->pulseTimes[k] < 2250){ // look for laser pulse
          pulse_time[j] = wf->pulseTimes[k];
          pulse_height[j] = (baseline[j] - wf->pulseCharges[k])/0.01;          
        }
      }
    }
    
   
    double time0 = wf0->mean;
    double time1 = wf1->mean;
    double cfd_time0 = wf0->sinw;
    double cfd_time1 = wf1->sinw;
    double time_diff = time0 - time1;

    for(int j = 0; j < 15; j++){

      
      double min_pe = (j+1.0) - 0.5;
      double max_pe = (j+1.0) + 0.5;
      if(pulse_height[0] < max_pe && pulse_height[1] < max_pe &&
         pulse_height[0] > min_pe && pulse_height[1] > min_pe){

        tdiffs[j]->Fill(time_diff);

        if(j == 0 and time_diff < -3 && i < 50000){
          std::cout << i << " Time difference: " << time0 << " - " << time1
                    << " : " << time0-time1 << std::endl;
          std::cout << "Charge: " << wf0->amp << " " << wf1->amp << " "
                    << pulse_time[0] << " " << pulse_time[1] << " "
                    << pulse_height[0] << " " << pulse_height[1] << " "
                    << std::endl;
          std::cout << "Ratio "
                    << wf0->amp/pulse_height[0] << " " 
                    << wf1->amp/pulse_height[1] << " " << std::endl;
        }
      }
    }
    

    // Only look at large pulses
    if(pulse_height[0] < 5.0 && pulse_height[1] < 5.0 &&
       pulse_height[0] > 3.0 && pulse_height[1] > 3.0       ){
      tdiff->Fill(time0-time1);
      tdiff2->Fill(cfd_time0-cfd_time1);

      if(i < 10000 && 0){
        std::cout << i << " Time difference: " << time0 << " - " << time1
                  << " : " << time0-time1 << std::endl;
        std::cout << "Charge: " << wf0->amp << " " << wf1->amp << " "
                  << pulse_time[0] << " " << pulse_time[1] << " "
                  << pulse_height[0] << " " << pulse_height[1] << " "
                  << std::endl;
        std::cout << "Ratio "
                  << wf0->amp/pulse_height[0] << " " 
                  << wf1->amp/pulse_height[1] << " " << std::endl;
        
      }
    }
     
  }

  TCanvas *c = new TCanvas("C");
  tdiff2->Draw();
  TF1 *gaus = new TF1("gaus","gaus",-1,3);
  tdiff2->Fit("gaus","R");


  gStyle->SetOptFit(1111);
  tdiff->Draw("SAME");
  tdiff->Fit("gaus","R");

  c->SaveAs("tdiff.png");

  TGraph *gr = new TGraph();
  
  
    
  for(int i = 0; i < 15; i++){

    char name[100];
    sprintf(name,"C%i",i);
    TCanvas *c = new TCanvas(name);
    tdiffs[i]->Draw();
    tdiffs[i]->Fit("gaus","R");
    gStyle->SetOptFit(1111);
    sprintf(name,"tdiff_%i.png",i+1);
    int x = i+1;
    gr->SetPoint(i,x,gaus->GetParameter(2));
    c->SaveAs(name);
    
    
  }

  TCanvas *c2 = new TCanvas("CC");
  gr->Draw("AP*");
  TF1 *sqrtf = new TF1("sqrtf","sqrt([0]/x + [1])",0.5,8);
  sqrtf->SetParameter(0,1);
  sqrtf->SetParameter(1,0.5);
  
  gr->Fit("sqrtf","R");
  gStyle->SetOptFit(1111);
  c2->SaveAs("tres_vs_ph.png");
  
  return 0;
}

