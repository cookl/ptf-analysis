#include "WaveformFitResult.hpp"
#include "ScanPoint.hpp"
#include "TFile.h"
#include "TCanvas.h"
#include "TH1.h"
#include "TAxis.h"
#include "TPad.h"
#include "TGraph.h"

#include <iostream>
#include <vector>
#include <fstream>


int main( int argc, char* argv[] ) {

  if ( argc != 3 ){
    std::cerr<<"Usage: ptf_ttree_analysis.app ptf_analysis.root run_number\n";
    exit(0);
  }

  // creating variables for afterpulse:  
  int lsrTotal = 0, afpTotal = 0, singleAfpCount = 0;
  int afpTimeThreshold;
  double wfValue, afpValue;
  bool afpCounted;

  // creating histogram file names:
  char fnmT[1024], fnmQ[1024], fnmQlsr[1024], fnmQafp[1024];
  sprintf(fnmT, "../mpmt_time-%s.png", argv[2]);
  sprintf(fnmQ, "../mpmt_charge-%s.png", argv[2]);
  sprintf(fnmQlsr, "../mpmt_charge-lsr-%s.png", argv[2]);
  sprintf(fnmQafp, "../mpmt_charge-afp-%s.png", argv[2]);

  // creating output file:
  char fnmOut[1024];
  sprintf(fnmOut, "../mpmt-afp-%s.txt", argv[2]);
  std::ofstream outFile;
  
  // arbritary values of a threshold and a time window
  afpTimeThreshold   = 2300;

  // opening the root file
  TFile * fin = new TFile( argv[1], "read" );

  // adding a canvas and a histogram
  auto canvas1 = new TCanvas("canvas1","",800,500);
  auto canvas2 = new TCanvas("canvas2","",800,500);
  auto canvas3 = new TCanvas("canvas3","",800,500);
  auto canvas4 = new TCanvas("canvas4","",800,500);
  auto canvas5 = new TCanvas("canvas5","",800,500);
  auto histT    = new TH1F("histT","",1000,0,8200*2.0*0.48828125);
  auto histQ    = new TH1F("histQ","",500,0,1000*1.5*0.48828125);
  auto histQlsr = new TH1F("histQlsr","",500,0.0,1000*1.5*0.48828125);
  auto histQafp = new TH1F("histQafp","",500,0.0,1000*1.5*0.48828125);

  // get the waveform fit TTree
  TTree * tt = (TTree*)fin->Get("ptfanalysis0");
  WaveformFitResult * wf = new WaveformFitResult;
  wf->SetBranchAddresses( tt );
  
  std::vector<double> lsrPulses;

  std::cout << "\n" << std::endl;
  std::cout << "LOOPING TREE, 1:" << std::endl;
  std::cout << "Performing global calculations." <<std::endl;
  for(int i = 0; i < tt->GetEntries(); i++){
    tt->GetEvent(i );
    // std::cout<<i<< ". number of found pulses="<<wf->numPulses << std::endl;
    if(wf->numPulses > 0){
      afpCounted = false;
      for(int i = 0; i < wf->numPulses; i++){
        // filling up the histogram
        histT->Fill(wf->pulseTimes[i]);
        histQ->Fill(wf->pulseCharges[i]*1000.0);

        // check if it is laser pulse. 
        // how many p.e. </to be implemented
        if (wf->pulseTimes[i]<afpTimeThreshold){
          histQlsr->Fill(wf->pulseCharges[i]*1000.0);
          lsrPulses.push_back(wf->pulseCharges[i]*1000.0);
          lsrTotal++;
        } else if (wf->pulseTimes[i]>=afpTimeThreshold){
          afpTotal++;
          histQafp->Fill(wf->pulseCharges[i]*1000.0);
          if (afpCounted == false) {
            singleAfpCount++;
            afpCounted = true;
          }
        }
      }             
    }
  }

  // fixing the pe height
  double peHeight = 9.52; //< mV 
  std::cout<<"pe value used: " << peHeight << std::endl;

  // finding maximum of laser pulses:
  double lsrMax = lsrPulses[0];
  for (int i=1; i<lsrPulses.size(); ++i){
    if (lsrMax < lsrPulses[i]) {
      lsrMax = lsrPulses[i];
    }
  }
  std::cout<<"The maximum pulse height is (mV):"<<lsrMax<<std::endl;

  double peDivision = 0.5; //< we want 0.5 pe resolution
  int N = int(ceil(lsrMax/peHeight/peDivision)); //< determining how many 0.5 pe regions we need

  // finding the most populated bin
  int maxY = histQlsr->GetBinContent(histQlsr->GetMaximumBin());
  double maxX2 = histQlsr->GetXaxis()->GetBinCenter(histQlsr->GetMaximumBin());
  std::cout<<"\nSUMMARY OF GLOBAL HISTOGRAMS"<<std::endl;
  std::cout<<"Population of highest bin: "<<maxY<<std::endl;
  std::cout<<"Position of highest bin (mV): "<<maxX2<<std::endl;

  // calculating global afterpulse rates
  std::cout<<"\nSUMMARY OF TOTALS"<<std::endl;
  std::cout<<"Lsr pulses: "<< afpTotal << std::endl;
  std::cout<<"Afp: "<< lsrTotal << std::endl;
  std::cout<<"Waveforms w/ afp: " <<singleAfpCount << std::endl;
  std::cout<<"Afp by lsr pulses: "<< double(afpTotal)/lsrTotal*100. << std::endl;
  std::cout<<"Waveforms w/ afp by lsr pulses:"<< double(singleAfpCount)/lsrTotal*100<<std::endl;


  // i always starts at 0
  // j always starts at 1

  // defining all pe regions we are looking for
  // and initializing the arrays
  double peRegions[N], afpRates[N];
  bool autoLsrHappened[N];
  int autoLsrCounter[N], autoAfpCounter[N];
  for (int j=1; j<=N; j++){
    peRegions[j] = peDivision*j;
    autoLsrHappened[j] = false;
    autoLsrCounter[j] = 0;
    autoAfpCounter[j] = 0;
  }

  std::cout << "\n" << std::endl;
  std::cout << "LOOPING TREE, 2:" << std::endl;
  std::cout << "Evaluating afterpulse rates by pe." <<std::endl;
  for(int i = 0; i < tt->GetEntries(); i++){
    tt->GetEvent(i );
    // std::cout<<i<< ". number of found pulses="<<wf->numPulses << std::endl;
    if(wf->numPulses > 0){
      

      for(int i = 0; i < wf->numPulses; i++){
        if (wf->pulseTimes[i]<afpTimeThreshold){
          
          // setting all lsr booleans to false at the start of every loop.
          for (int j=1; j<=N; j++) {
            autoLsrHappened[j] = false;
          }

          // transform charge value to p.e. height
          wfValue = wf->pulseCharges[i]*1000.0/peHeight;
          
          // classifying pulse heights as multiple of 0.5 p.e.
          // counting them 
          for (int j=1; j<=N; j++) {
            if (j==1){
              if (wfValue<=peRegions[1]) {
                autoLsrHappened[1] = true;
                autoLsrCounter[1] += 1;
              }
            } else if (j>1 && j<N) {
              if (wfValue>peRegions[j] && wfValue<=peRegions[j+1]){
                autoLsrHappened[j] = true;
                autoLsrCounter[j] += 1;
              }
            } else if (j==N) {
              if (wfValue>peRegions[N]) {
                autoLsrHappened[N] = true;
                autoLsrCounter[N] += 1;
              }
            }
          } 
        
          // counting the afterpulses generated by
          // each p.e. height
        } else if (wf->pulseTimes[i]>=afpTimeThreshold){
          afpValue = wf->pulseCharges[i];
          for (int j=1;j<=N;j++) {
            if (autoLsrHappened[j] == true) {
              autoAfpCounter[j] += 1;
              autoLsrHappened[j] = false;
            }
          }
        }
      }             
    }
  }

  // evaluating afterpulse rate.
  // as afp per lsr, and afp per total if necessary.
  for (int j=1; j<=N; j++){
    if (autoLsrCounter[j] != 0){
      afpRates[j] = double(autoAfpCounter[j])/double(autoLsrCounter[j]);
      std::cout<<peRegions[j]<<" p.e.: "<< afpRates[j]*100.0 << "%"<<std::endl; 
    }
  }


  // outFile.open(fnmOut);
  // for (int i=0; i<=15; i++) {
  //   // char line[1024];
  //   // sprintf(line,"             %f             %f",pe[i],afpAverages[i]);
  //   // outFile << line << std::endl;
  //   outFile << pe[i] << std::setw(18) << afpAverages[i] << std::endl;
  //   // std::cout<<pe[i]<<", "<<afpAverages[i]<<std::endl;
  // }


  // TIME histogram
  canvas1->cd();
  histT->SetTitle("Pulse times");
  histT->GetXaxis()->SetTitle("Time (ns)");
  histT->GetYaxis()->SetTitle("No. of events");
  histT->SetMaximum(1500);
  histT->Draw();
  canvas1->SaveAs(fnmT);

  // ALL PULSES histogram
  canvas2->cd();
  histQ->SetTitle("Total pulses");
  histQ->GetXaxis()->SetTitle("Pulse height (mV)");
  histQ->GetYaxis()->SetTitle("No. of events");
  histQ->GetXaxis()->SetRange(0,100);
  histQ->Draw();
  canvas2->SaveAs(fnmQ);

  // LASER PULSES histogram
  canvas3->cd();
  histQlsr->SetTitle("Laser pulses");
  histQlsr->GetXaxis()->SetTitle("Laser pulse height (mV)");
  histQlsr->GetYaxis()->SetTitle("No. of events");
  histQlsr->GetXaxis()->SetRange(0,100);
  histQlsr->Draw();
  canvas3->SaveAs(fnmQlsr);  

  // AFTER PULSES histogram
  canvas4->cd();
  histQafp->SetTitle("Afterpulses");
  histQafp->GetXaxis()->SetTitle("Afterpulse height (mV)");
  histQafp->GetYaxis()->SetTitle("No. of events");
  histQafp->GetXaxis()->SetRange(0,100);
  histQafp->Draw();
  canvas4->SaveAs(fnmQafp);

  // ALL-IN-ONE histogram
  // </Change the colors in the future.
  canvas5->cd();
  histQafp->SetStats(0);
  histQlsr->SetStats(0);
  histQlsr->SetTitle("Split pulses");
  histQlsr->GetXaxis()->SetTitle("Pulse height (mV)");
  histQlsr->GetYaxis()->SetTitle("No. of events");
  histQlsr->SetLineColor(3);
  histQlsr->SetMaximum(68000);
  histQlsr->Draw();
  histQafp->SetLineColor(4);
  histQafp->Draw("same");
  histQ->Draw("same");
  canvas5->SaveAs("../both.png");

  return 0;
}
