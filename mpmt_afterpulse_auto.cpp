#include "WaveformFitResult.hpp"
#include "ScanPoint.hpp"
#include "TFile.h"
#include "TCanvas.h"
#include "TH1.h"
#include "TF1.h"
#include "TPad.h"
#include "TLegend.h"
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
  double wfValue;
  bool afpCounted;

  // creating histogram file names:
  char fnmT[1024], fnmQ[1024], fnmQlsr[1024], fnmQafp[1024], fnmAfpRate[1024], fnmSingleAfpRate[1024], combinedPlot[1024];
  sprintf(fnmT, "../%s-time-hist.png", argv[2]);
  sprintf(fnmQ, "../%s-charge-hist.png", argv[2]);
  sprintf(fnmQlsr, "../%s-lsr-hist.png", argv[2]);
  sprintf(fnmQafp, "../%s-afp-hist.png", argv[2]);
  sprintf(fnmAfpRate, "../%s-afpRate.png", argv[2]);
  sprintf(fnmSingleAfpRate, "../%s-single-afpRate.png", argv[2]);
  sprintf(combinedPlot, "../%s-combined-afpRate.png", argv[2]);

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
  auto canvas6 = new TCanvas("canvas6","",1600,1000);
  auto canvas7 = new TCanvas("canvas7","",1600,1000);
  auto canvas8 = new TCanvas("canvas8","",1600,1000);
  auto histT    = new TH1F("histT","",1000,0,8200*2.0*0.48828125);
  auto histQ    = new TH1F("histQ","",500,0,1000*1.5*0.48828125);
  auto histQlsr = new TH1F("histQlsr","",500,0.0,1000*1.5*0.48828125);
  auto histQafp = new TH1F("histQafp","",500,0.0,1000*1.5*0.48828125);
  auto afpGraph = new TGraph();
  auto afpSingleGraph = new TGraph();

  // get the waveform fit TTree
  TTree * tt = (TTree*)fin->Get("ptfanalysis1");
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
  for (int i=1; i<int(lsrPulses.size()); ++i){
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
  double peRegions[N], afpRates[N], singleEventAfpRates[N];
  bool autoLsrHappened[N];
  int autoLsrCounter[N], autoAfpCounter[N], autoSingleAfpEventConter[N];
  for (int j=1; j<=N; j++){
    peRegions[j] = peDivision*j;
    autoLsrHappened[j] = false;
    autoLsrCounter[j] = 0;
    autoAfpCounter[j] = 0;
    autoSingleAfpEventConter[j] = 0;
  }

  std::cout << "\n" << std::endl;
  std::cout << "LOOPING TREE, 2:" << std::endl;
  std::cout << "Evaluating afterpulse rates by pe." <<std::endl;
  for(int i = 0; i < tt->GetEntries(); i++){
    tt->GetEvent(i );
    // std::cout<<i<< ". number of found pulses="<<wf->numPulses << std::endl;
    if(wf->numPulses > 0){
      
      bool foundAfpEvent = false;
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
          for (int j=1;j<=N;j++) {
            if (autoLsrHappened[j] == true) {

              // counting all p.e. in the waveforms:
              autoAfpCounter[j] += 1;
              //autoLsrHappened[j] = false;

              // counting whenever an afp happens:
              // whenever autoAfpCounter == 1, then it means an afp happened
              // if autoAfpCounter > 1, we do not count the pulse
              if (!foundAfpEvent) {
                autoSingleAfpEventConter[j]+=1;
              }

              foundAfpEvent = true;
            }
          }
        }
      }             
    }
  }

// evaluating afterpulse rate.
// as afp per lsr, and afp per total if necessary.
int nonZeroAfp=0;
for (int j=1; j<=N; j++){
  if (autoLsrCounter[j] != 0 && autoAfpCounter[j] != 0 ){
    ++nonZeroAfp;
    afpRates[j] = double(autoAfpCounter[j])/double(autoLsrCounter[j]);
    singleEventAfpRates[j] = double(autoSingleAfpEventConter[j])/double(autoLsrCounter[j]);
    std::cout<<j<<": "<<peRegions[j]<<" p.e.: "<< afpRates[j]*100.0 << "% and " <<singleEventAfpRates[j]*100.<<"%"<<std::endl; 
  }
}

// Drawing the graphs of afpRate by lsr height
int upToWhatPoint=43; // </ choose this after seeing the graphs a first time
// int upToWhatPoint = nonZeroAfp; // </ choose this if generating the graphs for the first time.
canvas6->cd();
afpGraph->Set(upToWhatPoint);  
for (int i=0; i<upToWhatPoint; i++){
  afpGraph->SetPoint(i, peRegions[i+1], 100.*afpRates[i+1]);
}
afpGraph->SetTitle("Afterpulse rate");
afpGraph->GetXaxis()->SetTitle("Laser pulse height (p.e.)");
afpGraph->GetYaxis()->SetTitle("Percentage \%");
afpGraph->SetMarkerStyle(20);
afpGraph->SetLineColor(0);
afpGraph->Draw();

TF1* polyFitting = new TF1("poly","[0] + [1]*x^[2]",0.0,peRegions[upToWhatPoint]);
polyFitting->SetParameters(1.,1.,0.5);
polyFitting->SetLineColor(11);
polyFitting->SetLineStyle(1);
afpGraph->Fit(polyFitting,"","",0.0,peRegions[upToWhatPoint+1]);

canvas6->SaveAs(fnmAfpRate);

// Drawing the graphs of afp single event rate by lsr height
canvas7->cd();
afpSingleGraph->Set(upToWhatPoint);  
for (int i=0; i<upToWhatPoint; i++){
  afpSingleGraph->SetPoint(i, peRegions[i+1], 100.0*singleEventAfpRates[i+1]);
}
afpSingleGraph->SetTitle("Single afterpulse rate");
afpSingleGraph->GetXaxis()->SetTitle("Laser pulse height (p.e.)");
afpSingleGraph->GetYaxis()->SetTitle("Percentage \%");
//afpSingleGraph->SetMaximum(0.005);
afpSingleGraph->SetLineColor(0);
// afpSingleGraph->GetYaxis()->SetRange(0,2.5);
afpSingleGraph->SetMarkerStyle(20);
afpSingleGraph->Draw();

// TF1* expFitting = new TF1("exp","[0]*e^([1]*x)",1.0,peRegions[upToWhatPoint]);
// expFitting->SetParameters(1.,1.);
// expFitting->SetLineColor(11);
// expFitting->SetLineStyle(1);
// afpSingleGraph->Fit(expFitting,"","",1.0,peRegions[upToWhatPoint+1]);

canvas7->SaveAs(fnmSingleAfpRate);


 canvas8->cd();
afpGraph->SetTitle("Afterpulse rate");
afpGraph->GetXaxis()->SetTitle("Laser pulse height (p.e.)");
afpGraph->GetYaxis()->SetTitle("Percentage \%");
afpGraph->SetMaximum(180);
afpGraph->Draw();


afpSingleGraph->SetMarkerStyle(21);
afpSingleGraph->SetMarkerColor(2);
afpSingleGraph->Draw("*");

TLegend *leg = new TLegend(0.6,0.7,0.89,0.89);
 leg->AddEntry(afpGraph,"All afterpulse hits");
 leg->AddEntry(afpSingleGraph,"Single afterpulse hit per event");
 leg->Draw("SAME");

canvas8->SaveAs(combinedPlot);



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
