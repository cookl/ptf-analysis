#include "WaveformFitResult.hpp"
#include "ScanPoint.hpp"
#include "TFile.h"
#include "TCanvas.h"
#include "TH1.h"
#include "TF1.h"
#include "TPad.h"
#include "TLegend.h"
#include "TGraph.h"
#include "TGraphErrors.h"

#include <iostream>
#include <vector>
#include <fstream>
#include <math.h>

int main( int argc, char* argv[] ) {

  if ( argc != 4 ){
    std::cerr<<"Usage: ptf_ttree_analysis.app ptf_analysis.root run_number ptfanalysis0\n";
    exit(0);
  }

  // creating variables for afterpulse:  
  int lsrTotal = 0, afpTotal = 0, singleAfpCount = 0;
  int afpTimeThreshold, lsrHeightInt, lsrAreaInt, lsrInt;
  double wfValue, area;
  bool afpCounted;

  // creating histogram file names:
  char fnmT[1024], fnmQ[1024], fnmQlsr[1024], fnmQafp[1024], fnmAfpRate[1024], fnmSingleAfpRate[1024], combinedPlot[1024];
  sprintf(fnmT, "../%s-%s-time-hist.png", argv[2], argv[3]);
  sprintf(fnmQ, "../%s-%s-charge-hist.png", argv[2], argv[3]);
  sprintf(fnmQlsr, "../%s-%s-lsr-hist.png", argv[2], argv[3]);
  sprintf(fnmQafp, "../%s-%s-afp-hist.png", argv[2], argv[3]);
  sprintf(fnmAfpRate, "../%s-%s-afpRate.png", argv[2], argv[3]);
  sprintf(fnmSingleAfpRate, "../%s-%s-single-afpRate.png", argv[2], argv[3]);
  sprintf(combinedPlot, "../%s-%s-combined-afpRate.png", argv[2], argv[3]);

  // creating output file:
  char fnmOut[1024];
  sprintf(fnmOut, "../mpmt-afp-%s.txt", argv[2]);
  std::ofstream outFile;
  
  // arbritary values of a threshold and a time window
  afpTimeThreshold   = 2300;

  // how many sigma under the gaussian should be integrated around the mean?
  double n = 2.0;

  // opening the root file
  TFile * fin = new TFile( argv[1], "read" );

  // adding two options of canvas size
  auto canvas1 = new TCanvas("canvas1","",800,500);
  auto canvas2 = new TCanvas("canvas6","",1600,1000);

  // initializing all histograms and graphs
  auto histT    = new TH1F("histT","",1000,0,8200*2.0*0.48828125);
  auto histQ    = new TH1F("histQ","",500,0,1000*1.5*0.48828125);
  auto histQlsr = new TH1F("histQlsr","",500,0.0,1000*1.5*0.48828125);
  auto histQafp = new TH1F("histQafp","",500,0.0,1000*1.5*0.48828125);
  auto histA    = new TH1F("histA","",500,0.0,1000*1.5*0.48828125);
  auto histAlsr = new TH1F("histAlsr","",500,0.0,1000*1.5*0.48828125);
  auto histAafp = new TH1F("histAafp","",500,0.0,1000*1.5*0.48828125);
  auto afpGraph = new TGraphErrors();
  auto afpSingleGraph = new TGraphErrors();

  // get the waveform fit TTree
  TTree * tt = (TTree*)fin->Get(argv[3]);
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

        area = sqrt(2.*M_PI)*(wf->amp*1000.)*(wf->sigma)*erf(n/sqrt(2) );
        histA->Fill(area);
        
        // check if it is laser pulse. 
        // how many p.e. </to be implemented
        if (wf->pulseTimes[i]<afpTimeThreshold){
          histAlsr->Fill(area);
          histQlsr->Fill(wf->pulseCharges[i]*1000.0);
          lsrTotal++;
        } else if (wf->pulseTimes[i]>=afpTimeThreshold){
          afpTotal++;
          histAafp->Fill(area);
          histQafp->Fill(wf->pulseCharges[i]*1000.0);
          if (afpCounted == false) {
            singleAfpCount++;
            afpCounted = true;
          }
        }
      }             
    }
  }


  // finding the most populated bin for height
  int maxYh = histQlsr->GetBinContent(histQlsr->GetMaximumBin());
  double maxXh = histQlsr->GetXaxis()->GetBinCenter(histQlsr->GetMaximumBin());

  int maxYc = histQlsr->GetBinContent(histAlsr->GetMaximumBin());
  double maxXc = histQlsr->GetXaxis()->GetBinCenter(histAlsr->GetMaximumBin());

  std::cout<<"\nSUMMARY OF TOTAL HISTOGRAMS"<<std::endl;
  std::cout<<"Bin with highest population (mV) at: "<<maxXh<<std::endl;
  std::cout<<"Bin with highest population (mV.ns) at: "<<maxXc<<std::endl;


  // calculating global afterpulse rates
  std::cout<<"\nSUMMARY OF TOTALS"<<std::endl;
  std::cout<<"Lsr pulses: "<< afpTotal << std::endl;
  std::cout<<"Afp: "<< lsrTotal << std::endl;
  std::cout<<"Waveforms w/ afp: " <<singleAfpCount << std::endl;
  std::cout<<"Afp by lsr pulses: "<< double(afpTotal)/lsrTotal*100. << std::endl;
  std::cout<<"Waveforms w/ afp by lsr pulses:"<< double(singleAfpCount)/lsrTotal*100<<std::endl;

  // fixing the pe height and charge
  // double peHeight = 9.52; //< mV 
  double peHeight = maxXh;
  double peArea   = maxXc;
  std::cout<<"pe height used (mV): " << peHeight << std::endl;
  std::cout<<"pe charge used (mv.ns): " <<peArea <<std::endl;

  // i always starts at 0
  // j always starts at 1

  // defining all pe regions we are looking for
  // and initializing the arrays
  int N = 2000;
  double peDivision = 0.5;
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

          // retrieve pulse amplitude
          wfValue = wf->pulseCharges[i]*1000.0/peHeight;
          area = sqrt(2.*M_PI)*(wf->amp*1000.)*(wf->sigma)*erf(n/sqrt(2) );
          
          // transform charge value to p.e. height or p.e. area
          lsrHeightInt = ceil(wf->pulseCharges[0]*1000.0/peHeight/peDivision);
          lsrAreaInt = ceil(area/peArea/peDivision);

          std::cout<<lsrHeightInt<<"   "<<lsrAreaInt<<std::endl;

          // choose one of the above
          // lsrInt = lsrHeightInt;
          lsrInt = lsrAreaInt;

          // count the lsr pulse
          autoLsrCounter[lsrInt] += 1;          
          autoLsrHappened[lsrInt] = true;
         
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

// calculating afterpulse rate.
// as afp per lsr, and afp per total if necessary.
double errorAfpCounts[N], errorSingleAfpCounts[N], errorLsrCounts[N];
double afpRatesError[N], afpSingleEventRatesError[N];
for (int j=1; j<=N; j++){

  // error on pulse counts:
  errorAfpCounts[j] = sqrt(autoAfpCounter[j]);
  errorSingleAfpCounts[j] = sqrt(autoSingleAfpEventConter[j]);
  errorLsrCounts[j] = sqrt(autoLsrCounter[j]);

  // afp rates:
  afpRates[j] = double(autoAfpCounter[j])/double(autoLsrCounter[j]);
  singleEventAfpRates[j] = double(autoSingleAfpEventConter[j])/double(autoLsrCounter[j]);

  // error in afp rates (added in quadrature):
  afpRatesError[j] = afpRates[j]*sqrt( pow(errorAfpCounts[j]/autoAfpCounter[j],2) + pow(errorLsrCounts[j]/autoLsrCounter[j],2)  );
  afpSingleEventRatesError[j] = autoSingleAfpEventConter[j]*sqrt( pow(errorSingleAfpCounts[j]/autoLsrCounter[j],2) + pow(errorLsrCounts[j]/autoLsrCounter[j],2)  );
  
  // printing the rate and error:
  std::cout<<j<<": "<<peRegions[j]<<" p.e.: "<< afpRates[j] << "+/- " <<afpRatesError[j]<<std::endl; 
}


// Drawing the graphs of afpRate by lsr height
int upToWhatPoint=43;

canvas2->cd(1);
afpGraph->Set(upToWhatPoint);  
for (int i=0; i<upToWhatPoint; i++){
  // all afp events:
  afpGraph->SetPoint(i, peRegions[i+1], 100.*afpRates[i+1]);
  afpGraph->SetPointError(i, 0.0, afpRatesError[i+1]);
  
  // single afp event:
  afpSingleGraph->SetPoint(i, peRegions[i+1], 100.0*singleEventAfpRates[i+1]);
  afpSingleGraph->SetPointError(i, 0.0, afpSingleEventRatesError[i+1]);
}

afpGraph->SetTitle("Afterpulse rate");
afpGraph->GetXaxis()->SetTitle("Laser pulse height (p.e.)");
afpGraph->GetYaxis()->SetTitle("Percentage \%");
afpGraph->SetMarkerStyle(20);
afpGraph->SetLineColor(3);
afpGraph->Draw();

TF1* polyFitting = new TF1("poly","[0] + [1]*x^[2]",0.0,peRegions[upToWhatPoint]);
polyFitting->SetParameters(1.,1.,0.5);
polyFitting->SetLineColor(11);
polyFitting->SetLineStyle(1);
afpGraph->Fit(polyFitting,"","",0.0,peRegions[upToWhatPoint+1]);

afpSingleGraph->SetLineColor(1);
afpSingleGraph->SetMarkerStyle(20);
afpSingleGraph->Draw();

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

canvas2->SaveAs(combinedPlot);



  // TIME histogram
  canvas1->cd(1);
  histT->SetTitle("Pulse times");
  histT->GetXaxis()->SetTitle("Time (ns)");
  histT->GetYaxis()->SetTitle("No. of events");
  histT->SetMaximum(1500);
  histT->Draw();
  canvas1->SaveAs(fnmT);

  // ALL PULSES histogram
  canvas2->cd(2);
  histQ->SetTitle("Total pulses");
  histQ->GetXaxis()->SetTitle("Pulse height (mV)");
  histQ->GetYaxis()->SetTitle("No. of events");
  histQ->GetXaxis()->SetRange(0,100);
  histQ->Draw();
  canvas2->SaveAs(fnmQ);

  // LASER PULSES histogram
  canvas2->cd(3);
  histAlsr->SetTitle("Laser pulses");
  histAlsr->GetXaxis()->SetTitle("Laser pulse area (mV.ns)");
  histAlsr->GetYaxis()->SetTitle("No. of events");
  histAlsr->GetXaxis()->SetRange(0,100);
  histAlsr->Draw();
  canvas2->SaveAs(fnmQlsr);  

  // AFTER PULSES histogram
  canvas2->cd(4);
  histQafp->SetTitle("Afterpulses");
  histQafp->GetXaxis()->SetTitle("Afterpulse height (mV)");
  histQafp->GetYaxis()->SetTitle("No. of events");
  histQafp->GetXaxis()->SetRange(0,100);
  histQafp->Draw();
  canvas2->SaveAs(fnmQafp);

  // ALL-IN-ONE histogram
  // </Change the colors in the future.
  canvas2->cd(5);
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
  canvas2->SaveAs("../both.png");

  return 0;
}
