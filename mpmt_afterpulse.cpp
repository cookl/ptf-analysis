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
  int afpTimeThreshold, lsrHeightInt=0, lsrAreaInt=0, lsrInt = 0;
  double wfValue, area;
  bool afpCounted;

  // creating histogram file names:
  char fnmT[1024], fnmQ[1024], fnmQlsr[1024]; 
  char fnmQafp[1024], fnmAfpRate[1024], fnmSingleAfpRate[1024];
  char combinedPlot[1024], fnmAlsr[1024], fnmAafp[1024];
  char fnmA[1024]; 
  sprintf(fnmT, "../%s-%s-hist-time.png", argv[2], argv[3]);
  sprintf(fnmQ, "../%s-%s-hist-amps-all.png", argv[2], argv[3]);
  sprintf(fnmQlsr, "../%s-%s-hist-amps-lsr.png", argv[2], argv[3]);
  sprintf(fnmQafp, "../%s-%s-hist-amps-afp.png", argv[2], argv[3]);
  sprintf(fnmAfpRate, "../%s-%s-plot-amps-afpRate.png", argv[2], argv[3]);
  sprintf(fnmSingleAfpRate, "../%s-%s-plot-amps-afpSingleRate.png", argv[2], argv[3]);
  sprintf(combinedPlot, "../%s-%s-plot-amps-combined.png", argv[2], argv[3]);
  sprintf(fnmA, "../%s-%s-hist-area-all.png", argv[2], argv[3]);
  sprintf(fnmAlsr, "../%s-%s-hist-area-lsr.png", argv[2], argv[3]);
  sprintf(fnmAafp, "../%s-%s-hist-area-afp.png", argv[2], argv[3]);

  // creating output file:
  char fnmOut[1024];
  sprintf(fnmOut, "../mpmt-afp-%s.txt", argv[2]);
  std::ofstream outFile;
  
  // arbritary values of a threshold and a time window
  afpTimeThreshold   = 2300;

  // how many sigma under the gaussian should be integrated around the mean?
  double spread = 2.0;

  // opening the root file
  TFile * fin = new TFile( argv[1], "read" );

  // adding two options of canvas size
  auto canvas1 = new TCanvas("canvas1","",900,900);
  auto canvas2 = new TCanvas("canvas6","",2500,1000);

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
        // filling up the histograms        
        histT->Fill(wf->pulseTimes[i]);  // Time in ns
        histQ->Fill(wf->pulseCharges[i]*1000.0); // < Multiply by 1000.0 to convert to mV
        histA->Fill(wf->pulseArea); // Area in mV.ns
        
        // check if it is laser pulse. 
        // how many p.e. </to be implemented
        if (wf->pulseTimes[i]<afpTimeThreshold){
          histAlsr->Fill(wf->pulseArea);
          histQlsr->Fill(wf->pulseCharges[i]*1000.0);
          lsrTotal++;
        } else if (wf->pulseTimes[i]>=afpTimeThreshold){
          afpTotal++;
          histAafp->Fill(wf->pulseArea);
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
  double peHeight = 8.00; //< mV 
  // double peHeight = maxXh;
  // double peArea   = maxXc;
  double peArea = 260.0;
  std::cout<<"pe height used (mV): " << peHeight << std::endl;
  std::cout<<"pe charge used (mv.ns): " <<peArea <<std::endl;

  // i always starts at 0
  // j always starts at 1

  // defining all pe regions we are looking for
  // and initializing the arrays
  int lsrMax = 0, N = 5000;
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
  std::cout << "Calculating afterpulse rates per pe." <<std::endl;

  std::vector<double> all_areas;

  for(int i = 0; i < tt->GetEntries(); i++){
    tt->GetEvent(i );
    
    if(wf->numPulses > 0){
      
      bool foundAfpEvent = false;
      for(int i = 0; i < wf->numPulses; i++){
        if (wf->pulseTimes[i]<afpTimeThreshold){
          
          // setting all lsr booleans to false at the start of every loop.
          for (int j=1; j<=N; j++) {
            autoLsrHappened[j] = false;
          }

          // retrieve amplitude and area information and
          // transform pulse values to multiples of p.e. height and p.e. area
          lsrHeightInt = ceil(wf->pulseCharges[0]*1000.0/peHeight/peDivision);
          lsrAreaInt = ceil(wf->pulseArea/peArea/peDivision);
          all_areas.push_back(wf->pulseArea);


          // CHOOSE ONE OF THE ABOVE
          lsrInt = lsrAreaInt;          

          // count the lsr pulse
          autoLsrCounter[lsrInt] += 1;          
          autoLsrHappened[lsrInt] = true;

          // finding the maximum lsrInt
          if (lsrMax < lsrInt) {
            lsrMax = lsrInt;
          }

          // counting the afterpulses generated by
          // each p.e. height
        } else if (wf->pulseTimes[i]>=afpTimeThreshold){
          for (int j=1;j<=N;j++) {
            if (autoLsrHappened[j] == true) {

               autoAfpCounter[j] += 1;
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
int M = lsrMax+1;
double errorAfpCounts[M], errorSingleAfpCounts[M], errorLsrCounts[M];
double afpRatesError[M], afpSingleEventRatesError[M];
for (int j=1; j<lsrMax; j++){

  // error on pulse counts:
  errorAfpCounts[j] = sqrt(autoAfpCounter[j]);
  errorSingleAfpCounts[j] = sqrt(autoSingleAfpEventConter[j]);
  errorLsrCounts[j] = sqrt(autoLsrCounter[j]);

  // afp rates:
  afpRates[j] = double(autoAfpCounter[j])/double(autoLsrCounter[j]);
  singleEventAfpRates[j] = double(autoSingleAfpEventConter[j])/double(autoLsrCounter[j]);

  // error in afp rates, added in quadrature:
  afpRatesError[j] = afpRates[j]*sqrt( pow(errorAfpCounts[j]/autoAfpCounter[j],2) + pow(errorLsrCounts[j]/autoLsrCounter[j],2)  );
  afpSingleEventRatesError[j] = singleEventAfpRates[j]*sqrt( pow(errorSingleAfpCounts[j]/autoLsrCounter[j],2) + pow(errorLsrCounts[j]/autoLsrCounter[j],2)  );
  
  // printing the rate and error:
  std::cout<<j<<": "<<peRegions[j]<<" p.e.: "<< afpRates[j]*100. << "+/- " <<afpRatesError[j]*100.<<std::endl; 
}

std::cout<<"\nROOT's output"<<std::endl;

// Drawing the graphs of afpRate by lsr height
int upToWhatPoint=10;

canvas2->cd(1);
char title1[1204];
sprintf(title1, "Afterpulse rate, %s, %s", argv[2], argv[3]);

afpGraph->Set(upToWhatPoint);  
for (int i=0; i<upToWhatPoint; i++){
  
  // all afp events:
  afpGraph->SetPoint(i, peRegions[i+1], 100.*afpRates[i+1]);
  afpGraph->SetPointError(i, 0.0, afpRatesError[i+1]*100.0);
  
  // single afp event:
  afpSingleGraph->SetPoint(i, peRegions[i+1], 100.0*singleEventAfpRates[i+1]);
  afpSingleGraph->SetPointError(i, 0.0, afpSingleEventRatesError[i+1]*100.0);
}

// afpGraph->SetTitle("Afterpulse rate");
afpGraph->GetXaxis()->SetTitle("Laser pulse area (p.e.)");
afpGraph->GetYaxis()->SetTitle("Percentage \%");
afpGraph->SetMarkerStyle(20);
afpGraph->SetLineColor(1);
afpGraph->Draw("AP");

afpSingleGraph->SetMarkerStyle(21);
afpSingleGraph->SetMarkerColor(3);
afpSingleGraph->Draw("P");

// TF1* polyFitting = new TF1("poly","[0] + [1]*x^[2]",0.0,peRegions[upToWhatPoint]);
// polyFitting->SetParameters(1.,1.);
// polyFitting->SetLineColor(11);
// polyFitting->SetLineStyle(1);
// afpGraph->Fit(polyFitting,"","",0.0,peRegions[upToWhatPoint+1]);

// TF1* polyFitting2 = new TF1("poly","[0] + [1]*x^[2]",0.0,peRegions[upToWhatPoint]);
// polyFitting2->SetParameters(1.,1.);
// polyFitting2->SetLineColor(11);
// polyFitting2->SetLineStyle(2);
// afpSingleGraph->Fit(polyFitting,"","",0.0,peRegions[upToWhatPoint+1]);

TLegend *leg = new TLegend(0.2,0.7,0.4,0.79);
 leg->AddEntry(afpGraph,"All afterpulse hits");
 leg->AddEntry(afpSingleGraph,"Single afterpulse hit per event");
 leg->Draw("SAME");

canvas2->SaveAs(combinedPlot);



  // TIME histogram
  canvas2->cd(2);
  char title2[1024];
  sprintf(title2, "Pulse times, %s, %s", argv[2], argv[3]);
  histT->SetTitle(title2);
  histT->GetXaxis()->SetTitle("Time (ns)");
  histT->GetYaxis()->SetTitle("No. of events");
  histT->SetMaximum(400);
  histT->SetLineColor(1);
  histT->Draw();
  canvas2->SaveAs(fnmT);

  // ALL PULSES histogram
  canvas2->cd(3);
  char title3[1024];
  sprintf(title3, "All, amps, %s, %s", argv[2], argv[3]);  
  histQ->SetTitle(title3);
  histQ->GetXaxis()->SetTitle("Pulse height (mV)");
  histQ->GetYaxis()->SetTitle("No. of events");
  histQ->GetXaxis()->SetRange(0,100);
  histQ->SetLineColor(1);
  histQ->Draw();
  canvas2->SaveAs(fnmQ);

  // LASER PULSES histogram
  canvas2->cd(4);
  char title4[1024];
  sprintf(title4, "Lsr, amps, %s, %s", argv[2], argv[3]);
  histQlsr->SetTitle(title4);
  histQlsr->GetXaxis()->SetTitle("Pulse height (mV)");
  histQlsr->GetYaxis()->SetTitle("No. of events");
  histQlsr->GetXaxis()->SetRangeUser(0,60);
  histQlsr->SetLineColor(1);
  histQlsr->SetLineWidth(2);
  histQlsr->Draw();
  canvas2->SaveAs(fnmQlsr);  

  // AFTER PULSES histogram
  canvas2->cd(5);
  char title5[1024];
  sprintf(title5, "Afp, amps, %s, %s", argv[2], argv[3]);
  histQafp->SetTitle(title5);
  histQafp->GetXaxis()->SetTitle("Pulse height (mV)");
  histQafp->GetYaxis()->SetTitle("No. of events");
  histQafp->GetXaxis()->SetRange(0,100);
  histQafp->SetLineColor(1);
  histQafp->Draw();
  canvas2->SaveAs(fnmQafp);


  // ALL PULSES histogram
  canvas2->cd(6);
  char title6[1024];
  sprintf(title6, "All, area, %s, %s", argv[2], argv[3]);  
  histA->SetTitle(title6);
  histA->GetXaxis()->SetTitle("Pulse area (mV.ns)");
  histA->GetYaxis()->SetTitle("No. of events");
  histA->SetMaximum(4000);
  histA->GetXaxis()->SetRangeUser(0,700);
  histA->SetLineColor(1);
  histA->Draw();
  canvas2->SaveAs(fnmA);

  // LASER PULSES histogram
  canvas2->cd(7);
  char title7[1024];
  sprintf(title7, "Lsr, area, %s, %s", argv[2], argv[3]);
  histAlsr->SetTitle(title7);
  histAlsr->GetXaxis()->SetTitle("Pulse area (mV.ns)");
  histAlsr->GetYaxis()->SetTitle("No. of events");
  histAlsr->GetXaxis()->SetRange(0,700);
  histAlsr->SetMaximum(2000);
  histAlsr->SetLineColor(1);
  histAlsr->Draw();
  canvas2->SaveAs(fnmAlsr);  

  // AFTER PULSES histogram
  canvas2->cd(8);
  char title8[1024];
  sprintf(title8, "Afp, area, %s, %s", argv[2], argv[3]);
  histAafp->SetTitle(title8);
  histAafp->GetXaxis()->SetTitle("Pulse area (mV.ns)");
  histAafp->GetYaxis()->SetTitle("No. of events");
  histAafp->GetXaxis()->SetRange(0,700);
  histAafp->Draw();
  canvas2->SaveAs(fnmAafp);
  return 0;
}
