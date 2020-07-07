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
  int lsrCounter1 = 0, lsrCounter2 = 0, lsrCounter3 = 0, lsrCounter4 = 0, lsrCounter5 = 0, lsrCounter6 = 0, lsrCounter7 = 0; 
  int lsrCounter8 = 0, lsrCounter9 = 0, lsrCounter10 = 0, lsrCounter11 = 0, lsrCounter12 = 0, lsrCounter13 = 0, lsrCounter14 = 0; 
  int lsrCounter15 = 0, lsrCounter16 = 0;
  
  int lsrTotal = 0, afpTotal = 0;

  int afpCounter1 = 0, afpCounter2 = 0, afpCounter3 = 0, afpCounter4 = 0, afpCounter5 = 0, afpCounter6 = 0, afpCounter7 = 0; 
  int afpCounter8 = 0, afpCounter9 = 0, afpCounter10 = 0, afpCounter11 = 0, afpCounter12 = 0, afpCounter13 = 0, afpCounter14 = 0; 
  int afpCounter15 = 0, afpCounter16 = 0;
  
  int afpTimeThreshold;
  
  double wfValue;

  bool lsr1Happened, lsr2Happened, lsr3Happened, lsr4Happened, lsr5Happened, lsr6Happened, lsr7Happened;
  bool lsr8Happened, lsr9Happened, lsr10Happened, lsr11Happened, lsr12Happened, lsr13Happened, lsr14Happened;
  bool lsr15Happened, lsr16Happened;

  int afpPeCounter1 = 0, afpPeCounter2 = 0, afpPeCounter3 = 0, afpPeCounter4 = 0, afpPeCounter5 = 0, afpPeCounter6 = 0; 
  int afpPeCounter7 = 0, afpPeCounter8 = 0, afpPeCounter9 = 0, afpPeCounter10 = 0, afpPeCounter11 = 0, afpPeCounter12 = 0;
  int afpPeCounter13 = 0, afpPeCounter14 = 0, afpPeCounter15 = 0, afpPeCounter16 = 0;

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
  auto histT    = new TH1F("histT","",1000,0,8200*2.0*0.48828125);
  auto histQ    = new TH1F("histQ","",500,0,1000*1.5*0.48828125);
  auto histQlsr = new TH1F("histQlsr","",500,0.0,1000*1.5*0.48828125);
  auto histQafp = new TH1F("histQafp","",500,0.0,1000*1.5*0.48828125);

  // get the waveform fit TTree
  TTree * tt = (TTree*)fin->Get("ptfanalysis0");
  WaveformFitResult * wf = new WaveformFitResult;
  wf->SetBranchAddresses( tt );
  
  // if (argv[1]=="247") {
  std::cout << "\n" << std::endl;
  std::cout << "LOOPING TREE, 1:" << std::endl;
  std::cout << "Creating histograms and finding total afp rate."<<std::endl;
  for(int i = 0; i < tt->GetEntries(); i++){
    tt->GetEvent(i );
    // std::cout<<i<< ". number of found pulses="<<wf->numPulses << std::endl;
    if(wf->numPulses > 0){
      for(int i = 0; i < wf->numPulses; i++){
        // filling up the histogram
        histT->Fill(wf->pulseTimes[i]);
        histQ->Fill(wf->pulseCharges[i]*1000.0);

        // check if it is laser pulse. 
        // how many p.e. </to be implemented
        if (wf->pulseTimes[i]<afpTimeThreshold){
          histQlsr->Fill(wf->pulseCharges[i]*1000.0);
          lsrTotal++;
        } else if (wf->pulseTimes[i]>=afpTimeThreshold){
          afpTotal++;
          histQafp->Fill(wf->pulseCharges[i]*1000.0);
        }
      }             
    }
  }

  std::cout<<"Total afp rate: "<< double(afpTotal)/lsrTotal*100. << std::endl;

  // // finding the most populated bin
  // int maxY = histQlsr->GetBinContent(histQlsr->GetMaximumBin());
  // double maxX = histQlsr->GetXaxis()->GetBinCenter(histQlsr->GetMaximumBin());
  // std::cout<<"Population of highest bin: "<<maxY<<std::endl;
  // std::cout<<"Position of highest bin (mV): "<<maxX<<std::endl;

  // } else {
  double maxX = 9.52; //< mV 
  // }

  std::cout<<"pe value used: " << maxX << std::endl;

  // start loop #2
  std::cout << "\n" << std::endl;
  std::cout << "LOOPING TREE, 2:" << std::endl;
  std::cout << "Creating histograms and finding maximum."<<std::endl;

  for(int i = 0; i < tt->GetEntries(); i++){
    tt->GetEvent(i );
    // finding pulses on the waveform
    if(wf->numPulses > 0){
      for(int i = 0; i < wf->numPulses; i++){

        // classifying laser pulses as multiples of pe
        if (wf->pulseTimes[i]<afpTimeThreshold){

          lsr1Happened=false;
          lsr2Happened=false;
          lsr3Happened=false;
          lsr4Happened=false;
          lsr5Happened=false;
          lsr6Happened=false;
          lsr7Happened=false;
          lsr8Happened=false;
          lsr9Happened=false;
          lsr10Happened=false;
          lsr11Happened=false;
          lsr12Happened=false;
          lsr13Happened=false;
          lsr14Happened=false;
          lsr15Happened=false;
          lsr16Happened=false;

          wfValue = (wf->pulseCharges[i]*1000)/maxX;
          if (wfValue <=0.5) {
            lsrCounter1++;
            lsr1Happened = true;
          } else if (wfValue >0.5 && wfValue <=1.0) {
            lsrCounter2++;
            lsr2Happened = true;
          } else if (wfValue >1.0 && wfValue <=1.5) {
            lsrCounter3++;
            lsr3Happened = true;
          } else if (wfValue >1.5 && wfValue <=2.0) {
            lsrCounter4++;
            lsr4Happened = true;
          } else if (wfValue >2.0 && wfValue <=2.5) {
            lsrCounter5++;
            lsr5Happened = true;
          } else if (wfValue >2.5 && wfValue <=3.0) {
            lsrCounter6++;
            lsr6Happened = true;
          } else if (wfValue >3.0 && wfValue <=3.5){
            lsrCounter7++;
            lsr7Happened = true;
          } else if (wfValue >3.5 && wfValue <=4.0){
            lsrCounter8++;
            lsr8Happened = true;
          } else if (wfValue >4.0 && wfValue <=4.5){
            lsrCounter9++;
            lsr9Happened = true;
          } else if (wfValue >4.5 && wfValue <=5.0){
            lsrCounter10++;
            lsr10Happened = true;
          } else if (wfValue >5.0 && wfValue <=5.5){
            lsrCounter11++;
            lsr11Happened = true;
          } else if (wfValue >5.5 && wfValue <=6.0){
            lsrCounter12++;
            lsr12Happened = true;
          } else if (wfValue >6.0 && wfValue <=6.5){
            lsrCounter13++;
            lsr13Happened = true;
          } else if (wfValue >6.5 && wfValue <=7.0){
            lsrCounter14++;
            lsr14Happened = true;
          } else if (wfValue >7.0 && wfValue <=7.5){
            lsrCounter15++;
            lsr15Happened = true;
          } else if (wfValue >7.5){
            lsrCounter16++;
            lsr16Happened = true;
          } 
        
        } 
        
        // now the afterpulses:
        if (wf->pulseTimes[i]>=afpTimeThreshold){

          double afpValue = wf->pulseCharges[i]*1000/maxX;

          if (lsr1Happened) {
            afpCounter1++;
            lsr1Happened = false;
            if (afpValue > 0.9 && afpValue < 1.1) {
            afpPeCounter1++;
            }
          } else if (lsr2Happened) {
            afpCounter2++;
            lsr2Happened = false;
            if (afpValue > 0.9 && afpValue < 1.1) {
            afpPeCounter2++;
            }
          } else if (lsr3Happened) {
            afpCounter3++;
            lsr3Happened = false;
            if (afpValue > 0.9 && afpValue < 1.1) {
            afpPeCounter3++;
            }
          } else if (lsr4Happened) {
            afpCounter4++;
            lsr4Happened = false;
            if (afpValue > 0.9 && afpValue < 1.1) {
            afpPeCounter4++;
            }
          } else if (lsr5Happened) {
            afpCounter5++;
            lsr5Happened = false;
            if (afpValue > 0.9 && afpValue < 1.1) {
            afpPeCounter5++;
            }
          } else if (lsr6Happened) {
            afpCounter6++;
            lsr6Happened = false;
            if (afpValue > 0.9 && afpValue < 1.1) {
            afpPeCounter6++;
            }
          } else if (lsr7Happened) {
            afpCounter7++;
            lsr7Happened = false;
            if (afpValue > 0.9 && afpValue < 1.1) {
            afpPeCounter7++;
            }
          } else if (lsr8Happened) {
            afpCounter8++;
            lsr8Happened = false;
            if (afpValue > 0.9 && afpValue < 1.1) {
            afpPeCounter8++;
            }
          } else if (lsr9Happened) {
            afpCounter9++;
            lsr9Happened = false;
            if (afpValue > 0.9 && afpValue < 1.1) {
            afpPeCounter9++;
            }
          } else if (lsr10Happened) {
            afpCounter10++;
            lsr10Happened = false;
            if (afpValue > 0.9 && afpValue < 1.1) {
            afpPeCounter10++;
            }
          }else if (lsr11Happened) {
            afpCounter11++;
            lsr11Happened = false;
            if (afpValue > 0.9 && afpValue < 1.1) {
            afpPeCounter11++;
            }
          }else if (lsr12Happened) {
            afpCounter12++;
            lsr12Happened = false;
            if (afpValue > 0.9 && afpValue < 1.1) {
            afpPeCounter12++;
            }
          }else if (lsr13Happened) {
            afpCounter13++;
            lsr13Happened = false;
            if (afpValue > 0.9 && afpValue < 1.1) {
            afpPeCounter13++;
            }
          }else if (lsr14Happened) {
            afpCounter14++;
            lsr14Happened = false;
            if (afpValue > 0.9 && afpValue < 1.1) {
            afpPeCounter14++;
            }
          }else if (lsr15Happened) {
            afpCounter15++;
            lsr15Happened = false;
            if (afpValue > 0.9 && afpValue < 1.1) {
            afpPeCounter15++;
            }
          }else if (lsr16Happened) {
            afpCounter16++;
            lsr16Happened = false;
            if (afpValue > 0.9 && afpValue < 1.1) {
            afpPeCounter16++;
            }
            std::cout<<afpValue<<", "<<wfValue<<std::endl;
          }
        }
      }             
    }
  }


  // estimating afterpulse rate.
  // as afp per lsr, and afp per total if necessary.
  std::vector<double> afpAverages;
  std::vector<double> pe;
  std::vector<double> afpPeFraction;

  std::cout<<"pe Information:"<<std::endl;
  for (int i=0; i<=15; i++) {
    pe.push_back(double(0.5*(i+1)));
    std::cout<<pe[i]<<", ";
  }
  std::cout<<std::endl;

  afpAverages.push_back(double(afpCounter1)/double(lsrCounter1));
  afpAverages.push_back(double(afpCounter2)/double(lsrCounter2));
  afpAverages.push_back(double(afpCounter3)/double(lsrCounter3));
  afpAverages.push_back(double(afpCounter4)/double(lsrCounter4));
  afpAverages.push_back(double(afpCounter5)/double(lsrCounter5));
  afpAverages.push_back(double(afpCounter6)/double(lsrCounter6));
  afpAverages.push_back(double(afpCounter7)/double(lsrCounter7));
  afpAverages.push_back(double(afpCounter8)/double(lsrCounter8));
  afpAverages.push_back(double(afpCounter9)/double(lsrCounter9));
  afpAverages.push_back(double(afpCounter10)/double(lsrCounter10));
  afpAverages.push_back(double(afpCounter11)/double(lsrCounter11));
  afpAverages.push_back(double(afpCounter12)/double(lsrCounter12));
  afpAverages.push_back(double(afpCounter13)/double(lsrCounter13));
  afpAverages.push_back(double(afpCounter14)/double(lsrCounter14));
  afpAverages.push_back(double(afpCounter15)/double(lsrCounter15));
  afpAverages.push_back(double(afpCounter16)/double(lsrCounter16));

  afpPeFraction.push_back(double(afpPeCounter1)/double(lsrCounter1));
  afpPeFraction.push_back(double(afpPeCounter2)/double(lsrCounter2));
  afpPeFraction.push_back(double(afpPeCounter3)/double(lsrCounter3));
  afpPeFraction.push_back(double(afpPeCounter4)/double(lsrCounter4));
  afpPeFraction.push_back(double(afpPeCounter5)/double(lsrCounter5));
  afpPeFraction.push_back(double(afpPeCounter6)/double(lsrCounter6));
  afpPeFraction.push_back(double(afpPeCounter7)/double(lsrCounter7));
  afpPeFraction.push_back(double(afpPeCounter8)/double(lsrCounter8));
  afpPeFraction.push_back(double(afpPeCounter9)/double(lsrCounter9));
  afpPeFraction.push_back(double(afpPeCounter10)/double(lsrCounter10));
  afpPeFraction.push_back(double(afpPeCounter11)/double(lsrCounter11));
  afpPeFraction.push_back(double(afpPeCounter12)/double(lsrCounter12));
  afpPeFraction.push_back(double(afpPeCounter13)/double(lsrCounter13));
  afpPeFraction.push_back(double(afpPeCounter14)/double(lsrCounter14));
  afpPeFraction.push_back(double(afpPeCounter15)/double(lsrCounter15));
  afpPeFraction.push_back(double(afpPeCounter16)/double(lsrCounter16));
  
  std::cout<<"afp Rate:"<<std::endl;
  for (auto afpRate : afpAverages) {
    std::cout<<afpRate*100.0<<",";
  }
  std::cout<<std::endl;

  std::cout<<"One pe fraction:"<<std::endl;
  for (auto afpPe : afpPeFraction) {
    std::cout<<afpPe*100.0<<", ";
  }
  std::cout<<std::endl;

  // outFile.open(fnmOut);
  // for (int i=0; i<=15; i++) {
  //   // char line[1024];
  //   // sprintf(line,"             %f             %f",pe[i],afpAverages[i]);
  //   // outFile << line << std::endl;
  //   outFile << pe[i] << std::setw(18) << afpAverages[i] << std::endl;
  //   // std::cout<<pe[i]<<", "<<afpAverages[i]<<std::endl;
  // }




  std::cout<<"lsrCounter1: "<<lsrCounter1<<std::endl;
  std::cout<<"lsrCounter2: "<<lsrCounter2<<std::endl;
  std::cout<<"lsrCounter3: "<<lsrCounter3<<std::endl;
  std::cout<<"lsrCounter4: "<<lsrCounter4<<std::endl;
  std::cout<<"lsrCounter5: "<<lsrCounter5<<std::endl;
  std::cout<<"lsrCounter6: "<<lsrCounter6<<std::endl;
  std::cout<<"lsrCounter7: "<<lsrCounter7<<std::endl;
  std::cout<<"lsrCounter8: "<<lsrCounter8<<std::endl;
  std::cout<<"lsrCounter9: "<<lsrCounter9<<std::endl;
  std::cout<<"lsrCounter10: "<<lsrCounter10<<std::endl;
  std::cout<<"lsrCounter11: "<<lsrCounter11<<std::endl;
  std::cout<<"lsrCounter12: "<<lsrCounter12<<std::endl;
  std::cout<<"lsrCounter13: "<<lsrCounter13<<std::endl;
  std::cout<<"lsrCounter14: "<<lsrCounter14<<std::endl;
  std::cout<<"lsrCounter15: "<<lsrCounter15<<std::endl;
  std::cout<<"lsrCounter16: "<<lsrCounter16<<std::endl;

  std::cout<<"afpCounter1: "<<afpCounter1<<std::endl;
  std::cout<<"afpCounter2: "<<afpCounter2<<std::endl;
  std::cout<<"afpCounter3: "<<afpCounter3<<std::endl;
  std::cout<<"afpCounter4: "<<afpCounter4<<std::endl;
  std::cout<<"afpCounter5: "<<afpCounter5<<std::endl;
  std::cout<<"afpCounter6: "<<afpCounter6<<std::endl;
  std::cout<<"afpCounter7: "<<afpCounter7<<std::endl;
  std::cout<<"afpCounter8: "<<afpCounter8<<std::endl;
  std::cout<<"afpCounter9: "<<afpCounter9<<std::endl;
  std::cout<<"afpCounter10: "<<afpCounter10<<std::endl;
  std::cout<<"afpCounter11: "<<afpCounter11<<std::endl;
  std::cout<<"afpCounter12: "<<afpCounter12<<std::endl;
  std::cout<<"afpCounter13: "<<afpCounter13<<std::endl;
  std::cout<<"afpCounter14: "<<afpCounter14<<std::endl;
  std::cout<<"afpCounter15: "<<afpCounter15<<std::endl;
  std::cout<<"afpCounter16: "<<afpCounter16<<std::endl;

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
  histQ->GetXaxis()->SetTitle("Height (mV)");
  histQ->GetYaxis()->SetTitle("No. of events");
  histQ->GetXaxis()->SetRange(0,100);
  histQ->Draw();
  canvas2->SaveAs(fnmQ);

  // LASER PULSES histogram
  canvas3->cd();
  histQlsr->SetTitle("Laser pulses");
  histQlsr->GetXaxis()->SetTitle("Height (mV)");
  histQlsr->GetYaxis()->SetTitle("No. of events");
  histQlsr->GetXaxis()->SetRange(0,100);
  histQlsr->Draw();
  canvas3->SaveAs(fnmQlsr);  

  // AFTER PULSES histogram
  canvas4->cd();
  histQafp->SetTitle("Afterpulses");
  histQafp->GetXaxis()->SetTitle("Pulse (mV)");
  histQafp->GetYaxis()->SetTitle("No. of events");
  histQafp->GetXaxis()->SetRange(0,100);
  histQafp->Draw();
  canvas4->SaveAs(fnmQafp);
  return 0;
}
