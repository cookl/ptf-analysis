#include "WaveformFitResult.hpp"
#include "ScanPoint.hpp"
#include "TFile.h"
#include "TCanvas.h"
#include "TH1.h"
#include "TAxis.h"

#include <iostream>
#include <vector>


int main( int argc, char* argv[] ) {

  if ( argc != 3 ){
    std::cerr<<"Usage: ptf_ttree_analysis.app ptf_analysis.root run_number\n";
    exit(0);
  }

  // creating variables for afterpulse:
  int lsrCounter1 = 0, lsrCounter2 = 0, lsrCounter3 = 0, lsrCounter4 = 0, lsrCounter5 = 0, lsrCounter6 = 0; 
  int afpCounter1 = 0, afpCounter2 = 0, afpCounter3 = 0, afpCounter4 = 0, afpCounter5 = 0, afpCounter6 = 0; 
  int afpTimeThreshold;
  double afpAverageLsr1, afpAverageLsr2, afpAverageLsr3, afpAverageLsr4, afpAverageLsr5, afpAverageLsr6;
  bool lsr1Happened, lsr2Happened, lsr3Happened, lsr4Happened, lsr5Happened, lsr6Happened;

  // creating histogram file names:
  char fnmT[1024], fnmQ[1024], fnmQlsr[1024], fnmQafp[1024];
  sprintf(fnmT, "../mpmt_time-%s.pdf", argv[2]);
  sprintf(fnmQ, "../mpmt_charge-%s.pdf", argv[2]);
  sprintf(fnmQlsr, "../mpmt_charge-lsr-%s.pdf", argv[2]);
  sprintf(fnmQafp, "../mpmt_charge-afp-%s.pdf", argv[2]);

  // arbritary values of a threshold and a time window
  afpTimeThreshold   = 2300;

  // opening the root file
  TFile * fin = new TFile( argv[1], "read" );

  // adding a canvas and a histogram
  auto canvas1 = new TCanvas("canvas1","",500,500);
  auto canvas2 = new TCanvas("canvas2","",500,500);
  auto canvas3 = new TCanvas("canvas3","",500,500);
  auto canvas4 = new TCanvas("canvas4","",500,500);
  auto histT    = new TH1F("histT","",1000,0,8200);
  auto histQ    = new TH1F("histQ","",100,0.85,1.0);
  auto histQlsr = new TH1F("histQlsr","",100,0.85,1.0);
  auto histQafp = new TH1F("histQafp","",100,0.85,1.0);

  // get the waveform fit TTree
  TTree * tt = (TTree*)fin->Get("ptfanalysis0");
  WaveformFitResult * wf = new WaveformFitResult;
  wf->SetBranchAddresses( tt );
  
  std::cout << "Looping tree " << std::endl;
  // Loop the first scan point and print something
  for(int i = 0; i < tt->GetEntries(); i++){
    tt->GetEvent(i );
    // std::cout<<i<< ". number of found pulses="<<wf->numPulses << std::endl;
    lsr1Happened = false;
    lsr2Happened = false;
    lsr3Happened = false;
    lsr4Happened = false;
    lsr5Happened = false;
    lsr6Happened = false;
    if(wf->numPulses > 0){
      for(int i = 0; i < wf->numPulses; i++){
        // filling up the histogram
        histT->Fill(wf->pulseTimes[i]);
        histQ->Fill(wf->pulseCharges[i]);

        // check if it is laser pulse. 
        // how many p.e. </to be implemented
        if (wf->pulseTimes[i]<afpTimeThreshold){
          histQlsr->Fill(wf->pulseCharges[i]);
          //0.9, 0.92, 0.94, 0.96, 0.98, 1.00
          double wfvalue = wf->pulseCharges[i];
          if (wfvalue<=0.9){
            lsrCounter1++;
            lsr1Happened = true;
          } else if (wfvalue>0.9 && wfvalue<=0.92) {
            lsrCounter2++;
            lsr2Happened = true;
          } else if (wfvalue>0.92 && wfvalue<=0.94) {
            lsrCounter3++;
            lsr3Happened = true;
          } else if (wfvalue>0.94 && wfvalue<=0.96) {
            lsrCounter4++;
            lsr4Happened = true;
          } else if (wfvalue>0.96 && wfvalue<=0.98) {
            lsrCounter5++;
            lsr5Happened = true;
          } else if (wfvalue>0.98) {
            lsrCounter6++;
            lsr6Happened = true;
          }
        }

        // simple time and charge threshold. 
        // if pulse time is higher than afpTimeThreshold,
        // and pulse charge is lower than afpChargeThreshold, 
        // then it is counted as an afterpulse.
        if(wf->pulseTimes[i]>=afpTimeThreshold){
          histQafp->Fill(wf->pulseCharges[i]);
          if (lsr1Happened) {
            afpCounter1++;
          } else if (lsr2Happened) {
            afpCounter2++;
          } else if (lsr3Happened) {
            afpCounter3++;
          } else if (lsr4Happened) {
            afpCounter4++;
          } else if (lsr5Happened) {
            afpCounter5++;
          } else if (lsr6Happened) {
            afpCounter6++;
          }
        }
        
        // terminal output
        // std::cout <<" pulse " << i << " has time = "<<wf->pulseTimes[i]
        //           << " ns " <<  std::endl;
                  
      }

    }
  
  }

  // std::cout<<histQlsr->GetXaxis()->GetBinCenter(27)<<std::endl;
  // std::cout<<histQlsr->GetBinContent(27)<<std::endl;

  // estimating afterpulse rate.
  // as afp per lsr, and afp per total if necessary.
  afpAverageLsr1   = double(afpCounter1)/double(lsrCounter1);
  afpAverageLsr2   = double(afpCounter2)/double(lsrCounter2);
  afpAverageLsr3   = double(afpCounter3)/double(lsrCounter3);
  afpAverageLsr4   = double(afpCounter4)/double(lsrCounter4);
  afpAverageLsr5   = double(afpCounter5)/double(lsrCounter5);
  afpAverageLsr6   = double(afpCounter6)/double(lsrCounter6);
  // std::cout<<"Laser pulses: " << lsrCounter << std::endl;
  // std::cout<<"Afterpulses: " << afpCounter << std::endl;
  std::cout<< afpAverageLsr1*100.<< " percent." << std::endl;
  std::cout<< afpAverageLsr2*100.<< " percent." << std::endl;
  std::cout<< afpAverageLsr3*100.<< " percent." << std::endl;
  std::cout<< afpAverageLsr4*100.<< " percent." << std::endl;
  std::cout<< afpAverageLsr5*100.<< " percent." << std::endl;
  std::cout<< afpAverageLsr6*100.<< " percent." << std::endl;


  std::cout<<"lsrCounter1: "<<lsrCounter1<<std::endl;
  std::cout<<"lsrCounter2: "<<lsrCounter2<<std::endl;
  std::cout<<"lsrCounter3: "<<lsrCounter3<<std::endl;
  std::cout<<"lsrCounter4: "<<lsrCounter4<<std::endl;
  std::cout<<"lsrCounter5: "<<lsrCounter5<<std::endl;
  std::cout<<"lsrCounter6: "<<lsrCounter6<<std::endl;

  std::cout<<"afpCounter1: "<<afpCounter1<<std::endl;
  std::cout<<"afpCounter2: "<<afpCounter2<<std::endl;
  std::cout<<"afpCounter3: "<<afpCounter3<<std::endl;
  std::cout<<"afpCounter4: "<<afpCounter4<<std::endl;
  std::cout<<"afpCounter5: "<<afpCounter5<<std::endl;
  std::cout<<"afpCounter6: "<<afpCounter6<<std::endl;

  // TIME histogram
  canvas1->cd();
  histT->SetTitle("Pulse times");
  histT->GetXaxis()->SetTitle("Time (ns)");
  histT->GetYaxis()->SetTitle("No. of events");
  histT->GetYaxis()->SetRange(0,100);
  histT->Draw();
  canvas1->SaveAs(fnmT);

  // ALL PULSES histogram
  canvas2->cd();
  histQ->SetTitle("Total pulses, charge");
  histQ->GetXaxis()->SetTitle("Charge");
  histQ->GetYaxis()->SetTitle("No. of events");
  histQ->GetXaxis()->SetRange(0,100);
  histQ->Draw();
  canvas2->SaveAs(fnmQ);

  // LASER PULSES histogram
  canvas3->cd();
  histQlsr->SetTitle("Laser pulses, charge");
  histQlsr->GetXaxis()->SetTitle("Charge");
  histQlsr->GetYaxis()->SetTitle("No. of events");
  histQlsr->GetXaxis()->SetRange(0,100);
  histQlsr->Draw();
  canvas3->SaveAs(fnmQlsr);  
  canvas3->SaveAs("../Laser charges.png");  

  // AFTER PULSES histogram
  canvas4->cd();
  histQafp->SetTitle("Afterpulses, charge");
  histQafp->GetXaxis()->SetTitle("Charge");
  histQafp->GetYaxis()->SetTitle("No. of events");
  histQafp->GetXaxis()->SetRange(0,100);
  histQafp->Draw();
  canvas4->SaveAs(fnmQafp);
  return 0;
}
