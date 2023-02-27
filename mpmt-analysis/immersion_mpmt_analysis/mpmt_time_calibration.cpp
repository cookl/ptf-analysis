// Hassan Elshabasy mPMT Time Calibration
// helshabasy@triumf.ca
// 02/17/23


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
#include "TProfile.h"

#include <iostream>
#include <vector>
#include <string>
#include <cstring>

#include <algorithm>
#include <numeric>
//using namespace std;


class Pulse {
    public:
		int event;
		int channel;
        double height;
		double midtime;
		double mintime;
        Pulse(int aEvent, int aChannel, double aHeight, double aMidTime, double aMinTime){
            event = aEvent;
			channel = aChannel;
            height = aHeight;
			midtime = aMidTime;	
			mintime = aMinTime;
        }
};

bool notSame(Pulse pulse1, Pulse pulse2){
	if((pulse1.channel == pulse2.channel) && (pulse1.midtime == pulse2.midtime)){
		return false;
	} else {
		return true;
	}
}


int main( int argc, char* argv[] ) {

    if ( argc != 2 ){
        std::cerr<<"Usage: ptf_ttree_analysis.app ptf_analysis.root\n";
        exit(0);
    }

    TFile * fin = new TFile( argv[1], "read" );
    TTree * tt[20];
    WaveformFitResult * wf[20];
	TH1F * h1[20]; // pulse times per channel
	TH1F * h2[20]; // pulse time differences relative to channel 1
	
    for(int j = 0; j < 20; j++){ // To initialize TTree and branches
	
		char branch_name[100];
		sprintf(branch_name,"ptfanalysis%i",j);
		tt[j] = (TTree*)fin->Get(branch_name);
			   
		wf[j] = new WaveformFitResult;
		if(tt[j]) wf[j]->SetBranchAddresses( tt[j] );
		
		h1[j] = new TH1F("pulse_time", "Pulse Time",200,2300,2400);
		h2[j] = new TH1F("pulse_time_difference", "Pulse Time Difference",50,0,50);
    }
	
	int myChannels[] = {0, 1, 2, 3, 4, 5, 6, 7, 10, 11, 12, 13, 14, 15, 16, 17};
	
	double myAverages[20];
	std::vector<double> myDifferences[20] = {};
	
	for(int i = 0; i < tt[0]->GetEntries()-1; i++){ // loop over events
		
		std::cout << "Event: " << i << std::endl;
		
		double sameEvent[17];
		
		for(int j : myChannels){ // loop over channels
			tt[j]->GetEvent(i);
			
			std::cout << "Channel: " << j << " Number of Pulses: " << wf[j]->numPulses << std::endl;
			
			std::vector<Pulse> sameChannel = {};
			
			for(int k = 0; k < wf[j]->numPulses; k++){ // loop over each pulse
				
				// pulse height
				double pulse_height = wf[j]->pulseCharges[k]*1000.0;
				
				// pulse midtime
				double pulse_midtime = wf[j]->pulseTimeErr[k];
				
				// pulse mintime
				double pulse_mintime = wf[j]->pulseTimes[k];
				
				Pulse pulse(i, j, pulse_height, pulse_midtime, pulse_mintime);
				sameChannel.push_back(pulse);
				
			}
			
			if(sameChannel.size() == 0){
				continue;
			}
			
			Pulse LED = sameChannel[0];
			for(Pulse pulse : sameChannel){
				if(pulse.height > LED.height){
					LED = pulse;
				}
			}
			
			std::cout << "Mid Time: " <<  LED.midtime << " Min Time: " << LED.mintime << " Pulse Height: " << LED.height << std::endl;
			h1[j]->Fill(LED.midtime);
			
			sameEvent[j] = LED.midtime;
		}
		
		for(int j : myChannels){
			if(abs(sameEvent[j] - sameEvent[1]) <= 50){
				myDifferences[j].push_back(abs(sameEvent[j] - sameEvent[1]));
			}
			std::cout << "Channel " << j << " Time Difference: " << abs(sameEvent[j] - sameEvent[1]) << std::endl;
			h2[j]->Fill(abs(sameEvent[j] - sameEvent[1]));
		}
	}
	
	for(int j : myChannels){
		myAverages[j] = std::accumulate(myDifferences[j].begin(), myDifferences[j].end(), 0.0) / myDifferences[j].size();
		std::cout << "Channel " << j << " Average Time Difference: " << myAverages[j] << std::endl;
	}
	
	// Print Pulse Time Histograms
	for(int j : myChannels){
		TCanvas * c1 = new TCanvas("C1");
		
		char hist_title[100];
		sprintf(hist_title,"Channel %d Pulse Times",j);
		h1[j]->SetTitle(hist_title);
		
		h1[j]->Draw();
		h1[j]->GetXaxis()->SetTitle("Pulse Times (ns)");
		h1[j]->GetYaxis()->SetTitle("Number of events");
		gStyle->SetOptFit(11);
		
		char png_name[100];
		sprintf(png_name,"mpmt_pulse_time_%d_LED.png",j);
		c1->SaveAs(png_name);
	}
	
	// Print Pulse Time Difference Histograms
	for(int j : myChannels){
		TCanvas * c2 = new TCanvas("C2");
		
		char hist_title[100];
		sprintf(hist_title,"Channel %d Pulse Time Differences",j);
		h2[j]->SetTitle(hist_title);
		
		h2[j]->Draw();
		h2[j]->GetXaxis()->SetTitle("Pulse Times (ns)");
		h2[j]->GetYaxis()->SetTitle("Number of events");
		gStyle->SetOptFit(11);
		
		char png_name[100];
		sprintf(png_name,"mpmt_pulse_time_difference_%d_LED.png",j);
		c2->SaveAs(png_name);
	}
	
	std::cout << "Number of Events: " << tt[0]->GetEntries() << std::endl;
	
	fin->Close();
	return 0;
}
	
