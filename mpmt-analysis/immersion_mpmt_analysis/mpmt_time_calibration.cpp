// Hassan Elshabasy mPMT Time Calibration
// helshabasy@triumf.ca
// 04/26/23


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
		double charge;
		double CFDtime;
		double time;
		double fittedtime;
		//double fittedtime;
        Pulse(int aEvent, int aChannel, double aHeight, double aCharge, double aCFDTime, double aTime, double aFittedTime){
            event = aEvent;
			channel = aChannel;
            height = aHeight;
			charge = aCharge;
			CFDtime = aCFDTime;	
			time = aTime;
			fittedtime = aFittedTime;
        }
};

bool notSame(Pulse pulse1, Pulse pulse2){
	if((pulse1.channel == pulse2.channel) && (pulse1.CFDtime == pulse2.CFDtime)){
		return false;
	} else {
		return true;
	}
}

int time_diff = 100;

bool isCoincidentEvent(int x, std::vector<Pulse> lst){ // Checks if lst is a x-fold coincident event containing a pulse in channel 2
	
	std::vector<Pulse> coincident_pulses = {}; // time_cut
	for(Pulse pulse_check : lst){
		for(Pulse pulse : lst){
			if((abs(pulse_check.time - pulse.time) < time_diff) && (notSame(pulse_check, pulse))){
				coincident_pulses.push_back(pulse_check);
				break;
			}
		}
	}
	
	std::vector<int> coincident_channels = {};
	for(Pulse pulse : coincident_pulses){ // find how many channels involved
		coincident_channels.push_back(pulse.channel);
	}
	
	std::sort(coincident_channels.begin(), coincident_channels.end());
	coincident_channels.erase(std::unique(coincident_channels.begin(), coincident_channels.end()), coincident_channels.end());
	
	bool isFound = (std::find(coincident_channels.begin(), coincident_channels.end(), 2) != coincident_channels.end());
	
	if((coincident_channels.size() >= x) && isFound){
		return true;
	} else {
		return false;
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
	TH1F * h1[20]; // pulse times CFD per channel
	TH1F * h2[20]; // pulse time CFD differences relative to channel 2
	TH1F * h5[20]; // pule time difference between CFDtime and time
	TH1F * h6[20]; // original pulse times per channel
	TH1F * h8[20]; // pulse time difference relative to mean per channel
	TH1F * h9[20]; // pulse fitted times per channel
	TH1F * h11[20]; // pulse time difference between CFDtime and fittedtime
	TH1F * h12[20]; // pulse time fitted differences relative to channel 2
	TH1F * h14[20]; // pulse charge histograms
	
    for(int j = 0; j < 20; j++){ // To initialize TTree and branches
	
		char branch_name[100];
		sprintf(branch_name,"ptfanalysis%i",j);
		tt[j] = (TTree*)fin->Get(branch_name);
			   
		wf[j] = new WaveformFitResult;
		if(tt[j]) wf[j]->SetBranchAddresses( tt[j] );
		
		h1[j] = new TH1F("pulse_time_CFD", "Pulse Time CFD",200,0,8000);
		//h2[j] = new TH1F("pulse_time_CFD_difference", "Pulse Time CFD Difference",100,-25,25); // relative to channel 2 
		h5[j] = new TH1F("pulse_time_cfd-time", "Pulse Time Difference Between CFDtime and time", 100, -25, 25);
		h6[j] = new TH1F("pulse_time_original", "Pulse Time Original",200,2300,2400);
		h8[j] = new TH1F("pulse_time_difference_mean_channel", "Pulse Time Difference Relative to Mean per Channel", 100, -25, 25);
		h9[j] = new TH1F("pulse_time_fitted", "Pulse Time Fitted",200,2300,2400);
		h11[j] = new TH1F("pulse_time_cfd-fitted", "Pulse Time Difference Between CFDTime and FittedTime", 100, 20, 35);
		h12[j] = new TH1F("pulse_time_fitted_difference", "Pulse Time Fitted Difference",100,-25,25); // relative to channel 2
		h14[j] = new TH1F("pulse_charge", "Pulse Charge",100,0,30);
		
		// Zooming in on Different Sections of the Histogram Depending on Channel Number
		switch(j){
		case 0:
			h2[j] = new TH1F("pulse_time_CFD_difference", "Pulse Time CFD Difference",100,-3,1); // relative to channel 2
			break;
		case 1:
			h2[j] = new TH1F("pulse_time_CFD_difference", "Pulse Time CFD Difference",100,-2,1); // relative to channel 2
			break;
		case 2:
			h2[j] = new TH1F("pulse_time_CFD_difference", "Pulse Time CFD Difference",100,-25,25); // relative to channel 2
			break;
		case 3:
			h2[j] = new TH1F("pulse_time_CFD_difference", "Pulse Time CFD Difference",100,-3,1); // relative to channel 2
			break;
		case 4:
			h2[j] = new TH1F("pulse_time_CFD_difference", "Pulse Time CFD Difference",100,-1,2); // relative to channel 2
			break;
		case 5:
			h2[j] = new TH1F("pulse_time_CFD_difference", "Pulse Time CFD Difference",100,-4,0); // relative to channel 2
			break;
		case 6:
			h2[j] = new TH1F("pulse_time_CFD_difference", "Pulse Time CFD Difference",100,-4,1); // relative to channel 2
			break;
		case 7:
			h2[j] = new TH1F("pulse_time_CFD_difference", "Pulse Time CFD Difference",100,-1,3); // relative to channel 2
			break;
		case 10:
			h2[j] = new TH1F("pulse_time_CFD_difference", "Pulse Time CFD Difference",100,-20,-14); // relative to channel 2
			break;
		case 11:
			h2[j] = new TH1F("pulse_time_CFD_difference", "Pulse Time CFD Difference",100,-20,-16); // relative to channel 2
			break;
		case 12:
			h2[j] = new TH1F("pulse_time_CFD_difference", "Pulse Time CFD Difference",100,-11,-5); // relative to channel 2
			break;
		case 13:
			h2[j] = new TH1F("pulse_time_CFD_difference", "Pulse Time CFD Difference",100,-14,-11); // relative to channel 2
			break;
		case 14:
			h2[j] = new TH1F("pulse_time_CFD_difference", "Pulse Time CFD Difference",100,-8,-3); // relative to channel 2
			break;
		case 15:
			h2[j] = new TH1F("pulse_time_CFD_difference", "Pulse Time CFD Difference",100,-9,-3); // relative to channel 2
			break;
		case 16:
			h2[j] = new TH1F("pulse_time_CFD_difference", "Pulse Time CFD Difference",100,-10,-5); // relative to channel 2
			break;
		case 17:
			h2[j] = new TH1F("pulse_time_CFD_difference", "Pulse Time CFD Difference",100,-13,-6); // relative to channel 2
			break;
		}
    }
	
	TH1F * h3 = new TH1F("pulse_time_difference_mean", "Pulse Time Difference Relative to Mean", 100, -25, 25);
	
	int myChannels[] = {0, 1, 2, 3, 4, 5, 6, 7, 10, 11, 12, 13, 14, 15, 16, 17};
	int nonSaturatingChannels[] = {1, 2, 3, 4, 5, 7, 12, 13, 14}; // Channels that don't saturate when using internal LED

	double myTimeDiffMeansCFD[20];
	double myTimeDiffSTDCFD[20];
	
	double myTimeDiffMeansFitted[20];
	double myTimeDiffSTDFitted[20];
	
	double myChargeMeans[20];
	std::vector<double> myCharges[20] = {};
	
	for(int i = 0; i < tt[0]->GetEntries()-1; i++){ // loop over events
		
		std::cout << "Event: " << i << std::endl;
		
		double sameEventCFD[20];
		double sameEventFitted[20];
		std::vector<Pulse> sameEventPulse;
		
		for(int j : myChannels){ // loop over channels
			tt[j]->GetEvent(i);
			
			std::cout << "Channel: " << j << " Number of Pulses: " << wf[j]->numPulses << std::endl;
			
			std::vector<Pulse> sameChannel = {};
			
			for(int k = 0; k < wf[j]->numPulses; k++){ // loop over each pulse
				
				// pulse height
				double pulse_height = wf[j]->pulseCharges[k]*1000.0;
				
				// pulse charge
				double pulse_charge = pulse_height/7.7;
				
				// pulse CFDtime
				double pulse_CFDtime = wf[j]->pulseTimesCFD[k];
				
				// pulse time
				double pulse_time = wf[j]->pulseTimes[k];
				
				//pulse fittedTime
				double pulse_fittedtime = wf[j]->mean;
				
				Pulse pulse(i, j, pulse_height, pulse_charge, pulse_CFDtime, pulse_time, pulse_fittedtime);
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
			
			std::cout << "CFD Time: " <<  LED.CFDtime << " Time: " << LED.time << " Pulse Height: " << LED.height << std::endl;
			sameEventPulse.push_back(LED);
		}
		
		if(sameEventPulse.size() == 0){
			continue;
		}
		
		if(isCoincidentEvent(2, sameEventPulse)){
			for(Pulse pulse : sameEventPulse){
				//if(pulse.channel == 1 || pulse.channel == 2 || pulse.channel == 3 || pulse.channel == 4 || pulse.channel == 6 ||
				   //pulse.channel == 7 || pulse.channel == 12 || pulse.channel == 13 || pulse.channel == 14){
					h1[pulse.channel]->Fill(pulse.CFDtime);
					h14[pulse.channel]->Fill(pulse.charge);
					h5[pulse.channel]->Fill(pulse.CFDtime - pulse.time);
					h6[pulse.channel]->Fill(pulse.time);
					h9[pulse.channel]->Fill(pulse.fittedtime);
				    h11[pulse.channel]->Fill(pulse.CFDtime - pulse.fittedtime);
					sameEventCFD[pulse.channel] = pulse.CFDtime;
					sameEventFitted[pulse.channel] = pulse.fittedtime;
					myCharges[pulse.channel].push_back(pulse.charge);
				//} Introduce commented section into code when analyzing a run taken with the internal LED
			}
		} else {
			continue;
		}
		
		std::vector<int> coincident_channels = {};
		for(Pulse pulse : sameEventPulse){ // find how many channels involved
			coincident_channels.push_back(pulse.channel);
		}
		std::sort(coincident_channels.begin(), coincident_channels.end());
		coincident_channels.erase(std::unique(coincident_channels.begin(), coincident_channels.end()), coincident_channels.end());
		
		
		for(int j : coincident_channels){ // time difference calculation relative to channel 2
			std::cout << "Channel " << j << " CFD Time Difference: " << sameEventCFD[j] - sameEventCFD[2] << " Fitted Time Difference "
					  << sameEventFitted[j] - sameEventFitted[2] << std::endl;
			h2[j]->Fill(sameEventCFD[j] - sameEventCFD[2]);
			h12[j]->Fill(sameEventFitted[j] - sameEventFitted[2]);
		}
		
		std::vector<double> lst_CFD = {};
		for(Pulse pulse : sameEventPulse){
			lst_CFD.push_back(pulse.CFDtime);
		}
		double meanCFD = std::accumulate(lst_CFD.begin(), lst_CFD.end(), 0.0) / lst_CFD.size();
		for(Pulse pulse : sameEventPulse){
			h3->Fill(pulse.CFDtime - meanCFD);
			h8[pulse.channel]->Fill(pulse.CFDtime - meanCFD);
		}
	}
	
	for(int j : myChannels){ // average CFD time difference per channel
		myTimeDiffMeansCFD[j] = h2[j]->GetMean();
		std::cout << "Channel " << j << " Mean CFD Time Difference: " << myTimeDiffMeansCFD[j] << " ns" << std::endl;
		myTimeDiffSTDCFD[j] = h2[j]->GetStdDev();
		std::cout << "Channel " << j << " CFD Time Difference STD: " << myTimeDiffSTDCFD[j] << std::endl;
	}
	
	for(int j : myChannels){ // average Fitted time difference per channel
		myTimeDiffMeansFitted[j] = h12[j]->GetMean();
		std::cout << "Channel " << j << " Mean Fitted Time Difference: " << myTimeDiffMeansFitted[j] << " ns" << std::endl;
		myTimeDiffSTDFitted[j] = h12[j]->GetStdDev();
		std::cout << "Channel " << j << " Fitted Time Difference STD: " << myTimeDiffSTDFitted[j] << std::endl;
		myChargeMeans[j] = std::accumulate(myCharges[j].begin(), myCharges[j].end(), 0.0) / myCharges[j].size();
		std::cout << "Channel " << j << " Mean Pulse Charge: " << myChargeMeans[j] << " PE" << std::endl;
	}
	
	
	// Print Pulse Time CFD Histograms
	for(int j : myChannels){
		TCanvas * c1 = new TCanvas("C1");
		
		char hist_title[100];
		sprintf(hist_title,"Channel %d CFD Pulse Times",j);
		h1[j]->SetTitle(hist_title);
		
		h1[j]->Draw();
		h1[j]->GetXaxis()->SetTitle("Pulse Times (ns)");
		h1[j]->GetYaxis()->SetTitle("Number of events");
		gStyle->SetOptFit(11);
		
		char png_name[100];
		sprintf(png_name,"pulse_time_CFD_%d_LED.png",j);
		c1->SaveAs(png_name);
	}
	
	// Print Pulse Time Difference Relative to Channel 2 Histograms
	for(int j : myChannels){
		TCanvas * c2 = new TCanvas("C2");
		
		char hist_title[100];
		sprintf(hist_title,"Channel %d CFD Pulse Time Differences",j);
		h2[j]->SetTitle(hist_title);
		
		h2[j]->Draw();
		h2[j]->GetXaxis()->SetTitle("Pulse Times (ns)");
		h2[j]->GetYaxis()->SetTitle("Number of events");
		gStyle->SetOptFit(11);
		
		char png_name[100];
		sprintf(png_name,"pulse_CFD_time_difference_CHANNEL2_%d_LED.png",j);
		c2->SaveAs(png_name);
	}
	
	// Print Pulse Time Differences Relative to Mean
	TCanvas * c3 = new TCanvas("C3");
	h3->Draw();
	h3->GetXaxis()->SetTitle("Time Difference Relative to Mean");
	h3->GetYaxis()->SetTitle("Number of Hits");
	//gStyle->SetOptStat(kFALSE);
	TLegend *leg3 = new TLegend(0.7,0.7,0.9,0.9);
	for(int j : myChannels){
		h8[j]->Draw("SAME");
		h8[j]->SetLineColor(j);
		
		char entry_title[100];
		sprintf(entry_title,"Channel %d Pulse Time Difference",j);
		leg3->AddEntry(h8[j],entry_title,"l");	
	}
	leg3->Draw();
	gStyle->SetOptFit(11);
	c3->SaveAs("pulse_time_difference_MEAN.png");
	
	// Print Standard Deviation vs Pulse Charge Graph
	TCanvas * c4 = new TCanvas("C4");
	TGraph * g4 = new TGraph(20, myChargeMeans, myTimeDiffSTDCFD);
	g4->SetTitle("Standard Deviation vs Mean Pulse Charge; Mean Pulse Charge (PE); Standard Deviation");
	gStyle->SetOptFit(11);
	g4->SetMarkerColor(4);
	g4->Draw("ap*");
	c4->SaveAs("std_deviation_vs_mean_pulse_charge_LED.png");
	
	// Print Pulse Time CFD - time
	for(int j : myChannels){
		TCanvas * c5 = new TCanvas("C5");
		
		char hist_title[100];
		sprintf(hist_title,"Channel %d Pulse Time Differences: CFD - time",j);
		h5[j]->SetTitle(hist_title);
		
		h5[j]->Draw();
		h5[j]->GetXaxis()->SetTitle("Pulse Time Difference (ns)");
		h5[j]->GetYaxis()->SetTitle("Number of Events");
		gStyle->SetOptFit(11);
		
		char png_name[100];
		sprintf(png_name,"pulse_time_difference_CFD-time_%d_LED.png",j);
		c5->SaveAs(png_name);
	}
	
	// Print Pulse Original Time Histograms
	for(int j : myChannels){
		TCanvas * c6 = new TCanvas("C6");
		
		char hist_title[100];
		sprintf(hist_title,"Channel %d Original Pulse Times",j);
		h6[j]->SetTitle(hist_title);
		
		h6[j]->Draw();
		h6[j]->GetXaxis()->SetTitle("Pulse Times (ns)");
		h6[j]->GetYaxis()->SetTitle("Number of events");
		gStyle->SetOptFit(11);
		
		char png_name[100];
		sprintf(png_name,"pulse_time_ORIGINAL_%d_LED.png",j);
		c6->SaveAs(png_name);
	}
	
	// Print Pulse Time CFD vs Original
	for(int j : myChannels){
		TCanvas * c7 = new TCanvas("C7","C7");
		//gStyle->SetOptStat(kFALSE);
		h1[j]->Draw();
		h6[j]->Draw("SAME");
		h6[j]->SetLineColor(kRed);
		TLegend *leg7 = new TLegend(0.7,0.7,0.9,0.9);
		leg7->AddEntry(h1[j],"Pulse Time CFD","l");
		leg7->AddEntry(h6[j],"Pulse Time Original","l");
		leg7->Draw();
		h1[j]->GetXaxis()->SetTitle("Pulse Times (ns)");
		h1[j]->GetYaxis()->SetTitle("Number of Events");
		gStyle->SetOptFit(11);
		
		char hist_title[100];
		sprintf(hist_title,"Channel %d Pulse Time CFD vs Original",j);
		h1[j]->SetTitle(hist_title);
		
		char png_name[100];
		sprintf(png_name,"pulse_time_CFDvsOriginal_%d.png",j);
		c7->SaveAs(png_name);
	}
	
	// Print Pulse Fitted Time Histograms (LED Runs Only)
	for(int j : myChannels){
		TCanvas * c9 = new TCanvas("C9");
		
		char hist_title[100];
		sprintf(hist_title,"Channel %d Fitted Pulse Times",j);
		h9[j]->SetTitle(hist_title);
		
		h9[j]->Draw();
		h9[j]->GetXaxis()->SetTitle("Pulse Times (ns)");
		h9[j]->GetYaxis()->SetTitle("Number of events");
		gStyle->SetOptFit(11);
		
		char png_name[100];
		sprintf(png_name,"pulse_time_FITTED_%d_LED.png",j);
		c9->SaveAs(png_name);
	}
	
	//Print Pulse Time CFD vs Fitted Time
	for(int j : myChannels){
		TCanvas * c10 = new TCanvas("C10","C10");
		gStyle->SetOptStat(kFALSE);
		h1[j]->Draw();
		h9[j]->Draw("SAME");
		h9[j]->SetLineColor(kRed);
		TLegend *leg10 = new TLegend(0.7,0.7,0.9,0.9);
		leg10->AddEntry(h1[j],"Pulse Time CFD","l");
		leg10->AddEntry(h9[j],"Pulse Time Fitted","l");
		leg10->Draw();
		h1[j]->GetXaxis()->SetTitle("Pulse Times (ns)");
		h1[j]->GetYaxis()->SetTitle("Number of Events");
		gStyle->SetOptFit(11);
		
		char hist_title[100];
		sprintf(hist_title,"Channel %d Pulse Time CFD vs Fitted",j);
		h1[j]->SetTitle(hist_title);
		
		char png_name[100];
		sprintf(png_name,"pulse_time_CFDvsFitted_%d.png",j);
		c10->SaveAs(png_name);
	}
	
	// Print Pulse Time CFD - Fitted Time
	for(int j : myChannels){
		TCanvas * c11 = new TCanvas("C11");
		
		char hist_title[100];
		sprintf(hist_title,"Channel %d Pulse Time Differences: CFD - Fitted Time",j);
		h11[j]->SetTitle(hist_title);
		
		h11[j]->Draw();
		h11[j]->GetXaxis()->SetTitle("Pulse Time Difference (ns)");
		h11[j]->GetYaxis()->SetTitle("Number of Events");
		gStyle->SetOptFit(11);
		
		char png_name[100];
		sprintf(png_name,"pulse_time_difference_CFD-fitted_%d.png",j);
		c11->SaveAs(png_name);
	}
	
	// Print Pulse Fitted Time Difference Relative to Channel 2
	for(int j : myChannels){
		TCanvas * c12 = new TCanvas("C12");
		
		char hist_title[100];
		sprintf(hist_title,"Channel %d Fitted Pulse Time Differences",j);
		h12[j]->SetTitle(hist_title);
		
		h12[j]->Draw();
		h12[j]->GetXaxis()->SetTitle("Pulse Times (ns)");
		h12[j]->GetYaxis()->SetTitle("Number of events");
		gStyle->SetOptFit(11);
		
		char png_name[100];
		sprintf(png_name,"pulse_FITTED_time_difference_CHANNEL2_%d_LED.png",j);
		c12->SaveAs(png_name);
	}
	
	// Print Pulse Time CFD Differences vs Fitted Differences
	for(int j : myChannels){
		TCanvas * c13 = new TCanvas("C13","C13");
		//gStyle->SetOptStat(kFALSE);
		h2[j]->Draw();
		h12[j]->Draw("SAME");
		h12[j]->SetLineColor(kRed);
		TLegend *leg13 = new TLegend(0.7,0.7,0.9,0.9);
		leg13->AddEntry(h2[j],"Pulse Time CFD Differences","l");
		leg13->AddEntry(h12[j],"Pulse Time Fitted Differences","l");
		leg13->Draw();
		h2[j]->GetXaxis()->SetTitle("Pulse Time Differences (ns)");
		h2[j]->GetYaxis()->SetTitle("Number of Events");
		gStyle->SetOptFit(11);
		
		char hist_title[100];
		sprintf(hist_title,"Channel %d Pulse Time CFD vs Fitted Differences",j);
		h2[j]->SetTitle(hist_title);
		
		char png_name[100];
		sprintf(png_name,"pulse_time_CFDvsFitted_DIFFERENCES_%d.png",j);
		c13->SaveAs(png_name);
	}
	
	// Print Pulse Time CFD Histograms
	for(int j : myChannels){
		TCanvas * c14 = new TCanvas("C14");
		
		char hist_title[100];
		sprintf(hist_title,"Channel %d Pulse Charges",j);
		h14[j]->SetTitle(hist_title);
		
		h14[j]->Draw();
		h14[j]->GetXaxis()->SetTitle("Pulse Charge (PE)");
		h14[j]->GetYaxis()->SetTitle("Number of events");
		gStyle->SetOptFit(11);
		
		char png_name[100];
		sprintf(png_name,"pulse_charges_%d_LED.png",j);
		c14->SaveAs(png_name);
	}
	
	std::cout << "Number of Events: " << tt[0]->GetEntries() << std::endl;
	
	fin->Close();
	return 0;
}
	
