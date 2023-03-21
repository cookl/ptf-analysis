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
		double charge;
		double CFDtime;
		double time;
        Pulse(int aEvent, int aChannel, double aHeight, double aCharge, double aCFDTime, double aTime){
            event = aEvent;
			channel = aChannel;
            height = aHeight;
			charge = aCharge;
			CFDtime = aCFDTime;	
			time = aTime;
        }
};

bool notSame(Pulse pulse1, Pulse pulse2){
	if((pulse1.channel == pulse2.channel) && (pulse1.CFDtime == pulse2.CFDtime)){
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
	TH1F * h1[20]; // pulse times CFD per channel
	TH1F * h2[20]; // pulse time differences relative to channel 1
	TH1F * h5[20]; // pule time difference between CFDtime and time
	TH1F * h6[20]; // original pulse times per channel
	TH1F * h8[20]; // pulse time difference relative to mean per channel
	
    for(int j = 0; j < 20; j++){ // To initialize TTree and branches
	
		char branch_name[100];
		sprintf(branch_name,"ptfanalysis%i",j);
		tt[j] = (TTree*)fin->Get(branch_name);
			   
		wf[j] = new WaveformFitResult;
		if(tt[j]) wf[j]->SetBranchAddresses( tt[j] );
		
		h1[j] = new TH1F("pulse_time_CFD", "Pulse Time CFD",200,2300,2400);
		h2[j] = new TH1F("pulse_time_difference", "Pulse Time Difference",100,-25,25); // relative to channel 1
		h5[j] = new TH1F("pulse_time_cfd-time", "Pulse Time Difference Between CFDtime and time", 100, -25, 25);
		h6[j] = new TH1F("pulse_time_original", "Pulse Time Original",200,2300,2400);
		h8[j] = new TH1F("pulse_time_difference_mean_channel", "Pulse Time Difference Relative to Mean per Channel", 100, -25, 25);
		
    }
	
	TH1F * h3 = new TH1F("pulse_time_difference_mean", "Pulse Time Difference Relative to Mean", 100, -25, 25);
	
	int myChannels[] = {1, 2, 3, 4, 6, 7, 12, 13, 14};
	
	// all active channels {0, 1, 2, 3, 4, 5, 6, 7, 10, 11, 12, 13, 14, 15, 16, 17};
	// none-saturated channels {1, 2, 3, 4, 6, 7, 12, 13, 14};
	
	double myTimeDiffMeans[20];
	
	double myTimeDiffSTD[20];
	
	double myChargeMeans[20];
	std::vector<double> myCharges[20] = {};
	
	for(int i = 0; i < tt[0]->GetEntries()-1; i++){ // loop over events
		
		std::cout << "Event: " << i << std::endl;
		
		double sameEventCFD[20];
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
				
				Pulse pulse(i, j, pulse_height, pulse_charge, pulse_CFDtime, pulse_time);
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
			h1[j]->Fill(LED.CFDtime);
			h5[j]->Fill(LED.CFDtime - LED.time);
			h6[j]->Fill(LED.time);
			
			sameEventCFD[j] = LED.CFDtime;
			sameEventPulse.push_back(LED);
			myCharges[j].push_back(LED.charge);
		}
		
		if(sameEventPulse.size() == 0){
			continue;
		}
		
		for(int j : myChannels){ // time difference calculation relative to channel 1
			std::cout << "Channel " << j << " Time Difference: " << sameEventCFD[j] - sameEventCFD[1] << std::endl;
			h2[j]->Fill(sameEventCFD[j] - sameEventCFD[1]);
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
	
	for(int j : myChannels){ // average time difference per channel
		myTimeDiffMeans[j] = h2[j]->GetMean();
		std::cout << "Channel " << j << " Mean Time Difference: " << myTimeDiffMeans[j] << " ns" << std::endl;
		myTimeDiffSTD[j] = h2[j]->GetStdDev();
		std::cout << "Channel " << j << " Time Difference STD: " << myTimeDiffSTD[j] << std::endl;
		myChargeMeans[j] = std::accumulate(myCharges[j].begin(), myCharges[j].end(), 0.0) / myCharges[j].size();
		std::cout << "Channel " << j << " Mean Pulse Charge: " << myChargeMeans[j] << " PE" << std::endl;
	}
	
	
	// Print Pulse Time Histograms
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
	
	// Print Pulse Time Difference Relative to Channel 1 Histograms
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
		sprintf(png_name,"pulse_time_difference_CHANNEL1_%d_LED.png",j);
		c2->SaveAs(png_name);
	}
	
	// Print Pulse Time Differences Mean
	TCanvas * c3 = new TCanvas("C3");
	h3->Draw();
	h3->GetXaxis()->SetTitle("Time Difference Relative to Mean");
	h3->GetYaxis()->SetTitle("Number of Hits");
	gStyle->SetOptStat(kFALSE);
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
	TGraph * g4 = new TGraph(20, myChargeMeans, myTimeDiffSTD);
	g4->SetTitle("Standard Deviation vs Mean Pulse Charge; Standard Deviation; Mean Pulse Charge (PE)");
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
		gStyle->SetOptStat(kFALSE);
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
	
	
	
		
	/*
	std::unique_ptr<TFile> myFile_978( TFile::Open("run_978_LED.root", "RECREATE") );
	myFile_978->WriteObject(h2[1], "h21_978");
	myFile_978->WriteObject(h2[2], "h22_978");
	myFile_978->WriteObject(h2[3], "h23_978");
	myFile_978->WriteObject(h2[4], "h24_978");
	myFile_978->WriteObject(h2[6], "h26_978");
	myFile_978->WriteObject(h2[7], "h27_978");
	myFile_978->WriteObject(h2[12], "h212_978");
	myFile_978->WriteObject(h2[13], "h213_978");
	myFile_978->WriteObject(h2[14], "h214_978");
	std::cout << "File run_978_LED.root is saved" << std::endl;*/
	/*
	std::unique_ptr<TFile> myFile_980( TFile::Open("run_980_LED.root", "RECREATE") );
	myFile_980->WriteObject(h2[1], "h21_980");
	myFile_980->WriteObject(h2[2], "h22_980");
	myFile_980->WriteObject(h2[3], "h23_980");
	myFile_980->WriteObject(h2[4], "h24_980");
	myFile_980->WriteObject(h2[6], "h26_980");
	myFile_980->WriteObject(h2[7], "h27_980");
	myFile_980->WriteObject(h2[12], "h212_980");
	myFile_980->WriteObject(h2[13], "h213_980");
	myFile_980->WriteObject(h2[14], "h214_980");
	std::cout << "File run_980_LED.root is saved" << std::endl;*/
	
	/*
	TFile *f_978 = new TFile("run_978_LED.root");
	TFile *f_980 = new TFile("run_980_LED.root");
	TH1F* h_978[100];
	TH1F* h_980[100];
	for(int j : myChannels){
		char hist_978[100];
		sprintf(hist_978,"h2%d_978",j);
		h_978[j] = (TH1F*)f_978->Get(hist_978);
		
		char hist_980[100];
		sprintf(hist_980,"h2%d_980",j);
		h_980[j] = (TH1F*)f_980->Get(hist_980);
		
		TCanvas *c5 = new TCanvas("C5","C5");
		gStyle->SetOptStat(kFALSE);
		h_978[j]->Draw();
		h_980[j]->Draw("SAME");
		h_980[j]->SetLineColor(kRed);
		TLegend *leg5 = new TLegend(0.7,0.7,0.9,0.9);
		leg5->AddEntry(h_978[j],"Run 978: No Water","l");
		leg5->AddEntry(h_980[j],"Run 980: Extra Water","l");
		leg5->Draw();
		h_978[j]->GetXaxis()->SetTitle("Time Difference (ns)");
		h_978[j]->GetYaxis()->SetTitle("Number of Events");
		gStyle->SetOptFit(11);
		
		char hist_title[100];
		sprintf(hist_title,"Channel %d Pulse Time Differences",j);
		h_978[j]->SetTitle(hist_title);
		
		char png_name[100];
		sprintf(png_name,"pulse_time_difference_COMBINED_%d.png",j);
		c5->SaveAs(png_name);
	}*/
	
	std::cout << "Number of Events: " << tt[0]->GetEntries() << std::endl;
	
	fin->Close();
	return 0;
}
	
