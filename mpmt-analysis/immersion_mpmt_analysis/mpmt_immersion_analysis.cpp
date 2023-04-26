// Hassan Elshabasy Immersion mPMT Analysis
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
		double time;
        Pulse(int aEvent, int aChannel, double aHeight, double aCharge, double aTime){
            event = aEvent;
			channel = aChannel;
            height = aHeight;
			charge = aCharge;
			time = aTime;	
        }
};

bool notSame(Pulse pulse1, Pulse pulse2){
	if((pulse1.channel == pulse2.channel) && (pulse1.time == pulse2.time)){
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
	TH1F * h1[20]; //pulse heights hist
	TH1F * h2[20]; //pulse charges hist
    
	double liveTime = 0;
	
	TH1F *h3 = new TH1F("hit_coincidences_2", "Hit Coincidences 2 Channels",200,0,17);
	TH1F *h4 = new TH1F("hit_coincidences_4", "Hit Coincidences 4 Channels",200,0,17);
	TH1F *h5 = new TH1F("coincident_heights_2", "Coincidences Heights 2 Channels",200,0,200*0.48828125);
	TH1F *h6 = new TH1F("coincident_heights_4", "Coincidences Heights 4 Channels",200,0,200*0.48828125);
	TH1F *h7 = new TH1F("coincidence_charges_max_2", "Max Coincidence Charges 2 Channels",100,0,30*0.48828125);
	TH1F *h8 = new TH1F("coincidence_charges_max_4", "Max Coincidence Charges 4 Channels",100,0,200*0.48828125);
	TH1F *h9 = new TH1F("coincidence_charges_sum_2", "Sum Coincidence Charges 2 Channels",100,0,30*0.48828125);
	TH1F *h10 = new TH1F("coincidence_charges_sum_4", "Sum Coincidence Charges 4 Channels",100,0,200*0.48828125);
	TH1F *h11 = new TH1F("coincidence_time_difference_2", "Coincidence Time Differences 2 Channels",200,0,100*0.48828125);
	TH1F *h12 = new TH1F("coincidence_time_difference_4", "Coincidence Time Differences 4 Channels",200,0,100*0.48828125);
	TH1F *h13 = new TH1F("coincident_times_2", "Coincident Times 2 Channels",100,0,17000*0.48828125);
	TH1F *h14 = new TH1F("coincident_times_4", "Coincident Times 4 Channels",100,0,17000*0.48828125);
	TH1F *h15 = new TH1F("coincident_channels_2", "Coincident Channels 2",200,0,20);
	TH1F *h16 = new TH1F("coincident_channels_4", "Coincident Channels 4",200,0,20);
	TH1F *h17 = new TH1F("coincidence_charges_ratio_2", "Ratio Coincidence Charges 2 Channels",100,0,2*0.48828125);
	TH1F *h18 = new TH1F("coincidence_charges_ratio_4", "Ratio Coincidence Charges 4 Channels",100,0,2*0.48828125);
	
	TH2F *h19 = new TH2F("coincident_charges_sumvratio_2", "Sum vs Ratio Coincident Charges 2 Channels",
	50,0,30*0.48828125,50,0,2*0.48828125);
	TH2F *h20 = new TH2F("coincident_charges_sumvratio_4", "Sum vs Ratio Coincident Charges 4 Channels",
	50,0,200*0.48828125,50,0,2*0.48828125);
	
	int twoNum = 0; // number of 2-fold events
	int fourNum = 0; // number of 4-fold events
	
	int time_diff = 100;
	
    for(int j = 0; j < 20; j++){ // To initialize TTree and branches
	
		char branch_name[100];
		sprintf(branch_name,"ptfanalysis%i",j);
		tt[j] = (TTree*)fin->Get(branch_name);
			   
		wf[j] = new WaveformFitResult;
		if(tt[j]) wf[j]->SetBranchAddresses( tt[j] );
		
		h1[j] = new TH1F("pulse_height", "Pulse Height",200,0,200*0.48828125);
		h2[j] = new TH1F("pulse_charge", "Pulse Charge",200,0,200*0.48828125/7.7);
    }
	
	int myChannels[16] = {0, 1, 2, 3, 4, 5, 6, 7, 10, 11, 12, 13, 14, 15, 16, 17};
	
	for(int i = 0; i < tt[0]->GetEntries()-1; i++){ // loop over events
		liveTime += 8192*pow(10,-9); // in seconds
		
		std::vector<Pulse> myPulses = {};
		std::vector<Pulse> coincident_pulses = {};
		std::vector<int> coincident_channels = {};
		std::vector<Pulse> coincident_pulses_2 = {};
		std::vector<Pulse> coincident_pulses_4 = {};
		std::vector<double> coincident_charges_2 = {};
		std::vector<double> coincident_charges_4 = {};
		std::vector<double> coincident_times_2 = {};
		std::vector<double> coincident_times_4 = {};
		
		for(int j : myChannels){ // loop over channels
			tt[j]->GetEvent(i);
			
			for(int k = 0; k < wf[j]->numPulses; k++){ // loop over each pulse
				
				std::cout << "Pulses found: " << wf[j]->numPulses << " Event: " << i << " Channel: " << j << std::endl;
				
				// pulse height
				double pulse_height = wf[j]->pulseCharges[k]*1000.0;
				h1[j]->Fill(pulse_height);
				
				// pulse charge in photoelectron assume 7.7mV
				double pulse_charge = pulse_height / 7.7;
				h2[j]->Fill(pulse_charge);
				
				// pulse time
				double pulse_time = wf[j]->pulseTimesCFD[k];
				
				// Time Calibration (not yet completed)
				/*
				switch(j){
				case 0:
					pulse_time -= -1.4766;
					break;
				case 1:
					pulse_time -= 0.0;
					break;
				case 2:
					pulse_time -= -0.245132;
					break;
				case 3:
					pulse_time -= -0.1279973;
					break;
				case 4:
					pulse_time -= 1.57659;
					break;
				case 5:
					pulse_time -= -0.83196;
					break;
				case 6:
					pulse_time -= -0.994029;
					break;
				case 7:
					pulse_time -= 1.46915;
					break;
				case 10:
					pulse_time -= 15.5602;
					break;
				case 11:
					pulse_time -= 9.65071;
					break;
				case 12:
					pulse_time -= 8.23973;
					break;
				case 13:
					pulse_time -= 3.2517;
					break;
				case 14:
					pulse_time -= 11.2103;
					break;
				case 15:
					pulse_time -= 5.525648;
					break;
				case 16:
					pulse_time -= 8.21802;
					break;
				case 17:
					pulse_time -= -2.47894;
					break;
				}*/
				
				Pulse pulse(i, j, pulse_height, pulse_charge, pulse_time);
				myPulses.push_back(pulse);
				
			}
		}
		
		if(myPulses.size() < 2){
			continue;
		}
		
		for(Pulse pulse_check : myPulses){
			for(Pulse pulse : myPulses){
				if((abs(pulse_check.time - pulse.time) < time_diff) && (notSame(pulse_check, pulse))){
					coincident_pulses.push_back(pulse_check);
					std::cout << "Channel: " << pulse_check.channel << " Time: " << pulse_check.time << std::endl;
					break;
				}
			}
		}
			
		for(Pulse pulse : coincident_pulses){ // find how many channels involved
			coincident_channels.push_back(pulse.channel);
		}
		
		std::sort(coincident_channels.begin(), coincident_channels.end());
		coincident_channels.erase(std::unique(coincident_channels.begin(), coincident_channels.end()), coincident_channels.end());
		
		if(coincident_channels.size() >= 2){// 2-fold coincident events
			coincident_pulses_2 = coincident_pulses;
			h3->Fill(coincident_channels.size());
		} else {
			continue;
		}
		
		for(Pulse pulse : coincident_pulses_2){ // 2-fold
			coincident_charges_2.push_back(pulse.charge); // analyze 2-fold pulse charges
			coincident_times_2.push_back(pulse.time); // analyze 2-fold pulse times
			h5->Fill(pulse.height);
			h13->Fill(pulse.time);
			h15->Fill(pulse.channel);
		}
		
		double max_2 = *max_element(coincident_charges_2.begin(), coincident_charges_2.end());
		double sum_2 = std::accumulate(coincident_charges_2.begin(), coincident_charges_2.end(), 0.0);
		double ratio_2 = max_2 / sum_2;
		h7->Fill(max_2); // max charge
		h9->Fill(sum_2); // sum charge
		h17->Fill(ratio_2);
		h19->Fill(sum_2,ratio_2);
		
		std::sort(coincident_times_2.begin(), coincident_times_2.end());

		while(coincident_times_2.size() >= 2){
			if(coincident_times_2[1] - coincident_times_2[0] <= time_diff){
				h11->Fill(coincident_times_2[1] - coincident_times_2[0]);
				coincident_times_2.erase(coincident_times_2.begin()+1);
			} else {
				coincident_times_2.erase(coincident_times_2.begin());
			}
		}
		
		std::cout << "2-fold event " << i << std::endl;
		twoNum++;
		
		if(coincident_channels.size() >= 4){// 4-fold coincident events
			coincident_pulses_4 = coincident_pulses;
			h4->Fill(coincident_channels.size());
		}
		
		if(coincident_pulses_4.size() == 0){
			continue;
		}
		
		for(Pulse pulse : coincident_pulses_4){ // 4-fold
			coincident_charges_4.push_back(pulse.charge); // analyze 4-fold pulse charges
			coincident_times_4.push_back(pulse.time); // analyze 4-fold pulse times
			h6->Fill(pulse.height);
			h14->Fill(pulse.time);
			h16->Fill(pulse.channel);
		}
		
		double max_4 = *max_element(coincident_charges_4.begin(), coincident_charges_4.end());
		double sum_4 = std::accumulate(coincident_charges_4.begin(), coincident_charges_4.end(), 0.0);
		double ratio_4 = max_4 / sum_4;
		h8->Fill(max_4); // max charge
		h10->Fill(sum_4); // sum charge
		h18->Fill(ratio_4);
		h20->Fill(sum_4,ratio_4);
		
		std::sort(coincident_times_4.begin(), coincident_times_4.end());

		// TIME CUT
		while(coincident_times_4.size() >= 2){
			if(coincident_times_4[1] - coincident_times_4[0] <= time_diff){
				h12->Fill(coincident_times_4[1] - coincident_times_4[0]);
				coincident_times_4.erase(coincident_times_4.begin()+1);
			} else {
				coincident_times_4.erase(coincident_times_4.begin());
			}
		}
		
		std::cout << "4-fold event " << i << std::endl;
		fourNum++;		
	}
	
	
	// Print Pulse Height Histograms
	for(int j : myChannels){
		TCanvas * c1 = new TCanvas("C1");
		
		char hist_title[100];
		sprintf(hist_title,"Channel %i Pulse Height",j);
		h1[j]->SetTitle(hist_title);
		
		h1[j]->Draw();
		h1[j]->GetXaxis()->SetTitle("Pulse height (mV)");
		h1[j]->GetYaxis()->SetTitle("Number of events");
		h1[j]->Fit("gaus","","",4,22);
		gStyle->SetOptFit(11);
		
		char png_name[100];
		sprintf(png_name,"mpmt_pulse_height_%i.png",j);
		c1->SaveAs(png_name);
	}

	// Print Pulse Charges Histograms
	for(int j : myChannels){
		TCanvas * c2 = new TCanvas("C2");
		
		char hist_title[100];
		sprintf(hist_title,"Channel %i Pulse Charge",j);
		h2[j]->SetTitle(hist_title);
		
		h2[j]->Draw();
		h2[j]->GetXaxis()->SetTitle("Pulse Charge (photoelectrons)");
		h2[j]->GetYaxis()->SetTitle("Number of events");
		h2[j]->Fit("gaus","","",4/7.7,22/7.7);
		gStyle->SetOptFit(11);
		
		char png_name[100];
		sprintf(png_name,"mpmt_pulse_charge_%i.png",j);
		c2->SaveAs(png_name);
	}*/
	
	// Print Hit Coincident Histograms
	TCanvas *c3 = new TCanvas("C3");
	h3->Draw();
	h3->GetXaxis()->SetTitle("Number of Channels Hit");
	h3->GetYaxis()->SetTitle("Number of Hits");
	h3->Fit("gaus","","",1.9,17);
	gStyle->SetOptFit(11);
	c3->SaveAs("hit_coincidences_2.png");
	
	TCanvas *c4 = new TCanvas("C4");
	h4->Draw();
	h4->GetXaxis()->SetTitle("Number of Channels Hit");
	h4->GetYaxis()->SetTitle("Number of Hits");
	h4->Fit("gaus","","",3.9,17);
	gStyle->SetOptFit(11);
	c4->SaveAs("hit_coincidences_4.png");
	
	// Print coincidence height
	TCanvas *c5 = new TCanvas("C5");
	h5->Draw();
	h5->GetXaxis()->SetTitle("Pulse Height (mV)");
	h5->GetYaxis()->SetTitle("Number of Hits");
	//h5->Fit("gaus","","",4,50);
	gStyle->SetOptFit(11);
	c5->SaveAs("coincidence_heights_2.png");
	
	TCanvas *c6 = new TCanvas("C6");
	h6->Draw();
	h6->GetXaxis()->SetTitle("Pulse Height (mV)");
	h6->GetYaxis()->SetTitle("Number of Hits");
	//h6->Fit("gaus","","",4,50);
	gStyle->SetOptFit(11);
	c6->SaveAs("coincidence_heights_4.png");
	
	// Print Max Coincident Charge Histogram
	TCanvas *c7 = new TCanvas("C7");
	h7->Draw();
	h7->GetXaxis()->SetTitle("Pulse Charge (PE)");
	h7->GetYaxis()->SetTitle("Number of Hits");
	h7->Fit("gaus","","",0.5,3);
	gStyle->SetOptFit(11);
	c7->SaveAs("coincidence_charges_MAX_2.png");
	
	TCanvas *c8 = new TCanvas("C8");
	h8->Draw();
	h8->GetXaxis()->SetTitle("Pulse Charge (PE)");
	h8->GetYaxis()->SetTitle("Number of Hits");
	//h8->Fit("gaus","","",0.5,3);
	gStyle->SetOptFit(11);
	c8->SaveAs("coincidence_charges_MAX_4.png");
	
	// Print Sum Coincident Charge Histogram
	TCanvas *c9 = new TCanvas("C9");
	h9->Draw();
	h9->GetXaxis()->SetTitle("Pulse Charge (PE)");
	h9->GetYaxis()->SetTitle("Number of Hits");
	h9->Fit("gaus","","",1,6);
	gStyle->SetOptFit(11);
	c9->SaveAs("coincidence_charges_SUM_2.png");
	
	TCanvas *c10 = new TCanvas("C10");
	h10->Draw();
	h10->GetXaxis()->SetTitle("Pulse Charge (PE)");
	h10->GetYaxis()->SetTitle("Number of Hits");
	//h10->Fit("gaus","","",0.5,3);
	gStyle->SetOptFit(11);
	c10->SaveAs("coincidence_charges_SUM_4.png");
	
	// Print coincidence time difference
	TCanvas *c11 = new TCanvas("C11");
	h11->Draw();
	h11->GetXaxis()->SetTitle("Pulse Time Differences (ns)");
	h11->GetYaxis()->SetTitle("Number of Hits");
	gStyle->SetOptFit(11);
	//gPad->SetLogy();
	c11->SaveAs("coincidence_time_differences_2.png");
	
	TCanvas *c12 = new TCanvas("C12");
	h12->Draw();
	h12->GetXaxis()->SetTitle("Pulse Time Differences (ns)");
	h12->GetYaxis()->SetTitle("Number of Hits");
	gStyle->SetOptFit(11);
	//gPad->SetLogy();
	c12->SaveAs("coincidence_time_differences_4.png");
	
	// Print Coincidence Times Historgrams
	TCanvas *c13 = new TCanvas("C13");
	h13->Draw();
	h13->GetXaxis()->SetTitle("Pulse Time (ns)");
	h13->GetYaxis()->SetTitle("Number of Hits");
	gStyle->SetOptFit(11);
	c13->SaveAs("coincident_times_2.png");
	
	TCanvas *c14 = new TCanvas("C14");
	h14->Draw();
	h14->GetXaxis()->SetTitle("Pulse Time (ns)");
	h14->GetYaxis()->SetTitle("Number of Hits");
	gStyle->SetOptFit(11);
	c14->SaveAs("coincident_times_4.png");
	
	// Print Coincident Channels Histograms
	TCanvas *c15 = new TCanvas("C15");
	h15->Draw();
	h15->GetXaxis()->SetTitle("Channels Number");
	h15->GetYaxis()->SetTitle("Number of Coincident Hits");
	//h15->Fit("gaus","","",1.9,17);
	gStyle->SetOptFit(11);
	c15->SaveAs("coincident_channels_2.png");
	
	TCanvas *c16 = new TCanvas("C16");
	h16->Draw();
	h16->GetXaxis()->SetTitle("Channels Number");
	h16->GetYaxis()->SetTitle("Number of Coincident Hits");
	//h16->Fit("gaus","","",1.9,17);
	gStyle->SetOptFit(11);
	c16->SaveAs("coincident_channels_4.png");
	
	TCanvas *c17 = new TCanvas("C17");
	h17->Draw();
	h17->GetXaxis()->SetTitle("Max/Sum");
	h17->GetYaxis()->SetTitle("Number of Occurance");
	//h17->Fit("gaus","","",1,6);
	gStyle->SetOptFit(11);
	c17->SaveAs("coincidence_charges_RATIO_2.png");
	
	TCanvas *c18 = new TCanvas("C18");
	h18->Draw();
	h18->GetXaxis()->SetTitle("Max/Sum");
	h18->GetYaxis()->SetTitle("Number of Occurance");
	//h18->Fit("gaus","","",0.5,3);
	gStyle->SetOptFit(11);
	c18->SaveAs("coincidence_charges_RATIO_4.png");
	
	// Print Sum vs Ratio Charge Histograms
	TCanvas *c19 = new TCanvas("C19");
	h19->Draw("colz");
	h19->SetStats(0);
	h19->SetContour(1000);
	h19->GetXaxis()->SetTitle("Sum Charge (photoelectrons)");
	h19->GetYaxis()->SetTitle("Charge Ratio (Max/Sum)");
	//h19->Fit("gaus","","",0.5,3);
	gStyle->SetOptFit(11);
	gStyle->SetPalette(kRainBow);
	c19->SaveAs("coincidence_charges_SUMVRATIO_2.png");
	
	TCanvas *c20 = new TCanvas("C20");
	h20->Draw("colz");
	h20->SetStats(0);
	h20->SetContour(1000);
	h20->GetXaxis()->SetTitle("Sum Charge (photoelectrons)");
	h20->GetYaxis()->SetTitle("Charge Ratio (Max/Sum)");
	//h20->Fit("gaus","","",0.5,3);
	gStyle->SetOptFit(11);
	gStyle->SetPalette(kRainBow);
	c20->SaveAs("coincidence_charges_SUMVRATIO_4.png");
	
	
	std::cout << "Number of Events: " << tt[0]->GetEntries() << std::endl;
	std::cout << "Live Time: " << liveTime << " seconds" << std::endl;
	
	std::cout << "Number of 2-fold coincident events: " << twoNum << std::endl;
	std::cout << "Rate of 2-fold coincident events: " << twoNum / liveTime << " hits/second" << std::endl;
	
	std::cout << "Number of 4-fold coincident events: " << fourNum << std::endl;
	std::cout << "Rate of 4-fold coincident events: " << fourNum / liveTime << " hits/second" << std::endl;
	
	
	fin->Close();
	return 0;
}
	