// Ashley/Thomas simple pulse height histogram
// 2020-01-20


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
//using namespace std;

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
	
	std::vector<int> hit_events = {}; //events that contain pulses
	std::vector<std::vector<double>> hit_heights = {}; //pulse heights of all hits
	std::vector<std::vector<double>> hit_charges = {}; //pulse charges of all hits
	std::vector<std::vector<double>> hit_times = {}; //pulse times of all hits
	std::vector<std::vector<double>> coincident_times_2 = {}; //2-fold pulse times of coincident hits
	std::vector<std::vector<double>> coincident_times_4 = {}; //4-fold pulse times of coincident hits
    
	double liveTime = 0;
	
    for(int j = 0; j < 20; j++){ // To initialize TTree and branchs
	
		char branch_name[100];
		sprintf(branch_name,"ptfanalysis%i",j);
		tt[j] = (TTree*)fin->Get(branch_name);
			   
		wf[j] = new WaveformFitResult;
		if(tt[j]) wf[j]->SetBranchAddresses( tt[j] );  // A bit of problem
		
		h1[j] = new TH1F("pulse_height", "Pulse Height",200,0,200*0.48828125);
		h2[j] = new TH1F("pulse_charge", "Pulse Charge",200,0,200*0.48828125/7.7);
    }
	
	int myChannels[16] = {0, 1, 2, 3, 4, 5, 6, 7, 10, 11, 12, 13, 14, 15, 16, 17};
	
	for(int i = 0; i < tt[0]->GetEntries()-1; i++){ // loop over events
		
		liveTime += 8192*pow(10,-9); //in seconds
		
		for(int j : myChannels){ // loop over channels
			tt[j]->GetEvent(i);
			
			if(wf[j]->numPulses > 0){
				hit_events.push_back(i);
			}
			
			for(int k = 0; k < wf[j]->numPulses; k++){ // loop over each pulse
				
				std::cout << "Number of pulses found: " << wf[j]->numPulses << 
				" Event Number: " << i << " Channel Number: " << j << std::endl;
				
				std::cout << "Pulse " << k << " has pulse height " << wf[j]->pulseCharges[k]*1000.0 
				<< "mV " << wf[j]->pulseTimes[k] << "ns" 
				<< std::endl;
				
				// pulse height
				double pulse_height = wf[j]->pulseCharges[k]*1000.0;
				h1[j]->Fill(pulse_height);
				
				// pulse charge in photoelectron assume 7.7mV
				double pulse_charge = pulse_height / 7.7;
				h2[j]->Fill(pulse_charge);
				
				// pulse time
				double pulse_time = wf[j]->pulseTimes[k];
				
				if(wf[j]->numPulses != 0){
					double ind = i;
					hit_heights.push_back(std::vector<double> {ind, pulse_height});
					hit_charges.push_back(std::vector<double> {ind, pulse_charge});
					hit_times.push_back(std::vector<double> {ind, pulse_time});
				}
				
			}
			
		}
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
	}
	
	//************************************************************************************************************
	// Fill hit coincidences histograms
	TH1F *h3 = new TH1F("hit_coincidences_2", "Hit Coincidences 2 Channels",200,0,17);
	TH1F *h4 = new TH1F("hit_coincidences_4", "Hit Coincidences 4 Channels",200,0,17);
	std::vector<int> coincident_events_2 = {};
	std::vector<int> coincident_events_4 = {};
	
	int count = 1;
	while(hit_events.size() > 1){
		if(hit_events[0] == hit_events[1]){
			count++;
		} else {
			if(count >= 2){
				h3->Fill(count);
				std::cout << "Event " << hit_events[0] << " had " << count << " coincidences" << std::endl;
				coincident_events_2.push_back(hit_events[0]);
			}
			if(count >= 4){
				h4->Fill(count);
				std::cout << "Event " << hit_events[0] << " had " << count << " coincidences" << std::endl;
				coincident_events_4.push_back(hit_events[0]);
			}
			count = 1;
		}
		hit_events.erase(hit_events.begin());
	}
	
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
	
	//************************************************************************************************************
	// Fill Coincident Height Histogram
	TH1F *h5 = new TH1F("coincident_heights_2", "Coincidences Heights 2 Channels",200,0,200*0.48828125);
	TH1F *h6 = new TH1F("coincident_heights_4", "Coincidences Heights 4 Channels",200,0,200*0.48828125);
	
	for(int event : coincident_events_2){ //2-fold Coincident Heights
		for(std::vector<double> pair : hit_heights){
			if(event == pair[0]){
				std::cout << "Coincidence " << pair[0] << " has " << pair[1] << " pulse height" << std::endl;
				h5->Fill(pair[1]);
			}
		}
	}
	
	for(int event : coincident_events_4){ //4-fold Coincident Heights
		for(std::vector<double> pair : hit_heights){
			if(event == pair[0]){
				std::cout << "Coincidence " << pair[0] << " has " << pair[1] << " pulse height" << std::endl;
				h6->Fill(pair[1]);
			}
		}
	}
    
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
	
	//************************************************************************************************************
	// Fill Max Coincident Charge Histogram
	TH1F *h7 = new TH1F("coincidence_charges_max_2", "Max Coincidence Charges 2 Channels",200,0,30*0.48828125);
	TH1F *h8 = new TH1F("coincidence_charges_max_4", "Max Coincidence Charges 4 Channels",200,0,200*0.48828125);
	
	double max_2 = 0;
	for(int event : coincident_events_2){ //2-fold max coincident charge
		max_2 = 0;
		for(std::vector<double> pair : hit_charges){
			if(event == pair[0] and pair[1] > max_2){
				max_2 = pair[1];
			}
		}
		h7->Fill(max_2);
		std::cout << "Coincidence " << event << " has charge max " << max_2 << std::endl;
	}
	
	double max_4 = 0;
	for(int event : coincident_events_4){ //4-fold max coincident charge
		max_4 = 0;
		for(std::vector<double> pair : hit_charges){
			if(event == pair[0] and pair[1] > max_4){
				max_4 = pair[1];
			}
		}
		h8->Fill(max_4);
		std::cout << "Coincidence " << event << " has charge max " << max_4 << std::endl;
	}
	
	// Print Max Coincident Charge Histogram
	TCanvas *c7 = new TCanvas("C7");
	h7->Draw();
	h7->GetXaxis()->SetTitle("Pulse Charge (photoelectrons)");
	h7->GetYaxis()->SetTitle("Number of Hits");
	h7->Fit("gaus","","",0.5,3);
	gStyle->SetOptFit(11);
	c7->SaveAs("coincidence_charges_MAX_2.png");
	
	TCanvas *c8 = new TCanvas("C8");
	h8->Draw();
	h8->GetXaxis()->SetTitle("Pulse Charge (photoelectrons)");
	h8->GetYaxis()->SetTitle("Number of Hits");
	//h8->Fit("gaus","","",0.5,3);
	gStyle->SetOptFit(11);
	c8->SaveAs("coincidence_charges_MAX_4.png");
	
	//************************************************************************************************************
	// Fill Sum Coincident Charge Histogram
	TH1F *h9 = new TH1F("coincidence_charges_sum_2", "Sum Coincidence Charges 2 Channels",200,0,30*0.48828125);
	TH1F *h10 = new TH1F("coincidence_charges_sum_4", "Sum Coincidence Charges 4 Channels",200,0,200*0.48828125);
	
	double sum_2 = 0;
	for(int event : coincident_events_2){ //2-fold sum coincident charge
		sum_2 = 0;
		for(std::vector<double> pair : hit_charges){
			if(event == pair[0]){
				sum_2 += pair[1];
			}
		}
		h9->Fill(sum_2);
		std::cout << "Coincidence " << event << " has charge sum " << sum_2 << std::endl;
	}
	
	double sum_4 = 0;
	for(int event : coincident_events_4){ //4-fold sum coincident charge
		sum_4 = 0;
		for(std::vector<double> pair : hit_charges){
			if(event == pair[0]){
				sum_4 += pair[1];
			}
		}
		h10->Fill(sum_4);
		std::cout << "Coincidence " << event << " has charge sum " << sum_4 << std::endl;
	}
	
	// Print Sum Coincident Charge Histogram
	TCanvas *c9 = new TCanvas("C9");
	h9->Draw();
	h9->GetXaxis()->SetTitle("Pulse Charge (photoelectrons)");
	h9->GetYaxis()->SetTitle("Number of Hits");
	h9->Fit("gaus","","",1,6);
	gStyle->SetOptFit(11);
	c9->SaveAs("coincidence_charges_SUM_2.png");
	
	TCanvas *c10 = new TCanvas("C10");
	h10->Draw();
	h10->GetXaxis()->SetTitle("Pulse Charge (photoelectrons)");
	h10->GetYaxis()->SetTitle("Number of Hits");
	//h10->Fit("gaus","","",0.5,3);
	gStyle->SetOptFit(11);
	c10->SaveAs("coincidence_charges_SUM_4.png");
	
	//************************************************************************************************************
	// Fill coincidence time difference histogram
	TH1F *h11 = new TH1F("coincidence_time_difference_2", "Coincidence Time Differences 2 (Log Scale)",200,0,16000*0.48828125);
	TH1F *h12 = new TH1F("coincidence_time_difference_4", "Coincidence Time Differences 4 (Log Scale)",200,0,16000*0.48828125);
	TH1F *h13 = new TH1F("coincident_times_2", "Coincident Times 2",200,0,17000*0.48828125);
	TH1F *h14 = new TH1F("coincident_times_4", "Coincident Times 4",200,0,17000*0.48828125);
	
	
	for(int event : coincident_events_2){// 2-fold coincident time difference
		for(std::vector<double> pair : hit_times){
			if(event == pair[0]){
				std::cout << "Coincidence event " << pair[0] << " has a hit at time " << pair[1] << std::endl;
				h13->Fill(pair[1]);
				coincident_times_2.push_back(pair);
			}
		}
	}
	std::sort(coincident_times_2.begin(), coincident_times_2.end());
	while(coincident_times_2.size() >= 2){
		if(coincident_times_2[1][0] == coincident_times_2[0][0]){
			h11->Fill(coincident_times_2[1][1] - coincident_times_2[0][1]);
			std::cout << "Coincidence Number " << coincident_times_2[0][0] << " has time difference "
					  << coincident_times_2[1][1] - coincident_times_2[0][1] << std::endl;
			coincident_times_2.erase(coincident_times_2.begin()+1);
		} else {
			coincident_times_2.erase(coincident_times_2.begin());
		}
	}
	
	for(int event : coincident_events_4){// 4-fold coincident time difference
		for(std::vector<double> pair : hit_times){
			if(event == pair[0]){
				std::cout << "Coincidence event " << pair[0] << " has a hit at time " << pair[1] << std::endl;
				h14->Fill(pair[1]);
				coincident_times_4.push_back(pair);
			}
		}
	}
	std::sort(coincident_times_4.begin(), coincident_times_4.end());
	while(coincident_times_4.size() >= 2){
		if(coincident_times_4[1][0] == coincident_times_4[0][0]){
			h12->Fill(coincident_times_4[1][1] - coincident_times_4[0][1]);
			std::cout << "Coincidence Number " << coincident_times_4[0][0] << " has time difference "
					  << coincident_times_4[1][1] - coincident_times_4[0][1] << std::endl;
			coincident_times_4.erase(coincident_times_4.begin()+1);
		} else {
			coincident_times_4.erase(coincident_times_4.begin());
		}
	}
	
	// Print coincidence time difference
	TCanvas *c11 = new TCanvas("C11");
	h11->Draw();
	h11->GetXaxis()->SetTitle("Pulse Time Differences (ns)");
	h11->GetYaxis()->SetTitle("Number of Hits");
	gStyle->SetOptFit(11);
	gPad->SetLogy();
	c11->SaveAs("coincidence_time_differences_LOG_2.png");
	
	TCanvas *c12 = new TCanvas("C12");
	h12->Draw();
	h12->GetXaxis()->SetTitle("Pulse Time Differences (ns)");
	h12->GetYaxis()->SetTitle("Number of Hits");
	gStyle->SetOptFit(11);
	gPad->SetLogy();
	c12->SaveAs("coincidence_time_differences_LOG_4.png");
	
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
	

	//************************************************************************************************************
	// Fill Coincident Channels Histogram
	TH1F *h15 = new TH1F("coincident_channels_2", "Coincident Channels 2",200,0,20);
	TH1F *h16 = new TH1F("coincident_channels_4", "Coincident Channels 4",200,0,20);
	
	for(int i : coincident_events_2){ //2-fold channel numbers of coincident hits
		for(int j : myChannels){
			tt[j]->GetEvent(i);
			for(int k = 1; k <= wf[j]->numPulses; k++){
				h15->Fill(j);
				std::cout << "Coincident event " << i << " has a pulse in channel " << j << std::endl;
			}
		}
	}
	for(int i : coincident_events_4){ //4-fold channel numbers of coincident hits
		for(int j : myChannels){
			tt[j]->GetEvent(i);
			for(int k = 1; k <= wf[j]->numPulses; k++){
				h16->Fill(j);
				std::cout << "Coincident event " << i << " has a pulse in channel " << j << std::endl;
			}
		}
	}
	
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
	
	
	std::cout << "Live Time: " << liveTime << " seconds" << std::endl;
	
	std::cout << "Number of 2-fold coincident events: " << coincident_events_2.size() << std::endl;
	std::cout << "Rate of 2-fold coincident events: " << coincident_events_2.size() / liveTime << " hits/second" << std::endl;
	
	std::cout << "Number of 4-fold coincident events: " << coincident_events_4.size() << std::endl;
	std::cout << "Rate of 4-fold coincident events: " << coincident_events_4.size() / liveTime << " hits/second" << std::endl;
	
	
	fin->Close();
	return 0;
}
	
