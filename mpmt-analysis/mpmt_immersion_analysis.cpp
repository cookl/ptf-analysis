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
    TTree * tt0;    
    WaveformFitResult * wf0;
	
	
	// Contains the event numbers of all hits in all channels
	std::vector<int> hit_events = {};
	// Contains the event numbers of all hits that happen in 2 or more events
	std::vector<int> coincidence_events = {};
	// Contains the pulse heights of all events in all channels
	std::vector<std::vector<double>> hit_heights = {};
	// Contains the pulse charges of all events in all channels
	std::vector<std::vector<double>> hit_charges = {};
	
	// Contains the time of occurance of all hits in all channels
	std::vector<std::vector<double>> hit_times = {};
	// Contains the time of occurance of all coincidence hits
	std::vector<double> coincidence_times = {};
	
	
	
	int myChannels[16] = {0, 1, 2, 3, 4, 5, 6, 7, 10, 11, 12, 13, 14, 15, 16, 17};
	
	for(int channel : myChannels){
	
		TH1F *h1 = new TH1F("pulse_height", "Pulse Height",200,0,200*0.48828125);
		TH1F *h3 = new TH1F("pulse_charge", "Pulse Charge",200,0,200*0.48828125/7.7);

		// Get the waveform fit TTree for channel 0
		std::cout << "TTree " << channel << std::endl;
		
		switch(channel){
			case 0:
				tt0 = (TTree*)fin->Get("ptfanalysis0");
				break;
			case 1:
				tt0 = (TTree*)fin->Get("ptfanalysis1");
				break;
			case 2:
				tt0 = (TTree*)fin->Get("ptfanalysis2");
				break;
			case 3:
				tt0 = (TTree*)fin->Get("ptfanalysis3");
				break;
			case 4:
				tt0 = (TTree*)fin->Get("ptfanalysis4");
				break;
			case 5:
				tt0 = (TTree*)fin->Get("ptfanalysis5");
				break;
			case 6:
				tt0 = (TTree*)fin->Get("ptfanalysis6");
				break;
			case 7:
				tt0 = (TTree*)fin->Get("ptfanalysis7");
				break;
			case 10:
				tt0 = (TTree*)fin->Get("ptfanalysis10");
				break;
			case 11:
				tt0 = (TTree*)fin->Get("ptfanalysis11");
				break;
			case 12:
				tt0 = (TTree*)fin->Get("ptfanalysis12");
				break;
			case 13:
				tt0 = (TTree*)fin->Get("ptfanalysis13");
				break;
			case 14:
				tt0 = (TTree*)fin->Get("ptfanalysis14");
				break;
			case 15:
				tt0 = (TTree*)fin->Get("ptfanalysis15");
				break;
			case 16:
				tt0 = (TTree*)fin->Get("ptfanalysis16");
				break;
			case 17:
				tt0 = (TTree*)fin->Get("ptfanalysis17");
				break;
		}
		
		wf0 = new WaveformFitResult;
		if(tt0) wf0->SetBranchAddresses( tt0 );


		std::cout << "Analyzing " << tt0->GetEntries() << " waveforms" << std::endl;
		// Loop over each waveform:
		for(int i = 0; i < tt0->GetEntries()-1; i++){
			tt0->GetEvent(i);

			if(wf0->numPulses != 0){
				hit_events.push_back(i);
			}
		
			// Loop over each pulse:
			for(int k = 0; k < wf0->numPulses; k++){
				
				std::cout << "Number of pulses found: " << wf0->numPulses << 
				" Event Number " << i << std::endl;

				std::cout << "Pulse " << k << " has pulse height " << wf0->pulseCharges[k]*1000.0 
				<< "mV " << wf0->pulseTimes[k] << "ns" 
				<< std::endl;
				
				// pulse height
				double pulse_height = wf0->pulseCharges[k]*1000.0;
				h1->Fill(pulse_height);
				
				// pulse charge in photoelectron assume 7.7mV
				double pulse_charge = pulse_height / 7.7;
				h3->Fill(pulse_charge);
				
				// pulse time
				double pulse_time = wf0->pulseTimes[k];
				
				if(wf0->numPulses != 0){
					double ind = i;
					hit_heights.push_back(std::vector<double> {ind, pulse_height});
					hit_charges.push_back(std::vector<double> {ind, pulse_charge});
					hit_times.push_back(std::vector<double> {ind, pulse_time});
					
					
				
				
				}
			}			
		}
		
		// Print total pulse height
		//TCanvas *c1 = new TCanvas("C1");
		//h1->Draw();
		//h1->GetXaxis()->SetTitle("Pulse height (mV)");
		//h1->GetYaxis()->SetTitle("Number of events");
		//h1->Fit("gaus","","",4,22);
		//gStyle->SetOptFit(11);
		
		// Print total pulse charge
		TCanvas *c3 = new TCanvas("C3");
		h3->Draw();
		h3->GetXaxis()->SetTitle("Pulse Charge (photoelectrons)");
		h3->GetYaxis()->SetTitle("Number of Events");
		h3->Fit("gaus","","",4/7.7,22/7.7);
		gStyle->SetOptFit(11);

		switch(channel){
			case 0:
				h1->SetTitle("Channel 0 Pulse Height");
				//c1->SaveAs("mpmt_pulse_height_0.png");
				h3->SetTitle("Channel 0 Pulse Charge");
				c3->SaveAs("mpmt_pulse_charge_0.png");
				break;
			case 1:
				h1->SetTitle("Channel 1 Pulse Height");
				//c1->SaveAs("mpmt_pulse_height_1.png");
				h3->SetTitle("Channel 1 Pulse Charge");
				c3->SaveAs("mpmt_pulse_charge_1.png");
				break;
			case 2:
				h1->SetTitle("Channel 2 Pulse Height");
				//c1->SaveAs("mpmt_pulse_height_2.png");
				h3->SetTitle("Channel 2 Pulse Charge");
				c3->SaveAs("mpmt_pulse_charge_2.png");
				break;
			case 3:
				h1->SetTitle("Channel 3 Pulse Height");
				//c1->SaveAs("mpmt_pulse_height_3.png");
				h3->SetTitle("Channel 3 Pulse Charge");
				c3->SaveAs("mpmt_pulse_charge_3.png");
				break;
			case 4:
				h1->SetTitle("Channel 4 Pulse Height");
				//c1->SaveAs("mpmt_pulse_height_4.png");
				h3->SetTitle("Channel 4 Pulse Charge");
				c3->SaveAs("mpmt_pulse_charge_4.png");
				break;
			case 5:
				h1->SetTitle("Channel 5 Pulse Height");
				//c1->SaveAs("mpmt_pulse_height_5.png");
				h3->SetTitle("Channel 5 Pulse Charge");
				c3->SaveAs("mpmt_pulse_charge_5.png");
				break;
			case 6:
				h1->SetTitle("Channel 6 Pulse Height");
				//c1->SaveAs("mpmt_pulse_height_6.png");
				h3->SetTitle("Channel 6 Pulse Charge");
				c3->SaveAs("mpmt_pulse_charge_6.png");
				break;
			case 7:
				h1->SetTitle("Channel 7 Pulse Height");
				//c1->SaveAs("mpmt_pulse_height_7.png");
				h3->SetTitle("Channel 7 Pulse Charge");
				c3->SaveAs("mpmt_pulse_charge_7.png");
				break;
			case 10:
				h1->SetTitle("Channel 10 Pulse Height");
				//c1->SaveAs("mpmt_pulse_height_10.png");
				h3->SetTitle("Channel 10 Pulse Charge");
				c3->SaveAs("mpmt_pulse_charge_10.png");
				break;
			case 11:
				h1->SetTitle("Channel 11 Pulse Height");
				//c1->SaveAs("mpmt_pulse_height_11.png");
				h3->SetTitle("Channel 11 Pulse Charge");
				c3->SaveAs("mpmt_pulse_charge_11.png");
				break;
			case 12:
				h1->SetTitle("Channel 12 Pulse Height");
				//c1->SaveAs("mpmt_pulse_height_12.png");
				h3->SetTitle("Channel 12 Pulse Charge");
				c3->SaveAs("mpmt_pulse_charge_12.png");
				break;
			case 13:
				h1->SetTitle("Channel 13 Pulse Height");
				//c1->SaveAs("mpmt_pulse_height_13.png");
				h3->SetTitle("Channel 13 Pulse Charge");
				c3->SaveAs("mpmt_pulse_charge_13.png");
				break;
			case 14:
				h1->SetTitle("Channel 14 Pulse Height");
				//c1->SaveAs("mpmt_pulse_height_14.png");
				h3->SetTitle("Channel 14 Pulse Charge");
				c3->SaveAs("mpmt_pulse_charge_14.png");
				break;
			case 15:
				h1->SetTitle("Channel 15 Pulse Height");
				//c1->SaveAs("mpmt_pulse_height_15.png");
				h3->SetTitle("Channel 15 Pulse Charge");
				c3->SaveAs("mpmt_pulse_charge_15.png");
				break;
			case 16:
				h1->SetTitle("Channel 16 Pulse Height");
				//c1->SaveAs("mpmt_pulse_height_16.png");
				h3->SetTitle("Channel 16 Pulse Charge");
				c3->SaveAs("mpmt_pulse_charge_16.png");
				break;
			case 17:
				h1->SetTitle("Channel 17 Pulse Height");
				//c1->SaveAs("mpmt_pulse_height_17.png");
				h3->SetTitle("Channel 17 Pulse Charge");
				c3->SaveAs("mpmt_pulse_charge_17.png");
				break;
		}
	
	}
	TH1F *h2 = new TH1F("hit_coincidences", "Hit Coincidences",200,0,17);
	
	// Fill hit coincidences histogram
	std::sort(hit_events.begin(), hit_events.end());
	int count = 1;
	while(hit_events.size() > 1){
		if(hit_events[0] == hit_events[1]){
			count++;
		} else {
			if(count >= 4){
				h2->Fill(count);
				coincidence_events.push_back(hit_events[0]);
				std::cout << "Event " << hit_events[0] << " had " << count << " coincidences" << std::endl;
			}
			count = 1;
		}
		hit_events.erase(hit_events.begin());
	}

	// Print hit coincidence
	TCanvas *c2 = new TCanvas("C2");
	h2->Draw();
	h2->GetXaxis()->SetTitle("Number of Channels Hit");
	h2->GetYaxis()->SetTitle("Number of Events");
	h2->Fit("gaus","","",4,17);
	gStyle->SetOptFit(11);
	c2->SaveAs("hit_coincidences.png");
	
	//***********************************************************************************************************
	
	// Fill coincidence height histogram
	TH1F *h4 = new TH1F("coincidence_heights", "Coincidences Heights",200,0,200*0.48828125);
	
	for(int entry : coincidence_events){
		for(std::vector<double> pair : hit_heights){
			if(entry == pair[0]){
				std::cout << "Coincidence " << pair[0] << " has " << pair[1] << " pulse height" << std::endl;
				h4->Fill(pair[1]);
			}
		}
	}
    
	// Print coincidence height
	TCanvas *c4 = new TCanvas("C4");
	h4->Draw();
	h4->GetXaxis()->SetTitle("Pulse Height (mV)");
	h4->GetYaxis()->SetTitle("Number of Events");
	h4->Fit("gaus","","",4,50);
	gStyle->SetOptFit(11);
	c4->SaveAs("coincidence_heights.png");
	
	//***********************************************************************************************************
	
	// Fill coincidence charge histogram
	TH1F *h5 = new TH1F("coincidence_charges", "Sum Coincidence Charges",200,0,1600*0.48828125/7.7);
	double sum = 0;
	
	for(int entry : coincidence_events){
		for(std::vector<double> pair : hit_charges){
			if(entry == pair[0]){
				sum += pair[1];
			}
		}
		std::cout << "Coincidence " << entry << " has charge sum " << sum << std::endl;
		h5->Fill(sum);
		sum = 0;		
	}
	
	// Print coincidence charge
	TCanvas *c5 = new TCanvas("C5");
	h5->Draw();
	h5->GetXaxis()->SetTitle("Pulse Charge (photoelectrons)");
	h5->GetYaxis()->SetTitle("Number of Events");
	h5->Fit("gaus","","",4/7.7,20);
	gStyle->SetOptFit(11);
	c5->SaveAs("coincidence_charges_SUM.png");
	
	//***********************************************************************************************************
	
	// Fill coincidence time histogram
	TH1F *h6 = new TH1F("coincidence_times", "Coincidence Time Differences",200,0,20000*0.48828125);
	
	for(int entry : coincidence_events){
		for(std::vector<double> pair : hit_times){
			if(entry == pair[0]){
				std::cout << "Time: " << pair[1] << std::endl;
				coincidence_times.push_back(pair[1]);
			}
		}
	}
	
	std::sort(coincidence_times.begin(), coincidence_times.end());
	for(int time : coincidence_times){
		std::cout << "Time Difference: " << time - coincidence_times[0] << std::endl;
		h6->Fill(time - coincidence_times[0]);
	}
	
	// Print coincidence times
	TCanvas *c6 = new TCanvas("C6");
	h6->Draw();
	h6->GetXaxis()->SetTitle("Pulse Time (ns)");
	h6->GetYaxis()->SetTitle("Number of Events");
	h6->Fit("gaus","","",0,8500);
	gStyle->SetOptFit(11);
	c6->SaveAs("coincidence_times.png");
	
	
	
	
    fin->Close();
    return 0;
}
