#include "WaveformFitResult.hpp"
#include "ScanPoint.hpp"
#include "TFile.h"
#include "TStyle.h"
#include "TCanvas.h"
#include "TH1.h"
#include "TH2.h"
#include "TF1.h"
#include "TPad.h"
#include "TLegend.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TText.h"
#include "TFitResult.h"
#include <iostream>
#include <vector>
#include <map>
#include <fstream>
#include <math.h>
#include <string>
#include <tgmath.h>



// This programs opens the ROOT file
// containing the waveforms' information,
// such as number, time position and 
// height of pulses.
// It outputs a plot of average pulse
// height in PE vs time, for a specified
// interval of time (eg every minute)



void analysis(TTree* tt, WaveformFitResult* wf, const Double_t peHeight, const Double_t peDivision, const Double_t &afpTimeThreshold1, Double_t &afpTimeThreshold2, Double_t dcr, std::vector<Double_t> &vecT, std::vector<Double_t> &vecQ, std::vector<Double_t> &vecA, std::vector<Double_t> &vecDeltaT, std::vector<Double_t> &vecAlsr, std::vector<Double_t> &vecQlsr, std::vector<Double_t> &vecAafp, std::vector<Double_t> &vecQafp, std::map<Int_t, Double_t> &lsrCounts, std::map<Int_t, Double_t> &afpRate, std::map<Int_t, Double_t> &singleAfpRate, std::map<Int_t, Double_t> &multipleAfp_dcr, TH1F* hist) {
    
    
}


int main( int argc, char* argv[] ) {

    // PRE-RUN TEST
    // Only run analysis if all arguments are provided.
    if ( argc != 4 ){
        std::cerr<<"Usage: ptf_ttree_analysis.app ptf_analysis.root run_number ptfanalysis0\n";
    exit(0);
    }


    // Constants
    // Arbritary values of a time window:
    double afpTimeThreshold1 = 2180; //< ns
    double afpTimeThreshold2 = 2300; //< ns

    // P.E. values
    double peHeight   = 13.18; //< mV

    // Gaussian fit bounds
    double gausMin = 3;
    double gausMax = 7;

    // Length of time each average data point analyzes, in minutes
    double periodLen = 1;



    // Open the root file
    TFile * fin = new TFile( argv[1], "read" );
  
    // Get the waveform fit TTree
    TTree * tt = (TTree*)fin->Get(argv[3]);
    WaveformFitResult * wf = new WaveformFitResult;
    wf->SetBranchAddresses( tt );

    // Initialize some variables/histograms
    TH1F * hist = new TH1F("pulse_height","Pulse Height",200,0,200*0.48828125);
    float_t lsrPulseFlt;
    char filename[1024];
    auto canvas = new TCanvas("canvas","",1200,800);



    // Determine General Timestamp info
    // Find the time of the first event
    tt->GetEvent(0);
    int UStaTime = wf->evt_timestamp;

    // Find the time of the final event
    tt->GetEvent((tt->GetEntries())-1);
    int UEndTime = wf->evt_timestamp;

    // Determine the total number of minutes the root file ran for, in minutes.
    // Partial minutes at the end of the run are NOT included (i.e. a 90s run would have totMins = 1)
    int totMins = int(floor((UEndTime - UStaTime)/60));

    // Determine the total number of (full) time periods the root file ran for, in minutes
    int totPeriods = int(floor(totMins / periodLen));

    
    
    //Define the arrays for storing mean pulse height
    double means[totPeriods];
    double meanErrs[totPeriods];
    double mins[totPeriods];
    double minErrs[totPeriods];
    for(int i=0; i<totPeriods; i++){
        mins[i] = periodLen*i;
        minErrs[i] = 0;    // The minute-error is always 0 because we don't want horizontal error bars later
    }


    
    // Create a Histogram for each minute of data
    // Declare some variables
    int endLoop;
    int k = 0; 
    // Note: we set k = 0 BEFORE the looping because we don't want to reset it to 0 once we start counting
    // events until we're done running the program to avoid double counting
    int currentTime;

    // Loop over each time-period
    for(int period = 0; period < totPeriods; period++){
        // Reset a few things before making a new minute-hist
        hist->Reset();
        endLoop = 0;

        // Loop over each waveform in that time-period
        while(endLoop == 0){
            // Get event information from tree.
	    tt->GetEvent(k);

	    // Check if the event is actually in the next time-period
	    if((wf->evt_timestamp - UStaTime) >= ((period + 1)*periodLen*60)){
	      endLoop = 1;
	      k -= 1;    // This is so that this event will be included in the next period
	    } 
	    else {
	        // If there are pulses on waveform
	        if(wf->numPulses > 0){

	            if ((wf->pulseTimes[0]>=afpTimeThreshold1) && (wf->pulseTimes[0]<afpTimeThreshold2)) {
	                // Find the pulse height
		        lsrPulseFlt = wf->pulseCharges[0]*1000.0/peHeight;
		    }

		    // Loop over each pulse
		    for(int i = 0; i < wf->numPulses; i++){
		        // check if it is laser pulse (and not an afterpulse) 
		        if (wf->pulseTimes[i]>=afpTimeThreshold1 && wf->pulseTimes[i]<afpTimeThreshold2){
		            hist->Fill(lsrPulseFlt);
		        }
		    }
	        }
	        
	        k++;
	    }
	}
	
	// Format the histogram
	hist->Draw();
	hist->SetTitle("Laser Intensity");
	hist->GetXaxis()->SetTitle("Pulse height (PE)");
	hist->GetXaxis()->SetRangeUser(0, 20);
	hist->GetYaxis()->SetTitle("Number of events");
	TFitResultPtr fitRes = hist->Fit("gaus","QS","",gausMin,gausMax);
	gStyle->SetOptFit(11);
	// Un-comment the following 2 lines to save each avg-pulse-height histogram
	//sprintf(filename, "../images/%s-%s-pulse-height-hist-minute-%d.png", argv[2], argv[3], min);
	//canvas->SaveAs(filename);

	// Store the fit results
	Double_t mean = fitRes->Parameter(1);
	means[period] = mean;
	Double_t meanErr = fitRes->ParError(1);
	meanErrs[period] = meanErr;
	//std::cout << "Mean: "  << mean << std::endl;
    }



    // Plot the mean pulse height over time
    // Create the graph and canvas
    auto canvas2 = new TCanvas("canvas2","",1200,800);
    TGraphErrors * timeGraph = new TGraphErrors(totPeriods,mins,means,minErrs,meanErrs);
    
    //formatting and saving
    timeGraph->SetTitle("Laser Stability over Time");
    timeGraph->GetXaxis()->SetTitle("Time (min)");
    timeGraph->GetYaxis()->SetTitle("Mean Pulse Height (PE)");
    timeGraph->Draw("AP*");
    canvas2->Update();
    sprintf(filename, "../images/%s-%s-laser-pulse-stability.png", argv[2], argv[3]);
    canvas2->SaveAs(filename);

        



    return 0;
}

