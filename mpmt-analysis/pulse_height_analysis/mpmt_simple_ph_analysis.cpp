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
//using namespace std;

int main( int argc, char* argv[] ) {

    if ( argc != 2 ){
        std::cerr<<"Usage: ptf_ttree_analysis.app ptf_analysis.root\n";
        exit(0);
    }

    TH1F *h1 = new TH1F("pulse_height","Pulse Height",200,0,200*0.48828125);

    TFile * fin = new TFile( argv[1], "read" );
    TTree * tt0;    
    WaveformFitResult * wf0;

    // Get the waveform fit TTree for channel 0
    std::cout << "TTree 0" << std::endl;
    tt0 = (TTree*)fin->Get("ptfanalysis0");
    wf0 = new WaveformFitResult;
    if(tt0) wf0->SetBranchAddresses( tt0 );


    std::cout << "Analyzing " << tt0->GetEntries() << " waveforms" << std::endl;
    // Loop over each waveform:
    for(int i = 0; i < tt0->GetEntries()-1; i++){
        tt0->GetEvent(i);

	std::cout << "Number of pulses found: " << wf0->numPulses << std::endl;
        // Loop over each pulse:
        for(int k = 0; k < wf0->numPulses; k++){

	  std::cout << "Pulse " << k << " has pulse height " << wf0->pulseCharges[k]*1000.0 
		    << "mV " << wf0->pulseTimes[k] << "ns" 
		    << std::endl;
	  if(wf0->pulseTimes[k] > 2100 and wf0->pulseTimes[k] < 2400 and wf0->pulseCharges[k]*1000.0 > 4.5)
	    h1->Fill(wf0->pulseCharges[k]*1000.0);
        }
    }
    

    // Print total pulse height
    TCanvas *c3 = new TCanvas("C3");
    h1->Draw();
    h1->GetXaxis()->SetTitle("Pulse height (mV)");
    h1->GetYaxis()->SetTitle("Number of events");
    h1->Fit("gaus","","",8,22);
    gStyle->SetOptFit(11);
    c3->SaveAs("mpmt_pulse_height.png");

    
    fin->Close();
    return 0;
}
