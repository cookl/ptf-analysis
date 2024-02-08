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

    if ( argc != 3  ){
        std::cerr<<"Usage: ptf_ttree_analysis.app ptf_analysis.root run# \n";
        exit(0);
    }

    TH1F *h1[19]; 
    
    // = new TH1F("pulse_height","Pulse Height",200,0,200*0.48828125);

    TFile * fin = new TFile( argv[1], "read" );
    TTree * tt0;    
    WaveformFitResult * wf0;
    
    int run_num = std::atoi(argv[2]);


    for(int ch=0; ch<19; ch++){
        std::cout << "Pulse height channel " << ch << std::endl;
        std::string name = "pulse_height_ch" + std::to_string(ch);
        h1[ch] = new TH1F(name.c_str(),name.c_str(),2000,0,2000);
        // Get the waveform fit TTree for channel 0
        std::cout << "TTree 0" << std::endl;
        std::string ttree_name = "ptfanalysis"+std::to_string(ch);
        tt0 = (TTree*)fin->Get(ttree_name.c_str());
        wf0 = new WaveformFitResult;
        if(tt0) wf0->SetBranchAddresses( tt0 );
        std::cout << "Here"<< std::endl;

        // std::cout << "Analyzing " << tt0->GetEntries() << " waveforms" << std::endl;
        // Loop over each waveform:
        for(int i = 0; i < tt0->GetEntries()-1; i++){
            tt0->GetEvent(i);

        // std::cout << "Number of pulses found: " << wf0->numPulses << std::endl;
            // Loop over each pulse:
            for(int k = 0; k < wf0->numPulses; k++){

                if(wf0->pulseTimes[k] > 0. and wf0->pulseTimes[k] < 400 and wf0->pulseCharges[k]*1000.0 > 15.0 ){

                // if(wf0->pulseTimes[k] > 2100 and wf0->pulseTimes[k] < 2400 and wf0->pulseCharges[k]*1000.0 > 4.5){
                    h1[ch]->Fill(wf0->pulseCharges[k]*1000.0);
                }
            }
        }
        

        // Print total pulse height
        // TCanvas *c3 = new TCanvas("C3");
        // h1->Draw();
        // h1->GetXaxis()->SetTitle("Pulse height (mV)");
        // h1->GetYaxis()->SetTitle("Number of events");
        // h1->Fit("gaus","","",8,22);
        // gStyle->SetOptFit(11);
        // c3->SaveAs("mpmt_pulse_height.png");
    }
    

    char *outputGraphFileName = Form("OutputGraphsRun%d.root",run_num);
	TFile* outputGraphs = new TFile(outputGraphFileName,"RECREATE");
    outputGraphs->cd();
    for(int ch=0; ch<19; ch++){
        h1[ch]->Write();
    }
    outputGraphs->Close();
    fin->Close();

    return 0;
}
