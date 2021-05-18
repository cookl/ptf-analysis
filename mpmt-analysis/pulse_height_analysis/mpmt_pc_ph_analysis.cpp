#include "WaveformFitResult.hpp"
#include "ScanPoint.hpp"
#include "TFile.h"
#include "TCanvas.h"
#include "TH1D.h"
#include "TH2F.h"
#include "THStack.h"
#include "TGraph.h"
#include "TGaxis.h"
#include "TF1.h"
#include "TStyle.h"
#include "TLegend.h"
#include "TPaveStats.h"
#include "TProfile.h"

#include <iostream>
#include <vector>

int main( int argc, char* argv[] ) {

    if ( argc != 2 ){
        std::cerr<<"Usage: ptf_ttree_analysis.app ptf_analysis.root\n";
        exit(0);
    }

    TFile * fin = new TFile( argv[1], "read" );
    
    // Init arrays and variables for scatterplots
    double pc[150000];
    double ph[150000];
    int n=0;
    TH2F *h2 = new TH2F("pc_ph_hist","Histogram of pulse height and pulse charge",100,0,50,100,-30,70);

    // Get the waveform fit TTree
    std::cout << "TTree 0" << std::endl;
    TTree * tt0 = (TTree*)fin->Get("ptfanalysis0");
    WaveformFitResult * wf0 = new WaveformFitResult;
    if(tt0) wf0->SetBranchAddresses( tt0 );

    // For each fixed window
    for(int i = 0; i < tt0->GetEntries()-1; i++){
        tt0->GetEvent(i );
                
        // For each pulse found in each fixed window
        // Plot ph and corresponding pc
        for(int k = 0; k < wf0->numPulses; k++){
            auto pulse_time = wf0->pulseTimes[k];
            if (pulse_time>=2100 && pulse_time<=2400) {
                pc[n]=wf0->qsum*1000.0;
                ph[n]=wf0->pulseCharges[k]*1000.0;
                h2->Fill(wf0->qsum*1000.0,wf0->pulseCharges[k]*1000.0);
                n++;
            }
        }
  }
    
    TCanvas *c1 = new TCanvas("C1");
    TGraph *pc_ph = new TGraph(n,ph,pc);
    pc_ph->SetTitle("Pulse charge vs pulse height");
    pc_ph->GetXaxis()->SetRangeUser(0,50);
    pc_ph->GetYaxis()->SetRangeUser(-30,70);
    pc_ph->GetXaxis()->SetTitle("Pulse height (mV)");
    pc_ph->GetYaxis()->SetTitle("Pulse charge (mV * 8ns)");
    pc_ph->Draw("ap");
    c1->SaveAs("mpmt_pulse_charge_vs_height.png");
    
    TCanvas *c2 = new TCanvas("C2");
    h2->Draw("COLZ");
    c2->SaveAs("mpmt_pc_vs_ph_histogram.png");

    fin->Close();
    return 0;
}



