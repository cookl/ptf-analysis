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

    // Get the waveform fit TTree
    std::cout << "TTree 0" << std::endl;
    TTree * tt0 = (TTree*)fin->Get("ptfanalysis0");
    WaveformFitResult * wf0 = new WaveformFitResult;
    if(tt0) wf0->SetBranchAddresses( tt0 );

    // Create histogram
    /* NEED
     TO
     CHANGE
     HIGHER
     END
     OF
     RANGE*/
    TH1F *laser_pulse_charge    = new TH1F("pc","Laser Pulse Charge Scaled",224,-10,0.000488*224*1000);

    // Loop the first scan point and print something
    for(int i = 0; i < tt0->GetEntries()-1; i++){
        tt0->GetEvent(i );
        
        if (i<20) {
//            std::cout << "i: " << i << ", charge: " << wf0->qsum << ", baseline: " << wf0->qped<<  std::endl;
        }
        
        laser_pulse_charge->Fill(wf0->qsum * 1000.0); // Convert to mV
  }

    TCanvas *c1 = new TCanvas("C1");
    laser_pulse_charge->GetXaxis()->SetTitle("Pulse charge (mV)");
    laser_pulse_charge->GetYaxis()->SetTitle("Number of events");
    gPad->SetLogy();
    laser_pulse_charge->Draw();
    
    c1->SaveAs("mpmt_pulse_charge.png");
    fin->Close();
    return 0;
}



