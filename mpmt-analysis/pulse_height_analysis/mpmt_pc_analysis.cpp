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

    // Create histograms
    TH1F *laser_pulse_height     = new TH1F("pc","Laser Pulse Charge",200,0,0.000488*200*1000);

    // Loop the first scan point and print something
    for(int i = 0; i < tt0->GetEntries()-1; i++){
        tt0->GetEvent(i );

        // Loop over the pulses we found
        for(int k = 0; k < wf0->numPulses; k++){

//            auto pulse_time = wf0->pulseTimes[k];
//
//            if (pulse_time <= 2151) {
//                before_laser_height->Fill(wf0->pulseCharges[k] * 1000.0);
//                continue;
//            }
//
//            if (pulse_time >= 2254) {
//                after_laser_height->Fill(wf0->pulseCharges[k] * 1000.0);
//                continue;
//            }
//
//            laser_pulse_height->Fill(wf0->pulseCharges[k] * 1000.0);  // Convert to mV
            
            std::cout << wf0->qsum << std::endl;
        }
  }

    TCanvas *c1 = new TCanvas("C1");
    
//  scaled histograms

    TPad *pad1 = new TPad("pad1", "", 0,0,1,1);
    TPad *pad2 = new TPad("pad2", "",0,0,1,1);
    TPad *pad3 = new TPad("pad3", "",0,0,1,1);
    pad2->SetFillStyle(4000);
    pad3->SetFillStyle(4000);

    pad1->Draw();
    pad1->cd();
    laser_pulse_height->Draw();
    pad1->Update();
    TPaveStats *ps1 = (TPaveStats*)laser_pulse_height->GetListOfFunctions()->FindObject("stats");
    ps1->SetX1NDC(0.2);ps1->SetX2NDC(0.4);
    ps1->SetTextColor(kBlue);
    pad1->Modified();
    c1->cd();

    double ymin = 0;
    double ymax = 80;
    double dy = (ymax-ymin)/0.8;
    double xmin = 0;
    double xmax = 100;
    double dx = (xmax-xmin)/0.8;

    pad2->Range(xmin-0.1*dx,ymin-0.1*dy,xmax+0.1*dx,ymax+0.1*dy);
    pad2->Draw();
    pad2->cd();
    before_laser_height->SetLineColor(kRed);
    before_laser_height->Draw("][sames");
    pad2->Update();
    TPaveStats *ps2 = (TPaveStats*)before_laser_height->GetListOfFunctions()->FindObject("stats");
    ps2->SetX1NDC(0.45); ps2->SetX2NDC(0.65);
    ps2->SetTextColor(kRed);
    pad2->Modified();

    ymax = 600;
    dy = (ymax-ymin)/0.8;
    pad3->Range(xmin-0.1*dx,ymin-0.1*dy,xmax+0.1*dx,ymax+0.1*dy);
    pad3->Draw();
    pad3->cd();
    after_laser_height->SetLineColor(kMagenta);
    after_laser_height->Draw("][sames");
    pad3->Update();
    TPaveStats *ps3 = (TPaveStats*)after_laser_height->GetListOfFunctions()->FindObject("stats");
    ps3->SetX1NDC(0.7); ps3->SetX2NDC(0.9);
    ps3->SetTextColor(kMagenta);
    pad3->Modified();
    
//    unscaled histograms
//        before_laser_height->SetLineColor(kRed);
//        after_laser_height->SetLineColor(kMagenta);
//        after_laser_height->Draw();
//        before_laser_height->Draw("][sames");
//        laser_pulse_height->Draw("][sames");
//        TLegend *legend = new TLegend(0.5,0.5,0.9,0.7);
//        legend->SetHeader("Legend","C");
//        legend->AddEntry(laser_pulse_height,"pulse height from laser","l");
//        legend->AddEntry(before_laser_height,"pulse height before laser","l");
//        legend->AddEntry(after_laser_height,"pulse height after laser");
//        legend->Draw();
    
    c1->SaveAs("mpmt_pulse_height.png");
    fin->Close();
    return 0;
}



