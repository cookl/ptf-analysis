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
    // Note: bins quantized in units of 0.4883
    TH1F *laser_pc    = new TH1F("pc","Laser Pulse Charge",400,-200*0.4883*2,200*0.4883);  //200 bins total
    TH1F *laser_ph     = new TH1F("ph-laser","Laser Pulse Height",200,0,0.4883*200);
    TH1F *before_ph = new TH1F("ph-before","Pulse Height Before Laser",200,0,0.4883*200);
    TH1F *after_ph = new TH1F("ph-after", "Pulse Height After Laser",200,0,0.4883*200);
    TH1F *total_ph = new TH1F("ph-total", "Pulse Height", 200,0,0.4883*200);
    
    // Init arrays and variables for scatterplots
    double pc[150000];
    double ph[150000];
    int n=0;
    TH2F *h2 = new TH2F("pc_ph_hist","Histogram of pulse height and pulse charge",100,0,100*0.4883,120,-40*0.4883*2,80*0.4883*2);

    // peak-to-valley calculation variables
    double min_amp = 10000;
    double max_amp = 0;
    double peak_to_valley;
    
    // Fill histograms
    for(int i = 0; i < tt0->GetEntries()-1; i++){
        tt0->GetEvent(i );
                
        // pulse charges
        auto pulse_charge = wf0->qsum;
        laser_pc->Fill(pulse_charge * 1000.0); // Convert to mV
                
        // pulse heights before, during, and after laser
        for(int k = 0; k < wf0->numPulses; k++){
            auto pulse_time = wf0->pulseTimes[k];
            total_ph->Fill(wf0->pulseCharges[k]*1000.0);
            if (pulse_time <= 2180) {
                before_ph->Fill(wf0->pulseCharges[k] * 1000.0);
                continue;
            }
            if (pulse_time >= 2320) {
                after_ph->Fill(wf0->pulseCharges[k] * 1000.0);
                continue;
            }
            laser_ph->Fill(wf0->pulseCharges[k] * 1000.0);
            
            // plot ph vs pc
            pc[n]=wf0->qsum*1000.0;
            ph[n]=wf0->pulseCharges[k]*1000.0;
            h2->Fill(wf0->pulseCharges[k]*1000.0,wf0->qsum*1000.0);
            n++;
        }
  }
    
    // find peak-to-valley ratio
    // range depends on run : O
    for (auto pulse_charge=1.9532; pulse_charge<=22.4618; pulse_charge+=0.4883) {    // higher pulse  range: 2.4415-49.8066
        auto bin_num = (20*0.4883 + pulse_charge)/0.4883;
        auto pc_count = laser_pc->GetBinContent(bin_num);
        if (pc_count<min_amp) {                 // find min amp
            min_amp=pc_count;
            continue;
        }
        if (pc_count>max_amp) max_amp=pc_count; // find max amp
    }
    peak_to_valley = max_amp/min_amp;           // find peak-to-valley ratio

    std::cout << "Peak-to-valley ratio: " << max_amp << "/" << min_amp << " = " << peak_to_valley << std::endl;;

    // Print pulse charge histogram on a log-scale
    TCanvas *c1 = new TCanvas("C1");
    laser_pc->GetXaxis()->SetTitle("Pulse charge (mV * 8ns)");
    laser_pc->GetYaxis()->SetTitle("Number of events");
    gPad->SetLogy();    // comment this line to view linear-scale histogram
    laser_pc->Draw();
    c1->SaveAs("mpmt_pulse_charge.png");
    
    // Print pulse height on a scaled axis
    // pad1: ph during laser
    // pad2: ph before laser
    // pad3: ph after laser
    TCanvas *c2 = new TCanvas("C2");
    TPad *pad1 = new TPad("pad1", "", 0,0,1,1);
    TPad *pad2 = new TPad("pad2", "",0,0,1,1);
    TPad *pad3 = new TPad("pad3", "",0,0,1,1);
    pad2->SetFillStyle(4000);
    pad3->SetFillStyle(4000);

    pad1->Draw();
    pad1->cd();
    laser_ph->GetXaxis()->SetTitle("Pulse height (mV)");
    laser_ph->GetYaxis()->SetTitle("Number of events");
    laser_ph->Draw();
    pad1->Update();
    TPaveStats *ps1 = (TPaveStats*)laser_ph->GetListOfFunctions()->FindObject("stats");
    ps1->SetX1NDC(0.7);ps1->SetX2NDC(0.9);
    ps1->SetY1NDC(0.6);ps1->SetY2NDC(0.75);
    ps1->SetTextColor(kBlue);
    pad1->Modified();
    c2->cd();

    double ymin = 0;
    double ymax = before_ph->GetBinContent(before_ph->GetMaximumBin())+10;
    double dy = (ymax-ymin)/0.8;
    double xmin = 0;
    double xmax = 100;
    double dx = (xmax-xmin)/0.8;
    pad2->Range(xmin-0.1*dx,ymin-0.1*dy,xmax+0.1*dx,ymax+0.1*dy);
    pad2->Draw();
    pad2->cd();
    before_ph->SetLineColor(kRed);
    before_ph->Draw("][sames");
    pad2->Update();
    TPaveStats *ps2 = (TPaveStats*)before_ph->GetListOfFunctions()->FindObject("stats");
    ps2->SetX1NDC(0.45); ps2->SetX2NDC(0.65);
    ps2->SetTextColor(kRed);
    pad2->Modified();

    ymax = after_ph->GetBinContent(after_ph->GetMaximumBin())+10;
    dy = (ymax-ymin)/0.8;
    pad3->Range(xmin-0.1*dx,ymin-0.1*dy,xmax+0.1*dx,ymax+0.1*dy);
    pad3->Draw();
    pad3->cd();
    after_ph->SetLineColor(kMagenta);
    after_ph->Draw("][sames");
    pad3->Update();
    TPaveStats *ps3 = (TPaveStats*)after_ph->GetListOfFunctions()->FindObject("stats");
    ps3->SetX1NDC(0.7); ps3->SetX2NDC(0.9);
    ps3->SetTextColor(kMagenta);
    pad3->Modified();

////    Print pulse heights on an unscaled axis
//    laser_ph->Draw();
//    before_ph->SetLineColor(kRed);
//    after_ph->SetLineColor(kMagenta);
//    after_ph->Draw("][sames");
//    before_ph->Draw("][sames");
//    after_ph->GetXaxis()->SetTitle("Pulse height (mV)");
//    after_ph->GetYaxis()->SetTitle("Number of events");
//
//    TLegend *legend = new TLegend(0.5,0.5,0.9,0.7);
//    legend->SetHeader("Legend","C");
//    legend->AddEntry(laser_ph,"pulse height from laser","l");
//    legend->AddEntry(before_ph,"pulse height before laser","l");
//    legend->AddEntry(after_ph,"pulse height after laser");
//    legend->Draw();

    c2->SaveAs("mpmt_pulse_height_separated.png");
    
    // Print total pulse height
    TCanvas *c3 = new TCanvas("C3");
    total_ph->Draw();
    total_ph->GetXaxis()->SetTitle("Pulse height (mV)");
    total_ph->GetYaxis()->SetTitle("Number of events");
    c3->SaveAs("mpmt_pulse_height_total.png");
    
    // Print pulse charge vs pulse height
    TCanvas *c4 = new TCanvas("C4");
    TGraph *pc_ph = new TGraph(n,ph,pc);
    pc_ph->SetTitle("Pulse charge vs pulse height");
    pc_ph->GetXaxis()->SetRangeUser(0,40);
    pc_ph->GetYaxis()->SetRangeUser(-40,80);
    pc_ph->GetXaxis()->SetTitle("Pulse height (mV)");
    pc_ph->GetYaxis()->SetTitle("Pulse charge (mV * 8ns)");
    pc_ph->Draw("ap");
    c4->SaveAs("mpmt_pulse_charge_vs_height_scatter.png");
    
    TCanvas *c5 = new TCanvas("C5");
    h2->GetXaxis()->SetTitle("Pulse height (mV)");
    h2->GetYaxis()->SetTitle("Pulse charge (mV*ns)");
    h2->Draw("COLZ");
    c5->SaveAs("mpmt_pulse_charge_vs_height_hist.png");
    

    fin->Close();
    return 0;
}



