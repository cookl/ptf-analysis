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
#include <cstring>

double pc[1000000];
double ph[1000000];
double input_pulse[100000];
double pmt_pulse[1000000];

int main( int argc, char* argv[] ) {

    if ( argc != 2 ){
        std::cerr<<"Usage: ptf_ttree_analysis.app ptf_analysis.root\n";
        exit(0);
    }

    TFile * fin = new TFile( argv[1], "read" );

    // Get the waveform fit TTree for channel 1 and 2
    std::cout << "TTree 0" << std::endl;
    TTree * tt1 = (TTree*)fin->Get("ptfanalysis1");
    TTree * tt2 = (TTree*)fin->Get("ptfanalysis2");
    WaveformFitResult * wf1 = new WaveformFitResult;
    WaveformFitResult * wf2 = new WaveformFitResult;
    if(tt1) wf1->SetBranchAddresses( tt1 );
    if(tt2) wf2->SetBranchAddresses( tt2 );
    
    // Histograms for PMT signals
    // Note: bins quantized in units of 0.4883
    double x_low = 20;
    TH1F *laser_pc    = new TH1F("pc","Laser Pulse Charge",200,-1*x_low*0.4883,180*0.4883);  //200 bins total
    TH1F *laser_ph     = new TH1F("ph-laser","Laser Pulse Height",200,0,0.4883*200);
    TH1F *before_ph = new TH1F("ph-before","Pulse Height Before Laser",200,0,0.4883*200);
    TH1F *after_ph = new TH1F("ph-after", "Pulse Height After Laser",200,0,0.4883*200);
    TH1F *total_ph = new TH1F("ph-total", "Pulse Height", 200,0,0.4883*200);//123,820*0.4883,943*0.4883);//
    TH2F *pc_ph_hist = new TH2F("pc_ph_hist","Histogram of pulse height and pulse charge",100,0,100*0.4883,120,-20*0.4883*2,100*0.4883*2);
    TH1F *pmt_shift = new TH1F("pmt-pulse-shift", "PMT: Pulse height spread",12,275*8,287*8);
    
    // Histograms for input signal
    TH1F *input_shift = new TH1F("inputed-pulse-shift", "Input signal: Pulse height spread",13,262*8,275*8);
    TH2F *input_pmt_hist = new TH2F("input_pmt_hist", "Input vs PMT pulse time",13,262*8,275*8,12 ,275*8,287*8);//175,225*8,400*8,175,225*8,400*8);//,12,275*8,287*8);
    
    //test
//    TH1F *h1 = new TH1F("test", "Test", 200,0,)
    
    int n=0;

    // peak-to-valley calculation variables
    double min_amp = 10000;
    double max_amp = 0;
    double peak_to_valley;
    
    // Fill histograms
    // For each waveform:
    std::cout<< "ch1 num entries: " <<tt1->GetEntries()<< std::endl;
    std::cout<< "ch2 num entries: " <<tt2->GetEntries()<< std::endl;
    for(int i = 0; i < tt2->GetEntries()-1; i++){
        tt1->GetEvent(i);
        tt2->GetEvent(i );

//        std::cout << "ch1 num pulses: " <<wf1->numPulses << std::endl;
//        std::cout << "ch2 num pulses: " <<wf2->numPulses << std::endl;
        
        // pulse charge histogram
        auto pulse_charge = wf2->qsum;
        laser_pc->Fill(pulse_charge * 1000.0); // Convert to mV
                        
        // pulse heights before, during, and after laser
        // For each pulse:
        for(int k = 0; k < wf2->numPulses; k++){
            
//            std::cout << "numPulses: " << wf2->numPulses << "ph: " << wf2->pulseCharges[0] << std::endl;
            
            auto pmt_pulse_time = wf2->pulseTimes[k];
            total_ph->Fill(wf2->pulseCharges[k]*1000.0);
            if (pmt_pulse_time <= 2230) {  //2100
                before_ph->Fill(wf2->pulseCharges[k] * 1000.0);
                continue;
            }
            if (pmt_pulse_time >= 2270) { //3000
                after_ph->Fill(wf2->pulseCharges[k] * 1000.0);
                continue;
            }
            
            laser_ph->Fill(wf2->pulseCharges[k] * 1000.0);
            
            // pulse height spread/shift
            pmt_shift->Fill(pmt_pulse_time);
            input_shift->Fill(wf1->pulseTimes[0]);
            
            // plot time of pulse of input signal vs pmt signal
            
//            std::cout << "ch1 time: " << wf1->pulseTimes[0] << std::endl;
//            std::cout << "ch2 time: " << wf2->pulseTimes[k] << std::endl;
            
            input_pulse[n]=wf1->pulseTimes[0];
            pmt_pulse[n]=pmt_pulse_time;
            input_pmt_hist->Fill(wf1->pulseTimes[0],pmt_pulse_time);
            
            // plot ph vs pc
            pc[n]=wf2->qsum*1000.0;
            ph[n]=wf2->pulseCharges[k]*1000.0;
            pc_ph_hist->Fill(wf2->pulseCharges[k]*1000.0,wf2->qsum*1000.0);
            n++;
            
        }
  }
    
    // find peak-to-valley ratio
    // range depends on run : O
    for (auto pulse_charge=0.4883*8; pulse_charge<=50*0.4883; pulse_charge+=0.4883) {    // higher pulse  range: 6.8362-25.3916
        auto bin_num = (x_low*0.4883 + pulse_charge)/0.4883;
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
//    laser_pc->Fit("gaus","Q","C",-4,4);     // noise fit
    laser_pc->Fit("gaus","Q","C",10,30);      // p.e. fit
//    laser_pc->Fit("gaus","Q","C",10,40);
    gStyle->SetOptFit(11);
    laser_pc->SetMarkerStyle(6);
    gStyle->SetOptStat(11);
    laser_pc->Draw("P");
    c1->Update();
    
    TPaveStats *ps = (TPaveStats*)c1->GetPrimitive("stats");
    ps->SetName("peak-to-valley");
    TList *listOfLines = ps->GetListOfLines();

    std::string text = "Peak-to-valley   " + std::to_string(peak_to_valley);
    TLatex *myt = new TLatex(0,0, text.c_str());
    listOfLines->Add(myt);
    laser_pc->SetStats(0);
    c1->Modified();
    
    c1->SaveAs("mpmt_pulse_charge.png");
    
    // Print pulse height on a scaled axis
    // pad1: ph during laser
    // pad2: ph before laser
    // pad3: ph after laser
    TCanvas *c2 = new TCanvas("C2");
    gStyle->SetOptStat();
    gStyle->SetOptFit(0);
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
    TGraph *pc_ph_scatter = new TGraph(n,ph,pc);
    pc_ph_scatter->SetTitle("Pulse charge vs pulse height");
    pc_ph_scatter->GetXaxis()->SetRangeUser(0,40);
    pc_ph_scatter->GetYaxis()->SetRangeUser(-20,100);
    pc_ph_scatter->GetXaxis()->SetTitle("Pulse height (mV)");
    pc_ph_scatter->GetYaxis()->SetTitle("Pulse charge (mV * 8ns)");
    pc_ph_scatter->Draw("ap");
    c4->SaveAs("mpmt_pulse_charge_vs_height_scatter.png");
    
    TCanvas *c5 = new TCanvas("C5");
    pc_ph_hist->GetXaxis()->SetTitle("Pulse height (mV)");
    pc_ph_hist->GetYaxis()->SetTitle("Pulse charge (mV*ns)");
    pc_ph_hist->Draw("COLZ");
    c5->SaveAs("mpmt_pulse_charge_vs_height_hist.png");
    
    TCanvas *c6 = new TCanvas("C6");
    pmt_shift->GetXaxis()->SetTitle("Time (ns)");
    pmt_shift->GetYaxis()->SetTitle("Number of events");
    pmt_shift->Draw();
    c6->SaveAs("mpmt_pulse_shift.png");
    
    TCanvas *c7 = new TCanvas("C7");
    input_shift->GetXaxis()->SetTitle("Time (ns)");
    input_shift->GetYaxis()->SetTitle("Number of events");
    input_shift->Draw();
    c7->SaveAs("input_pulse_shift.png");
    
    TCanvas *c8 = new TCanvas("C8");
    input_pmt_hist->GetXaxis()->SetTitle("Input signal: time of minimum amplitude (ns)");
    input_pmt_hist->GetYaxis()->SetTitle("PMT signal: time of minimum amplitude (ns)");
    input_pmt_hist->Draw("COLZ");
//    input_pmt_hist->Print("all");
    c8->SaveAs("input_pmt_pulse_hist.png");
    
    TCanvas *c9 = new TCanvas("C9");
    TGraph *input_pmt_scatter = new TGraph(n,input_pulse,pmt_pulse);
    input_pmt_scatter->SetTitle("Input vs PMT pulse time");
    input_pmt_scatter->GetXaxis()->SetRangeUser(2100,2200);  //1800,3200
    input_pmt_scatter->GetYaxis()->SetRangeUser(2000,2600);
    input_pmt_scatter->GetXaxis()->SetTitle("Input signal: time of minimum amplitude (ns)");
    input_pmt_scatter->GetYaxis()->SetTitle("PMT signal: time of minimum amplitude (ns)");
    input_pmt_scatter->Draw("ap");
    input_pmt_scatter->Fit("pol1");
    gStyle->SetOptFit(11);
    c9->SaveAs("input_pmt_pulse_scatter.png");

    fin->Close();
    return 0;
}



