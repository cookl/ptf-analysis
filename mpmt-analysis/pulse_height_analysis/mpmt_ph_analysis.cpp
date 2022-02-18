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
using namespace std;

const double bin_unit = 0.4883;

double pc[1000000];
double ph[1000000];
double input_pulse[100000];
double pmt_pulse[1000000];

TH1F *pedestal;
TH1F *laser_pc;
TH1F *laser_ph;
TH1F *before_ph;
TH1F *after_ph;
TH1F *total_ph;
TH2F *pc_ph_hist;
TH1F *input_ph_shift;
TH1F *pmt_ph_shift;
TH2F *input_pmt_hist;

TTree * tt1;
TTree * tt2;

WaveformFitResult * wf1;
WaveformFitResult * wf2;

double peak_to_valley;
double ptv_min_amp;
double ptv_max_amp;

// Initialize histograms for output
// Note: bins quantized in units of bin_unit
void initHist(int pc_lower_range){
    pedestal = new TH1F("pedestal","Pedestal value per waveform",60,0.9990,1.0050);
    laser_pc = new TH1F("pc","Laser Pulse Charge",200,-1*pc_lower_range*bin_unit,180*bin_unit);
    laser_ph = new TH1F("ph-laser","Laser Pulse Height",200,0,bin_unit*200);
    before_ph = new TH1F("ph-before","Pulse Height Before Laser",200,0,bin_unit*200);
    after_ph = new TH1F("ph-after", "Pulse Height After Laser",200,0,bin_unit*200);
    total_ph = new TH1F("ph-total", "Pulse Height", 200,0,bin_unit*200);        //123,820*bin_unit,943*bin_unit);
    pc_ph_hist = new TH2F("pc_ph_hist","Histogram of pulse height and pulse charge",100,0,100*bin_unit,120,-20*2*bin_unit,100*2*bin_unit);
    pmt_ph_shift = new TH1F("pmt-pulse-shift", "PMT: Pulse height spread",12,275*8,287*8);
    input_ph_shift = new TH1F("inputed-pulse-shift", "Input signal: Pulse height spread",13,262*8,275*8);
    input_pmt_hist = new TH2F("input_pmt_hist", "Input vs PMT pulse time",13,262*8,275*8,12 ,275*8,287*8);        //175,225*8,400*8,175,225*8,400*8);//,12,275*8,287*8);
}

// Collect pulse heights and plot them on the related histograms
//  - k indicates the pulse_num in the current waveform
//  - pulse_low (ns) is the lower bin for the detected pulse
//  - pulse_high (ns) is the higher bin for the detected pulse
void collectPH(int k, int pulse_low, int pulse_high) {
    auto pulse_time = wf2->pulseTimes[k];
    auto pulse_height = wf2->pulseCharges[k];

    total_ph->Fill(pulse_height*1000.0);

    if (pulse_time <= pulse_low) {
        before_ph->Fill(pulse_height * 1000.0);
        return;
    }
    if (pulse_time >= pulse_high) {
        after_ph->Fill(pulse_height * 1000.0);
        return;
    }

    laser_ph->Fill(pulse_height * 1000.0);
    
    pmt_ph_shift->Fill(pulse_time);
    input_ph_shift->Fill(wf1->pulseTimes[0]);
}

// Calculate peak-to-valley ratio for pulse charge
// Range can change depending on the run settings
void peakToValley(int range_low, int range_high, int pc_lower_range) {
    int bin_low = range_low/bin_unit;
    int bin_high = range_high/bin_unit;
    for (auto pulse_charge=bin_unit*bin_low; pulse_charge<=bin_unit*bin_high; pulse_charge+=bin_unit) {
        auto bin_num = (pc_lower_range*bin_unit + pulse_charge)/bin_unit;
        auto pc_count = laser_pc->GetBinContent(bin_num);
        if (pc_count<ptv_min_amp) {                 // find min amp
            ptv_min_amp=pc_count;
            continue;
        }
        if (pc_count>ptv_max_amp) ptv_max_amp=pc_count; // find max amp
    }
    peak_to_valley = ptv_max_amp/ptv_min_amp;           // find peak-to-valley ratio
}


// Print plots
//  - pick "noise" or "p.e." for pulse charge fitting
//  - pick "scaled" or "unscaled" for pulse height plot settings
//  - event_num is number of pulses
void printPlots(string pc_fit_type, string ph_axis_type, int event_num) {
    // Print pulse charge histogram on a log-scale
    TCanvas *c1 = new TCanvas("C1");
    laser_pc->GetXaxis()->SetTitle("Pulse charge (mV * 8ns)");
    laser_pc->GetYaxis()->SetTitle("Number of events");
    gPad->SetLogy();    // comment this line to view linear-scale histogram
    if (pc_fit_type == "noise") laser_pc->Fit("gaus","Q","C",-4,4);
    if (pc_fit_type == "p.e.") laser_pc->Fit("gaus","Q","C",10,30);
    gStyle->SetOptFit(11);
    laser_pc->SetMarkerStyle(6);
    gStyle->SetOptStat(11);
    laser_pc->Draw("P");
    c1->Update();
    // Print peak-to-valley ratio on histogram
    TPaveStats *ps = (TPaveStats*)c1->GetPrimitive("stats");
    ps->SetName("peak-to-valley");
    TList *listOfLines = ps->GetListOfLines();
    std::string text = "Peak-to-valley   " + std::to_string(peak_to_valley);
    TLatex *myt = new TLatex(0,0, text.c_str());
    listOfLines->Add(myt);
    laser_pc->SetStats(0);
    c1->Modified();
    c1->SaveAs("mpmt_pc.png");
    
    // Print separated pulse height histogram
    TCanvas *c2 = new TCanvas("C2");
    if (ph_axis_type == "scaled") {
        // pad1: ph during laser
        // pad2: ph before laser
        // pad3: ph after laser
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
    }
    if (ph_axis_type == "unscaled") {
        laser_ph->Draw();
        before_ph->SetLineColor(kRed);
        after_ph->SetLineColor(kMagenta);
        after_ph->Draw("][sames");
        before_ph->Draw("][sames");
        after_ph->GetXaxis()->SetTitle("Pulse height (mV)");
        after_ph->GetYaxis()->SetTitle("Number of events");

        TLegend *legend = new TLegend(0.5,0.5,0.9,0.7);
        legend->SetHeader("Legend","C");
        legend->AddEntry(laser_ph,"pulse height from laser","l");
        legend->AddEntry(before_ph,"pulse height before laser","l");
        legend->AddEntry(after_ph,"pulse height after laser");
        legend->Draw();
    }
    c2->SaveAs("mpmt_pulse_height_separated.png");

    // Print total pulse height
    TCanvas *c3 = new TCanvas("C3");
    total_ph->Draw();
    total_ph->GetXaxis()->SetTitle("Pulse height (mV)");
    total_ph->GetYaxis()->SetTitle("Number of events");
    c3->SaveAs("mpmt_pulse_height_total.png");
    
    // Print pulse charge vs pulse height
    TCanvas *c4 = new TCanvas("C4");
    TGraph *pc_ph_scatter = new TGraph(event_num,ph,pc);
    pc_ph_scatter->SetTitle("Pulse charge vs pulse height");
    pc_ph_scatter->GetXaxis()->SetRangeUser(0,40);
    pc_ph_scatter->GetYaxis()->SetRangeUser(-20,100);
    pc_ph_scatter->GetXaxis()->SetTitle("Pulse height (mV)");
    pc_ph_scatter->GetYaxis()->SetTitle("Pulse charge (mV * 8ns)");
    pc_ph_scatter->Draw("ap");
    c4->SaveAs("mpmt_pc_vs_height_scatter.png");
    
    TCanvas *c5 = new TCanvas("C5");
    pc_ph_hist->GetXaxis()->SetTitle("Pulse height (mV)");
    pc_ph_hist->GetYaxis()->SetTitle("Pulse charge (mV*ns)");
    pc_ph_hist->Draw("COLZ");
    c5->SaveAs("mpmt_pc_vs_height_hist.png");
    
    TCanvas *c6 = new TCanvas("C6");
    pmt_ph_shift->GetXaxis()->SetTitle("Time (ns)");
    pmt_ph_shift->GetYaxis()->SetTitle("Number of events");
    pmt_ph_shift->Draw();
    c6->SaveAs("mpmt_pulse_shift.png");
    
    TCanvas *c7 = new TCanvas("C7");
    input_ph_shift->GetXaxis()->SetTitle("Time (ns)");
    input_ph_shift->GetYaxis()->SetTitle("Number of events");
    input_ph_shift->Draw();
    c7->SaveAs("input_pulse_shift.png");
    
    TCanvas *c8 = new TCanvas("C8");
    input_pmt_hist->GetXaxis()->SetTitle("Input signal: time of minimum amplitude (ns)");
    input_pmt_hist->GetYaxis()->SetTitle("PMT signal: time of minimum amplitude (ns)");
    input_pmt_hist->Draw("COLZ");
    c8->SaveAs("input_pmt_pulse_hist.png");
    
    TCanvas *c9 = new TCanvas("C9");
    TGraph *input_pmt_scatter = new TGraph(event_num,input_pulse,pmt_pulse);
    input_pmt_scatter->SetTitle("Input vs PMT pulse time");
    input_pmt_scatter->GetXaxis()->SetRangeUser(2100,2200);  //1800,3200
    input_pmt_scatter->GetYaxis()->SetRangeUser(2000,2600);
    input_pmt_scatter->GetXaxis()->SetTitle("Input signal: time of minimum amplitude (ns)");
    input_pmt_scatter->GetYaxis()->SetTitle("PMT signal: time of minimum amplitude (ns)");
    input_pmt_scatter->Draw("ap");
    c9->SaveAs("input_pmt_pulse_scatter.png");
    
    TCanvas *c10 = new TCanvas("C10");
    pedestal->SetStats(1);
    gStyle->SetOptStat(1111);
    pedestal->GetXaxis()->SetTitle("Pedestal per waveform (V)");
    pedestal->GetYaxis()->SetTitle("Number of events");
    pedestal->Draw();
    pedestal->Fit("gaus");
    gStyle->SetOptFit(11);
    c10->SaveAs("pedestals.png");
}

int main( int argc, char* argv[] ) {

    if ( argc != 2 ){
        std::cerr<<"Usage: ptf_ttree_analysis.app ptf_analysis.root\n";
        exit(0);
    }

    TFile * fin = new TFile( argv[1], "read" );

    // Get the waveform fit TTree for channel 1 and 2
    //  - Let tt1 be the inputted electrical signal (channel 1)
    //  - Let tt2 be the detected PMT signal (channel 2)
    std::cout << "TTree 0" << std::endl;
    tt1 = (TTree*)fin->Get("ptfanalysis1");
    tt2 = (TTree*)fin->Get("ptfanalysis2");
    wf1 = new WaveformFitResult;
    wf2 = new WaveformFitResult;
    if(tt1) wf1->SetBranchAddresses( tt1 );
    if(tt2) wf2->SetBranchAddresses( tt2 );
    
    // Initiliaze histograms
    double pc_x_low = 20;
    initHist(pc_x_low);

    // Reset peak-to-valley calculation variables
    ptv_min_amp = 10000;
    ptv_max_amp = 0;

    // Reset event number for array indexing
    int event_num=0;

    // For each waveform:
    for(int i = 0; i < tt2->GetEntries()-1; i++){
        tt1->GetEvent(i);
        tt2->GetEvent(i);

        // Collect pulse charge
        auto pulse_charge = wf2->qsum;
        laser_pc->Fill(pulse_charge * 1000.0); // Convert to mV
        
        // Collect pedestal
        pedestal->Fill(wf2->qped);
                        
        // For each pulse:
        for(int k = 0; k < wf2->numPulses; k++){
            // Collect pulse height
            int pulse_low = 2230;
            int pulse_high = 2270;
            collectPH(k,pulse_low,pulse_high);

            // Collect data for remaining analyses
            auto pulse_time = wf2->pulseTimes[k];
            if (pulse_time>=pulse_low && pulse_time<=pulse_high) {
                // Compare time of PMT pulse to time of input pulse
                input_pulse[event_num]=wf1->pulseTimes[0];
                pmt_pulse[event_num]=pulse_time;
                input_pmt_hist->Fill(wf1->pulseTimes[0],pulse_time);
                
                // Compare PMT ph to pc
                pc[event_num]=wf2->qsum*1000.0;
                ph[event_num]=wf2->pulseCharges[k]*1000.0;
                pc_ph_hist->Fill(wf2->pulseCharges[k]*1000.0,wf2->qsum*1000.0);
                event_num++;
            }
        }
    }
    
    // Find peak-to-valley ratio
    // Range (mV*8ns) depends on run
    int range_low = 3;
    int range_high = 26;
    peakToValley(range_low, range_high, pc_x_low);
    std::cout << "Peak-to-valley ratio: " << ptv_max_amp << "/" << ptv_min_amp << " = " << peak_to_valley << std::endl;;

    // Print plots
    //  - pick "noise" or "p.e." for pulse charge fitting
    //  - pick "scaled" or "unscaled" for pulse height plot settings
    printPlots("noise", "scaled", event_num);
    
    fin->Close();
    return 0;
}
