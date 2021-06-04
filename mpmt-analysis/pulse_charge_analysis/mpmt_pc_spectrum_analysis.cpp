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

TTree * tt2;
WaveformFitResult * wf2;


double peak_to_valley;
double ptv_min_amp;
double ptv_max_amp;

// Calculate peak-to-valley ratio for pulse charge
// Range can change depending on the run settings
void peakToValley(int range_low, int range_high, int pc_lower_range, TH1F* pc) {
    int bin_low = range_low/bin_unit;
    int bin_high = range_high/bin_unit;
    for (auto pulse_charge=bin_unit*bin_low; pulse_charge<=bin_unit*bin_high; pulse_charge+=bin_unit) {
        auto bin_num = (pc_lower_range*bin_unit + pulse_charge)/bin_unit;
        auto pc_count = pc->GetBinContent(bin_num);
        if (pc_count<ptv_min_amp) {                 // find min amp
            ptv_min_amp=pc_count;
            continue;
        }
        if (pc_count>ptv_max_amp) ptv_max_amp=pc_count; // find max amp
    }
    peak_to_valley = ptv_max_amp/ptv_min_amp;           // find peak-to-valley ratio
}

int main( int argc, char* argv[] ) {
    
    // Set up files
    TFile * files[5];
    files[0] = new TFile( "../../mpmt_Analysis_run0844.root" , "read" );
    files[1] = new TFile( "../../mpmt_Analysis_run0853.root" , "read" );
    files[2] = new TFile( "../../mpmt_Analysis_run0854.root" , "read" );
    files[3] = new TFile( "../../mpmt_Analysis_run0855.root" , "read" );
    files[4] = new TFile( "../../mpmt_Analysis_run0856.root" , "read" );
    
    // Set up voltages
    int voltages[5] = {1258,1275,1234,1307,1331};
    
    // Set up canvas
    TCanvas *c1 = new TCanvas("C1");
    int color[5] = {1,800,632,616,600}; //black, orange, red, magenta, blue
    
    // For each file:
    for (int v=0; v<4; v++) {
        
        // Get the waveform fit TTree
        tt2 = (TTree*)files[v]->Get("ptfanalysis2");
        wf2 = new WaveformFitResult;
        if(tt2) wf2->SetBranchAddresses( tt2 );
        
        // Initiliaze histograms
        double x_low = 20;
        string hist_name = to_string(voltages[v]) + "V";
        TH1F *pc = new TH1F(hist_name.c_str(),"Laser Pulse Charge",x_low+180,-1*x_low*bin_unit,180*bin_unit);
        
        // Reset peak-to-valley calculation variables
        ptv_min_amp = 10000;
        ptv_max_amp = 0;
        
        // For each waveform:
        for(int i = 0; i < tt2->GetEntries()-1; i++){
            tt2->GetEvent(i);

            // Collect pulse charge
            auto pulse_charge = wf2->qsum;
            pc->Fill(pulse_charge * 1000.0); // Convert to mV
        }
        
        // Find peak-to-valley ratio
        // Range (mV*8ns) depends on run
        int range_low = 3;
        int range_high = 26;
        peakToValley(range_low, range_high, x_low,pc);
        std::cout << "HV: "<< voltages[v] <<"V. Peak-to-valley ratio: " << ptv_max_amp << "/" << ptv_min_amp << " = " << peak_to_valley << std::endl;;
        
        // Draw histogram
        pc->SetLineColor(color[v]);
        gPad->SetLogy();
        gStyle->SetOptStat(11);
        pc->SetMarkerStyle(6);
        if(v==0) {
            pc->GetXaxis()->SetTitle("Pulse charge (mV * 8ns)");
            pc->GetYaxis()->SetTitle("Number of events");
        }
        pc->Draw("][sames");
        c1->Update();
        
        // Print peak-to-valley ratio on histogram
        TPaveStats *ps = (TPaveStats*)c1->GetPrimitive("stats");
        ps->SetName("peak-to-valley");
        TList *listOfLines = ps->GetListOfLines();
        std::string text = "Peak-to-valley   " + std::to_string(peak_to_valley);
        TLatex *myt = new TLatex(0,0, text.c_str());
        listOfLines->Add(myt);
        pc->SetStats(0);
        c1->Modified();
        
        files[v]->Close();
    }
    
    c1->SaveAs("pc_spectrum.png");
    return 0;
}
