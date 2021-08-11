#include "WaveformFitResult.hpp"
#include "ScanPoint.hpp"
#include "TCanvas.h"
#include "TFile.h"
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

#include <iostream>
#include <vector>
#include <string>
#include <cstring>
#include <math.h> 
using namespace std;

const double bin_unit = 0.4883;

TTree * tt2;
WaveformFitResult * wf2;

double peak_to_valley;
double ptv_min_amp;
double ptv_max_amp;

// Calculate peak-to-valley ratio for pulse charge
//  - calculated by finding minimum and maximum in given range
void peakToValley(int range_low, int range_high, int pc_lower_range, TH1F* pc) {
    int bin_low = range_low/bin_unit;
    int bin_high = range_high/bin_unit;
    for (auto pulse_charge=bin_unit*bin_low; pulse_charge<=bin_unit*bin_high; pulse_charge+=bin_unit) {
        auto bin_num = (pc_lower_range*bin_unit + pulse_charge)/bin_unit;
        auto pc_count = pc->GetBinContent(bin_num);
        if (pc_count<ptv_min_amp) {                     // find min amp
            ptv_min_amp=pc_count;
            continue;
        }
        if (pc_count>ptv_max_amp) ptv_max_amp=pc_count; // find max amp
    }
    peak_to_valley = ptv_max_amp/ptv_min_amp;           // find peak-to-valley ratio
}

// Define fit function for pulse charge distribution
Double_t fitf(Double_t *x, Double_t *p) {
    Double_t    Sped=0,
                Snoise=0,
                S1=0,
                Sn=0;

    Sped=(1-p[2])/(sqrt(2*M_PI)*p[1])*exp(-0.5*pow((x[0]-p[0])/(p[1]),2)-p[4]);
    if (x[0] > p[0]) Snoise=p[3]*p[2]*exp(-1*p[3]*(x[0]-p[0])-p[4]);
    S1=p[4]/(sqrt(2*M_PI)*p[5])*exp(-0.5*pow((x[0]-p[6]-p[0]-p[2]/p[3])/p[5],2)-p[4]);

    for (int n=2; n<=4; n++) {
        int fact = 1;
        for (int i=n; i>1; i--) {fact *= i;}
        Sn += pow(p[4],n)/(sqrt(2*M_PI*n)*p[5]*fact)*exp(pow((x[0]-n*p[6]-p[0]-p[2]/p[3])/p[5],2)/(-2*n)-p[4]);
    }
    
    return p[7]*(Sped+Snoise+S1+Sn);
}

Double_t Sped(Double_t *x, Double_t *p) {
    Double_t Sped = p[4]*(1-p[0])/(sqrt(2*M_PI)*p[1])*exp(-0.5*pow((x[0]-p[2])/(p[1]),2)-p[3]);
    return Sped;
}

Double_t Snoise(Double_t *x, Double_t *p) {
    Double_t Snoise = 0;
    if (x[0] > p[2]) Snoise=p[4]*p[1]*p[0]*exp(-1*p[1]*(x[0]-p[2])-p[3]);
    return Snoise;
}

Double_t S1(Double_t *x, Double_t *p) {
    Double_t S1 = p[3]*p[2]/(sqrt(2*M_PI)*p[0])*exp(-0.5*pow((x[0]-p[1]-p[4]-p[5])/p[0],2)-p[2]);

    return S1;
}

// "sig1","Q1","miu","N","n","Q0","Qsh"
// 0     ,1   ,2    ,3  ,4  ,5   ,6
Double_t Sn(Double_t *x, Double_t *p) {
    Double_t Sn = 0;
    // for (int n=3; n<=4; n++) {
    int fact = 1;
    for (int i=p[4]; i>1; i--) {
        fact *= i;
    }
    Sn += p[3]*pow(p[2],p[4])/(sqrt(2*M_PI*p[4])*p[0]*fact)*exp(pow((x[0]-p[4]*p[1]-p[5]-p[6])/p[0],2)/(-2*p[4])-p[2]);
    return Sn;
}


int main( int argc, char* argv[] ) {
    
    // Set up files
    TFile * files[8];
    files[0] = new TFile( "../../mpmt_Analysis_run0933.root" , "read" );    
    files[1] = new TFile( "../../mpmt_Analysis_run0917.root" , "read" );
    files[2] = new TFile( "../../mpmt_Analysis_run0914.root" , "read" );
    files[3] = new TFile( "../../mpmt_Analysis_run0912.root" , "read" );
    files[4] = new TFile( "../../mpmt_Analysis_run0909.root" , "read" );
    files[5] = new TFile( "../../mpmt_Analysis_run0910.root" , "read" );
    files[6] = new TFile( "../../mpmt_Analysis_run0911.root" , "read" );
    files[7] = new TFile( "../../mpmt_Analysis_run0916.root" , "read" );
    
    // Set up voltages
    // double voltages[7] = {1209,1234,1258,1275,1307,1331,1356};
    double voltages[8] = {1000,1225,1250,1275,1292,1325,1350,1375};
    double means[8];
    
    // Set up canvas
    TCanvas *c1 = new TCanvas("C1","C1",1600,1300);
    int color[8] = {1,2,4,6,8,9,46,874};
    int range_low[8]{3,4,4,4,4,5,6,7};
    int range_high[8]={4,18,22,24,27,33,36,40};
    
    // For each file:
    for (int v=0; v<1; v++) {
        
        // Get the waveform fit TTree
        tt2 = (TTree*)files[v]->Get("ptfanalysis2");
        wf2 = new WaveformFitResult;
        if(tt2) wf2->SetBranchAddresses( tt2 );
        
        // Initiliaze histograms
        double x_low = 20;
        double x_hi = 300;
        string hist_name = to_string((int)voltages[v]) + "V";
        TH1F *pc = new TH1F(hist_name.c_str(),"Miu=0.36, Pulse charge distribution at different HV",x_low+x_hi,-1*x_low*bin_unit,x_hi*bin_unit);
        
        // Reset peak-to-valley calculation variables
        ptv_min_amp = 10000;
        ptv_max_amp = 0;
        
        // For each waveform:
        for(int i = 0; i < tt2->GetEntries()-1; i++){
            tt2->GetEvent(i);
            // Collect pulse charge
            auto pulse_charge = wf2->qsum;
            pc->Fill(pulse_charge * 1000.0); // Convert to mV
            
            // if (pulse_charge>=-4 && pulse_charge<=4) miu+=pulse_charge;
            if (i==961240) break;
        }

        // Calculate miu
        double N = 0;
        for (int bin=12; bin<=28; bin++) N+=pc->GetBinContent(bin);
        double miu = -log(N/961240);
        
        // Find peak-to-valley ratio
        // Range (mV*8ns) depends on run
        peakToValley(range_low[v], range_high[v], x_low,pc);
        std::cout << "HV: "<< voltages[v] <<"V. Peak-to-valley ratio: " << ptv_max_amp << "/" << ptv_min_amp << " = " << peak_to_valley << std::endl;;

        // Draw histogram
        gPad->SetLogy();
        if(v==0) {
            pc->GetXaxis()->SetTitle("Pulse charge (mV * 8ns)");
            pc->GetYaxis()->SetTitle("Number of events");
        }
        pc->SetLineColor(color[v]);
        if (v==0) {
            pc->Draw();
        } else {
            pc->Draw("SAMES");
        }
        
        

        TF1 * S_n[6];
        // double init_par[4] = {0};
        double par[8] = { 0,      0,        0.057,    0.1457,       miu,    0,      0,      N};
                        //"Q0",   "sig0",   "W"  ,    "alpha",    "miu",  "sig1", "Q1",   "N"
                        //  0 ,     1   ,     2  ,    3      ,     4   ,  5     , 6   ,   7

        // PRE-FITS:
        // Rough fit for initial parameters
        TF1 *pe_fit = new TF1("pe_fit","gaus",range_low[v]*1.2,range_high[v]*3+v);
        pc->Fit("pe_fit","Q0R");
        par[5]= pe_fit->GetParameter(2);      //1 p.e. gauss rms
        par[6] = pe_fit->GetParameter(1);      //1 p.e. gauss mean
        TF1 *ped_fit = new TF1("ped_fit","gaus",-2,3);
        pc->Fit("ped_fit","Q0R");
        par[1] = ped_fit->GetParameter(2);   //ped gauss rms
        par[0] = ped_fit->GetParameter(1);   //ped gauss mean

        // Pedestal peak
        S_n[0] = new TF1("Sped",Sped,-2,2,5);
        S_n[0]->SetParNames("W","sig0","Q0","miu","N");
        S_n[0]->SetParameters(par[2],par[1],par[0],par[4],par[7]);
        S_n[0]->SetLineColor(1);
        pc->Fit("Sped","QR0");
        par[0] = S_n[0]->GetParameter("Q0");
        par[1] = S_n[0]->GetParameter("sig0");
        par[2] = S_n[0]->GetParameter("W");
        par[7] = S_n[0]->GetParameter("N");

        // 1 p.e. peak
        if (v<1) S_n[1] = new TF1("S1",S1,range_low[v],range_high[v]+1,6);
        else S_n[1] = new TF1("S1",S1,range_low[v]*3,range_high[v]+5,6);
        S_n[1]->SetParNames("sig1","Q1","miu","N","Q0","Qsh");
        S_n[1]->SetParameters(par[5],par[6],par[4],par[7],par[0],par[2]/par[3]);
        S_n[1]->SetLineColor(6);
        pc->Fit("S1", "RQ+");
        par[5] = S_n[1]->GetParameter("sig1");
        par[6] = S_n[1]->GetParameter("Q1");

        // TLine * ptv_hi = new TLine(range_low[v],0,range_low[v],100000);
        // ptv_hi->SetLineColor(color[v]);
        // ptv_hi->Draw("SAME");

        // TLine * ptv_low = new TLine(range_high[v]+1,0,range_high[v]+1,100000);
        // ptv_low->SetLineColor(color[v]);
        // ptv_low->Draw("SAME");

        // Overall fit
        TF1 *pc_f = new TF1("pc_f",fitf,-3,140,8);       //87
        pc_f->SetParNames("Q0","sig0","W","alpha","miu","sig1","Q1","N");
        pc_f->SetParameters(par[0],par[1],par[2],par[3],par[4],par[5],par[6],par[7]);
        // pc_f->SetParameter(0,par[0]);
        // pc_f->SetParameter(1,par[1]);
        // pc_f->FixParameter(2,par[2]);
        // pc_f->FixParameter(3,par[3]);
        // pc_f->SetParameter(4,par[4]);
        // pc_f->SetParameter(5,par[5]);
        // pc_f->SetParameter(6,par[6]);
        // pc_f->SetParameter(7,par[7]);
        pc_f->SetNpx(200);
        pc->Fit("pc_f","0R+");          
        
        // // POST-FITS:
        // // Uncomment to view separate fits
        // for (int save=0; save<8; save++) par[save] = pc_f->GetParameter(save); 

        // TF1 *Sped_test = new TF1("Sped-test",Sped,-3,3,5);
        // Sped_test->SetParNames("W","sig0","Q0","miu","N");
        // Sped_test->SetParameters(par[2],par[1],par[0],par[4],par[7]);
        // Sped_test->SetLineColor(1);
        // Sped_test->Draw("same");

        // TF1 * Snoise_test = new TF1("Snoise-test",Snoise,-2,250,5);
        // Snoise_test->SetParNames("W","alpha","Q0","miu","N");
        // Snoise_test->SetParameters(par[2],par[3],par[0],par[4],par[7]);
        // Snoise_test->SetLineColor(4);
        // Snoise_test->Draw("same");

        // TF1 * S1_test = new TF1("S1-test",S1,0,50,6);
        // S1_test->SetParNames("sig1","Q1","miu","N","Q0","Qsh");
        // S1_test->SetParameters(par[5],par[6],par[4],par[7],par[0],par[2]/par[3]);
        // S1_test->SetLineColor(6);
        // S1_test->Draw("same");

        //  // 2+ p.e. peak
        // for (int n=2; n<=10; n++) {
        //     string fitname_test = "S" + to_string(n) = "-test";
        //     TF1* Sn_test = new TF1(fitname_test.c_str(),Sn,4,250,7);
        //     Sn_test->SetParNames("sig1","Q1","miu","N","n","Q0","Qsh");
        //     Sn_test->SetParameters(par[5],par[6],par[4],par[7],n,par[0],par[2]/par[3]);
        //     Sn_test->SetLineColor(8);
        //     Sn_test->Draw("same");
        // }

        gPad->Update();
        gStyle->SetOptStat(11);        
        gStyle->SetOptFit();

        means[v]=pc_f->GetParameter("Q1");
       
        c1->Update();
        
        // Print peak-to-valley ratio on histogram
        TPaveStats *ps = (TPaveStats*)pc->GetListOfFunctions()->FindObject("stats");
        ps->SetName("peak-to-valley");
        TList *listOfLines = ps->GetListOfLines();
        string text = "Peak-to-valley   " + to_string(peak_to_valley).substr(0,4);
        TLatex *myt = new TLatex(0,0, text.c_str());
        listOfLines->Add(myt);
        pc->SetStats(0);
        ps->SetX1NDC(0.25+v*0.15);ps->SetX2NDC(0.40+v*0.15);
        ps->SetY1NDC(0.7);ps->SetY2NDC(0.9);
        if (v>=5) {
            ps->SetX1NDC(0.25+(v-2)*0.15);ps->SetX2NDC(0.40+(v-2)*0.15);
            ps->SetY1NDC(0.5);ps->SetY2NDC(0.7);
        }
        ps->SetTextColor(color[v]);

        c1->Modified();
        
//        files[v]->Close();            // including this deletes the histograms
    }
    
    c1->SaveAs("hv_pc_spectrum.png");
    
    TCanvas *c2 = new TCanvas("C2");
    TGraph *pc_mean = new TGraph(7,voltages,means);
    pc_mean->SetTitle("1PE peak charge at different HV");
    pc_mean->GetXaxis()->SetTitle("HV (V)");
    pc_mean->GetYaxis()->SetTitle("1PE peak charge (mV*8ns)");
    pc_mean->GetYaxis()->SetRangeUser(9,100);
    gPad->SetLogy();
    pc_mean->Draw("a*");
    c2->SaveAs("hv_pc_spectrum_mean.png");
    
    return 0;
}