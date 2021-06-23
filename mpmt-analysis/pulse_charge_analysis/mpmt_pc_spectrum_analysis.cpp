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
#include <math.h> 
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
        if (pc_count<ptv_min_amp) {                     // find min amp
            ptv_min_amp=pc_count;
            continue;
        }
        if (pc_count>ptv_max_amp) ptv_max_amp=pc_count; // find max amp
    }
    peak_to_valley = ptv_max_amp/ptv_min_amp;           // find peak-to-valley ratio
}

// // Define fit function for pulse charge distribution
// Double_t fitf(Double_t *x, Double_t *p) {
//     Double_t    Sped=0,
//                 Snoise=0,
//                 S1=0,
//                 Sn=0;

//     Sped=(1-p[2])/(sqrt(2*M_PI)*p[1])*exp(-0.5*pow((x[0]-p[0])/(p[1]),2)-p[4]);
//     if (x[0] > p[0]) Snoise=p[3]*p[2]*exp(-1*p[3]*(x[0]-p[0])-p[4]);
//     S1=p[4]/(sqrt(2*M_PI)*p[5])*exp(-0.5*pow((x[0]-p[6]-(p[2]/p[3]))/(p[5]),2)-p[4]);

//     for (int n=2; n<=4; n++) {
//         Sn+=pow(p[4],n)*exp(-1*p[4])/(sqrt(2*M_PI*n)*p[5])*exp(-0.5*pow((x[0]-p[0]-(n*p[6])-(p[2]/p[3]))/(p[5]),2)/n);
//     }
    
//     return p[7]*(Sped+Snoise+S1+Sn);
// }

Double_t Sped(Double_t *x, Double_t *p) {
    Double_t Sped = p[4]*(1-p[0])/(sqrt(2*M_PI)*p[1])*exp(-0.5*pow((x[0]-p[2])/(p[1]),2)-p[3]);
    return Sped;
}

Double_t Snoise(Double_t *x, Double_t *p) {
    Double_t Snoise = 0;
    if (x[0] > p[2]) Snoise=p[1]*p[0]*exp(-1*p[1]*(x[0]-p[2])-p[3]);    //p[4]*
    return Snoise;
}

Double_t S1(Double_t *x, Double_t *p) {
    // Double_t S1 = p[4]/(sqrt(2*M_PI)*p[2])*exp(-0.5*pow((x[0]-p[3]-(p[0]/p[1]))/p[2],2)-p[4]);
    Double_t S1 = p[3]*p[2]/(sqrt(2*M_PI)*p[0])*exp(-0.5*pow((x[0]-p[1])/p[0],2)-p[2]);

    return S1;
}

Double_t Sn(Double_t *x, Double_t *p) {
    Double_t Sn = 0;
    int fact = 1;
    int scale = 10;
    // for (int n=3; n<=4; n++) {
    fact *= p[4];
    scale = pow(scale,(p[4]-1));
    Sn += p[3]*pow(p[2],p[4])/(sqrt(2*M_PI*p[4])*p[0]*fact*scale)*exp(pow((x[0]-p[4]*p[1])/p[0],2)/(-2*p[4])-p[2]);
    // }
    
    return Sn;
}

// Double_t S3(Double_t *x, Double_t *p) {
//     Double_t S3 = p[3]*pow(p[2],3)/(sqrt(2*M_PI*3)*p[0]*2*3)*exp(-0.5*pow((x[0]-3*p[1])/p[0],2)-p[2])/3;
//     return S3;
// }

// Double_t fitf(Double_t *x, Double_t *p) {
//     return Sped(x,p)+Snoise(x,&p[5])+S1(x,&p[10])+Sn(x,&p[14]);
// }


int main( int argc, char* argv[] ) {
    
    // Set up files
    TFile * files[8];
    files[0] = new TFile( "../../mpmt_Analysis_run0861.root" , "read" );
    files[1] = new TFile( "../../mpmt_Analysis_run0854.root" , "read" );
    files[2] = new TFile( "../../mpmt_Analysis_run0853.root" , "read" );
    files[3] = new TFile( "../../mpmt_Analysis_run0857.root" , "read" );
    files[4] = new TFile( "../../mpmt_Analysis_run0855.root" , "read" );
    files[5] = new TFile( "../../mpmt_Analysis_run0856.root" , "read" );
    files[6] = new TFile( "../../mpmt_Analysis_run0862.root" , "read" );
    files[7] = new TFile( "../../mpmt_Analysis_run0863.root" , "read" );
    
    // Set up voltages
    double voltages[8] = {1209,1234,1258,1275,1307,1331,1356,1478};
    double means[8];
    
    // Set up canvas
    TCanvas *c1 = new TCanvas("C1","C1",1600,1300);
    int color[8] = {1,810,632,616,600,882,417,861}; //black, orange, red, magenta, blue, violet, green
    int range_high[8] = {13,15,18,21,22,24,26,28};
    
    // For each file:
    for (int v=0; v<1; v++) {       //7
        
        // Get the waveform fit TTree
        tt2 = (TTree*)files[v]->Get("ptfanalysis2");
        wf2 = new WaveformFitResult;
        if(tt2) wf2->SetBranchAddresses( tt2 );
        
        // Initiliaze histograms
        double x_low = 20;
        string hist_name = to_string((int)voltages[v]) + "V";
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
            
            if (i==961240) break;
        }
        
        // Find peak-to-valley ratio
        // Range (mV*8ns) depends on run
        int range_low = 3;
        peakToValley(range_low, range_high[v], x_low,pc);
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
        
        TF1 *pe_fit = new TF1("pe_fit","gaus",5,range_high[v]+4);
        pc->Fit("pe_fit","Q0R");
        TF1 *noise_fit = new TF1("noise_fit","gaus",-2,2);
        pc->Fit("noise_fit","Q0R+");

        // c1->Update();

        // TF1 *pc_f = new TF1("pc_f",fitf,-3,35,8);       //87
        // pc_f->SetParNames("Q0","sig0","W","alpha","miu","sig1","Q1","N");
        // pc_f->SetParameters(noise_fit->GetParameter(1),
        //                     noise_fit->GetParameter(2),
        //                     0.35,
        //                     0.06,
        //                     1,
        //                     pe_fit->GetParameter(2),
        //                     pe_fit->GetParameter(1),
        //                     2000);
        // pc_f->FixParameter(0,noise_fit->GetParameter(1));
        // pc_f->FixParameter(1,noise_fit->GetParameter(2));
        // pc_f->SetParLimits(2,0,0.6);
        // pc_f->SetParLimits(3,0,0.1);
        // pc_f->FixParameter(4,1);
        // pc_f->FixParameter(5,pe_fit->GetParameter(2));
        // pc_f->SetParLimits(6,pe_fit->GetParameter(1)-5,pe_fit->GetParameter(1)+5);
        // pc_f->SetParLimits(7,1000,3000);
        // pc->Fit("pc_f","","",5,16);  


        TF1 * S_ped = new TF1("Sped",Sped,-3,3.5,5);
        S_ped->SetParNames("W","sig0","Q0","miu","N");
        S_ped->SetParLimits(0,0,1);
        S_ped->SetParameter(1,noise_fit->GetParameter(2));
        S_ped->SetParameter(2,noise_fit->GetParameter(1));
        S_ped->SetParameter(3,1.4);
        S_ped->SetParameter(4,75000);
        S_ped->SetLineColor(2);
        pc->Fit("Sped","R+");

        // TF1 * S_noise = new TF1("Snoise",Snoise,-3,60,5);
        // S_noise->SetParNames("W","alpha","Q0","miu","N");   //,"N"
        // S_noise->SetParameter(0,S_ped->GetParameter("W"));
        // S_noise->SetParameter(1,0.04);
        // S_noise->FixParameter(2,S_ped->GetParameter("Q0"));
        // S_noise->SetParameter(3,-8.8);
        // S_noise->SetParameter(4,-0.46);  //,2000
        // S_noise->SetLineColor(7);
        // pc->Fit("Snoise","Q0R+");

        TF1 * S_1 = new TF1("S1",S1,3.5,range_high[v]+10,4);      //range_high[v]+3
        S_1->SetParNames("sig1","Q1","miu","N");
        S_1->SetParameter(0,pe_fit->GetParameter(2));
        S_1->SetParameter(1,pe_fit->GetParameter(1));
        S_1->SetParameter(2,1.3);
        S_1->SetParameter(3,75000);
        S_1->SetLineColor(6);
        pc->Fit("S1", "R+");

        TF1 * S_n[5];
        for (int n=2; n<=4; n++) {
            S_n[n] = new TF1("Sn",Sn,4,80,5);
            S_n[n]->SetParNames("sig1","Q1","miu","N","n");
            S_n[n]->SetParameter(0,S_1->GetParameter(0));
            S_n[n]->FixParameter(1,S_1->GetParameter(1));
            S_n[n]->SetParameter(2,1.3);
            S_n[n]->FixParameter(3,75000);
            S_n[n]->FixParameter(4,n);
            S_n[n]->SetLineColor(3);
            pc->Fit("Sn", "R+");
        }


//how do i know for sure if the separate fits
        // // TF1 * fitf = new TF1("total_fit", "Sped(0)+Snoise(5)+S1(10)+Sn(14)",-3,50);
        // TF1* fit_f = new TF1("fit_f",fitf,-3,80,18);
        // Double_t par[18];
        // S_ped->GetParameters(&par[0]);
        // S_noise->GetParameters(&par[5]);
        // S_1->GetParameters(&par[10]);
        // S_n->GetParameters(&par[14]);
        // fit_f->SetParameters(par);
        // pc->Fit("fit_f","R+");
        

        // S_ped->SetParNames("Q0","sig0","W","alpha","miu","sig1","Q1","N");
        // S_noise->SetParNames("Q0","sig0","W","alpha","miu","sig1","Q1","N");
        // S_ped->SetParameters(noise_fit->GetParameter(1),
        //                     noise_fit->GetParameter(2),
        //                     0.35,
        //                     0.06,
        //                     0.53,
        //                     pe_fit->GetParameter(2),
        //                     pe_fit->GetParameter(1),
        //                     2000);
        // // S_ped->SetParLimits(2,0,1);
        // // S_ped->SetParLimits(4,0,2);
        // S_noise->SetParameters(noise_fit->GetParameter(1),
        //                     noise_fit->GetParameter(2),
        //                     0.35,
        //                     0.06,
        //                     0.53,
        //                     pe_fit->GetParameter(2),
        //                     pe_fit->GetParameter(1),
        //                     2000);
            
        // try minimizing range
        
        gPad->Update();
        gStyle->SetOptStat(11);        
         gStyle->SetOptFit();
        // TPaveStats *fit_stats = (TPaveStats*)pc->GetListOfFunctions()->FindObject("stats");
        // fit_stats->SetOptStat();
        // pc->SetStats(1);


        means[v]=S_1->GetParameter(1);
       
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
        ps->SetY1NDC(0.65);ps->SetY2NDC(0.85);
        if (v>=5) {
            ps->SetX1NDC(0.25+(v-2)*0.15);ps->SetX2NDC(0.40+(v-2)*0.15);
            ps->SetY1NDC(0.45);ps->SetY2NDC(0.65);
        }
        ps->SetTextColor(color[v]);
        
    
        c1->Modified();
        
//        files[v]->Close();            // including this deletes the histograms
    }
    
    c1->SaveAs("pc_spectrum.png");
    
    TCanvas *c2 = new TCanvas("C2");
    TGraph *pc_mean = new TGraph(7,voltages,means);
    pc_mean->SetTitle("Pulse charge p.e. peak mean at different HV");
    pc_mean->GetXaxis()->SetTitle("HV (V)");
    pc_mean->GetYaxis()->SetTitle("Mean of pulse charge p.e. peak (mV*8ns)");
    gPad->SetLogy();
    pc_mean->Draw("a*");
    c2->SaveAs("pc_spectrum_mean.png");
    
    return 0;
}
