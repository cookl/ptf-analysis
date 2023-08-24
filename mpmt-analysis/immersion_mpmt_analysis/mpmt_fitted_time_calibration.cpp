// Hamza Dastgir mPMT Time Calibration
// hdastgir@triumf.ca
// 08/15/23


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
#include "TMath.h"

#include <iostream>
#include <vector>
#include <string>
#include <cstring>

#include <algorithm>
#include <numeric>
#include <cmath>
#include <map>
//using namespace std;

//This will return 3 types of histograms for 1 signal channel and the reference channel
//1:a fitted pulse CFD time difference histogram with FWHM calculations
//2:all pulse charges for signal pulses that are within 100ns of any reference channel pulse
//3:all pulse charges for only signal pulses that are within 100ns of the reference pulse (earliest pulse above charge cut)
//will also produce a scatter plot of FWHM vs mean pulse charge for the signal channel
//one can change the chosen channel using line 159 and the reference pulse charge cut using line 162


class Pulse {
    public:
        int event;
        int channel;
        double height;
        double charge;
        double CFDtime;
        double time;
        double fittedtime;
        Pulse(int aEvent, int aChannel, double aHeight, double aCharge, double aCFDTime, double aTime, double aFittedTime){
            event = aEvent;
            channel = aChannel;
            height = aHeight;
            charge = aCharge;
            CFDtime = aCFDTime;
            time = aTime;
            fittedtime = aFittedTime;
        }
};

bool notSame(Pulse pulse1, Pulse pulse2){ //determines if 2 pulses are the same
    if((pulse1.channel == pulse2.channel) && (pulse1.CFDtime == pulse2.CFDtime)){
        return false;
    } else {
        return true;
    }
}
double fAssymGauss(double *x, double *p)//skewed gaussian fit function
{
    double f=0,
                    y1=0,
                    y2=0,
                    stot=0,
                    n=0;
    y1 = x[0]-p[1];
    y1*=y1;
    y1/=2*p[2]*p[2];

    y2=x[0]-p[1];
    y2*=y2;
    y2/=2*p[3]*p[3];

    stot=p[2]+p[3];
    n=2.;
    n/=stot;
    n*=p[0];
    if( x[0]<p[1] ){ f+=TMath::Exp( -y1); }
    else{ f+=TMath::Exp( -y2 ); }

    return n*f;


}
int time_diff = 100;//time cut for filtering charges



std::vector<Pulse> isCoincidentEvent(std::vector<Pulse> ref, std::vector<Pulse> lst){ // checks lst against ref(reference pulse)
    std::vector<Pulse> coincident_pulses = {};                                        // to see if time diff<100ns, vector of pulses
    Pulse refpulse = ref[0];
    for(Pulse pulse : lst){
        if((std::abs(refpulse.time - pulse.time) < time_diff) && (notSame(refpulse, pulse))){
            coincident_pulses.push_back(pulse);
            break;
        }
    }
    return coincident_pulses;


}

std::vector<Pulse> timeCheck(std::vector<Pulse> lst1 , std::vector<Pulse> lst2){//checks lst1 against lst2 to see if time diff <100ns
    std::vector<Pulse> outPulses = {};                                          //returns pulses that pass
    for(Pulse pulse1 : lst1){
        std:: vector<Pulse> checkPulse = {};
        for(Pulse pulse2 : lst2){
            if(((pulse1.time - pulse2.time) < time_diff) && checkPulse.empty()){
             checkPulse.push_back(pulse1);
            } else {
                continue;
            }
        }
        if(checkPulse.empty()){
            continue;
        } else {
            outPulses.push_back(checkPulse[0]);
        }
    }
    return outPulses;
}

int main( int argc, char* argv[] ) {

    if ( argc != 2 ){
        std::cerr<<"Usage: ptf_ttree_analysis.app ptf_analysis.root\n";
        exit(0);
    }

    TFile * fin = new TFile( argv[1], "read" );
    TTree * tt[20];
    WaveformFitResult * wf[20];
    TH1F * h1[20]; //all pulse charges within time diff of any pulse within ref channel (not only ref pulse)
    TH1F * h2[20]; //pulse time cfd diff for pulses with charge > 1.5 PE
    TH1F * h14[20]; //pulse charges filtered to be within time_diff of reference pulse(earliest pulse above charge cut in ref channel)
    TH1F * h15[20]; //pulse time cfd diff for pulses with charge < 1.5 PE
    TH1F * h3[20]; //FWHM of time_diff vs pulse charge
    std::map<std::string, std::vector<int>> channel_change{
        //key -- letter: V = U-Vic driver, T= Tamadenshi Laser -- number = channel -- ex: V14 = U-Vic ch14
        //value -- vector<int> = {channel number, lower bound of hist, upper bound of hist}
        {"V7", {7,-3,6}}, {"V10", {10,-20,-11}}, {"V11", {11,-21,-12}}, {"V12", {12,-11,-2}},
        {"V13", {13,-15,-6}}, {"V14", {14,-8,1}}, {"V15", {15,-9,0}}, {"V16", {16,-11,-2}}, {"V17", {17,-14,-5}},
        {"T7", {7,-4,5}}, {"T10", {10,-13,-4}}, {"T11", {11,-14,-5}}, {"T12", {12,-11,-2}},
        {"T13", {13,-15,-6}}, {"T14", {14,-9,0}}, {"T15", {15,-9,0}}, {"T16", {16,-11,-2}}, {"T17", {17,-5,4}}

    };

    std::string chosen_channel = "V13"; //IMPORTANT use this to change the channel and bounds of the histograms in accordance to the map above
                                       //to add more drivers or channels, add to the map above

    int ref_charge_cut = 4; //charge cut for reference pulses, change this value to change ref pulse cut

    std::vector<int> bounds = channel_change.find(chosen_channel)->second; //finds value(vector<int>) using key, saves as bounds




    for(int j = 0; j < 20; j++){ // To initialize TTree and branches

        char branch_name[100];
        sprintf(branch_name,"ptfanalysis%i",j);
        tt[j] = (TTree*)fin->Get(branch_name);

        wf[j] = new WaveformFitResult;
        if(tt[j]) wf[j]->SetBranchAddresses( tt[j] );

        h14[j] = new TH1F("pulse_chargefil", "Pulse Charge Filtered",100,0,20);
        h1[j] = new TH1F("pulse_charge", "Pulse Charge",100,0,20);
        h3[j] = new TH1F("FWHM_vs_mean_pulse_charge", "FWHM vs Mean Pulse Charge",100,bounds[1],bounds[2]);




        // Zooming in on Different Sections of the Histogram Depending on Channel Number
        switch(j){
        case 0:
            h2[j] = new TH1F("pulse_time_CFD_difference_PE > 1.5", "Pulse Time CFD Difference",100,bounds[1],bounds[2]);
            break;
        case 1:
            h2[j] = new TH1F("pulse_time_CFD_difference_PE > 1.5", "Pulse Time CFD Difference",100,-2,1);
            break;
        }
        switch(j){
        case 0:
            h15[j] = new TH1F("pulse_time_CFD_difference_PE < 1.5", "Pulse Time CFD Difference",100,bounds[1],bounds[2]);
            break;
        case 1:
            h15[j] = new TH1F("pulse_time_CFD_difference_PE < 1.5", "Pulse Time CFD Difference",100,-2,1);
            break;

        }
    }

    int myChannels[] = {bounds[0], 2};

    int chansize = std::size(myChannels);

    for(int i = 0; i < tt[0]->GetEntries()-1; i++){ // loop over events

        std::cout << "Event: " << i << std::endl;

        std::vector<Pulse> sameEventPulse;
        std::vector<Pulse> refPulse;
        std::vector<Pulse> ch2Pulse;
        for(int j = 0; j < chansize; j++){ // loop over channels
            int chan = myChannels[j];
            tt[chan]->GetEvent(i);

            std::cout << "Channel: " << chan << " Number of Pulses: " << wf[chan]->numPulses << std::endl;

            std::vector<Pulse> sameChannel = {};
            for(int k = 0; k < wf[chan]->numPulses; k++){ // loop over each pulse

                // pulse height
                double pulse_height = wf[chan]->pulseCharges[k]*1000.0;

                // pulse charge
                double pulse_charge = pulse_height/7.7;

                // pulse CFDtime
                double pulse_CFDtime = wf[chan]->pulseTimesCFD[k];

                // pulse time
                double pulse_time = wf[chan]->pulseTimes[k];

                //pulse fittedTime
                double pulse_fittedtime = wf[chan]->mean;

                Pulse pulse(i, chan, pulse_height, pulse_charge, pulse_CFDtime, pulse_time, pulse_fittedtime);
                sameChannel.push_back(pulse);
                std::cout << "CFD Time: " <<  pulse_CFDtime << " Channel: " << chan << std::endl;

            }
            if(sameChannel.size() == 0){
                continue;
            }
            Pulse REF = sameChannel[0];

            if(chan == 2){
                for(Pulse pulse : sameChannel){
                    if((pulse.charge > ref_charge_cut)){// && pulse.charge < 1.5){ //implementing charge cut to find ref pulse above cut
                        REF = pulse;
                    }
                }
            }

            for(Pulse pulse : sameChannel){ //looping over ch2 (reference channel) to find earliest pulse above the charge cut
                if(pulse.channel == 2){
                    ch2Pulse.push_back(pulse);
                    if((pulse.time < REF.time) && (pulse.charge > ref_charge_cut)){ //&& pulse.charge < 1.5){
                        REF = pulse;
                    } else {
                        continue;
                    }
                }
                else{
                    sameEventPulse.push_back(pulse);
                }
            }

            if ((REF.channel == 2) && (REF.charge > ref_charge_cut)){// && REF.charge < 1.5){ //ensuring the final ref pulse is above charge cut
                refPulse.push_back(REF);
        }
        }
        std::vector<Pulse> reftimepulses = timeCheck(ch2Pulse,sameEventPulse); // finding pulses to fill in h1
        std::vector<Pulse> sigtimepulses = timeCheck(sameEventPulse,ch2Pulse);
        for(Pulse pulse : sigtimepulses){
            h1[0]->Fill(pulse.charge);
        }
        for(Pulse pulse : reftimepulses){
            h1[1]->Fill(pulse.charge);
        }


        if(refPulse.empty() || sameEventPulse.empty()){
            continue;
        } else {
            std::vector<Pulse> filteredPulse = isCoincidentEvent(refPulse,sameEventPulse);//filtering signal pulses using ref pulse
            Pulse ReferencePulse = refPulse[0];                                           //using time_diff
            h14[1]->Fill(ReferencePulse.charge);
            for(Pulse pulse : filteredPulse){ //filling h14 using filtered pulses, and filling h3 based on pulse charge
                h14[0]->Fill(pulse.charge);   //and filling h2 and h15 using pulse charges
                if((pulse.charge >= 0.5) && (pulse.charge < 1.5)){
                    h15[0]->Fill(pulse.CFDtime - ReferencePulse.CFDtime);//finding time difference between signal pulse and ref pulse
                    h3[0]->Fill(pulse.CFDtime - ReferencePulse.CFDtime); //and filling histograms h2,h3,h15 with time diff
                } else if((pulse.charge >= 1.5) && (pulse.charge < 2.5)){
                    h3[1]->Fill(pulse.CFDtime - ReferencePulse.CFDtime);
                    h2[0]->Fill(pulse.CFDtime - ReferencePulse.CFDtime);
                } else if((pulse.charge >= 2.5) && (pulse.charge < 3.5)){
                    h3[2]->Fill(pulse.CFDtime - ReferencePulse.CFDtime);
                    h2[0]->Fill(pulse.CFDtime - ReferencePulse.CFDtime);
                } else if((pulse.charge >= 3.5) && (pulse.charge < 4.5)){
                    h3[3]->Fill(pulse.CFDtime - ReferencePulse.CFDtime);
                    h2[0]->Fill(pulse.CFDtime - ReferencePulse.CFDtime);
                } else if((pulse.charge >= 4.5) && (pulse.charge < 5.5)){
                    h3[4]->Fill(pulse.CFDtime - ReferencePulse.CFDtime);
                    h2[0]->Fill(pulse.CFDtime - ReferencePulse.CFDtime);
                } else if((pulse.charge >= 5.5) && (pulse.charge < 6.5)){
                    h3[5]->Fill(pulse.CFDtime - ReferencePulse.CFDtime);
                    h2[0]->Fill(pulse.CFDtime - ReferencePulse.CFDtime);
                } else if((pulse.charge >= 6.5) && (pulse.charge < 7.5)){
                    h3[6]->Fill(pulse.CFDtime - ReferencePulse.CFDtime);
                    h2[0]->Fill(pulse.CFDtime - ReferencePulse.CFDtime);
                } else if((pulse.charge >= 7.5) && (pulse.charge < 8.5)){
                    h3[7]->Fill(pulse.CFDtime - ReferencePulse.CFDtime);
                    h2[0]->Fill(pulse.CFDtime - ReferencePulse.CFDtime);
                } else if((pulse.charge >= 8.5) && (pulse.charge < 9.5)){
                    h3[8]->Fill(pulse.CFDtime - ReferencePulse.CFDtime);
                    h2[0]->Fill(pulse.CFDtime - ReferencePulse.CFDtime);
                } else { continue;
                }

            }
        }
    }



    // Print Pulse Time CFD Histograms


    // Print Pulse Time Difference Relative to Channel 2 Histograms
    for(int j = 0; j < chansize; j++){ //overlapping h15 and h2 using the fit function to find the FWHM
        TCanvas * c2 = new TCanvas("C2");

        char hist_title[100];
        sprintf(hist_title,"Channel %d CFD Pulse Time Differences",myChannels[j]);
        h2[j]->SetTitle(hist_title);
        h15[j]->SetTitle(hist_title);
        TF1 *func = new TF1("fit",fAssymGauss,((h2[j]->GetMean()) - 2),((h2[j]->GetMean()) + 2),4);//fit for h2
        func->SetParameters(1,h2[j]->GetMean(),h2[j]->GetRMS(),h2[j]->GetRMS(),2);
        func->SetParNames("Constant","Mean_value","Sigma", "Sigma2");
        func->SetLineColor(kBlue);
        h2[j]->Fit("fit", "R");
        h2[j]->GetXaxis()->SetTitle("Pulse Times (ns)");
        h2[j]->GetYaxis()->SetTitle("Number of events");
        TF1 *func2 = new TF1("fit2",fAssymGauss,((h15[j]->GetMean()) - 2.1),((h15[j]->GetMean()) + 2.1),4);//fit for h15
        func2->SetParameters(1,h15[j]->GetMean(),h15[j]->GetRMS(),h15[j]->GetRMS(),2);
        func2->SetParNames("Constant","Mean_value","Sigma", "Sigma2");
        h15[j]->Draw();
        h15[j]->Fit("fit2", "R");
        double Sconst = func2->GetParameter(0);//h15 fit stats
        double Smean = func2->GetParameter(1);
        double Ssigma1 = func2->GetParameter(2);
        double Ssigma2 = func2->GetParameter(3);
        double Sscale = 2/(Ssigma1+Ssigma2)*Sconst;
        double Shalf_max = Sscale/2;
        double Sleft_fwhm = -sqrt(2)*Ssigma1*sqrt(-log(Shalf_max/Sscale)) + Smean;
        double Sright_fwhm = sqrt(2)*Ssigma2*sqrt(-log(Shalf_max/Sscale)) + Smean;
        double Sfwhm = Sright_fwhm - Sleft_fwhm;
        double Sx_val[2];
        double Sy_val[2];
        Sx_val[0] = Sleft_fwhm;
        Sx_val[1] = Sright_fwhm;
        Sy_val[0] = Shalf_max;
        Sy_val[1] = Shalf_max;
        double Lconst = func->GetParameter(0);
        double Lmean = func->GetParameter(1);//h2 fit stats
        double Lsigma1 = func->GetParameter(2);
        double Lsigma2 = func->GetParameter(3);
        double Lscale = 2/(Lsigma1+Lsigma2)*Lconst;
        double Lhalf_max = Lscale/2;
        double Lleft_fwhm = -sqrt(2)*Lsigma1*sqrt(-log(Lhalf_max/Lscale)) + Lmean;
        double Lright_fwhm = sqrt(2)*Lsigma2*sqrt(-log(Lhalf_max/Lscale)) + Lmean;
        double Lfwhm = Lright_fwhm - Lleft_fwhm;
        double Lx_val[2];
        double Ly_val[2];
        Lx_val[0] = Lleft_fwhm;
        Lx_val[1] = Lright_fwhm;
        Ly_val[0] = Lhalf_max;
        Ly_val[1] = Lhalf_max;


        h15[j]->GetMean();

        h2[j]->Draw();
        double h2_height = h2[j]->GetMaximum();
        h15[j]->Draw();
        double h15_height = h15[j]->GetMaximum();
        if(h15_height >= h2_height){
            h2[j]->Draw("SAMES");
        } else {
            h2[j]->Draw();
            h15[j]->Draw("SAMES");
        }

        h15[j]->SetLineColor(kRed);
        TGraph *SFWHM = new TGraph(2,Sx_val,Sy_val);
        SFWHM->SetLineColor(kRed);
        SFWHM->Draw("SAME");
        TGraph *LFWHM = new TGraph(2,Lx_val,Ly_val);
        LFWHM->SetLineColor(kBlue);
        LFWHM->Draw("SAME");
        char Sstat[100];
        sprintf(Sstat,"FWHM PE<1.5=%f", Sfwhm);
        char Lstat[100];
        sprintf(Lstat,"FWHM PE>1.5=%f", Lfwhm);
        TLegend *leg2 = new TLegend(0.1,0.7,0.3,0.9);
        leg2->AddEntry(h15[j],"Pulse Charge PE < 1.5","l");
        leg2->AddEntry(SFWHM,Sstat,"l");
        leg2->AddEntry(h2[j],"Pulse Charge PE > 1.5","l");
        leg2->AddEntry(LFWHM,Lstat,"l");
        leg2->Draw();
        h15[j]->GetXaxis()->SetTitle("Pulse Times (ns)");
        h15[j]->GetYaxis()->SetTitle("Number of events");
        gStyle->SetOptFit(11);
        gPad->Update();

        TPaveStats * st = (TPaveStats*)h2[j]->FindObject("stats"); //moving stat box so they dont overlap
        st->SetX1NDC(0.78); //new x start position
        st->SetX2NDC(.98); //new x end position
        st->SetY1NDC(0.5); //new y start position
        st->SetY2NDC(0.7); //new y end position
        TPaveStats * st2 = (TPaveStats*)h15[j]->FindObject("stats");
        st2->SetX1NDC(0.78); //new x start position
        st2->SetX2NDC(.98); //new x end position
        st2->SetY1NDC(0.7); //new y start position
        st2->SetY2NDC(0.9); //new x end position
        char png_name[100];
        sprintf(png_name,"pulse_CFD_time_difference_CHANNEL2_%d_LED.png",myChannels[j]);
        c2->SaveAs(png_name);
    }
    double fwhmvar[9] = {};
    // Print Pulse Time Difference Relative to Channel 2 Histograms
    for(int j = 0; j < 9; j++){//fitting h3 histograms to find FWHM for scatter plot
        TCanvas * c3 = new TCanvas("C3");

        char hist_title[100];
        sprintf(hist_title,"Channel %d CFD Pulse Time Differences",myChannels[j]);
        h3[j]->SetTitle(hist_title);
        TF1 *func = new TF1("fit",fAssymGauss,((h3[j]->GetMean()) - 2.2),((h3[j]->GetMean()) + 2.2),4);
        func->SetParameters(1,h3[j]->GetMean(),h3[j]->GetRMS(),h3[j]->GetRMS(),2);
        func->SetParNames("Constant","Mean_value","Sigma", "Sigma2");
        func->SetLineColor(kBlue);
        h3[j]->Fit("fit", "R");
        h3[j]->GetXaxis()->SetTitle("Pulse Times (ns)");
        h3[j]->GetYaxis()->SetTitle("Number of events");
        double Lconst = func->GetParameter(0);
        double Lmean = func->GetParameter(1);
        double Lsigma1 = func->GetParameter(2);
        double Lsigma2 = func->GetParameter(3);
        double Lscale = 2/(Lsigma1+Lsigma2)*Lconst;
        double Lhalf_max = Lscale/2;
        double Lleft_fwhm = -sqrt(2)*Lsigma1*sqrt(-log(Lhalf_max/Lscale)) + Lmean;
        double Lright_fwhm = sqrt(2)*Lsigma2*sqrt(-log(Lhalf_max/Lscale)) + Lmean;
        double Lfwhm = Lright_fwhm - Lleft_fwhm;
        if(Lfwhm > 5){
            fwhmvar[j] =0;
        }else{ fwhmvar[j] = Lfwhm;
        }
    }
    TCanvas * c4 = new TCanvas("C4");//making scatterplot from FWHM found from h3 histograms
    char chname[100];
    double charges[9] = {1,2,3,4,5,6,7,8,9};
    sprintf(chname,"FWHM vs Mean Pulse Charge ch %d",myChannels[0]);
    TGraph * g4 = new TGraph(9, charges, fwhmvar);
    g4->SetTitle(chname);
    g4->GetXaxis()->SetTitle("Mean Pulse Charge (PE)");
    g4->GetYaxis()->SetTitle("FWHM");
    gStyle->SetOptFit(11);
    g4->SetMarkerColor(4);
    g4->Draw("ap*");
    char png_names[100];
    sprintf(png_names,"ch%d_FWHM_vs_MeanPulseCharge.png",myChannels[0]);
    c4->SaveAs(png_names);


    // Print Pulse Charge Histograms
    for(int j = 0; j < chansize; j++){
        TCanvas * c14 = new TCanvas("C14");

        char hist_title[100];
        sprintf(hist_title,"Channel %d Pulse Charges Filtered",myChannels[j]);
        h14[j]->SetTitle(hist_title);
        h14[j]->Draw();
        h14[j]->GetXaxis()->SetTitle("Pulse Charge (PE)");
        h14[j]->GetYaxis()->SetTitle("Number of events");
        gStyle->SetOptFit(11);

        char png_name[100];
        sprintf(png_name,"pulse_charges_fil%d_LED.png",myChannels[j]);
        c14->SaveAs(png_name);
    }
    // Print Pulse Charge Histograms
    for(int j = 0; j < chansize; j++){
        TCanvas * c1 = new TCanvas("C1");

        char hist_title[100];
        sprintf(hist_title,"Channel %d Pulse Charges",myChannels[j]);
        h1[j]->SetTitle(hist_title);
        h1[j]->Draw();
        h1[j]->GetXaxis()->SetTitle("Pulse Charge (PE)");
        h1[j]->GetYaxis()->SetTitle("Number of events");
        gStyle->SetOptFit(11);

        char png_name[100];
        sprintf(png_name,"pulse_charges_%d_LED.png",myChannels[j]);
        c1->SaveAs(png_name);
    }
    std::cout << "Number of Events: " << tt[0]->GetEntries() << std::endl;

    fin->Close();
    return 0;
}

