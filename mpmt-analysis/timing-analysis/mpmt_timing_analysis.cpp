#include "WaveformFitResult.hpp"
#include "ScanPoint.hpp"
#include "TFile.h"
#include "TCanvas.h"
#include "TH1D.h"
#include "TH2F.h"
#include "TGraph.h"
#include "TF1.h"
#include "TStyle.h"
#include "TLegend.h"
#include "TProfile.h"

#include <iostream>
#include <vector>
double funcEMG(double *x, double *p){
  //Exponential gaussian used for fitting
  // p[0]: amplitude
  // p[1]: gaussian mu
  // p[2]: gaussian sig
  // p[3]: exponential decay constant
  // p[4]: baseline
  
 double y = p[4] + (p[0]/0.3)*(p[3]/2.)*exp((p[1]+p[2]*p[2]/p[3]/2.-x[0])/(p[3]))*
   TMath::Erfc((p[1]+p[2]*p[2]/p[3] -x[0])/sqrt(2.)/p[2]) *
   TMath::Sin(p[5] * (x[0]-p[6]));

 return y ;
}

double bessel(double *x, double *p){
  //Exponential gaussian used for fitting
  double xx = (x[0] - p[1]) * p[0];

  double y;
  if(x[0] < p[1]){
    y = p[3];
  }else{
    //y = p[3] + p[2] * (TMath::Sin(xx) / (xx * xx) - TMath::Cos(xx)/xx); 
    y = p[3] + p[2] / pow(xx,p[4]) *
      ((15/(xx*xx*xx) - 6/xx) * TMath::Sin(xx) / (xx) - ((15/(xx*xx) -1) *TMath::Cos(xx)/xx) ); 
  }

  return y ;
}


int main( int argc, char* argv[] ) {

  if ( argc != 2 ){
    std::cerr<<"Usage: ptf_ttree_analysis.app ptf_analysis.root\n";
    exit(0);
  }

  //TFile * fin = new TFile( argv[1], "update" );
  TFile * fin = new TFile( argv[1], "read" );


  // get the waveform fit TTree
  std::cout << "TTree 0" << std::endl;
  TTree * tt0 = (TTree*)fin->Get("ptfanalysis0");
  WaveformFitResult * wf0 = new WaveformFitResult;
  if(tt0) wf0->SetBranchAddresses( tt0 );

  std::cout << "TTree 1" << std::endl;
  TTree * tt1 = (TTree*)fin->Get("ptfanalysis1");
  WaveformFitResult * wf1 = new WaveformFitResult;
  wf1->SetBranchAddresses( tt1 );

  std::cout << "TTree 17" << std::endl;  
  TTree * tt17 = (TTree*)fin->Get("ptfanalysis17");
  if(!tt17) tt17 = (TTree*)fin->Get("ptfanalysis17");
  WaveformFitResult * wf17 = new WaveformFitResult;
  if(tt17)  wf17->SetBranchAddresses( tt17 );

  TTree * tt18 = (TTree*)fin->Get("ptfanalysis18");
  WaveformFitResult * wf18 = new WaveformFitResult;
  wf18->SetBranchAddresses( tt18 );
  std::cout << "TTree done" << std::endl;

  TTree * tt19 = (TTree*)fin->Get("ptfanalysis19");
  WaveformFitResult * wf19 = new WaveformFitResult;
  if(tt19) wf19->SetBranchAddresses( tt19 );
  std::cout << "TTree done" << std::endl;

  TH1F *tdiff = new TH1F("time diff","Ch 0 minus Ch 1 time difference",100,-5,1);
  TH1F *tdiff0 = new TH1F("time diff0","PMT0 Time relative to Trigger Time",200,316,326);
  TH1F *tdiff1 = new TH1F("time diff1","PMT1 Time relative to Trigger Time",100,80,90);
  //  TH1F *tdiff1 = new TH1F("time diff1","PMT1 Time relative to Trigger Time",200,70,80);
  TH1F *tdiff2 = new TH1F("time diff2","timediff2",800,-6,1);
  TH1F *tdiff_inj = new TH1F("time diff inj","Time difference injected pulses",200,-9.3,-8.8);


  TH1F *tdiff00 = new TH1F("time diff00","PMT1 Time relative to Trigger Time",50,325,333);
  TH1F *tdiff01 = new TH1F("time diff01","PMT1 Time relative to Trigger Time",50,325,333);
  TH1F *tdiff02 = new TH1F("time diff02","PMT1 Time relative to Trigger Time",50,325,333);
  TH1F *tdiff03 = new TH1F("time diff03","PMT1 Time relative to Trigger Time",50,325,333);

  TH1F *tdiff_phase[8];
  for (int i = 0; i < 8; i++){
    char name1[100];
    sprintf(name1,"time diff phase %i",i);
    tdiff_phase[i]  = new TH1F(name1,"PMT1 Time relative to Trigger Time - phase",100,91.5,93.5);
  }

  TH1F *tdiff_phase_inj[8];
  for (int i = 0; i < 8; i++){
    char name1[100];
    sprintf(name1,"time diff phase %i inj",i);
    tdiff_phase_inj[i]  = new TH1F(name1,"Time difference injected pulses - by phase",200,-9.3,-8.8);
  }


  TH2F *tdiff_vs_ph = new TH2F("tdiff_vs_ph","Time difference vs pulse height",40,0,0.000488/0.018*80,50,70, 80);
  TProfile *tdiff_vs_ph_prof = new TProfile("tdiff_vs_ph_prof", "Time difference vs pulse height - Profile", 20,0,0.000488/0.018*80);
  
  TH1F *ph[2];
  ph[0] = new TH1F("PH0","Pulse Heights ",120,0,0.000488/0.018*120);
  ph[1] = new TH1F("PH1","Pulse Heights",2000,0,0.000488/0.018*2000);
  std::cout << "Looping tree " << tt0->GetEntries() << " " << tt1->GetEntries() << std::endl;
  int total_hits0 = 0, success_fits0 = 0;
  int total_hits1 = 0, success_fits1 = 0;

  int total_pe = 0.0;
  int total_nohits = 0.0;
  // Loop the first scan point and print something 
  for(int i = 0; i < tt0->GetEntries()-1; i++){
    //for(int i = 0; i < 50000; i++){
    tt0->GetEvent(i );
    tt1->GetEvent(i );
    //if(tt16) tt16->GetEvent(i );
    if(tt17) tt17->GetEvent(i );
    if(tt18) tt18->GetEvent(i );
    if(tt19) tt19->GetEvent(i );

    // Find the pulses in list of pulses
    double pulse_time[2] = {-1,-1};
    double pulse_height[2] {-1,-1};

    for(int j =0; j < 2; j++){
      WaveformFitResult *wf;
      if(j==0) wf = wf0;
      if(j==1) wf = wf1;
      for(int k = 0; k < wf->numPulses; k++){
        //if(wf->pulseTimes[k] > 2420 && wf->pulseTimes[k] < 2480){ // look for laser pulse
        if(wf->pulseTimes[k] > 2020 && wf->pulseTimes[k] < 2240){ // look for laser pulse
          pulse_time[j] = wf->pulseTimes[k];
          if(j == 1 && 0) std::cout << "Found " << i << " " << wf->pulseTimes[k] << " "
                               << (wf->pulseCharges[k])/0.018 << std::endl;

          //pulse_height[j] = (baseline[j] - wf->pulseCharges[k])/0.01;
          //pulse_height[j] = (wf->pulseCharges[k])/0.018 ;
          pulse_height[j] = (wf->pulseCharges[k])/0.018 ;
          //pulse_height[j] = (wf->pulseCharges[k])/0.016;
          //          std::cout << "Pulse Charge: " << wf->pulseCharges[k] << std::endl;
        }
      }
    }
    ph[0]->Fill(pulse_height[0]);
    ph[1]->Fill(pulse_height[1]);

    if(pulse_height[1] < 0.5) total_nohits += 1.0;
    if(pulse_height[1] > 0.5 && pulse_height[1] < 1.5) total_pe += 1.0;
    if(pulse_height[1] > 1.5 && pulse_height[1] < 2.5) total_pe += 2.0;
    if(pulse_height[1] > 2.5 && pulse_height[1] < 3.5) total_pe += 3.0;
    if(pulse_height[1] > 3.5 && pulse_height[1] < 4.5) total_pe += 4.0;
    if(pulse_height[1] > 4.5 && pulse_height[1] < 5.5) total_pe += 5.0;
    if(pulse_height[1] > 5.5 && pulse_height[1] < 6.5) total_pe += 6.0;
    if(pulse_height[1] > 6.5 && pulse_height[1] < 7.5) total_pe += 7.0;
    if(pulse_height[1] > 7.5 && pulse_height[1] < 8.5) total_pe += 8.0;
    if(pulse_height[1] > 8.5 && pulse_height[1] < 9.5) total_pe += 9.0;
    if(pulse_height[1] > 9.5 && pulse_height[1] < 10.5) total_pe += 10.0;
   
    double time0 = wf0->mean;
    double time1 = wf1->mean;
    double time16 = 0.0;//wf16->mean;
    double time17 = wf17->mean;
    double time18 = wf18->mean;
    double time19 = wf19->mean;
    double cfd_time0 = 0.0;//wf16->sinw;
    double cfd_time1 = 0.0;//wf18->sinw;
    double time_diff = 0.0;//time16 - time18;

    if(1 &&  i < 1000)std::cout << i << " Times: " << time18 << " "
              << time1 << " "
              << time1 - time18 << " "
              << time1 - time18 -75.6<< " "
              << pulse_height[0] << " "
              << pulse_height[1] << " "
              << cfd_time0 << " " 
              << cfd_time1 << " " 
              << std::endl;

    // Ch 0 to trigger
    //    if(pulse_height[0] < 1.5 && pulse_height[0] > 0.5){  tdiff0->Fill(time0-time16 + 1.5); }

    //    if(pulse_height[0] < 1.5 && pulse_height[0] > 1.25){  tdiff0->Fill(time0-time16); }
    if(pulse_height[0] <4.5 && pulse_height[0] > 3.5){
      double mtdiff = time0-time18;
      tdiff0->Fill(mtdiff);
      total_hits0++;
      if(mtdiff > 312 && mtdiff < 330) {success_fits0++;}
    }
    
    double td = time1-time18;


    if(pulse_height[1] < 0.75 && pulse_height[1] > 0.5){  tdiff00->Fill(time1-time18); }
    if(pulse_height[1] < 1.0 && pulse_height[1] > 0.75){  tdiff01->Fill(time1-time18); }
    if(pulse_height[1] < 1.25 && pulse_height[1] > 1.0){  tdiff02->Fill(time1-time18); }
    if(pulse_height[1] < 1.5 && pulse_height[1] > 1.25){  tdiff03->Fill(time1-time18); }

    tdiff_vs_ph->Fill(pulse_height[1],td);
    if(td > 327.5 && td < 330){
      tdiff_vs_ph_prof->Fill(pulse_height[1],td);
    }// Ch 1 to trigger

    int phase = ((int)(time18-0.5)) % 8 ;
    if(pulse_height[1] < 32.5 && pulse_height[1] > 28.5){

      double extra = 0.0;
      if(phase == 0) extra = 0.16;
      if(phase == 1) extra = 0.14;
      if(phase == 2) extra = 0.09;
      if(phase == 3) extra = 0.11;
      if(phase == 4) extra = 0.05;
      if(phase == 5) extra = 0.00;
      if(phase == 6) extra = 0.07;
      if(phase == 7) extra = 0.18;

      double mtdiff = time1-time18;
      //double mtdiff = time1-time18 - extra;

      //      if(phase == 2 || phase==3)
      //if(phase != 6)
        tdiff1->Fill(mtdiff);

      tdiff_phase[phase]->Fill(mtdiff);

      total_hits1++;
      if(mtdiff > 60 && mtdiff < 340) {success_fits1++;}
    

    }

    // check injected pulse timing
    double tdiff_i = time1-time17;
    //std::cout << "Injected: " << tdiff_i << std::endl;
    tdiff_inj->Fill(tdiff_i);
    tdiff_phase_inj[phase]->Fill(tdiff_i);
      
    // ch 0 - ch 1
    if(pulse_height[0] < 200.5 && pulse_height[1] < 200.5 &&
       pulse_height[0] > -2 && pulse_height[1] > -2       ){
      tdiff->Fill(time16-time18);
      tdiff2->Fill(cfd_time0-cfd_time1);

      if(i < 1000 && 1 && (time16-time18) > -1.25 &&
         (time16-time18) < -1.15 
         ){
        std::cout << i << " Time difference: " << time16 << " - " << time18
                  << " : " << time16-time18 << std::endl;
        std::cout << "Charge: " << wf0->amp << " " << wf1->amp << " "
                  << pulse_time[0] << " " << pulse_time[1] << " "
                  << pulse_height[0] << " " << pulse_height[1] << " "
                  << std::endl;
        std::cout << "Ratio "
                  << wf0->amp/pulse_height[0] << " " 
                  << wf1->amp/pulse_height[1] << " " << std::endl;
        
      }
    }     
  } 

  std::cout << "Successful fit (chan 0) = " << success_fits0 << " / " << total_hits0
            << " : " << (((double) success_fits0) /((double)total_hits0) * 100.0) << "%" << std::endl;
  std::cout << "Successful fit (chan 1) = " << success_fits1 << " / " << total_hits1
            << " : " << (((double) success_fits1) /((double)total_hits1) * 100.0) << "%" << std::endl;

  std::cout << "Total events " << tt0->GetEntries() << " total PE " << total_pe << std::endl;
  std::cout << "Mean: " << total_pe/(double)tt0->GetEntries() << std::endl;
  std::cout << "P(0) meas: " << total_nohits/(double)tt0->GetEntries() << std::endl;
  double lambda = total_pe/(double)tt0->GetEntries();
  double calc_p0 = exp(-lambda);// https://en.wikipedia.org/wiki/Poisson_distribution#Poisson_Approximation
  std::cout << "P(0) calc: " << calc_p0 << std::endl;
  double better_lambda = - log(total_nohits/(double)tt0->GetEntries());
  std::cout << "Better lamba: " << better_lambda << std::endl;

    
    TCanvas *c = new TCanvas("C");

  tdiff->Draw();
  //TF1 *gaus = new TF1("gaus","gaus",-5,2);

  //tdiff->Fit("gaus","R");
  tdiff->SetLineColor(2);
  //tdiff2->Draw("SAME");
  //tdiff->Fit("gaus","R");
  tdiff->SetXTitle("time difference (ns)"); 

  gStyle->SetOptFit(1111);


  
  c->SaveAs("tdiff.png");

  TCanvas *c1 = new TCanvas("C1");
  tdiff0->Draw();
  //TF1 *gaus = new TF1("gaus","gaus",321.5,324);
  TF1 *gaus1 = new TF1("gaus1","gaus",319.5,321.5);
  tdiff0->Fit("gaus1","R");
  tdiff0->SetXTitle("time difference (ns)"); 

  gStyle->SetOptFit(1111);

  double mean0 = gaus1->GetParameter(0);
  double hm0 = mean0/2.0;
  double start0 = -1, end0 = -1;
  for (int i = 1; i < tdiff0->GetNbinsX(); i++){
    double tmp = tdiff0->GetBinContent(i);
    if(tmp > hm0 && start0 < 0) start0 = tdiff0->GetBinCenter(i);
    if(tmp < hm0 && start0 > 0 && end0 < 0) end0 = tdiff0->GetBinCenter(i);
    
  }
  double fwhm0 = end0-start0;
  std::cout << "mean0 : " << mean0 << " "
            << start0 << " "
            << end0 << " "
            << fwhm0 
            << std::endl;


  
  c1->SaveAs("tdiff_injected0.png");
  
  TCanvas *c2 = new TCanvas("C2");
  tdiff1->Draw();
  //TF1 *gaus = new TF1("gaus","gaus",321.5,324);
    TF1 *gaus2 = new TF1("gaus2","gaus",85.5,86.5);
  //  TF1 *gaus2 = new TF1("gaus2","gaus",74,75.6);
  tdiff1->Fit("gaus2","R");
  tdiff1->SetXTitle("time difference (ns)"); 

  gStyle->SetOptFit(1111);

  double mean1 = gaus2->GetParameter(0);
  double hm1 = mean1/2.0;
  double start1 = -1, end1 = -1;
  for (int i = 1; i < tdiff1->GetNbinsX(); i++){
    double tmp = tdiff1->GetBinContent(i);
    if(tmp > hm1 && start1 < 0) start1 = tdiff1->GetBinCenter(i);
    if(tmp < hm1 && start1 > 0 && end1 < 0) end1 = tdiff1->GetBinCenter(i);
    
  }
  double fwhm = end1-start1;
  std::cout << "mean1 : " << mean1 << " " << hm1 << " " 
            << start1 << " "
            << end1 << " "
            << fwhm
            << " " << tdiff1->GetRMS() << " " 
            << std::endl;
  

  c2->SaveAs("tdiff_injected1.png");



  TCanvas *c2inj = new TCanvas("C2inj");
  tdiff_inj->Draw();
  //TF1 *gaus = new TF1("gaus","gaus",321.5,324);
  TF1 *gaus3 = new TF1("gaus3","gaus",-9.1,-9.05);
  tdiff_inj->Fit("gaus3","R");
  tdiff_inj->SetXTitle("time difference (ns)"); 

  gStyle->SetOptFit(1111);

  double mean3 = gaus3->GetParameter(0);
  double hm3 = mean3/2.0;
  double start3 = -9999, end3 = -9999;
  for (int i = 1; i < tdiff_inj->GetNbinsX(); i++){
    double tmp = tdiff_inj->GetBinContent(i);
    if(tmp > hm3 && start3 < -1000) start3 = tdiff_inj->GetBinCenter(i);
    if(tmp < hm3 && start3 > -1000 && end3 < -1000) end3 = tdiff_inj->GetBinCenter(i);
    
  }
  double fwhm3 = end3-start3;
  std::cout << "mean (inj) : " << mean3 << " " << hm3 << " " 
            << start3 << " "
            << end3 << " "
            << fwhm3
            << " " << tdiff_inj->GetRMS() << " " 
            << std::endl;
  

  c2inj->SaveAs("tdiff_injected_injected.png");


  TCanvas *c3 = new TCanvas("C3");
  ph[1]->Draw();
  ph[1]->SetXTitle("Pulse height (PE)");

  //  ph[1]->Draw("SAME");
  //ph[1]->SetLineColor(2);

  TLegend *leg = new TLegend(0.5,0.7,0.79,0.89);
  leg->AddEntry(ph[0],"Channel 0");
  leg->AddEntry(ph[1],"Channel 1");
  leg->Draw("SAME");    
  
  c3->SaveAs("pulse_heights.png");


  TCanvas *c4 = new TCanvas("C4");

  char name[100];
  sprintf(name,"PMT1_NoWaveforms/hwf_%i;1",103);
  TH1 *ch0 = (TH1*) fin->Get(name);

  if(ch0){
    ch0->Draw();
    ch0->GetXaxis()->SetRangeUser(2000,2450);


    TF1 *ffitfunc1 = new TF1("mygauss1",bessel,2100,2300,5);
    

    
    ffitfunc1->SetParameter(0, 0.12);
    ffitfunc1->SetParameter(1, 2200 );
    ffitfunc1->SetParameter(2, -0.3 );
    ffitfunc1->FixParameter(3, 0.996);
    ffitfunc1->SetParameter(4, 1.4 );

       
    ch0->Fit( ffitfunc1, "", "", 2170,2266);

    
    
  }
  c4->SaveAs("waveform_example.png");


  TCanvas *c5 = new TCanvas("C5","Timing Diff",1000,800);
  tdiff03->Draw("HIST");
  tdiff03->SetLineColor(4);
  tdiff03->SetXTitle("Timing difference (ns)");
  
  tdiff00->Draw("SAMEHIST");
  tdiff01->Draw("SAMEHIST");
  tdiff01->SetLineColor(2);
 
  tdiff02->Draw("SAMEHIST");
  tdiff02->SetLineColor(3);

  tdiff01->Scale(tdiff00->GetEntries()/tdiff01->GetEntries());
  tdiff02->Scale(tdiff00->GetEntries()/tdiff02->GetEntries());
  tdiff03->Scale(tdiff00->GetEntries()/tdiff03->GetEntries());

  TLegend *leg4 = new TLegend(0.6,0.6,0.89,0.89);
  leg4->AddEntry(tdiff00,"0.5PE - 0.75PE");
  leg4->AddEntry(tdiff01,"0.75PE - 1.0PE");
  leg4->AddEntry(tdiff02,"1.0PE - 1.25PE");
  leg4->AddEntry(tdiff03,"1.25PE - 1.5PE");

  leg4->Draw("SAME");
  
  c5->SaveAs("tdiff_injected0_many.png");
  
  TCanvas *c6 = new TCanvas("C6","Timing Diff vs PH",1000,800);
  tdiff_vs_ph->Draw("COLZ");
  tdiff_vs_ph->SetXTitle("Pulse Height (PE)");
  tdiff_vs_ph->SetYTitle("Time Difference (ns)");
  //tdiff_vs_ph_prof->Draw("SAME");

  
  c6->SaveAs("tdiff_vs_ph_1pe.png");


    
  TCanvas *c7 = new TCanvas("C7","Timing Diff vs PH",1000,800);
  tdiff_vs_ph_prof->Draw();
  tdiff_vs_ph_prof->SetXTitle("Pulse Height (PE)");
  tdiff_vs_ph_prof->SetYTitle("Time Difference (ns)");
  tdiff_vs_ph_prof->GetYaxis()->SetRangeUser(328.6,328.9);

  

  c7->SaveAs("tdiff_vs_ph_1pe_prof.png");

  TCanvas *c8 = new TCanvas("C7","Timing Diff vs PH",1000,800);

  TLegend *legp = new TLegend(0.6,0.6,0.89,0.89);
  for(int i = 0; i < 8; i++){

    if(i == 0){
      tdiff_phase[0]->Draw();
      tdiff_phase[0]->GetYaxis()->SetRangeUser(0,tdiff_phase[0]->GetMaximum()*1.2);
    }else{
      tdiff_phase[i]->Scale((float)tdiff_phase[0]->GetEntries()/((float)tdiff_phase[i]->GetEntries()));
      tdiff_phase[i]->Draw("SAMEHIST");
      tdiff_phase[i]->SetLineColor(i+1);
    }
    char nnn[100];
    sprintf(nnn,"Phase %i",i);
    legp->AddEntry(tdiff_phase[i],nnn);

    std::cout << i << " " << tdiff_phase[i]->GetEntries()
              << " " << tdiff_phase[i]->GetMean()
              << " " << tdiff_phase[i]->GetRMS()
              << std::endl;
    
  }
  legp->Draw("SAME");

  c8->SaveAs("tdiff_vs_phase.png");


  TCanvas *c9 = new TCanvas("C7","Timing Diff vs PH",1000,800);

  TLegend *legp2 = new TLegend(0.6,0.6,0.89,0.89);
  for(int i = 0; i < 8; i++){

    if(i == 0){
      tdiff_phase_inj[0]->Draw();
      tdiff_phase_inj[0]->GetYaxis()->SetRangeUser(0,tdiff_phase_inj[0]->GetMaximum()*1.2);
    }else{
      tdiff_phase_inj[i]->Scale((float)tdiff_phase_inj[0]->GetEntries()/((float)tdiff_phase_inj[i]->GetEntries()));
      tdiff_phase_inj[i]->Draw("SAMEHIST");
      tdiff_phase_inj[i]->SetLineColor(i+1);
    }
    char nnn[100];
    sprintf(nnn,"Phase %i",i);
    legp2->AddEntry(tdiff_phase_inj[i],nnn);

    std::cout << i << " " << tdiff_phase[i]->GetEntries()
              << " " << tdiff_phase[i]->GetMean()
              << " " << tdiff_phase[i]->GetRMS()
              << std::endl;
    
  }
  legp2->Draw("SAME");

  c9->SaveAs("tdiff_vs_phase_inj.png");
  
  tdiff0->Write();
  fin->Close();
  return 0;
}

