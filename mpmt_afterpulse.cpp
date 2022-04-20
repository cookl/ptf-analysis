#include "WaveformFitResult.hpp"
#include "ScanPoint.hpp"
#include "TFile.h"
#include "TStyle.h"
#include "TCanvas.h"
#include "TH1.h"
#include "TH2.h"
#include "TF1.h"
#include "TPad.h"
#include "TLegend.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TText.h"
#include <iostream>
#include <vector>
#include <map>
#include <fstream>
#include <math.h>
#include <string>

//// This programs opens the ROOT file
//   containing the waveforms' information,
//   such as number, time position and 
//   height of pulses.
//   It outputs afterpulse rates and graphs, 
//   and histograms with pulse height and time 
//   distributions.

void analysis(TTree* tt, WaveformFitResult* wf, const Double_t peHeight, const Double_t peDivision, const Double_t &afpTimeThreshold1, Double_t &afpTimeThreshold2, Double_t dcr, std::vector<Double_t> &vecT, std::vector<Double_t> &vecQ, std::vector<Double_t> &vecA, std::vector<Double_t> &vecDeltaT, std::vector<Double_t> &vecAlsr, std::vector<Double_t> &vecQlsr, std::vector<Double_t> &vecAafp, std::vector<Double_t> &vecQafp, std::map<Int_t, Double_t> &lsrCounts, std::map<Int_t, Double_t> &afpRate, std::map<Int_t, Double_t> &singleAfpRate, std::map<Int_t, Double_t> &multipleAfp_dcr, TH1F* histPE);
void rates(std::map<Int_t, Double_t> &lsrCounts, std::map<Int_t, Double_t> &singleCounts, std::map<Int_t, Double_t> &multipleCounts, std::map<Int_t, Double_t> afpCounter_corrected, std::map<Int_t, Double_t> &singleRates, std::map<Int_t, Double_t> &multipleRates, std::map<Int_t, Double_t> &multipleRates_corrected, std::map<Int_t, Double_t> &singleErrors, std::map<Int_t, Double_t> &multipleErrors, std::map<Int_t, Double_t> &multipleRates_corrected_errors);
void dcr_estimate(TTree* tt, WaveformFitResult* wf, const Double_t early_threshold, const Double_t afpTimeThreshold1, Double_t &dcr_estimation, Double_t &dcr_error);
void drawHistogram(std::vector<double> data, TH1F* histogram, const char runNumber[], const Int_t y1_in, const Int_t y2_in, const Double_t x1_in, const Double_t x2_in, const Double_t percentage_x, const Double_t percentage_y);
// void draw2DHistogram();

int main( int argc, char* argv[] ) {

  //// PRE-RUN TEST
  //   Only run analysis if all arguments are provided.
  if ( argc != 5 ){
    std::cerr<<"Usage: ptf_ttree_analysis.app ptf_analysis.root run_number ptfanalysis0 channel_number\n";
    exit(0);
  }
  ////
  
  // Arbritary values of a time window:
  double afpTimeThreshold1 = 2180; //< ns
  double afpTimeThreshold2 = 2300; //< ns

  // P.E. values
  double peDivision = 0.5;
  double peHeight   = 9.261; //< mV 
  double peArea     = 260.0; // mV.ns ?

  // DCR parameters
  Double_t early_threshold = 20.;//<ns
  
  // Open ROOT file
  TFile * fin = new TFile( argv[1], "read" );
  
  // Get the waveform fit TTree
  TTree * tt = (TTree*)fin->Get(argv[3]);
  WaveformFitResult * wf = new WaveformFitResult;
  wf->SetBranchAddresses( tt );

  // histogram for pulse heights
  TH1F * histPE = new TH1F("pulse_height","Pulse Height",200,0,200*0.48828125);
  
  //// Vectors for storing pulse and afterpulse information
  std::vector<Double_t> vecT, vecQ, vecA, vecDeltaT, vecAlsr, vecQlsr, vecAafp, vecQafp;
  std::map<Int_t, Double_t> lsrCounts, afpCounts, singleAfpCounts, afpCounter_corrected;
  std::map<Int_t, Double_t> singleRates, multipleRates, multipleRates_corrected, singleErrors, multipleErrors, multipleErrors_corrected, multipleRates_corrected_error;

  // Estimate and apply dark count rate (DCR)
  Double_t dcr_estimation, dcr_estimation_error, dark_pulses_per_waveform;
  dcr_estimate(tt, wf, early_threshold, afpTimeThreshold1, dcr_estimation, dcr_estimation_error);

  // Find pulses, count and classify them according to p.e. level.
  analysis(tt, wf, peHeight, peDivision, afpTimeThreshold1, afpTimeThreshold2, dcr_estimation, vecT, vecQ, vecA, vecDeltaT, vecAlsr, vecQlsr, vecAafp, vecQafp, lsrCounts, afpCounts, singleAfpCounts, afpCounter_corrected, histPE);

  // Measure afterpulse rate according to p.e. level.
  rates(lsrCounts,singleAfpCounts,afpCounts,afpCounter_corrected,singleRates,multipleRates,multipleRates_corrected,singleErrors,multipleErrors,multipleErrors_corrected);
  
    // // Debugging
  std::map<Int_t, Double_t>::iterator it_lsr;
  for (it_lsr = lsrCounts.begin(); it_lsr != lsrCounts.end(); it_lsr++) {
    std::cout<<it_lsr->first<<"  "<<multipleRates_corrected.find(it_lsr->first)->second<<"  "<<multipleErrors_corrected.find(it_lsr->first)->second<<std::endl;
  }

  //
  std::cout<<"DCR:       "<<dcr_estimation<<" Hz."<<std::endl;
  std::cout<<"DCR error: "<<dcr_estimation_error<<" Hz."<<std::endl;
  std::cout<<"Dark pulses:       "<<dcr_estimation*1024*8*1.e-9<<" per waveform."<<std::endl;
  std::cout<<"Dark pulses error: "<<dcr_estimation_error*1024*8*1.e-9<<" per waveform."<<std::endl;

  ////////////////////////////////////////////////////
  // DRAWING STUFF
  
  // Useful variables:
  Int_t ymin, ymax, dirInt;
  Double_t xmin, xmax, titlex, titley;
  char filename[1024];

  auto canvas3 = new TCanvas("canvas3","",1200,800);

  // Change margins (left, right, bottom, top)
  canvas3->SetMargin(0.14, 0.03, 0.09, 0.04);
        
  
  // // HISTOGRAMS
  //    For the histograms, it was easy enough to split it
  //    in a separate function.
  canvas3->cd(0);
  TH1F* histTimeDistance = new TH1F("histDeltaTime","", 155, -100., 6100.);
  drawHistogram(vecDeltaT, histTimeDistance, argv[4], ymin = 0, ymax = 25000, xmin = 0., xmax = 0., titlex = 70., titley = 80.);
  sprintf(filename, "../images/%s-%s-hist-DeltaTime.png", argv[2], argv[3] );
  canvas3->SaveAs(filename);



  // // AFTERPULSE GRAPHS
  //
  canvas3->cd(1);
  TGraphErrors* graphSingle = new TGraphErrors();
  TGraphErrors* graphMultiple = new TGraphErrors();
  TGraphErrors* graphCorrected = new TGraphErrors();
 
  // Find data
  std::map<Int_t,Double_t>::iterator it;
  graphSingle->Set(singleRates.size()-1);  
  Int_t i = 0;
  for (it = singleRates.begin(); it != singleRates.end(); ++it){
    ++i;
    graphSingle->SetPoint(i-1, it->first*0.5, it->second*100.);
    graphSingle->SetPointError(i-1, 0., singleErrors.find(it->first)->second*100.);
    graphMultiple->SetPoint(i-1, it->first*peDivision, multipleRates.find(it->first)->second*100.);
    graphMultiple->SetPointError(i-1, 0., multipleErrors.find(it->first)->second*100.);
    graphCorrected->SetPoint(i-1, it->first*peDivision, multipleRates_corrected.find(it->first)->second*100.);
    graphCorrected->SetPointError(i-1, 0., multipleErrors_corrected.find(it->first)->second*100.);
    std::cout << "Multi error is " << multipleErrors.find(it->first)->second*100. << std::endl;
    std::cout << "Multi corrected error is " << multipleErrors_corrected.find(it->first)->second*100. << std::endl;
  }
  
  // Styles
  graphSingle->GetXaxis()->SetTitle("Laser pulse (p.e.)");
  graphSingle->GetYaxis()->SetTitle("Percentage (\%)");
  graphSingle->SetMarkerStyle(8);
  graphSingle->SetMarkerSize(1.4);
  graphSingle->SetMarkerColor(8);
  graphSingle->SetLineColor(1);
  graphSingle->GetYaxis()->SetRangeUser(0, 100);
  graphSingle->GetXaxis()->SetRangeUser(0, 15);
  graphMultiple->SetMarkerColor(2);
  graphMultiple->SetMarkerStyle(8);
  graphMultiple->SetMarkerSize(1.4);
  graphCorrected->SetMarkerColor(4);
  graphCorrected->SetMarkerStyle(8);
  graphCorrected->SetMarkerSize(1.4);

  // Legend
  TLegend *leg2 = new TLegend(0.16,0.7,0.6,0.85);
  leg2->AddEntry(graphMultiple,  "All events");
  leg2->AddEntry(graphCorrected, "All events (corrected)");
  leg2->AddEntry(graphSingle,    "Single event");
  leg2->SetTextSize(0.036);
  leg2->SetLineWidth(0);

  // Draw
  graphSingle->Draw("AP");
  graphMultiple->Draw("P");
  graphCorrected->Draw("P");
  leg2->Draw("SAME");

  // Multiple rate fitting
  std::cout<<"\nMultiple events:"<<std::endl;
  TF1* multipleFitting = new TF1("linear_multiple","[0] + [1]*x",0.0,multipleRates.end()->first);
  multipleFitting->SetParameters(1., 1., 1.);
  multipleFitting->SetLineColor(46);
  multipleFitting->SetLineStyle(2);
  graphMultiple->Fit(multipleFitting, "", "", 0.0, multipleRates.end()->first);
  Double_t intercept = multipleFitting->GetParameter(0);
  Double_t slope     = multipleFitting->GetParameter(1);

  // Multiple rate fitting
  std::cout<<"\nSingle events:"<<std::endl;
  TF1* singleFitting = new TF1("poly_single","[0] + [1]*x^[2]", 0.0, singleRates.end()->first);
  singleFitting->SetParameters(1.,1.);
  singleFitting->SetLineColor(30);
  singleFitting->SetLineStyle(2);
  graphSingle->Fit(singleFitting, "", "", 0.0, singleRates.end()->first);

  // Corrected rate fitting
  std::cout<<"\nMultiple events corrected for DCR:"<<std::endl;
  TF1* correctedFitting = new TF1("linear_corrected","[0] + [1]*x", multipleRates_corrected.begin()->first*peDivision, multipleRates_corrected.end()->first*peDivision);
  correctedFitting->SetParameters(intercept-dcr_estimation*1024*1.e-9*100, slope);
  std::cout<<"Estimated intercept: "<<intercept-dcr_estimation*1024*1.e-9*100<<std::endl;
  correctedFitting->SetParameters(0.0, 14.);
  correctedFitting->SetLineColor(38);
  correctedFitting->SetLineStyle(2);
  graphCorrected->Fit(correctedFitting, "", "", 0., 16);

  // Add a title to the graph
  char inGraphTitle[1024];
  sprintf(inGraphTitle, "Channel %s", argv[4] );
  auto *th1 = new TText(2,250,inGraphTitle);
  th1->SetTextAlign(11); th1->SetTextSize(0.040);
  th1->Draw();

  // Save file
  sprintf(filename, "../images/%s-%s-graph-afp-all.png", argv[2], argv[3] );
  canvas3->SaveAs(filename);


  // // 2D HISTOGRAMS
  // 

  // These are color properties sent to me by Patrick.
  gStyle->SetPalette(1,0); 
  const Int_t NRGBs = 5;
  const Int_t NCont = 255;
  Double_t stops[NRGBs] = { 0.00, 0.34, 0.61, 0.84, 1.00 }; // < Original
  Double_t red[NRGBs]   = { 0.00, 0.00, 0.87, 1.00, 0.51 };
  Double_t green[NRGBs] = { 0.00, 0.81, 1.00, 0.20, 0.00 };
  Double_t blue[NRGBs]  = { 0.51, 1.00, 0.12, 0.00, 0.00 };
  TColor::CreateGradientColorTable(NRGBs, stops, red, green, blue, NCont);
  gStyle->SetNumberContours(NCont);

  // // 2D HISTOGRAM
  //
  canvas3->cd(2);
  canvas3->SetMargin(0.10, 0.11, 0.09, 0.04);
  auto hist2D = new TH2F("hist2D","",500,1000,0,8300*2.0*0.48828125, 0,1000*0.08*0.48828125);
  for (int i=0; i<vecQ.size(); i++) {
    hist2D->Fill(vecT.at(i),vecQ.at(i));
  }
  
  // Set style
  hist2D->GetXaxis()->SetTitle("Time (ns)");
  hist2D->GetYaxis()->SetTitle("Pulse height (p.e.)");
  hist2D->SetFillColor(1);
  // hist2D->Draw("box");
  hist2D->SetStats(0);
  hist2D->GetXaxis()->SetRangeUser(1000.,8000.);
  // hist2D->GetYaxis()->SetRangeUser(0,20.);
  gPad->SetLogz(); // </ Change COLZ to a logarithmic scale

  hist2D->Draw("COLZ");
  
  // Save file
  sprintf(filename, "../images/%s-%s-hist2D-all.png", argv[2], argv[3] );
  canvas3->SaveAs(filename);

  //Stuff for the pulse height histogram
  TCanvas *cPE = new TCanvas("CPE");
  histPE->Draw();
  histPE->GetXaxis()->SetTitle("Pulse height (PE)");
  histPE->GetXaxis()->SetRangeUser(0, 20);
  histPE->GetYaxis()->SetTitle("Number of events");
  histPE->Fit("gaus","","",3,6.5);
  gStyle->SetOptFit(11);
  sprintf(filename, "../images/%s-%s-mpmt-pulse-height.png", argv[2], argv[3] );
  cPE->SaveAs(filename);
 
  return 0;
}

void analysis(TTree* tt, WaveformFitResult* wf, const Double_t peHeight, const Double_t peDivision, const Double_t &afpTimeThreshold1, Double_t &afpTimeThreshold2, Double_t dcr, std::vector<Double_t> &vecT, std::vector<Double_t> &vecQ, std::vector<Double_t> &vecA, std::vector<Double_t> &vecDeltaT, std::vector<Double_t> &vecAlsr, std::vector<Double_t> &vecQlsr, std::vector<Double_t> &vecAafp, std::vector<Double_t> &vecQafp, std::map<Int_t, Double_t> &lsrCounts, std::map<Int_t, Double_t> &afpRate, std::map<Int_t, Double_t> &singleAfpRate, std::map<Int_t, Double_t> &multipleAfp_dcr, TH1F* histPE) {

  // Defining variables
  bool afpCounter;
  Int_t lsrPulseInt;
  float_t lsrPulseFlt;

  for(int k = 0; k < tt->GetEntries(); k++){
      
    // Get event information from tree.
    tt->GetEvent(k);
    
    // If there are pulses on waveform
    if(wf->numPulses > 0){
      // std::cout << "." ; 
      // Start afterpulse counter at zero for this waveform.
      Int_t singleAfpCounter = 0, lsrCounter = 0;
      Double_t afpCounter = 0, afpCounter_corrected = 0;

      // Count number of pulses on the waveform,
      // considering only the waveforms whose first pulses 
      // are within the laser window.

    // if (wf->pulseTimes[0]<afpTimeThreshold2) { // < To consider all types of waveforms.
	  if ((wf->pulseTimes[0]>=afpTimeThreshold1) && (wf->pulseTimes[0]<afpTimeThreshold2)) { //< To consider only waveform with laser pulse as first signal.
        
        // Count the laser and find its height
        lsrCounter++;
	lsrPulseFlt = wf->pulseCharges[0]*1000.0/peHeight/peDivision;
	lsrPulseInt = round(lsrPulseFlt);
	
	// std::cout << "Laser=" << lsrPulseInt << " " << wf->numPulses ;
        // Count the number of afterpulses, if any
        afpCounter_corrected -= (dcr)*(1024*8)*1.e-9;
        if (wf->numPulses > 1) {
          singleAfpCounter++;
          afpCounter += Double_t(wf->numPulses)-1.;
          afpCounter_corrected += Double_t(wf->numPulses)-1.;
        }
	// std::cout << std::endl;
      }
      

      // Add laser-count information to map
      std::map<Int_t, Double_t>::iterator it_lsr = lsrCounts.find(lsrPulseInt);
      if (it_lsr != lsrCounts.end()) {
        it_lsr->second += (Double_t)lsrCounter;
      } else {
        lsrCounts.insert(std::make_pair(lsrPulseInt, (Double_t)lsrCounter));
      }

      // Add multiple-afterpulse information to a map
      std::map<Int_t, Double_t>::iterator it_multiple = afpRate.find(lsrPulseInt);
      if (it_multiple != afpRate.end()) {
        it_multiple->second += (Double_t)afpCounter;
      } else {
        afpRate.insert(std::make_pair(lsrPulseInt, (Double_t)afpCounter));
      }

      // Add single-afterpulse information to a map
      std::map<Int_t, Double_t>::iterator it_single   = singleAfpRate.find(lsrPulseInt);
      if (it_single != singleAfpRate.end()) {
        it_single->second += (Double_t)singleAfpCounter;
      } else {
        singleAfpRate.insert(std::make_pair(lsrPulseInt, (Double_t)singleAfpCounter));
      }

      // Adding pair or incrementing respective pair for DCR corrected waveforms
      std::map<Int_t, Double_t>::iterator it_corrected = multipleAfp_dcr.find(lsrPulseInt);
      if (it_corrected != multipleAfp_dcr.end()) {
        it_corrected->second += (Double_t)afpCounter_corrected;
      } else {
        multipleAfp_dcr.insert(std::make_pair(lsrPulseInt, (Double_t)afpCounter_corrected));
      }


      // Data for the histograms.
      for(int i = 0; i < wf->numPulses; i++){

        // Overall pulse distributions
        vecT.push_back(wf->pulseTimes[i]);
        vecQ.push_back(wf->pulseCharges[i]*1000./peHeight);
        vecA.push_back(wf->pulseArea);  

        // Collect pulse time with laser time at zero.
        if (wf->pulseTimes[i]-afpTimeThreshold1 > 0) {
          vecDeltaT.push_back(wf->pulseTimes[i]-afpTimeThreshold1);
        }
        
        // check if it is laser pulse. 
        if (wf->pulseTimes[i]>=afpTimeThreshold1 && wf->pulseTimes[i]<afpTimeThreshold2){
          vecAlsr.push_back(wf->pulseArea);
          vecQlsr.push_back(wf->pulseCharges[i]*1000.0/peHeight);
	  histPE->Fill(wf->pulseCharges[i]*1000.0/peHeight);
        
        // otherwise it is afterpulse
        } else if (wf->pulseTimes[i]>=afpTimeThreshold2){
          vecAafp.push_back(wf->pulseArea);
          vecQafp.push_back(wf->pulseCharges[i]*1000.0/peHeight);
        }
      }
    }             
  }
  // // Debugging
  // std::map<Int_t, Double_t>::iterator it_lsr = lsrCounts.find(lsrPulseInt);
  // for (it_lsr = lsrCounts.begin(); it_lsr != lsrCounts.end(); it_lsr++) {
  //   std::cout<<it_lsr->first<<" "<<multipleAfp_dcr.find(it_lsr->first)->second<<std::endl;
  // }
}

void rates(std::map<Int_t, Double_t> &lsrCounts, std::map<Int_t, Double_t> &singleCounts, std::map<Int_t, Double_t> &multipleCounts, std::map<Int_t, Double_t> afpCounter_corrected, std::map<Int_t, Double_t> &singleRates, std::map<Int_t, Double_t> &multipleRates, std::map<Int_t, Double_t> &multipleRates_corrected, std::map<Int_t, Double_t> &singleErrors, std::map<Int_t, Double_t> &multipleErrors, std::map<Int_t, Double_t> &multipleRates_corrected_errors) {
  
  // Finding the "lsrInt" iterator from lsrCounts
  std::map<Int_t, Double_t>::iterator it_lsr;

  // Calculate the rates, 
  // mapping their value to the lsrInt index.
  for (it_lsr = lsrCounts.begin(); it_lsr != lsrCounts.end(); ++it_lsr) {
    singleRates.insert(std::make_pair(it_lsr->first, (singleCounts.find(it_lsr->first)->second)/(it_lsr->second) ));
    multipleRates.insert(std::make_pair(it_lsr->first, (multipleCounts.find(it_lsr->first)->second)/(it_lsr->second) ));
    multipleRates_corrected.insert(std::make_pair(it_lsr->first, (afpCounter_corrected.find(it_lsr->first)->second)/(it_lsr->second) ));
    // std::cout<<singleCounts.find(it_lsr->first)->second<<"  "<<it_lsr->second<<"  "<<singleRates.find(it_lsr->first)->second <<std::endl;
    // std::cout<<it_lsr->first*0.5<<"  "<<100.*singleCounts.find(it_lsr->first)->second/(it_lsr->second)<<std::endl;
  }

  // Doing the same for the errors.
  for (it_lsr = lsrCounts.begin(); it_lsr != lsrCounts.end(); ++it_lsr) {

    // Errors in the count of pulses
    Double_t lsrCount_error      = sqrt(it_lsr->second);
    Double_t singleCount_error   = sqrt(singleCounts.find(it_lsr->first)->second);
    Double_t multipleCount_error = sqrt(multipleCounts.find(it_lsr->first)->second);
    Double_t afpCounter_corr_err = sqrt(afpCounter_corrected.find(it_lsr->first)->second);

    // The compound error of rates, in quadrature
    Double_t singleError    = singleRates.find(it_lsr->first)->second*sqrt(pow(singleCount_error/singleCounts.find(it_lsr->first)->second, 2) + pow(lsrCount_error/it_lsr->second, 2)); 
    Double_t multipleError  = multipleRates.find(it_lsr->first)->second*sqrt(pow(multipleCount_error/multipleCounts.find(it_lsr->first)->second, 2) + pow(lsrCount_error/it_lsr->second, 2)); 
    Double_t correctedError = multipleRates_corrected.find(it_lsr->first)->second*sqrt(pow(afpCounter_corr_err/afpCounter_corrected.find(it_lsr->first)->second, 2) + pow(lsrCount_error/it_lsr->second, 2) );

    // Add the errors to their maps
    singleErrors.insert(std::make_pair(it_lsr->first, singleError));
    multipleErrors.insert(std::make_pair(it_lsr->first, multipleError));
    multipleRates_corrected_errors.insert(std::make_pair(it_lsr->first, correctedError));

  }
}

void dcr_estimate(TTree* tt, WaveformFitResult* wf, const Double_t early_threshold, const Double_t afpTimeThreshold1, Double_t &dcr_estimation, Double_t &dcr_error){

  // Estimate the dark-count rate.
  // Count all dark pulses, 
  // the pulses in between early_threshold 
  // and afpTimeThreshold1,
  // divide it by a time window and 
  // find average DCR.

  Double_t timeWindow = afpTimeThreshold1*1.e-9;
  std::vector<Double_t> waveformDCRs, waveformDCRs_errors;

  std::cout << "Total number of waveforms is " << tt->GetEntries() << std::endl;

  // Analyze all waveforms and find the ones with dark pulses
  for (int k = 0; k < tt->GetEntries(); k++){
  
    // Find event
    tt->GetEvent(k);

    if (wf->numPulses > 0) {
    
      Int_t darkPulse_counter = 0;
      while (wf->pulseTimes[darkPulse_counter] >= 20. && wf->pulseTimes[darkPulse_counter] < afpTimeThreshold1) {
        darkPulse_counter++;
      }

      Double_t dcr_waveform = Double_t(darkPulse_counter)/timeWindow; 
      waveformDCRs.push_back(dcr_waveform);

      if (darkPulse_counter != 0) {
        Double_t dcr_waveform_error = dcr_waveform*sqrt(pow(8.e-9/timeWindow,2) + pow(sqrt(darkPulse_counter)/Double_t(darkPulse_counter),2));
        waveformDCRs_errors.push_back(dcr_waveform_error);
      } else if (darkPulse_counter == 0) {
        Double_t dcr_waveform_error = dcr_waveform*sqrt(pow(8.e-9/timeWindow,2) );
        waveformDCRs_errors.push_back(dcr_waveform_error);
      }
    }
  }

  // Find sum of all elements
  Double_t DCRsum = 0., DCRerrorsquared = 0.;
  for (int i = 0; i < waveformDCRs.size(); i++ ) {
    DCRsum += waveformDCRs.at(i);
    DCRerrorsquared += pow(waveformDCRs_errors.at(i),2);
  }  

  // Find average of DCR
  dcr_estimation = DCRsum/waveformDCRs.size();

  std::cout << "Total number of waveformDCRs is " << waveformDCRs.size() << std::endl;

  // Find dcr error
  dcr_error = sqrt(DCRerrorsquared)/waveformDCRs_errors.size();

  std::cout << "Total number of waveformDCRs_errors is " << waveformDCRs_errors.size() << std::endl;
}

void drawHistogram(std::vector<double> data, TH1F* histogram, const char runNumber[], const Int_t y1_in, const Int_t y2_in, const Double_t x1_in, const Double_t x2_in, const Double_t percentage_x, const Double_t percentage_y) {

  // This function will save the histogram to a file, 
  // following a same style, 
  // but changing the histogram, canvas size, user ranges.
  // The definition of the histogram is made at the start
  // of the program. 

  // THE USER RANGES ARE ONLY APPLIED IF x1 != x2 and y1 != y2, RESPECTIVELY. 
  // IF THE PAIRS ARE EQUAL, THE DEFAULT RANGES WILL REMAIN.
  Double_t x1, x2, y1, y2;

  // Fill histogram with data
  for (int i=0; i < data.size(); i++) {
    histogram->Fill(data.at(i));
  }
  
  // Bin width for y label
  float binWidth = histogram->GetBinWidth(0);
  std::cout<<binWidth<<std::endl;
  
  // Write y label
  char yLabel[1024];
  sprintf(yLabel, "No. of events/ %3.1f ns", binWidth);

  // Style
  histogram->SetFillColor(0);
  histogram->SetLineColor(1);
  histogram->SetLineWidth(2);
  
  // Ranges
  // If ranges are equal, the default value
  // will be used for title positioning.
  if (x1_in != x1_in) {
    x1 = x1_in;
    x2 = x2_in;
    histogram->GetXaxis()->SetRangeUser(x1, x2);
  } else {
    x1 = histogram->GetXaxis()->GetXmin();
    x2 = histogram->GetXaxis()->GetXmax();
  }
  if (y1 != y2) {
    y1 = y1_in;
    y2 = y2_in;
    histogram->GetYaxis()->SetRangeUser(y1, y2);
  } else {
    y1 = histogram->GetMaximum();
    y2 = histogram->GetMinimum();
  }

  // Change axes' labels
  histogram->GetYaxis()->SetTitle(yLabel);
  histogram->GetXaxis()->SetTitle("Time (ns)");
  
  // Remove stats box
  histogram->SetStats(0);

  // Calculate position of in-graph title from percentage
  Double_t position_x = percentage_x/100.*(x2-x1);
  Double_t position_y = percentage_y/100.*(y2-y1);

  // Create in-ghaph text
  char inGraphTitle[1024];
  sprintf(inGraphTitle, "Channel %s", runNumber );
  std::cout<<"\n"<<inGraphTitle<<std::endl;
  // Add in-graph text to canvas
  auto *th = new TText(position_x, position_y, inGraphTitle);
  th->SetTextAlign(11); 
  th->SetTextSize(0.04);
  
  // Draw
  histogram->Draw();
  th->Draw();

}
