#include "wrapper.hpp"
#include "MeanRMSCalc.hpp"
#include "ErrorBarAnalysis.hpp"
#include "WaveformFitResult.hpp"
#include "ScanPoint.hpp"
#include "PTFAnalysis.hpp"
#include "PTFQEAnalysis.hpp"
#include "Utilities.hpp"
#include "temperature_correction.hpp"



#include <string>
#include <iostream>
#include <ostream>
#include <fstream>
#include <sstream>
#include <exception>
#include <stdio.h>
#include <stdlib.h>
#include <vector>
#include <algorithm>
#include <cmath>
#include <unordered_set>

#include "TH1D.h"
#include "TH2D.h"
#include "TFile.h"
#include "TF1.h"
#include "TVirtualFitter.h"
#include "TFitter.h"
#include "TMath.h"
#include "TStyle.h"
#include "TDirectory.h"
#include "TFile.h"
#include "TH2D.h"
#include "TStyle.h"
#include "TCanvas.h"
#include "TGraph.h"//scatter plot
#include "TLegend.h"
#include "TGaxis.h"
#include "TMultiGraph.h"
#include "TFrame.h"
using namespace std;

double round(double var)
{
  double value = (int)(var * 100 + .5);
  return (double)value / 100;

}
static bool comparison (double i, double j){ return (fabs( i-j ) < 1e-5); }

int main( int argc, char* argv[] ) {

  if ( argc != 3 ){
    std::cerr<<"Usage: temperature_sensor_analysis.app ptf_analysis.root run_number\n";
    exit(0);
  }

  // Get utilities
  Utilities utils;

  // Set style
  utils.set_style();
	
  TFile * fin = new TFile( argv[1], "read" );

  // Opening the output root file
  string outname = string("temperature_sensor_") + argv[2] + ".root";
  TFile * fout = new TFile(outname.c_str(), "NEW");
  //Get the scanpoints information
  std::vector< ScanPoint > scanpoints = ReadScanPoints( fin );
  //Go fin the values for temperature and timing
  vector< double > temp = utils.get_bins( scanpoints, 'e' );
  vector< double > tempbins = utils.get_bins( scanpoints, 'b' );
  vector< double > timing = utils.get_bins( scanpoints, 't' );
  
  vector< double > xbins = utils.get_bins( scanpoints, 'x' ); 
  vector< double > ybins = utils.get_bins( scanpoints, 'y' );
  //Plots in order to use the fit
  TH2D * pmt0_qe = new TH2D("pmt0_qe", "Detection efficiency", xbins.size()-1, &xbins[0], ybins.size()-1, &ybins[0]);
  TH2D * pmt1_qe = (TH2D*) pmt0_qe->Clone("pmt1_qe");
  TH2D * pmt1_diff = (TH2D*) pmt0_qe->Clone("pmt1_diff");//Calculate the difference  

  TH2D * temp_plot = new TH2D("temp_plot", "Temperature ", xbins.size()-1, &xbins[0], ybins.size()-1, &ybins[0]);
  temp_plot->GetXaxis()->SetTitle("X position [m]");//Create normal plot for all the corrected efficiency
  temp_plot->GetYaxis()->SetTitle("Y position [m]"); 

  //Get detection efficiency from monitor PMT and SK PMT
  TTree * tt0 = (TTree*)fin->Get("ptfanalysis0");// how to create tree for the wave form                                                                                           
  WaveformFitResult * wf = new WaveformFitResult;                                                                                                                                  
  wf->SetBranchAddresses( tt0 );    
  vector< double > v_pmt0_qe;
  vector< double > QEBins;
  vector< double > count;
  for(unsigned int iscan=0; iscan<scanpoints.size(); iscan++){
    if (iscan%1000==0) std::cout<<"Filling PMT0 histograms for iscan = "<<iscan<<" / "<<scanpoints.size()<<std::endl;
    ScanPoint scanpoint = scanpoints[ iscan ];
    v_pmt0_qe.push_back( 0.0 ); // store the data of the efficiency                                                                                                                  
    //Loop over scanpoint                                                                                                                                                            
    int k=0;
    for ( unsigned iev = 0; iev < scanpoint.nentries(); ++iev ){
      tt0->GetEvent( scanpoint.get_entry() + iev );
      bool haswf = utils.HasWaveform( wf, 0 );//only use data that has a waveform                                                                                                
      v_pmt0_qe[iscan] += (double)haswf/(double)scanpoint.nentries();
      pmt0_qe->Fill(wf->x, wf->y, (double)haswf/(double)scanpoint.nentries());
     //qe_effect_h->Fill(temperature_ext2[iev], Time[iev], (double)haswf/(double)scanpoint.nentries());
      count.push_back((double)haswf/(double)scanpoint.nentries());
      k++;
    
    }
  }
  

  TTree * tt1 = (TTree*)fin->Get("ptfanalysis1");
  wf->SetBranchAddresses( tt1 );
  
  vector< double > v_pmt1_qe;
  //Loop through scanpoints                                                                                                                                                          
  for(unsigned int iscan=0; iscan<scanpoints.size(); iscan++){
    if (iscan%1000==0) std::cout<<"Filling PMT1 histograms for iscan = "<<iscan<<" / "<<scanpoints.size()<<std::endl;
    ScanPoint scanpoint = scanpoints[ iscan ];
    v_pmt1_qe.push_back( 0.0 );// what is exactly push back function ?                                                                                                               
    //Loop over scanpoint                                                                                                                                                            
    for ( unsigned iev = 0; iev < scanpoint.nentries(); ++iev ){ // unsigned just means positif value only                                                                           
      tt1->GetEvent( scanpoint.get_entry() + iev );
      bool haswf = utils.HasWaveform( wf, 1 );
      //temp_plot->Fill(wf->x, wf->y,temperature_ext2[iev]); 
      v_pmt1_qe[iscan] += (double)haswf/(double)scanpoint.nentries();
      pmt1_qe->Fill(wf->x, wf->y, (double)haswf/(double)scanpoint.nentries());
    }
  }



  //Classified the QE data
  QEBins=v_pmt0_qe;
  for (int i=0;i<QEBins.size();i++) {
    QEBins[i]=round(QEBins[i]);
  
    }
 sort(  QEBins.begin(),  QEBins.end() );
 QEBins.erase(unique(  QEBins.begin(),  QEBins.end() ,comparison),QEBins.end() );


 //Calculation of the fit
 double pmt1_qe_av = 0.496535;
 for(unsigned int iscan=0; iscan<scanpoints.size(); iscan++){//Start of the loop witn iscan                                                                                         
   ScanPoint scanpoint = scanpoints[ iscan ];
   //Calculate rolling average efficiency over 2*ntemps scan points                                                                                                                 
   int ntemps = 5;
   double pmt1_qe_av_rolling = 0.0;
   int nbins_rolling = 0;
   int imin = max(0, (int)iscan-ntemps);
   int imax = min((int)v_pmt1_qe.size(), (int)iscan+ntemps);
   for(int i=imin; i<imax; i++){
     if( v_pmt1_qe[i] < 0.01 ) continue;
     pmt1_qe_av_rolling += v_pmt1_qe[i];
     nbins_rolling++;
   }
   pmt1_qe_av_rolling /= (double)nbins_rolling;                                                                                                                                   

   
   pmt1_diff->Fill(scanpoint.x(), scanpoint.y(), pmt1_qe_av_rolling-pmt1_qe_av);
    } 
 

 g_pmt0_qe = (TH2D*) pmt0_qe->Clone("g_pmt0_qe");
 g_pmt1_diff = (TH2D*) pmt1_diff->Clone("g_pmt1_diff");
 g_pmt1_fit = (TH2D*) pmt1_diff->Clone("g_pmt1_fit");
 g_pmt1_fit->Reset();

 Double_t p0;
 Double_t p1;
 Double_t p2;
 Double_t p3;
 Double_t p4;
 Double_t err0;
 Double_t err1;
 Double_t err2;
 Double_t err3;
 Double_t err4;
 TVirtualFitter *gMinuit = TVirtualFitter::Fitter(0,g_nb);  //initialize TMinuit with a maximum of 3 params                                                                         
 Ifit(gMinuit);
 p0 = gMinuit->GetParameter(0);
 p1 = gMinuit->GetParameter(1);
 p2 = gMinuit->GetParameter(2);
 p3 = gMinuit->GetParameter(3);
 p4 = gMinuit->GetParameter(4);
 err0 = gMinuit->GetParError(0);
 err1 = gMinuit->GetParError(1);
 err2 = gMinuit->GetParError(2);
 err3 = gMinuit->GetParError(3);
 err4 = gMinuit->GetParError(4);


 double parameter[5];
 parameter[0]= p0;
 parameter[1]=p1;
 parameter[2]=p2;
 parameter[3]=p3 ;
 parameter[4]=p4 ;

 //Creates the histogram for QE and Temp
 TH2D * nbevents = new TH2D("nbevents", "Temperature correlation with QE_SK ",20,QEBins[0],QEBins.back() ,20,19,tempbins.back()  );
 nbevents->GetXaxis()->SetTitle("QE");
 nbevents->GetYaxis()->SetTitle("Temperature");
 TH2D * nbevents_MN = new TH2D("nbevents", "Temperature correlation with QE_MN ",20,0.3,0.6 ,20,19,tempbins.back()  );
 nbevents_MN->GetXaxis()->SetTitle("QE");
 nbevents_MN->GetYaxis()->SetTitle("Temperature");

 for (int j=0; j<=temp.size();j++){
     nbevents->Fill(v_pmt0_qe[j],temp[j]);
     nbevents_MN->Fill(v_pmt1_qe[j],temp[j]);

  }
  
 
  vector<double> QE_corr;
   for(int ix=1; ix<=pmt1_qe->GetNbinsX(); ix++){
     for(int iy=1; iy<=pmt1_qe->GetNbinsY(); iy++){
       //We choose if we plot the correction of the corrected value
       double pmt_SK=pmt0_qe->GetBinContent(ix, iy);
       double pmt_mn=pmt1_diff->GetBinContent(ix, iy);
       //QE_corr.push_back(func5P(pmt_SK, pmt_mn,parameter));
       QE_corr.push_back(t_model(pmt1_diff->GetBinContent(ix, iy),parameter));                                                                                                     
     }
   }

   TH2D * nbevents_corr = new TH2D("nbevents_corr", "Temperature correlation with QE_SK ",20,QEBins[0],QEBins.back() ,20,19,tempbins.back()  );
   nbevents_corr->GetXaxis()->SetTitle("QE");//Create normal plot for all the corrected efficiency                                                                                                                                                                                                                                                                
   nbevents_corr->GetYaxis()->SetTitle("Temperature");
   
   for (int j=0; j<=temp.size();j++){
     nbevents_corr->Fill(QE_corr[j],temp[j]);

   }


   nbevents->SetDirectory(fout );
   TCanvas* c1 = new TCanvas("canvas");
   string plotname;
   nbevents->Draw("colz0");
   
   plotname = string("~/Desktop/ptf-analysis/temp_his_2_SK")+argv[3]+".pdf";
   c1->SaveAs(plotname.c_str(),"pdf");
   nbevents_MN->Draw("colz0");
   plotname = string("~/Desktop/ptf-analysis/temp_his_2_MN")+argv[3]+".pdf";
   c1->SaveAs(plotname.c_str(),"pdf");
   nbevents_corr->Draw("colz0");
   plotname = string("~/Desktop/ptf-analysis/temp_his_2_corr")+argv[3]+".pdf";
   c1->SaveAs(plotname.c_str(),"pdf");
   

   //Create the scatter plots
   /*
  
  <vector> string variables;
  for (int i=0;i<=variables.size();i++){
  
  }





    */
   TGraph*QE_SK = new TGraph(temp.size(),&timing[0],&v_pmt0_qe[0]);
   TGraph*QE_MN = new TGraph(temp.size(),&timing[0],&v_pmt1_qe[0]);
   TGraph*QE_correction = new TGraph(temp.size(),&timing[0],&QE_corr[0]);
   TGraph* Temper = new TGraph(temp.size(),&timing[0],&temp[0]);

   QE_SK->SetMarkerColor(2);
   QE_MN->SetMarkerColor(3);
   Temper->SetMarkerColor(4);
   QE_correction->SetMarkerColor(6);

   QE_SK->SetTitle("QE_SK");
   QE_MN->SetTitle("QE_MN");
   QE_correction->SetTitle("QE_corrected");

   QE_MN->SetMarkerStyle(2);
   QE_correction->SetMarkerStyle(2);
   Temper->SetMarkerStyle(2);


   QE_MN->GetYaxis()->SetRangeUser(0,1.2);
   QE_SK->GetYaxis()->SetRangeUser(0,1.2);
   QE_correction->GetYaxis()->SetRangeUser(0,1.2);
   QE_SK->SetMarkerStyle(2);
   //We want to overlay to pad in order to have to different scale
   TPad*pad1 = new TPad("QE","QE",0,0,1,1);
   pad1->SetGrid();
   pad1->Draw();
   pad1->cd();

   QE_MN->Draw("ap1");                                                                                                                                                       
   QE_SK->Draw("p1 same"); 
   QE_correction->Draw("p1 same");


   QE_SK->GetYaxis()->SetTitle("PMT QE");
   QE_MN->GetYaxis()->SetTitle("PMT QE");


   c1->cd();
   //Creating the Second pad
   TPad *overlay = new TPad("Temperature","Temperature",0,0,1,1);
   
   overlay->SetFillStyle(4000);
   overlay->SetFillColor(0);
   overlay->SetFrameFillStyle(4000);
   overlay->Draw();
   overlay->cd();
   
   Temper->GetXaxis()->SetTitle("Time");
   Temper->GetYaxis()->SetTitle("Temperature");
   Temper->Draw("APY+");//Options to plot on the other axis

   //auto leg = new TLegend(0.2,0.5,0.5,0.7);
   // leg->SetFillColor(2);
   //  leg->AddEntry("temperature_data","Int","p");
   //  leg->AddEntry("temperature_data_2","ext_1","p");
   //  leg->AddEntry("temperature_data_3","ext_2","p");
   //  leg->Draw();
   //c1->BuildLegend();
   gPad->Modified();
   gPad->Update();
   plotname = string("~/Desktop/ptf-analysis/qe_effect")+argv[3]+".pdf";
   c1->SaveAs(plotname.c_str(),"pdf");
   
   fout->Write();
   fout->Close();
  } 
