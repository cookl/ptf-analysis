#include "WaveformFitResult.hpp"
#include "ScanPoint.hpp"
#include "Utilities.hpp"
#include "FindCircle.hpp"
#include "TDirectory.h"
#include "TFile.h"
#include "TH2D.h"
#include "TStyle.h"
#include "TCanvas.h"
#include "TGraph.h"//scatter plot
#include "TMinuit.h"//minimisation 
#include "TVirtualFitter.h"
#include <iostream>
#include <vector>
#include <algorithm>
#include <string>
#include <fstream>

using namespace std;


//0-> SK PMT, 1-> Monitor PMT
TH2D* g_pmt0_qe;
TH2D* g_pmt1_diff;
TH2D* g_pmt1_fit;
int g_nb =5;
//_____________________________________________________________________________________
//Function to correct the temperature


Double_t func5P(float qe_sk,float qe_mon_diff, Double_t *par)//Double_t is a data type from root                                                                                      
                                                                                                                                                                                   
{
   Double_t t_correction = qe_sk * ( 1+par[0] * qe_mon_diff+par[1] * pow(qe_mon_diff, 2) + par[2] * pow(qe_mon_diff, 3)
				     + par[3] * pow(qe_mon_diff, 4)+par[4]* pow(qe_mon_diff, 5) );
  return t_correction;
}

Double_t t_model(float qe_mon_diff2, Double_t *par2)//Double_t is a data type from root                                                                                                                                                                                                                                                              
{
   Double_t t_corr= ( 1+par2[0] * qe_mon_diff2+par2[1] * pow(qe_mon_diff2, 2) + par2[2] * pow(qe_mon_diff2, 3)
                                     + par2[3] * pow(qe_mon_diff2, 4)+par2[4]* pow(qe_mon_diff2, 5) );
  return t_corr;
}




//Function to minimize
void fcn5P(Int_t &npar, Double_t *gin, Double_t &f, Double_t *par, Int_t iflag)
{
  double metric =0;
  double metric2 =0;
  for(int ix=1; ix<=g_pmt0_qe->GetNbinsX(); ix++){
    for(int iy=2; iy<=g_pmt0_qe->GetNbinsY(); iy++){
      if( g_pmt0_qe->GetBinContent(ix, iy)<1e-10)continue ;
        /*
	g_pmt1_fit->SetBinContent(ix,iy,0);
	metric += fabs( func5P(g_pmt0_qe->GetBinContent(ix, iy), g_pmt1_diff->GetBinContent(ix, iy), par)
                      - func5P(g_pmt0_qe->GetBinContent(ix, iy-1), g_pmt1_diff->GetBinContent(ix, iy-1), par) );
	metric2 += fabs( pow((func5P(g_pmt0_qe->GetBinContent(ix, iy), g_pmt1_diff->GetBinContent(ix, iy), par)
                           - func5P(g_pmt0_qe->GetBinContent(ix, iy-1), g_pmt1_diff->GetBinContent(ix, iy-1), par)),2));
	continue; //don't fill the bin
	*/
      if( g_pmt0_qe->GetBinContent(ix, iy-1)<1e-10 )continue ;
	/*
	g_pmt1_fit->SetBinContent(ix,iy,0);
        metric += fabs( func5P(g_pmt0_qe->GetBinContent(ix, iy), g_pmt1_diff->GetBinContent(ix, iy), par)
                      - func5P(g_pmt0_qe->GetBinContent(ix, iy-1), g_pmt1_diff->GetBinContent(ix, iy-1), par) );
	metric2 += fabs( pow((func5P(g_pmt0_qe->GetBinContent(ix, iy), g_pmt1_diff->GetBinContent(ix, iy), par)
                           - func5P(g_pmt0_qe->GetBinContent(ix, iy-1), g_pmt1_diff->GetBinContent(ix, iy-1), par)),2));
	continue;
	*/
     
      metric += fabs( func5P(g_pmt0_qe->GetBinContent(ix, iy), g_pmt1_diff->GetBinContent(ix, iy), par)
		      - func5P(g_pmt0_qe->GetBinContent(ix, iy-1), g_pmt1_diff->GetBinContent(ix, iy-1), par) );
      metric2 += fabs( pow((func5P(g_pmt0_qe->GetBinContent(ix, iy), g_pmt1_diff->GetBinContent(ix, iy), par)
			   - func5P(g_pmt0_qe->GetBinContent(ix, iy-1), g_pmt1_diff->GetBinContent(ix, iy-1), par)),2));
      g_pmt1_fit->SetBinContent(ix,iy,func5P(g_pmt0_qe->GetBinContent(ix, iy), g_pmt1_diff->GetBinContent(ix, iy), par)
				- func5P(g_pmt0_qe->GetBinContent(ix, iy-1), g_pmt1_diff->GetBinContent(ix, iy-1), par));        
                     
    }
  }
  std::cout<<"metric "<<metric2<<std::endl;
  f=metric2*10;
  //f = pow(metric,2);
}

//Fit
void Ifit(TVirtualFitter*& gMinuit)
{
  
  
  gMinuit->SetFCN(fcn5P); //Metric to minimize

  Double_t arglist[2];//errors+ parameter value                                                                                                                                     
  // Unused Int_t ierflg = 0;
  
  arglist[0] = 0.0001;
  gMinuit->ExecuteCommand("SET ERR", arglist ,1);//We want to minimize the chi^2                                                                                                     
  arglist[0] = 0;
  gMinuit->ExecuteCommand("SET PRINT", arglist, 1);

  // Set starting values and step sizes for parameters                                                                                                                                 
  static Double_t vstart[5] = {0., 0., 0.,0.,0.}; // only 3 parameter                                                                                                               
  static Double_t step[5] = {0.01,  0.01 , 0.01,0.01,0.01};
  gMinuit->SetParameter(0, "a1", vstart[0], step[0], 0, 0);//number/name/step/limit, lower bound and higher bound, 0 and 0, nothing differen
  gMinuit->SetParameter(1, "a2", vstart[1], step[1], 0, 0);
  gMinuit->SetParameter(2, "a3", vstart[2], step[2], 0, 0);
  gMinuit->SetParameter(3, "a4", vstart[3], step[3], 0, 0);
  gMinuit->SetParameter(4, "a5", vstart[4], step[4], 0, 0);
  // Now ready forarglist[0] = 0.0001;
  
  
  for (int i=4; (i>g_nb-1);i--){
   gMinuit->FixParameter(i);
  } 
   
  arglist[0] = 5000;
  arglist[1] = 10;
  gMinuit->ExecuteCommand("MIGRAD", arglist ,2);

  // Print results                                                                                                                                                                     
  //Double_t amin, edm, errdef;
  //Int_t nvpar, nparx, icstat;
  //gMinuit->mnstat(amin, edm, errdef, nvpar, nparx, icstat);
  //gMinuit->mnprin(1, amin);

}


//_______________________________________________________________________
int main( int argc, char* argv[] ) {

  if ( argc != 3 ){
    std::cerr<<"Usage: ptf_qe_analysis.app ptf_analysis.root run_number\n";
    exit(0);
  }

  // Get utilities
  Utilities utils;

  // Set style
  utils.set_style();

  TFile * fin = new TFile( argv[1], "read" );

  // Opening the output root file
  string outname = string("ptf_qe_analysis_run0") + argv[2] + ".root";
  TFile * fout = new TFile(outname.c_str(), "NEW");
  //Get the scanpoints information
  std::vector< ScanPoint > scanpoints = ReadScanPoints( fin );

  vector< double > xbins = utils.get_bins( scanpoints, 'x' ); // how to use the template correctly
  vector< double > ybins = utils.get_bins( scanpoints, 'y' );

  //Create detection efficiency histogram
  TH2D * pmt0_qe = new TH2D("pmt0_qe", "Detection efficiency", xbins.size()-1, &xbins[0], ybins.size()-1, &ybins[0]);
  pmt0_qe->GetXaxis()->SetTitle("X position [m]");//Create normal plot for all the corrected efficiency
  pmt0_qe->GetYaxis()->SetTitle("Y position [m]");
  TH2D * pmt1_qe = (TH2D*) pmt0_qe->Clone("pmt1_qe");
  TH2D * pmt1_diff = (TH2D*) pmt0_qe->Clone("pmt1_diff");//Calculate the difference
  TH2D * temp_corr = (TH2D*) pmt0_qe->Clone("temp_corr");
  TH2D * temp_corr_min = (TH2D*) pmt0_qe->Clone("temp_corr_min");
  TH2D * pmt_correlation_h = (TH2D*) pmt0_qe->Clone("temp_corr_min");
  TH2D * model_h = (TH2D*) pmt0_qe->Clone("model_h");
  TH2D * pmt0_qe_corr = (TH2D*) pmt0_qe->Clone("pmt0_qe_corr");
  temp_corr->SetTitle("Difference from the rolling average monitor PMT");
  temp_corr_min->SetTitle("Detection efficiency (Corrected), using Minuit");
  pmt1_diff->SetTitle("Difference from the rolling average");
  pmt0_qe_corr->SetTitle("Detection efficiency (Corrected)");
  pmt_correlation_h->SetTitle("Correlation between the PMT");
  model_h->SetTitle("Modelization of the temperature fluctuation");
  //______________________________________________________________________________________________________________
  // get the waveform fit TTree for PMT0
  TTree * tt0 = (TTree*)fin->Get("ptfanalysis0");// how to create tree for the wave form
  WaveformFitResult * wf = new WaveformFitResult;
  wf->SetBranchAddresses( tt0 );

  // Vector to store the efficiencies
  // Used to calculate the correction below
  vector< double > v_pmt0_qe; 

  //Loop through scanpoints
  for(unsigned int iscan=0; iscan<scanpoints.size(); iscan++){
    if (iscan%1000==0) std::cout<<"Filling PMT0 histograms for iscan = "<<iscan<<" / "<<scanpoints.size()<<std::endl;
    ScanPoint scanpoint = scanpoints[ iscan ];
    v_pmt0_qe.push_back( 0.0 ); // store the data of the efficiency
    //Loop over scanpoint
    for ( unsigned iev = 0; iev < scanpoint.nentries(); ++iev ){
      tt0->GetEvent( scanpoint.get_entry() + iev );
      bool haswf = utils.HasWaveform( wf, 0 );//only use data that has a waveform
      pmt0_qe->Fill(wf->x, wf->y, (double)haswf/(double)scanpoint.nentries()); //
      v_pmt0_qe[iscan] += (double)haswf/(double)scanpoint.nentries();
    }
  }

  //________________________________________________________________________________________________________________
  // Get the waveform fit TTree for PMT1
  TTree * tt1 = (TTree*)fin->Get("ptfanalysis1");
  wf->SetBranchAddresses( tt1 );

  // Vector to store the efficiencies
  // Used to calculate the correction below
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
      pmt1_qe->Fill(wf->x, wf->y, (double)haswf/(double)scanpoint.nentries());
      v_pmt1_qe[iscan] += (double)haswf/(double)scanpoint.nentries();

    }
  }
  //__________________________________________________________________________________________________________________
  //Calculate temperature correction

  //First calculate average from histogram
  double pmt1_qe_av = 0.0;
  vector <double> pmt1_qe_tot;   
  int n_filled = 0;
  for(int ix=1; ix<=pmt1_qe->GetNbinsX(); ix++){
    for(int iy=1; iy<=pmt1_qe->GetNbinsY(); iy++){
      if( pmt1_qe->GetBinContent(ix, iy) < 1e-10 ) continue; // do the operation as long as the statement is true
      pmt1_qe_av += pmt1_qe->GetBinContent(ix, iy);
      pmt1_qe_tot.push_back(pmt1_qe->GetBinContent(ix, iy));// added line
      n_filled++;
    }
  }
  pmt1_qe_av /= (double)n_filled;// add all bins value and then divide by total -> gives the mean -> <QesK>
  cout << "PMT1 average qe: " << pmt1_qe_av << endl; //1- Calculate average value of the PMT
  //Change average to be from a single reference scan

  //_______________________________________________________________________________________________________
  //Addition to the script

  //Vector containing all the efficiency for pmt0
  vector <double> pmt0_qe_tot;
  for(int ix=1; ix<=pmt0_qe->GetNbinsX(); ix++){
    for(int iy=1; iy<=pmt0_qe->GetNbinsY(); iy++){
      if( pmt0_qe->GetBinContent(ix, iy) < 1e-10 ) continue; // do the operation as long as the statement is true                                               
      pmt0_qe_tot.push_back(pmt0_qe->GetBinContent(ix, iy));
    }

  }
  /*
  //Calculating the differential
  for(int ix=1; ix<=pmt1_qe->GetNbinsX(); ix++){
    for(int iy=1; iy<=pmt1_qe->GetNbinsY(); iy++){
      if( pmt1_qe->GetBinContent(ix, iy) < 1e-10 ) continue;
      double bincontent=pmt1_qe->GetBinContent(ix, iy);
      pmt1_diff->SetBinContent(ix, iy, bincontent-pmt1_qe_av);
    }
  }
  */
  


  //_________________________________________________________________________________________________________

  //Now create correction histogram
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
    pmt1_qe_av_rolling /= (double)nbins_rolling;//
    
    //Divide rolling average by overall average and store in histogram
    //pmt1_qe_av_rolling /= pmt1_qe_av;
    pmt1_diff->Fill(scanpoint.x(), scanpoint.y(), pmt1_qe_av_rolling-pmt1_qe_av);
    //pmt1_qe_av_rolling /= pmt1_qe_av; 
    temp_corr->Fill(scanpoint.x(), scanpoint.y(), pmt1_qe_av_rolling);
  } // create the corrected histogram with temperature
  //Create corrected QE
  
  pmt0_qe_corr->Divide(pmt0_qe, temp_corr);
  
  //Remove data outside circle
  TH2D* pmt0_qe_corr_grad;
  Circle_st circ = find_circle_max_grad( pmt0_qe_corr, pmt0_qe_corr_grad, 0.5 );
  zero_outside_circle( pmt0_qe_corr, circ );
  zero_outside_circle( pmt0_qe, circ );

  //_________________________________________________________________________________________________________

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
  //fcn()
    //std::cout<<"a1 "<<metric<<std::endl;
  //std::cout<< "a1:"<<p0 << std::end;
  //____________________________________________________________________________
  //Now calculate the temperature correction using the value of the fit and compute histogram

  double parameter[5];
  parameter[0]= p0;
  parameter[1]=p1;
  parameter[2]=p2;
  parameter[3]=p3 ;
  parameter[4]=p4 ;

  ofstream myfile("parameter.txt",ios::app);
  myfile <<argv[2]<<" "<< p0<<" "<< p1<< " "<<p2<< " "<<p3<< " "<<p4<<endl;
  myfile<<argv[2]<< " "<<err0<<" "<< err1<< " "<<err2<< " "<<err3<< " "<<err4<<endl;
  myfile.close();

  
  for(int ix=1; ix<=pmt1_qe->GetNbinsX(); ix++){
    for(int iy=1; iy<=pmt1_qe->GetNbinsY(); iy++){
      double pmt_SK=pmt0_qe->GetBinContent(ix, iy);
      double pmt_mn=pmt1_diff->GetBinContent(ix, iy);
       temp_corr_min->SetBinContent(ix, iy, func5P(pmt_SK, pmt_mn,parameter));
       
    }  
  }
  
  for(int ix=1; ix<=pmt1_diff->GetNbinsX(); ix++){
    for(int iy=1; iy<=pmt1_diff->GetNbinsY(); iy++){
      model_h->SetBinContent(ix, iy,(t_model(pmt1_diff->GetBinContent(ix, iy),parameter)));

   }
  }
 
  double metric_roll;
  double metric_mn;
  for(int ix=1; ix<=pmt0_qe_corr->GetNbinsX(); ix++){
    for(int iy=2; iy<=pmt0_qe_corr->GetNbinsY(); iy++){
      if( pmt0_qe_corr->GetBinContent(ix, iy) < 1e-10 ) continue;
      if( pmt0_qe_corr->GetBinContent(ix, iy-1) < 1e-10 ) continue;
      metric_roll+=pow(pmt0_qe_corr->GetBinContent(ix, iy)-pmt0_qe_corr->GetBinContent(ix, iy-1),2);
      metric_mn+=pow(temp_corr_min->GetBinContent(ix, iy)-temp_corr_min->GetBinContent(ix, iy-1),2);

    }

  }

  
  
  std::cout<<"Metric of the older correction "<<metric_roll<<endl;
  std::cout<<"Metric of the new correction "<<metric_mn<<endl;
  
  //_________________________________________________________________________________________________________
  //Scatter plot to look at the correlation between the 2 PMT
  TGraph* pmt_correlation=new TGraph(pmt1_qe_tot.size(),&pmt0_qe_tot[0],&pmt1_qe_tot[0]);
  pmt_correlation->GetXaxis()->SetTitle("PMT efficiency");//Create normal plot for all the corrected efficiency                                                                     
  pmt_correlation->GetYaxis()->SetTitle("Monitor PMT efficiency");
  pmt_correlation->SetTitle("Correlation between PMT");
  pmt_correlation->SetMarkerStyle(1);
  pmt_correlation->SetMarkerColor(2);
  pmt_correlation->GetXaxis()->SetRangeUser(0.1,0.3);
  pmt_correlation->GetYaxis()->SetRangeUser(0.4,0.6);
  for (int ix=1;ix<=pmt0_qe_tot.size();ix++){
    for(int iy=1; iy<=pmt0_qe_tot.size(); iy++){
      pmt_correlation_h->Fill(pmt0_qe->GetBinContent(ix, iy),pmt1_qe->GetBinContent(ix, iy));     
    }
  }
  pmt_correlation_h->GetXaxis()->SetRangeUser(0.1,0.3);
  pmt_correlation_h->GetYaxis()->SetRangeUser(0.4,0.6);


  //_________________________________________________________________________________________________________
  //outfile->cd();
  //xpmt_correlation->SetDirectory(fout);
  
  pmt0_qe->SetDirectory(fout );
  pmt1_qe->SetDirectory( fout );
  pmt1_diff->SetDirectory( fout );
  temp_corr->SetDirectory( fout );
  temp_corr_min->SetDirectory( fout );
  pmt0_qe_corr->SetDirectory( fout );
  pmt_correlation_h->SetDirectory( fout );
  g_pmt1_fit->SetDirectory( fout );
  //model_h->SetDirectory( fout );

  //Set plot ranges
  int run = stoi(argv[2]);
  if( run<4525 ){
    pmt0_qe->SetMinimum(0.0);
    temp_corr_min->SetMinimum(0);
    pmt0_qe_corr->SetMinimum(0.0);
    pmt0_qe->SetMaximum(0.22);
    pmt0_qe_corr->SetMaximum(0.22);
    temp_corr_min->SetMaximum(0.22);
  }
  else{
    pmt0_qe->SetMinimum(0.0);
    pmt0_qe_corr->SetMinimum(0.0);
    pmt0_qe->SetMaximum(0.3);
    pmt0_qe_corr->SetMaximum(0.3);
    temp_corr_min->SetMaximum(0.3);
  }
  pmt1_qe->SetMinimum(0.3);
  pmt1_qe->SetMaximum(0.6);
  //temp_corr->SetMinimum(0.7);
  //temp_corr->SetMaximum(1.2);
  //temp_corr->SetMinimum(1.25);
  //temp_corr->SetMaximum(1.45);
  pmt_correlation_h->SetMaximum(70);
  pmt1_diff->SetMinimum(-0.15);
  pmt1_diff->SetMaximum(0.15);
  g_pmt1_fit->SetMinimum(-0.1);
  g_pmt1_fit->SetMaximum(0.1);
  //Make plots
  TCanvas* c = new TCanvas("canvas");
  string plotname;
  g_pmt0_qe->Draw("colz0");
  plotname = string("~/projects/def-pdeperio/vgousy/ptf/ptf-analysis-2/plot_dir/1par/ptf_qe_analysis_run0")+argv[2]+"_pmt0.pdf";
  c->SaveAs(plotname.c_str(),"pdf");
  //TCanvas* c2 = new TCanvas();
  pmt1_qe->Draw("colz0");
  gPad->Modified();
  gPad->Update();
  plotname = string("~/projects/def-pdeperio/vgousy/ptf/ptf-analysis-2/plot_dir/1par/ptf_qe_analysis_run0")+argv[2]+"_pmt1.pdf";
  c->SaveAs(plotname.c_str(),"pdf");
  //TCanvas* c3 = new TCanvas();
  temp_corr->Draw("colz0");
  gPad->Modified();
  gPad->Update();
  plotname = string("~/projects/def-pdeperio/vgousy/ptf/ptf-analysis-2/plot_dir/1par/ptf_qe_analysis_run0")+argv[2]+"_tempcorr.pdf";
  c->SaveAs(plotname.c_str(),"pdf");
  //TCanvas* c4 = new TCanvas("c4");
  pmt0_qe_corr->Draw("colz0");
  gPad->Modified();
  gPad->Update();
  plotname = string("~/projects/def-pdeperio/vgousy/ptf/ptf-analysis-2/plot_dir/1par/ptf_qe_analysis_run0")+argv[2]+"_pmt0corr.pdf";
  c->SaveAs(plotname.c_str(),"pdf");
  //TCanvas* c5 = new TCanvas("c5"); 
  pmt1_diff->Draw("colz0");
  gPad->Modified();
  gPad->Update();
  plotname = string("~/projects/def-pdeperio/vgousy/ptf/ptf-analysis-2/plot_dir/1par/ptf_qe_analysis_run0")+argv[2]+"_diff.pdf";
  c->SaveAs(plotname.c_str(),"pdf");


  //TCanvas* c6 = new TCanvas("c6");                                                                                                                                                 
  pmt_correlation->Draw("ap");
  gPad->Modified();
  gPad->Update();
  plotname = string("~/projects/def-pdeperio/vgousy/ptf/ptf-analysis-2/plot_dir/1par/ptf_qe_analysis_run0")+argv[2]+"_corr.pdf";
  c->SaveAs(plotname.c_str(),"pdf");


  //TCanvas* c6 = new TCanvas("c6");   
  temp_corr_min->Draw("colz0");
  gPad->Modified();
  gPad->Update();
  plotname = string("~/projects/def-pdeperio/vgousy/ptf/ptf-analysis-2/plot_dir/1par/ptf_qe_analysis_run0")+argv[2]+"_tempcorr_min.pdf";
  c->SaveAs(plotname.c_str(),"pdf");
  
  //TCanvas* c6 = new TCanvas("c6");                                                                                                                                                 
  pmt_correlation_h->Draw("colz0");
  gPad->Modified();
  gPad->Update();
  plotname = string("~/projects/def-pdeperio/vgousy/ptf/ptf-analysis-2/plot_dir/1par/ptf_qe_analysis_run0")+argv[2]+"_pmtcorrelation.pdf";
  c->SaveAs(plotname.c_str(),"pdf");


   //TCanvas* c6 = new TCanvas("c6");                                                                                                                                               
                                                                                                                                                                                     
  g_pmt1_fit->Draw("colz0");
  gPad->Modified();
  gPad->Update();
  plotname = string("~/projects/def-pdeperio/vgousy/ptf/ptf-analysis-2/plot_dir/1par/ptf_qe_analysis_run0")+argv[2]+"fit.pdf";
  c->SaveAs(plotname.c_str(),"pdf");

  model_h->Draw("colz0");
  gPad->Modified();
  gPad->Update();
  plotname = string("~/projects/def-pdeperio/vgousy/ptf/ptf-analysis-2/plot_dir/1par/ptf_qe_analysis_run0")+argv[2]+"model.pdf";
  c->SaveAs(plotname.c_str(),"pdf");

  //Write and close output file
  fout->Write();
  fout->Close();

  cout << "Done" << endl; 

  return 0;
}



