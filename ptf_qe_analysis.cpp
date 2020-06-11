#include "WaveformFitResult.hpp"
#include "ScanPoint.hpp"
#include "Utilities.hpp"
#include "FindCircle.hpp"
#include "TFile.h"
#include "TH2D.h"
#include "TStyle.h"
#include "TCanvas.h"
#include "TArc.h"

#include <iostream>
#include <vector>
#include <algorithm>
#include <string>

using namespace std;

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

  // get the scanpoints information
  std::vector< ScanPoint > scanpoints = ReadScanPoints( fin );

  vector< double > xbins = utils.get_bins( scanpoints, 'x' );
  vector< double > ybins = utils.get_bins( scanpoints, 'y' );

  //Create detection efficiency histogram
  TH2D * pmt0_qe = new TH2D("pmt0_qe", "Detection efficiency", xbins.size()-1, &xbins[0], ybins.size()-1, &ybins[0]);
  pmt0_qe->GetXaxis()->SetTitle("X position [m]");
  pmt0_qe->GetYaxis()->SetTitle("Y position [m]");
  TH2D * pmt1_qe = (TH2D*) pmt0_qe->Clone("pmt1_qe");
  TH2D * temp_corr = (TH2D*) pmt0_qe->Clone("temp_corr");
  TH2D * pmt0_qe_corr = (TH2D*) pmt0_qe->Clone("pmt0_qe_corr");
  temp_corr->SetTitle("Temperature correction");
  pmt0_qe_corr->SetTitle("Detection efficiency (corrected)");

  // get the waveform fit TTree for PMT0
  TTree * tt0 = (TTree*)fin->Get("ptfanalysis0");
  WaveformFitResult * wf = new WaveformFitResult;
  wf->SetBranchAddresses( tt0 );
  
  // Vector to store the efficiencies
  // Used to calculate the correction below
  vector< double > v_pmt0_qe;
  
  //Loop through scanpoints
  for(unsigned int iscan=0; iscan<scanpoints.size(); iscan++){
    if (iscan%1000==0) std::cout<<"Filling PMT0 histograms for iscan = "<<iscan<<" / "<<scanpoints.size()<<std::endl;
    ScanPoint scanpoint = scanpoints[ iscan ];
    v_pmt0_qe.push_back( 0.0 );
    //Loop over scanpoint
    for ( unsigned iev = 0; iev < scanpoint.nentries(); ++iev ){
      tt0->GetEvent( scanpoint.get_entry() + iev );
      bool haswf = utils.HasWaveform( wf, 0 );
      pmt0_qe->Fill(wf->x, wf->y, (double)haswf/(double)scanpoint.nentries());
      v_pmt0_qe[iscan] += (double)haswf/(double)scanpoint.nentries();
    }
  }

  // get the waveform fit TTree for PMT0
  TTree * tt1 = (TTree*)fin->Get("ptfanalysis1");
  wf->SetBranchAddresses( tt1 );
  
  // Vector to store the efficiencies
  // Used to calculate the correction below
  vector< double > v_pmt1_qe;
  
  //Loop through scanpoints
  for(unsigned int iscan=0; iscan<scanpoints.size(); iscan++){
    if (iscan%1000==0) std::cout<<"Filling PMT1 histograms for iscan = "<<iscan<<" / "<<scanpoints.size()<<std::endl;
    ScanPoint scanpoint = scanpoints[ iscan ];
    v_pmt1_qe.push_back( 0.0 );
    //Loop over scanpoint
    for ( unsigned iev = 0; iev < scanpoint.nentries(); ++iev ){
      tt1->GetEvent( scanpoint.get_entry() + iev );
      bool haswf = utils.HasWaveform( wf, 1 );
      pmt1_qe->Fill(wf->x, wf->y, (double)haswf/(double)scanpoint.nentries());
      v_pmt1_qe[iscan] += (double)haswf/(double)scanpoint.nentries();
    }
  }

  //Calculate temperature correction

  //First calculate average from histogram
  double pmt1_qe_av = 0.0;
  int n_filled = 0;
  for(int ix=1; ix<=pmt1_qe->GetNbinsX(); ix++){
    for(int iy=1; iy<=pmt1_qe->GetNbinsY(); iy++){
      if( pmt1_qe->GetBinContent(ix, iy) < 1e-10 ) continue;
      pmt1_qe_av += pmt1_qe->GetBinContent(ix, iy);
      n_filled++;
    }
  }
  pmt1_qe_av /= (double)n_filled;
  cout << "PMT1 average qe: " << pmt1_qe_av << endl;
  //Change average to be from a single reference scan
  //Results will then be comparable
  pmt1_qe_av = 0.496535;

  //Now create correction histogram
  for(unsigned int iscan=0; iscan<scanpoints.size(); iscan++){
    ScanPoint scanpoint = scanpoints[ iscan ];
    //Calculate rolling average efficiency over 2*ntemps+1 scan points
    int ntemps = 5;
    double pmt1_qe_av_rolling = 0.0;
    int nbins_rolling = 0;
    int imin = max(0, (int)iscan-ntemps);
    int imax = min((int)v_pmt1_qe.size()-1, (int)iscan+ntemps);
    for(int i=imin; i<=imax; i++){
      if( v_pmt1_qe[i] < 0.01 ) continue;
      pmt1_qe_av_rolling += v_pmt1_qe[i];
      nbins_rolling++;
    }
    pmt1_qe_av_rolling /= (double)nbins_rolling;
    //Divide rolling average by overall average and store in histogram
    pmt1_qe_av_rolling /= pmt1_qe_av;
    temp_corr->Fill(scanpoint.x(), scanpoint.y(), pmt1_qe_av_rolling);
  }

  //Create corrected QE
  pmt0_qe_corr->Divide(pmt0_qe, temp_corr);
  
  //Remove data outside circle
  TH2D* pmt0_qe_corr_grad;
  Circle_st circ = find_circle_max_grad( pmt0_qe_corr, pmt0_qe_corr_grad, 0.8 );
  zero_outside_circle( pmt0_qe_corr, circ );

  //outfile->cd();
  pmt0_qe->SetDirectory( fout );
  pmt1_qe->SetDirectory( fout );
  temp_corr->SetDirectory( fout );
  pmt0_qe_corr->SetDirectory( fout );

  //Set plot ranges
  int run = stoi(argv[2]);
  if( run<4525 ){
    pmt0_qe->SetMinimum(0.0);
    pmt0_qe_corr->SetMinimum(0.0);
    pmt0_qe->SetMaximum(0.22);
    pmt0_qe_corr->SetMaximum(0.22);
  }
  else{
    pmt0_qe->SetMinimum(0.0);
    pmt0_qe_corr->SetMinimum(0.0);
    pmt0_qe->SetMaximum(0.3);
    pmt0_qe_corr->SetMaximum(0.3);
  }
  pmt1_qe->SetMinimum(0.3);
  pmt1_qe->SetMaximum(0.6);
  temp_corr->SetMinimum(0.7);
  temp_corr->SetMaximum(1.2);
  //temp_corr->SetMinimum(1.25);
  //temp_corr->SetMaximum(1.45);

  //Make plots
  TCanvas* c = new TCanvas("canvas");
  string plotname;
  pmt0_qe->Draw("colz0");
  plotname = string("ptf_qe_analysis_run0")+argv[2]+"_pmt0.pdf";
  c->SaveAs(plotname.c_str(),"pdf");
  //TCanvas* c2 = new TCanvas();
  pmt1_qe->Draw("colz0");
  gPad->Modified();
  gPad->Update();
  plotname = string("ptf_qe_analysis_run0")+argv[2]+"_pmt1.pdf";
  c->SaveAs(plotname.c_str(),"pdf");
  //TCanvas* c3 = new TCanvas();
  temp_corr->Draw("colz0");
  gPad->Modified();
  gPad->Update();
  plotname = string("ptf_qe_analysis_run0")+argv[2]+"_tempcorr.pdf";
  c->SaveAs(plotname.c_str(),"pdf");
  //TCanvas* c4 = new TCanvas("c4");
  pmt0_qe_corr->Draw("colz0");
  gPad->Modified();
  gPad->Update();
  plotname = string("ptf_qe_analysis_run0")+argv[2]+"_pmt0corr.pdf";
  c->SaveAs(plotname.c_str(),"pdf");
  TArc *arc1 = new TArc(circ.xc, circ.yc, circ.r );
  arc1->SetLineColor(kRed);
  arc1->SetLineWidth(3);
  arc1->SetFillStyle(0);
  arc1->Draw();
  pmt0_qe_corr_grad->Draw("colz0");
  gPad->Modified();
  gPad->Update();
  plotname = string("ptf_qe_analysis_run0")+argv[2]+"_pmt0corr_grad.pdf";
  c->SaveAs(plotname.c_str(),"pdf");

  //Write and close output file
  fout->Write();
  fout->Close();

  cout << "Done" << endl; 

  return 0;
}
