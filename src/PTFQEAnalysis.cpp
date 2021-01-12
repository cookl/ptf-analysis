#include "PTFQEAnalysis.hpp"

#include "ScanPoint.hpp"
#include "WaveformFitResult.hpp"

#include "TStyle.h"

using namespace std;

PTFQEAnalysis::PTFQEAnalysis( TFile *fout, PTFAnalysis *ptfanalysis ){

  // Turn off stats box
  gStyle->SetOptStat(0);

  vector< double > xbins = ptfanalysis->get_bins( 'x' );
  vector< double > ybins = ptfanalysis->get_bins( 'y' );

  //Create detection efficiency histogram
  pmt0_qe = new TH2D("pmt0_qe", "Detection efficiency", xbins.size()-1, &xbins[0], ybins.size()-1, &ybins[0]);
  pmt0_qe->GetXaxis()->SetTitle("X position [m]");
  pmt0_qe->GetYaxis()->SetTitle("Y position [m]");

  //Loop through events and fill histogram

  //Method 1
  //Loop through ptf_tree
  //for(unsigned int entry=0; entry<ptf_tree->GetEntries(); entry++){
  //  if( entry>10 ) continue;
  //  //Loop through waveforms
  //  ptf_tree->GetEvent( entry );
  //  cout << fitresult->x << "," << fitresult->y << "," << fitresult->z << " haswf: " << fitresult->haswf << " nwaves: " << fitresult->nwaves << endl;
  //}

  //Method 2
  //Loop through scanpoints
  vector< ScanPoint > scanpoints = ptfanalysis->get_scanpoints();
  for(unsigned int iscan=0; iscan<scanpoints.size(); iscan++){
    ScanPoint scanpoint = scanpoints[ iscan ];
    long long entry = scanpoint.get_entry();
    //Loop through fit results
    for(unsigned long long iwavfm=0; iwavfm<scanpoint.nentries(); iwavfm++){
      entry += iwavfm;
      const WaveformFitResult fit_result = ptfanalysis->get_fitresult(iscan, iwavfm);
      //if( iwavfm>9 ) continue;
      //cout << fitresult->x << "," << fitresult->y << "," << fitresult->z << " haswf: " << fitresult->haswf << " nwaves: " << fitresult->nwaves << " scan_entries: " << scanpoint.nentries() << endl;
      pmt0_qe->Fill(fit_result.x, fit_result.y, (double)fit_result.haswf/(double)scanpoint.nentries());
    }
  }

  //fout->cd();
  pmt0_qe->SetDirectory( fout );

}

PTFQEAnalysis::PTFQEAnalysis( TFile *fout, PTFAnalysis *ptfanalysis0, PTFAnalysis *ptfanalysis1 ){

  // Turn off stats box
  gStyle->SetOptStat(0);

  vector< double > xbins = ptfanalysis0->get_bins( 'x' );
  vector< double > ybins = ptfanalysis0->get_bins( 'y' );

  //Create detection efficiency histogram
  pmt0_qe = new TH2D("pmt0_qe", "Detection efficiency", xbins.size()-1, &xbins[0], ybins.size()-1, &ybins[0]);
  pmt0_qe->GetXaxis()->SetTitle("X position [m]");
  pmt0_qe->GetYaxis()->SetTitle("Y position [m]");
  pmt1_qe = (TH2D*) pmt0_qe->Clone("pmt1_qe");
  temp_corr = (TH2D*) pmt0_qe->Clone("temp_corr");
  pmt0_qe_corr = (TH2D*) pmt0_qe->Clone("pmt0_qe_corr");
  temp_corr->SetTitle("Temperature correction");
  pmt0_qe_corr->SetTitle("Detection efficiency (corrected)");
  
  // Vector to store the efficiencies
  // Used to calculate the correction below
  vector< double > v_pmt0_qe;
  
  //Loop through scanpoints
  vector< ScanPoint > scanpoints = ptfanalysis0->get_scanpoints();
  for(unsigned int iscan=0; iscan<scanpoints.size(); iscan++){
    ScanPoint scanpoint = scanpoints[ iscan ];
    v_pmt0_qe.push_back( 0.0 );
    //Loop over scanpoint
    long long entry = scanpoint.get_entry();
    for ( unsigned long long iwavfm = 0; iwavfm < scanpoint.nentries(); ++iwavfm ){
      entry += iwavfm;
      const WaveformFitResult wf = ptfanalysis0->get_fitresult( iscan, iwavfm);
      pmt0_qe->Fill(wf.x, wf.y, (double)wf.haswf/(double)scanpoint.nentries());
      v_pmt0_qe[iscan] += (double)wf.haswf/(double)scanpoint.nentries();
    }
  }
  
  // Vector to store the efficiencies
  // Used to calculate the correction below
  vector< double > v_pmt1_qe;
  
  //Loop through scanpoints
  for(unsigned int iscan=0; iscan<scanpoints.size(); iscan++){
    ScanPoint scanpoint = scanpoints[ iscan ];
    v_pmt1_qe.push_back( 0.0 );
    //Loop over scanpoint
    long long entry = scanpoint.get_entry();
    for ( unsigned iwavfm = 0; iwavfm < scanpoint.nentries(); ++iwavfm ){
      entry += iwavfm;
      const WaveformFitResult wf = ptfanalysis1->get_fitresult( iscan, iwavfm);
      pmt1_qe->Fill(wf.x, wf.y, (double)wf.haswf/(double)scanpoint.nentries());
      v_pmt1_qe[iscan] += (double)wf.haswf/(double)scanpoint.nentries();
    }
  }

  //Calculate temperature correction

  //First calculate average from histogram
  double pmt1_qe_av = 0.0;
  int n_filled = 0;
  for(int ix=1 ix<=pmt1_qe->GetNbinsX(); ix++){
    for(int iy=1; iy<=pmt1_qe->GetNbinsY(); iy++){
      if( pmt1_qe->GetBinContent(ix, iy) < 1e-10 ) continue;
      pmt1_qe_av += pmt1_qe->GetBinContent(ix, iy);
      n_filled++;
    }
  }
  pmt1_qe_av /= (double)n_filled;

  //Now create correction histogram
  for(unsigned int iscan=0; iscan<scanpoints.size(); iscan++){
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
    //Divide rolling average by overall average and store in histogram
    pmt1_qe_av_rolling /= pmt1_qe_av;
    temp_corr->Fill(scanpoint.x(), scanpoint.y(), pmt1_qe_av_rolling);
  }

  //Create corrected QE
  pmt0_qe_corr->Divide(pmt0_qe, temp_corr);

  //fout->cd();
  pmt0_qe->SetDirectory( fout );
  pmt1_qe->SetDirectory( fout );
  temp_corr->SetDirectory( fout );
  pmt0_qe_corr->SetDirectory( fout );

}
