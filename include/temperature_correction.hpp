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
