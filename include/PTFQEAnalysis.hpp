#ifndef __PTFQEANALYSIS__
#define __PTFQEANALYSIS__

#include "TFile.h"
#include "TH2D.h"

#include "PTFAnalysis.hpp"

using namespace std;

class PTFQEAnalysis {
public:
  PTFQEAnalysis(){}; // default constructor
  PTFQEAnalysis( TFile *fout, PTFAnalysis *ptfanalysis );
  PTFQEAnalysis( TFile *fout, PTFAnalysis *ptfanalysis0, PTFAnalysis *ptfanalysis1 );
  ~PTFQEAnalysis(){
    delete pmt0_qe;
    delete pmt1_qe;
    delete temp_corr;
    delete pmt0_qe_corr;
  }

private:
  TH2D *pmt0_qe{nullptr};
  TH2D *pmt1_qe{nullptr};
  TH2D *temp_corr{nullptr};
  TH2D *pmt0_qe_corr{nullptr};

};

#endif // __PTFQEANALYSIS__
