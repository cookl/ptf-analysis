#ifndef TBaselineSing_h
#define TBaselineSing_h

#include <vector>
#include <array>
#include <iostream>
#include "TFile.h"
#include "TTree.h"

// Singleton for retrieving and returning various
// values that are stored in the 'Settings' tree for the mPMT mainboard ROOT files
class BrbSettingsTree
{

protected:
  BrbSettingsTree()
  {
  }

  static BrbSettingsTree* singleton_;

  std::vector<double> fBaselines;

  std::vector<double> fHV_setpoints;

  

public:

  BrbSettingsTree(BrbSettingsTree &other) = delete;

  void operator=(const BrbSettingsTree &) = delete;

  static BrbSettingsTree *Get();

  // Return the baselines we stored
  double GetBaseline(int i){
    if(i < 0 || i >= fBaselines.size()) return 0;
    return fBaselines[i];
  }

  // Return the HV setpoints
  double GetHV(int i){
    if(i < 0 || i >= fHV_setpoints.size()) return 0;
    return fHV_setpoints[i];
  }

  // Load the settings tree from the ROOT TTree called 'Settings'
  int LoadSettingsTree(TFile *file);


};




#endif
