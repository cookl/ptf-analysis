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
    baselines = std::vector<double>(20,0);
  }

  static BrbSettingsTree* singleton_;

  std::vector<double> baselines;

  std::vector<double> hv_setpoints;

  

public:

  BrbSettingsTree(BrbSettingsTree &other) = delete;

  void operator=(const BrbSettingsTree &) = delete;

  static BrbSettingsTree *Get();

  // Return the baselines we stored
  double GetBaseline(int i){
    if(i < 0 || i >= baselines.size()) return 0;
    return baselines[i];
  }

  // Return the HV setpoints
  double GetHV(int i){
    if(i < 0 || i >= hv_setpoints.size()) return 0;
    return hv_setpoints[i];
  }

  // Load the settings tree from the ROOT TTree called 'Settings'
  int LoadSettingsTree(TFile *file);


};




#endif
