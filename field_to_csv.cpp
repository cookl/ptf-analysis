#include "wrapper.hpp"

#include <string>
#include <iostream>
#include <ostream>
#include <fstream>
#include <unordered_set>


using namespace std;


int main(int argc, char** argv) {
  if (argc != 2) {
    cerr << "Enter the run number.";
  }

  const string run_no = argv[1];
  string root_f = "/neut/datasrv2a/jmgwalker/ptf/ptf-analysis-2/develop/ptf-analysis-2/data/out_run0" + run_no + ".root";
  string csv_f  = "/neut/datasrv2a/jmgwalker/ptf/ptf-analysis-2/develop/ptf-analysis-2/field_to_csv/out_run0" + run_no + ".csv";

  vector<int> phidgets = {0, 1, 3, 4};
  vector<PTF::PMTChannel> activeChannels = {};
  //   {0, 3},
  //   {1, 4},
  //   {2, 5},
  //   {3, 6},
  //   {4, 7}, 
  //   {5, 8},
  //   {6, 9},
  //   {7, 10}
  // };
  PTF::Wrapper wrapper = PTF::Wrapper(16384, 70, activeChannels, phidgets);

  unordered_set<int> skipLines = {};// {962,1923,2884,5240,6201,9611,10572,11533,12494,13455,15811,16771};

  wrapper.openFile(root_f, "scan_tree");
  ofstream csv(csv_f);

  cerr << "Num entries: " << wrapper.getNumEntries() << endl << endl;

  csv << "X0,Y0,Z0, Bx0,By0,Bz0, Bx1,By1,Bz1, Bx3,By3,Bz3, Bx4,By4,Bz4" << endl;

  uint32_t lines = 0;
  const uint32_t freq = 100;
  
  for (unsigned int i = 0; i < wrapper.getNumEntries(); i++) {
    // cerr << "Entry " << i;
    if (i % freq == 0 || i == wrapper.getNumEntries() - 1) {
      cerr << "Entry " << i << "/" << wrapper.getNumEntries() << "\u001b[34;1m (" << (((double)i)/wrapper.getNumEntries()*100) << "%)\u001b[0m\033[K";
      if (skipLines.find(i) != skipLines.end()) {
        cerr << "\u001b[31;1m Skipping...\u001b[0m\r";
        continue;
      } else {
        cerr << "\r";
      }
    }

    if (skipLines.find(i) != skipLines.end()) continue;

    lines++;
    wrapper.setCurrentEntry(i);

    auto location = wrapper.getDataForCurrentEntry(PTF::Gantry0);

    csv << location.x << "," << location.y << "," << location.z << ",";

    for (int phidget : phidgets) {

      auto reading  = wrapper.getReadingForPhidget(phidget);

      csv << reading.Bx[0] << "," << reading.By[0] << "," << reading.Bz[0];
      
      if (phidget != 4) {
        csv << ",";
      }
    }

    csv << endl;
  }

  cerr << endl << "Done. Wrote " << lines << " lines.";
  csv.close();
}

