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
  string root_f = "/neut/data19/vincent/ptf-analysis-2/out_run0" + run_no + ".root";
 string csv_f  = "/neut/data19/vincent/ptf-analysis_test/ptf-analysis/acc" + run_no + ".csv";

  vector<int> phidgets = {0, 1};
  vector<PTF::PMT> activePMTs = {};
  //   {0, 3},                                                                                                                                                                                                                                
  //   {1, 4},                                                                                                                                                                                                                                
  //   {2, 5},                                                                                                                                                                                                                                
  //   {3, 6},                                                                                                                                                                                                                                
  //   {4, 7},                                                                                                                                                                                                                                
  //   {5, 8},                                                                                                                                                                                                                                
  //   {6, 9},                                                                                                                                                                                                                                
  //   {7, 10}                                                                                                                                                                                                                                
  // ;                                                                                                                                                                                                
  vector<PTF::Gantry> gantries = {PTF::Gantry0, PTF::Gantry1};                                       
  Wrapper wrapper = Wrapper(16384, 70, activePMTs, phidgets, gantries, PTF_CAEN_V1730);

  unordered_set<int> skipLines = {};// {962,1923,2884,5240,6201,9611,10572,11533,12494,13455,15811,16771};                                                                                                                                    

  wrapper.openFile(root_f, "scan_tree");
  ofstream csv(csv_f);

  cerr << "Num entries: " << wrapper.getNumEntries() << endl << endl;

  csv << "time,phid_0_accx,phid_0_accy,phid_0_accz,phid_1_accx,phid_1_accy,phid_1_accz,accx,accy,accz, " << endl;

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
  
  
  auto time_before=wrapper.getReadingTime();
  wrapper.setCurrentEntry(0);
  auto time_after=wrapper.getReadingTime();

  csv <<time_before.time_c-time_after.time_c << "," ;

  

    for (int phidget : phidgets) {

	       auto reading  = wrapper.getReadingForPhidget(phidget);

	       csv << reading.Ax[0] << "," << reading.Ay[0] << "," << reading.Az[0];
		   
		   if (phidget != 100) {
		   		csv << ",";
		                  }
		   if (phidget == 100) {
		       auto acceleration  = wrapper.getReadingAcceleration();
			    csv <<  acceleration.acc_x << "," << reading.acc_y << "," << reading.acc_z;
			   		                  }

		}				  
	     csv << endl;
	   }
   

	   cerr << endl << "Done. WroteZEGOAL" << lines << " lines.";
	   csv.close();
	 }



		
