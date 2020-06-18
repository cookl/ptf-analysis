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

  vector<int> phidgets = {0, 1, 3, 4};
  vector<PTF::PMTChannel> activeChannels = {};
 
  PTF::Wrapper wrapper = PTF::Wrapper(16384, 70, activeChannels, phidgets);

  unordered_set<int> skipLines = {};// {962,1923,2884,5240,6201,9611,10572,11533,12494,13455,15811,16771};

  wrapper.openFile(root_f, "scan_tree");

  cerr << "Num entries: " << wrapper.getNumEntries() << endl << endl;
  
  for (unsigned int i = 0; i < wrapper.getNumEntries(); i++) {

    auto reading_T = wrapper.getDataForCurrentEntry(PTF::Gantry0);


      auto reading  = wrapper.getReadingTemperature(PTF::T);

      cerr << reading.int_1<<reading.ext_1<<reading,ext_2<<endl;
      
  }
}

