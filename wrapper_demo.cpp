#include <vector>

#include "wrapper.hpp"


using namespace std;


int main(void) {
  // decide which PMTs we'd like
  vector<PTF::PMT> activePMTs = {
    {0,0,PTF::Hamamatsu_R3600_PMT} // this is saying we want pmt #0, which is on channel 0 and is the SK PMT.
  };

  // decide which phidgets we'd like to read
  vector<int> phidgets = {1, 3, 4};

  // decide which gantries we'd like to include
  vector<PTF::Gantry> gantries = {PTF::Gantry0, PTF::Gantry1};

  // initialize the wrapper
  auto wrapper = PTF::Wrapper(
    6000, // the maximum number of samples
    70, // the size of one sample
    activePMTs,
    phidgets,
    gantries
  );

  // now we can open our file
  wrapper.openFile("/path/to/file.root");

  cout << "There are " << wrapper.getNumEntries() << " entries." << endl;

  // we can iterate over all the entries

  for (size_t i = 0; i < wrapper.getNumEntries(); i++) {
    wrapper.setCurrentEntry(i);

    // get data from phidget 3
    PhidgetReading phidgetReading = wrapper.getReadingForPhidget(3);

    // get data from gantry 1
    GantryData gantryData = wrapper.getDataForCurrentEntry(PTF::Gantry1);

    // see how many samples there are for the current entry
    auto numSamples = wrapper.getNumSamples();

    for (size_t sample = 0; sample < numSamples; sample++) {
      // Gets a pointer to the data for PMT 1 for this sample
      // It's an array with the length of one sample, set above, in this case 34.
      double* data = getPmtSample(1, sample);
      // do something with data
    }
  }

  // wrapper is automatically deallocated when it goes out of scope here, and its destructor cleans up memory
}
