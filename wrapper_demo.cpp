#include <vector>

#include "wrapper.hpp"


using namespace std;


int main(void) {
  // decide which channels we'd like
  vector<PTF::PMTChannel> channels = {
    {1, 3} // this is saying we want pmt #1, which is on channel 3.
  };

  // decide which phidgets we'd like to read
  vector<int> phidgets = {1, 3, 4};

  // initialize the wrapper
  auto wrapper = PTF::Wrapper(
    16384, // the maximum number of samples
    34, // the size of one sample
    channels,
    phidgets
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
    GantryPos gantryData = wrapper.getDataForCurrentEntry(PTF::Gantry1);

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
