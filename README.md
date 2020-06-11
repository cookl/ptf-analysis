# PTF Analysis

This is a simple C++ library for reading the ROOT files produced from the PTF. It handles loading the files and provides a simple way to access the data.

It can be build with `make wrapper.o` for a debug build (uses `-g3 -Og`), or with `make wrapper.o RELEASE=TRUE` for an optimized build (uses `-g -O2`). You can then link your program with `g++ -o myprog myprog.cpp wrapper.o {flags here} -I$(ROOTSYS)/include/root -L$(ROOTSYS)/lib/root -lCore -lHist -lRIO -lTree -lGpad`. `clang++` works as well, but make sure that ROOT, `wrapper.o` and your program are all compiled with the same compiler.

Here's a simple example of how you might use it:

```c++
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
```

## Data Types

Here is a brief overview of the data types you'll use (all in "wrapper.hpp", in `namespace PTF`):

### Gantry (`enum Gantry`)

```c++
enum Gantry {
  Gantry0,
  Gantry1
};
```

Just an enum for which gantry you want to reference.


### PMT Channel (`struct PMTChannel`)

```c++
typedef struct PMTChannel {
  int pmt;
  int channel;
} PMTChannel;
```

Represents a pair used to map PMTs to their channel. `int pmt` is the PMT's number, and `int channel` is the channel used to read the info.


### Phidget Reading (`struct PhidgetReading`)

```c++
typedef struct PhidgetReading {
  double Bx[10];
  double By[10];
  double Bz[10];
} PhidgetReading;
```

Data read for the phidget. Mostly you'll probably just want index 0 of each, which is the magnetic field.


### Ganty Position (`struct GantryPos`)

```c++
typedef struct GantryPos {
  double x;
  double y;
  double z;
  double theta;
  double phi;
} GantryPos;
```

Contains the position information for a gantry.


## Methods of `PTF::Wrapper`


Here are the methods of `PTF::Wrapper`:

- `Wrapper(size_t maxSamples, size_t sampleSize, const std::vector<PMTChannel>& activeChannels, const std::vector<int>& phidgets)`
    - Constructs a wrapper object and prepares to read the given channels and phidgets.
- `Wrapper(size_t maxSamples, size_t sampleSize, const std::vector<PMTChannel>& activeChannels, const std::vector<int>& phidgets, const std::string& fileName, const std::string& treeName = "scan_tree")`
    - Constructs a wrapper object like above, but immediately opens a file and loads a scan tree ("scan_tree" by default).
- `void openFile(const std::string& fileName, const std::string& treeName = "scan_tree")`
    - Opens a file, as described above.
- `bool isFileOpen() const`
    - Returns `true` if the wrapper currently has a file loaded, `false` otherwise.
- `void closeFile()`
    - Closes the current file, if one is open. Does nothing if there is no file open.
- `int getChannelForPmt(int pmt) const`
    - Gets the channel for the given PMT. Returns -1 if it's not found.
- `int getPmtForChannel(int channel) const`
    - Does the inverse of the above, also returning -1 if not found.
- `size_t getCurrentEntry() const`
    - Gets the current entry. Throws `NoFileIsOpen` if no file is open.
- `size_t getNumEntries() const`
    - Gets the total number of entries. Throws `NoFileIsOpen` if no file is open.
- `void setCurrentEntry(size_t entry)`
    - Sets the current entry. Throws `NoFileIsOpen` if no file is open, and `EntryOutOfRange` if the entry is too large.
- `size_t getNumSamples() const`
    - Gets the number of samples for the current entry. Throws `NoFileIsOpen` if no file is open.
- `double* getPmtSample(int pmt, size_t sample) const`
    - Returns a pointer to an array of length `sampleSize` which is the sample for the PMT on the current entry. Throws `SampleOutOfRange` if the sample is too large, `InvalidPMT` if the PMT can't be found and `NoFileIsOpen` if no file is open.
- `int getSampleLength() const`
    - Returns `sampleSize`.
- `GantryPos getDataForCurrentEntry(Gantry whichGantry) const`
    - Gets gantry position info for a given gantry. Throws if no file is open.
- `PhidgetReading getReadingFOrPhidget(int phidget) const`
    - Gets the phidget data for the given phidget and current entry. Throws `InvalidPhidget` if the phidget wasn't registered, and `NoFileIsOpen` if no file is open.
