#ifndef __PTF_WRAPPER__
#define __PTF_WRAPPER__

#include <climits>
#include <vector>
#include <array>
#include <string>
#include <utility>
#include <exception>
#include <algorithm>
#include <unordered_map>
#include <iostream>
#include <iostream>
#include <fstream>


#include "TROOT.h"
#include "TFile.h"
#include "TTree.h"
#include "TBranch.h"

#include "config.hpp"

/// Classes to to help with reading in PTF data
/// PTF::PmtChannel           holds pmt number and channel number
///                           pmt seems to be an arbitrary number of user 
///                           channel is the PMT channel in the data (digitizer channel)
///
/// PTF::PhigetReading        Holds magnetic field readings
///
/// PTF::GantryData           Holds gantry data
///
/// PTF::PMTSet               Holds PMT data as pointer to array of doubles, and branch from input
///
/// PTF::PhigetSet            Holds one phidget reading, and branches from input
///
/// PTF::GantrySet            Holds one GantryData object, and branches from input
///
/// PTF::Wrapper              Main helper class for accessing the PTF data
///                           See comments on public methods for help

// using namespace std;
// using namespace boost;


namespace PTF {

enum Gantry {
  Gantry0 = 0,
  Gantry1 = 1
};


enum PMTType {
  Hamamatsu_R3600_PMT = 0,
  PTF_Monitor_PMT = 1,
  PTF_Receiver_PMT = 2,
  mPMT_REV0_PMT = 3,
  mPMT_Monitor_PMT = 4,
  Reference = 5
};



struct PMT {

  int pmt;
  int channel;
  PMTType type;
};


struct PhidgetReading {
  double Bx[150];
  double By[150];
  double Bz[150];
  double Ax[150];
  double Ay[150];
  double Az[150];
};

}


struct PhidgetSet {
    PTF::PhidgetReading data;
    TBranch*       branchX{nullptr};
    TBranch*       branchY{nullptr};
    TBranch*       branchZ{nullptr};
    TBranch*       branchaccx{nullptr};
    TBranch*       branchaccy{nullptr};
    TBranch*       branchaccz{nullptr};
  };




struct GantryData {

  double x;
  double y;
  double z;
  double theta;
  double phi;
};

struct Temperature_r {
  //double int_1;
  //double ext_1;
  double ext_2;
};

struct Phidget_acce00 {
  double acc_x;
  double acc_y;
  double acc_z;
};


struct Timing {
   Double_t time_c;
};

struct PMTSet {
  int      channel;
  PTF::PMTType  type;
  double*  data{nullptr};
  TBranch* branch{nullptr};

};

//struct PhidgetSet {
//  PTF::PhidgetReading data;
//  TBranch*       branchX{nullptr};
//  TBranch*       branchY{nullptr};
//  TBranch*       branchZ{nullptr};
//};

struct GantrySet {
  PTF::Gantry     gantry;
  GantryData data;
  TBranch*   branchX{nullptr};
  TBranch*   branchY{nullptr};
  TBranch*   branchZ{nullptr};
  TBranch*   branchTheta{nullptr};
  TBranch*   branchPhi{nullptr};
};


enum DigitizerModel {
  PTF_CAEN_V1730 = 0,
  mPMT_DIGITIZER = 1
};


struct Digitizer {
  DigitizerModel model;
  int samplingRate;
  double fullScaleRange;
  int resolution;
};




struct Wrapper {
	Wrapper(unsigned long long maxSamples, unsigned long long sampleSize, const std::vector<PTF::PMT>& activePMTs, const std::vector<int>& phidgets, const std::vector<PTF::Gantry>& gantries, DigitizerModel digi);
    Wrapper(unsigned long long maxSamples, unsigned long long sampleSize, const std::vector<PTF::PMT>& activePMTs, const std::vector<int>& phidgets, const std::vector<PTF::Gantry>& gantries, DigitizerModel digi, const std::string& fileName, const std::string& treeName = "scan_tree");
	~Wrapper();


public:
  // Public interface
  
  // Opens the selected file, and loads the first entry
  void openFile(const std::string& fileName, const std::string& treeName = "scan_tree");
  bool isFileOpen() const;
  // Closes the currently open file and deletes the tree.
  // Does nothing if the file is already closed 
  void closeFile();

  // Load the BRB settings tree; returns -1 on failure
  int LoadBrbSettingsTree();
  
  // Returns -1 on not found
  int getChannelForPmt(int pmt) const;
  // Ditto
  int getPmtForChannel(int channel) const;

  // Functions for getting/setting data entry
  // throws `NoFileIsOpen` if no file is open.
  unsigned long long getCurrentEntry() const;
  // ditto
  unsigned long long getNumEntries() const;
  void setCurrentEntry(unsigned long long entry);  // Throws exception on invalid entry

  /* Reading data */

  // Throws on file not open
  unsigned long long getNumSamples() const;

  // Gets the data for a given sample on the current
  // Throws on invalid sample or file not open
  double* getPmtSample(int pmt, unsigned long long sample) const;

  // Returns the length of the samples
  int getSampleLength() const;

  GantryData getDataForCurrentEntry(PTF::Gantry whichGantry) const;

  PTF::PhidgetReading getReadingForPhidget(int phidget) const;
  
  Temperature_r getReadingTemperature() const;
  

  //Phidget_acce00 getReadingAcceleration() const;


  
  Timing getReadingTime() const;

  Digitizer getDigitizerSettings() const {return digiData;}

private:
  TFile* file{0};
  TTree* tree{0};
  unsigned long long maxSamples;
  unsigned long long sampleSize;
  unsigned long long entry{ULONG_MAX};

  // data
  std::unordered_map<int, PMTSet*>     pmtData;
  std::unordered_map<int, PhidgetSet*> phidgetData;
  std::unordered_map<int, GantrySet*>  gantryData;
  Digitizer digiData;
  

  Temperature_r Temp;
  Timing ti;
  //*Phidget_acce00 ACC;
  
  unsigned long long    numEntries;
  unsigned long long    numSamples;

  /* Private methods */

  // Gets the data pointer for the specified pmt
  // Returns nullptr if not found
  double* getDataForPmt(int pmt) const;

  // Sets the pointers in the tree to the newly opened tree
  // Returns false on failure, true on success
  bool setDataPointers();

  // Sets all branch pointers to nullptr
  // Returns false on failure, true on success
  bool unsetDataPointers();

};


namespace Exceptions {
  class FileDoesNotExist : public std::runtime_error {
  public:
    FileDoesNotExist(const std::string& msg) : runtime_error(msg) {}
  };

  class InvalidTreeName : public std::runtime_error {
  public:
    InvalidTreeName(const std::string& msg) : runtime_error(msg) {}
  };

  class NoFileIsOpen : public std::runtime_error {
  public:
    NoFileIsOpen() : runtime_error("No file is open.") {}
  };

  class EntryOutOfRange : public std::runtime_error {
  public:
    EntryOutOfRange() : runtime_error("Entry is out of range.") {}
  };

  class SampleOutOfRange : public std::runtime_error {
  public:
    SampleOutOfRange() : runtime_error("Sample is out of range.") {}
  };

  class InvalidPhidget : public std::runtime_error {
  public:
    InvalidPhidget() : runtime_error("No such phidget.") {}
  };

  class InvalidPMT : public std::runtime_error {
  public:
    InvalidPMT() : runtime_error("No such PMT.") {}
  };

  class InvalidGantry : public std::runtime_error {
  public:
    InvalidGantry() : runtime_error("No such gantry.") {}
  };

  class DataPointerError : public std::runtime_error {
  public:
    DataPointerError() : runtime_error("Error while setting data pointers.") {}
  };

  class CSVFileError : public std::runtime_error {
  public:
    CSVFileError() : runtime_error("Error while trying to open CSV file.") {}
  };
} // end namespace Exceptions


 // end namespace PTF


// template <typename T, typename Ts>
// bool has(Ts _variant) {
//   if (boost::get<T>(&_variant)) {
//     return true;
//   } else {
//     return false;
//   }
// }


#endif // __PTF_WRAPPER__
