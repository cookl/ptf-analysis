# PTF Analysis

This is the codebase for analysing the ROOT files produced from the (PMT test facility) PTF. It handles the loading of files, accessing the data (e.g. PMT waveform samples and Phidget readings), fitting of waveforms to produce a tree of fitted parameters, and analyses of the waveform fits (including charge, timing, and efficiency measurements).

## Table of Contents

### 1. [Directory Layout](#directory_layout)
### 2. [Installation](#installation)
### 3. [Getting the PTF data](#data)
### 4. [Usage](#usage)

## Directory Layout <a id="directory_layout"></a>

```bash
.
+-- bin/                                # Location of compiled executables
+-- include/                            # Header files
+-- macros/                             # ROOT macros to produce plots from the output of the analysis executables
+-- magnetic-field/                     # Python scripts to process the output of field_to_csv
+-- obj/                                # Location for compiled .o files
+-- ptf_bfield/                         # Standalone analysis of the predicted PTF magnetic field
+-- src/                                # Source files
+-- Makefile                            # Makefile to build executables
+-- field_to_csv.cpp
+-- ptf.config.dat                      # A configuration file that sets analysis options
+-- ptf_analysis.cpp
+-- ptf_charge_analysis.cpp
+-- ptf_field_analysis.cpp
+-- ptf_qe_analysis.cpp
+-- ptf_timing_analysis.cpp
+-- ptf_ttree_analysis.cpp
```

## Installation <a id="installation"></a>

To download the repository use:

`git clone https://<username>@bitbucket.org/ttriumfdaq/ptf-analysis-2.git`

## Getting the PTF data <a id="data"></a>

The PTF MIDAS DAQ produces output files. These can be converted to a ROOT TTree with the `rootana/libAnalyzer/analyzer_convert_ptf_scan_to_rootTree.cxx` script in the `ptf-online-converter` repository. This step it typically completed automatically on the PTF machine when a scan completes. 

The MIDAS files are located here on the PTF machine:  
`~/online/data/`

The ROOT trees are located here on the PTF machine:  
`~/online/rootfiles/`

## Usage <a id="usage"></a>

To compile the code run `make`. To build new analyses add them to the `Makefile` following the example of the existing analyses.

The `ptf_analysis` executable fits the PMT waveforms and produces a ROOT file that contains a TTree with the fitted parameter values. The fitted parameter values can then be analysed by the `ptf_charge_analysis`, `ptf_qe_analysis` and `ptf_timing_analysis` executables. The command to run the code from the root directory is:  
`./bin/ptf_analysis.dat filename.root run_number config_file`  
The `run_number` argument is to produce an output file with a name specific to the run.  

The `ptf_ttree_analysis` executable is a demonstration of how the TTree produced by `ptf_analysis` could be accessed. The command to run the code from the root directory is:  
`./bin/ptf_ttree_analysis.app ptf_analysis.root`

The `ptf_charge_analysis` executable reads the fitted waveforms from `ptf_analysis` and computes the charge of the events. The command to run the code from the root directory is:  
`./bin/ptf_charge_analysis.app ptf_analysis.root run_number [T/F/I]`  
Where the T/F/I is for True to do/not do circle fit to find PMT, I to cut inside circle (default T).  
The `run_number` argument is to produce an output file with a name specific to the run.  

The `ptf_qe_analysis` executable reads the fitted waveforms from `ptf_analysis` and calculates the detection efficiency for the PMT. The command to run the code from the root directory is:  
`./bin/ptf_qe_analysis.app ptf_analysis.root run_number`  
The `run_number` argument is to produce an output file with a name specific to the run.  

The `ptf_timing_analysis` executable reads the fitted waveforms from `ptf_analysis` and calculates the timing response for the PMT. The command to run the code from the root directory is:  
`./bin/ptf_timing_analysis.app ptf_analysis.root run_number`  
The `run_number` argument is to produce an output file with a name specific to the run.  

The `ptf_field_analysis` executable reads the data from Phidget04 which is fixed inside the Helmholtz coils and plots its magnetic field values as the scan progresses. This provides an indication of the field stability over the course of a run. The command to run the script from the root directory is:  
`./bin/ptf_field_analysis.app /data/directory run_number`  
The `run_number` argument is to produce an output file with a name specific to the run.  

The `field_to_csv` analysis reads the magnetic field values from the Phidgets and outputs them to a csv file for analysis by the python scripts in the `magnetic-field` directory.  

The `mpmt_analysis` executable fits the mPMT PMT waveforms and produces a ROOT file that contains a TTree with the fitted parameter values. The fitted parameter values can then be analysed by the further executables. The command to run the code from the root directory is:  
`./bin/mpmt_analysis.dat filename.root run_number mpmt.config.dat`  
The `run_number` argument is to produce an output file with a name specific to the run.  


## The different classes

These are the classes that are most important to understand.

```bash
+-- wrapper               Handles the loading of files and accessing the data
+-- Configuration         Loads options from configuration file
+-- PTFErrorBarAnalysis   For calculating the error bar size to use on the waveforms
+-- PTFAnalysis           For doing analysis of all of the waveforms, and keep track of scan points, stores results in TTree
+-- WaveformFitResult     Structure to hold one waveform fit result
+-- ScanPoint             Holds location of scan point, first entry number in TTree of scan point, and number of waveforms
```

## The wrapper class

It handles loading the files and provides a simple way to access the data. A simple example of how you might use it can be found in `wrapper_demo.cpp`.

## Data Types

Here is a brief overview of the data types you'll use (all in "wrapper.hpp", in `namespace PTF`):

### Gantry (`enum Gantry`)

```c++
enum Gantry {
  Gantry0 = 0,
  Gantry1 = 1
};
```

Just an enum for which gantry you want to reference.

### PMTType (`enum PMTType`)

```c++
enum PMTType {
  Hamamatsu_R3600_PMT = 0,
  PTF_Monitor_PMT = 1,
  PTF_Receiver_PMT = 2,
  mPMT_REV0_PMT = 3,
  mPMT_Monitor_PMT = 4
};
```

An enum for the PMT type


### PMT (`struct PMT`)

```c++
typedef struct PMT {
  int pmt;
  int channel;
  PMTType type;
} PMTChannel;
```

A structure which maps PMTs to their channel and type. `int pmt` is the PMT's number, `int channel` is the digitizer channel used to read the info, and `PMTType type` is the type of PMT.


### Phidget Reading (`struct PhidgetReading`)

```c++
typedef struct PhidgetReading {
  double Bx[10];
  double By[10];
  double Bz[10];
} PhidgetReading;
```

Data read for the phidget. Mostly you'll probably just want index 0 of each, which is the magnetic field.


### Ganty Data (`struct GantryData`)

```c++
typedef struct GantryData {
  double x;
  double y;
  double z;
  double theta;
  double phi;
} GantryData;
```

Contains the position information for a gantry.


### DigitizerModel (`enum DigitizerModel`)

```c++
enum DigitizerModel {
  PTF_CAEN_V1730 = 0,
  mPMT_DIGITIZER = 1
};
```

An enum for the digitizer type

### Digitizer Settings (`struct Digitizer`)

```c++
typedef struct Digitizer {
  DigitizerModel model;
  int samplingRate;
  double fullScaleRange;
  int resolution;
};
```

A structure to store the digitizer settings. Units are MS/s (mega samples per second) for the sampling rate, Vpp (Voltage peak-to-peak) for the full scale range, and bits for the resolution.


## Methods of `PTF::Wrapper`


Here are the methods of `PTF::Wrapper`:

- `Wrapper(size_t maxSamples, size_t sampleSize, const std::vector<PMT>& activePMTs, const std::vector<int>& phidgets, const std::vector<Gantry>& gantries, DigitizerModel digi)`
    - Constructs a wrapper object and prepares to read the given PMTs and phidgets.
- `Wrapper(size_t maxSamples, size_t sampleSize, const std::vector<PMT>& activePMTs, const std::vector<int>& phidgets, const std::vector<Gantry>& gantries, DigitizerModel digi, const std::string& fileName, const std::string& treeName = "scan_tree")`
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
- `GantryData getDataForCurrentEntry(Gantry whichGantry) const`
    - Gets gantry position info for a given gantry. Throws if no file is open.
- `PhidgetReading getReadingFOrPhidget(int phidget) const`
    - Gets the phidget data for the given phidget and current entry. Throws `InvalidPhidget` if the phidget wasn't registered, and `NoFileIsOpen` if no file is open.
- `Digitizer getDigitizerSettings() const`
    - Gets the digitizer settings.
