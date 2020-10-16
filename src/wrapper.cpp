#include "wrapper.hpp"


using namespace std;
using namespace PTF;


Wrapper::Wrapper(unsigned long long maxSamples, unsigned long long sampleSize, const vector<PMT>& activePMTs, const vector<int>& phidgets, const vector<Gantry>& gantries, DigitizerModel digi)
  : maxSamples(maxSamples), sampleSize(sampleSize)
{
  for (auto pmt : activePMTs) {
    double* data = new double[maxSamples * sampleSize];
    PMTSet* pmtSet  = new PMTSet();
    pmtSet->channel = pmt.channel;
    pmtSet->type = pmt.type;
    pmtSet->data    = data;
    pmtData[pmt.pmt] = pmtSet;
  }
  for (auto phidget : phidgets) {
    PhidgetSet* pSet = new PhidgetSet();
    phidgetData[phidget] = pSet;
  }
  for (auto gantry : gantries) {
    GantrySet* gSet = new GantrySet();
    gSet->gantry = gantry;
    gantryData[gantry] = gSet;
  }
  digiData.model = digi;
  switch( digi ) {
    case PTF_CAEN_V1730:
      digiData.samplingRate = PTF_CAEN_V1730_SAMPLE_RATE;
      digiData.fullScaleRange = PTF_CAEN_V1730_FULL_SCALE_RANGE;
      digiData.resolution = PTF_CAEN_V1730_RESOLUTION;
      break;
    case mPMT_DIGITIZER:
      digiData.samplingRate = mPMT_DIGITIZER_SAMPLE_RATE;
      digiData.fullScaleRange = mPMT_DIGITIZER_FULL_SCALE_RANGE;
      digiData.resolution = mPMT_DIGITIZER_RESOLUTION;
      break;
  }
}


Wrapper::Wrapper(unsigned long long maxSamples, unsigned long long sampleSize, const vector<PMT>& activePMTs, const vector<int>& phidgets, const vector<Gantry>& gantries, DigitizerModel digi, const string& fileName, const string& treeName)
  : Wrapper(maxSamples, sampleSize, activePMTs, phidgets, gantries, digi) {
  openFile(fileName, treeName);
}


Wrapper::~Wrapper() {
  return;
  // for some reason doing cleanup below causes
  // program to crash as exiting?
  for (auto pmt : pmtData) {
    delete[] pmt.second->data;
    delete pmt.second;
  }
  
  for (auto phidget : phidgetData) {
    delete phidget.second;
  }

  for (auto gantry : gantryData) {
    delete gantry.second;
  }

  if (tree) {
    unsetDataPointers();
    delete tree;
    tree = nullptr;
  }
  if (file) {
    file->Close();
    delete file;
    file = nullptr;
  }
}


/* Private functions */


double* Wrapper::getDataForPmt(int pmt) const {
  // Could be replaced with binary search, but probably list is small enough to not matter
  auto res = pmtData.find(pmt);

  if (res != pmtData.end()) {
    return res->second->data;
  }
  else {
    return nullptr;
  }
}


bool Wrapper::setDataPointers() {
  if (tree == nullptr || file == nullptr)
    return false;
  
  // Set PMT branches
  char branchName[64];
  for (auto pmt : pmtData) {
    snprintf(branchName, 64, PMT_CHANNEL_FORMAT, pmt.second->channel);
    pmt.second->branch = nullptr;
    pmt.second->branch = tree->GetBranch(branchName);
    if (pmt.second->branch == nullptr) {
      return false;
    }
    pmt.second->branch->SetAddress(pmt.second->data);
  }

  // Set phidget branches
  for (auto phidget : phidgetData) {
    snprintf(branchName, 64, PHIDGET_FORMAT_X, phidget.first);
    phidget.second->branchX = nullptr;
    phidget.second->branchX = tree->GetBranch(branchName);
    phidget.second->branchX->SetAddress(&phidget.second->data.Bx);

    snprintf(branchName, 64, PHIDGET_FORMAT_Y, phidget.first);
    phidget.second->branchY = nullptr;
    phidget.second->branchY = tree->GetBranch(branchName);
    phidget.second->branchY->SetAddress(&phidget.second->data.By);

    snprintf(branchName, 64, PHIDGET_FORMAT_Z, phidget.first);
    phidget.second->branchZ = nullptr;
    phidget.second->branchZ = tree->GetBranch(branchName);
    phidget.second->branchZ->SetAddress(&phidget.second->data.Bz);

    if (phidget.second->branchX == nullptr
        || phidget.second->branchY == nullptr
        || phidget.second->branchZ == nullptr) {
      return false;
    }
  }

  // Set gantry branches
  for (auto gantry : gantryData) {
    snprintf(branchName, 64, GANTRY_FORMAT_X, (int)gantry.second->gantry);
    gantry.second->branchX = nullptr;
    gantry.second->branchX = tree->GetBranch(branchName);
    gantry.second->branchX->SetAddress(&gantry.second->data.x);

    snprintf(branchName, 64, GANTRY_FORMAT_Y, (int)gantry.second->gantry);
    gantry.second->branchY = nullptr;
    gantry.second->branchY = tree->GetBranch(branchName);
    gantry.second->branchY->SetAddress(&gantry.second->data.y);

    snprintf(branchName, 64, GANTRY_FORMAT_Z, (int)gantry.second->gantry);
    gantry.second->branchZ = nullptr;
    gantry.second->branchZ = tree->GetBranch(branchName);
    gantry.second->branchZ->SetAddress(&gantry.second->data.z);

    snprintf(branchName, 64, GANTRY_FORMAT_THETA, (int)gantry.second->gantry);
    gantry.second->branchTheta = nullptr;
    gantry.second->branchTheta = tree->GetBranch(branchName);
    gantry.second->branchTheta->SetAddress(&gantry.second->data.theta);

    snprintf(branchName, 64, GANTRY_FORMAT_PHI, (int)gantry.second->gantry);
    gantry.second->branchPhi = nullptr;
    gantry.second->branchPhi = tree->GetBranch(branchName);
    gantry.second->branchPhi->SetAddress(&gantry.second->data.phi);

    if (gantry.second->branchX == nullptr
        || gantry.second->branchY == nullptr
        || gantry.second->branchZ == nullptr
        || gantry.second->branchTheta == nullptr
        || gantry.second->branchPhi == nullptr) {
      return false;
    }
  }

  TBranch* brNumSamples = tree->GetBranch("num_points");

  brNumSamples->SetAddress(&numSamples);

  return true;
}


bool Wrapper::unsetDataPointers() {
  if (tree == nullptr || file == nullptr)
    return false;
  
  for (auto pmt : pmtData) {
    pmt.second->branch = nullptr;
  }

  for (auto phidget : phidgetData) {
    phidget.second->branchX = nullptr;
    phidget.second->branchY= nullptr;
    phidget.second->branchZ = nullptr;
  }

  for (auto gantry : gantryData) {
    gantry.second->branchX = nullptr;
    gantry.second->branchY= nullptr;
    gantry.second->branchZ = nullptr;
    gantry.second->branchTheta = nullptr;
    gantry.second->branchPhi = nullptr;
  }

  return true;
}


/* Public functions */


void Wrapper::openFile(const string& fileName, const string& treeName) {
  file = new TFile(fileName.c_str(), "READ");

  if (!file->IsOpen()) {
    delete file;
    file = nullptr;
    throw new Exceptions::FileDoesNotExist(fileName);
  }

  tree = nullptr;
  file->GetObject(treeName.c_str(), tree);

  if (!tree) {
    throw new Exceptions::InvalidTreeName(treeName);
  }

  auto res = setDataPointers();

  if (!res) {
    throw new Exceptions::DataPointerError();
  }

  numEntries = tree->GetEntries();

  tree->GetEntry(0);
  entry = 0;
}


bool Wrapper::isFileOpen() const {
  return file && tree;
}


void Wrapper::closeFile() {
  if (tree) {
    unsetDataPointers();
    delete tree;
    tree = nullptr;
  }
  if (file) {
    file->Close();
    delete file;
    file = nullptr;
  }
  entry = UINT32_MAX;
}


// Open the settings tree to retrieve information:
PTF::brbReadings Wrapper::GetSettings(const string& fileName, const string& settingsTreeName) {

  // Declare a struct to contain the settings
  PTF::brbReadings brbSettings;

  // Define variables to store the data.
  int channelMask;        // "channel_mask/Int_t"
  double baseline[20];    // "CalcBaseline[20]/Double_t"
  double hvSetPoints[20]; // "HVsetpoints[20]/Double_t"
  double hvReadback[20];  // "HVreadback[20]/Double_t"
  double hvCurrent[20];   // "HVcurrent[20]/Double_t"

  

  // Open the root file if it is not open
  // file = new TFile(fileName.c_str(), "READ");

  // Get tree
  auto settingsTree = new TTree();
  file->GetObject(settingsTreeName.c_str(), settingsTree);

  // Get the branches
  auto dataMask = settingsTree->GetBranch("channel_mask");
  auto dataBaseline = settingsTree->GetBranch("CalcBaseline");
  auto dataHVsetPoints = settingsTree->GetBranch("HVsetpoints");
  auto dataHVreadBack = settingsTree->GetBranch("HVreadback");
  auto dataHVcurrent = settingsTree->GetBranch("HVcurrent");

  // Point data to address
  dataMask->SetAddress(&channelMask);
  dataBaseline->SetAddress(&baseline);
  dataHVsetPoints->SetAddress(&hvSetPoints);
  dataHVreadBack->SetAddress(&hvReadback);
  dataHVcurrent->SetAddress(&hvCurrent);

  // Get entry
  settingsTree->GetEntry(0);

  brbSettings.channelMask = channelMask; 
  for (int k =0; k < 20; k++) {
    brbSettings.baseline[k] = baseline[k];   
    brbSettings.hvSetPoints[k] = hvSetPoints[k];
    brbSettings.hvReadback[k] = hvReadback[k];
    brbSettings.hvCurrent[k] = hvCurrent[k]; 
  }

  // Return the settings
  return brbSettings;
}


int Wrapper::getChannelForPmt(int pmt) const {
  // Could be replaced with binary search, but probably list is small enough to not matter
  auto res = pmtData.find(pmt);

  if (res != pmtData.end()) {
    return res->second->channel;
  }
  else {
    return -1;
  }
}


int Wrapper::getPmtForChannel(int channel) const {
  for (auto pmt : pmtData) {
    if (pmt.second->channel == channel) {
      return pmt.first;
    }
  }
  return -1;
}


unsigned long long Wrapper::getCurrentEntry() const {
  if (!isFileOpen()) {
    throw new Exceptions::NoFileIsOpen();
  }
  return entry;
}


unsigned long long Wrapper::getNumEntries() const {
  if (!isFileOpen()) {
    throw new Exceptions::NoFileIsOpen();
  }
  return numEntries;
}


void Wrapper::setCurrentEntry(unsigned long long entry) {
  if (!isFileOpen()) {
    throw new Exceptions::NoFileIsOpen();
  }
  if (entry >= numEntries) {
    throw new Exceptions::EntryOutOfRange();
  }

  this->tree->GetEntry(entry);
  this->entry = entry;
}


unsigned long long Wrapper::getNumSamples() const {
  if (!isFileOpen()) {
    throw new Exceptions::NoFileIsOpen();
  }
  return numSamples;
}


double* Wrapper::getPmtSample(int pmt, unsigned long long sample) const {
  if (sample > numSamples) {
    throw new Exceptions::SampleOutOfRange();
  }
  auto pmtData = getDataForPmt(pmt);
  if (pmtData == nullptr) {
    throw new Exceptions::InvalidPMT();
  }
  return pmtData + (sample * sampleSize);
}


int Wrapper::getSampleLength() const {
  return sampleSize;
}


GantryData Wrapper::getDataForCurrentEntry(Gantry whichGantry) const {
  if (!isFileOpen()) {
    throw new Exceptions::NoFileIsOpen();
  }
  auto res = gantryData.find(whichGantry);
  if (res == gantryData.end()) {
    throw new Exceptions::InvalidGantry();
  }
  else {
    return res->second->data;
  }
}


PhidgetReading Wrapper::getReadingForPhidget(int phidget) const {
  if (!isFileOpen()) {
    throw new Exceptions::NoFileIsOpen();
  }
  auto res = phidgetData.find(phidget);
  if (res == phidgetData.end()) {
    throw new Exceptions::InvalidPhidget();
  }
  else {
    return res->second->data;
  }
}
