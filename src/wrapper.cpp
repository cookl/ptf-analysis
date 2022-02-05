#include "wrapper.hpp"
#include "BrbSettingsTree.hxx"

using namespace std;
using namespace PTF;


//constructor of the class


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
  if (tree == nullptr || file == nullptr){
    cout << "false tree or file ppointer" << endl;
    return false;
}
  
  // Set PMT branches
  char branchName[64];
  for (auto pmt : pmtData) {
    snprintf(branchName, 64, PMT_CHANNEL_FORMAT, pmt.second->channel);
    pmt.second->branch = nullptr;
    pmt.second->branch = tree->GetBranch(branchName);
    if (pmt.second->branch == nullptr) {
      cout << "False second branch pointer " << branchName << endl;   
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
	
    snprintf(branchName, 64, PHIDGET_FORMAT_ACCX, phidget.first);
    phidget.second->branchX = nullptr;
    phidget.second->branchX = tree->GetBranch(branchName);
    phidget.second->branchX->SetAddress(&phidget.second->data.Ax);

    snprintf(branchName, 64, PHIDGET_FORMAT_ACCY, phidget.first);
    phidget.second->branchY = nullptr;
    phidget.second->branchY = tree->GetBranch(branchName);
    phidget.second->branchY->SetAddress(&phidget.second->data.Ay);

    snprintf(branchName, 64, PHIDGET_FORMAT_ACCZ, phidget.first);
    phidget.second->branchZ = nullptr;
    phidget.second->branchZ = tree->GetBranch(branchName);
    phidget.second->branchZ->SetAddress(&phidget.second->data.Az);

    if (phidget.second->branchX == nullptr
        || phidget.second->branchY == nullptr
        || phidget.second->branchZ == nullptr) {
      cout << "False branch xyz" << endl; 
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
  TBranch
    //*T_int = tree->GetBranch("int_temp"),//, *T_ext1 = tree->GetBranch("ext1_temp")
    *T_ext2 = tree->GetBranch("ext2_temp");
    
   // *braNumSamples = tree->GetBranch("num_points");
  //T_int->SetAddress(&Temp.int_1);
  //T_ext1->SetAddress(&Temp.ext_1);

  // Make sure this branch exists first
  if(T_ext2) T_ext2->SetAddress(&Temp.ext_2);


  //braNumSamples->SetAddress(&numSamples);
  TBranch
    *Time_1=tree->GetBranch("timestamp");

  // Make sure this branch exists first
  if(Time_1) Time_1->SetAddress(&ti.time_c);
	
   // TBranch
   //   *ACC_x= tree->GetBranch("gantry0_x"), *g0Y = tree->GetBranch("gantry0_y"), *g0Z = tree->GetBranch("gantry0_z"),
   //     *ACC_y = tree->GetBranch("gantry0_rot"), *g0Phi = tree->GetBranch("gantry0_tilt"),
   //   *ACC_z = tree->GetBranch("gantry1_x"), *g1Y = tree->GetBranch("gantry1_y"), *g1Z = tree->GetBranch("gantry1_z"),
    //    *g1Theta = tree->GetBranch("gantry1_rot"), *g1Phi = tree->GetBranch("gantry1_tilt");

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
    std::cout << "Error, ttreename not valid: " << treeName.c_str() << std::endl;
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


int Wrapper::LoadBrbSettingsTree(){

  return BrbSettingsTree::Get()->LoadSettingsTree(file);
  
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

Temperature_r Wrapper::getReadingTemperature() const {
  if (!isFileOpen()) {
    throw new Exceptions::NoFileIsOpen();
  }
  return Temp;
}
Timing Wrapper::getReadingTime() const {
 if (!isFileOpen()) {
    throw new Exceptions::NoFileIsOpen();
  }
  return ti;

}

 
