#ifndef __ODB_LAYER
#define __ODB_LAYER

#include "midas.h"

#include <unordered_map>
#include <vector>

#include "degauss.hxx"


struct PhidgetReading {
  double acc0;
  double acc1;
  double acc2;
  double b0;
  double b1;
  double b2;
  double mag;
  double tilt;
  double etime;
  double uetime;

  static PhidgetReading all_nan();
};


struct ODBInterface {
  ODBInterface(const vector<uint32_t>& phidgets);
  ~ODBInterface();

  double get_coil_voltage(Coil coil);
  void   set_coil_voltage(Coil coil, double voltage);

  PhidgetReading get_phidget_reading(uint32_t phidget);

private:

  // stuff for connecting to ODB here
  HNDLE
    hDB, hkeyclient,
    coil_v_read, coil_v_write, coil_u_read;
  char host_name[256], exp_name[256];

  // cache the ODB keys for phidgets
  std::unordered_map<uint32_t, HNDLE> phidget_handle_cache;
};




#endif 
