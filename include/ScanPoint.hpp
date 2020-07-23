#ifndef __SCANPOINT__
#define __SCANPOINT__

#include <vector>
#include <iostream>
#include "TFile.h"

/// Holds x,y,z of scan point
/// First TTree entry number of this scan point
/// Number of waveforms at this scan point
class ScanPoint {
public:
  ScanPoint(); // default constructor
  void set_xyz( double x, double y, double z,double time_1,double t_ext2 ){ fX=x; fY=y; fZ=z;fT_ext2=t_ext2;
    fTime_1=time_1; }
  void set_entry( unsigned long long entry ){ fEntry = entry; ++fEntries; }
  
  ScanPoint( double x, double y, double z, double time_1,double t_ext2,
	     unsigned long long entry, unsigned long long entries = 0 );

  ScanPoint & operator++(){ ++fEntries; return *this; }// prefix increment operator
  
  void get_xyz( double& xx, double &yy, double&zz, double&tti,double&tt_ext2 ) const { xx=fX; yy=fY; zz=fZ;tti=fTime_1;tt_ext2=fT_ext2; }
  double x() const { return fX; }
  double y() const { return fY; }
  double z() const { return fZ; }
  double time_1() const {return fTime_1;}
  //double t_int1()	const	{return fT_int1;}
  //double t_ext1()	const	{return fT_ext1;}
  //double t_ext2()	const	{return fT_ext2;}
  unsigned long long get_entry() const { return fEntry; }
  unsigned long long nentries() const { return fEntries; }

private:
  double fX, fY, fZ,fTime_1,fT_ext2;            // x, y, z location of scan point
  unsigned long long fEntry;    // first entry in TTree of this scan pt
  unsigned long long fEntries;  // number of entries (waveforms) in TTree for this scan pt
};


// Write a vector of ScanPoint to TTree
void WriteScanPoints( const std::vector< ScanPoint >& scanpoints );

// Read a TTree data into vector of ScanPoints
std::vector< ScanPoint > ReadScanPoints( TFile * fin );

/// Print scanpoint
std::ostream& operator<<( std::ostream& os, const ScanPoint& sp );

#endif // __SCANPOINT__
