#include "ScanPoint.hpp"
#include "TTree.h"
#include "TBranch.h"

ScanPoint::ScanPoint(){
  fX=0.; fY=0.; fZ=0.;
  fEntry=0;
  fEntries=0; // haven't set an entry yet
}

ScanPoint::ScanPoint( double x, double y, double z, unsigned long long entry, unsigned long long entries ){
  fX = x;  fY = y;  fZ = z;
  fEntry = entry;
  fEntries = entries;
}

void WriteScanPoints( const std::vector< ScanPoint > & scanpoints ){
  float X,Y,Z;
  unsigned long long Entry,Entries;

  TTree * tt = new TTree( "scanpoints", "scanpoints" );
  tt->Branch( "X",       &X,       "X/F" );
  tt->Branch( "Y",       &Y,       "Y/F" );
  tt->Branch( "Z",       &Z,       "Z/F" );
  tt->Branch( "Entry",   &Entry,   "Entry/l" );
  tt->Branch( "Entries", &Entries, "Entries/l" );

  for ( const ScanPoint& sp : scanpoints ){
    X=sp.x();
    Y=sp.y();
    Z=sp.z();
    Entry=sp.get_entry();
    Entries=sp.nentries();
    std::cout<<"Filling TTree with "<<sp<<std::endl;
    tt->Fill();
  }
  tt->Write();
}

std::vector< ScanPoint > ReadScanPoints( TFile * fin ){
  std::vector< ScanPoint > result;
  float X,Y,Z;
  unsigned long long Entry,Entries;

  TTree * tt = (TTree*)fin->Get("scanpoints");
  if ( tt == nullptr ) {
    std::cerr<<"ReadScanPoints: could not find TTree named scanpoints"<<std::endl;
    return result;
  }
  tt->SetBranchAddress( "X", &X );
  tt->SetBranchAddress( "Y", &Y );
  tt->SetBranchAddress( "Z", &Z );
  tt->SetBranchAddress( "Entry", &Entry );
  tt->SetBranchAddress( "Entries", &Entries );

  unsigned long long n = tt->GetEntries();
  for ( unsigned long long i = 0 ; i < n ; ++i ){
    tt->GetEvent( i );
    result.push_back( ScanPoint( X, Y, Z, Entry, Entries ) );
  }
  return result;
}

std::ostream& operator<<( std::ostream& os, const ScanPoint& sp ){
  os << "(X,Y,Z)=( "<<sp.x()<<", "<<sp.y()<<", "<<sp.z()
     <<") First entry="<<sp.get_entry()
     <<" Entries ="<<sp.nentries()
     << std::endl;
  return os;
}
