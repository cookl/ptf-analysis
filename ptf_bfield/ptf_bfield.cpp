/*
 * ptf_bfield.cpp
 * Program to calculate magnetic field around PMT test facility.
 *
 * Author: Blair Jamieson (Oct. 2019)
 */

#include <iostream>
#include <sstream>
#include <fstream>
#include <iomanip>
#include <vector>
#include <string>
#include <cmath>
#include "TStyle.h"
#include "TVector3.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TCanvas.h"
#include "TFile.h"
#include "TLine.h"
#include "t2kstyle.h"
#include "TGraph2D.h"

using std::cout;
using std::endl;

//const double    pi          = std::acos(-1.0);
const double   mu0_over_4pi = 1.0e-7;

// Earth's field in this problem is (303, 0, -366  ) mgauss
// where -y is vertical direction, +x is direction of B-field in horizontal plane
const TVector3 BEarth{ 0.0, 0.0, 0.0 }; // Tesla

// Use root TVector3 to represent vectors

// represent wire element by vector of length three
// all lenghts in meters!
struct WireElement{
  WireElement();
  // constructors take vector dl and location rp
  WireElement( TVector3 adl, TVector3 arp ) : dl( adl ), rp( arp ){;}

  // get back info about element
  TVector3 get_dl() const { return dl; }
  TVector3 get_rp() const { return rp; }
 
  // method to get dl cross ( r-r' )
  void get_dl_x_r( const TVector3 & ar, TVector3 & dlxr ) const;

  // method to get dl cross ( r-r' ) / |r-r'|^3
  void get_dl_x_r3( const TVector3 & ar, TVector3 & dlxr ) const;
 
private:
  TVector3 dl; // Direction
  TVector3 rp; // Location of current element  
};

typedef std::vector< WireElement > Wire;

// method to get dl cross ( r-r' )
void WireElement::get_dl_x_r( const TVector3 & ar, TVector3 & dlxr ) const{
  TVector3 rr =ar - rp;
  dlxr = dl.Cross( rr );
}

// method to get dl cross ( r-r' ) / |r-r'|^3
void WireElement::get_dl_x_r3( const TVector3 & ar, TVector3 & dlxr ) const{
    TVector3 rr =ar - rp;
    dlxr = dl.Cross( rr );
    dlxr *= 1.0/( pow( rr.Mag(), 3 ) );
}

// class to calculate magnetic field at a point.
// initialize class with:
//   1) current in amps
//   2) std::vector< WireElement > to represent locations and directions of current
// methods then exist to get field calculated at:
//   1) a point (TVector3)
//   2) a vector of points ( std::vector< TVector3 > )
//   3) as a TH1D along a line
//   4) as a TH2D on one of planes
struct BiotSavart{
  BiotSavart( double aI, const Wire& ww ) :
    I(aI), w( ww ) { ; }

  BiotSavart() : I(0.) { }

  // setters
  void set_I( double aI ){ I=aI; }
  void set_wire( const Wire& ww ){ w = ww; }

  // field at point
  TVector3 get_B( const TVector3 &p ) const;

  // field at vector of points
  std::vector< TVector3 > get_B( const std::vector< TVector3 > & p ) const;

  // Build and fill a 1D histogram of |B| along line
  // ceneter at point loc and going in direction dir (dir must be normalized)
  TH1D* get_B( const TVector3 &  loc, const TVector3 & dir,
	       const int nbinsx, const double xmin, const double xmax ) const;

  // Build and fill 2D histogram of B in plane
  // Define plane by a center location and two (presumed to be orthogonal) vectors as TVector3
  // assumes xdir and ydir are normalized vectors
  // Define region of histogram as coordinates relative to center location
  TH2D* get_B( const TVector3 & c, const TVector3 & xdir, const TVector3 & ydir,
	       const int nbinsx, const double xmin, const double xmax,
	       const int nbinsy, const double ymin, const double ymax ) const;

  const Wire & get_Wire() const { return w; }
  double get_I() const { return I; }
  
private:
  double I;
  Wire w;
};

// field at point
TVector3 BiotSavart::get_B( const TVector3& p ) const{
  TVector3 B(0.,0.,0.);
  TVector3 dB(0.,0.,0.);
  for ( auto elem : w ){
    elem.get_dl_x_r3( p, dB );
    B += dB;
  }
  B *= mu0_over_4pi * I;
  return B;
}

// field at vector of points
std::vector< TVector3 > BiotSavart::get_B( const std::vector< TVector3 > & vp ) const{
  std::vector< TVector3 > B;
  for ( const auto p : vp ){
    B.push_back( get_B( p ) );
  }
  return B;
}

// Build and fill a 1D histogram of B along line
// center at point loc and going in direction dir (dir must be normalized)
TH1D* BiotSavart::get_B( const TVector3 &  loc, const TVector3 & dir,
			 const int nbinsx, const double xmin, const double xmax ) const {
  static int count=0;
  count++;
  std::string hname{"get_B1D"};
  hname += std::to_string( count ); 
  TH1D* h = new TH1D( hname.c_str(), " ; loc (m) ; B ( Tesla )", nbinsx, xmin, xmax );
  for ( int i=1; i<=nbinsx; ++i){
    double bc = h->GetBinCenter( i );
    TVector3 p = loc + bc*dir;
    TVector3 B = get_B( p );
    h->SetBinContent( i, B.Mag() );
  }
  return h;
}

// Build and fill 2D histogram of B in plane
// Define plane by a center location and two (presumed to be orthogonal) vectors as TVector3
// assumes xdir and ydir are normalized vectors
// Define region of histogram as coordinates relative to center location
TH2D* BiotSavart::get_B( const TVector3 & c, const TVector3 & xdir, const TVector3 & ydir,
			 const int nbinsx, const double xmin, const double xmax,
			 const int nbinsy, const double ymin, const double ymax ) const {
  static int count=0;
  count++;
  std::string hname{"get_B2D"};
  hname += std::to_string( count ); 
  TH2D* h = new TH2D( hname.c_str(), " ; loc (m) ; loc (m); B ( Tesla )",
		      nbinsx, xmin, xmax, nbinsy, ymin, ymax );
  TAxis * xax = h->GetXaxis();
  TAxis * yax = h->GetYaxis();
  for ( int i=1; i<=nbinsx; ++i){
    for ( int j=1; j<=nbinsy; ++j){
      double x = xax->GetBinCenter( i );
      double y = yax->GetBinCenter( j );
      TVector3 p = c + x * xdir + y * ydir;
      TVector3 B = get_B( p );
      h->SetBinContent( i, j, B.Mag() );
    }
  }
  return h;
}


// PTF coil dimensions
const double   x_length = 1.96; // meters
const double   y_length = 2.04; // meters
const double   z_length = 1.89; // meters

// Coil currents:
const int      num_coils = 6;
const int      coil_turns[ num_coils ]    = {  52,  52,   52,   52,  260,  260 };
const double   coil_voltages[ num_coils ] = { 6.6, 9.3, 2.47, 1.58, 1.74, 1.26 }; 
const double   coil_resistances[ num_coils]= { 3.8, 3.8,  3.8,  3.9, 11.1, 11.1 };


/// PTFCoils
/// Builds 6 coils with dimensions of length x_length, y_length, z_length
/// Spacing of the coils is base on Helmholtz geometry (spaced by 1/2 length)
/// Can independently set current of each coil.
/// Default currents are from coil_turns * coil_voltages / coil_resistances
///
/// Numbering for the current loops is:
/// 1) Bottom coil compensating in z direction (vertical)
/// 2) Top coil compensating in z direction
/// 3) Left coil compensating in x direction (toward M11)
/// 4) Right coil compensating in x direction (away from M11)
/// 5) Front coil compensating in y direction (further into room)
/// 6) Back coil compensating in y direction (closer to door)
///
struct PTFCoils{
  PTFCoils();   
  // setters
  void        set_voltage( int iwire, double V );

  // getters
  const Wire& get_wire_loop( int iwire ) const { return bs[iwire].get_Wire(); }
  double      get_current( int iwire ) const { return bs[iwire].get_I(); }
  TVector3    get_B( const TVector3 & p ) const;

  friend void plot_coils( const PTFCoils& ptfc );  
protected:

  BiotSavart  bs[num_coils];
};


std::ostream& operator<<( std::ostream& os, const PTFCoils& ptfc ){
  for (int i=0; i<num_coils ; ++i ) {
    os << "Coil " << i
       <<" Current = " << ptfc.get_current( i )
       <<" turns=" << coil_turns[i]
       << " Current per turn = "<<ptfc.get_current( i ) / coil_turns[i]
       << endl;
  }
  return os;
}


/// addtowire
/// Helper function to add to an existing Wire w
/// Adds wire elements in steps from from to to in npts steps
void addtowire( Wire& w, const TVector3& from, const TVector3& to, int npts ){
  TVector3 dlvec = to - from;
  double length = dlvec.Mag();
  double dl = length / npts;
  dlvec *= dl/length; // dl scaled to step size
  TVector3 steploc = from - 0.5*dlvec;
  for ( int step=0; step<npts; ++step ){
    steploc += dlvec;
    w.push_back( WireElement( dlvec, steploc ) );
  }
}
  

PTFCoils::PTFCoils(){
  // build the coils
  const int ndl = 100;
  Wire coils[num_coils];

  // coil 1 (a z-coil so lengths of wires are all z_length):
  //         at z = -z_length/4 
  //         in an offset xy-plane
  // wire at +x (going in +y) 
  {
    TVector3 from( +z_length/2, -z_length/2, -z_length/4 );
    TVector3 to  ( +z_length/2, +z_length/2, -z_length/4 );
    addtowire( coils[0], from, to, ndl );
  }
  // wire at +y (going in -x)
  {
    TVector3 from( +z_length/2, +z_length/2, -z_length/4 );
    TVector3 to  ( -z_length/2, +z_length/2, -z_length/4 );
    addtowire( coils[0], from, to, ndl );
  }
  // wire at -x (going in -y) 
  {
    TVector3 from( -z_length/2, +z_length/2, -z_length/4 );
    TVector3 to  ( -z_length/2, -z_length/2, -z_length/4 );
    addtowire( coils[0], from, to, ndl );
  }
  // wire at -y (going in +x)
  {
    TVector3 from( -z_length/2, -z_length/2, -z_length/4 );
    TVector3 to  ( +z_length/2, -z_length/2, -z_length/4 );
    addtowire( coils[0], from, to, ndl );
  }


  // coil 2 (a z-coil so lengths of wires are all z_length):
  //         at z = +z_length/4 
  //         in an offset xy-plane
  // wire at +x (going in +y) 
  {
    TVector3 from( +z_length/2, -z_length/2, +z_length/4 );
    TVector3 to  ( +z_length/2, +z_length/2, +z_length/4 );
    addtowire( coils[1], from, to, ndl );
  }
  // wire at +y (going in -x)
  {
    TVector3 from( +z_length/2, +z_length/2, +z_length/4 );
    TVector3 to  ( -z_length/2, +z_length/2, +z_length/4 );
    addtowire( coils[1], from, to, ndl );
  }
  // wire at -x (going in -y) 
  {
    TVector3 from( -z_length/2, +z_length/2, +z_length/4 );
    TVector3 to  ( -z_length/2, -z_length/2, +z_length/4 );
    addtowire( coils[1], from, to, ndl );
  }
  // wire at -y (going in +x)
  {
    TVector3 from( -z_length/2, -z_length/2, +z_length/4 );
    TVector3 to  ( +z_length/2, -z_length/2, +z_length/4 );
    addtowire( coils[1], from, to, ndl );
  }


  // coil 3 (a x-coil so lengths of wires are all x_length):
  //         at x = -x_length/4 
  //         in an offset yz-plane
  // wire at +y (going in +z) 
  {
    TVector3 from( -x_length/4, +x_length/2, -x_length/2 );
    TVector3 to  ( -x_length/4, +x_length/2, +x_length/2 );
    addtowire( coils[2], from, to, ndl );
  }
  // wire at +z (going in -y)
  {
    TVector3 from( -x_length/4, +x_length/2, +x_length/2 );
    TVector3 to  ( -x_length/4, -x_length/2, +x_length/2 );
    addtowire( coils[2], from, to, ndl );
  }
  // wire at -y (going in -z) 
  {
    TVector3 from( -x_length/4, -x_length/2, +x_length/2 );
    TVector3 to  ( -x_length/4, -x_length/2, -x_length/2 );
    addtowire( coils[2], from, to, ndl );
  }
  // wire at -z (going in +y)
  {
    TVector3 from( -x_length/4, -x_length/2, -x_length/2 );
    TVector3 to  ( -x_length/4, +x_length/2, -x_length/2 );
    addtowire( coils[2], from, to, ndl );
  }

  // coil 4 (a x-coil so lengths of wires are all x_length):
  //         at x = +x_length/4 
  //         in an offset yz-plane
  // wire at +y (going in +z) 
  {
    TVector3 from( +x_length/4, +x_length/2, -x_length/2 );
    TVector3 to  ( +x_length/4, +x_length/2, +x_length/2 );
    addtowire( coils[3], from, to, ndl );
  }
  // wire at +z (going in -y)
  {
    TVector3 from( +x_length/4, +x_length/2, +x_length/2 );
    TVector3 to  ( +x_length/4, -x_length/2, +x_length/2 );
    addtowire( coils[3], from, to, ndl );
  }
  // wire at -y (going in -z) 
  {
    TVector3 from( +x_length/4, -x_length/2, +x_length/2 );
    TVector3 to  ( +x_length/4, -x_length/2, -x_length/2 );
    addtowire( coils[3], from, to, ndl );
  }
  // wire at -z (going in +y)
  {
    TVector3 from( +x_length/4, -x_length/2, -x_length/2 );
    TVector3 to  ( +x_length/4, +x_length/2, -x_length/2 );
    addtowire( coils[3], from, to, ndl );
  }

  // coil 5 (a y-coil so lengths of wires are all y_length):
  //         at y = -y_length/4 
  //         in an offset zx-plane
  // wire at +z (going in +x) 
  {
    TVector3 from( -y_length/2, -y_length/4, +y_length/2 );
    TVector3 to  ( +y_length/2, -y_length/4, +y_length/2 );
    addtowire( coils[4], from, to, ndl );
  }
  // wire at +x (going in -z)
  {
    TVector3 from( +y_length/2, -y_length/4, +y_length/2 );
    TVector3 to  ( +y_length/2, -y_length/4, -y_length/2 );
    addtowire( coils[4], from, to, ndl );
  }
  // wire at -z (going in -x) 
  {
    TVector3 from( +y_length/2, -y_length/4, -y_length/2 );
    TVector3 to  ( -y_length/2, -y_length/4, -y_length/2 );
    addtowire( coils[4], from, to, ndl );
  }
  // wire at -x (going in +z)
  {
    TVector3 from( -y_length/2, -y_length/4, -y_length/2 );
    TVector3 to  ( -y_length/2, -y_length/4, +y_length/2 );
    addtowire( coils[4], from, to, ndl );
  }

  // coil 6 (a y-coil so lengths of wires are all y_length):
  //         at y = +y_length/4 
  //         in an offset zx-plane
  // wire at +z (going in +x) 
  {
    TVector3 from( -y_length/2, +y_length/4, +y_length/2 );
    TVector3 to  ( +y_length/2, +y_length/4, +y_length/2 );
    addtowire( coils[5], from, to, ndl );
  }
  // wire at +x (going in -z)
  {
    TVector3 from( +y_length/2, +y_length/4, +y_length/2 );
    TVector3 to  ( +y_length/2, +y_length/4, -y_length/2 );
    addtowire( coils[5], from, to, ndl );
  }
  // wire at -z (going in -x) 
  {
    TVector3 from( +y_length/2, +y_length/4, -y_length/2 );
    TVector3 to  ( -y_length/2, +y_length/4, -y_length/2 );
    addtowire( coils[5], from, to, ndl );
  }
  // wire at -x (going in +z)
  {
    TVector3 from( -y_length/2, +y_length/4, -y_length/2 );
    TVector3 to  ( -y_length/2, +y_length/4, +y_length/2 );
    addtowire( coils[5], from, to, ndl );
  }

  /// all the wires are made
  /// now setup BiotSavart objects
  for (int i=0; i< num_coils; ++i ){
    double current = coil_voltages[i] / coil_resistances[i] * coil_turns[i];
    bs[i].set_I( current );
    bs[i].set_wire( coils[i] );
  }
}

void PTFCoils::set_voltage( int iwire, double V){
  double current = V / coil_resistances[iwire] * coil_turns[iwire];
  bs[iwire].set_I( current );
}

TVector3 PTFCoils::get_B( const TVector3 & p ) const{
  TVector3 btot(0.,0.,0.);
  for (int icoil=0; icoil<num_coils; ++icoil){
    btot += bs[icoil].get_B( p );
  }
  return btot;
}


void plot_coils( const PTFCoils& ptfc ) {
  std::vector<double> x,y,z;
  for (int i=0; i<num_coils; ++i ){
    const Wire& w = ptfc.get_wire_loop( i );
    for ( const WireElement& we : w ){
      TVector3 loc = we.get_rp();
      x.push_back( loc.X() );
      y.push_back( loc.Y() );
      z.push_back( loc.Z() );
    }
  }
  TGraph2D * tg2d = new TGraph2D( x.size(), &x[0], &y[0], &z[0] );
  tg2d->SetName("tg_ptfcoilloc");
  tg2d->SetTitle("PTF Coil Locations; X(m); Y(m); Z(m)");
  tg2d->Write();
}


void biotsavart(){
  t2kstyle();
  gStyle->SetPadRightMargin(0.2);
  gStyle->SetPadLeftMargin(0.3);
  gStyle->SetPadBottomMargin(0.2);

  
  // Output file
  TFile * fout = new TFile("ptf_bfield.root","recreate");

  // Build PTF coils
  PTFCoils ptfc;

  // Change voltages on coils
  if (0){
    ptfc.set_voltage( 0, 5.0  );  // plot per 5 amp
    ptfc.set_voltage( 1, 5.0  );  // plot per 5 amp
    ptfc.set_voltage( 2, 5.0  );  // plot per 5 amp
    ptfc.set_voltage( 3, 5.0  );  // plot per 5 amp
    ptfc.set_voltage( 4, 5.0  );  // plot per 5 amp
    ptfc.set_voltage( 5, 5.0  );  // plot per 5 amp
  }
  
  plot_coils( ptfc ); // plot locations of coils

  cout << ptfc << endl;

  // Make histograms of magnetic field
  std::string tag = "ptfc";
  std::ostringstream os;

  // make histograms of magnetic field from PTF coils
  os.str(""); os.clear(); os << tag << "_hxyBx";
  TH2D* hxyBx = new TH2D(os.str().c_str()," Bx (gauss); x (m); y (m) ",100,-x_length/2, x_length/2, 100, -y_length/2, y_length/2 );
  os.str(""); os.clear(); os << tag << "_hxyBy";
  TH2D* hxyBy = new TH2D(os.str().c_str()," By (gauss); x (m); y (m) ",100,-x_length/2, x_length/2, 100, -y_length/2, y_length/2 );
  os.str(""); os.clear(); os << tag << "_hxyBz";
  TH2D* hxyBz = new TH2D(os.str().c_str()," Bz (gauss); x (m); y (m) ",100,-x_length/2, x_length/2, 100, -y_length/2, y_length/2 );
  os.str(""); os.clear(); os << tag << "_hxyB";
  TH2D* hxyB = new TH2D(os.str().c_str()," B (gauss); x (m); y (m) "  ,100,-x_length/2, x_length/2, 100, -y_length/2, y_length/2 );

  for ( int ix=1; ix<=100; ++ix ){
    for ( int iy=1; iy<=100; ++iy ){
      double x = hxyBx->GetXaxis()->GetBinCenter( ix );
      double y = hxyBx->GetYaxis()->GetBinCenter( iy );
      TVector3 loc = TVector3{ x, y, 0.0 } ;
      TVector3 btot = ptfc.get_B( loc ) + BEarth;
      btot *= 1.0e4;//convert to gauss
      hxyBx->SetBinContent( ix, iy, btot.X() );
      hxyBy->SetBinContent( ix, iy, btot.Y() );
      hxyBz->SetBinContent( ix, iy, btot.Z() );
      hxyB->SetBinContent(  ix, iy, btot.Mag() );
    }
  }

  hxyBx->SetMinimum( -2.0 ); hxyBx->UseCurrentStyle();
  hxyBx->SetMaximum(  2.0 ); hxyBx->UseCurrentStyle();
  hxyBy->SetMinimum( -2.0 ); hxyBy->UseCurrentStyle();
  hxyBy->SetMaximum(  2.0 ); hxyBy->UseCurrentStyle();
  hxyBz->SetMinimum( -2.0 ); hxyBz->UseCurrentStyle();
  hxyBz->SetMaximum(  2.0 ); hxyBz->UseCurrentStyle();
  hxyB ->SetMinimum(  0.0 ); hxyB ->UseCurrentStyle();
  hxyB ->SetMaximum(  2.0 ); hxyB ->UseCurrentStyle();
  
  TCanvas * cxy = new TCanvas("cxy","cxy",800,600);
  cxy->Divide(2,2);
  cxy->cd(1);  hxyBx->Draw("colz");
  cxy->cd(2);  hxyBy->Draw("colz");
  cxy->cd(3);  hxyBz->Draw("colz");
  cxy->cd(4);  hxyB->Draw("colz");
  cxy->Write();
  
  os.str(""); os.clear(); os << tag << "_hyzBx";
  TH2D* hyzBx = new TH2D(os.str().c_str()," Bx (gauss); y (m); z (m) ",100, -y_length/2, y_length/2, 100, -z_length/2, z_length/2 );
  os.str(""); os.clear(); os << tag << "_hyzBy";
  TH2D* hyzBy = new TH2D(os.str().c_str()," By (gauss); y (m); z (m) ",100, -y_length/2, y_length/2, 100, -z_length/2, z_length/2 );
  os.str(""); os.clear(); os << tag << "_hyzBz";
  TH2D* hyzBz = new TH2D(os.str().c_str()," Bz (gauss); y (m); z (m) ",100, -y_length/2, y_length/2, 100, -z_length/2, z_length/2 );
  os.str(""); os.clear(); os << tag << "_hyzB";
  TH2D* hyzB = new TH2D(os.str().c_str()," B (gauss); y (m); z (m) "  ,100, -y_length/2, y_length/2, 100, -z_length/2, z_length/2 );

  for ( int iy=1; iy<=100; ++iy ){
    for ( int iz=1; iz<=100; ++iz ){
      double y = hyzBx->GetXaxis()->GetBinCenter( iy );
      double z = hyzBx->GetYaxis()->GetBinCenter( iz );
      TVector3 loc = TVector3{ 0.0, y, z } ;
      TVector3 btot = ptfc.get_B( loc ) + BEarth;
      btot *= 1.0e4;//convert to gauss

      hyzBx->SetBinContent( iy, iz, btot.X() );
      hyzBy->SetBinContent( iy, iz, btot.Y() );
      hyzBz->SetBinContent( iy, iz, btot.Z() );
      hyzB->SetBinContent( iy, iz, btot.Mag() );
    }
  }


  hyzBx->SetMinimum( -2.0 ); hyzBx->UseCurrentStyle();
  hyzBx->SetMaximum(  2.0 ); hyzBx->UseCurrentStyle();
  hyzBy->SetMinimum( -2.0 ); hyzBy->UseCurrentStyle();
  hyzBy->SetMaximum(  2.0 ); hyzBy->UseCurrentStyle();
  hyzBz->SetMinimum( -2.0 ); hyzBz->UseCurrentStyle();
  hyzBz->SetMaximum(  2.0 ); hyzBz->UseCurrentStyle();
  hyzB ->SetMinimum(  0.0 ); hyzB ->UseCurrentStyle();
  hyzB ->SetMaximum(  2.0 ); hyzB ->UseCurrentStyle();
  
  TCanvas * cyz = new TCanvas("cyz","cyz",800,600);
  cyz->Divide(2,2);
  cyz->cd(1);  hyzBx->Draw("colz");
  cyz->cd(2);  hyzBy->Draw("colz");
  cyz->cd(3);  hyzBz->Draw("colz");
  cyz->cd(4);  hyzB->Draw("colz");
  cyz->Write();
  
  os.str(""); os.clear(); os << tag << "_hzxBx";
  TH2D* hzxBx = new TH2D(os.str().c_str()," Bx (gauss); z (m); x (m) ",100, -z_length/2, z_length/2, 100, -x_length/2, x_length/2 );
  os.str(""); os.clear(); os << tag << "_hzxBy";
  TH2D* hzxBy = new TH2D(os.str().c_str()," By (gauss); z (m); x (m) ",100, -z_length/2, z_length/2, 100, -x_length/2, x_length/2 );
  os.str(""); os.clear(); os << tag << "_hzxBz";
  TH2D* hzxBz = new TH2D(os.str().c_str()," Bz (gauss); z (m); x (m) ",100, -z_length/2, z_length/2, 100, -x_length/2, x_length/2 );
  os.str(""); os.clear(); os << tag << "_hzxB";
  TH2D* hzxB = new TH2D(os.str().c_str()," B (gauss); z (m); x (m) "  ,100, -z_length/2, z_length/2, 100, -x_length/2, x_length/2 );

  for ( int iz=1; iz<=100; ++iz ){
    for ( int ix=1; ix<=100; ++ix ){
      double z = hzxBx->GetXaxis()->GetBinCenter( iz );
      double x = hzxBx->GetYaxis()->GetBinCenter( ix );
      TVector3 loc = TVector3{ x, 0, z } ;
      TVector3 btot = ptfc.get_B( loc ) + BEarth;
      btot *= 1.0e4;//convert to gauss

      hzxBx->SetBinContent( iz, ix, btot.X() );
      hzxBy->SetBinContent( iz, ix, btot.Y() );
      hzxBz->SetBinContent( iz, ix, btot.Z() );
      hzxB->SetBinContent( iz, ix, btot.Mag() );
    }
  }


  hzxBx->SetMinimum( -2.0 ); hzxBx->UseCurrentStyle();
  hzxBx->SetMaximum(  2.0 ); hzxBx->UseCurrentStyle();
  hzxBy->SetMinimum( -2.0 ); hzxBy->UseCurrentStyle();
  hzxBy->SetMaximum(  2.0 ); hzxBy->UseCurrentStyle();
  hzxBz->SetMinimum( -2.0 ); hzxBz->UseCurrentStyle();
  hzxBz->SetMaximum(  2.0 ); hzxBz->UseCurrentStyle();
  hzxB ->SetMinimum(  0.0 ); hzxB ->UseCurrentStyle();
  hzxB ->SetMaximum(  2.0 ); hzxB ->UseCurrentStyle();
  
  TCanvas * czx = new TCanvas("czx","czx",800,600);
  czx->Divide(2,2);
  czx->cd(1);  hzxBx->Draw("colz");
  czx->cd(2);  hzxBy->Draw("colz");
  czx->cd(3);  hzxBz->Draw("colz");
  czx->cd(4);  hzxB->Draw("colz");
  czx->Write();
  
  
  fout->Write();
  fout->Close();
}


int main(){

  biotsavart();
  

  return 0;
}
  

