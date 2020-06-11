/// A simple class to do Hough transform on a circle
/// Author:  Blair Jamieson (Oct. 2019)
#ifndef __HOUGH__
#define __HOUGH__

#include "XYPoint.hpp"

#include <vector>
#include <utility>

#include "TH2D.h"
#include "TDirectory.h"

enum   HoughResultType {HoughCircle, HoughUnusedPoints };

/// Define structure type to store results in
/// A single of Hough transform is a vector of data points
struct HoughResult {
  HoughResult( float r, xypoint xy, float pkval ) :
    rc( r ), xyc( xy ), peakval(pkval), type(HoughUnusedPoints) { }
  
  std::vector< xypoint > data;      // the points for this result
  float                  rc;        // circle radius
  xypoint                xyc;       // circle center
  float                  peakval;   // peak height in hough space
  HoughResultType        type;      // result type
  unsigned               rbin{0};   // radial bin number
};
/// Vector of Hough transform results
typedef std::vector< HoughResult >  HoughResults ; 

/// number of circles found
unsigned num_circles( const HoughResults & hrs );     

/// print results
std::ostream& operator<<( std::ostream& os, const HoughResults& hrs );
std::ostream& operator<<( std::ostream& os, const HoughResult& hr );


/// A simple class to do Hough transform on a circle
/// Example usage, using default radius and bins range.
/// Data is provided as a std::vector< xypoint > data
///
/// CircleHough h;
/// HoughResults res = h.find_circles( data );
///
/// cout<<"Number of circles found is "<<res.num_circles()<<endl;
/// // Loop over circles and do a fit to the hits
/// for (int ires=0; ires<res.data.size(); ++ires ){
///     if ( res[ires].type == HoughCircle ){
///         circle_fit( res[ires].data );
///
/// Parameters that can be tweaked to tune the transform:
///   minhits:   Minimum number of hits to keep searching for circles
///              (default value 3)
///   threshold: Minimum peak height in transform space to call a result
///              (default value 10)
class CircleHough {
 public:
  /// Constructor sets parameters of the Hough Transform
  /// can set number of bins in radius, center x, and center y
  /// along with range in each of those values
  //  CircleHough( unsigned nbins_radius =  33, double rmin=3.0, double  rmax=99.0,
  CircleHough( unsigned nbins_radius =  50, double rmin=10.0, double  rmax=60.0,
	       unsigned nbins_x      = 166, double xmin=-249, double xmax=249,
	       unsigned nbins_y      = 166, double ymin=-249, double ymax=249);

  ~CircleHough();

  /// Tunable parameters
  /// Minimum hits to be called a circle
  /// default is 3 hits
  void set_minhits( int n ){ minhits = n; }
  /// Minimum height in hough space to be called a circle
  /// default is 5
  void set_threshold( int n ){ threshold = n; }
  /// Scaling factor for distance from circle for hit to be part of circle
  /// Value of 1.0 is default, uses sqrt(dx*dx+dy*dy) where dx, dy are
  /// the xc,yc bin sizes.
  void set_distance_factor( float dist_factor ){ drscaling=dist_factor; }
  
  /// Main method to find circles
  const HoughResults& find_circles( const std::vector< xypoint >& data );

  /// Get the histogram of the transformed data
  std::vector< TH2D* > get_transform() { return fTransformed; }
  std::vector< std::pair< double, double > > get_rbins() { return fRbins; }

 private:
  /// Parameters
  int minhits{3};
  int threshold{5};
  float drscaling{1.0};
  
  /// Store results
  HoughResults fresults;

  // Histogram used to store the transformed data
  // bin ranges as pairs of doubles for each bin
  std::vector< std::pair< double, double > > fRbins;
  // one XY histogram per bin
  std::vector< TH2D* > fTransformed;

  // current directory to store results
  TDirectory * houghdir;
  
  // helper methods

  // Fill a HoughResult, pass it a copy of all the hits that
  // and any that don't get used in the HoughResult are returned
  // as unused hits
  HoughResult find_maximum( std::vector< xypoint >& unused_hits );

  // Reset the hough transformed histograms and fill them
  // with the data passed
  void hough_transform( const std::vector< xypoint >& data );

  // Make a clone of one of hough transformed histograms for
  // a candidate circle.  Pass it the candidate number, and
  // histogram pointer
  int nfind_circles{0};//call count for find_circles
  void save_hough_histo( unsigned num, TH2D* histo );

  // make a histogram of the candidate circle, and draw
  // the hough estimate for radius and center
  void plot_candidate( unsigned num, const HoughResult & hr );
}; /// end of CircleHough



#endif // __HOUGH__
