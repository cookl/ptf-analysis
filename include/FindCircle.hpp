#ifndef __FINDCIRCLE__
#define __FINDCIRCLE__

#include "TH2D.h"

/// Simple structure to hold circle parameters
struct Circle_st {
  double r; // radius
  double xc; // center x coord
  double yc; // center y coord
  
  // return true if point x, y is inside circle
  bool is_inside( double x, double y ) const;
  
  Circle_st( double ar, double axc, double ayc ) : r(ar), xc(axc), yc(ayc) { }
  Circle_st() : r(1.0), xc(0.0), yc(0.0) { }
};

/// Take an input TH2D histogram, and compute the gradient (largest
/// change in value between 8 nearest neighbours) and fill that in the
/// gradhist that is built in this function.  Set any gradient values
/// that are below cut_below * maximum_gradient to zero, then fit to
/// circle. Returns Circle_st representing the best fit circle.
Circle_st find_circle_max_grad( const TH2D* in, TH2D*& grad, double cut_below = 0.8 );

/// Function sets any bins outside of circle to zero
void zero_outside_circle( TH2D* in, const Circle_st& c );

/// Function sets any bins inside of circle to zero
void zero_inside_circle( TH2D* in, const Circle_st& c );

#endif // __FINDCIRCLE__

