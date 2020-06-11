#ifndef __HOUGHDISPLAY__
#define __HOUGHDISPLAY__

#include "Hough.hpp"

/// Event display for HoughResults
/// Grabs hit pixels and hough results from HoughResults object
/// Grabs dimensions of display from CircleHough binning
void hough_display( CircleHough & ch,
		    const HoughResults & hcr );

#endif // __HOUGHDISPLAY__
