/* A basic R wrapper for the Munkres implementation of the Hungarian
   algorithm (John Weaver) */

// To be compiled by R CMD SHLIB munkres_Rwrap.cpp munkres.cpp
#include "munkres.h"

extern "C" {
  void munkres_Rwrap(int *nrows, int *ncols, double *dist) {
    Matrix<double> m(*nrows,*ncols);
    for ( int row = 0 ; row < *nrows ; row++ ) {
      for ( int col = 0 ; col < *ncols ; col++ ) {
	m(row,col) = dist[row+col*(*nrows)];
      }
    }
    // Apply Munkres algorithm to matrix.
    Munkres h;
    h.solve(m);
    for ( int row = 0 ; row < *nrows ; row++ ) {
      for ( int col = 0 ; col < *ncols ; col++ ) {
	dist[row+col*(*nrows)] = m(row,col);
      }
    }
  }
}


