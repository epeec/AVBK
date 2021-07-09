/*
 * Gather-Scatter Mini-App
 * To compile : g++ -I. -I/Users/mohanamuraly/NutsCFD_sandbox/NutsCFD/external/DIST/include -L/Users/mohanamuraly/NutsCFD_sandbox/NutsCFD/external/DIST/lib main.cpp -std=c++11 -lhdf5
 *
*/

#include <hdf5.h>
#include <iostream>
#include <string>
#include <vector>
#include <cassert>
#include <sstream>

#include "backend.h"
#include "kernel.cpp"

int main(int nargs, char *args[] ) {
  if( nargs != 2 )
    abort();

  // Read mesh data using C API call
  read_data(args[1], 20);
  int nnodes, ncells, npedges;
  int *triangles, *pedges;
  double *x, *y;
  get_data_ptr(&nnodes, &ncells, &npedges, &x, &y, &triangles, &pedges);

  // Allocate data for kernel
  double *vol_e = new double[ncells];
  double *vol_n = new double[nnodes];
 
  // Run the kernel
  element_volumes(nnodes, ncells, x, y, triangles, vol_e, vol_n);
  apply_periodic_bc(npedges, pedges, vol_n);

//  for( size_t i=0; i<nnodes; ++i)
//    std::cout << "Node " << i + 1 << " Volume = " << vol_n[i] << "\n"; 

  // Clean up the memory
  delete [] vol_e;
  delete [] vol_n;
  free_data();
  return 0;
}

