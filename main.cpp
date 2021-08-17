/*
 *  @file main.cpp
 *  @brief Gather-Scatter Mini-App
 *  @author Pavanakumar Mohanamuraly
 */
#include "backend.h"
#include "kernel.h"
#include <iostream>

/**
 *
 * @param nargs
 * @param args
 * @return
 */
int main(int nargs, char *args[]) {

  if (nargs != 2) {
    std::cout << "Need just one argument <input-hip-hdf5-mesh-file>\n";
    abort();
  }

  // Read mesh data using C API call
  read_data(args[1], 20);
  int nnodes, ncells, npedges;
  int *triangles, *pedges;
  double *x, *y;
  get_data_ptr(&nnodes, &ncells, &npedges, &x, &y, &triangles, &pedges);

  // Allocate data for kernel
  double *vol_e = new double[ncells];
  double *vol_n = new double[nnodes];
  double *work = new double[2 * npedges];

  for (unsigned i = 0; i < ncells; ++i)
    vol_e[i] = 0.0;
  for (unsigned i = 0; i < nnodes; ++i)
    vol_n[i] = 0.0;
  for (unsigned i = 0; i < 2 * npedges; ++i)
    work[i] = 0.0;

  // Run the main kernel
  run_kernel(ncells, npedges, x, y, triangles, pedges, work, vol_e, vol_n);

  for (size_t i = 0; i < nnodes; ++i)
    std::cout << "   (" << i << "): " << vol_n[i] << ",\n";

  // Clean up the memory
  free_data();
  delete[] vol_n;
  delete[] vol_e;
  delete[] work;
  return 0;
}
